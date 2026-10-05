#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <thread>
#include <vector>

// KS statistic of one scored point, as inst/perl/hic_merge_ks.pl computes it (see
// shaman_merge_ks_cpp). x1: the point's k round()ed observed kNN distances (ascending), x2: its
// k_exp expected ones. Their last elements are overwritten (the perl sentinels). Writes the
// truncated min/max to v1/v2. Returns 0, or 1/2 when x2/x1 fails the perl monotonicity check.
static inline int ks_row(std::vector<double> &x1, std::vector<double> &x2, double *v1, double *v2) {
    const int n1 = (int)x1.size() - 1; // perl: $n1 = $#x1
    const int n2 = (int)x2.size() - 1;
    const double c1 = 1.0 / n1;
    const double c2 = 1.0 / n2;

    // perl: $x1[$n1] = $x2[$n2]+1; $x2[$n2] = $x1[$n2]+1;
    // $x1[$n2] is undef (0) when n2 > n1. The loop never reads index n1 / n2, so these
    // only matter for the monotonicity check below when k or k_exp is 2 or 3.
    x1[n1] = x2[n2] + 1;
    x2[n2] = (n2 <= n1 ? x1[n2] : 0) + 1;
    // perl reads index 2 even when it is past the end (undef, i.e. 0)
    const double x1_2 = n1 >= 2 ? x1[2] : 0;
    const double x2_2 = n2 >= 2 ? x2[2] : 0;
    if (x2_2 < x2[1]) {
        return 1;
    }
    if (x1_2 < x1[1]) {
        return 2;
    }

    int i1 = 1, i2 = 1;
    double d = 0, min_d = 0, max_d = 0;
    while (i1 < n1 && i2 < n2) {
        if (x1[i1] < x2[i2]) {
            d += c1;
            i1++;
        } else {
            d -= c2;
            i2++;
        }
        if (d > max_d) {
            max_d = d;
        } else if (d < min_d) {
            min_d = d;
        }
    }
    // perl: int($x*1000)/1000. int() truncates toward zero and yields an integer, so a
    // small negative min gives +0, not -0.
    *v1 = (double)(long long)(min_d * 1000) / 1000;
    *v2 = (double)(long long)(max_d * 1000) / 1000;
    return 0;
}

// In-memory port of inst/perl/hic_merge_ks.pl, which shaman used to run on the kNN distance
// matrices written to text files. Output is identical to the perl script's output as read back
// with data.table::fread, including the script's quirks:
// - walks start at column 1 (column 0 is skipped on both sides),
// - steps are 1/(k-1) and the last column is never compared (loop runs while i < k-1),
// - the walk stops as soon as either side runs out,
// - ties go to the expected side (strict <),
// - min/max are truncated toward zero to 3 decimals.
//
// o_dist, e_dist: round()ed distance matrices (one row per scored point, k and k_exp columns).
// Returns list(V1 = truncated min, V2 = truncated max), one value per row.
// [[Rcpp::export]]
Rcpp::List shaman_merge_ks_cpp(Rcpp::NumericMatrix o_dist, Rcpp::NumericMatrix e_dist) {
    const R_xlen_t nrow = o_dist.nrow();
    if (e_dist.nrow() != nrow) {
        // the perl script skips one missing line of its second file and dies on the second
        Rcpp::stop("observed and expected distance matrices have different numbers of rows (%d, %d)",
                   (long)nrow, (long)e_dist.nrow());
    }
    const int n1 = o_dist.ncol() - 1;
    const int n2 = e_dist.ncol() - 1;
    if (n1 < 1 || n2 < 1) {
        Rcpp::stop("Illegal division by zero: k and k_exp must be at least 2");
    }
    const double *o = o_dist.begin();
    const double *e = e_dist.begin();
    std::vector<double> x1(n1 + 1), x2(n2 + 1);
    Rcpp::NumericVector v1(nrow), v2(nrow);

    for (R_xlen_t r = 0; r < nrow; r++) {
        for (int j = 0; j <= n1; j++) x1[j] = o[r + (R_xlen_t)j * nrow];
        for (int j = 0; j <= n2; j++) x2[j] = e[r + (R_xlen_t)j * nrow];
        const int err = ks_row(x1, x2, &v1[r], &v2[r]);
        if (err) {
            Rcpp::stop("non monotonic distance sequence at %s! (row %d)", err == 1 ? "x2" : "x1", (long)(r + 1));
        }
    }
    return Rcpp::List::create(Rcpp::Named("V1") = v1, Rcpp::Named("V2") = v2);
}

// Exact 2D k-nearest-neighbor distances, a replacement for RANN::nn2(data, query, k)$nn.dist
// (exact search, eps = 0) that does not build the n x k index and distance matrices and can use
// threads.
//
// The values are identical to RANN's: both compute a squared distance as dx*dx + dy*dy in double
// (ANN adds the squared coordinate differences in dimension order) and report its square root, and
// the k smallest values of a multiset do not depend on which of several tied points are picked.
// Queries are independent, so splitting them across threads changes nothing.
namespace {

struct Pt {
    double x, y;
};

class KdTree {
  public:
    KdTree(const double *x, const double *y, R_xlen_t n) : m_pts(n) {
        for (R_xlen_t i = 0; i < n; i++) m_pts[i] = Pt{x[i], y[i]};
        if (n > 0) build(0, n);
    }

    // squared distances of the k nearest points to (qx, qy), ascending, into out[0..k)
    void knn(double qx, double qy, int k, std::vector<double> &heap, std::vector<int> &stack, double *out) const {
        heap.clear();
        stack.clear();
        stack.push_back(0);
        while (!stack.empty()) {
            const Node &nd = m_nodes[stack.back()];
            stack.pop_back();
            // no point in this box can replace the current k-th distance
            if ((int)heap.size() == k && box_dist2(nd, qx, qy) >= heap.front()) continue;
            if (nd.left < 0) {
                for (R_xlen_t i = nd.lo; i < nd.hi; i++) {
                    const double dx = qx - m_pts[i].x;
                    const double dy = qy - m_pts[i].y;
                    const double d2 = dx * dx + dy * dy;
                    if ((int)heap.size() < k) {
                        heap.push_back(d2);
                        std::push_heap(heap.begin(), heap.end());
                    } else if (d2 < heap.front()) {
                        std::pop_heap(heap.begin(), heap.end());
                        heap.back() = d2;
                        std::push_heap(heap.begin(), heap.end());
                    }
                }
                continue;
            }
            // visit the nearer child first (pushed last)
            if (box_dist2(m_nodes[nd.left], qx, qy) <= box_dist2(m_nodes[nd.right], qx, qy)) {
                stack.push_back(nd.right);
                stack.push_back(nd.left);
            } else {
                stack.push_back(nd.left);
                stack.push_back(nd.right);
            }
        }
        std::sort_heap(heap.begin(), heap.end());
        std::copy(heap.begin(), heap.end(), out);
    }

  private:
    struct Node {
        double x0, x1, y0, y1; // bounding box of the node's points
        R_xlen_t lo, hi;       // the node's points are m_pts[lo, hi)
        int left, right;       // children, -1 for a leaf
    };
    static const R_xlen_t LEAF_SIZE = 16;

    std::vector<Pt> m_pts;
    std::vector<Node> m_nodes;

    static double box_dist2(const Node &nd, double qx, double qy) {
        const double dx = qx < nd.x0 ? nd.x0 - qx : (qx > nd.x1 ? qx - nd.x1 : 0);
        const double dy = qy < nd.y0 ? nd.y0 - qy : (qy > nd.y1 ? qy - nd.y1 : 0);
        return dx * dx + dy * dy;
    }

    int build(R_xlen_t lo, R_xlen_t hi) {
        Node nd;
        nd.lo = lo;
        nd.hi = hi;
        nd.left = nd.right = -1;
        nd.x0 = nd.x1 = m_pts[lo].x;
        nd.y0 = nd.y1 = m_pts[lo].y;
        for (R_xlen_t i = lo + 1; i < hi; i++) {
            nd.x0 = std::min(nd.x0, m_pts[i].x);
            nd.x1 = std::max(nd.x1, m_pts[i].x);
            nd.y0 = std::min(nd.y0, m_pts[i].y);
            nd.y1 = std::max(nd.y1, m_pts[i].y);
        }
        const int id = (int)m_nodes.size();
        m_nodes.push_back(nd);
        if (hi - lo <= LEAF_SIZE || (nd.x0 == nd.x1 && nd.y0 == nd.y1)) return id;
        // split the wider side at the median
        const R_xlen_t mid = lo + (hi - lo) / 2;
        if (nd.x1 - nd.x0 >= nd.y1 - nd.y0) {
            std::nth_element(m_pts.begin() + lo, m_pts.begin() + mid, m_pts.begin() + hi,
                             [](const Pt &a, const Pt &b) { return a.x < b.x; });
        } else {
            std::nth_element(m_pts.begin() + lo, m_pts.begin() + mid, m_pts.begin() + hi,
                             [](const Pt &a, const Pt &b) { return a.y < b.y; });
        }
        const int left = build(lo, mid);
        const int right = build(mid, hi);
        m_nodes[id].left = left;
        m_nodes[id].right = right;
        return id;
    }
};

// f(i) for i in [0, n), on `threads` threads (no R API calls allowed in f)
template <class F> void parallel_for(R_xlen_t n, int threads, F f) {
    const R_xlen_t block = 1024;
    std::atomic<R_xlen_t> next(0);
    auto work = [&]() {
        for (R_xlen_t b = next.fetch_add(block); b < n; b = next.fetch_add(block)) {
            const R_xlen_t e = std::min(n, b + block);
            for (R_xlen_t i = b; i < e; i++) f(i);
        }
    };
    if (threads <= 1) {
        work();
        return;
    }
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) pool.emplace_back(work);
    for (auto &th : pool) th.join();
}

void check_knn_args(R_xlen_t n, int k) {
    if (k < 1) Rcpp::stop("k must be at least 1");
    if (k > n) Rcpp::stop("Cannot find more nearest neighbours than there are points"); // RANN::nn2's check
}

} // namespace

// shaman_merge_ks_cpp(round(RANN::nn2(cbind(obs_x, obs_y), pts, k)$nn.dist),
//                     round(RANN::nn2(cbind(exp_x, exp_y), pts, k_exp)$nn.dist))
// with pts = cbind(pts_x, pts_y), computed one point at a time on `threads` threads.
// [[Rcpp::export]]
Rcpp::List shaman_knn_ks_cpp(Rcpp::NumericVector obs_x, Rcpp::NumericVector obs_y, Rcpp::NumericVector exp_x,
                             Rcpp::NumericVector exp_y, Rcpp::NumericVector pts_x, Rcpp::NumericVector pts_y, int k,
                             int k_exp, int threads) {
    check_knn_args(exp_x.size(), k_exp);
    check_knn_args(obs_x.size(), k);
    if (k < 2 || k_exp < 2) {
        Rcpp::stop("Illegal division by zero: k and k_exp must be at least 2");
    }
    const KdTree e_tree(exp_x.begin(), exp_y.begin(), exp_x.size());
    const KdTree o_tree(obs_x.begin(), obs_y.begin(), obs_x.size());
    const R_xlen_t n = pts_x.size();
    const double *px = pts_x.begin();
    const double *py = pts_y.begin();
    std::vector<double> v1(n), v2(n);
    std::vector<char> err(n);

    parallel_for(n, threads, [&](R_xlen_t i) {
        thread_local std::vector<double> heap, x1, x2;
        thread_local std::vector<int> stack;
        x1.resize(k);
        x2.resize(k_exp);
        o_tree.knn(px[i], py[i], k, heap, stack, x1.data());
        e_tree.knn(px[i], py[i], k_exp, heap, stack, x2.data());
        for (double &d : x1) d = std::nearbyint(std::sqrt(d)); // R's round(x) is nearbyint(x)
        for (double &d : x2) d = std::nearbyint(std::sqrt(d));
        err[i] = ks_row(x1, x2, &v1[i], &v2[i]);
    });
    // only possible with k or k_exp of 2 or 3, where the perl sentinels reach the check
    for (R_xlen_t i = 0; i < n; i++) {
        if (err[i]) Rcpp::stop("non monotonic distance sequence at %s! (row %d)", err[i] == 1 ? "x2" : "x1", (long)(i + 1));
    }
    return Rcpp::List::create(Rcpp::Named("V1") = Rcpp::wrap(v1), Rcpp::Named("V2") = Rcpp::wrap(v2));
}

// RANN::nn2(cbind(data_x, data_y), cbind(query_x, query_y), k)$nn.dist computed with the tree
// above; for testing it against RANN.
// [[Rcpp::export]]
Rcpp::NumericMatrix shaman_knn_dist_cpp(Rcpp::NumericVector data_x, Rcpp::NumericVector data_y,
                                        Rcpp::NumericVector query_x, Rcpp::NumericVector query_y, int k,
                                        int threads = 1) {
    check_knn_args(data_x.size(), k);
    const KdTree tree(data_x.begin(), data_y.begin(), data_x.size());
    const R_xlen_t n = query_x.size();
    const double *qx = query_x.begin();
    const double *qy = query_y.begin();
    std::vector<double> d(n * (R_xlen_t)k);
    parallel_for(n, threads, [&](R_xlen_t i) {
        thread_local std::vector<double> heap;
        thread_local std::vector<int> stack;
        tree.knn(qx[i], qy[i], k, heap, stack, &d[i * k]);
    });
    Rcpp::NumericMatrix out(n, k);
    for (R_xlen_t i = 0; i < n; i++)
        for (int j = 0; j < k; j++) out(i, j) = std::sqrt(d[i * k + j]);
    return out;
}
