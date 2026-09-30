#include <Rcpp.h>
#include <vector>

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
    const int n1 = o_dist.ncol() - 1; // perl: $n1 = $#x1
    const int n2 = e_dist.ncol() - 1;
    if (n1 < 1 || n2 < 1) {
        Rcpp::stop("Illegal division by zero: k and k_exp must be at least 2");
    }
    const double c1 = 1.0 / n1;
    const double c2 = 1.0 / n2;
    const double *o = o_dist.begin();
    const double *e = e_dist.begin();
    std::vector<double> x1(n1 + 1), x2(n2 + 1);
    Rcpp::NumericVector v1(nrow), v2(nrow);

    for (R_xlen_t r = 0; r < nrow; r++) {
        for (int j = 0; j <= n1; j++) x1[j] = o[r + (R_xlen_t)j * nrow];
        for (int j = 0; j <= n2; j++) x2[j] = e[r + (R_xlen_t)j * nrow];

        // perl: $x1[$n1] = $x2[$n2]+1; $x2[$n2] = $x1[$n2]+1;
        // $x1[$n2] is undef (0) when n2 > n1. The loop never reads index n1 / n2, so these
        // only matter for the monotonicity check below when k or k_exp is 2 or 3.
        x1[n1] = x2[n2] + 1;
        x2[n2] = (n2 <= n1 ? x1[n2] : 0) + 1;
        // perl reads index 2 even when it is past the end (undef, i.e. 0)
        const double x1_2 = n1 >= 2 ? x1[2] : 0;
        const double x2_2 = n2 >= 2 ? x2[2] : 0;
        if (x2_2 < x2[1]) {
            Rcpp::stop("non monotonic distance sequence at x2! (row %d)", (long)(r + 1));
        }
        if (x1_2 < x1[1]) {
            Rcpp::stop("non monotonic distance sequence at x1! (row %d)", (long)(r + 1));
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
        v1[r] = (double)(long long)(min_d * 1000) / 1000;
        v2[r] = (double)(long long)(max_d * 1000) / 1000;
    }
    return Rcpp::List::create(Rcpp::Named("V1") = v1, Rcpp::Named("V2") = v2);
}
