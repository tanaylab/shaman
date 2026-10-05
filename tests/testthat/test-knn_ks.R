# shaman_knn_ks_cpp() replaces RANN::nn2() + round() + shaman_merge_ks_cpp() in the score step.
# The kNN distances must be identical to RANN's, not just close.

# contacts-like integer coordinates: clustered near the diagonal, large values, duplicated points
contacts <- function(n, dup = 1000) {
    x <- round(runif(n, 38e6, 47e6))
    y <- x + round(rexp(n, 1 / 2e5)) * sample(c(-1, 1), n, replace = TRUE)
    list(x = c(x, rep(x[1:dup], 5)), y = c(y, rep(y[1:dup], 5)))
}

rann_dist <- function(dx, dy, qx, qy, k) RANN::nn2(cbind(dx, dy), cbind(qx, qy), k = k)$nn.dist

test_that("kNN distances are identical to RANN::nn2", {
    set.seed(1)
    p <- contacts(2e5)
    q <- sample(length(p$x), 2e4)
    for (k in c(1, 7, 100, 200)) {
        for (threads in c(1, 3)) {
            expect_identical(shaman_knn_dist_cpp(p$x, p$y, p$x[q], p$y[q], k, threads), rann_dist(p$x, p$y, p$x[q], p$y[q], k))
        }
    }
    # queries that are not data points
    expect_identical(shaman_knn_dist_cpp(p$x, p$y, p$x[q] + 0.5, p$y[q] - 3, 50), rann_dist(p$x, p$y, p$x[q] + 0.5, p$y[q] - 3, 50))
    # a small grid: every distance is tied many times; k = all points
    gx <- rep(0:20, 21)
    gy <- rep(0:20, each = 21)
    expect_identical(shaman_knn_dist_cpp(gx, gy, gx, gy, 50), rann_dist(gx, gy, gx, gy, 50))
    expect_identical(shaman_knn_dist_cpp(gx, gy, gx, gy, length(gx)), rann_dist(gx, gy, gx, gy, length(gx)))
    # all points identical
    expect_identical(shaman_knn_dist_cpp(rep(5, 100), rep(9, 100), c(5, 0), c(9, 0), 30), rann_dist(rep(5, 100), rep(9, 100), c(5, 0), c(9, 0), 30))
})

test_that("KS scores are identical to the RANN + round + shaman_merge_ks_cpp path", {
    set.seed(2)
    o <- contacts(1e5)
    e <- list(x = c(o$x, o$x + round(rnorm(length(o$x), 0, 5e3))), y = c(o$y, o$y))
    p <- sample(length(o$x), 5000)
    for (kk in list(c(100, 200), c(100, 137), c(7, 5), c(2, 2), c(3, 2))) {
        ref <- tryCatch(shaman_merge_ks_cpp(
            round(rann_dist(o$x, o$y, o$x[p], o$y[p], kk[1])),
            round(rann_dist(e$x, e$y, o$x[p], o$y[p], kk[2]))
        ), error = conditionMessage)
        for (threads in c(1, 4)) {
            # including the error for k of 2 or 3, where the perl sentinels fail the monotonicity check
            new <- tryCatch(shaman_knn_ks_cpp(o$x, o$y, e$x, e$y, o$x[p], o$y[p], kk[1], kk[2], threads), error = conditionMessage)
            expect_identical(new, ref)
        }
    }
    ref <- shaman_merge_ks_cpp(round(rann_dist(o$x, o$y, o$x[p], o$y[p], 100)), round(rann_dist(e$x, e$y, o$x[p], o$y[p], 200)))
    # both signs of the score occur
    expect_true(any(ref$V1 < 0) && any(ref$V2 > 0))
})

test_that("errors like RANN when k exceeds the number of points", {
    expect_error(shaman_knn_ks_cpp(1:50, 1:50, 1:150, 1:150, 1:3, 1:3, 100, 200, 1), "Cannot find more nearest neighbours than there are points")
    expect_error(RANN::nn2(cbind(1:150, 1:150), cbind(1:3, 1:3), k = 200), "Cannot find more nearest neighbours than there are points")
})
