# shaman_merge_ks_cpp() replaces the perl KS step (inst/perl/hic_merge_ks.pl). These tests run the old
# path (round, fwrite, perl, fread) on the same distance matrices and require identical output.

perl_ks <- function(o, e) {
    pl <- system.file("perl", "hic_merge_ks.pl", package = "shaman")
    fo <- tempfile("knn_o_")
    fe <- tempfile("knn_e_")
    fk <- tempfile("knn_ks_")
    ferr <- tempfile("knn_err_")
    on.exit(unlink(c(fo, fe, fk, ferr)))
    data.table::fwrite(as.data.frame(round(o)), fo, sep = "\t", col.names = F, quote = F, row.names = F)
    data.table::fwrite(as.data.frame(round(e)), fe, sep = "\t", col.names = F, quote = F, row.names = F)
    status <- system(sprintf("perl %s %s %s >%s 2>%s", pl, fo, fe, fk, ferr))
    lines <- length(readLines(fk))
    s_ks <- if (lines > 0) as.data.frame(data.table::fread(fk)) else NULL
    list(status = status, lines = lines, s_ks = s_ks, err = readLines(ferr))
}

ks_score <- function(s_ks) 100 * ifelse(-s_ks$V1 < s_ks$V2, s_ks$V2, s_ks$V1)

expect_same_as_perl <- function(o, e) {
    p <- perl_ks(o, e)
    if (p$status != 0) {
        # perl died at row p$lines + 1; the C++ version must stop at the same row with the same message
        err <- tryCatch(shaman_merge_ks_cpp(round(o), round(e)), error = conditionMessage)
        expect_equal(err, sprintf("%s (row %d)", p$err[1], p$lines + 1))
        return(invisible(p))
    }
    s_ks <- shaman_merge_ks_cpp(round(o), round(e))
    expect_equal(p$lines, nrow(o))
    expect_true(identical(s_ks$V1, as.numeric(p$s_ks$V1), num.eq = FALSE))
    expect_true(identical(s_ks$V2, as.numeric(p$s_ks$V2), num.eq = FALSE))
    expect_true(identical(ks_score(s_ks), ks_score(p$s_ks), num.eq = FALSE))
    invisible(p)
}

# rows sorted ascending like RANN::nn2 distances; self distance 0 in column 1 of the observed matrix
knn_mat <- function(n, k, values, self = FALSE) {
    m <- matrix(sample(values, n * k, replace = TRUE), n, k)
    m <- t(apply(m, 1, sort))
    if (self) m[, 1] <- 0
    m
}

skip_if_no_perl <- function() skip_if(Sys.which("perl") == "", "perl not available")

test_that("random distances match perl, k = 100, k_exp = 200", {
    skip_if_no_perl()
    set.seed(1)
    o <- knn_mat(3000, 100, runif(1e5, 0, 5e4), self = TRUE)
    e <- knn_mat(3000, 200, runif(1e5, 0, 5e4))
    p <- expect_same_as_perl(o, e)
    # the data must exercise both signs of the score
    expect_true(any(p$s_ks$V1 < 0) && any(p$s_ks$V2 > 0))
})

test_that("heavy ties and exact .5 distances match perl", {
    skip_if_no_perl()
    set.seed(2)
    o <- knn_mat(3000, 100, seq(0, 30, by = 0.5), self = TRUE)
    e <- knn_mat(3000, 200, seq(0, 30, by = 0.5))
    expect_same_as_perl(o, e)
    # coarse values: many exact ties across the two sides
    o <- knn_mat(2000, 100, 0:5, self = TRUE)
    e <- knn_mat(2000, 200, 0:5)
    expect_same_as_perl(o, e)
})

test_that("k_exp != 2k matches perl", {
    skip_if_no_perl()
    set.seed(3)
    for (kk in list(c(100, 137), c(100, 100), c(50, 30), c(7, 5), c(5, 7), c(20, 400))) {
        o <- knn_mat(1500, kk[1], runif(1e4, 0, 1e4), self = TRUE)
        e <- knn_mat(1500, kk[2], runif(1e4, 0, 1e4))
        expect_same_as_perl(o, e)
    }
})

test_that("one side running out first matches perl", {
    skip_if_no_perl()
    set.seed(4)
    near <- knn_mat(500, 100, runif(1e3, 0, 10))
    far <- knn_mat(500, 200, runif(1e3, 100, 200))
    expect_same_as_perl(near, far)
    near <- knn_mat(500, 200, runif(1e3, 0, 10))
    far <- knn_mat(500, 100, runif(1e3, 100, 200))
    expect_same_as_perl(far, near)
    # mixed within one matrix
    expect_same_as_perl(rbind(near[, 1:100], far), rbind(far, near[, 1:100]))
})

test_that("constant and identical rows match perl", {
    skip_if_no_perl()
    set.seed(5)
    expect_same_as_perl(matrix(7, 10, 100), matrix(7, 10, 200))
    expect_same_as_perl(matrix(0, 10, 100), matrix(0, 10, 200))
    o <- knn_mat(500, 100, runif(1e3, 0, 1e3))
    expect_same_as_perl(o, o)
    expect_same_as_perl(o, o[, c(1:100, 1:100)])
})

test_that("large distances written in scientific notation by fwrite match perl", {
    skip_if_no_perl()
    set.seed(6)
    big <- c(1e5, 2e5, 1e6, 1.2e6, 3e6, 1e7, 1.5e8, 2e8)
    o <- knn_mat(1000, 100, c(big, runif(200, 0, 3e8)), self = TRUE)
    e <- knn_mat(1000, 200, c(big, runif(200, 0, 3e8)))
    expect_same_as_perl(o, e)
})

test_that("single row matches perl", {
    skip_if_no_perl()
    set.seed(7)
    expect_same_as_perl(knn_mat(1, 100, runif(100, 0, 100), self = TRUE), knn_mat(1, 200, runif(200, 0, 100)))
})

test_that("small k, where the perl sentinels reach the monotonicity checks, matches perl", {
    skip_if_no_perl()
    set.seed(8)
    for (kk in list(c(3, 3), c(3, 5), c(5, 3), c(4, 3), c(3, 4), c(3, 2), c(2, 3), c(2, 2))) {
        for (i in 1:20) {
            o <- knn_mat(30, kk[1], runif(50, 0, 20))
            e <- knn_mat(30, kk[2], runif(50, 0, 20))
            expect_same_as_perl(o, e)
        }
    }
})

test_that("errors where the perl path fails", {
    skip_if_no_perl()
    set.seed(9)
    o <- knn_mat(100, 100, runif(1e3, 0, 1e3))
    e <- knn_mat(100, 200, runif(1e3, 0, 1e3))
    # non monotonic rows: perl dies at the first one
    o_bad <- o
    o_bad[37, 2:3] <- c(500, 100)
    expect_equal(perl_ks(o_bad, e)$lines, 36)
    expect_error(shaman_merge_ks_cpp(round(o_bad), round(e)), "non monotonic distance sequence at x1! \\(row 37\\)")
    e_bad <- e
    e_bad[37, 2:3] <- c(500, 100)
    expect_error(shaman_merge_ks_cpp(round(o), round(e_bad)), "non monotonic distance sequence at x2! \\(row 37\\)")
    # k = 1: perl dies on 1/0
    expect_true(perl_ks(o[, 1, drop = FALSE], e)$status != 0)
    expect_error(shaman_merge_ks_cpp(round(o[, 1, drop = FALSE]), round(e)), "division by zero")
    # row count mismatch: the old code failed when assigning the score column
    expect_error(shaman_merge_ks_cpp(round(o), round(e[-1, ])), "different numbers of rows")
})

test_that("every possible 3-decimal output survives perl printing and fread unchanged", {
    skip_if_no_perl()
    # the statistic is int(x * 1000) / 1000 with |x| < 1; the C++ code returns that value directly,
    # the perl path printed it and read it back with fread
    f <- tempfile()
    on.exit(unlink(f))
    system(sprintf("perl -e 'for my $j (-1000..1000) { my $x = $j / 1000; print \"$x\\t$j\\n\" }' > %s", f))
    back <- data.table::fread(f)
    expect_equal(back$V2, -1000:1000)
    expect_true(identical(as.numeric(back$V1), (-1000:1000) / 1000, num.eq = FALSE))
})
