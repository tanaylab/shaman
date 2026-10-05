## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission. Possibly misspelled words in DESCRIPTION: Mendelson, et al. (the cited
  authors) and enrichments are correct; Aparametric is the second A of the package name.
* installed size is 5.2Mb, sub-directory libs 3.5Mb: debug information of the compiled code; the
  stripped library is about 0.2Mb.

## Test environments

* local: AlmaLinux 8.10, R 4.4.1, gcc 13.3 (R's default flags, and -Wall -pedantic) and clang 18
* local: clang 18 with -fsanitize=address,undefined (examples, tests)
* GitHub Actions: macOS (arm64), R release; Ubuntu, R devel, release and oldrel-1
* GitHub Actions: R-hub's clang-asan container (R-devel with AddressSanitizer and
  UndefinedBehaviorSanitizer)

## Notes

* shaman is `OS_type: unix`: it depends on misha, which is unix-only.
* Examples and tests run on a small example database (0.5Mb of data in inst/extdata, built in
  `tempdir()`). The full example database (about 100Mb) is downloaded only when asked for with
  `shaman_get_test_track_db(full = TRUE)`, never in examples, tests or vignettes.
* No `\dontrun{}`. The gene and ideogram tracks of `shaman_plot_tracks_and_annotations()` and
  `shaman_plot_map_score_with_annotations()` are downloaded from the UCSC genome browser, so those
  examples run only in interactive sessions (`if (interactive())`). The versions without them run in
  `\donttest{}`, because loading Gviz (Suggests, Bioconductor) takes several seconds.
* Examples, tests and vignettes run without the Suggests packages Gviz and foreach (checked with
  only the hard dependencies, testthat, knitr and rmarkdown installed).
* Examples and tests use at most 2 cores.
