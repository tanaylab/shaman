## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission. Possibly misspelled words in DESCRIPTION: Mendelson (an author's
  surname), et al. (in the citation) and enrichments are correct; Aparametric is the second A of the
  package name.
* installed size is about 5Mb, sub-directory libs about 3.5Mb: debug information of the compiled
  code; the stripped library is about 0.25Mb.

## Test environments

* local: AlmaLinux 8.10, R 4.4.1, gcc 13.3 with -Wall -pedantic
* GitHub Actions: macOS (arm64), R release; Ubuntu, R devel, release and oldrel-1
* GitHub Actions: R-hub's clang-asan container (R-devel with AddressSanitizer and
  UndefinedBehaviorSanitizer)

## Notes

* shaman is `OS_type: unix` because it depends on misha, which is unix-only.
* Examples, tests and vignettes use a small example database bundled with the package (built in
  `tempdir()`) and download nothing.
* The examples that download gene and ideogram tracks from UCSC run only in interactive sessions;
  the versions without them are in `\donttest{}` because loading Gviz (Suggests) takes several
  seconds.
