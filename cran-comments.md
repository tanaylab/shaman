## R CMD check results

0 errors | 0 warnings | 2 notes

* Days since last update: this update fixes the warning CRAN asked to correct before 2026-11-06.
* installed size is about 5Mb, sub-directory libs about 3.5Mb: debug information of the compiled
  code; the stripped library is about 0.25Mb.

## This update

Fixes the `-Wstringop-overflow` warning of g++ 16 in src/ContactShuffler.cpp
(r-devel-linux-x86_64-fedora-gcc). With g++ 16.2.1 (the same version as on that machine) and the
same flags, the package compiles without warnings.

## Test environments

* local: AlmaLinux 8.10, R 4.4.1, gcc 13.3 with -Wall -pedantic
* R-hub's gcc16 container: Fedora 44, g++ 16.2.1, at -O2 and -O3
* GitHub Actions: macOS (arm64), R release; Ubuntu, R devel, release and oldrel-1
* GitHub Actions: R-hub's clang-asan container (R-devel with AddressSanitizer and
  UndefinedBehaviorSanitizer)

## Notes

* shaman is `OS_type: unix` because it depends on misha, which is unix-only.
