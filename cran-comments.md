## R CMD check results

0 errors | 0 warnings | 2 notes

* Days since last update: this update fixes the warning CRAN asked to correct before 2026-11-06.
* installed size is about 5Mb: debug information of the compiled code.

## This update

Fixes the g++ 16 `-Wstringop-overflow` warning in src/ContactShuffler.cpp. With g++ 16.2.0 and
16.2.1 and CRAN's flags, the package now compiles without warnings.
