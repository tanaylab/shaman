# shaman on CRAN: plan

Working plan for getting shaman onto CRAN. It lives on the local branch `cran-readiness` (worktree
`~/src/shaman-cran`) and is excluded from the package build (`.Rbuildignore`). Nothing here has been
pushed, uploaded or submitted anywhere.

Status keys: **done** (committed on this branch), **decision** (blocked on an open decision),
**macOS** (needs a macOS build to verify), **upstream** (needs an open PR to merge first),
**not done** (nice-to-have left for later).

## Base

`cran-readiness` starts at `origin/master` and merges the open PRs in the recommended order
(`git merge --no-ff`):

| PR | branch | head merged |
|---|---|---|
| - | origin/master | 5be5ce2f514c38f53cbda2497b73bca9bf3bd126 |
| #6 | fix-sge-mode | 09d39f36dfa942d66a7becbe25048c41babd4706 |
| #7 | fix-test-db-extraction | c4f63c8ec074180af8c794242df784896f5014e0 |
| #9 | perf-shuffle-fast | 13665cf27f733edf2b58faf5f440b46229283daf |
| #8 | perf-ks-in-memory | 17172b28994ed5fcce5fd0eeadeba5a9953e9a72 |
| #10 | perf-score-fast | 810556a07267a2ab0d68940d150fe198292f7d7a |
| #11 | release-2.1.0 | f171f43df1dfb05412064a3fa8af9da27279cf4a |
| #12 | perf-shuffle-uniq | 53e9f481a335e24f596d4220da75b22dfde89171 |
| #13 | fix-smooth-vector-overflow | c8a331805b73314e31d6ce86d7be33a32787be7f |
| #14 | feat-shuffle-seed | 949ff41b00f5b906f4a718f07cab8d422db83fa3 |

Integration base (the #14 merge): `ca5a8b64991d44f8b49dcf98e119d589420e3812`. Everything after it
(`git log ca5a8b6..cran-readiness`) is the CRAN work, as small independent commits.

Merge conflicts, both in #14 and both as its PR describes: in `R/shaman.R`
(`shaman_shuffle_hic_mat_for_track()`), #12's in-memory `sort_uniq` block is kept without the
`ret <- 1` line that #14 removes; in `src/RcppExports.cpp`, the registration table, regenerated
with `Rcpp::compileAttributes()` (Rcpp 1.0.13; the only change from the pre-#14 file is the `seed`
argument). #6-#12 merge to the same tree as the first build (e48958d); #13 merges cleanly.

History: the first build of this branch (base e48958d, without #13 and #14) is kept as the local
tag `pre-rebuild/cran-readiness` (6e9f54b). Its CRAN commits were replayed with
`git cherry-pick`. Conflicts: the docs commit (now 7509d36) next to #14's `@param seed` lines,
resolved by keeping #14's `@return` for `shaman_shuffle_hic_mat_for_track()` (the seed, or NA)
and adding the `@return` for `shaman_shuffle_hic_track()`, then regenerating the Rd files. The old
2975f17 ("returns 0 instead of failing when sort_uniq is FALSE") was dropped: #14 initializes
`ret <- NA`, which fixes the same failure. Everything else applied as is, including the macOS
commit next to #13's reordered initializer list.

Rebuilding again after the PRs merge on origin: branch from the new `origin/master`, merge the PRs
that are still open in the order above, then `git cherry-pick ca5a8b6..cran-readiness`.

git-lfs is not installed here, so `inst/trackdb.tar.gz` is the 134-byte LFS pointer in this
checkout. Merges and commits ran with `-c filter.lfs.*=` and `core.hooksPath=/dev/null`.

## Decisions

1. **Maintainer:** Aviezer Lifshitz <aviezer.lifshitz@weizmann.ac.il>, `aut, cre` (final). Netta
   Mendelson Cohen stays `aut`. The maintainer answers CRAN email, including the automatic notices
   when a check starts failing.
2. **Full 104 MB example DB:** not in the package (`.Rbuildignore`); examples and tests use the small
   bundled one. `shaman_get_test_track_db(full = TRUE)` downloads it on request and checks its md5
   (88552541e7bf346ef50187f1e42fc25b, the same file as `inst/trackdb.tar.gz`). Decided location:
   the lab's public S3 bucket, `https://misha-genome.s3.eu-west-1.amazonaws.com/shaman/trackdb.tar.gz`.
   Done: the URL in `.shaman_get_full_test_track_db()` (`R/params.R`) and the wording in `README.Rmd`/`README.md`,
   the vignette, the Rd details and NEWS point to that bucket. The object is public (HTTP 200, 108,400,483 bytes,
   md5 88552541e7bf346ef50187f1e42fc25b, byte-identical to the LFS tarball).
3. **#9 ships** in the CRAN release, so the macOS `std::pmr` fallback (3841bd2) stays. Not yet
   verified on macOS (see CI below).
4. **Gviz stays in Suggests** (36f5e53, `\donttest` in 81ba85e). The README installation section
   (e381a52) says how to install it from Bioconductor; `.shaman_check_gviz()`'s error gives the same
   `BiocManager::install("Gviz")` line.
5. Names approved: `full`, `inst/extdata/hoxd.tsv.xz`, `tempdir()/shaman_test_db`,
   `.shaman_get_full_test_track_db()`, `.shaman_check_gviz()`, and from #14 `seed` (argument and track
   attribute), the `seed = N` line, `.shaman_check_seed()`.
6. Still for the authors: the licence stays `GPL` (unversioned, which CRAN accepts); a method
   reference with a DOI in Description needs the right citation.

## Blockers

| # | item | status | commits |
|---|---|---|---|
| 1 | 104 MB example DB in the package; examples downloaded it into the user cache | done; S3 URL pending (decision 2) | aa7e16d, efaf775, 04f9458 |
| 2 | examples failed (misha not attached, DB download, UCSC network, all cores, options not restored) | done | 8b59cbc, 81ba85e |
| 3 | NAMESPACE imports (utils, grDevices, graphics, stats), Rd `\usage` mismatches, `\value` | done | 26692d4, 7509d36 |
| 4 | compiled code: `std::cerr`/`std::cout`, terminating `ASSERT` | done | 85287f6 |
| 4b | compiled code: `-Wreorder` (install WARNING with `-Wall`), `smooth_vector` heap overflow (ASAN) | done in #13 (merged in the base) | - |
| 5 | macOS build of #9 (`std::pmr` needs macOS 14 at run time; `MADV_HUGEPAGE` was already guarded by #9) | done in code; needs macOS to verify | 3841bd2 |
| 6 | DESCRIPTION (Title, Description typo, Date, Remotes, URL, OS_type, parallel, maintainer) | done | 3c73f0d |
| 7 | policy: `.GlobalEnv` assign, `options()`/`par()` restored with `on.exit`, `T`/`F`, <= 2 threads | done | a705aa5, 35c082a |

Details:

- **Example DB (1).** `shaman_get_test_track_db()` builds, on first call in a session, a
  one-chromosome misha DB in `tempdir()/shaman_test_db` from `inst/extdata/hoxd.tsv.xz` (497 KB):
  the `hic_obs`, `hic_exp` and `hic_score` contacts of chr2:176.5e06-177e06 of the full example DB
  (46,735 / 88,264 / 45,150 contacts with start1 < start2, stored in both orientations as in the
  full DB). Layout as in PR #7: `chrom_sizes.txt` with `2`, an empty `seq/chr2.seq`. Building
  takes 0.6 s; later calls reuse it. It restores the previous misha root (`gdb.info()`, misha >=
  5.3). The window has 4 convergent pairs of the bundled `ctcf_forward`/`ctcf_reverse` sites
  100-500 kb apart, used by the feature grid examples. The data file was made with:

  ```r
  library(misha)
  gsetroot(shaman_get_test_track_db(full = TRUE))
  reg <- gintervals.2d(2, 176.5e6, 177e6, 2, 176.5e6, 177e6)
  d <- do.call(rbind, lapply(c("hic_obs", "hic_exp", "hic_score"), function(track) {
      x <- gextract(track, reg, colnames = "value")
      x <- x[x$start1 < x$start2, ]
      x <- x[order(x$start1, x$start2), ]
      # the score track holds floats of one-decimal scores
      data.frame(track = track, start1 = x$start1, start2 = x$start2, value = round(x$value, 1))
  }))
  write.table(d, xzfile("inst/extdata/hoxd.tsv.xz", compression = 9), sep = "\t", quote = FALSE, row.names = FALSE)
  ```
- **Examples (2).** Each example attaches misha, works in chr2:176.5e06-177e06, uses at most 2 jobs
  (`max_jobs = 2`) and restores the options it sets. The UCSC gene/ideogram versions of the two
  Gviz functions are `\dontrun` (network); runnable versions without them are in `\donttest`
  because loading Gviz takes 8-10 s here. Scoring examples use smaller
  focus intervals, and `k = 20` for the score track, to stay under 5 s. `shaman_get_test_track_db()`
  got an example.
- **Also fixed on the way:** undefined variables in error messages of `shaman_score_hic_track()`
  and `init_params()` (26692d4); the vignette placeholder title (efaf775); stale
  `inst/doc` and `vignettes/shaman-package.html` removed (efaf775); the unused 1.5 MB
  `inst/extdata/example_data.txt` removed (21a2e02); `@docType package` replaced by `"_PACKAGE"`
  (roxygen2 7.3.2 regenerated all Rd files; defaults now print as `2e+06` etc.); NEWS (b105cc0).
- **#14 (`seed`).** Covered by 37253f2: the two shuffle examples pass `seed = 1`, the track example
  prints the `seed` attribute, `shaman_shuffle_hic_track()`'s `\value` mentions it, and NEWS has
  #14's line plus the new return value of `shaman_shuffle_hic_mat_for_track()` (the seed or NA,
  was 0 or 1). No new threads. The `seed = N` line goes through `Rcpp::Rcerr`, that is `REprintf`,
  R's console error stream, which is what Writing R Extensions asks compiled code to use; the
  compiled-code check flags only `std::cout`/`std::cerr`, `printf` and the like. I left it there
  rather than making it a `message()`: the seed of a time-seeded run is chosen inside the C++
  code, next to the other progress lines it prints the same way, and it is also returned.
- **macOS (5).** On `__APPLE__` the grid cells are plain `std::vector`s (no pool, no huge pages);
  elsewhere the code is unchanged. Verified here: the Linux objects are identical before and after
  the commit (`objdump -d` and `.rodata`/`.data` of all 7 objects, g++ 13.3 -O2); the fallback
  path, forced on Linux by turning the `#ifndef __APPLE__` into `#if 0` in a scratch copy, builds
  without `<memory_resource>` and gives byte-identical results (below). Not verified: an actual
  macOS build. That needs mac-builder or R-hub, which are external submissions and need your OK.
  From the headers I expect the rest of the C++17 code (`std::to_chars` for integers,
  `__builtin_prefetch`) to build with Apple clang, but that is inference.

## Done after the decisions (2026-10-01)

- **Random 1-10 s wait** (4abf34d): only SGE shuffle jobs wait (`Sys.sleep(sample(1:10, 1))` in the job
  command, staggering their start). Direct calls of `shaman_shuffle_hic_mat_for_track()` and
  multi-core mode no longer wait, and the `shaman_shuffle_hic_track()` example (now with `seed = 1`)
  runs in about a second outside `\donttest`.
- **doMC** (e90da78): multi-core mode used `doMC::registerDoMC()` + `plyr::ddply(.parallel = TRUE)`,
  which left doMC registered as the user's foreach backend. Now `parallel::mclapply()` (what doMC
  runs underneath, with the same defaults: prescheduled, `mc.set.seed = TRUE`); an error in a job
  still stops the run, as it did (checked: plyr + doMC stops with "task 2 failed"). doMC left
  Imports; parallel (base R) is in Imports. Chosen over saving and restoring the foreach backend,
  which needs foreach internals.
- **Full DB download** (04f9458): a failed download stops with the URL in the message, and the md5
  is checked.
- **README** (e381a52): CRAN and GitHub installation, Gviz from Bioconductor; the `biocLite`, GenomeInfoDb,
  remotes requirement and old tarball instructions are gone. `README.md` re-knit (`--wrap=none`,
  no smart quotes, so only the changed section differs).
- **Tests** (79bd8d8, `tests/testthat/test-example_db.R`, about 4 s, at most 2 processes): the small
  DB has the expected contacts; a shuffle keeps every contact end (the marginal coverage, exactly),
  the same seed gives the same file and the misha options are restored; `shaman_score_hic_mat()`
  gives a score in [-100, 100] for each focus point; the multi-core track functions leave the foreach
  backend as it was. foreach is in Suggests for that test.
- **pkgdown** (1788f7e): `url` added, all exported topics and the datasets in the reference index,
  the two config helpers `@keywords internal`; `pkgdown::check_pkgdown()` passes and the site builds
  locally (Bootstrap 3 is reported as deprecated).
- **cran-comments.md** (49483e1, in `.Rbuildignore`).
- **CI** (503c642): `.github/workflows/R-CMD-check.yaml` (as misha: macOS arm64 with R release,
  Ubuntu with R devel, release and oldrel-1; r-lib/actions installs misha from CRAN and Gviz from
  Bioconductor; a macOS step prints `sw_vers` and the number of `std::pmr` symbols in the built
  library), `sanitizers.yaml` (R CMD check in R-hub's clang-asan container, failing on any ASan or
  UBSan report), `pkgdown.yaml` (builds the site as an artifact, no Pages deployment). Triggers: push to
  master and cran-readiness, pull requests to master. misha's `style.yaml` was not copied: it commits
  restyled code back to the branch by itself. `actionlint` passes on all three. **Not run:**
  pushing `cran-readiness` to GitHub was refused by this session's permission check (reason given:
  "out-of-place publication"), so no workflow has run and macOS is still unverified. Needs a push
  by you, or your direct go-ahead.

## Full example: the usage article as an end-to-end test (2026-10-01)

The usage article now runs on the full example database. Its code is in
`vignettes/shaman-package.Rmd.orig`; `vignettes/precompute.R` knits it into
`vignettes/shaman-package.Rmd` (code, output, and the figures `vignettes/shaman-package-*.png`)
and then runs the checks (below). The shipped `.Rmd` has no code to run, so R CMD check and
pkgdown only render it. Both the `.orig` and the script are in `.Rbuildignore`. To regenerate by
hand: `Rscript vignettes/precompute.R` from the package root, with the branch installed (it
downloads the database into a new temporary directory unless `R_USER_CACHE_DIR` is set). Commit
the `.Rmd` and figures again when they change. Names (approved): `vignettes/precompute.R`,
`vignettes/shaman-package.Rmd.orig`, `.github/workflows/full-example.yaml`, the artifact
`full-example`, the cache keys `shaman-full-db-<md5>` and
`shaman-full-example-tracks-<image>-<R version>-<hash>-<md5>`, the figure files (knitr's
`shaman-package-<chunk label>-1.png`: `score-map`, `ctcf-grid`). The track names `hic_obs_shuffle`
and `hic_score_new` come from the old article.

**Measured** on n111 (all8.q, AMD EPYC 9384X), `taskset -c 0-3`, 16G requested, misha's options set
to what it picks on a 4-core 16GB machine (`gmax.processes = 2` etc.). Memory: peak of the summed PSS
of all processes (the forked workers share pages, so summed RSS overcounts); `/usr/bin/time` gives
only the largest single process. Disk: peak of the database, work dir and temp dir.

| step | wall | memory |
|---|---|---|
| download (S3) and extract the database | 11 s | 0.1 GB |
| `shaman_shuffle_hic_track()`, `max_jobs = 4`, `seed = 1` | 5 min 14 s | 1.0 GB |
| `shaman_score_hic_track()`, `max_jobs = 4` | 2 min 59 s | 2.3 GB |
| score map of chr2:175e06-178e06 (929,416 points) | 20 s | 1.4 GB |
| feature grid and its figure | 22 s | 1.5 GB |
| comparison with `hic_score` | 7 s | 0.4 GB |
| local mode, chr2:176.5e06-177e06 (not timed with the monitor) | 41 s | - |
| `precompute.R`, all of the above plus the `hic_exp` check, fresh download (first version) | 10 min 33 s | 3.1 GB, 2.2 GB disk |
| `precompute.R`, final version (adds the chr2 rerun and the coverage check), fresh download | 12 min 21 s | 3.0 GB, 2.2 GB disk |

**CI, not local**: well under the ~45 min / 12 GB budget of a 4-vCPU, 16 GB runner. A runner is
slower per core than this node; by how much is a guess (up to 2x), which still leaves a wide margin.
`full-example.yaml` runs `precompute.R` weekly (Monday 04:00 UTC), by hand, and from `pkgdown.yaml`
on pushes to master (before every deployment) and to cran-readiness (to test it on GitHub before
merging; no deployment from there). It caches the downloaded database (key: md5) and the two
computed tracks (key: runner image and R version, `src/`, `R/`, the `.orig`, `precompute.R`, md5;
the article holds the shuffle and score parameters; the image and R version because the seeded
rerun is compared with the cached shuffle, and another compiler or libm could change it), saved with `actions/cache/save` right after the
download, before the article adds tracks. Restored tracks are picked up by `gsetroot(rescan = TRUE)`
(misha's `.db.cache` otherwise hides them), and the article skips the shuffle and score when their
tracks exist. It uploads the regenerated `.Rmd` and figures; `pkgdown.yaml` builds the site from them on master
and cran-readiness (deploying from master only), and builds nothing if the run fails. PR builds use
the committed article. `actionlint` passes. **Not run on GitHub** (not pushed).

**Figure changes.** Both figures use the database's `hic_exp` (decided after the cut-out finding
below): the score map draws scores of `hic_obs` against `hic_exp` from `shaman_score_hic_mat()`
(r = 0.9996 with the bundled `hic_score` the old article plotted, so it is the 2017 map), and the
feature grid uses `exp_track_nm = "hic_exp"`. The article's shuffle and score still run, as the
pipeline test; the figures do not depend on them. Score map: no gene track and ideogram (UCSC; `Gviz::UcscTrack()` fails here: rtracklayer's UCSC
table browser query gets an empty response, "missing value where TRUE/FALSE needed" in `htmlParse`;
`IdeogramTrack()` works but needs the network); no `K562.k27ac`/`rna-seq` tracks (not in the
database); the CTCF annotations are the package datasets `ctcf_forward`/`ctcf_reverse` (data/*.rda;
the old code passed their names as strings, which are not interval sets of the database). Feature
grid: now on the example database (where the old one came from is not recorded; its file name is
`hic.K562.ela_k562.ctcf_neg.ctcf_pos.1k.png`), with the two panels its caption always described (log10 obs, log2 obs/exp; the old one had only the
second). The old call passed `grid` instead of `list(grid)`. Local mode now runs too (it prints a summary
of its scores; its map was dropped, since its expected is a shuffle of the cut-out). The article had
`warning = FALSE` on every chunk; now only the shuffle chunk has it: `shaman_shuffle_hic_track()` warns
"1 full chrom files were not shuffled" for chrY, which has no contacts in the database (and a run with
cached tracks skips the shuffle, so the article would differ between fresh and cached runs).

**Scores.** On the 464,708 points of chr2:175e06-178e06 in `hic_score` (2017):

| score | Pearson r with `hic_score` |
|---|---|
| `hic_score_new` (shuffle seed 1) | 0.378 (Spearman 0.303) |
| the same with seed 2 | 0.377; seed 1 vs seed 2: 0.979 |
| today's scoring of `hic_obs` against the bundled `hic_exp` | 0.9996 (track route and `shaman_score_hic_mat()`; 1.3% of points differ by > 0.1) |
| local mode on chr2:176.5e06-177e06 (45,150 points) | 0.891 |

So the scoring reproduces 2017; the low correlation is the expected track. Cause: the database is a
cut-out. Outside the hoxd window every contact is 200kb-2Mb long (min 200,001, max 1,999,957), around
convergent CTCF pairs (the convergent grid collects 3.16M contacts, the divergent one 0.55M). The
shuffle keeps each whole chromosome's distance distribution, so in the window it puts too few
expected contacts at < 100kb and too many at 200kb-1Mb: mean new score +53 at < 20kb (2017: -6), -33
at 200-500kb (2017: +7). Per distance bin in the window, the new expected has 0.9x the observed
contacts below 10kb and 3.7x at 200-500kb; the bundled `hic_exp` has 2.0-2.2x up to 500kb (1.5-1.9x
beyond). On chr2 (contacts with start1 < start2), `hic_exp` has 2.0x the observed inside the window
but 0.9x outside it, while the new shuffle, which keeps every contact end, has 2.0x and 2.1x: so
`hic_exp` is not a shuffle of this `hic_obs` (my guess: it was cut, like `hic_obs`, from a shuffle of
the complete data). Feature grid, log2 obs/exp at the center: 0.21 with the
new expected, 0.56 with `hic_exp`. The article has a short note on this.

**Checks in `precompute.R`** (all must pass): (1) scoring `hic_obs` against `hic_exp` gives r >= 0.999
with `hic_score` (0.9996); (2) chr2 shuffled again with the seed recorded in `hic_obs_shuffle` gives
the same contacts and counts (2,706,104 contacts, seed 12); (3) the shuffle keeps the marginal
coverage on every chromosome: per position, the shuffled contacts are exactly twice the contact ends
it shuffled (once on unshuffled chromosomes); checked to fail with `hic_exp` in place of the shuffle.
`hic_score_new` vs `hic_score` (0.378) is reported without a threshold.

**Check** (`R CMD check --as-cran` of 7c6bd32, as in Check results below): 0 ERRORs, 0 WARNINGs, the
same 3 NOTEs, no ggplot2 deprecation in the log. The tarball is 1.64 MB (1.77 MB before this work) and
the installed size 5.1 MB, with only `libs` over 1 MB. The shipped vignette rendered in 1.4 s (b8fa34a).

**ggplot2**: `shaman_plot_map.R` passed `size` to `element_line()`, deprecated since ggplot2 3.4.0
(a warning in every session's first map plot); now `linewidth`, with `ggplot2 (>= 3.4.0)` in
DESCRIPTION (a90038c). The 9 test PNGs of the map functions are byte-identical before and after.

## Nice-to-haves

| item | status |
|---|---|
| replace the other `system()` calls (`echo` header, `rm <work_dir><track>*` glob) with R code, quote paths | not done (fine on unix) |
| dead code (`Parser`, `GenomeGridLog`), `-pthread` in `src/Makevars`, `-Wsign-compare` and `-Wpedantic` "extra ;" warnings (not counted by R CMD check) | not done |
| pkgdown: Bootstrap 5, and deploying the site (there is a `gh-pages` branch) | not done, needs your OK for deployment |
| installed size: the 5.2 MB NOTE is 3.5 MB of `-g` debug info in `shaman.so` (209 KB stripped) | nothing to do |

## Check results

Run on cran-readiness 503c642 (base ca5a8b6, with #13 and #14; the earlier run on 37253f2 gave the
same results). R 4.4.1, misha 5.11.23 (CRAN),
`R CMD build` + `R CMD check --as-cran` of a `git archive` of the branch, with
`_R_CHECK_THINGS_IN_OTHER_DIRS_` and `_R_CHECK_THINGS_IN_TEMP_DIR_` on, a scratch HOME/TMPDIR,
pinned to 8 cores. TeX (TinyTeX), qpdf and tidy were installed into the scratchpad, so the PDF
manual, HTML manual and PDF size checks ran. For comparison, the first integration base e48958d
(`--no-examples`, no manual, no qpdf yet): 4 WARNINGs (one of them the missing qpdf), 9 NOTEs.

| build | result |
|---|---|
| gcc 13.3, R's default flags | 0 ERRORs, 0 WARNINGs, 3 NOTEs (1-3 below) |
| gcc 13.3, `-Wall -pedantic` (CRAN's Linux warning flags) | 0 ERRORs, 0 WARNINGs, 4 NOTEs (1-4 below); no significant compiler warnings |
| clang 18 + libstdc++, `-Wall -pedantic` | 0 ERRORs, 0 WARNINGs, 2 NOTEs (1, 3 below) |
| clang 18 + libc++ (conda), first build | builds; examples OK; tests crash in the first `Rcpp::stop()` - a toolchain problem here, a 3-line Rcpp function that calls `Rcpp::stop()` crashes the same way with this libc++ and works with libstdc++. Not re-run. |

Before #13 was merged (first build), `-Wall` gave 1 WARNING for `-Wreorder`; #13 fixed it.

The NOTEs:

1. New submission, maintainer. Expected.
2. Installed size 5.2 MB, of which `libs` 3.5 MB: `-g` debug info (the stripped `.so` is 209 KB).
   CRAN's Linux builds use `-g` too, so this NOTE is likely there as well; it is common for
   compiled packages. (clang's debug info is smaller and stays under 5 MB.)
3. "unable to verify current time": this machine cannot reach the time server. Environment.
4. "new files in some other directories": `/tmp/tmp*wandb-media`, `/tmp/tmp*wandb-artifacts`.
   Directories owned by you that some other process on the node creates and removes (they came
   and went between my looks; I did not pin down which process). wandb is a Python library,
   shaman has no Python, and the check ran with its own TMPDIR. Environment; it shows up in some
   runs and not others.

Sanitizers: clang 18, `-fsanitize=address,undefined`, libstdc++, run with the ASAN runtime
preloaded, on 503c642: all 16 example files including `\donttest`, the 17 tests, and the
assessment's synthetic shuffle + score script (one shuffle with seed 7, one time-seeded): no
AddressSanitizer or UndefinedBehaviorSanitizer report. (The first build, without #13, hit the
`smooth_vector` heap overflow that #13 fixes.)

## Results unchanged

Compared: the integration base ca5a8b6 ("base"), this branch at 503c642 ("new": after the sleep
move and the switch to `parallel::mclapply`), and 503c642 with the macOS branch forced on Linux by
turning `#ifndef __APPLE__` into `#if 0` ("fallback"). The same comparison on 37253f2 gave the same
results.
All three were installed with R's default flags (g++ 13.3, -g -O2) and run by the same script on
a hard-linked copy of the same database. `time()` was fixed with an `LD_PRELOAD` shim
(`FAKE_TIME=1700000000` unless noted), so time-seeded shuffles (`seed = NULL`) seed the same way
in every run. Every output file was compared with `cmp`; md5 manifests of all outputs are kept in
the scratchpad (`ident/manifests`).

Outputs of the script, on each database:

- `shaman_score_hic_mat()` (k = 100; and k = 50 with `k_exp = NA`), `shaman_score_hic_points()`,
  `shaman_kk_norm()`: returned objects (`saveRDS`).
- `shaman_score_hic_mat_for_track()`, one vectorized call for 3 matrices of a row: the 3 `.score`
  files and the return value.
- `shaman_score_hic_track()` in multi-core mode (2 jobs): `gextract()` of the new score track and
  the track's files.
- `shaman_shuffle_and_score_hic_mat()` (shuffle = 20, grid_step_iter = 10): the shuffled file and
  the scores.
- `shaman_shuffle_hic_track()` in multi-core mode (2 jobs, shuffle = 2, grid_step_iter = 1): every
  chromosome's `.full_chrom_shuffled` and `.uniq` files, the imported track's files, its
  `gextract()` and its `seed` attribute.
- `shaman_generate_feature_grid()` and `shaman_plot_feature_grid()` (matrices), `shaman_score_pal()`,
  and PNGs of `shaman_gplot_map()`, `shaman_gplot_map_score()` (plus its layer data),
  `shaman_plot_feature_grid()` and `shaman_plot_tracks_and_annotations()` without the UCSC tracks.

| comparison | full example DB (114 files) | small example DB (25 files) |
|---|---|---|
| base vs new, `seed = NULL` | all identical | all identical |
| base vs fallback, `seed = NULL` | all identical | all identical |
| base vs new, `seed = 7` | all identical | all identical |
| new `seed = 7` twice, `FAKE_TIME` 1700000000 and 1800000000 | all identical | all identical |
| new, `seed = NULL` vs `seed = 7` (the seed takes effect) | 73 differ (the shuffle outputs) | 7 differ |

The full example DB is hg19 with 4.6M contacts (tarball sha256 a2aba695..., the LFS oid); the
small one was built by the new `shaman_get_test_track_db()`, and the same copy was used for all
runs. The `seed` attributes were `chr1:7 chr10:8 ...` for seed 7 and `chr1:1700000000 ...` under
the shim.

The macOS commit (3841bd2) also leaves the Linux object code unchanged: `objdump -d` and
`.rodata`/`.data` of all 7 objects are identical between 36f5e53 and 3841bd2 (g++ 13.3 -O2).

The small DB itself: `gextract()` of its three tracks over chr2:176.5e06-177e06 equals
`gextract()` of the full DB over the same region (every row and value, including the float
scores), apart from the chromosome factor levels.

This compares the branch with the merged PRs. That the PRs match stock shaman is what the PRs
report; it is not re-checked here.

## Likely questions from a CRAN reviewer

- `\dontrun` in the two Gviz examples: justified by the network (UCSC); the assessment also saw
  `Gviz::UcscTrack()` fail with UCSC reachable.
- The C++ code prints progress (a few lines per iteration) with no way to turn it off; some
  reviewers ask for a `verbose` switch.
- `shaman_generate_feature_grid()` writes interval sets into the misha DB (cached
  `giterator.intervals()` results); this is not documented.
- A reference (DOI) for the method in Description.

## Remaining effort (estimate)

- Your decisions 1-4: minutes each; dropping the full DB or keeping Gviz in Imports is one commit.
- When the PRs merge on origin, rebuild this branch on the new master (cherry-pick
  ca5a8b6..cran-readiness), re-run the check and the identity script: 1-2 h.
- macOS: one mac-builder or R-hub run (needs your OK), plus fixes if it fails: 0.5 day.
- `cran-comments.md`, final check, submission, and a round of CRAN reviewer comments: 0.5-1 day
  spread over the review wait.
