# shaman on CRAN: plan

Working plan for getting shaman onto CRAN. It lives on the local branch `cran-readiness` (worktree
`~/src/shaman-cran`) and is excluded from the package build (`.Rbuildignore`). Nothing here has been
pushed or submitted anywhere.

Status keys: **done** (committed on this branch), **decision** (blocked on an open decision),
**macOS** (needs a macOS build to verify), **upstream** (needs an open PR to merge first), **todo**.

## Base

`cran-readiness` starts at `origin/master` and merges the open PRs in the recommended order
(`git merge --no-ff`, no conflicts):

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

Integration base (after the last merge): `e48958d1d6c329825b20ee99bad13e4af74e263c`. Everything
after it is the CRAN work.

Not merged, because they were not on origin when this branch was made (2026-09-30): the
`smooth_vector` heap overflow + `-Wreorder` fix (branch `fix-smooth-vector-overflow`) and the
shuffler `seed` argument (branch `feat-shuffle-seed`, stacked on #9). Both are being written by
other agents; this branch does not re-implement them.

Rebuilding after the PRs merge: `git checkout -b cran-readiness-2 origin/master`, merge whatever
PRs are still open in the order above plus the two above, then `git cherry-pick e48958d..cran-readiness`
(the CRAN commits are small and independent). Expected conflicts: `R/shaman.R` examples and
`src/hic_matrix_shuffler.cpp` with `feat-shuffle-seed`, `src/ContactShuffler.cpp` with
`fix-smooth-vector-overflow`.

git-lfs is not installed here, so `inst/trackdb.tar.gz` is the 134-byte LFS pointer in this
checkout. All git commands used `-c filter.lfs.*=` and `core.hooksPath=/dev/null`.

## Open decisions

1. **Maintainer.** Set to Aviezer Lifshitz <aviezer.lifshitz@weizmann.ac.il> (`cre`), tentative
   ("I will be the maintainer, I think"). Netta Mendelson Cohen stays `aut`. The maintainer must
   answer CRAN email, including the automatic check-failure notices.
2. **Where the full 104 MB example DB goes** (GitHub release, Zenodo, or dropped). Either way the
   package no longer ships it, and examples/tests use a tiny bundled example (see below).
   Implemented: the recommended option, a download on explicit request only. The URL is one constant
   in `R/params.R`; it still points at the git-lfs file on master, which keeps working as long as
   the file stays in the repository. To swap: change the URL (release/Zenodo), or, to drop it,
   revert the commit that adds the download argument.
3. **Does #9 ship in the CRAN release?** Assumed yes (recommended). That makes the macOS/`std::pmr`
   item a blocker; its fix is one separate commit. If #9 does not ship, drop the #9 merge and
   that commit.
4. **Gviz to Suggests.** Implemented (recommended) as one separate commit; `git revert` it to keep
   Gviz in Imports.
5. Also for the user: the license stays `GPL` (unversioned, which CRAN accepts); a versioned licence
   (`GPL-3`, ...) is the authors' call. A method reference with a DOI in Description is a
   nice-to-have that needs the right citation from the authors.

## Blockers

| # | item | status | effort |
|---|---|---|---|
| 1 | 104 MB example DB in the package; examples download it into the user cache | todo | 1-1.5 d |
| 2 | examples fail (misha not attached; the DB needs a download; UCSC/Gviz network) | todo | 0.5 d |
| 3 | NAMESPACE imports, undocumented/mismatched Rd arguments, missing `\value` | todo | 0.5 d |
| 4 | compiled code: `std::cerr`/`std::cout`, terminating `ASSERT` | todo | 1-2 h |
| 4b | compiled code: `-Wreorder`, `smooth_vector` overflow | upstream (fix-smooth-vector-overflow) | - |
| 5 | macOS build of #9 (`std::pmr` availability) | todo, then macOS | 0.5 d |
| 6 | DESCRIPTION (Title, Description, Date, Remotes, URL, OS_type, parallel, maintainer) | todo | 1 h |
| 7 | policy: `.GlobalEnv` assign, `options()`/`par()` restored, `T`/`F`, <= 2 threads in tests/examples | todo | 2-3 h |

## Nice-to-haves

| item | status |
|---|---|
| seed the C++ RNG from R (`seed` argument) | upstream (feat-shuffle-seed) |
| replace `system()` calls (`echo`, `rm` glob, `sleep`) with R code, quote paths | todo |
| replace doMC with `parallel::mclapply` (registerDoMC changes the user's backend) | not planned |
| delete dead code (`Parser`, `GenomeGridLog`) | not planned |
| drop `-pthread` from `src/Makevars` | not planned |
| sign-compare warnings | not planned |

## Check results

To be filled in.

## Results unchanged

To be filled in.
