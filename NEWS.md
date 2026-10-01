# shaman 2.1.0

* `shaman_shuffle_hic_track()` is 6-8x faster; results are unchanged for the same seed.
* The `sort_uniq` step of the shuffle (used by `shaman_shuffle_hic_track()`) counts contacts in memory instead of with `sort | uniq -c`: 44 to 3.5 min on chr1 of a Hi-C dataset. `.uniq` lines are now ordered by start1, start2.
* Scoring no longer needs perl or temporary files and is about 2x faster on dense regions; scores are unchanged.
* Scoring needs much less memory (17-25GB instead of 42-74GB for a dense 5Mb matrix) and can use several threads with `options(shaman.score.threads = N)`; scores are unchanged.
* `shaman_score_hic_mat_for_track()` accepts vectors of matrices; matrices in the same row share one read of each track.
* `shaman_shuffle_hic_track()`, `shaman_shuffle_hic_mat_for_track()` and `shaman_shuffle_and_score_hic_mat()` take a `seed` argument for reproducible shuffles. `shaman_shuffle_hic_track()` stores each chromosome's seed, time-based ones included, in the track attribute `seed`; `shaman_shuffle_hic_mat_for_track()` returns the seed it used.
* SGE mode (`options(shaman.sge_support = 1)`) works again with current misha.
* `shaman_shuffle_hic_track()` no longer fails after importing the track because the `shaman.debug` option was missing.
* The shuffler no longer writes one value past the end of a buffer when smoothing the decay curve; results are unchanged.
* `library(shaman)` no longer tries to unpack the example database; `shaman_get_test_track_db()` unpacks it into the user cache on first use (downloading it if needed), and current misha can open it.
* shaman now requires R >= 4.0.0 and a C++17 compiler (GCC >= 9).
