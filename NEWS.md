# shaman 2.1.0

* `shaman_shuffle_hic_track()` is 6-8x faster; results are unchanged for the same seed.
* Scoring no longer needs perl or temporary files and is about 2x faster on dense regions; scores are unchanged.
* Scoring needs much less memory (17-25GB instead of 42-74GB for a dense 5Mb matrix) and can use several threads with `options(shaman.score.threads = N)`; scores are unchanged.
* `shaman_score_hic_mat_for_track()` accepts vectors of matrices; matrices in the same row share one read of each track.
* SGE mode (`options(shaman.sge_support = 1)`) works again with current misha.
* `shaman_shuffle_hic_track()` no longer fails after importing the track because the `shaman.debug` option was missing.
* `library(shaman)` no longer tries to unpack the example database; `shaman_get_test_track_db()` unpacks it into the user cache on first use (downloading it if needed), and current misha can open it.
* shaman now requires R >= 4.0.0 and a C++17 compiler (GCC >= 9).
