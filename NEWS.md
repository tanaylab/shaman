# shaman 2.1.0

* `shaman_shuffle_hic_track()` is about 6x faster; results are unchanged for the same seed.
* The `sort_uniq` step of the shuffle (used by `shaman_shuffle_hic_track()`) counts contacts in memory instead of with `sort | uniq -c`: 44 to 3.5 min on chr1 of a Hi-C dataset.
* Scoring no longer needs perl or temporary files and is about 2x faster on dense regions; scores are unchanged.
* Scoring needs less memory (peak virtual memory 16-36GB instead of 48-77GB with shaman 2.0 for dense 5Mb matrices) and can use several threads with `options(shaman.score.threads = N)`; scores are unchanged.
* `shaman_score_hic_mat_for_track()` accepts vectors of matrices; matrices in the same row share one read of each track.
* `shaman_shuffle_hic_track()`, `shaman_shuffle_hic_mat_for_track()` and `shaman_shuffle_and_score_hic_mat()` take a `seed` argument for reproducible shuffles. `shaman_shuffle_hic_track()` stores each chromosome's seed, time-based ones included, in the track attribute `seed`; `shaman_shuffle_hic_mat_for_track()` returns the seed it used (NA when it did not shuffle) instead of 0 or 1, and no longer fails with "object 'ret' not found" when `sort_uniq = FALSE` and the matrix is not shuffled; `shaman_shuffle_and_score_hic_mat()` returns the seed as `seed`.
* SGE mode (`options(shaman.sge_support = 1)`) works again with current misha.
* `shaman_shuffle_hic_track()` no longer fails after importing the track because the `shaman.debug` option was missing.
* The shuffler no longer writes one value past the end of a buffer when smoothing the decay curve; results are unchanged.
* The shuffler no longer kills the R session (division by zero) on fewer than 10,000 contacts: a small chromosome or contig in `shaman_shuffle_hic_track()`, or a sparse region in `shaman_shuffle_and_score_hic_mat()`.
* The shuffler stops with an error when it cannot write its output (e.g. a full disk) instead of leaving a truncated file.
* `shaman_score_hic_track()` stops with an error naming the matrices that still have no score after 3 rounds, instead of resubmitting them forever; the finished scores stay in `work_dir`, so a rerun computes only the missing ones.
* A matrix whose expected has fewer contacts than `k_exp` gets no scores (as with fewer than `k`) instead of failing with "Cannot find more nearest neighbours than there are points".
* `library(shaman)` no longer tries to unpack the example database; `shaman_get_test_track_db()` unpacks it into the user cache on first use (downloading it from the lab's public S3 bucket and checking its md5 if needed), and current misha can open it.
* The unused `shaman.ks_pl` option was removed from `shaman.conf`; configuration files that still set it load as before.
* shaman now requires R >= 4.0.0 and a C++17 compiler.
* `shaman_get_test_track_db()` returns a small example database (chr2:176.5e06-177e06 of the previous one), built in the session's temporary directory; `shaman_get_test_track_db(full = TRUE)` returns the full one, downloaded from the lab's public S3 bucket on first use. The package no longer includes the 104MB database, and the examples run on the small one.
* Gviz is suggested instead of required; `shaman_plot_tracks_and_annotations()` and `shaman_plot_map_score_with_annotations()` need it.
* The shuffle, score and feature grid functions restore the misha options they set (`gmultitasking`, `gmax.data.size`), and `shaman_plot_feature_grid()` restores `par()`.
* `shaman_shuffle_hic_track()`, `shaman_shuffle_hic_mat_for_track()` and `shaman_shuffle_and_score_hic_mat()` take a `seed` argument; `shaman_shuffle_hic_track()` records each chromosome's seed, time-based ones included, in the track attribute `seed`. `shaman_shuffle_hic_mat_for_track()` returns the seed (NA when it did not shuffle) instead of 0 or 1, and no longer fails with "object 'ret' not found" when `sort_uniq = FALSE` and the matrix is not shuffled.
* The shuffler writes its progress through R's console instead of directly to stderr.
* The usage article runs on the full example database, and its figures are made by the code it shows (precomputed, so the vignette builds without it).
