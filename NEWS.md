# HPAclusteR 1.1.0

This release removes the Python and Seurat dependencies that made the package
hard to install and broke continuous integration, adds the package's first test
suite, and fixes the bugs that suite uncovered.

## Breaking changes

* `hc_pca()` now stores a plain list in `uns$pca` with elements `scores`,
  `loadings`, `sdev`, `r2cum` and `method`, instead of a `pcaMethods::pcaRes`
  S4 object. Code that reached into the old object needs updating:
  `pca@R2cum` becomes `pca$r2cum`, `pcaMethods::sDev(pca)` becomes `pca$sdev`,
  and `pca@scores` becomes `pca$scores`.

* `hc_snn()` no longer takes a `similarity` argument. It was passed to Seurat as
  an Annoy metric, but Annoy is only used when neighbours are searched in a
  coordinate space; on the precomputed-distance path this package uses, the
  argument never had any effect. The metric is the one chosen in
  `hc_distance()`.

* `hc_auto_cluster()` replaces `snn_similarity` with `snn_prune`, which is
  forwarded to `hc_snn()`.

* UMAP embeddings are no longer numerically identical to those from earlier
  versions. The input, the algorithm and every parameter are unchanged: the
  embedding is still optimised from the same shared nearest neighbour graph,
  with the same `n_epochs`, `learning_rate`, `min_dist`, `spread`,
  `repulsion_strength`, `negative_sample_rate` and spectral initialisation, and
  `uwot` and `umap-learn` derive the same `a`/`b` curve parameters from
  `spread` and `min_dist` (agreeing to seven decimal places). The two
  implementations differ only in their random number generators and SGD
  internals.

  Measured on the bundled example data, embeddings from the two implementations
  agree with a Spearman correlation of 0.90 across all pairwise distances and a
  25-nearest-neighbour overlap of 0.77 -- which is within the 0.75-0.77 range
  that `uwot` produces against *itself* across different seeds. In other words,
  switching implementation moves the embedding no more than changing the seed
  does. Embeddings remain reproducible for a given `seed` within this version.

* Community detection now uses `igraph` rather than Seurat's bundled Louvain
  implementation. On the bundled example data the consensus clusterings from
  the two implementations agree with an adjusted Rand index of 0.96, but
  individual cluster labels and counts can shift slightly. If a resolution now
  falls outside the accepted range, adjust `resolution` or the new `min_k` and
  `max_k` arguments.

## Installation and continuous integration

* **Removed the Python dependency.** `hc_umap()` now uses
  `uwot::optimize_graph_layout()` instead of driving the Python `umap-learn`
  package through `reticulate` and `Seurat::RunUMAP()`. No Python interpreter,
  no `umap-learn`, and no pinned `numpy` version is required, and the
  workflow step that pinned `numpy` to 1.26.4 is gone.

* **Removed Seurat.** `hc_snn()` builds the shared nearest neighbour graph
  directly; the implementation reproduces `Seurat::FindNeighbors()` on a
  precomputed distance matrix exactly, and is verified against it. Clustering
  moved to `igraph`.

* **Removed the `concaveman` dependency**, and with it the `V8` and `sf`
  packages and the GDAL, GEOS and PROJ system libraries they need. Cluster
  outlines now come from a native implementation of the same algorithm.

* **`pcaMethods` moved from Imports to Suggests**, so a default installation no
  longer pulls in Bioconductor. It is only needed for
  `hc_pca(na_action = "nipals")`.

* **Dropped `factoextra`, `matrixStats`, `MASS`, `fpc`, `GGally`, `network` and
  `sna`**, replaced by base R or existing dependencies. `ggrepel` and
  `patchwork` moved from Suggests to Imports, since exported functions need
  them; `hc_classify()` and `hc_cluster_compare()` now work on a fresh install.

* Minimum versions are declared for every dependency, and the minimum R version
  is now 4.2.0 rather than 4.5.0, which excluded most users unnecessarily.

* `R-CMD-check` now runs on Windows in addition to macOS and Ubuntu, and covers
  R release, devel and oldrel. Both workflows can be triggered manually.

* The vignette's annotation section is skipped rather than failing when the
  Bioconductor annotation stack is unavailable.

## Bug fixes

* `hc_cluster_consensus()` can be run again on an object that already carries a
  clustering. It previously failed with "Join columns in `x` must be present in
  the data", because the new clustering collided with the existing `cluster`
  column. Stale `cluster` and `cluster_colors` columns are now replaced.

* `check_annotation_alignment()` returned `NULL` instead of `FALSE` when an
  object was misaligned, so `AnnDatR$validate()` failed with "argument is of
  length zero" rather than reporting the misalignment.

* TRRUST was read as if it had a header row. The file does not, so its first
  record was silently dropped and the column names were taken from whichever
  transcription factor sorted first. Reactome and TRRUST column names are now
  supplied explicitly.

* `hc_cluster_hulls()` read gene identifiers from a hard-coded `ensembl_id`
  column, and now uses the object's `obs_names_col`. `hc_umap()` had the same
  problem when adding its coordinate columns.

* `hc_cluster_hulls()` no longer emits degenerate one- and two-vertex polygons
  for very small landmasses, which could not be drawn.

* `hc_cluster_consensus()` reads gene names from `obs` rather than through the
  `obs_names` active binding. R6 active bindings do not survive serialisation,
  so an object restored from disk could report stale names.

* `hc_umap()` no longer duplicates `UMAP1` and `UMAP2` columns when re-run on
  an object that already has them.

* Failed annotation database downloads are now reported with a warning naming
  the databases and URLs, instead of silently yielding a shorter database list.

* Several `dplyr` joins had no `by` argument, producing "Joining with ..."
  messages and silently depending on incidental column names.

* Fixed the deprecated `size` argument of `element_line()` and a deprecated
  external-vector `tidyselect` usage.

## Missing values and performance

* `hc_pca()` gained an `na_action` argument. Missing values previously made
  `pcaMethods::pca()` fall back to the iterative NIPALS algorithm without any
  indication, which is orders of magnitude slower than an SVD. `hc_pca()` now
  reports how many values are missing and, by default, imputes them with their
  column mean so the fast path is taken. `na_action = "nipals"` reproduces the
  old behaviour exactly, `"omit"` drops affected genes, and `"error"` stops.

* On the bundled example data, the pipeline is substantially faster:

  | Step | Before | After |
  | --- | --- | --- |
  | `hc_pca()` | 1.40 s | 0.02 s |
  | `hc_cluster_consensus()` | 5.73 s | 2.75 s |
  | `hc_umap()` | 6.19 s | 0.79 s |
  | `hc_cluster_hulls()` | 2.46 s | 0.74 s |
  | `hc_cluster_stability()` | 4.19 s | 0.60 s |

* `hc_cluster_consensus()` builds the graph once and reuses it across seeds,
  rather than constructing a Seurat object for every one of the 100 runs.

* `hc_cluster_stability()` extracts the per-seed labels once instead of
  re-selecting and deframing the table inside each of the seed pairs.

* `hc_classify()` computes per-gene detection counts once instead of filtering
  the full long table once per gene, and preallocates its result vectors.

* `build_annotation_terms_tibble()` reads `proteinatlas.tsv` once and shares it
  across formatters. Ten of the supported databases are derived from that one
  file, which was previously re-read for each of them.

* `get_density()` builds its long-form grid directly rather than reshaping and
  re-joining an n-by-n matrix.

## Annotation databases

* All four source URLs were verified to resolve and serve the expected content.

* `hpa_version` now defaults to `NULL`, meaning the current Human Protein Atlas
  release, rather than being pinned to version 24. Passing a number still
  selects an archived release host.

* Database downloads use a 30 minute timeout; R's 60 second default was not
  enough for files of this size.

## New features

* `hc_cluster_consensus()` gained `min_k` and `max_k` arguments. The accepted
  range for the median number of communities was previously hard-coded to
  30-110, and produced a misleading error message when a resolution fell
  outside it.

* `hc_pca()`, `hc_kaisers_rule()` and `hc_auto_cluster()` gained a `verbose`
  argument, and progress reporting uses `message()` rather than `print()` so it
  can be suppressed.

* `hc_umap()` exposes `min_dist` and `spread`.

## Testing

* Added the package's first test suite: 401 expectations in 126 test blocks
  across 13 files, covering every exported function, the internal numerical and
  geometric helpers, the annotation database formatters, and end-to-end runs of
  the documented pipeline. The replacements for `pcaMethods`, `factoextra`,
  `MASS::kde2d`, `fpc::dbscan` and `Seurat::FindNeighbors` were each validated
  against the original implementation before those dependencies were dropped:
  scaling, distances, the density grid, DBSCAN and the SNN graph all reproduce
  the originals exactly.

* `R CMD check --as-cran` passes with no errors or warnings, and also passes
  with `_R_CHECK_DEPENDS_ONLY_=true`, which confirms the package works with only
  its declared dependencies installed.


# HPAclusteR 1.0.0

* Initial release.
