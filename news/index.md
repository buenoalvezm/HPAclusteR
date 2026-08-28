# Changelog

## HPAclusteR 1.1.0

### Breaking changes

- [`hc_pca()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_pca.md)
  now stores PCA results as a list rather than a
  [`pcaMethods::pcaRes`](https://rdrr.io/pkg/pcaMethods/man/pcaRes.html)
  object. Use `$scores`, `$loadings`, `$sdev`, `$r2cum`, and `$method`
  to access the results.

- [`hc_snn()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_snn.md)
  no longer has a `similarity` argument. The distance metric is
  determined by
  [`hc_distance()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_distance.md).

- [`hc_auto_cluster()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_auto_cluster.md)
  now uses `snn_prune` instead of `snn_similarity`.

- UMAP embeddings may differ from previous versions due to the switch
  from `umap-learn` to `uwot`. Results remain reproducible for a given
  `seed`. The two implementations showed a Spearman correlation of 0.90
  for pairwise distances and 0.77 overlap among 25-nearest neighbours,
  comparable to the variation between uwot runs with different seeds
  (~0.75).

- Community detection now uses `igraph`. Cluster labels and counts may
  differ slightly from previous versions. The resulting clusterings
  showed an Adjusted Rand Index (ARI) of 0.96 against the previous
  implementation on the bundled example data, although individual
  cluster labels and counts may differ slightly.

### New features

- [`hc_pca()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_pca.md)
  gains `na_action` for handling missing values. The default uses
  column-mean imputation; `"nipals"`, `"omit"`, and `"error"` are also
  available.

- [`hc_cluster_consensus()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_consensus.md)
  gains `min_k` and `max_k` to control the accepted range of cluster
  numbers.

- [`hc_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_umap.md)
  exposes the `min_dist` and `spread` parameters.

- [`hc_pca()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_pca.md),
  [`hc_kaisers_rule()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_kaisers_rule.md),
  and
  [`hc_auto_cluster()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_auto_cluster.md)
  gain a `verbose` argument.

- `hpa_version` now defaults to `NULL`, using the current Human Protein
  Atlas release. A numeric version can still be used to select an
  archived release.

### Dependencies

- Removed the Python, Seurat, `concaveman`, and associated system
  dependencies.

- `pcaMethods` is now a suggested dependency and is only required for
  `hc_pca(na_action = "nipals")`.

- Removed several unused dependencies, including `factoextra`,
  `matrixStats`, `MASS`, `fpc`, `GGally`, `network`, and `sna`.

- Minimum R version is now 4.2.0.

### Bug fixes

- [`hc_cluster_consensus()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_consensus.md)
  can now be rerun on objects that already contain clustering results.

- Fixed annotation alignment validation returning `NULL` instead of
  `FALSE`.

- Fixed loss of the first TRRUST record when reading annotation data.

- Fixed use of hard-coded gene identifier columns in
  [`hc_cluster_hulls()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_hulls.md)
  and
  [`hc_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_umap.md).

- Fixed degenerate polygons produced by
  [`hc_cluster_hulls()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_hulls.md)
  for very small clusters.

- Fixed stale gene names after serialising and restoring objects.

- Fixed duplicated UMAP coordinates when rerunning
  [`hc_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_umap.md).

- Annotation database download failures now produce informative
  warnings.

- Fixed several joins that relied on implicit column matching.

- Fixed `hc_annotate(dbs = "GO")` failing with “must have n \>= 2
  objects to cluster”. `rrvgo` silently drops GO terms it cannot
  resolve, so a group could reach the reduction step with fewer than two
  usable terms. Such groups are now skipped with a message.

- Fixed GO simplification failing when no group could be reduced.

- Fixed infinite term scores arising from adjusted p-values that
  underflow to zero.

- [`reduce_go_terms()`](https://buenoalvezm.github.io/HPAclusteR/reference/reduce_go_terms.md)
  now reports missing input columns and absent packages instead of
  failing with an empty message.

- Fixed deprecated `ggplot2` and `tidyselect` usage.

### Performance

- PCA, consensus clustering, UMAP, cluster hulls, and cluster stability
  are substantially faster.

- Consensus clustering now reuses the nearest-neighbour graph across
  runs.

- Annotation data are loaded once and reused across annotation
  databases.

| Step | Before | After |
|----|----|----|
| [`hc_pca()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_pca.md) | 1.40 s | 0.02 s |
| [`hc_cluster_consensus()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_consensus.md) | 5.73 s | 2.75 s |
| [`hc_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_umap.md) | 6.19 s | 0.79 s |
| [`hc_cluster_hulls()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_hulls.md) | 2.46 s | 0.74 s |
| [`hc_cluster_stability()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_stability.md) | 4.19 s | 0.60 s |

### Annotation databases

- Annotation downloads now allow up to 30 minutes for large files.

- The annotation workflow handles unavailable optional Bioconductor
  packages without failing the package examples or vignette.

## HPAclusteR 1.0.0

- Initial release.
