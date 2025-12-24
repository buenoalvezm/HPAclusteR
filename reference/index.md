# Package index

## Utilities

General utility functions.

- [`hc_initialize()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_initialize.md)
  : Initialize AnnDatR object for gene clustering (genes as samples)

## Workflows

Functions for complete, one-call workflows.

- [`hc_auto_cluster()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_auto_cluster.md)
  : One-call clustering pipeline

## Clustering Pipeline

Functions for gene clustering.

- [`hc_pca()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_pca.md)
  : Perform PCA on AnnDatR object
- [`hc_kaisers_rule()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_kaisers_rule.md)
  : Determine number of components using Kaiser's rule
- [`hc_distance()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_distance.md)
  : Compute distance matrix from PCA scores
- [`hc_snn()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_snn.md)
  : Compute Shared Nearest Neighbors (SNN) Graph
- [`hc_cluster_consensus()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_consensus.md)
  : Consensus clustering of genes
- [`hc_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_umap.md)
  : Create UMAP embeddings from SNN graph
- [`hc_cluster_hulls()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_hulls.md)
  : Calculate UMAP cluster hulls

## Clustering Evaluation

Functions for evaluating clustering results.

- [`hc_cluster_stability()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_stability.md)
  : Calculate ARI scores for multiple clustering seeds
- [`hc_cluster_compare()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_cluster_compare.md)
  : Compare clusters using hypergeometric test

## Cluster Annotation

Functions for annotating clusters with functional databases.

- [`hc_annotate()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_annotate.md)
  : Master annotation pipeline: download, load, and run all enrichments
- [`hc_classify()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_classify.md)
  : Perform sample category gene classification (HPA logic,
  sample-agnostic)

## Visualization

Functions for visualizing clustering results.

- [`hc_plot_umap()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_plot_umap.md)
  : Plot UMAP with clusters and hulls
- [`hc_plot_expression()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_plot_expression.md)
  : Plot per-cluster heatmaps for scaled and z-score expression
- [`theme_hc()`](https://buenoalvezm.github.io/HPAclusteR/reference/theme_hc.md)
  : A simple theme for visualizations

## Built in datasets

Datasets included with the package for demonstration and testing
purposes.

- [`example_adata`](https://buenoalvezm.github.io/HPAclusteR/reference/example_adata.md)
  : Human Protein Atlas Tissue Data (transposed and subsampled)
