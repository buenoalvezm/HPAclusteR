# Compute Shared Nearest Neighbors (SNN) Graph

`hc_snn()` computes a Shared Nearest Neighbors (SNN) graph from the
distance matrix stored in the AnnDatR object.

## Usage

``` r
hc_snn(AnnDatR, neighbors = 20, prune = 1/15, verbose = TRUE)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with distance matrix results.

- neighbors:

  Number of nearest neighbors to consider (default is 20).

- prune:

  Pruning threshold for the SNN graph (default is 1/15). Edges with a
  Jaccard weight below this value are removed.

- verbose:

  Logical indicating whether to print progress messages (default is
  `TRUE`).

## Value

An AnnDatR object with the nearest neighbour and shared nearest
neighbour graphs stored in `uns$neighbors`.

## Details

Neighbours are read directly from the distance matrix computed by
[`hc_distance()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_distance.md),
so the metric is the one chosen there. Earlier versions accepted a
`similarity` argument that was passed to `Seurat` as an Annoy metric; it
never had any effect on this code path, because Annoy is only used when
neighbours are searched in a coordinate space rather than taken from a
precomputed distance matrix. The argument has been removed.

## Examples

``` r
# Compute SNN graph on an AnnDatR object after distance calculation
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN graph from the provided distance matrix.
adata_res$uns$neighbors$snn[1:5, 1:5]
#> 5 x 5 sparse Matrix of class "dgCMatrix"
#>                 ENSG00000002745 ENSG00000004660 ENSG00000006047 ENSG00000006059
#> ENSG00000002745               1               .               .               .
#> ENSG00000004660               .               1               .               .
#> ENSG00000006047               .               .               1               .
#> ENSG00000006059               .               .               .               1
#> ENSG00000006453               .               .               .               .
#>                 ENSG00000006453
#> ENSG00000002745               .
#> ENSG00000004660               .
#> ENSG00000006047               .
#> ENSG00000006059               .
#> ENSG00000006453               1
```
