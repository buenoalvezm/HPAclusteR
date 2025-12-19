# Compute Shared Nearest Neighbors (SNN) Graph

`hc_snn()` computes a Shared Nearest Neighbors (SNN) graph from the
distance matrix stored in the AnnDatR object.

## Usage

``` r
hc_snn(
  AnnDatR,
  neighbors = 20,
  prune = 1/15,
  similarity = "euclidean",
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with distance matrix results.

- neighbors:

  Number of nearest neighbors to consider (default is 20).

- prune:

  Pruning threshold for the SNN graph (default is 1/15).

- similarity:

  Distance metric to use for nearest neighbor search (default is
  "euclidean"). Other options include "cosine", "manhattan", and
  "hamming".

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

SNN graph stored within the AnnDatR object.

## Examples

``` r
# Compute SNN graph on an AnnDatR object after distance calculation
adata_res <- hc_pca(example_adata, components = 40)
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res$uns$neighbors$snn
#> A Graph object containing 1000 cells
```
