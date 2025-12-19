# Create UMAP embeddings from SNN graph

`hc_umap()` generates UMAP embeddings based on the Shared Nearest
Neighbors (SNN) graph stored in the AnnDatR object.

## Usage

``` r
hc_umap(AnnDatR, n_epochs = NULL, seed = 42, verbose = TRUE)
```

## Arguments

- AnnDatR:

  An AnnDatR object containing the data with SNN graph results.

- n_epochs:

  Number of epochs for UMAP optimization. If NULL, it will be set to 200
  for datasets with more than 10,000 cells and 500 otherwise (default is
  NULL).

- seed:

  Random seed for reproducibility (default is 42).

- verbose:

  Logical indicating whether to print progress messages (default is
  TRUE).

## Value

UMAP embeddings stored within the AnnDatR object.

## Examples

``` r
# Generate UMAP embeddings on an AnnDatR object after SNN graph computation
adata_res <- hc_pca(example_adata, components = 40)
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res <- hc_umap(adata_res)
head(adata_res$obsm$X_umap)
#>                      UMAP_1     UMAP_2
#> ENSG00000000003  0.11184286 -0.7571628
#> ENSG00000000005 -0.02018156 -0.8654167
#> ENSG00000000419 -0.28973365 -0.0678750
#> ENSG00000000457 -0.04218644  0.2172408
#> ENSG00000000460  0.09176493  0.4276099
#> ENSG00000000938  0.80212647  0.8295610
```
