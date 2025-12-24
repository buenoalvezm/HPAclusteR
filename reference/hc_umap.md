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
#> ENSG00000002745  0.39367986  0.9791806
#> ENSG00000004660 -0.39131860  0.9349804
#> ENSG00000006047  0.35270161 -0.8999294
#> ENSG00000006059  0.01820977  0.6333006
#> ENSG00000006453 -0.06702215  0.7121772
#> ENSG00000006740 -0.41821617  0.9093925
```
