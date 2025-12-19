# Clusters genes using specified method and parameters

Clusters genes using specified method and parameters

## Usage

``` r
cluster_genes(
  genes,
  neighbors,
  method = "louvain",
  resolution = 1,
  seed = seed
)
```

## Arguments

- genes:

  Gene names

- neighbors:

  Neighbor graph object

- method:

  Clustering method: "louvain" or "leiden"

- resolution:

  Resolution parameter for clustering

- seed:

  Random seed for clustering

## Value

Data frame with gene cluster assignments
