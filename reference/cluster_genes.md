# Clusters genes using specified method and parameters

Community detection on the shared nearest neighbour graph. Both
algorithms optimise modularity at the requested resolution and both
depend on the random number generator, so different seeds give different
partitions – which is what makes the consensus step meaningful.

## Usage

``` r
cluster_genes(
  genes,
  graph,
  method = c("louvain", "leiden"),
  resolution = 1,
  seed = 42
)
```

## Arguments

- genes:

  Gene names, in the row order of the graph.

- graph:

  An `igraph` graph built from the SNN adjacency matrix.

- method:

  Clustering method: `"louvain"` or `"leiden"`.

- resolution:

  Resolution parameter for clustering.

- seed:

  Random seed for clustering.

## Value

Data frame with gene cluster assignments, numbered from zero.
