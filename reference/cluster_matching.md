# Cluster matching algorithm

`cluster_matching()` implements an algorithm to match clusters from two
datasets based on significant overlaps, ensuring that there are no
many-to-many matches, while still possible to have one-to-many or
many-to-one matches.

## Usage

``` r
cluster_matching(cluster_comparison, total_genes)
```

## Arguments

- cluster_comparison:

  A tibble containing pairwise cluster comparison results, including
  columns: cluster_A, cluster_B, n_overlap, adj_p_val.

- total_genes:

  Total number of genes considered in the comparison.

## Value

A list containing:

- matched_clusters: A tibble of matched clusters adhering to the
  matching rules.

- match_score: A numeric score representing the total normalized overlap
  percentage.
