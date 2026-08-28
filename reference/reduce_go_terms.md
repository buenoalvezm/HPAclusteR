# Reduce GO terms

Reduce GO terms

## Usage

``` r
reduce_go_terms(go_enrichment, threshold = 0.7, verbose = TRUE)
```

## Arguments

- go_enrichment:

  Data frame with GO enrichment results (must have Cluster ID, Term ID,
  Database, Adjusted P-value)

- threshold:

  Similarity threshold for reduction (default: 0.7)

- verbose:

  Logical indicating whether to report groups that could not be reduced
  (default is TRUE).

## Value

List with two elements:

- combined: Data frame with original and simplified GO terms

- reducedTerms: Data frame with reduced GO terms information
