# Formatting functions for annotation databases (excluding GO)

Each function takes `(db, hpa, cluster_genes)` and returns a tibble with
the standardised columns `ensg_id`, `term` and `term_id`. `hpa` is the
already loaded Human Protein Atlas table, shared by every formatter that
needs it.

## Usage

``` r
database_formatting_functions
```
