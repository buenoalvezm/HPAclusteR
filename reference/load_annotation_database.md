# Load and format an annotation database

Load and format an annotation database

## Usage

``` r
load_annotation_database(db_id, db_files, AnnDatR)
```

## Arguments

- db_id:

  The ID of the database to load (e.g., "reactome").

- db_files:

  Named list mapping database IDs to file paths.

- AnnDatR:

  AnnDatR object containing clustering results.

## Value

A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
