# Format an annotation database (excluding GO)

Format an annotation database (excluding GO)

## Usage

``` r
format_annotation_database(raw_db, db_id, db_files, cluster_genes)
```

## Arguments

- raw_db:

  A tibble containing the raw database data.

- db_id:

  The ID of the database to format (e.g., "secretome_location").

- db_files:

  Named list mapping database IDs to file paths.

- cluster_genes:

  Character vector of Ensembl IDs to keep.

## Value

A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
