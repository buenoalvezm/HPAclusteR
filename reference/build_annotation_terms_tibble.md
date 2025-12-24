# Build annotation terms tibble from all databases in db_files

Build annotation terms tibble from all databases in db_files

## Usage

``` r
build_annotation_terms_tibble(db_files, AnnDatR)
```

## Arguments

- db_files:

  Named list mapping database IDs to file paths (as returned by
  get_annot_dbs()).

- AnnDatR:

  AnnDatR object containing clustering results.

## Value

A tibble with columns: db_id, ensg_id, term, term_id.
