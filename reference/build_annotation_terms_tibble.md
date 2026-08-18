# Build annotation terms tibble from all databases in db_files

Build annotation terms tibble from all databases in db_files

## Usage

``` r
build_annotation_terms_tibble(db_files, AnnDatR)
```

## Arguments

- db_files:

  Named list mapping database IDs to file paths (as returned by
  [`get_annot_dbs()`](https://buenoalvezm.github.io/HPAclusteR/reference/get_annot_dbs.md)).

- AnnDatR:

  AnnDatR object containing clustering results.

## Value

A tibble with columns: `db_id`, `ensg_id`, `term`, `term_id`.

## Details

Eight of the supported databases are different columns of the same Human
Protein Atlas table, and two more need it to map gene symbols to Ensembl
IDs. That file is several hundred megabytes, so it is read once here and
shared across all formatters rather than re-read for each database.
