# Map Ensembl IDs in a column to gene symbols and add as a new column

Map Ensembl IDs in a column to gene symbols and add as a new column

## Usage

``` r
map_ensembl_to_symbol(df, geneid_col = "Gene IDs", new_col = "Gene names")
```

## Arguments

- df:

  Data frame with a column of Ensembl IDs (slash-separated, default:
  "Gene IDs")

- geneid_col:

  Name of the column with Ensembl IDs (default: "Gene IDs")

- new_col:

  Name for the new gene symbol column (default: "Gene names")

## Value

Data frame with an added column of gene symbols (slash-separated)
