# Map Entrez IDs in a column back to Ensembl IDs

Map Entrez IDs in a column back to Ensembl IDs

## Usage

``` r
map_entrez_to_ensembl(df, geneid_col = "Gene IDs")
```

## Arguments

- df:

  Data frame with a column named 'Gene IDs' (slash-separated Entrez IDs)

- geneid_col:

  Name of the column with Entrez IDs (default: "Gene IDs")

## Value

Data frame with an added column of Ensembl IDs (slash-separated)
