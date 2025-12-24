# Initialize AnnDatR object for gene clustering (genes as samples)

`hc_initialize()` takes three data frames: an expression matrix, gene
metadata, and sample metadata. The expression matrix must have genes as
rows and samples as columns, with the first column as gene IDs. The gene
metadata (`obs`) must have a column matching the gene IDs in the
expression matrix. The sample metadata (`var`) must have a column
matching the sample IDs (column names) in the expression matrix.

## Usage

``` r
hc_initialize(
  X,
  obs,
  var,
  obs_names_col = "ensembl_id",
  var_names_col = "sample_id"
)
```

## Arguments

- X:

  Expression matrix data frame (genes as rows, samples as columns, first
  column = gene IDs)

- obs:

  Gene metadata data frame (must contain a column matching gene IDs)

- var:

  Sample metadata data frame (must contain a column matching sample IDs)

- obs_names_col:

  Name of the column in gene metadata to use as gene IDs (default:
  "ensembl_id")

- var_names_col:

  Name of the column in sample metadata to use as sample IDs (default:
  "sample_id")

## Value

AnnDatR object with genes as samples, ready for clustering
