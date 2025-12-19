# Human Protein Atlas Tissue Data (transposed)

Human Protein Atlas tissue bulk-RNA expression data in AnnDatR format
transposed.

## Usage

``` r
adata_t
```

## Format

An AnnDatR object with the following structure:

- X:

  A data frame with 78 rows and 20,163 columns representing the main
  data matrix (e.g., gene expression values).

- obs:

  A data frame with 20,162 rows and 11 columns containing metadata for
  observations (e.g., genes).

- obs_names_col:

  A string specifying the column used for observation names
  (`ensembl_id`).

- var:

  A data frame with 78 rows and 7 columns containing metadata for
  variables (e.g., samples).

- var_names_col:

  A string specifying the column used for variable names (`sample_id`).

- uns:

  A list for storing unstructured metadata (currently `NULL`).

- obsm:

  A list for storing multidimensional annotations for observations
  (currently `NULL`).

- raw:

  A slot for storing raw data (currently `NULL`).

## Source

<https://www.science.org/doi/10.1126/science.1260419>

## Examples

``` r
adata_t
#> <AnnDatR>
#>   Public:
#>     X: data.frame
#>     clone: function (deep = FALSE) 
#>     filter_obs: function (column, string, negative = FALSE) 
#>     filter_var: function (column, string, negative = FALSE) 
#>     initialize: function (prefix_name = NULL, var_names = "ensembl_ids", layer = NULL, 
#>     layers: list
#>     obs: spec_tbl_df, tbl_df, tbl, data.frame
#>     obs_names: function () 
#>     obs_names_col: ensembl_id
#>     obsm: NULL
#>     raw: NULL
#>     uns: NULL
#>     var: spec_tbl_df, tbl_df, tbl, data.frame
#>     var_names: function () 
#>     var_names_col: sample_id
```
