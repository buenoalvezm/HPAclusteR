# Human Protein Atlas Tissue Data (transposed and subsampled)

Human Protein Atlas tissue bulk-RNA expression data in AnnDatR format
transposed and subsampled (981 selected genes).

## Usage

``` r
example_adata
```

## Format

An AnnDatR object with the following structure:

- X:

  A data frame with 78 rows and 982 columns representing the main data
  matrix (e.g., gene expression values).

- obs:

  A data frame with 981 rows and 11 columns containing metadata for
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
example_adata
#> AnnDatR object with n_obs x n_vars = 981 x 78
#>     obs: ensembl_id, gene_symbol, chromosome_name, gene_biotype, start_position, end_position, length
#>     var: sample_id, tissue_name, consensus_tissue, organ_name, color, colors_esk, colors_esj
```
