# AnnDatR Class

A flexible, Anndata-like container for transcriptomics data.

## Public fields

- `X`:

  Primary expression data matrix (features x samples).

- `var`:

  Feature/Gene metadata (data.frame).

- `obs`:

  Cell/Sample metadata (data.frame).

- `var_names_col`:

  Column name in `var` used as feature names.

- `layers`:

  List of alternative expression data matrices (features x samples).

- `obs_names_col`:

  Column name in `obs` used as sample names.

- `uns`:

  Unstructured data (list).

- `obsm`:

  Multidimensional annotations of observations (list).

- `raw`:

  Raw unprocessed data (AnnDatR object or NULL).

## Active bindings

- `n_obs`:

  Number of observations (samples/cells)

- `n_vars`:

  Number of variables (features/genes)

- `shape`:

  Dimensions of the expression matrix (n_obs, n_vars)

- `obs_names`:

  Sample/cell names from obs metadata

- `var_names`:

  Feature/gene names from var metadata

## Methods

### Public methods

- [`AnnDatR$new()`](#method-AnnDatR-new)

- [`AnnDatR$validate()`](#method-AnnDatR-validate)

- [`AnnDatR$clone()`](#method-AnnDatR-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    AnnDatR$new(
      prefix_name = NULL,
      var_names = "ensembl_ids",
      layer = NULL,
      file_dir = "./",
      X = NULL,
      var = NULL,
      obs = NULL
    )

#### Arguments

- `prefix_name`:

  Optional prefix for loading data from files.

- `var_names`:

  Column name in `var` to use as feature names. Defaults to
  "ensembl_ids".

- `layer`:

  Optional layer to load as X.

- `file_dir`:

  Directory containing the data files. Defaults to "./".

- `X`:

  Optional expression matrix (features x samples).

- `var`:

  Optional feature metadata data.frame.

- `obs`:

  Optional sample metadata data.frame. Validates the consistency and
  integrity of the AnnDatR object.

------------------------------------------------------------------------

### Method `validate()`

#### Usage

    AnnDatR$validate()

#### Details

Checks that core matrices (X, obs, var, layers) share correct dimensions
and names.

#### Returns

The object itself (invisibly) if validation passes, otherwise stops with
an error.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    AnnDatR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
