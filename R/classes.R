#' Check AnnDatR Alignment
#' Verifies that the expression matrix (X) rows and columns are
#' correctly ordered and aligned with the gene (var) and cell (obs) metadata tables.
#'
#' @param adata AnnDatR object.
#'
#' @return Logical TRUE if alignment is perfect, FALSE otherwise.
#' @keywords internal
check_annotation_alignment <- function(adata) {
    x <- adata$X
    obs <- adata$obs
    var <- adata$var

    # Extract column names (excluding first column) and associated obs labels
    x_colnames <- colnames(x)[-1]
    obs_names <- as.character(obs[[1]])

    # Extract row identifiers from first column and var
    x_varorder <- x[[1]]
    var_order <- var[[1]]

    col_match <- identical(x_colnames, obs_names)
    row_match <- identical(x_varorder, var_order)

    if (col_match && row_match) {
        message(
            "Column names and row identifiers match and are in the same order."
        )
        return(TRUE)
    }
}

#' AnnDatR Class
#'
#' A flexible, Anndata-like container for transcriptomics data.
#'
#' @importFrom R6 R6Class
#'
#' @field X Primary expression data matrix (features x samples).
#' @field var Feature/Gene metadata (data.frame).
#' @field obs Cell/Sample metadata (data.frame).
#' @field var_names_col Column name in `var` used as feature names.
#' @field layers List of alternative expression data matrices (features x samples).
#' @field obs_names_col Column name in `obs` used as sample names.
#' @field uns Unstructured data (list).
#' @field obsm Multidimensional annotations of observations (list).
#' @field raw Raw unprocessed data (AnnDatR object or NULL).
#'
#' @keywords internal
AnnDatR <- R6::R6Class(
    "AnnDatR",
    public = list(
        X = NULL,
        var = NULL,
        obs = NULL,
        var_names_col = NULL,
        layers = list(),
        obs_names_col = NULL,
        uns = list(),
        obsm = list(),
        raw = NULL,

        #' Create a new AnnDatR object
        #'
        #' @param prefix_name Optional prefix for loading data from files.
        #' @param var_names Column name in `var` to use as feature names. Defaults to "ensembl_ids".
        #' @param layer Optional layer to load as X.
        #' @param file_dir Directory containing the data files. Defaults to "./".
        #' @param X Optional expression matrix (features x samples).
        #' @param var Optional feature metadata data.frame.
        #' @param obs Optional sample metadata data.frame.
        initialize = function(
            prefix_name = NULL,
            var_names = "ensembl_ids",
            layer = NULL,
            file_dir = "./",
            X = NULL,
            var = NULL,
            obs = NULL
        ) {
            if (!is.null(X) && !is.null(var) && !is.null(obs)) {
                # Direct assignment
                self$X <- X
                self$var <- var
                self$obs <- obs
                self$var_names_col <- var_names
                # Attempt to guess obs_names_col if not provided or standard
                self$obs_names_col <- colnames(self$obs)[1]
            } else if (!is.null(prefix_name)) {
                # Load data
                data <- read_anndatr_tsv(
                    prefix_name,
                    var_names,
                    layer,
                    file_dir
                )

                if (is.null(layer)) {
                    self$X <- data[['X']]
                } else {
                    self$X <- data[[layer]]
                    self$layers[[layer]] <- data[[layer]]
                }

                self$var <- data$var
                self$obs <- data$obs
                self$var_names_col <- var_names
                self$obs_names_col <- colnames(self$obs)[1]
            } else {
                stop(
                    "Provide either 'prefix_name' to load data or 'X', 'var', 'obs' for direct assignment."
                )
            }

            # Initialize lists if null
            if (is.null(self$uns)) {
                self$uns <- list()
            }
            if (is.null(self$obsm)) {
                self$obsm <- list()
            }
            if (is.null(self$layers)) {
                self$layers <- list()
            }

            self$validate()
        },

        #' Validates the consistency and integrity of the AnnDatR object.
        #'
        #' @details Checks that core matrices (X, obs, var, layers) share correct dimensions and names.
        #' @return The object itself (invisibly) if validation passes, otherwise stops with an error.
        validate = function() {
            errors <- character(0)

            # --- 1. Check for data presence ---
            if (is.null(self$X) || is.null(self$obs) || is.null(self$var)) {
                errors <- c(errors, "X, obs, or var matrix/data frame is NULL.")
            }

            # Stop early if essential components are missing
            if (length(errors) > 0) {
                stop(
                    "AnnDatR Validation Failed (Missing Core Data):\n",
                    paste(errors, collapse = "\n")
                )
            }

            # --- 2. Check Data Types ---
            if (
                !inherits(self$obs, "data.frame") &&
                    !inherits(self$obs, "tbl_df")
            ) {
                errors <- c(errors, "obs must be a data frame or tibble.")
            }
            if (
                !inherits(self$var, "data.frame") &&
                    !inherits(self$var, "tbl_df")
            ) {
                errors <- c(errors, "var must be a data frame or tibble.")
            }

            # --- 3. Check Dimensional Alignment (The most critical check) ---

            # X should have (n_vars + 1) rows and (n_obs + 1) columns
            # (since both obs/var IDs are included as a column in X, assuming data is stored wide/long-style tibble)
            expected_n_vars <- nrow(self$var)
            expected_n_obs <- nrow(self$obs)

            # Check number of rows in X (genes/features)
            if (nrow(self$X) != expected_n_vars) {
                errors <- c(
                    errors,
                    paste0(
                        "X has ",
                        nrow(self$X),
                        " rows, but var has ",
                        expected_n_vars,
                        " rows (features). They must match."
                    )
                )
            }

            # Check number of columns in X (samples/cells + 1 ID column)
            if (ncol(self$X) != expected_n_obs + 1) {
                errors <- c(
                    errors,
                    paste0(
                        "X has ",
                        ncol(self$X),
                        " columns, but obs has ",
                        expected_n_obs,
                        " rows. Expected ",
                        expected_n_obs + 1,
                        " columns (cells + 1 feature ID col)."
                    )
                )
            }

            # --- 4. Check ID Column Presence ---
            if (!(self$obs_names_col %in% colnames(self$obs))) {
                errors <- c(
                    errors,
                    paste0(
                        "obs is missing the specified observation ID column: '",
                        self$obs_names_col,
                        "'."
                    )
                )
            }
            if (!(self$var_names_col %in% colnames(self$var))) {
                errors <- c(
                    errors,
                    paste0(
                        "var is missing the specified variable ID column: '",
                        self$var_names_col,
                        "'."
                    )
                )
            }

            # --- 5. Check Layers Alignment ---
            for (layer_name in names(self$layers)) {
                layer_dim <- dim(self$layers[[layer_name]])
                if (
                    layer_dim[1] != expected_n_vars ||
                        layer_dim[2] != expected_n_obs + 1
                ) {
                    errors <- c(
                        errors,
                        paste0(
                            "Layer '",
                            layer_name,
                            "' has incorrect dimensions: ",
                            layer_dim[1],
                            "x",
                            layer_dim[2],
                            ". Expected ",
                            expected_n_vars,
                            "x",
                            expected_n_obs + 1,
                            " (vars x obs+1)."
                        )
                    )
                }
            }

            # --- 6. Check Alignment Order (NEW STEP) ---
            if (!check_annotation_alignment(self)) {
                errors <- c(
                    errors,
                    "X/layers data alignment failed. Row/column names are mismatched or misordered with obs/var tables."
                )
            }

            # --- Final Report ---
            if (length(errors) > 0) {
                stop(
                    "AnnDatR Validation Failed:\n",
                    paste(errors, collapse = "\n")
                )
            }

            invisible(self)
        }
    ),

    active = list(
        #' @field n_obs Number of observations (samples/cells)
        n_obs = function() {
            nrow(self$obs)
        },
        #' @field n_vars Number of variables (features/genes)
        n_vars = function() {
            nrow(self$var)
        },
        #' @field shape Dimensions of the expression matrix (n_obs, n_vars)
        shape = function() {
            c(self$n_obs, self$n_vars)
        },
        #' @field obs_names Sample/cell names from obs metadata
        obs_names = function(value) {
            if (missing(value)) {
                return(self$obs[[self$obs_names_col]])
            } else {
                stop(
                    "Setting obs_names directly is not yet supported. Modify 'obs' instead."
                )
            }
        },
        #' @field var_names Feature/gene names from var metadata
        var_names = function(value) {
            if (missing(value)) {
                return(self$var[[self$var_names_col]])
            } else {
                stop(
                    "Setting var_names directly is not yet supported. Modify 'var' instead."
                )
            }
        }
    )
)
