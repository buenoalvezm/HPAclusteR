# Shared fixtures. The clustering pipeline is deterministic but not free, so
# each stage is built once and cached for the rest of the suite.

.fixture_cache <- new.env(parent = emptyenv())

cached <- function(name, build) {
  if (is.null(.fixture_cache[[name]])) {
    .fixture_cache[[name]] <- build()
  }
  .fixture_cache[[name]]
}

adata_pca <- function() {
  cached("pca", function() {
    suppressMessages(hc_pca(example_adata, components = 40, verbose = FALSE))
  })
}

adata_snn <- function() {
  cached("snn", function() {
    suppressMessages(hc_snn(
      hc_distance(adata_pca(), components = 20),
      neighbors = 15,
      verbose = FALSE
    ))
  })
}

adata_clustered <- function() {
  cached("clustered", function() {
    suppressMessages(hc_cluster_consensus(
      adata_snn(),
      resolution = 8,
      n_seeds = 20,
      verbose = FALSE
    ))
  })
}

adata_umap <- function() {
  cached("umap", function() {
    suppressMessages(hc_umap(adata_clustered(), verbose = FALSE))
  })
}

#' A tiny synthetic AnnDatR with three well-separated gene groups.
#'
#' Built through `hc_initialize()`, so it exercises the documented entry point
#' and matches the gene-as-observation layout the pipeline expects.
toy_adata <- function(n_genes = 30L, n_samples = 9L, seed = 1L) {
  set.seed(seed)
  gene_ids <- sprintf("ENSG%08d", seq_len(n_genes))
  sample_ids <- sprintf("S%02d", seq_len(n_samples))

  group <- rep(seq_len(3L), length.out = n_genes)
  block_size <- n_samples %/% 3L

  expression <- t(vapply(
    seq_len(n_genes),
    function(i) {
      values <- stats::runif(n_samples, 0, 1)
      block <- seq(
        (group[i] - 1L) * block_size + 1L,
        length.out = block_size
      )
      values[block] <- values[block] + 50
      values
    },
    numeric(n_samples)
  ))

  X <- tibble::as_tibble(
    stats::setNames(as.data.frame(expression), sample_ids)
  )
  X <- dplyr::bind_cols(tibble::tibble(ensembl_id = gene_ids), X)

  obs <- tibble::tibble(
    ensembl_id = gene_ids,
    gene_symbol = paste0("G", seq_len(n_genes))
  )
  var <- tibble::tibble(
    sample_id = sample_ids,
    tissue_name = rep(c("brain", "liver", "blood"), each = block_size)
  )

  hc_initialize(X, obs = obs, var = var)
}
