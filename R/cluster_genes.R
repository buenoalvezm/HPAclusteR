#' Cluster genes
#'
#' @param distance The distance matrix.
#' @param k The number of clusters.
#' @param method The clustering method.
#' @param genes The genes to cluster.
#' @param id An identifier for the clustering run.
#' @param r
#'
#' @return A list containing the clustering results, time taken, and an identifier.
#' @export
#'
#' @examples
#' # cluster_genes(distance = distance_matrix, k = 10, method = "kmeans", genes = c("gene1", "gene2"))

cluster_genes <-
  function(distance,
           k = 10,
           method,
           genes,
           id = NULL,
           r = 10) {

    start_time <- Sys.time()

    if (method %in% c("kmeans", "pam", "clara", "fanny", "hclust", "agnes", "diana")) {
      res <- eclust(distance,
                    FUNcluster = m,
                    k = k,
                    nboot = 500)
    }

    if (method == "fastkmeans") {
      centroids <- KMeans_arma(
        distance %>% as.matrix(),
        clusters = k,
        n_iter = 10,
        seed_mode = "random_subset"
      )
      res <- predict_KMeans(distance %>% as.matrix(), centroids)
      class(res) <- "numeric"
    }

    if (method == "medoids") {
      res <-  fastkmed(distdata = distance,
                       ncluster = k,
                       iterate = 100)
      res <- as.numeric(res$cluster)
    }


    if (method == "louvain" | method == "leiden") {
      if (method == "louvain") {
        alg = 1
      }

      if (method == "leiden") {
        alg = 4
      }

      distance <- distance %>%
        as.matrix() %>%
        set_colnames(genes) %>%
        set_rownames(genes) %>%
        as.dist()

      louv <-
        CreateSeuratObject(assay = "Exp",
                           counts = t(data.frame(row.names = genes, c(
                             1:length(genes)
                           ))))
      neighbors <-
        FindNeighbors(
          dist,
          k.param = 20,
          return.neighbor = FALSE,
          compute.SNN = TRUE,
          prune.SNN = 1 / 15,
          nn.method = "annoy",
          #Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
          n.trees = 50,
          annoy.metric = "euclidean",
          ####mirar
          nn.eps = 0,
          verbose = TRUE,
          force.recalc = FALSE,
          l2.norm = FALSE,
          cache.index = FALSE
        )

      louv@graphs$Exp_snn <- neighbors$snn

      louv <-
        FindClusters(
          louv,
          graph.name = "Exp_snn",
          resolution = r,
          algorithm = alg
        )

      col <- paste("Exp_snn_res.", r, sep = "")
      res <- as.numeric(as.character(louv@meta.data[, col]))
      k <-
        tail(sort(as.numeric(as.character(
          louv@meta.data[, col]
        ))), n = 1)
    }

    if (method == "dbscan") {
      res_dbscan <-
        distance %>%
        dbscan::dbscan(eps = 100, minPts = 5)
      res <- as.numeric(res_dbscan$cluster)

    }

    if (method == "SOM") {
      res_som <-
        som(distance %>% as.matrix(), grid = somgrid(10, k / 10, "hexagonal"))
      res <- res_som$unit.classif
    }

    if (method == "fuzzyc") {
      res_fuzzy <- FKM(distance %>% as.matrix(), k = k)
      res <- as.numeric(res_fuzzy$clus[, 1])
    }


    if (method == "fasthclust") {
      res <- flashClust(distance, method = "ward", members = NULL)
      res <-
        res %>% cutree(k) %>% dplyr::as_tibble(rownames = "gene") %>%
        column_to_rownames("gene")
    }

    end_time <- Sys.time()
    total_time <- end_time - start_time
    print(total_time)

    res <- data.frame(gene = genes, value = res)


    return(list(
      cluster = res,
      time = total_time,
      id = paste(id, method, as.character(k))
    ))
  }
