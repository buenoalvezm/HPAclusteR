# Master annotation pipeline: download, load, and run all enrichments

`hc_annotate()` is a master function that performs functional enrichment
analysis. It downloads necessary annotation databases, formats them, and
runs enrichment analyses for Gene Ontology (GO), KEGG pathways, and
custom databases (Human Protein Atlas, Reactome, TRRUST, PanglaoDB)
based on the clustering results stored in the AnnDatR object. It also
generates visualizations such as treemaps for GO terms and bubblemaps
for all enrichment results.

## Usage

``` r
hc_annotate(
  AnnDatR,
  dbs = c("GO", "KEGG", "Others"),
  universe = NULL,
  simplify_go = TRUE,
  db_loc = "databases",
  hpa_version = NULL,
  verbose = TRUE
)
```

## Arguments

- AnnDatR:

  AnnDatR object containing clustering results

- dbs:

  Character vector of databases to use for enrichment ("GO", "KEGG",
  "Others") (default: all)

- universe:

  Character vector of background genes (default: NULL, all genes in
  clustering data)

- simplify_go:

  Logical, group semantically similar GO terms and plot Tree map

- db_loc:

  Directory to store annotation databases (default: "databases")

- hpa_version:

  Human Protein Atlas release to download. `NULL` (default) uses the
  current release; a number such as `24` pins an archived release.

- verbose:

  Logical, print progress messages (default: TRUE)

## Value

A named list with elements:

- enrichment: Data frame with combined enrichment results from all
  selected databases

- treemaps: List of ggplot objects with GO treemaps per cluster (if GO
  selected)

- bubblemap_go: ggplot object with GO enrichment bubblemap (if GO
  selected)

- bubblemap_kegg: ggplot object with KEGG enrichment bubblemap (if KEGG
  selected)

- bubblemap_others: ggplot object with custom database enrichment
  bubblemap (if Others selected)

## Details

In case some of the databases did not download correctly, please rerun
the function. If that does not help, please manually download the
databases, place them in the specified db_loc directory and unzip the
zipped files. Links:

- Human Protein Atlas (current release):
  <https://www.proteinatlas.org/download/proteinatlas.tsv.zip>

- Human Protein Atlas (archived release `n`):
  `https://vn.proteinatlas.org/download/proteinatlas.tsv.zip`

- Reactome:
  <https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt>

- TRRUST:
  <https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv>

- PanglaoDB:
  <https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz>

Or just run the annotation only using KEGG and GO by setting dbs =
c("KEGG", "GO").

This function will require time depending on the size of the dataset and
the number of clusters, especially if GO is included. Set verbose = TRUE
to monitor the progress.

For the bubblemaps, the odds ratio is capped at 50 for better
visualization. Also, only the top 150 terms (by mean odds ratio across
clusters) are shown.

## Examples

``` r
# \donttest{
# Run clustering pipeline
adata_res <- hc_pca(example_adata, components = 40)
#> 111 of 76518 values (0.15%) are missing; using na_action = 'impute'.
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Building SNN graph from the provided distance matrix.
adata_res <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
#> Running louvain clustering at resolution 8 across 20 seeds.
#> Median number of communities across seeds: 34 (range 32-36).
#> Iteration: 0 *** value: 950.475
#> Iteration: 1 *** value: 273.082
#> Iteration: 2 *** value: 101.189
#> Iteration: 3 *** value: 92.0887
#> Iteration: 4 *** value: 90.7377
#> Iteration: 5 *** value: 86.5719
#> Iteration: 6 *** value: 86.5719
#> Minimum: 86.5719

# Enrichment against KEGG. This queries the KEGG web service, so it needs
# network access and takes a while.
enrichment_results <- hc_annotate(adata_res, dbs = "KEGG")
#> 'select()' returned 1:many mapping between keys and columns
#> Warning: 0.82% of input gene IDs are fail to map...
#> No KEGG enrichment result for cluster: 22
#> No KEGG enrichment result for cluster: 4
#> No KEGG enrichment result for cluster: 8
#> No KEGG enrichment result for cluster: 25
#> No KEGG enrichment result for cluster: 9
#> No KEGG enrichment result for cluster: 30
#> No KEGG enrichment result for cluster: 13
#> No KEGG enrichment result for cluster: 21
#> No KEGG enrichment result for cluster: 18
#> No KEGG enrichment result for cluster: 6
#> No KEGG enrichment result for cluster: 20
#> No KEGG enrichment result for cluster: 23
#> No KEGG enrichment result for cluster: 26
#> No KEGG enrichment result for cluster: 11
#> Warning: qvalue::qvalue() failed, returning NA for qvalue. Error: missing values and NaN's not allowed if 'na.rm' is FALSE
#> No KEGG enrichment result for cluster: 10
#> Warning: no package '--> No gene can be mapped....' was found
#> Warning: no package '--> Expected input gene ID: ' was found
#> Warning: no package '--> return NULL...' was found
#> No KEGG enrichment result for cluster: 1
#> No KEGG enrichment result for cluster: 29
#> No KEGG enrichment result for cluster: 15
#> No KEGG enrichment result for cluster: 32
#> No KEGG enrichment result for cluster: 5
#> No KEGG enrichment result for cluster: 14
#> 'select()' returned 1:many mapping between keys and columns
#> KEGG enrichment done.
#> 'select()' returned 1:many mapping between keys and columns
head(enrichment_results$enrichment)
#> # A tibble: 6 × 10
#>   `Cluster ID` Database      `Term ID` Term          GeneRatio BgRatio `P-value`
#>   <chr>        <chr>         <chr>     <chr>         <chr>     <chr>       <dbl>
#> 1 12           KEGG pathways hsa04382  Cornified en… 7/14      40/352   2.77e- 4
#> 2 16           KEGG pathways hsa04382  Cornified en… 3/4       40/352   5.05e- 3
#> 3 16           KEGG pathways hsa04519  Cadherin sig… 2/4       24/352   2.46e- 2
#> 4 17           KEGG pathways hsa04382  Cornified en… 15/19     40/352   1.31e-12
#> 5 19           KEGG pathways hsa04024  cAMP signali… 3/5       11/352   2.21e- 4
#> 6 19           KEGG pathways hsa04014  Ras signalin… 3/5       16/352   7.34e- 4
#> # ℹ 3 more variables: `Adjusted P-value` <dbl>, `Gene IDs` <chr>,
#> #   `Gene names` <chr>
# }
```
