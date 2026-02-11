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
  db_loc = "databases",
  hpa_version = 24,
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

- db_loc:

  Directory to store annotation databases (default: "databases")

- hpa_version:

  Version of the Human Protein Atlas to download (default: 24)

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

- Human Protein Atlas version 24:
  https://v%d.proteinatlas.org/download/proteinatlas.tsv.zip

- Reactome:
  https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt

- Trrust: https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv

- Panglao:
  https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz

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
# Run clustering pipeline
adata_res <- hc_pca(example_adata, components = 40)
adata_res <- hc_distance(adata_res, components = 20)
adata_res <- hc_snn(adata_res, neighbors = 15)
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by ‘spam’
#> Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
#> Also defined by ‘spam’
#> Building SNN based on a provided distance matrix
#> Computing SNN
adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#> Iteration: 0 *** value: 948.279
#> Iteration: 1 *** value: 70.5063
#> Iteration: 2 *** value: 22.246
#> Iteration: 3 *** value: 22.2447
#> Iteration: 4 *** value: 22.2447
#> Minimum: 22.2447
#> Joining with `by = join_by(cons_cluster)`

# Run annotation pipeline
enrichment_results <- hc_annotate(adata_res, dbs = "KEGG")
#> 
#> 
#> 'select()' returned 1:many mapping between keys and columns
#> Warning: 0.82% of input gene IDs are fail to map...
#> No KEGG enrichment result for cluster: 13
#> No KEGG enrichment result for cluster: 22
#> No KEGG enrichment result for cluster: 12
#> No KEGG enrichment result for cluster: 19
#> No KEGG enrichment result for cluster: 26
#> No KEGG enrichment result for cluster: 9
#> No KEGG enrichment result for cluster: 31
#> No KEGG enrichment result for cluster: 23
#> No KEGG enrichment result for cluster: 24
#> No KEGG enrichment result for cluster: 8
#> No KEGG enrichment result for cluster: 28
#> No KEGG enrichment result for cluster: 10
#> No KEGG enrichment result for cluster: 1
#> No KEGG enrichment result for cluster: 20
#> No KEGG enrichment result for cluster: 6
#> No KEGG enrichment result for cluster: 25
#> No KEGG enrichment result for cluster: 4
#> No KEGG enrichment result for cluster: 11
#> No KEGG enrichment result for cluster: 18
#> No KEGG enrichment result for cluster: 30
#> 'select()' returned 1:many mapping between keys and columns
#> KEGG enrichment done.
#> 'select()' returned 1:many mapping between keys and columns
head(enrichment_results$enrichment)
#> # A tibble: 6 × 10
#>   `Cluster ID` Database      `Term ID` Term          GeneRatio BgRatio `P-value`
#>   <chr>        <chr>         <chr>     <chr>         <chr>     <chr>       <dbl>
#> 1 14           KEGG pathways hsa03040  Spliceosome   8/27      13/350    4.30e-7
#> 2 15           KEGG pathways hsa01100  Metabolic pa… 6/11      67/350    8.29e-3
#> 3 16           KEGG pathways hsa05020  Prion disease 5/9       20/350    3.96e-5
#> 4 16           KEGG pathways hsa04540  Gap junction  4/9       10/350    4.01e-5
#> 5 16           KEGG pathways hsa05130  Pathogenic E… 4/9       13/350    1.32e-4
#> 6 16           KEGG pathways hsa04814  Motor protei… 4/9       14/350    1.83e-4
#> # ℹ 3 more variables: `Adjusted P-value` <dbl>, `Gene IDs` <chr>,
#> #   `Gene names` <chr>
```
