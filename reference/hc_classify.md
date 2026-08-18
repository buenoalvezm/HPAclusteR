# Perform sample category gene classification (HPA logic, sample-agnostic)

`hc_classify` classifies genes based on their expression patterns across
different sample categories using a logic similar to that employed by
the Human Protein Atlas (HPA). The function processes each cluster of
genes identified in the AnnDatR object and categorizes them into
specificity categories such as "Enriched", "Group enriched", "Enhanced",
and "Low specificity". It also generates pie charts to visualize the
distribution of these categories within each cluster.

## Usage

``` r
hc_classify(
  AnnDatR,
  sample_categories,
  cutoff_detected = 1,
  lowest_exp_limit = 0.1,
  group_enriched_limit = NULL,
  fold = 4
)
```

## Arguments

- AnnDatR:

  AnnDatR object containing gene expression data and sample annotations

- sample_categories:

  Name of the column in AnnDatR\$var containing sample categories

- cutoff_detected:

  Expression value threshold to consider a gene as detected (default: 1)

- lowest_exp_limit:

  Minimum expression value to consider for calculations (default: 0.1)

- group_enriched_limit:

  Maximum number of sample categories for group enriched classification
  (default: NULL, auto-calculated)

- fold:

  Fold change threshold for classification (default: 4)

## Value

A tibble with gene classification results including ENSG, specificity
category, specific sample categories, tau value, and distribution
category

## Examples

``` r
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

# Classify genes based on sample categories
gene_classification <- hc_classify(adata_res, "tissue_name")
head(gene_classification$classification[["1"]])
#> # A tibble: 6 × 5
#>   ENSG            spec_category  spec_sample_categories   tau dist_category   
#>   <chr>           <chr>          <chr>                  <dbl> <chr>           
#> 1 ENSG00000105428 Enriched       testis                  0.99 Detected in some
#> 2 ENSG00000124449 Enriched       testis                  0.98 Detected in some
#> 3 ENSG00000153498 Group enriched testis                  0.98 Detected in some
#> 4 ENSG00000163958 Group enriched testis                  0.97 Detected in some
#> 5 ENSG00000166329 Group enriched testis                  0.99 Detected in some
#> 6 ENSG00000173699 Group enriched testis                  0.97 Detected in some
gene_classification$pie_charts[["1"]]
```
