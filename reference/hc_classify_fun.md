# Perform sample category gene classification (HPA logic, sample-agnostic)

Perform sample category gene classification (HPA logic, sample-agnostic)

## Usage

``` r
hc_classify_fun(
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
