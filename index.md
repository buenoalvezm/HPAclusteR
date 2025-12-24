# HPAclusteR

`HPAclusteR` is an R package developed by the **Human Protein Atlas** to
perform **gene clustering** from transcriptomics data. The package
provides a complete pipeline for clustering genes, visualizing
clustering results, and annotating gene clusters with functional
databases such as **Gene Ontology** and others. It is designed to
facilitate the analysis of transcriptomics data and help researchers
identify biologically meaningful gene clusters.

> Note: This package is built to work with `AnnDatR` objects, which are
> designed to handle transcriptomics data efficiently. The package that
> conceptualizes this data structure can be found
> [here](https://github.com/emiliosk/AnnDatR).

## Features

### Gene Clustering Pipeline:

- Perform PCA for dimensionality reduction.
- Calculate distances between genes.
- Construct Shared Nearest Neighbor (SNN) graphs.
- Perform consensus clustering to identify robust gene clusters.
- Generate UMAP embeddings for visualization.
- Create UMAP cluster hulls to highlight cluster boundaries.

### Annotate Gene Clusters with Functional Databases:

- Gene Ontology (GO)
- KEGG Pathways
- Humnan Protein Atlas (HPA)
  - Tissue Expression
  - Single Cell Type Expression
  - Blood Expression
  - Brain Expression
  - Subcellular Location
  - Secretome Location
  - Protein Class
- Reactome Pathways
- PanglaoDB Cell Markers
- Trrusted Cell Markers

### Ready-to-use Visualization Functions:

- Plot UMAP embeddings with clusters and/or hulls.
- Plot gene expression per cluster across samples.
- Plot annotation results with treemaps and buble heatmaps.
- Plot cluster comparison results with network graphs.

## Installation

You can install the development version of `HPAclusteR` from GitHub:

``` r
# Install devtools if not already installed
install.packages("devtools")

# Install HPAclusteR
devtools::install_github("buenoalvezm/HPAclusteR")
```

## Usage

Start right away with the
[`hc_auto_cluster()`](https://buenoalvezm.github.io/HPAclusteR/reference/hc_auto_cluster.md)
function, which performs the complete gene clustering pipeline on an
`AnnDatR` object.

``` r
library(HPAclusteR)

# Example input: AnnDatR object with transcriptomics data
adata_res <- hc_auto_cluster(example_adata)
```

## Issues and Support

If you encounter any bugs or you want to recommend new features and
changes to existing ones, please open a [new
issue](https://github.com/buenoalvezm/HPAclusteR/issues) on our GitHub
repository.

## Contact

For any questions or further information, please contact us at
<k.antono@outlook.com>.
