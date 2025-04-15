# scIntegrationMetrics

A collection of metrics to evaluate scRNA-seq data integration quality, including LISI, per-cell type LISI (CiLISI), and silhouette average width (ASW). This package provides convenient utilities to assess how well datasets are integrated, in terms of batch mixing and cell type separation.

This package builds upon [immunogenomics/LISI](https://github.com/immunogenomics/LISI) and adds new features:

-   Per-cell type metrics (e.g., CiLISI, ASW by cell type)

-   Normalization options for easier comparison across datasets

-   Seamless compatibility with Seurat objects

-   Functions to compute metrics in subsets (e.g., per cell type or region)

For a detailed description and motivation for the new metrics, see the sections below and refer to [Andreatta et al. (2024)](https://doi.org/10.1038/s41467-024-45240-z).

### Installation

``` r
install.packages("remotes")
remotes::install_github("carmonalab/scIntegrationMetrics")
```

### Usage

To illustrate how to calculate these integration metrics, we will start from a collection of pancreas datasets sequenced with different technologies and distributed with [SeuratData](https://github.com/satijalab/seurat-data).

First, install the data:

``` r
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)

InstallData("panc8")
data("panc8")

# Enforce Seurat object version compatibility, if needed
panc8 <- UpdateSeuratObject(panc8)

panc8 <- panc8 |> NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |> RunPCA(npcs=20)
```

We can calculate metrics for batch mixing and cell type separation on this unintegrated object, by specifying the metadata columns that contain batch and cell type information ("tech" and "celltype" respectively):

``` r
library(scIntegrationMetrics)

metrics <- getIntegrationMetrics(panc8, meta.label = "celltype",
                                 meta.batch = "tech",
                                 iLISI_perplexity = 20)

unlist(metrics)
```

```         
iLISI              1.08611688
norm_iLISI         0.02152922
CiLISI             0.02364108
CiLISI_means       0.08793179
norm_cLISI         0.96961856
norm_cLISI_means   0.87557344
celltype_ASW       0.15838012
celltype_ASW_means 0.21731279
```

We can apply an integration method such as [STACAS](https://github.com/carmonalab/STACAS) to mitigate batch effects:

``` r
remotes::install_github("carmonalab/STACAS")
library(STACAS)
panc8.list <- SplitObject(panc8, split.by = "tech")

panc8.stacas <- Run.STACAS(panc8.list)
```

We can then calculate metrics after integration, and verify whether batch mixing and cell type separation were improved upon integration:

``` r
metrics.integrated <- getIntegrationMetrics(panc8.stacas, meta.label = "celltype",
                                 meta.batch = "tech",
                                 iLISI_perplexity = 20)
unlist(metrics.integrated)
```

```         
iLISI              2.0815327
norm_iLISI         0.2703832
CiLISI             0.2705296
CiLISI_means       0.2886549
norm_cLISI         0.9861910
norm_cLISI_means   0.8945726
celltype_ASW       0.2872263
celltype_ASW_means 0.3309478
```

### :bar_chart: Metrics implemented in the package

-   ***iLISI***: Local Inverse Simpson's Index (LISI) for batch mixing, as defined by [Korsunsky et al. Nat Methods (2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884693); it quantifies the effective number of datasets in a local neighborhood, thereby measuring batch mixing. Higher values indicate better batch mixing.

-   ***norm_iLISI***: iLISI normalized between 0 and 1, for comparibility.

-   ***CiLISI***: Cell-type aware version of iLISI. iLISI computed separately for each cell type or cluster, normalized between 0 and 1, and averaged across all cells (global mean). By default, CiLISI is calculated only for groups with at least 10 cells and 2 distinct batch labels (configurable).

-   ***CiLISI_means***: As CiLISI, but returns mean of per-group CiLISI values (i.e., average of the means per group). instead of a global average.

-   ***norm_cLISI***: LISI is calculated on the cell type labels, quantifying the effective number of cell types in a local neighborhood. We calculate this quantity as `norm_cLISI = 1 - normalized cell-type LISI`, to vary between 0 (bad cell type separation) to 1 (perfect cell type separation)

-   ***norm_cLISI_means***: As norm_cLISI, but using mean of means per cell type instead of global mean.

-   ***celltype_ASW***: Average Silhouette Width by celltype: it quantifies distances of cells of the same type compared to the distances to cells of other types.

-   ***celltype_ASW_means***: As celltype_ASW, but using mean of means per cell type instead of global mean


### :newspaper_roll: NEWS

A recent independent benchmarking study [Rautenstrauch & Ohler. bioRxiv (2025)](https://doi.org/10.1101/2025.01.21.634098) demonstrates that **CiLISI** as defined in this package, is among the top-performing metrics for evaluating batch effect removal — particularly in the presence of nested batch effects. Unlike silhouette-based metrics, which often fail in such scenarios, **CiLISI** shows robust and discriminative performance across both simulated and real-world datasets.


### References

-   Andreatta, M., Herault, L., Gueguen, P., Gfeller, D., Berenstein, A. J., & Carmona, S. J. (2024). **Semi-supervised integration of single-cell transcriptomics data.** *Nature Communications, 15(1), 1-13*

-   Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... & Raychaudhuri, S. (2019). **Fast, sensitive and accurate integration of single-cell data with Harmony.** *Nature methods, 16(12), 1289-1296*
