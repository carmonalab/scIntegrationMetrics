## scIntegrationMetrics

A collection of metrics to evaluate scRNA-seq data integration quality, including LISI, per-cell type LISI, and silhouette average width (ASW).

This package is based on [immunogenomics/LISI](https://github.com/immunogenomics/LISI) by [Korsunsky et al](https://www.nature.com/articles/s41592-019-0619-0), adds new metrics and can be applied directly on Seurat objects. For a description and motivation for the new metrics please refer to [Andreatta et al. (2023) *bioRxiv preprint*](https://www.biorxiv.org/content/10.1101/2023.07.07.548105v1).


### Installation

```r
    install.packages("remotes")
    library(remotes)
    remotes::install_github("carmonalab/scIntegrationMetrics")
```

### Usage

To illustrate how to calculate these integration metrics, we will start from a collection of pancreas datasets distributed with SeuratData.

First, install the data:
```r
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)

InstallData("panc8")
data("panc8")

panc8 <- panc8 |> NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |> RunPCA(npcs=20)

head(panc8[[]])
```

```
        orig.ident nCount_RNA nFeature_RNA   tech replicate assigned_cluster celltype dataset
D101_5        D101   4615.810         1986 celseq    celseq             <NA>    gamma  celseq
D101_7        D101  29001.563         4209 celseq    celseq             <NA>   acinar  celseq
D101_10       D101   6707.857         2408 celseq    celseq             <NA>    alpha  celseq
D101_13       D101   8797.224         2964 celseq    celseq             <NA>    delta  celseq
D101_14       D101   5032.558         2264 celseq    celseq             <NA>     beta  celseq
D101_17       D101  13474.866         3982 celseq    celseq             <NA>   ductal  celseq
```


This object consists of 8 pancreas datasets ("replicate") across five technologies ("tech"); cell types were manually annotated and stored in the "celltype" metadata column.

We can calculate metrics for batch mixing (in terms of tech used) and cell type separation on this unintegrated object:
```r
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
```r
remotes::install_github("carmonalab/STACAS")
library(STACAS)
panc8.list <- SplitObject(panc8, split.by = "tech")

panc8.stacas <- Run.STACAS(panc8.list)
```

Upon integration, metrics for both batch mixing and cell type separation were improved.
```r
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


### Metrics implemented in the package


* ***iLISI***: Local Inverse Simpson's Index (LISI) for batch mixing, as defined by [Korsunsky et al. Nat Methods (2019)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884693); it quantifies the effective number of datasets in a local neighborhood, thereby measuring batch mixing.

* ***norm_iLISI***: iLISI normalized between 0 and 1.

* ***CiLISI***: iLISI is computed separately for each cell type, and normalized between 0 and 1.

* ***CiLISI_means***: As CiLISI, but returns the mean of means per cell type instead of global mean.

* ***norm_cLISI***: LISI is calculated on the cell type labels, quantifying the effective number of cell types in a local neighborhood. We calculate this quantity as `norm_cLISI = 1 - normalized cell-type LISI`, to vary between 0 (bad cell type separation) to 1 (perfect cell type separation)

* ***norm_cLISI_means***: As norm_cLISI, but using mean of means per cell type instead of global mean.

* ***celltype_ASW***: Average Silhouette Width by celltype: it quantifies distances of cells of the same type compared to the distances to cells of other types.

* ***celltype_ASW_means***: As celltype_ASW, but using mean of means per cell type instead of global mean

### References

* Andreatta, M., Herault, L., Gueguen, P., Gfeller, D., Berenstein, A. J., & Carmona, S. J. (2023). **Semi-supervised integration of single-cell transcriptomics data.** *bioRxiv, 2023-07*

* Korsunsky, I., Millard, N., Fan, J., Slowikowski, K., Zhang, F., Wei, K., ... & Raychaudhuri, S. (2019). **Fast, sensitive and accurate integration of single-cell data with Harmony.** *Nature methods, 16(12), 1289-1296*
