 #' Compute Local Inverse Simpson's Index (LISI)
#' 
#' Use this function to compute LISI scores of one or more labels.
#' NOTE: Forked from https://github.com/immunogenomics/LISI
#' 
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell. 
#' @param label_colnames Which variables to compute LISI for. 
#' @param perplexity The effective number of each cell's neighbors.
#' @param nn_eps Error bound for nearest neighbor search with \code{RANN:nn2()}.
#' Default of 0.0 implies exact nearest neighbor search.
#' 
#' @return A data frame of LISI values. Each row is a cell and each
#' column is a different label variable. 
#' 
#' @importFrom RANN nn2
#' @export 
#' 
#' @examples
#' 
#' ## Example with 400 cells. 
#' library(scIntegrationMetrics)
#' library(dplyr)
#' library(tidyr)
#' library(ggplot2)
#' library(magrittr)
#' 
#' head(scIntegrationMetrics::meta_data)
#' 
#' ## Let's color cells by labels. For label 1, there are mixed and non-mixed
#' ## groups. For label 2, all cells are well mixed. 
#' scIntegrationMetrics::X %>% 
#'   cbind(scIntegrationMetrics::meta_data) %>% 
#'   sample_frac(1L, FALSE) %>% 
#'   gather(key, val, label1, label2) %>% 
#'   ggplot(aes(X1, X2, color = val)) +
#'     geom_point(shape = 21) + 
#'     facet_wrap(~key)
#'   
#' ## Now to compute and plot the LISI values for each label. 
#' lisi_res <- compute_lisi(scIntegrationMetrics::X, scIntegrationMetrics::meta_data, c('label1', 'label2'))
#' head(lisi_res)
#' 
#scIntegrationMetrics::X %>% 
#'   cbind(lisi_res) %>% 
#'   sample_frac(1L, FALSE) %>% 
#'   gather(key, lisi_value, label1, label2) %>% 
#'   ggplot(aes(X1, X2, color = lisi_value)) +
#'     geom_point(shape = 21) + 
#'     facet_wrap(~key)
#' 
compute_lisi <- function(
  X,
  meta_data,
  label_colnames,
  perplexity = 30,
  nn_eps = 0
) {
  N <- nrow(meta_data)
  
  if (perplexity >= N) {
    perplexity <- N-1
    warning(sprintf("Cannot calculate more neighbors than there are points. Setting perplexity to %i", perplexity))
  }
  
  dknn <- nn2(X, k = perplexity + 1, eps = nn_eps)
  lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
  lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message('LISI: Cannot compute LISI on missing values')      
      return(rep(NA, N))
    } else {
      ## don't count yourself in your neighborhood
      dknn$nn.idx <- dknn$nn.idx[, 2:ncol(dknn$nn.idx)]
      dknn$nn.dists <- dknn$nn.dists[, 2:ncol(dknn$nn.dists)]
      labels <- as.integer(factor(labels)) - 1
      n_batches <- length(unique(labels))
      simpson <- compute_simpson_index(
        t(dknn$nn.dists),
        t(dknn$nn.idx) - 1, 
        labels,
        n_batches,
        perplexity
      )
      return(1 / simpson)
    }
  }))
  lisi_df <- as.data.frame(lisi_df)  
  colnames(lisi_df) <- label_colnames
  row.names(lisi_df) <- row.names(meta_data)
  return(lisi_df)
}


#' Compute Silhouette coefficient
#' 
#' Use this function to compute silhouette scores of one or more labels. 
#' 
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell. 
#' @param label_colnames Which variables to compute silhouettes for. 
#' 
#' @return A data frame of silhouette values. Each row is a cell and each
#' column is a different label variable. 
#' 
#' @importFrom cluster silhouette
#' @importFrom vegan vegdist 
#' @export 
compute_silhouette <- function(X, meta_data, label_colnames ) {
  
  N <- nrow(meta_data)
  sil_df <- data.frame(matrix(NA, N, length(label_colnames)))
  sil_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message(paste("ASW: Cannot compute silhouette on missing values.","Skipping",label_colname))
      return(rep(NA, N))
    } else if(sum(table(labels)>0) < 2){
      message(paste("ASW: Cannot compute silhouette without at least 2 label levels.","Skipping",label_colname))
      return(rep(NA, N))
    }
    else {
      dists <- vegdist(X, method="euclidean")
      labels.num <- as.numeric(as.factor(labels))
      sil <- as.data.frame(silhouette(labels.num, dists))
      return(sil$sil_width)
    }
  }))
  sil_df <- as.data.frame(sil_df)
  colnames(sil_df) <- label_colnames
  row.names(sil_df) <- row.names(meta_data)
  
  return(sil_df)
  
}

#' Single-level means
#' 
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell. 
#' @param label_colnames Which variables to compute averages for. 
#' @param level Which level to consider
#' 
#' @return A vector of means
#' 
#' @export 
compute_mean_singleLevel <- function(res, meta_data, label_colnames, level) {
  
  means <- lapply(label_colnames, function(x)
    mean(res[,x][meta_data[,x] == level])
  )
  names(means) <- label_colnames
  return(unlist(means))
  
}


#' Compute Local Inverse Simpson's Index (LISI) split by groups
#' 
#' @param X A matrix with cells (rows) and features (columns).
#' @param meta_data A data frame with one row per cell. 
#' @param label_colnames Which variables to compute averages for.
#' @param normalize Normalize LISI between 0 and 1 
#' @param split_by_colname Which variable levels use to split data
#' @param metricsLabels Which cell types to evaluate for LISI calculation
#' 
#' @return A list of data frames of LISI values (one per split_by_colname level). Each row is a cell and each
#' column is a different label variable. 
#' 
#' @export compute_lisi_splitBy
compute_lisi_splitBy <- function (X, meta_data,
                                  label_colnames,
                                  split_by_colname,
                                  metricsLabels=NULL,
                                  normalize=TRUE, ...){
  
  if (is.null(metricsLabels))
    metricsLabels <- unique(na.omit(meta_data[[meta.label]]))
  
  X.list <- split.data.frame(X,meta_data[,split_by_colname])
  meta_data.list <- split.data.frame(meta_data,meta_data[,split_by_colname])
  
  #on selected categories, and excluding NAs
  X.list <- X.list[metricsLabels]
  meta_data.list <- meta_data.list[metricsLabels]
  
  names <- names(X.list)
  
  X.list.list <- lapply(seq_along(X.list), function(i) {
    x <- X.list[[i]]
    m <- meta_data.list[[i]]
    n <- names[i]
    
    x.lisi <- compute_lisi(x, m, label_colnames=label_colnames, ... )
    if(normalize){
      label_colnames_levels <- apply(m[,label_colnames,drop=F],2,
                                     function(w) length(unique(w)))
      x.lisi <- data.frame(t(t(x.lisi-1)/(label_colnames_levels-1)))
    }
    return(x.lisi)
  })
  
  names(X.list.list) <- names
  
  return(X.list.list)
}



#' Integration metrics for single-cell data
#' 
#' Computes multiple integration metrics from a Seurat object, including several flavors of Local Inverse Simpson's Index and
#' Average Silhouette Coefficient - see details below.
#' 
#' @param object A seurat object with dimensionality reduction and meta data
#' @param meta.label Which meta data column contains cluster/celltype labels 
#' @param meta.batch Which meta data column contains batch 
#' @param method.reduction Dimensionality reduction for calculation of metrics, e.g. 'pca'
#' @param ndim Number of dimensions in methods.reduction (if NULL use all dimensions)
#' @param iLISI_perplexity Number of cells in neighborhood for integration LISI metrics (iLISI)
#' @param cLISI_perplexity Number of cells in neighborhood for cluster LISI metrics (cLISI). By default,
#'     it is calculated as 2 * average number of cells per label per dataset, not counting zero-size clusters.
#' @param metrics one or more of 'iLISI', 'norm_iLISI',
#'     'CiLISI', 'CiLISI_means','norm_cLISI', 'norm_cLISI_means',
#'     'celltype_ASW', 'celltype_ASW_means';
#'     or 'NULL' to calculate them all
#' @param metricsLabels which cluster/celltype labels to consider to summarize metrics, by default use all
#' @details The following metrics have been implemented:
#' \itemize{
#'  \item{"iLISI"}{Integration LISI: defined by \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884693/}{Korsunsky et al. Nat Methods (2019)},
#'  it quantifies the effective number of datasets in a local neighborhood, thereby
#'  measuring batch mixing}
#'  \item{"norm_iLISI"}{iLISI normalized between 0 and 1}
#'  \item{"CiLISI"}{Per-cell type iLISI: iLISI is computed separately for each cell type,
#'  and normalized between 0 and 1}
#'  \item{"CiLISI_means"}{As CiLISI, but returns the mean of means per cell type instead of global mean}
#'  \item{"norm_cLISI"}{Normalized cell-type LISI: defined by \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6884693/}{Korsunsky et al. Nat Methods (2019)},
#'  it quantifies the effective number of cell types in a local neighborhood. Here we define norm_cLISI as: 1 - normalized celltype LISI,
#'  to vary between 0 (bad cell type separation) to 1 (good cell type separation)}
#'  \item{"norm_cLISI_means"}{As norm_cLISI, but using mean of means per cell type instead of global mean}
#'  \item{"celltype_ASW"}{Cell type Average Silhoette Coefficient: It quantifies distances
#'   of cells of the same type compared to the distances to cells of other types.}
#'  \item{"celltype_ASW_means"}{As celltype_ASW, but using mean of means per cell type instead of global mean}
#' }
#' @return A list of mean values for each metric 
#' 
#' @export
getIntegrationMetrics <- function(object,
           metrics=NULL,
           meta.label,
           meta.batch,
           iLISI_perplexity=30,
           cLISI_perplexity="default",
           method.reduction="pca",
           ndim=NULL,
           metricsLabels = NULL) {
    # check input parameters
    metricsAvailable <-
      c("iLISI",
        "norm_iLISI",
        "CiLISI",
        "CiLISI_means",
        "norm_cLISI",
        "norm_cLISI_means", 
        "celltype_ASW",
        "celltype_ASW_means")
    
    if(is.null(metrics)) {metrics <- metricsAvailable }
    if (!all(metrics %in% metricsAvailable)) {
      metrics_vector <- paste(c(metricsAvailable), collapse=",")
      str <- sprintf("'metrics' is unknown. Please define one or more of %s, or 'metrics=NULL' to calculate them all", metrics_vector)
      stop(str)
    }
    
    if (!is(object, "Seurat")) {
      stop("Error: 'object' must be a 'Seurat' object.")
    }
    
    integrationMetrics <- list()
    
    #Exclude cells with NA labels (most metrics cannot account for NAs)
    meta <- object@meta.data[,c(meta.label, meta.batch)]
    ncells <- nrow(meta)
    notNA.cells <- rownames(meta)[!is.na(meta[,meta.label]) & !is.na(meta[,meta.batch])]
    object <- subset(object, cells = notNA.cells)
    meta <- meta[notNA.cells,]
    
    n.rem <- ncells - length(notNA.cells)
    if (n.rem > 0) {
      message(sprintf("Found labels with NA value. Excluding %i cells from calculation of metrics", n.rem))
    }
    
    if (is.null(metricsLabels))
      metricsLabels <- unique(na.omit(object@meta.data[[meta.label]]))
    
    message(paste("Cell type labels:", paste(metricsLabels, collapse = ",")))
    
    metricsLabels_logic <-
      object@meta.data[[meta.label]] %in% metricsLabels
    batchNames <- unique(object@meta.data[[meta.batch]])
    
    message(paste("Batches:", paste(batchNames, collapse = ",")))
    
    #get embeddings
    emb <- Embeddings(object, reduction = method.reduction)
    if (is.null(emb)) {
      ndim <- ncol(emb)
    } else if (ndim > ncol(emb)) {
      ndim <- ncol(emb)
    }
    emb <- emb[,1:ndim]

    #Determine default cLISI_perplexity
    if (cLISI_perplexity == "default") {
      t <- table(meta[,meta.label], meta[,meta.batch])
      labsize_means <- apply(t, 1, function(x){mean(x[x>0])})
      cLISI_perplexity <- round(2 * mean(labsize_means))
      mess <- sprintf("Setting default cLISI_perplexity to %0.f", cLISI_perplexity)
      message(mess)
    }
    
    #Integration LISI
    if (any(c("iLISI", "norm_iLISI") %in% metrics)) {
      lisi.this <-
        compute_lisi(
          X = emb,
          meta_data = object@meta.data,
          label_colnames = meta.batch,
          perplexity = iLISI_perplexity
        )[[1]]
      
      if ("iLISI" %in% metrics) {
        integrationMetrics[["iLISI"]] <-
          mean(lisi.this[metricsLabels_logic])
      }
      
      if ("norm_iLISI" %in% metrics) {
        lisi.this.normalized <- (lisi.this - 1) / (length(batchNames) - 1)
        integrationMetrics[["norm_iLISI"]] <-
          mean(lisi.this.normalized[metricsLabels_logic])
      }
    }
    
    #Integration LISI per celltype
    if (any(c(
      "CiLISI",
      "CiLISI_means"
    ) %in% metrics)) {
      lisi_splitByCelltype <-
        compute_lisi_splitBy(
          X = emb,
          meta_data = object@meta.data,
          label_colnames = meta.batch,
          perplexity = iLISI_perplexity,
          split_by_colname = meta.label,
          metricsLabels = metricsLabels,
          normalize = T
        )
      
      if ("CiLISI" %in% metrics) {
        integrationMetrics[["CiLISI"]] <-
          mean(unlist(lisi_splitByCelltype))
      }
      if ("CiLISI_means" %in% metrics) {
        classMeans <- sapply(lisi_splitByCelltype, function(x) mean(x[,1]))
        message("CiLISI: ", paste(names(lisi_splitByCelltype), round(classMeans, 2), " "))
        integrationMetrics[["CiLISI_means"]] <-
          mean(classMeans)
      }
    }
    
    #Cluster/celltype LISI
    if (any(c("norm_cLISI", "norm_cLISI_means") %in% metrics)) {
      lisi.this <-
        compute_lisi(
          X = emb,
          meta_data = object@meta.data,
          label_colnames = meta.label,
          perplexity = cLISI_perplexity
        )[[1]]
      
      lisi.this.normalized <-
        (lisi.this - 1) / (length(metricsLabels) - 1)
      
      if ("norm_cLISI" %in% metrics) {
        integrationMetrics[["norm_cLISI"]] <-
          1-mean(lisi.this.normalized[metricsLabels_logic])
      }
      
      if ("norm_cLISI_means" %in% metrics) {
        integrationMetrics[["norm_cLISI_means"]] <-
          1-mean(tapply(lisi.this.normalized, object@meta.data[[meta.label]], mean)[metricsLabels])
      }
    }
    
    # silhouette
    if (any(c("celltype_ASW", "celltype_ASW_means") %in% metrics)) {
      sil.this <-
        compute_silhouette(
          X = emb,
          meta_data = object@meta.data,
          label_colnames = meta.label
        )[[1]]
      
      if ("celltype_ASW" %in% metrics) {
        integrationMetrics[["celltype_ASW"]] <-
          mean(sil.this[metricsLabels_logic])
      }
      if ("celltype_ASW_means" %in% metrics) {
        integrationMetrics[["celltype_ASW_means"]] <-
          mean(tapply(sil.this, object@meta.data[[meta.label]], mean)[metricsLabels])
      }
    }
    return(integrationMetrics)
  }


