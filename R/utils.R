 #' Compute Local Inverse Simpson's Index (LISI)
#' 
#' Use this function to compute LISI scores of one or more labels.
#' NOTE: Function from https://github.com/immunogenomics/LISI repository 
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
#' 
#' @return A list of data frames of LISI values (one per split_by_colname level). Each row is a cell and each
#' column is a different label variable. 
#' 
#' @export compute_lisi_splitBy
compute_lisi_splitBy <- function (X, meta_data,
                                  label_colnames,
                                  split_by_colname,
                                  normalize=TRUE, ...){
  
  X.list <- split.data.frame(X,meta_data[,split_by_colname])
  meta_data.list <- split.data.frame(meta_data,meta_data[,split_by_colname])
  
  names <- names(X.list)
  
  X.list.list <- lapply(seq_along(X.list), function(i) {
    x <- X.list[[i]]
    m <- meta_data.list[[i]]
    n <- names[i]
    
    x.lisi <- compute_lisi(x, m, label_colnames=label_colnames, ... )
    if(normalize){
      label_colnames_levels <- apply(m[,label_colnames,drop=F],2,
                                     function(x) length(unique(x)))
      x.lisi <- data.frame(t(t(x.lisi-1)/(label_colnames_levels-1)))
    }
    message("LISI splitBy: Processing group ", n)
    return(x.lisi)
  })
  
  #names(X.list.list) <- unique(meta_data[,split_by_colname]) # NOTE: make sure the order of split list and unique() is always the same
  names(X.list.list) <- names
  
  return(X.list.list)
}



#' Compute multiple integration metrics from a Seurat object
#' 
#' @param object A seurat object with dimensionality reduction and meta data
#' @param meta.label Which meta data column contains cluster/celltype labels 
#' @param meta.batch Which meta data column contains batch 
#' @param method.reduction reduction method to consider, eg 'pca'
#' @param metrics one or more of 'batch_LISI', 'batch_nLISI', 'batch_nLISI_means', 'batch_nLISI_perCellType', 'batch_nLISI_perCellType_means','1-celltype_nLISI', '1-celltype_nLISI_means', 'celltype_ASW', 'celltype_ASW_means', or leave 'NULL' to calculate them all
#' @param metricsLabels which cluster/celltype labels to consider to summarize metrics, by default use all
#' @return A list of mean values for each metric 
#' 
#' @export
getIntegrationMetrics <- function(object,
           metrics=NULL,
           meta.label,
           meta.batch,
           lisi_perplexity=30,
           method.reduction="pca",
           metricsLabels = NULL) {
    # check input parameters
    metricsAvailable <-
      c(
        "batch_LISI",
        "batch_nLISI",
        "batch_nLISI_means",
        "batch_nLISI_perCellType",
        "batch_nLISI_perCellType_means",
        "1-celltype_nLISI",
        "1-celltype_nLISI_means",
        "celltype_ASW",
        "celltype_ASW_means"
      )
    
    if(is.null(metrics)) {metrics <- metricsAvailable }
    if (!all(metrics %in% metricsAvailable)) {
      stop(
        "Error: 'metrics' is unknown. Please define one or more of 'batch_LISI', 'batch_nLISI', 'batch_nLISI_means', 'batch_nLISI_perCellType', 'batch_nLISI_perCellType_means','1-celltype_nLISI', '1-celltype_nLISI_means', 'celltype_ASW', 'celltype_ASW_means', or leave 'NULL' to calculate them all"
      )
    }
    
    if (!is(object, "Seurat")) {
      stop("Error: 'object' must be a 'Seurat' object.")
    }
    
    integrationMetrics <- list()
    
    if (is.null(metricsLabels))
      metricsLabels <- unique(object@meta.data[[meta.label]])
    
    message(paste("Cell type labels:", paste(metricsLabels, collapse = ",")))
    
    metricsLabels_logic <-
      object@meta.data[[meta.label]] %in% metricsLabels
    batchNames <- unique(object@meta.data[[meta.batch]])
    
    message(paste("Batches:", paste(batchNames, collapse = ",")))
    
    
    #batch lisi
    if (any(c("batch_LISI", "batch_nLISI", "batch_nLISI_means") %in% metrics)) {
      lisi.this <-
        compute_lisi(
          object@reductions[[method.reduction]]@cell.embeddings,
          meta_data = object@meta.data,
          label_colnames = meta.batch,
          perplexity = lisi_perplexity
        )[[1]]
      
      if ("batch_LISI" %in% metrics) {
        integrationMetrics[["batch_LISI"]] <-
          mean(lisi.this[metricsLabels_logic])
        
      }
      
      if ("batch_nLISI" %in% metrics) {
        lisi.this.normalized <- (lisi.this - 1) / (length(batchNames) - 1)
        
        integrationMetrics[["batch_nLISI"]] <-
          mean(lisi.this.normalized[metricsLabels_logic])
        
        if ("batch_nLISI_means" %in% metrics) {
          integrationMetrics[["batch_nLISI_means"]] <-
            mean(tapply(lisi.this.normalized, object@meta.data[[meta.label]], mean)[metricsLabels]) #means per cell type
          
        }
      }
    }
    
    #batch lisi per celltype
    
    
    if (any(c(
      "batch_nLISI_perCellType",
      "batch_nLISI_perCellType_means"
    ) %in% metrics)) {
      lisi_splitByCelltype <-
        compute_lisi_splitBy(
          object@reductions[[method.reduction]]@cell.embeddings,
          meta_data = object@meta.data,
          label_colnames = meta.batch,
          perplexity = lisi_perplexity,
          split_by_colname = meta.label,
          normalize = T
        )
      
      if ("batch_nLISI_perCellType" %in% metrics) {
        integrationMetrics[["batch_nLISI_perCellType"]] <-
          mean(unlist(lisi_splitByCelltype)[metricsLabels_logic])
      }
      if ("batch_nLISI_perCellType_means" %in% metrics) {
        classMeans <- sapply(lisi_splitByCelltype, function(x) mean(x[,1]))[metricsLabels] # only considering the first `label_colnames`
        message("batch_nLISI_perCellType: ", paste(names(lisi_splitByCelltype), round(classMeans, 2), " "))
        integrationMetrics[["batch_nLISI_perCellType_means"]] <-
          mean(classMeans)
      }
    }
    
    #cluster/celltype lisi
    if (any(c("1-celltype_nLISI", "1-celltype_nLISI_means") %in% metrics)) {
      lisi.this <-
        compute_lisi(
          object@reductions[[method.reduction]]@cell.embeddings,
          meta_data = object@meta.data,
          label_colnames = meta.label,
          perplexity = lisi_perplexity
        )[[1]]
      
      lisi.this.normalized <-
        (lisi.this - 1) / (length(metricsLabels) - 1)
      
      if ("1-celltype_nLISI" %in% metrics) {
        integrationMetrics[["1-celltype_nLISI"]] <-
          1-mean(lisi.this.normalized[metricsLabels_logic])
      }
      
      if ("1-celltype_nLISI_means" %in% metrics) {
        integrationMetrics[["1-celltype_nLISI_means"]] <-
          1-mean(tapply(lisi.this.normalized, object@meta.data[[meta.label]], mean)[metricsLabels])
      }
    }
    
    # silhouette
    if (any(c("celltype_ASW", "celltype_ASW_means") %in% metrics)) {
      sil.this <-
        compute_silhouette(
          object@reductions[[method.reduction]]@cell.embeddings,
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


