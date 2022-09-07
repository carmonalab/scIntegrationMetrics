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
  X, meta_data, label_colnames, perplexity = 30, nn_eps = 0
) {
  N <- nrow(meta_data)
  dknn <- nn2(X, k = perplexity * 3, eps = nn_eps)
  lisi_df <- data.frame(matrix(NA, N, length(label_colnames)))
  lisi_df <- Reduce(cbind, lapply(label_colnames, function(label_colname) {
    labels <- data.frame(meta_data)[, label_colname, drop = TRUE]
    if (any(is.na(labels))) {
      message('Cannot compute LISI on missing values')      
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
      message(paste("Cannot compute silhouette on missing values.","Skipping",label_colname))
      return(rep(NA, N))
    } else if(sum(table(labels)>0) < 2){
      message(paste("Cannot compute silhouette without at least 2 label levels.","Skipping",label_colname))
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
#' @param split_by_colname Which variable levels use to split data
#' 
#' @return A list of data frames of LISI values (one per split_by_colname level). Each row is a cell and each
#' column is a different label variable. 
#' 
#' @export 
compute_lisi_splitBy <- function (X, meta_data, label_colnames, split_by_colname, normalize=T, ...){
  
  X.list <- split.data.frame(X,meta_data[,split_by_colname])
  meta_data.list <- split.data.frame(meta_data,meta_data[,split_by_colname])
  
  X.list.list <- lapply(seq_along(X.list), function(i) {
    x <- X.list[[i]]
    m <- meta_data.list[[i]]
    x.lisi <- compute_lisi(x, m, label_colnames=label_colnames, ... )
    if(normalize){
      label_colnames_levels <- apply(m[,label_colnames,drop=F],2,function(x) length(unique(x)))
      x.lisi <- t(t(x.lisi-1)*(label_colnames_levels-1))
    }
    message("Processing group ",unique(m[,split_by_colname]))
    return(x.lisi)
  })
  names(X.list.list) <- unique(meta_data[,split_by_colname])
  return(X.list.list)
}


