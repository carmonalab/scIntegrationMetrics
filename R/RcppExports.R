# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Compute the Local Inverse Simpson Index (LISI)
#' 
#' @param D Distance matrix of K nearest neighbors.
#' @param knn_idx Adjacency matrix of K nearest neighbors.
#' @param batch_labels A categorical variable.
#' @param n_batches The number of categories in the categorical variable.
#' @param perplexity The effective number of neighbors around each cell.
#' @param tol Stop when the score converges to this tolerance.
#' @export
compute_simpson_index <- function(D, knn_idx, batch_labels, n_batches, perplexity = 15, tol = 1e-5) {
    .Call('_scIntegrationMetrics_compute_simpson_index', PACKAGE = 'scIntegrationMetrics', D, knn_idx, batch_labels, n_batches, perplexity, tol)
}

