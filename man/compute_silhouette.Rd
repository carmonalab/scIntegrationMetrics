% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_silhouette}
\alias{compute_silhouette}
\title{Compute Silhouette coefficient}
\usage{
compute_silhouette(X, meta_data, label_colnames)
}
\arguments{
\item{X}{A matrix with cells (rows) and features (columns).}

\item{meta_data}{A data frame with one row per cell.}

\item{label_colnames}{Which variables to compute silhouettes for.}
}
\value{
A data frame of silhouette values. Each row is a cell and each
column is a different label variable.
}
\description{
Use this function to compute silhouette scores of one or more labels.
}
