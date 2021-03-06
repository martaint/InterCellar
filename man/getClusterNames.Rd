% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_cluster_verse.R
\name{getClusterNames}
\alias{getClusterNames}
\title{Get clusters names from initial input data}
\usage{
getClusterNames(input.data)
}
\arguments{
\item{input.data}{preprocessed input data}
}
\value{
named list of clusters
}
\description{
Get clusters names from initial input data
}
\examples{
data(input.data)
cluster_list <- getClusterNames(input.data)
}
