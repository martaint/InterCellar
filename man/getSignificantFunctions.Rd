% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_int_pair_modules.R
\name{getSignificantFunctions}
\alias{getSignificantFunctions}
\title{Calculate significant function per intpair module}
\usage{
getSignificantFunctions(
  subGenePairs_func_mat,
  gpModules_assign,
  rank.terms,
  input_maxPval
)
}
\arguments{
\item{subGenePairs_func_mat}{subset of binary mat}

\item{gpModules_assign}{assignment of intpairs to modules}

\item{rank.terms}{table of ranked functions}

\item{input_maxPval}{threshold of significance}
}
\value{
table with significant functions
}
\description{
Calculate significant function per intpair module
}
