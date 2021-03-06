% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_upload.R
\name{read.CPDBv2}
\alias{read.CPDBv2}
\title{Read output from CellPhoneDB v2.}
\usage{
read.CPDBv2(folder)
}
\arguments{
\item{folder}{folder containing output}
}
\value{
input.data which is the pre-processed object with annotated L-R pairs
}
\description{
Output is a folder containing 4 .txt files 
- deconvoluted.txt: containing list of single genes and their mean 
expression in each cluster (not considered);
- means.txt: containing list of interacting pairs with info regarding L/R, 
annotation strategy and mean value of all pairs over cluster couples.
- pvalues.txt: same as means, but containing pvalue of each pair, for each 
cluster couple.
- significant_means.txt: only means of those pairs that have pvalue < 0.05. 
Has one more column:rank. 
If the statistical analysis is not run, the folder would contain only 
deconvoluted and means
}
