% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_function_verse.R
\name{annotateGO}
\alias{annotateGO}
\title{Perform GO annotation of input data}
\usage{
annotateGO(
  input_select_ensembl,
  input_go_evidence_exclude,
  input_go_sources_checkbox,
  input.data
)
}
\arguments{
\item{input_select_ensembl}{ensembl version selected by user}

\item{input_go_evidence_exclude}{evidence codes to exclude by user}

\item{input_go_sources_checkbox}{GO sources to use by user}

\item{input.data}{preprocessed input data}
}
\value{
GO_annotation
}
\description{
Perform GO annotation of input data
}
