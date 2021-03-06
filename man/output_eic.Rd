% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_eic.R
\name{output_eic}
\alias{output_eic}
\title{output_eic}
\usage{
output_eic(
  path = ".",
  query_sample_name,
  query_peak_name,
  polarity = c("positive", "negative"),
  threads = 6
)
}
\arguments{
\item{path}{Work directory.}

\item{query_sample_name}{sample names}

\item{query_peak_name}{Feature names}

\item{polarity}{The polarity of data, "positive"or "negative".}

\item{threads}{Number of threads.}
}
\value{
EICs.
}
\description{
Output EIC for some peaks in some samples.
}
\author{
Xiaotao Shen
\email{shenxt1990@163.com}
}
