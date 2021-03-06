% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated_functions.R
\name{getWorklist}
\alias{getWorklist}
\title{getWorklist}
\usage{
getWorklist(
  table.name = "batch.xlsx",
  instrument = c("Thermo", "Agilent", "AB"),
  each.mode.number = 32,
  randommethod = c("no", "position", "injection"),
  samplenumber = NULL,
  QCstep = 8,
  conditionQCnumber = 8,
  qc.index.from = 1,
  dir = "D:\\\\Liang\\\\data\\\\PS4U\\\\HILIC\\\\batch3\\\\",
  method.path = "D:\\\\Liang\\\\Method\\\\urine\\\\HILIC\\\\",
  ms1.method.pos = "ZIC-HILIC_MS_pos",
  ms1.method.neg = "ZIC-HILIC_MS_neg",
  ms2.method.pos = c("ZIC-HILIC_MSMS_pos_NCE25", "ZIC-HILIC_MSMS_pos_NCE25",
    "ZIC-HILIC_MSMS_pos_NCE25", "ZIC-HILIC_MSMS_pos_NCE25"),
  ms2.method.neg = c("ZIC-HILIC_MSMS_neg_NCE25", "ZIC-HILIC_MSMS_neg_NCE25",
    "ZIC-HILIC_MSMS_neg_NCE25", "ZIC-HILIC_MSMS_neg_NCE25"),
  path = "."
)
}
\arguments{
\item{table.name}{Sample name.}

\item{instrument}{instrument.}

\item{each.mode.number}{each.mode.number.}

\item{randommethod}{See xcms.}

\item{samplenumber}{See xcms.}

\item{QCstep}{See xcms.}

\item{conditionQCnumber}{See xcms.}

\item{qc.index.from}{Number of threads.}

\item{dir}{dir.}

\item{method.path}{method.path.}

\item{ms1.method.pos}{ms1.method.pos.}

\item{ms1.method.neg}{ms1.method.neg.}

\item{ms2.method.pos}{ms2.method.pos.}

\item{ms2.method.neg}{ms2.method.neg.}

\item{path}{The working directory.}
}
\value{
Peak table.
}
\description{
Generate sample and worklist.
\lifecycle{deprecated}
}
\author{
Xiaotao Shen
\email{shenxt1990@163.com}
}
