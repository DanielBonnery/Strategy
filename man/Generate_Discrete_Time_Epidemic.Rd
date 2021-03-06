% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_data.R
\name{Generate_Discrete_Time_Epidemic}
\alias{Generate_Discrete_Time_Epidemic}
\title{Generate epidemic}
\usage{
Generate_Discrete_Time_Epidemic(
  U,
  TT,
  .distriskhalf = 5 * 10^(-4),
  jumprisk = 10^-6,
  delta = 0.05
)
}
\arguments{
\item{U}{a data.frame}

\item{TT}{an integer}

\item{.distriskhalf}{a positive number(default 5*10^(-4))}

\item{jumprisk}{=10^-6 a positive number}

\item{delta}{=0.05 a positive number}
}
\description{
Generate epidemic
}
\examples{
.distriskhalf=5*10^(-4);jumprisk=10^-6;delta=0.05; TT=10
UE<-Generate_Discrete_Time_Epidemic(U,3)
}
