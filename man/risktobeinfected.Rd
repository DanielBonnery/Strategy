% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Generate_data.R
\name{risktobeinfected}
\alias{risktobeinfected}
\title{Computes the risk to be infected}
\usage{
risktobeinfected(
  U,
  closedistances = NULL,
  sicks,
  new.sicks = NULL,
  .distriskhalf = 5 * 10^(-4),
  jumprisk = 10^-6,
  delta = 0.01,
  previouslyexposed = c(),
  previousrisk = NULL
)
}
\description{
Computes the risk to be infected
}
\examples{
y=rep("Sane",nrow(U));y[sample(length(y),10)]<-"sick"
jumprisk=10^-6 
.distriskhalf=10^-6
}
