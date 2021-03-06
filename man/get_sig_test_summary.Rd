% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_test_summary-class.R
\name{get_sig_test_summary}
\alias{get_sig_test_summary}
\title{get_sig_test_summary}
\usage{
get_sig_test_summary(fit, testable_features = row.names(fit), contrast,
  p_threshold = 0.05)
}
\arguments{
\item{fit}{A DGEGLM from a edgeR model fit.}

\item{testable_features}{Which rows of the DGEGLM should be
considered?}

\item{contrast}{A matrix or vector defining the contrasts over
which the LRT test should summarise.}

\item{p_threshold}{A single float given the FDR threshold under
which a given feature is considered 'significant'.}
}
\description{
get_sig_test_summary
}
