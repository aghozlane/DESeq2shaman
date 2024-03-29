% Generated by roxygen2 (4.0.2): do not edit by hand
\name{estimateBetaPriorVar}
\alias{estimateBetaPriorVar}
\alias{estimateMLEForBetaPriorVar}
\title{Steps for estimating the beta prior variance}
\usage{
estimateBetaPriorVar(object, betaPriorMethod = c("weighted", "quantile"),
  upperQuantile = 0.05)

estimateMLEForBetaPriorVar(object, maxit = 100, useOptim = TRUE,
  useQR = TRUE)
}
\arguments{
\item{object}{a DESeqDataSet}

\item{betaPriorMethod}{the method for calculating the beta prior variance,
either "quanitle" or "weighted":
"quantile" matches a normal distribution using the upper quantile of the finite MLE betas.
"weighted" matches a normal distribution using the upper quantile, but weighting by the variance of the MLE betas.}

\item{upperQuantile}{the upper quantile to be used for the
"quantile" or "weighted" method of beta prior variance estimation}

\item{maxit}{as defined in \code{link{nbinomWaldTest}}}

\item{useOptim}{as defined in \code{link{nbinomWaldTest}}}

\item{useQR}{as defined in \code{link{nbinomWaldTest}}}
}
\value{
for \code{estimateMLEForBetaPriorVar}, a DESeqDataSet, with the
necessary information stored in order to calculate the prior variance.
for \code{estimateBetaPriorVar}, the vector of variances for the prior
on the betas in the \code{\link{DESeq}} GLM
}
\description{
These lower-level functions are called within \code{\link{DESeq}} or \code{\link{nbinomWaldTest}}.
End users should use those higher-level function instead.
For advanced users to use this function, first run \code{estimateMLEForBetaPriorVar}
and then run \code{estimateBetaPriorVar}.
}

