\name{IGG-package}
\alias{IGG-package}
\alias{IGG}
\docType{package}
\title{Inverse Gamma-Gamma}
\description{
This package contains the functions, \code{igg} and \code{igg.normalmeans}, for implementing Bayesian linear regression, sparse normal means estimation, and variable selection with the inverse gamma-gamma (IGG) prior, as introduced by Bai and Ghosh (2018) <arXiv:1710.04369>
}

\details{

The DESCRIPTION file:
\packageDESCRIPTION{Inverse Gamma-Gamma}
\packageIndices{IGG}

This package implements the IGG model for sparse Bayesian linear regression and the normal means problem. Our package performs both estimation and model selection. The \code{igg} and \code{igg.normalmeans} functions also returns the endpoints of the credible intervals (i.e. the 2.5th and 97.5th percentiles) for every single parameter of interest, so that uncertainty quantification can be assessed.

}
\author{
Ray Bai and Malay Ghosh

Maintainer: Ray Bai <raybai07@ufl.edu>
}

\references{
Bai, R. and Ghosh, M. (2018). "The Inverse Gamma-Gamma Prior for Optimal Posterior Contraction and Multiple Hypothesis Testing." Submitted, arXiv:1711.07635.
}