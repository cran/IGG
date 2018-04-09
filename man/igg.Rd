\name{igg}
\alias{igg}

\title{Inverse Gamma-Gamma Regression}

\description{This function provides a fully Bayesian approach for obtaining a sparse estimate of the \eqn{p \times 1} vector, \eqn{\beta} in the univariate linear regression model, 

\deqn{y = X \beta + \epsilon,} 

where \eqn{\epsilon \sim N_n (0, \sigma^2 I_n)}. This is achieved by placing the inverse gamma-gamma (IGG) prior on the coefficients of \eqn{\beta}. In the case where \eqn{p > n}, we utilize an efficient sampler from Bhattacharya et al. (2016) to reduce the computational cost of sampling from the full conditional density of \eqn{\beta} to \eqn{O(n^2 p)}.
}

\usage{
igg(X, y, a=NA, b=NA, sigma2=NA, max.steps=10000, burnin=5000)}

\arguments{
  \item{X}{\eqn{n \times p} design matrix.}

  \item{y}{\eqn{n \times 1} response vector.}

  \item{a}{The parameter for \eqn{IG(a,1)}. If not specified (\code{a=NA}), defaults to \eqn{1/2+1/p}. User may specify a value for \code{a} between 0 and 1.}

  \item{b}{The parameter for \eqn{G(b,1)}. If not specified (\code{b=NA}), defaults to \eqn{1/p}. User may specify a value for \code{b} between 0 and 1.}
  
  \item{sigma2}{The variance parameter. If the user does not specify this (\code{sigma2=NA}), the Gibbs sampler will estimate this using Jeffreys prior. If \eqn{\sigma^2} is known or estimated separately (e.g. through empirical Bayes), the user may also specify it. }

\item{max.steps}{The total number of iterations to run in the Gibbs sampler. Defaults to 10,000.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to 5,000.}

}

\details{
The function performs sparse estimation of \eqn{\beta} in the standard linear regression model and variable selection from the \eqn{p} covariates. Variable selection is performed by assessing the 95 percent marginal posterior credible intervals. The lower and upper endpoints of the 95 percent credible intervals for each of the \eqn{p} covariates are also returned so that the user may assess uncertainty quantification. The full model is:

\deqn{Y | (X, \beta ) \sim N_n(X \beta, \sigma^2 I_n),}
\deqn{\beta_i | ( \lambda_i, \xi_i, \sigma^2) \sim N(0, \sigma^2 \lambda_i \xi_i), i = 1, ..., p,}
\deqn{\lambda_i \sim IG(a,1), i = 1, ..., p,}
\deqn{\xi_i \sim G(b,1), i = 1, ..., p,}
\deqn{\sigma^2 \propto 1/\sigma^2.}

If \eqn{\sigma^2} is known or estimated separately, the Gibbs sampler does not sample from the full conditional for \eqn{\sigma^2}.
}

\value{The function returns a list containing the following components:
\item{beta.hat}{ The posterior mean estimate of \eqn{\beta}. }
\item{beta.med}{ The posterior median estimate of \eqn{\beta}. }
\item{beta.intervals}{ The lower and upper endpoints of the 95 percent credible intervals for all \eqn{p} estimates in \eqn{\beta}. }
\item{igg.classifications}{ A \eqn{p}-dimensional binary vector with "1" if the covariate is selected and "0" if it is deemed irrelevant.}
}

\references{
Bai, R. and Ghosh, M. (2018). "The Inverse Gamma-Gamma Prior for Optimal Posterior Contraction and Multiple Hypothesis Testing." Submitted, arXiv:1711.07635.

Bhattacharya, A., Chakraborty, A., and Mallick, B.K. (2016). "Fast Sampling with Gaussian Scale Mixture Priors in High-Dimensional Regression." \emph{Biometrika}, \bold{69}(2): 447-457.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

###################### 
# Load diabetes data #
######################
data(diabetes)
attach(diabetes)
X <- scale(diabetes$x)
y <- scale(diabetes$y)

################################
# Fit the IGG regression model #
################################
igg.model <- igg(X=X, y=y, max.steps=5000, burnin=2500)

##############################
# Posterior median estimates #
##############################
igg.model$beta.med

###########################################
# 95 percent posterior credible intervals #
###########################################
igg.model$beta.intervals

######################
# Variable selection #
######################
igg.model$igg.classifications

}
