\name{igg.normalmeans}
\alias{igg.normalmeans}

\title{Normal Means Estimation and Classification with the IGG Prior}

\description{This function provides a fully Bayesian approach for obtaining a sparse estimate of \eqn{\theta = (\theta_1, ..., \theta_n)} in the normal means problem,

\deqn{X_i = \theta_i + \epsilon_i,} 

where \eqn{\epsilon_i \sim N(0, \sigma^2)}. This is achieved by placing the inverse gamma-gamma (IGG) prior on the individual \eqn{\theta_i}'s. Variable selection can also be performed by using either thresholding the shrinkage factor in the posterior mean, or by examining the marginal 95 percent credible intervals.
}

\usage{
igg.normalmeans(x, a=NA, b=NA, sigma2=NA, var.select = c("threshold", "intervals"), 
max.steps=10000, burnin=5000)}

\arguments{
  \item{x}{an \eqn{n \times 1} multivariate normal vector.}

  \item{a}{The parameter for \eqn{IG(a,1)}. If not specified (\code{a=NA}), defaults to \eqn{1/2+1/n}. User may specify a value for \code{a} between 0 and 1.}

  \item{b}{The parameter for \eqn{G(b,1)}. If not specified (\code{b=NA}), defaults to \eqn{1/n}. User may specify a value for \code{b} between 0 and 1.}
  
  \item{sigma2}{The variance parameter. If the user does not specify this (\code{sigma2=NA}), the Gibbs sampler will estimate this using Jeffreys prior. If \eqn{sigma^2} is known or estimated separately (e.g. through empirical Bayes), the user may also specify it. }
  
  \item{var.select}{The method of variable selection. \code{threshold} selects variables by thresholding the shrinkage factor in the posterior mean. \code{intervals} will classify entries of \eqn{x} as either signals (\eqn{x_i \neq 0}) or as noise (\eqn{x_i = 0}) by examining the 95 percent marginal credible intervals. }

\item{max.steps}{The total number of iterations to run in the Gibbs sampler. Defaults to 10,000.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to 5,000.}

}

\details{
The function performs sparse estimation of \eqn{\theta = (\theta_1, ..., \theta_n)} in normal means problem. The full model is:

\deqn{X | \theta \sim N_n( \theta, \sigma^2 I_n),}
\deqn{\theta_i | (\lambda_i,\xi_i, \sigma^2) \sim N(0, \sigma^2 \lambda_i \xi_i), i = 1, ..., n,}
\deqn{\lambda_i \sim IG(a,1), i = 1, ..., n,}
\deqn{\xi_i \sim G(b,1), i = 1, ..., n,}
\deqn{\sigma^2 \propto 1/\sigma^2.}

If \eqn{\sigma^2} is known or estimated separately, the Gibbs sampler does not sample from the full conditional for \eqn{\sigma^2}.

As described in Bai and Ghosh (2018), the posterior mean is of the form \eqn{[E(1-\kappa_i | X_i)] X_i, i = 1, ..., n.} The "threshold" method for variable selection is to classify \eqn{\theta_i} as signal (\eqn{\theta_i \neq 0}) if 

\deqn{E(1 - \kappa_i | X_i) > 1/2,}

and to classify \eqn{\theta_i} as noise (\eqn{\theta_i = 0}) if

\deqn{E(1 - \kappa_i |X_i) \leq 1/2.}

}

\value{The function returns a list containing the following components:
\item{theta.hat}{ The posterior mean estimate of \eqn{\theta}. }
\item{theta.med}{ The posterior median estimate of \eqn{\theta}. }
\item{theta.intervals}{ The lower and upper endpoints of the 95 percent credible intervals for all \eqn{n} components of \eqn{\theta}. }
\item{igg.classifications}{ An \eqn{n}-dimensional binary vector with "1" if the covariate is selected and "0" if it is deemed irrelevant.}
}

\references{
Bai, R. and Ghosh, M. (2018). "The Inverse Gamma-Gamma Prior for Optimal Posterior Contraction and Multiple Hypothesis Testing." Submitted, arXiv:1711.07635.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{
###################################################
###################################################
## Example on synthetic data.                    ## 
## 5 percent of entries in a sparse vector theta ##
## are set equal to signal value A =7.           ##
###################################################
###################################################

n <- 100
sparsity.level <- 5
A <- 7

# Initialize theta vector of all zeros
theta.true <- rep(0,n)

# Set (sparsity.level)% of them to be A
q <- floor(n*(sparsity.level/100))


# Pick random indices of theta.true to equal A
signal.indices <- sample(1:n, size=q, replace=FALSE)

#######################
# Generate true theta #
#######################
theta.true[signal.indices] <- A

####################################################### 
# Generate data X by corrupting theta.true with noise #
#######################################################
X <- theta.true + rnorm(n,0,1)

##########################
# Run the IGG model on X #
##########################
# For optimal performance, should set max.steps=10,000 and burnin=5000.

igg.model <- igg.normalmeans(X, var.select="threshold", max.steps=2000, burnin=1000)

################################ 
# Calculate mean squared error #
################################
igg.mse <- sum((igg.model$theta.med-theta.true)^2)/n
igg.mse

# To compute misclassification probability and false discovery rate
true.classifications <- rep(0,n)
signal.indicies <- which(theta.true != 0)
true.classifications[signal.indices] <- 1 
igg.classifications <- igg.model$igg.classifications

false.pos <- length(which(igg.classifications != 0 & true.classifications == 0))
num.pos <- length(which(igg.classifications != 0))
false.neg <- length(which(igg.classifications == 0 & true.classifications != 0))

########################################################
# Calculate FDR and misclassification probability (MP) #
########################################################
igg.fdr <- false.pos/num.pos
igg.fdr

igg.mp <- (false.pos+false.neg)/n
igg.mp

\dontrun{
#
#
######################################################
######################################################
## Prostate cancer data analysis.                   ##
## Running this code will allow you to reproduce    ##
## the results in Section 6 of Bai and Ghosh (2018) ##
######################################################
######################################################

# Load the data
data(singh2002)
attach(singh2002)

# Only look at the gene expression data.
# First 50 rows are the cancer patients,
# and the last 52 rows are the control subjects.d

prostate.data <- singh2002$x

# Perform 2-sample t-test and obtain z-scores
n <- ncol(prostate.data)
test.stat <- rep(NA,n)
z.scores <- rep(NA, n)

# Fill in the vectors
for(i in 1:n){
  test.stat[i] <- t.test(prostate.data[51:102,i], 
                        prostate.data[1:50,i])$statistic
  z.scores[i] <- qnorm(pt(test.stat[i],100))
}


#######################################
# Apply IGG model on the z-scores.    #
# Here sigma2 is known with sigma2= 1 #
#######################################

igg.model <- igg.normalmeans(z.scores, sigma2=1, var.select="threshold")

##########################################
# How many genes flagged as significant? #
##########################################
num.sig <- sum(igg.model$igg.classifications != 0)
num.sig

####################################################### 
# Estimated effect size for 10 most significant genes #
#######################################################
most.sig <- c(610,1720,332,364,914,3940,4546,1068,579,4331)
igg.model$theta.hat[most.sig]
}

}
