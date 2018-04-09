###########################################
# Function to implement Gibbs sampler for #
# Regression Using IGG                    #
###########################################

######################
# FUNCTION ARGUMENTS #
######################
# X = design matrix (n*p)
# y = response vector (n*1)
# a = parameteter for IG(a,1). If not specified, defaults to (1/2+1/p).
#     User should specify a value between 0 and 1
# b = parameter for G(b,1). If not specified, defaults to 1/p.
#      User should specify a value between 0 and 1
# sigma2 = variance term
#      User may specify this. If no value is given, then we use 
#      Jeffrey's prior to update sigma2
# max.steps = # of times to run MCMC
# burnin = # of samples in burn-in

##################
# RETURNS A LIST #
##################
# beta.hat = posterior mean estimate for beta (p*1 vector)
# beta.med = posterior median estimate for beta (p*1 vector)
# beta.intervals = 95% posterior credible intervals for each coefficient
# igg.classifications = binary vector of classifications, according to 
#                       method of examining 95% credible intervals 

igg = function(X, y, a=NA, b=NA, sigma2=NA, max.steps = 10000, burnin=5000) {
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max.steps)
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  
  # Extract dimensions n and p
  n <- nrow(X)
	p <- ncol(X)
	
	# if a = NA, set it equal to 1/2+1/p
	if(is.na(a))
	  a <- 1/2+1/p	  
	# if user specified a, must be between 0 and 1.
	if ( (a<0) | (a>1) )
	  stop("ERROR: a should be between 0 and 1. \n")
	
	# if b = NA, set it equal to 1/p
	if(is.na(b))
	  b <- 1/p
	# if user specified b, must be between 0 and 1
	if ( (b<0) | (b>1))
	  stop("ERROR: b should be between 0 and 1. \n")
	
	# Time saving
	XtX <- crossprod(X, X)	
	Xty <- crossprod(X, y)
	# List to hold the draws of B
	beta.samples <- rep(list(rep(NA,p)), max.steps)
	
	############################
	# Initial guesses for beta #
	############################
	# Take the LASSO estimator with lambda chosen through CV
	lambda.min <- cv.glmnet(X,y)$lambda.min
	beta <- as.vector(glmnet(X,y,lambda=lambda.min)$beta)
	
	############################
	# Initial guess for sigma2 #
	############################
	# If sigma2 was specified by the user, then ignore this.
	# Otherwise initialize using the variance of the residuals
	# of the LASSO estimator.
	if(is.na(sigma2)){
	    resid <- y - crossprod(t(X), beta)
      sigma2 <- as.double(var(resid)) # Initial guess for Sigma
      sigma2.known <- FALSE
	} else { 
	    sigma2.known <- TRUE
  }
	#####################
	# Initial guess for #
	# lambda_i and xi_i # 
	#####################
	lambda <- rep(1,p) # Initial guesses for lambda_i's
	xi <- rep(1,p)     # Initial guesses for xi_i's
  
  ###########################
	# Start the Gibbs sampler #
  ###########################
	j <- 0
	while (j < max.steps) {
		j <- j + 1

		if (j %% 100 == 0) {
			cat("Iteration:", j, "\n")
		}
		
		###############
		# Sample beta #
		###############
		if (p <= n){
		  ridge <- XtX + diag(1/(lambda*xi))
		  inv.ridge <- chol2inv(chol(ridge))
		  # Conditional posterior mean
		  cond.M <- crossprod(t(inv.ridge), Xty) 
		  # Draw from MVN density
		  beta <- mvrnorm(1, cond.M, sigma2*inv.ridge)
		  
		} else if (p>n){
		  # Use the more efficient sampling algorithm if p>n
		  
		  # Set D, Phi, and alpha
		  Phi <- X/sqrt(sigma2)
		  D <- sigma2*diag(lambda*xi)
		  alpha <- y/sqrt(sigma2)
		  
		  # Draw u ~ MVN(0, D)
		  u <- mvrnorm(1, rep(0,p), D)
		  
		  # Draw delta ~ MVN(0, I_n)
		  delta <- mvrnorm(1, rep(0,n), diag(n))
		  
		  # Set v = Phi%*%u + delta
		  v <- crossprod(t(Phi), u) + delta
		  
		  # Solve for w
		  w <- chol2inv(chol(Phi%*%D%*%t(Phi)+diag(n)))%*%(alpha-v)
		  
		  # Draw beta from conditional distribution
		  beta <- u + D%*%t(Phi)%*%w
		}
		
		################################
		# Sample lambda_i's and xi_i's #
		################################
		
		# Shape and scale for lambda_i's
		ig.shape <- a + 1/2
		ig.scale <- beta^2/(2*sigma2*xi)+1
		
		# Parameters for xi_i's
		u <- b - 1/2
		w <- 2
		
		for (i in 1:p){
		  lambda[i] <- rigamma(1, ig.shape, ig.scale[i])
		  # To prevent chi parameter from collapsing to zero
		  v <- max((beta[i]^2)/(sigma2*lambda[i]),.Machine$double.eps) 
		  xi[i] <- rgig(1, lambda=u, chi=v, psi=w)
		}
		
		############################		
		# Sample sigma2 if unknown #
		############################
		if(!sigma2.known){
		  ig.shape <- (n+p)/2
		  resid <- sum((y-crossprod(t(X),beta))^2)
		  scaleterm.2 <- t(beta) %*% diag(1/(lambda*xi))%*%beta
		  sigma2 <- rigamma(1, ig.shape, (resid+scaleterm.2)/2)
		}
		
		# Save the most recent estimate of beta to the list
		beta.samples[[j]] <- beta
	}
	
  ###################
  # Discard burn-in #
  ###################
  beta.samples <- tail(beta.samples,max.steps-burnin)
  
  #######################################
  # Extract the posterior mean, median, #
  # 2.5th, and 97.5th percentiles       #
  #######################################
  beta.sample.mat <- simplify2array(beta.samples)
  # Posterior mean
	beta.hat <- rowMeans(beta.sample.mat)
	# Posterior median
	beta.med <- apply(beta.sample.mat, 1, function(x) median(x))
	# endpoints of 95% posterior credible intervals
	beta.intervals <- apply(beta.sample.mat, 1, function(x) quantile(x, prob=c(.025,.975)))

  ##############################
  # Perform variable selection #
  ##############################
	# Initialize vector of binary entries: 0 for inactive variable b_i, 1 for active b_j
	igg.classifications <- rep(0,p)
	
	# Find the active covariates
	for(k in 1:p){
	  if(beta.intervals[1,k] < 0 && beta.intervals[2,k] < 0)
	      igg.classifications[k] <- 1
	   else if(beta.intervals[1,k] > 0 && beta.intervals[2,k] > 0)
	      igg.classifications[k] <- 1
	 }
	
  ######################################
  # Return list of beta.est, beta.med, #
	# lower.endpoints, upper.endpoints,  #
	# and igg.classifications            #
  ######################################
	# beta.hat = posterior mean point estimator
	# beta.med = posterior median point estimator
	# beta.intervals = endpoints of 95% posterior credible intervals
	# igg.classifications = variable selection
	
  igg.output <- list(beta.hat = beta.hat,
                     beta.med = beta.med,
                     beta.intervals = beta.intervals,
                     igg.classifications = igg.classifications)
  # Return list
	return(igg.output)
}
