###########################################
# Function to implement IGG Gibbs Sampler #
###########################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
# a = parameteter for IG(a,1). If not specified, defaults to (1/2+1/n).
#     User should specify a value between 0 and 1
# b = parameter for G(b,1). If not specified, defaults to 1/n.
#     User should specify a value between 0 and 1
# sigma2 = variance term
#      User may specify this. If no value is given, then we use 
#      Jeffrey's prior to update sigma2
# var.select = method for selecting variables
#      Can select "threshold" for thresholding rule. Defaults to this 
#      Can also select based on marginal credible intervals
# max.steps = # of iterations to run MCMC
# burnin = # of samples in burn-in

##################
# RETURNS A LIST #
##################
# theta.hat = posterior mean estimate for theta
# theta.med = posterior median estimate for theta
# lower.endpoints = 2.5th percentile of posterior density
# upper.endpoints = 97.5th percentile of posterior density
# theta.classifications = binary vector of classifications, according to 
#                         classification method chosen by user


igg.normalmeans = function(x, a=NA, b=NA, sigma2=NA, 
                           var.select = c("threshold","intervals"), 
                           max.steps = 10000, burnin=5000) {
  
  if (burnin > max.steps)
    stop("Burnin cannot be greater than # of iterations.")
  if (length(x) == 0)
    stop("Please enter a vector length greater than 0.")
  
  # Number of samples in noisy vector
  n <- length(x)

	# Lists to hold draws of theta_i's and kappa_i's
	theta.samples <- rep(list(rep(NA,n)), max.steps)
	kappa.samples <- rep(list(rep(NA,n)), max.steps)

	# if a = NA, set it equal to 1/2+1/n
	if(is.na(a))
	  a <- 1/2+1/n
	# if user specified a, must be between 0 and 1.
	if ( (a<0) | (a>1) )
	  stop("ERROR: a should be between 0 and 1. \n")
	
	# if b = NA, set it equal to 1/n
	if(is.na(b))
	  b <- 1/n
	# if user specified b, must be between 0 and 1
	if ( (b<0) | (b>1))
	  stop("ERROR: b should be between 0 and 1. \n")
	
	
	############################
	# Initial guess for sigma2 #
	############################
	# If sigma2 was specified by the user, then ignore this.
	# Otherwise initialize using the variance of the residuals
	# of the LASSO estimator.
	if(is.na(sigma2)){
	  sigma2 <- var(x)
	  sigma2.known <- FALSE
	} else { 
	  sigma2.known <- TRUE
	}
	
	###########################
	# Initial guess for theta #
	###########################
	theta <- rep(mean(x), n)
	
	###################################
	# Initial guesses for lambda, tau #
	###################################
	lambda <- rep(1, n)
	xi <- rep(1, n)
	kappa <- 1/(1+lambda*xi)

	# Start the Gibbs sampler
	j <- 0
	while (j < max.steps) {
		j <- j + 1

		###########
		# Counter #
		###########
		if (j %% 100 == 0) 
		   cat("Iteration:", j, "\n")
		
		####################
		# Sample theta_i's #
		####################
	  theta.means <- (1-kappa)*x   # conditional means for thetas
	  theta.sds <- sqrt(sigma2*(1-kappa))   # conditional variances for thetas
	  theta <- rnorm(n, mean=theta.means, sd=theta.sds)  # sample thetas
	  
	  # Store theta samples
	  theta.samples[[j]] <- theta
	  
	  ################################
	  # Sample lambda_i's and xi_i's #
	  ################################
	  
	  # Rate and scale parameters for Inverse Gamma
	  alpha <- a + 1/2
	  beta <- theta^2/(2*sigma2*xi)+1
	 
	  # Parameters for Generalized Inverse Gaussian
    u <- b - 1/2
	  w <- 2
	  
	  for (i in 1:n){
	    lambda[i] <- rigamma(1, alpha, beta[i])
	    v <- max((theta[i]^2)/(sigma2*lambda[i]),.Machine$double.eps) # To prevent chi parameter from collapsing to zero
	    xi[i] <- rgig(1,lambda=u, chi=v, psi=w)
	  }
	  
	  # Store kappa samples
	  kappa <- 1/(1+lambda*xi)    # shrinkage factors
	  kappa.samples[[j]] <- kappa
	  
	  ############################		
	  # Sample sigma2 if unknown #
	  ############################
	  if(!sigma2.known){
	    resid <- sum((x-theta)^2)
	    scaleterm.2 <- t(theta)%*%diag(1/(lambda*xi))%*%theta
	    sigma2 <- rigamma(1, n, (resid+scaleterm.2)/2)
	  }
  }
	
	###################
	# Discard burn-in #
	################### 
	theta.samples <- tail(theta.samples,max.steps-burnin)
	kappa.samples <- tail(kappa.samples,max.steps-burnin)  

	#######################################
	# Extract the posterior mean, median, #
	# 2.5th, and 97.5th percentiles       #
	#######################################
	theta.sample.mat <- simplify2array(theta.samples)
	# Posterior mean
	theta.hat <- rowMeans(theta.sample.mat)
	# Posterior median
	theta.med <- apply(theta.sample.mat, 1, function(x) median(x))
	# endpoints of 95% posterior credible intervals
	theta.intervals <- apply(theta.sample.mat, 1, function(x) quantile(x, prob=c(.025,.975)))
	
	##############################
	# Perform variable selection #
	##############################
	# Initialize vector of binary entries: 0 for inactive variable b_i, 1 for active b_j
	igg.classifications <- rep(0,n)
	
	if(var.select=="intervals"){
	  # Find the active covariates
	  for(k in 1:n){
	    if(theta.intervals[1,k] < 0 && theta.intervals[2,k] < 0){
	      igg.classifications[k] <- 1
	    }
	    else if(theta.intervals[1,k] > 0 && theta.intervals[2,k] > 0){
	      igg.classifications[k] <- 1
	    }
	  }
	}
	  
	if(var.select=="threshold"){
	  # Estimate the shrinkage factor kappa_i's from the MCMC samples
	  kappa.sample.mat <- simplify2array(kappa.samples)
	  kappa.estimates <- rowMeans(kappa.sample.mat)

	  # Return indices of the signals according to our classification rule
	  signal.indices <- which((1-kappa.estimates)>=0.5)
	  # Reset classified signals as 1
	  igg.classifications[signal.indices] <- 1
	 }

	######################################
	# Return list of beta.est, beta.med, #
	# lower.endpoints, upper.endpoints,  #
	# and igg.classifications            #
	######################################
	# beta.hat = posterior mean point estimator
	# beta.med = posterior median point estimator
	# igg.intervals = endpoints of 95% posterior credible intervals
	# igg.classifications = selected variables
	
	igg.output <- list(theta.hat = theta.hat,
	                   theta.med = theta.med,
	                   theta.intervals = theta.intervals,
	                   igg.classifications = igg.classifications)
	# Return list
	return(igg.output)
}