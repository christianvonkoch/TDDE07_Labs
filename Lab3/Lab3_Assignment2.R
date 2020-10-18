## Assignment 2: Consider the following Poisson regression model yi given Beta ~ Poisson(exp(xiT*Beta)), i=1,...,n
## where yi isthecountforthe ithobservationinthesampleand xi isthe p-dimensional vector with covariate observations
## for the ith observation. Use the data set eBayNumberOfBidderData.dat. This dataset contains observations from 1000
## eBay auctions of coins. The response variable is nBids and records the number of bids in each auction. The 
## remaining variables are features/covariates (x): 

## a) Obtain the maximum likelihood estimator of Beta in the Poisson regression model for the eBay data 
## [Hint: glm.R, don't forget that glm() adds its own intercept so don't input the covariate Const]. Which 
## covariates are significant?

# Read data
ebay = read.table("ebayNumberOfBidderData.dat", header=TRUE)
data = ebay[, -2]

# Create model
model = glm(nBids~., family="poisson", data=data)
print(model$coefficients)
summary(model)

## Conclusion: The covariates that are significant are VerifyID, Sealed, MajBlem, LogBook, MinBidShare.

## b) Let's now do a Bayesian analysis of the Poisson regression. Let the prior be Beta~N(0,100*(XTX)^(-1)) where X
## is the n x p covariate matrix. This is a commonly used prior which is called Zellner's g-prior. Assume first that
## the posterior density is approximately multivariate normal: Beta given y ~ N(Beta~, Jy(Beta~)^(-1)) where Beta~
## is the posterior mode and Jy(Beta~) is the negative Hessian at the posterior mode. Beta~ and J can be obtained
## by numerical optimization (optim.R) exactly like you already did for the logistic regression in Lab2. 

library(mvtnorm)
# Defining constants
X = as.matrix(ebay[,2:ncol(ebay)])
Y = ebay[,1]
nFeatures = dim(X)[2]
covNames=names(ebay[,2:ncol(ebay)])

# Constructing prior
mu_prior = rep(0,nFeatures)
sigma_prior = 100*solve(t(X)%*%X) 

# Defining function for returning the log posterior
logPostPoisson = function(beta, Y, X, mu, sigma) {
  n=length(Y)
  XBeta=beta%*%t(X)
  # Defining loglikelihood
  logLike <- sum(-log(factorial(Y))+XBeta*Y-exp(XBeta))
  # Defining prior
  prior=dmvnorm(beta, mean=mu, sigma=sigma, log=TRUE)
  # Adding loglikelihood and logprior together. Since it is log both of them are added instead of multiplied
  return(logLike + prior)
}

# Defining initial values to be passed on to the optimizer
set.seed(12345)
initVals = rnorm(dim(X)[2])

# Finding the optimized betavector
optimResult = optim(initVals, logPostPoisson, Y=Y, X=X, mu=mu_prior, sigma=sigma_prior, method=c("BFGS"),
                    control=list(fnscale=-1), hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = -solve(optimResult$hessian)
names(postMode) = covNames
approx_PostStd = sqrt(diag(postCov))
names(approx_PostStd) = covNames
print("The posterior mode is:")
print(postMode)
print("The approximated standard deviations are:")
print(approx_PostStd)

## Conclusion: Through optimization we have obtained the optimal betavector as well as the hessian evaluated at the posterior
## mode. 

## c) Now, let's simulate from the actual posterior of beta using the Metropolis algorithm and compare with the
## approximate results in b). Program a general function that uses the Metropolis algorithm to generate random
## draws from an arbitrary posterior density. In order to show that it is a general function for any model, I will
## denote the vector of model parameters by theta. Let the proposal density be the multivariate normal density
## mentioned in Lecture 8 (random walk Metropolis): Theta_p given Theta(i-1) ~ N(Theta(i-1), c*Cov) where 
## Cov = Jy(Beta~)^(-1) obtained in b). The value c is a tuning parameter and should be an input to your Metropolis
## function. The user of your Metropolis function should be able to supply her own posterior density function, not
## necessarily for the Poisson regression, and still be able to use your Metropolis function. This is not so
## straightforward, unless you have come across function objects in R and the triple dot (...) wildcard argument.
## I have posted a note (HowToCodeRWM.pdf) on the course web page that describes how to do this in R. Now, use your
## new Metropolis function to sample from the posterior of beta in the Poisson regression for the eBay dataset. 
## Assess MCMC convergence by graphical methods. 

# Defining function for sampling through metropolis hastings
RVMSampler = function(previousVal, postCov, c, myFunction, ...) {
  proposalVal=rmvnorm(1, mean=previousVal, sigma=c*postCov)
  alpha=min(1, exp(myFunction(proposalVal,...)-myFunction(previousVal, ...)))
  u=runif(1)
  if(u < alpha) {
    return(proposalVal)
  } else {
    return(previousVal)
  }
}

nDraws=5000
beta_matrix = matrix(0, nDraws, ncol(X))
# Setting initial values of beta to same initVals as in the optimizer (taken randomly from normal distrib)
beta_matrix[1,]=initVals
c=0.5
set.seed(12345)

# Performing metropolis hastings
for(i in 1:nDraws) {
  if(i<nDraws) {
    beta_matrix[i+1,]=RVMSampler(beta_matrix[i,], postCov, c, logPostPoisson, Y, X, mu_prior, sigma_prior)
  }
}

iter=seq(1,nDraws,1)
par(mfrow=c(3,3))
for (i in 1:9) {
  plot(iter, beta_matrix[,i], type="l", main=paste("Convergence plot for covariate", covNames[i]),
       ylab=covNames[i])
}
par(mfrow=c(1,1), new=FALSE)

# Calculating distinct rows and dividing by total rows to get average acceptance probability
avg_alpha=dim(beta_matrix[!duplicated(beta_matrix),])[1]/dim(beta_matrix)[1]

## Conclusion: As seen in the convergence plots the covariates oscillate around the same value which was found in the previous
## problem where the optimal beta values were found through optimization. Since the variable c should be chosen 
## in a way to acquire an average acceptance rate of approximately 25-30%, the average acceptance rate were 
## calculated to approximately 33 % which is deemed to be sufficiently satisfying. 

## d) Use the MCMC draws from c) to simulate from the predictive distribution of the number of bidders in a new
## auction with the characteristics below. Plot the predictive distribution. What is the probability of no bidders
## in this new auction? Use vector x=c(1,1,1,1,0,0,0,1,0.5)

obs_X=c(1,1,1,1,0,0,0,1,0.5)
# Removing first 1000 rows since they are before the start of the convergence
approx_post_beta=beta_matrix[1001:nrow(beta_matrix),]
mean_vector=exp(approx_post_beta%*%obs_X)
set.seed(12345)
pred_distrib_bidder=rpois(10000, mean_vector)
barplot(table(pred_distrib_bidder),
        main="Histogram of the predictive distribution of no. of bidders", xlab="No. of bidders")
# Calculating the probability of no bidders with the given characteristics
prob_noBidders=sum(pred_distrib_bidder==0)/length(pred_distrib_bidder)
print(prob_noBidders)

## Conclusion: As seen in the predictive distribution the majority of cases given the specified characteristics, will result in
## either 0 or 1 bidder with the probability decreasing for additional bidders. The calculated probability for
## no bidder is 0.3581.

