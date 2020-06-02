# Reading the data from file
library(MASS)
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)

library(mvtnorm)

# Defining a function that simulates from the scaled inverse Chi-square distribution
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

BayesLinReg <- function(y, X, mu_0, Omega_0, v_0, sigma2_0, nIter){
  # Direct sampling from a Gaussian linear regression with conjugate prior:
  #
  # beta | sigma2 ~ N(mu_0, sigma2*inv(Omega_0))
  # sigma2 ~ Inv-Chi2(v_0,sigma2_0)
  # 
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   y - n-by-1 vector with response data observations
  #   X - n-by-nCovs matrix with covariates, first column should be ones if you want an intercept.
  #   mu_0 - prior mean for beta
  #   Omega_0  - prior precision matrix for beta
  #   v_0      - degrees of freedom in the prior for sigma2
  #   sigma2_0 - location ("best guess") in the prior for sigma2
  #   nIter - Number of samples from the posterior (iterations)
  #
  # OUTPUTS:
  #   results$betaSample     - Posterior sample of beta.     nIter-by-nCovs matrix
  #   results$sigma2Sample   - Posterior sample of sigma2.   nIter-by-1 vector
  
  # Compute posterior hyperparameters
  n = length(y) # Number of observations
  nCovs = dim(X)[2] # Number of covariates
  XX = t(X)%*%X
  betaHat <- solve(XX,t(X)%*%y)
  Omega_n = XX + Omega_0
  mu_n = solve(Omega_n,XX%*%betaHat+Omega_0%*%mu_0)
  v_n = v_0 + n
  sigma2_n = as.numeric((v_0*sigma2_0 + ( t(y)%*%y + t(mu_0)%*%Omega_0%*%mu_0 - t(mu_n)%*%Omega_n%*%mu_n))/v_n)
  invOmega_n = solve(Omega_n)
  
  # The actual sampling
  sigma2Sample = rep(NA, nIter)
  betaSample = matrix(NA, nIter, nCovs)
  for (i in 1:nIter){
    
    # Simulate from p(sigma2 | y, X)
    sigma2 = rScaledInvChi2(n=1, df = v_n, scale = sigma2_n)
    sigma2Sample[i] = sigma2
    
    # Simulate from p(beta | sigma2, y, X)
    beta_ = rmvnorm(n=1, mean = mu_n, sigma = sigma2*invOmega_n)
    betaSample[i,] = beta_
    
  }
  return(results = list(sigma2Sample = sigma2Sample, betaSample=betaSample))
}

mu_0=rep(0, ncol(X))
omega_0=0.01*diag(ncol(X))
v_0=1
sigma2_0=36
nIter=5000

bayes_lin_results=BayesLinReg(y, X, mu_0, omega_0, v_0, sigma2_0, nIter)
# Under quadratic loss, posterior mean is point estimate
beta_estimates=rep(0,ncol(X))
beta_credIntervals=matrix(0, ncol(X), 2)
for (i in 1:ncol(X)) {
  beta_estimates[i]=mean(bayes_lin_results$betaSample[,i])
  beta_credIntervals[i,]=quantile(bayes_lin_results$betaSample[,i], c(0.025, 0.975))
}
sigma_estimate=mean(bayes_lin_results$sigma2Sample)
sigma_credInterval=quantile(bayes_lin_results$sigma2Sample, c(0.025, 0.975))
rownames(beta_credIntervals)=covNames
beta_credIntervals[which(rownames(beta_credIntervals)=="rm"), ]

## Interpretation: for one unit increase of rooms the hosing prices will rise between 3991,475 and 5009,826 dollars
## with 95 % posterior probability. 

## b) Owner of house 381 is considering selling their house. Bought house for 10400

old_obs=as.vector(X[381,])
new_obs=old_obs
new_obs[2]=10
pred_draw=rep(0,nIter)
for (i in 1:nIter) {
  pred_draw[i]=bayes_lin_results$betaSample[i,]%*%new_obs+rnorm(1, mean=0,
                                                                sd=sqrt(bayes_lin_results$sigma2Sample[i]))
}
pred_mean=mean(pred_draw)
hist(pred_draw, breaks=50)
quantile(posterior_prices, c(0.025, 0.975))
sum(pred_draw>=30)/nIter

## c) See paper.
