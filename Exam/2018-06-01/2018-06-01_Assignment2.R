## Linear regression model for fish with 3 covariates. 

# Reading the data from file
load(file = 'fish.RData')

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

## a) Plot marginal posterior for each param

y=as.matrix(subset(fish, select="length"))
X=as.matrix(fish[,2:ncol(fish)])
covNames=colnames(X)
mu_0=rep(0,3)
omega_0=0.01*diag(1,3)
v_0=1
sigma2_0=100^2
nIter=5000

linPost=BayesLinReg(y, X, mu_0, omega_0, v_0, sigma2_0, nIter)
betaPost=linPost$betaSample
colnames(betaPost)=covNames
sigma2Post=linPost$sigma2Sample
par(mfrow=c(2,2))
for (i in 1:ncol(betaPost)) {
  hist(betaPost[,i], xlab=paste("Beta",i,sep=""), main=paste("Marginal posterior distribution of beta", i, sep=""))
}
hist(sigma2Post, xlab=expression(sigma), main="Marginal posterior distribution of sigma2")

par(mfrow=c(1,1))

## Construct 90 % equal tail interval for beta1 and interpret it.

quantile(subset(betaPost, select="age"), probs=c(0.05, 0.95))

## It can be concluded that when the age of the fish increases with one unit the length of the fish increases 
## with approximately between 2.284 and 2.960 mm with 90 % posterior probability. 

## d) New experiment fish has been grown in water tank with water temp 30 degrees celsius. Newborn fish have
## have been inserted into the tank at two time points, 30 days ago and 100 days ago. Equal amount of fish
## in the two different ages. You pick up fish randomly from water tank. Do bayesian analysis (using sim methods)
## to determine predictive distrib of the length of the picked up fish. 

x1=c(1,30,30)
x2=c(1,100,30)
x_pred=rep(0,nIter)
for (i in 1:nIter) {
  prob=runif(1)
  if(prob>0.5) {
    x_pred[i]=betaPost[i,]%*%x1+rnorm(1, mean=0, sd=sqrt(sigmaPost[i]))
  } else {
    x_pred[i]=betaPost[i,]%*%x2+rnorm(1, mean=0, sd=sqrt(sigmaPost[i]))
  }
}
hist(x_pred, main="Histogram of predictive distribution of length of fish",
     xlab="Length in mm", freq=FALSE, breaks=50)


