## a) Use conjugate priors, standard normal and invchisq and use BayesLinReg to simulate 5000 draws from posterior
## distrib

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

mu_0=rep(0,ncol(X))
omega_0=1/100*diag(ncol(X))
v_0=1
sigma2_0=36
nIter=5000
post_distrib=BayesLinReg(y,X, mu_0, omega_0, v_0, sigma2_0, nIter)
post_beta=post_distrib$betaSample
colnames(post_beta)=covNames
lstat_post=subset(post_beta, select="lstat")
par(mfrow=c(1,1))
plot(density(lstat_post), main="Posterior density of lstat", lwd=2)
credInterval=quantile(lstat_post, probs=c(0.05, 0.95))
abline(v=credInterval[1], col="grey", lwd=3, lty=3)
abline(v=credInterval[2], col="grey", lwd=3, lty=3)

# Since posterior of beta is the student t-distrib the distrib is symmetric and therefore HPD interval is the same
# as equal tail interval

new_obs=X[9,]
names(new_obs)=covNames
new_obs_2=new_obs
new_obs_2[which(names(new_obs)=="lstat")]=new_obs_2[which(names(new_obs)=="lstat")]*0.7
post_sigma2=post_distrib$sigma2Sample
pred_price1=post_beta%*%new_obs+rnorm(nIter, mean=0, sd=sqrt(post_sigma2))
pred_price2=post_beta%*%new_obs_2+rnorm(nIter, mean=0, sd=sqrt(post_sigma2))
hist(pred_price1, breaks=50, main="Histogram of predicted price before change")
hist(pred_price2, breaks=50, main="Histogran of predicted price after change")
pred_price_house9=post_beta[,14]*(new_obs[14]*0.7-new_obs[14])
mean(pred_price_house9)
quantile(pred_price_house9, probs=c(0.025, 0.975))

# For a house like number 9 it will increase the house price with high posterior probability. 
