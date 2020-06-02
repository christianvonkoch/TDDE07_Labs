## a) Use BayesLinReg to sim 5000 draws from posterior distrib of all coeff coefficients. Summarize posterior
## with point estimate under quadratic loss function and 95 % equal tail intervals. Interpret cred intervals for
## regression coefficient on nitrogen oxides concentration.

###############################
########## Problem 1 ########## 
############################### 

# Reading the data from file
library(MASS)
BostonHousing = Boston
y = BostonHousing$medv
X = cbind(1,BostonHousing[,1:13]) # Adding a column of ones for the intercept
names(X)[1] <- "intercept"
covNames <- names(X)
y <- as.numeric(y)
X <- as.matrix(X)
XNewHouse <- c(1,0.03,40,1.5,0,0.5,6,30,5,3,300,17,390,4)

if(length((grep("mvtnorm",installed.packages()[,1])))==0)
  install.packages("mvtnorm")
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
omega_0=1/10^2*diag(ncol(X))
v_0=1
sigma2_0=5^2
nIter=5000
linPost=BayesLinReg(y, X, mu_0, omega_0, v_0, sigma2_0, nIter)
betaPost=linPost$betaSample
sigma2Post=linPost$sigma2Sample
results=matrix(0,ncol(X)+1,3)
results_names=covNames
results_names=append(results_names, "sigma2")
rownames(results)=results_names
colnames(results)=c("Point estimator", "2,5%", "97,5%")
for (i in 1:ncol(X)) {
  results[i,1]=mean(betaPost[,i])
  results[i,-1]=quantile(betaPost[,i], probs=c(0.025, 0.975))
}
results[(ncol(X)+1),1]=mean(sigma2Post)
results[(ncol(X)+1),-1]=quantile(sigma2Post, probs=c(0.025, 0.975))
results

## b) Kernel density estimates. Compute posterior mode and HPD 90 % for sigma2

sigma2_kernel=density(sigma2Post)
sigma2_kernel.df=data.frame(sigma2=sigma2_kernel$x, density=sigma2_kernel$y)
sigma2_kernel.df=sigma2_kernel.df[order(-sigma2_kernel.df[,2]),]
index=dim(sigma2_kernel.df)[1]
sigma2_kernel.df$density=cumsum(sigma2_kernel.df$density)/sum(sigma2_kernel.df$density)
sigma2Cred=sigma2_kernel.df[sigma2_kernel.df$density<0.9,]
credInterval=c(min(sigma2Cred$sigma2), max(sigma2Cred$sigma2))
sigma2Mode=sigma2_kernel.df[1,]$sigma2

plot(sigma2_kernel, type="l", lwd=2, main="Kernel density estimate of sigma2", xlab=expression(sigma^2))
abline(v=sigma2Mode, col="red", lwd=1, lty=2)
abline(v=credInterval[1], col="grey", lwd=1, lty=3)
abline(v=credInterval[2], col="grey", lwd=1, lty=3)
legend("topright", legend=c("Kernel density estimate", "Posterior mode", "90 % HPD Interval"), lty=c(1,2,3),
       lwd=c(2,1,1), col=c("black", "red", "grey"))

## c) Construction company planning to build a new house with covariates given in XNewHouse. Cost is 20000 dollars
## and the company is planning to sell the house when finished. Do Bayesian analysis to determine how probable
## it is that the company will make money (that the house will sell for more than 20000 dollars).

XNewHouse <- c(1,0.03,40,1.5,0,0.5,6,30,5,3,300,17,390,4)
profitVec=rep(0,nIter)
for (i in 1:nIter) {
  profitVec[i]=-20+betaPost[i]%*%XNewHouse+rnorm(1, mean=0, sd=sqrt(sigma2Post[i]))
}
hist(profitVec)
probProfit=sum(profitVec>0)/nIter
print(probProfit)
quantile(profitVec, probs=c(0.025, 0.975))

## Very probable that the company will make a profit since 98.82 % of the posterior draws are above zero. Negative
## values are also not present in the 95 % equal tail interval which also indicates that the company will make
## a profit. 
