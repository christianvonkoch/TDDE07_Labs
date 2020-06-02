###############################
########## Problem 3 ########## 
############################### 

# Reading the cars data from file
load("cars.RData")

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

## a) Linear regression problem with given dataset. Use Mattias function to derive joint posterior. 
## i) Plot marginal distributions of each param
## ii) Compute point estimates for each regression coefficient assuming loss function
## iii) Construct 95 % equal tail probability intervals for each parameter and interpret them


y=cars$mpg
x=as.matrix(cars[2:ncol(cars)])
mu_0=c(0,0,0,0)
omega_0=0.01*diag(x=4)
nu_0=1
sigma_sq_0=36
jointPostDistrib=BayesLinReg(y, x, mu_0, omega_0, nu_0, sigma_sq_0, 1000)
hist(jointPostDistrib$sigma2Sample, breaks=10, main=paste("Marginal distribution of", expression(sigma^2)),
     xlab=expression(sigma^2))
par(mfrow=c(2,2))
for(i in 1:4) {
  hist(jointPostDistrib$betaSample[,i], breaks=10, main=paste("Marginal distribution of ", expression(beta), i, 
                                                              sep=""), xlab=paste(expression(beta),i, sep=""))
}
title("Marginal distributions of the different betavalues", line=-1, outer=TRUE)
par(mfrow=c(1,1))

# Linear loss function is posterior median
median(jointPostDistrib$sigma2Sample)
median(jointPostDistrib$betaSample[,1])
median(jointPostDistrib$betaSample[,2])
median(jointPostDistrib$betaSample[,3])
median(jointPostDistrib$betaSample[,4])

# Prediction intervals for each param
quantile(jointPostDistrib$sigma2Sample, c(0.025, 0.975))
quantile(jointPostDistrib$betaSample[,1], c(0.025, 0.975))
quantile(jointPostDistrib$betaSample[,2], c(0.025, 0.975))
quantile(jointPostDistrib$betaSample[,3], c(0.025, 0.975))
quantile(jointPostDistrib$betaSample[,4], c(0.025, 0.975))

## Answer: Interpretation of the credible interval for weight [-4.759964, -1.531457]. A one unit increase of weight
## lowers the amount of miles per gallon between -4.759964 and -1.531457 with 95 % posterior probability. 

## b) Investigate if effect on mpg is different in cars with six cylinders compared to cars with 8 cylinders

hist(jointPostDistrib$betaSample[,4]-jointPostDistrib$betaSample[,3], 50)
quantile(jointPostDistrib$betaSample[,4]-jointPostDistrib$betaSample[,3], c(0.025, 0.975))

## Answer: Since 0 is present in interval we can not say that there is a difference between 8 and 6 cylinders
## with 95 % posterior probability.

## c) Compute by simulation predictive distrib for a new car 4 cylinders and weight=3.5

new_x=c(1,3.5,0,0)
pred_y=rep(0,nIter)
for (i in 1:nIter) {
  pred_y[i]=sum(new_x*jointPostDistrib$betaSample[i,])+rnorm(1,sd=sqrt(jointPostDistrib$sigma2Sample[i]))
}
hist(pred_y, breaks=40, freq=FALSE)




