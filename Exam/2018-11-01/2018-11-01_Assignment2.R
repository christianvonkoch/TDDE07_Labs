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

data=fish
y=data[,1]
X=as.matrix(data[,-1])
X=cbind(X, subset(X, select=c("age", "temp"))^2)
X=cbind(X, X[,2]*X[,3])
covNames=names(data[,-1])
covNames=append(covNames, c("age^2", "temp^2", "age*temp"))
mu_0=rep(0, ncol(X))
omega_0=0.01*diag(ncol(X))
v_0=1
sigma2_0=10000
nIter=5000

bayes_lin_results=BayesLinReg(y, X, mu_0, omega_0, v_0, sigma2_0, nIter)

## a) Compute posterior mean and 95 % equal tail credible intervals for all beta-params

results=matrix(0,ncol(X)+1,3)
for (i in 1:ncol(X)) {
  results[i,1]=mean(bayes_lin_results$betaSample[,i])
  results[i,-1]=quantile(bayes_lin_results$betaSample[,i], probs=c(0.025, 0.975))
}
results[ncol(X)+1,1]=mean(bayes_lin_results$sigma2Sample)
results[ncol(X)+1,-1]=quantile(bayes_lin_results$sigma2Sample, probs=c(0.025, 0.975))
covNames=append(covNames, "sigma2")
rownames(results)=covNames
colnames(results)=c("Posterior mean", "2,5%", "97,5%")
results

## b) Compute the posterior mean and posterior median of the noise standard deviation theta

median(bayes_lin_results$sigma2Sample)
results[ncol(X)+1,1]

## Results shown above

## c) First eleven datapoints come from watertank with 25 degrees celsius. Produce scatter plot of these datapoints
## with length and age on the two axes. Overlay a curve for the posterior mean of the regression curve with respect
## to age. 

betaMatrix=bayes_lin_results$betaSample
tempData=data[1:11,]
plot(tempData$length, tempData$age, main="Plot of data with 25 degrees temperature in tank", xlab="Length",
     ylab="Age", col="blue")
ageGrid=seq(0,160, 0.01)
credInt=matrix(0,length(ageGrid),2)
fAgePostMean=rep(0,length(ageGrid))
fAge=rep(0,nIter)
count=1
for (a in ageGrid) {
  fAge=betaMatrix%*%c(1,a,25,a^2, 25^2, a*25)
  fAgePostMean[count]=mean(fAge)
  credInt[count,]=quantile(fAge, probs=c(0.025, 0.975))
  count=count+1
}
lines(fAgePostMean, ageGrid, type="l", lwd=2, col="red")
lines(credInt[,1], ageGrid, col="grey", lty=2)
lines(credInt[,2], ageGrid, col="grey", lty=2)

## d) Assume that you want to make predictions for fish in a new water tank with a temperature of 15 degrees celsius
## , which is lower than any of the temperatures in the original data set. Discuss how the current data set
## might be a problem regarding this matter and how the prior could be changed to control this problem. 

## Since we have a model with high order terms the risk for overfitting is bigger than with a simpler model.
## When testing the model on data which are far away from the data used to fit the model, it is of high risk
## that the model might perform badly. To reduce the risk of overfitting one can use a betaprior close to zero
## to force many of the covariates to become zero. To further reduce the risk for overfitting one can increase the
## values in the diagonal of omega_0 to larger values to further reduce the variance. This should be done for the 
## covariates and not for the intercept so the place [1,1] in the omega_0 matrix can remain the same but the other
## values can be increased. 
