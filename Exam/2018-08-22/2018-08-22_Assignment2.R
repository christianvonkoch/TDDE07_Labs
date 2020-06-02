## a) Use numeric optimization to approximate joint posterior of beta.

library(mvtnorm)

data=titanic
y=data$survived
X=as.matrix(data[,-1])


nFeatures = dim(X)[2]
covNames=names(data[,2:ncol(data)])

# Constructing prior
tau=50
mu_prior = rep(0,nFeatures)
sigma_prior = tau^2*diag(nFeatures) 

logPostLogistic = function(beta, Y, X, mu, sigma) {
  nFeat = length(beta)
  XBeta=X%*%beta
  # Defining loglikelihood
  logLike = sum(Y*XBeta-log(1+exp(XBeta)))
  if (abs(logLike) == Inf) logLik = -20000; # Likelihood is not finite, stear the optimizer away from here!
  # Defining prior
  prior = dmvnorm(beta, mean=mu, sigma=sigma, log=TRUE)
  # Adding loglikelihood and logprior together. Since it is log both of them are added instead of multiplied
  return(logLike + prior)
}

# Defining initial values to be passed on to the optimizer
set.seed(12345)
initVals = rnorm(dim(X)[2])

# Finding the optimized betavector
optimResult = optim(initVals, logPostLogistic, Y=y, X=X, mu=mu_prior, sigma=sigma_prior, method=c("BFGS"),
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

par(mfrow=c(2,2))
for (i in 2:nFeatures) {
  grid=seq(postMode[i]-3*approx_PostStd[i], postMode[i]+3*approx_PostStd[i], length=1000)
  plot(grid, dnorm(grid, mean=postMode[i], sd=approx_PostStd[i]), 
       main=paste("Marginal posterior of", covNames[i]), xlab=covNames[i], ylab="Density", type="l", lwd=2)
}

## b) Compute posterior probability that the adult feature is smaller than 0

prob=pnorm(0, mean=postMode[which(covNames=="adult")], sd=approx_PostStd[which(covNames=="adult")])
prob

## The interpretation of the probability 0.76 is that we can say with approximately 76 % posterior probability
## that being an adult contributed negatively in regards to survival on the titanic. If you were an adult 
## you were more probable to die than if not.

## c) A first class adult woman and a third class adult man are together during the disaster.
## Compute predictive probability that the woman survives but the man dies. 

man=c(1,1,1,0,0)
woman=c(1,1,0,1,0)
nDraws=5000
results=matrix(0,nDraws, nFeatures)
for (i in 1:nDraws) {
  results[i,]=rmvnorm(1, mean=postMode, sigma=postCov)
}

manPred=results%*%man
womanPred=results%*%woman
manSim=rbinom(nDraws, 1, exp(manPred)/(1+exp(manPred)))
womanSim=rbinom(nDraws, 1, exp(womanPred)/(1+exp(womanPred)))
final=ifelse(womanSim == 1 & manSim ==0, 1, 0)
mean(final)

## Reasonable. Do simulation of param, use that param in new obs likelihood. Check probability.

