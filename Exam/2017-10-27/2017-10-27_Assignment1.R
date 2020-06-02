## a) Likelihood: Beta symmetric, prior, expon(1). Plot posterior distrib.

thetaGrid=seq(0.01, 15, length=1000)
data=yProp
lambda=1

logPriorExp = function(theta, lambda) {
  return(dexp(theta, rate=lambda, log=TRUE))
}

logPosterior = function(x, theta, lambda) {
  prior=logPriorExp(theta, lambda)
  likelihood=sum(dbeta(x, theta, theta, log=TRUE))
  return(likelihood+prior)
}

theta_post=sapply(thetaGrid, logPosterior, x=data, lambda=lambda)
theta_post_norm=1/((15-0.01)/1000)*exp(theta_post)/sum(exp(theta_post))
plot(thetaGrid, theta_post_norm, type="l", lwd=2, xlab=expression(theta), ylab="Posterior density")

# Zero to 1 loss means posterior mode is the optimal point estimator

index=which(theta_post_norm==max(theta_post_norm))
opt_theta=thetaGrid[index]
print(opt_theta)

## Optimal theta is around 4.481491

## b) Theta1 and theta2 are independent apriori. Plot joint posterior distrib

logPosteriorMult = function(theta, x, lambda) {
  theta1=theta[1]
  theta2=theta[2]
  prior1=logPriorExp(theta1, lambda)
  prior2=logPriorExp(theta2, lambda)
  likelihood=sum(dbeta(x, theta1, theta2, log=TRUE))
  return(likelihood+prior1+prior2)
}

# Defining initial values to be passed on to the optimizer
initVal = c(1,1)

# Finding the optimized betavector
optimResult = optim(initVal, logPosteriorMult, x=data, lambda=1, method=c("L-BFGS-B"),
                    control=list(fnscale=-1), lower=c(0.01,0.01), upper=c(Inf, Inf), hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = -solve(optimResult$hessian)
names(postMode)=c("Theta1", "Theta2")
rownames(postCov)=c("Theta1", "Theta2")
colnames(postCov)=c("Theta1", "Theta2")
print("The posterior mode is:")
print(postMode)
print("The approximated standard deviation is:")
print(postCov)

## c) Discuss how a Bayesian can determine if the symmetric model in 1a) or the non-symmetric model in 1b) 
## is most appropriate for this data. No need to compute anything here, just discuss.

## By calculating marginal likelihood for each model and check which has the highest. One can also calculate
## the bayes factor or the posterior model probabilities and choose the model with the highest probability.
