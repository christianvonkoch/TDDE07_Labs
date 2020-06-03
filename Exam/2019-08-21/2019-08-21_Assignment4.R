## Aircraft incidents assumed to be independent, follow negative binomial distrib. Assume joint prior 
## 1/phi^2
## a) Simulate from posterior using Metropolis algorithm. Denote theta=c(mu, phi) and use as proposal dens
## the multivariate normal density (random walk metropolis).

# Load airline incidents data
load(file = 'incidents.RData')
data=incidents$incidents
library(mvtnorm)

nIter=1000
burnIn=50
theta_0=c(200,20)
c=0.1
postCov=diag(c(100,5))

# Defining function for sampling through metropolishastings
RVMSampler = function(previousVal, postCov, c, myFunction, ...) {
  proposalVal=rmvnorm(1, mean=previousVal, sigma=c*postCov)
  proposalVal[proposalVal<=0]=1e-6
  alpha=min(1, exp(myFunction(proposalVal,...)-myFunction(previousVal, ...)))
  u=runif(1)
  if(u < alpha) {
    return(list(theta=proposalVal, acceptProb=alpha))
  } else {
    return(list(theta=previousVal, acceptProb=alpha))
  }
}

logPrior = function(phi) {
  return(-2*log(phi))
}

logLike <- function(param, x){
  theta1 = param[1]
  theta2 = param[2]
  logPost = sum(logdNegBin(x, theta1, theta2))  - 2*log(theta2)
  return(logPost)
}

logPost = function(theta, data) {
  log_Prior=logPrior(theta[2])
  log_Like=logLike(theta, data)
  return(log_Prior+log_Like)
}

post_matrix = matrix(0, nIter+burnIn, 2)
# Setting initial values of beta to same initVals as in the optimizer (taken randomly from normal distrib)
post_matrix[1,]=theta_0
accProb=rep(0, nIter)
set.seed(12345)

for(i in 1:(nIter+burnIn)) {
  if(i<(nIter+burnIn)) {
    draw=RVMSampler(post_matrix[i,], postCov, c, logPost, data)
    post_matrix[i+1,]=draw$theta
    accProb[i+1]=draw$acceptProb
  }
}

iter=seq(1,nIter+burnIn,1)
plot(iter[-(1:burnIn)], post_matrix[-(1:burnIn),1], type="l", lwd=1, col="grey", main="Traceplot of mu in RVM",
     xlab=expression(mu), ylab="Value")
plot(iter[-(1:burnIn)], post_matrix[-(1:burnIn),2], type="l", lwd=1, col="grey", main="Traceplot of phi in RVM",
     xlab=expression(phi), ylab="Value")
mean(accProb)

## This MCMC sampler is not efficient since it moves very slowly and is therefore probably not exploring
## the whole posterior distribution.We can also see that the acceptance probability for this algorithm
## is around 84,4 % and it should be around 30 %. Once could tune the c param to lower the acceptance probability.
## One example is to increase c to a value of 3 which would yield in approximately 30 % acceptance rate. 

## b) Instead simulate from posterior using metropolis hastings. 

c=0.8

MHSampler = function(previousVal, postCov, c, myFunction, ...) {
  proposalVal_mu=rgamma(1, c*previousVal[1], c)
  proposalVal_phi=rgamma(1, c*previousVal[2], c)
  proposalVal=c(proposalVal_mu, proposalVal_phi)
  proposalVal[proposalVal<=0]=1e-6
  alpha=min(1, exp(myFunction(proposalVal,...)-myFunction(previousVal, ...)+
                     dgamma(previousVal[1], c*proposalVal[1], c)+dgamma(previousVal[2],c*proposalVal[2],c)-
                     dgamma(proposalVal[1], c*previousVal[1], c)-dgamma(proposalVal[2], c*proposalVal[2],c)))
  u=runif(1)
  if(u < alpha) {
    return(list(theta=proposalVal, acceptProb=alpha))
  } else {
    return(list(theta=previousVal, acceptProb=alpha))
  }
}

post_matrix2 = matrix(0, nIter+burnIn, 2)
theta_0=c(200,20)
post_matrix2[1,]=theta_0
accProb2=rep(0, nIter)
set.seed(12345)

for(i in 1:(nIter+burnIn)) {
  if(i<(nIter+burnIn)) {
    draw=MHSampler(post_matrix2[i,], postCov, c, logPost, data)
    post_matrix2[i+1,]=draw$theta
    accProb2[i+1]=draw$acceptProb
  }
}

plot(iter[-(1:burnIn)], post_matrix2[-(1:burnIn),1], type="l", lwd=1, col="grey", main="Traceplot of mu in MH",
     xlab=expression(mu), ylab="Value")
plot(iter[-(1:burnIn)], post_matrix2[-(1:burnIn),2], type="l", lwd=1, col="grey", main="Traceplot of phi in MH",
     xlab=expression(phi), ylab="Value")
mean(accProb2)

## The new algorithm seems to rapidly explore the posterior which is good. The acceptance probability is also lower
## around 30 % which also indicates that this algorithm is better than the previous one. 
