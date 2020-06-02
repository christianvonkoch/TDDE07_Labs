## Using dataset ebay. Data describing the number of bids of 100 ebay auctions. 
## a) Assume bimodial model with N=50. Use prior (theta-1)^2.
## Write function in R that computes (unnormalized) log posterior density function of theta. Use function to 
## plot normalized posterior density function of theta on the interval [0,1] with at least 1000 grid points. 
## Report the (approximate) value of the posterior mode based on the computed values needed for the plot.

N=50
n=100
data=ebay

logPrior = function(theta) {
  return(log((theta-1)^2))
}

# The likelihood is proportional to the Beta(sum(data)+1, Nn-sum(data)+1) density
logLike = function(data, theta, N=50) {
  return(dbeta(theta, sum(data)+1, N*length(data)-sum(data)+1, log=TRUE))
}

logPost = function(data, theta, N=50) {
  log_Prior=logPrior(theta)
  log_Like=logLike(data, theta, N)
  return(log_Prior+log_Like)
}

thetaGrid=seq(0,1,0.001)
post_theta=exp(sapply(thetaGrid, logPost, data=data, N=50))
post_theta_norm=1/0.001*post_theta/sum(post_theta)
plot(thetaGrid, post_theta_norm, type="l", lwd=2, main="Approximated posterior density of theta", 
     xlab=expression(theta), ylab="Density")
postMode=thetaGrid[which(post_theta_norm==max(post_theta_norm))]
print(postMode)     
abline(v=postMode, col="red", lty=2)
title(sub="Black = Density, Red = Posterior mode")

## b) Use supplied function GibbsMixPoisin file ExamData to do Gibbs sampling for a mixture of Poissons model
## where each data pointis modeled as independent with density given.

set.seed(100)
K=2
nIter=500
xGrid=seq(min(data), max(data))
results=GibbsMixPois(ebay, K, alpha=1, alphaGamma=1, betaGamma=1, xGrid=xGrid, nIter)
post_theta=results$thetaSample
cum_mean=matrix(0,nIter,2)
theta1_cumsum=cumsum(post_theta[,1])
theta2_cumsum=cumsum(post_theta[,2])
for (i in 1:nIter) {
  cum_mean[i,]=c(theta1_cumsum[i]/i, theta2_cumsum[i]/i)
}
par(mfrow=c(2,1))
plot(seq(1,500), post_theta[,1], xlab="No. of bids", ylab=expression(theta), main="Trace plot", type="l")
title(line=3, main="Convergence of sampler for theta1")
plot(seq(1,500), cum_mean[,1], xlab="No. of bids", ylab="Cumulative mean", type="l", main="Cumulative means")
plot(seq(1,500), post_theta[,2], xlab="No. of bids", ylab=expression(theta), main="Trajectory over theta2", type="l")
title(line=3, main="Convergence of sampler for theta2")
plot(seq(1,500), cum_mean[,2], xlab="No. of bids", ylab="Cumulative mean", type="l", main="Cumulative means")

## According to plots we should choose burnin = 50 approximately. 

## c) Use graphical methods to investigate if mixture of Poissons with K=2 fits data well.

data_norm=as.vector(bidsCounts/sum(bidsCounts))
postMean=results$mixDensMean
par(mfrow=c(1,1))
plot(xGrid, data_norm, xlab="No. of bids", ylab="Density", main="Fitted models", type="o", lwd=2, ylim=c(0,0.3))
lines(xGrid, postMean, col="red", lwd=1, type="o", lty=2)
lines(xGrid, dbinom(xGrid, 50, postMode), col="blue", lwd=1, type="o", lty=2)
legend("topright", col=c("black", "red", "blue"), 
       legend=c("Data", "Posterior mean of mixture model", "Binomial model"), 
       lty=c(1,2,2), lwd=c(2,1,1), pch=c("o", "o", "o"))

## I would recommend mixture of poissons since it fits data better. The binomial model is clearly worse 
## specifically when we look at the number of 0 bids which is higher in the data.


