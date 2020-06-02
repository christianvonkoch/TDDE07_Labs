## c) Calc unnormalized posterior and plot normalized posterior. Gamma prior and indep likelihoods.

gridWidth=0.01
thetaGrid=seq(0,2,gridWidth)
xData <- c(1.888, 2.954, 0.364, 0.349, 1.090, 7.237)
yData <- c(-1.246, -1.139, -0.358, -1.308, -0.930, -0.157, -0.111, -0.635)
alpha=3
beta=2

logPosteriorX = function(theta, alpha, beta) {
  return(dgamma(theta, alpha, beta, log=TRUE))
}

likeY = function(y, theta) {
  return(-3*sum(log(1+(1/5)*(y-log(theta))^2)))
}

logPosterior = function(theta, alpha, beta, xDat, yDat) {
  likelihoodY=likeY(yDat, theta)
  logPostX=logPosteriorX(theta, length(xDat+3), sum(xDat)+2)
  return(likelihoodY+logPostX)
}

post_theta=sapply(thetaGrid, logPosterior, alpha=alpha, beta=beta, xDat=xData, yDat=yData)
post_theta_norm=1/gridWidth*exp(post_theta)/sum(exp(post_theta))
plot(thetaGrid, post_theta_norm, type="l", lwd=2, main="Posterior of theta", xlab=expression(theta),
     ylab="Density")


