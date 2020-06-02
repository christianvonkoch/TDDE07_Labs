## a) Plot the posterior distribution of theta

riceData <- c(1.556, 1.861, 3.135, 1.311, 1.877, 0.622, 3.219, 0.768, 2.358, 2.056)

# Random number generator for the Rice distribution
rRice <-function(n = 1, theta = 1, psi = 1){
  x <- rnorm(n = n, mean = 0, sd = sqrt(psi))
  y <- rnorm(n = n, mean = theta, sd = sqrt(psi))
  return(sqrt(x^2+y^2))
}

# Function for calculating the log posterior distrib with theta prior set to 1
logPosterior = function(data, theta, psi) {
  bessel_factor=1
  for (i in data) {
    bessel_factor=bessel_factor*besselI(i*theta/psi, nu=0)
  }
  post=-log(psi)-1/(2*psi)*sum(data^2+theta^2)+log(bessel_factor)
  return(post+0) # If prior is assumed to be constant we set the prior to 1 which in log scale yields 0
}

gridWidth=0.01
theta_grid=seq(0,3,gridWidth)
posterior_distrib_log=sapply(theta_grid, logPosterior, data=riceData, psi=1)
posterior_distrib_norm=1/gridWidth*exp(posterior_distrib_log)/sum(exp(posterior_distrib_log))
sum(posterior_distrib_norm)
plot(theta_grid, posterior_distrib_norm, xlab=expression(theta), ylab="Density", main="Posterior density of theta",
     type="l", lwd=2)

## b) Use numerical optimization to obtain a normal approx. of the posterior distrib of theta. Overlay curve 
## from a) with the approximated normal distribution

# Defining initial values to be passed on to the optimizer
set.seed(12345)
initVal = rnorm(1, mean=0, sd=1)

# Finding the optimized betavector
optimResult = optim(initVal, logPosterior, data=riceData, psi=1, method=c("L-BFGS-B"),
                    control=list(fnscale=-1), lower=0, hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = as.numeric(-solve(optimResult$hessian))
print("The posterior mode is:")
print(postMode)
print("The approximated standard deviation is:")
print(postCov)
lines(theta_grid, dnorm(theta_grid, mean=postMode, sd=sqrt(postCov)), col="red", lwd=2)
legend(x = 1.8, y = 1, legend = c("True posterior", "Approximate posterior"), 
       col = c("black","red"), lty = c(1,1), lwd = c(2,2), cex = 0.8)

## Answer: Not perfect approx but fairly good. 

## c) Simulate distrib for new observation using normal approx in b)

nDraws=5000
set.seed(12345)
theta=rnorm(nDraws, mean=postMode, sd=sqrt(postCov))
pred_distrib=c()
for (i in theta) {
  pred_distrib=c(pred_distrib, rRice(theta=i))
}

hist(pred_distrib, breaks=100, xlab="Index", main="Predictive density of new obs")