## a) Plot posterior density of theta, with normal prior and cauchy distrib as likelihood

# Reading the data vector yVect from file
load(file = 'CauchyData.RData')
cauchydata=yVect

dCauchy <- function(x, theta = 0, gamma = 1){
  return(dens = (1/(pi*gamma))*(1/(1+((x-theta)/gamma)^2)))
}

dlognormal <- function(x, mu, sigma2){
  return(dens = (1/(sqrt(2*pi*sigma2)*x))*exp((-1/(2*sigma2))*(log(x)-mu)^2))
}

logPrior_theta = function(theta, mu, sigma_sq) {
  return(dnorm(theta, mean=mu, sd=sqrt(sigma_sq), log=TRUE))
}

logPosterior = function(data, mu, sigma_sq, theta=0, gamma=1) {
  prior=logPrior(theta, mu, sigma_sq)
  likelihood=dCauchy(data, theta, gamma)
  likelihood=sum(log(likelihood))
  return(likelihood + prior)
}

mu=0
sigma_sq=100
gamma=1
gridWidth=0.01
theta_grid=seq(0,8,gridWidth)
posterior_distrib=sapply(theta_grid, logPosterior, data=cauchydata, mu=mu, sigma_sq=sigma_sq, gamma=1)
posterior_distrib=1/gridWidth*exp(posterior_distrib)/sum(exp(posterior_distrib))
plot(theta_grid, posterior_distrib, type="l", lwd=2, main="Posterior density for theta", xlab=expression(theta),
     ylab="Density")

## b) gamma is unknown with prior lognormal. 

set.seed(12345)
initVal = c(0,0)

logJointPosterior = function(joint, data, mu, sigma_sq) {
  prior_theta=logPrior(joint[1], mu, sigma_sq)
  prior_gamma=log(dlognormal(joint[2], mu, sigma_sq))
  likelihood=dCauchy(data, joint[1], joint[2])
  likelihood=sum(log(likelihood))
  return(likelihood + prior_theta + prior_gamma)
}

# Finding the optimized theta and gamma
optimResult = optim(initVal, logJointPosterior, data=cauchydata, mu=mu, sigma_sq=sigma_sq, method=c("L-BFGS-B"),
                    control=list(fnscale=-1), lower=c(-Inf, 0.001), upper=c(Inf, Inf), hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = -solve(optimResult$hessian)
names(postMode)=c("Theta", "Gamma")
print("The posterior mode is:")
print(postMode)
print("The approximated standard deviation is:")
print(postCov)

## c) Use normal approx in 1b) to obtain marginal posterior for the 99 % percentile of the caucy distrib
## theta + gamma * tan(pi(0.99-0.5))

library(rmvnorm)
normal_approx=rmvnorm(5000, mean=postMode, sigma=postCov)
cauchy_distrib=normal_approx[,1]+normal_approx[,2]*tan(pi*(0.99-0.5))
hist(cauchy_distrib, breaks=50, main="Marginal distribution of special case of caucby", xlab="Function value")
