### Assignment 1

## Assignment 1: Let y1,...,yn be Bernoulli distributed with parameter theta. Assume that you have obtained a sample
## with s=5 successes in n=20 trials. Assume a Beta(alpha0, beta0) prior for theta and let alpha0=beta0=2.

## a) Draw random numbers from the posterior theta given y ~ Beta(alpha0+s, beta0+f) and verify graphically that the
## posterior mean and standard deviation converges to the true values as the number of random draws grows large. 

set.seed(12345)
alpha0=2
beta0=2
s=5
f=15
n=20

# Function for calculating the mean of a beta-distribution
calcMeanBeta = function(alpha, beta) {
  return(alpha/(alpha+beta))
}

# Function for calculating the standard deviation of a beta-distribution
calcStdDevBeta = function(alpha, beta) {
  return(sqrt(alpha*beta/((alpha+beta)^2*(alpha+beta+1))))
}

# Function for calculating the mean squared error of drawed data
calcMSE = function(n, mean, data){
  return(sqrt(1/(n-1)*sum((data-mean)^2)))
}

# Function for drawing random values from the betadistribution
drawBetaValues = function(n, alpha, beta) {
  return(rbeta(n, alpha, beta))
}

MeanOfPosterior = calcMeanBeta(alpha0+s, beta0+f)
StdOfPosterior = calcStdDevBeta(alpha0+s, beta0+f)
nVector = seq(1, 5000, 1)
meanVector=c()
stdVector=c()
for (i in nVector) {
  set.seed(12345)
  betaValues= drawBetaValues(i, alpha0+s, beta0+f)
  meanVector=c(meanVector, mean(betaValues))
  stdVector=c(stdVector, calcMSE(i, mean(betaValues), betaValues))
}
plot(nVector, meanVector, main="Plot of how the mean converges with respect to number of draws",
     xlab="Number of draws", ylab="Mean")
plot(nVector, stdVector, main="Plot of how the standard deviation converges with respect to the number of draws",
     xlab="Number of draws", ylab="Standard deviation")

## As seen in the plot the posterior mean as well as the posterior standard deviation converges towards its true
## value of approx 0.29 and 0.09 respectively as the number of randow draws grows large.

## b) Use simulation (nDraws=10000) to compute the posterior probability Pr(theta>0.3 given y) and compare with
## with the exact value

trueProb=1-pbeta(0.3, alpha0+s, beta0+f)
set.seed(12345)
draw10000=rbeta(10000, alpha0+s, beta0+f)
probHat=sum(draw10000>0.3)/10000
print(trueProb)
print(probHat)

## As seen in the results from both calculations the probHat is very close to the true probability from the beta
## distribution. As the number of draws increases the approximated probability will converge towards the true
## value.

## c) Compute the posterior distribution of the log-odds phi= log(theta/(1-theta)) by simulation (nDraws=10000)

phi=log(draw10000/(1-draw10000))
hist(phi, breaks=20, main="Distribution of the log-odds")
plot(density(phi), main="Density function of phi")

### Assignment 2

## Assignment 1: Let y1,...,yn be Bernoulli distributed with parameter theta. Assume that you have obtained a sample
## with s=5 successes in n=20 trials. Assume a Beta(alpha0, beta0) prior for theta and let alpha0=beta0=2.

## a) Draw random numbers from the posterior theta given y ~ Beta(alpha0+s, beta0+f) and verify graphically that the
## posterior mean and standard deviation converges to the true values as the number of random draws grows large. 

set.seed(12345)
alpha0=2
beta0=2
s=5
f=15
n=20

# Function for calculating the mean of a beta-distribution
calcMeanBeta = function(alpha, beta) {
  return(alpha/(alpha+beta))
}

# Function for calculating the standard deviation of a beta-distribution
calcStdDevBeta = function(alpha, beta) {
  return(sqrt(alpha*beta/((alpha+beta)^2*(alpha+beta+1))))
}

# Function for calculating the mean squared error of drawed data
calcMSE = function(n, mean, data){
  return(sqrt(1/(n-1)*sum((data-mean)^2)))
}

# Function for drawing random values from the betadistribution
drawBetaValues = function(n, alpha, beta) {
  return(rbeta(n, alpha, beta))
}

MeanOfPosterior = calcMeanBeta(alpha0+s, beta0+f)
StdOfPosterior = calcStdDevBeta(alpha0+s, beta0+f)
nVector = seq(1, 5000, 1)
meanVector=c()
stdVector=c()
for (i in nVector) {
  set.seed(12345)
  betaValues= drawBetaValues(i, alpha0+s, beta0+f)
  meanVector=c(meanVector, mean(betaValues))
  stdVector=c(stdVector, calcMSE(i, mean(betaValues), betaValues))
}
plot(nVector, meanVector, main="Plot of how the mean converges with respect to number of draws",
     xlab="Number of draws", ylab="Mean")
plot(nVector, stdVector, main="Plot of how the standard deviation converges with respect to the number of draws",
     xlab="Number of draws", ylab="Standard deviation")

## As seen in the plot the posterior mean as well as the posterior standard deviation converges towards its true
## value of approx 0.29 and 0.09 respectively as the number of randow draws grows large.

## b) Use simulation (nDraws=10000) to compute the posterior probability Pr(theta>0.3 given y) and compare with
## with the exact value

trueProb=1-pbeta(0.3, alpha0+s, beta0+f)
set.seed(12345)
draw10000=rbeta(10000, alpha0+s, beta0+f)
probHat=sum(draw10000>0.3)/10000
print(trueProb)
print(probHat)

## As seen in the results from both calculations the probHat is very close to the true probability from the beta
## distribution. As the number of draws increases the approximated probability will converge towards the true
## value.

## c) Compute the posterior distribution of the log-odds phi= log(theta/(1-theta)) by simulation (nDraws=10000)

phi=log(draw10000/(1-draw10000))
hist(phi, breaks=20, main="Distribution of the log-odds")
plot(density(phi), main="Density function of phi")

### Assignment 3

## 3.a) Bayesian inference for the concentration parameter in the von Mises distribution. This exercise is concerned
## with directional data. The point is to show you that the posterior distribution for somewhat weird models can be
## obtained by plotting it over a grid of values. The data points are observed wind directions at a given location on
## ten different days. The data are recorded in degrees: (40, 303, 326, 285, 296, 314, 20, 308, 299, 296) where North
## is located at zero degrees (see Figure 1 on the next page, where the angles are measured clockwise). To fit with 
## Wikipedias description of probability distributions for circular data we convert the data into radians -pi<=y<=pi.
## The 10 observations in radians are (-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02).
## Assume that these data points are independent observations following the von Mises distribution
## p(y given my,k) = exp(k*cos(y-u))/(2*pi*I0(k)), -pi<=y<=pi, where I0(k) is the modified Bessel function of the 
## first kind of order zero (see ?besselI in R). The parameter my (-pi<=my<=pi) is the mean direction and k>0 is
## called the concentration parameter. Large k gives a small variance around my, and vice versa. Assume that my is
## known to be 2.39. Let K ~ Exponential(Lambda=1) a priori, where lambda is the rate parameter of the exponential
## distribution (so that the mean is 1/lambda).

## a) Plot the posterior distribution of k for the wind direction data over a fine grid of k values. 

data_radian=c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
my=2.39
lambda=1

# Function for computing the vonMisesDistrib for a given dataset
vonMisesDistrib = function(kappa, data, my){
  likelihood=1
  for (i in data) {
    likelihood=likelihood*exp(kappa*cos(i-my))/(2*pi*besselI(kappa, 0))
  }
  return(likelihood)
}

# Function for computing the exponential distribution
exponDistrib = function(data, lambda) {
  return(1/lambda*exp(-1/lambda*data))
}

kappa_values=seq(0,10,0.01)

# Function for computing the posterior distribution
posteriorDistrib = function(kappa, lambda, data, my) {
  likelihood=vonMisesDistrib(kappa, data, my)
  prior=exponDistrib(kappa, lambda)
  return(likelihood*prior)
}

posteriorLikelihood=posteriorDistrib(kappa_values, lambda, data_radian, my)
plot(kappa_values, posteriorLikelihood, xlab="Kappa", ylab="Likelihood",
     main="Posterior likelihood for different kappavalues", type="l", col="blue")

## As seen in the plot the likelihood of the posterior peaks between 2 and 4 and then dies off for larger
## kappa-values.

## b) Find the (approximate) posterior mode of k from the information in a).

# Puts likelihood values with corresponding kappa-values to be able to retrieve the kappa-value corresponding to
## the highest likelihood (mode)
posterior.df=data.frame(kappa=kappa_values, likelihood=posteriorLikelihood)
posteriorMode=subset(posterior.df, likelihood==max(likelihood), kappa)
print(posteriorMode$kappa)

## The approximated posterior mode is found to be 2.12.









