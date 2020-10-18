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
     xlab="Number of draws", ylab="Mean", type="l")
abline(h=MeanOfPosterior, col="red")
plot(nVector, stdVector, main="Plot of how the standard deviation converges with respect to the number of draws",
     xlab="Number of draws", ylab="Standard deviation", type="l")
abline(h=StdOfPosterior, col="red")
## Conclusion: As seen in the plot the posterior mean as well as the posterior standard deviation converges towards its true
## value of approx 0.29 and 0.09 respectively as the number of randow draws grows large.

## b) Use simulation (nDraws=10000) to compute the posterior probability Pr(theta>0.3 given y) and compare with
## with the exact value

trueProb=1-pbeta(0.3, alpha0+s, beta0+f)
set.seed(12345)
draw10000=rbeta(10000, alpha0+s, beta0+f)
probHat=sum(draw10000>0.3)/10000
print(trueProb)
print(probHat)

## Conclusion: As seen in the results from both calculations the probHat is very close to the true probability from the beta
## distribution. As the number of draws increases the approximated probability will converge towards the true
## value.

## c) Compute the posterior distribution of the log-odds phi= log(theta/(1-theta)) by simulation (nDraws=10000)

phi=log(draw10000/(1-draw10000))
hist(phi, breaks=20, main="Distribution of the log-odds")
plot(density(phi), main="Density function of phi")

  



