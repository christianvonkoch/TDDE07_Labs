## a) Data is normally distributed, assume non-informative prior. Simulate 1000 draws from the predictive distribution
## of the maximal weight in a given future week and plot them.

par(mfrow=c(1,1))
data=c(1690, 1790, 1760, 1750)
n=length(data)
sigma2=50^2

# Prior assumed to be constant. This yields posterior distribution of mu~N(mean(data), simga2/n) according to L2

postDistrib=rnorm(1000, mean=mean(data), sd=sqrt(sigma2/n))
predDistrib=rnorm(1000, mean=postDistrib, sd=sqrt(sigma2))
hist(predDistrib, breaks=100, main="Approximated predictive distribution of mu", xlab=expression(mu), freq=FALSE)

# To check if reasonable the real predictive distribution is plotted. We know from L4 that the predictive distribution
# of new obs is distributed N(mean(data), sigma2*(1+1/n))

grid=seq(1500,1900)
lines(grid, dnorm(grid, mean=mean(data), sd=sqrt(sigma2*(1+1/n))), col="red")

## Since the histogram follows the real distribution well it was performed correctly. 

## b) Use simulation to approximate the expected number of weeks out of the coming 52 weeks in which the maximal 
## weight will exceed 1850 kg, based on the predictive distribution. 

nDraws=1000
weekMatrix=matrix(0,52,nDraws)
for (i in 1:nDraws) {
  weekMatrix[,i]=rnorm(52, mean=mean(data), sd=sqrt(sigma2*(1+1/n)))
}
countWeeks=colSums(weekMatrix>1850)
barplot(table(countWeeks), main="Approximated predictive distribution", xlab="No. of weeks")
mean(countWeeks)

## Important here to simulate the number of predictive draws taken, in this case 52 samples. Then sum the no of
## observations in each sample which satisfies the condition. This then becomes the predictive distribution.
## We can then take the mean out of this sample to obtain the expected number of weeks. 

## c) The weight that the escalator can hold at any given time is given by 1000log(a), a is the build cost.
## If the weight is exceeded the excalator breaks and has to be repaired. Loss function for shopping mall is
## L(a, theta) = a+n(a,theta) where n(a,theta) is the no. of weeks out of the 52 in which the escalator breaks.
## Compute the optimal build cost (a) using Bayesian approach.

# Want to maximize the negative loss function.

countOfBreak=function(a, countMatrix) {
  return(colSums(countMatrix>1000*log(a)))
}

utilityFunction = function(a, n) {
  return(-(a+n))
}

aGrid=seq(0,20, 0.001)
utility=rep(0,length(aGrid))
for(i in 1:length(aGrid)) {
  counts=countOfBreak(aGrid[i], weekMatrix)
  utility[i]=mean(utilityFunction(aGrid[i], counts))
}

plot(aGrid, utility, type="l", lwd=2, main="Utility function")
aOpt=aGrid[which(utility==max(utility))]
points(aOpt, max(utility), col="red", cex=2, lwd=2)
aOpt

## 6.74 yields maximum utility and equivalently minimum loss.
