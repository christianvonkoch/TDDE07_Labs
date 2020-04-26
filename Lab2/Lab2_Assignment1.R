## Assignment 1: The dataset TempLinkoping.txt contains daily average tamperatures (in Celcius degrees) at 
## Malmslatt, Linkoping over the course of the year 2018. The response variable is temp and the covariate is 
## time=(the number of days since beginning of year)/365
## You're task is to perform a Bayesian analysis of a quadratic regression
## temp=beta0+beta1*time+beta2*time^2+epsilon, epsilon~N(0,sigma^2)

## a)  Determining the prior distribution of the model parameters. Use the conjugate prior for the linear 
## regression model. Your task is to set the prior hyperparameters my0, omega0, v0 and sigma0^2 to sensible 
## values. Start with my0=(-10,100,-100)T, omega0=0.01*I3, v0=4 and sigma0^2=1.  0 = 1. Check if this prior
## agrees with your prior opinions by simulating draws from the joint prior of all parameters and for every draw
## compute the regression curve. This gives a collection of regression curves, one for each draw from the prior.
## Do the collection of curves look reasonable? If not, change the prior hyperparameters until the collection
## of prior regression curves agrees with your prior beliefs about the regression curve. [Hint: the R package
## mvtnorm will be handy. And use your Inv-chisquared simulator from Lab1.

# Read file
temp = read.table("TempLinkoping.txt", header=TRUE)

## install.packages("mvtnorm")
library(mvtnorm)
# Defining the parameters for the prior distribution
# Switched to beta0=0 since it seems more reasonable and -10 seems too low.
my0=c(0,100,-100)
omega0=0.01*diag(3)
v0=4
sigma0_sq=1
omega0Inv=solve(omega0)

calcRegr = function(betaMatrix, row, x) {
  return(betaMatrix[row,1]+betaMatrix[row,2]*x+betaMatrix[row,3]*x^2)
}

drawBeta = function(my, sigma_sq, omegaInv) {
  return(rmvnorm(1, mean=my, sigma=sigma_sq*omegaInv))
}

nDraws=1000
set.seed(12345)
drawX=rchisq(nDraws, v0)
sigma_sq=(v0)*sigma0_sq/drawX
betaMatrix=matrix(0,nDraws,3)
plot.new()
plot.window(xlim=c(0,1), ylim=c(-50,50))
axis(side=1)
axis(side=2)
set.seed(12345)
for (i in 1:nDraws) {
  betaMatrix[i,]=drawBeta(my0, sigma_sq[i], omega0Inv)
  lines(temp$time, calcRegr(betaMatrix, i, temp$time), col=rgb(0,0,0,0.2))
}
title(main="Temps depending on different times for different simulated models", xlab="Time", ylab="Temp")

## The collection of curves look reasonable and in line with our prior beliefs. The temperature rises during the
## summer months and stays low in the beginning and the end of the year respectively.However, the value of -10
## were switched to 0 since it seems more reasonable with a measurement of the temperature 0 on the 1st of
## January than a measurement of -10.

## b) Write a program that simulates from the joint posterior distribution of beta0, beta1, beta2 and sigma^2.
## Plot the marginal posteriors of each parameter as a histogram. Also produce another figure with a scatter plot
## of the temperature data and overlay a curve for the posterior median of the regression function
## f(time)=beta0+beta1*time+beta2*time^2, computed for every value of time. Also overlay curves for the lower
## 2.5% and upper 97.5% posterior credible interval for f(time). That is, compute the 95% equal tail posterior
## probability intervals for every value of time and then connect the lower and upper limits of the interval by
## curves. Does the interval bands contain most of the data points? Should they?

# Calculating the parameters for the posterior distribution
v_n=v0+length(temp$temp)
X=cbind(1, temp$time, temp$time^2)
Y=temp$temp
beta_hat=solve(t(X)%*%X)%*%t(X)%*%Y
my_n=solve(t(X)%*%X+omega0)%*%(t(X)%*%X%*%beta_hat+omega0%*%my0)
omega_n=t(X)%*%X+omega0
omega_n_Inv=solve(omega_n)
sigma_sq_n=(v0*sigma0_sq+(t(Y)%*%Y+t(my0)%*%omega0%*%my0-t(my_n)%*%omega_n%*%my_n))/v_n

# Simulate the joint posterior
sigma_sq_post=(v_n)*c(sigma_sq_n)/drawX
betaMatrix_post=matrix(0,nDraws,3)
response_post_temp=matrix(0,nDraws,length(temp$time))
for (i in 1:nDraws) {
  betaMatrix_post[i,]=drawBeta(my_n, sigma_sq_post[i], omega_n_Inv)
}
hist(betaMatrix_post[,1], breaks=100, main="Marginal posterior for beta0")
hist(betaMatrix_post[,2], breaks=100, main="Marginal posterior for beta1")
hist(betaMatrix_post[,3], breaks=100, main="Marginal posterior for beta2")

plot(temp$time, Y, main="Plot of the temp data for different times", col="blue", 
     xlab="Time coefficient", ylab="Temp")
for (i in 1:nDraws) {
  betaTemp=sapply(temp$time, calcRegr, betaMatrix=betaMatrix_post, row=i)
  response_post_temp[i,]=betaTemp
}

response_post=c()
credInterval=matrix(0, length(temp$time), 2)
for (i in 1:length(temp$time)) {
  sortedTemp=sort(response_post_temp[,i])
  response_post=c(response_post, (sortedTemp[500]+sortedTemp[501])/2)
  credInterval[i,]=c(sortedTemp[(0.025*length(sortedTemp)+1)], sortedTemp[(0.975*length(sortedTemp))])
}

lines(temp$time, response_post)
lines(temp$time, credInterval[,1], lty=21, col="gray")
lines(temp$time, credInterval[,2], lty=21, col="gray")
title(sub="Grey = 95 % credible intervals, Black = Median")

## The interval bands contain most of the data points. They should contain most of the data points if the model
## is accurate in terms of describing the reality. In this case, it seems like the model has captured most of
## the data points which means that the model describes the reality fairly well. 

## c) It is of interest to locate the time with the highest expected temperature (that is, the time where
## f(time) is maximal). Let's call this value xtilde. Use the simulations in b) to simulate from posterior
## distribution of xtilde. [Hint: The regression curve is quadratic. You can find a simple formula for xtilde
## given beta0, beta1 and beta2]

calcMaxTemp = function(betaMatrix, row) {
  return(-betaMatrix[row,2]/(2*betaMatrix[row,3]))
}

time_max_temp=c()
for (i in 1:nDraws) {
  time_max_temp=c(time_max_temp, calcMaxTemp(betaMatrix_post, i))
}

hist(time_max_temp, breaks=1000, xlim=c(0,1), main="Frequency of max temperatures simulated from xtilde",
     xlab="Temperature")

## As seen in the histogram the derived highest temperature from the simulated models is mostly present in late
## june which seems reasonable if applying to Motala in Sweden where the temperature is the highest during the 
## summer time. 

## d) Say now that you want to estimate a polynomial model of order 7, but you suspect that higher order terms
## may not be needed, and you worry about overfitting. Suggest a suitable prior that mitigates this potential
## problem. You do not need to compute the posterior, just write down your prior. [Hint: the task is to specify
## my0 and omega0 in a smart way.]

## A suitable prior for this task would be to set my0 to 0 since you want most of the coefficients close to zero
## to obtain increased shrinkage. You would also want to set omega0 to Lambda*IdentityMatrix. This would mean
## that for larger values of lambda more and more of the beta values would be close to zero since the spread of
## the distribution of the beta values would decrease. In this case, where there is a worry about overfitting,
## it might be a good idea to choose a large lambda to decrease the spread of the beta values and increase the
## probability that most of the beta values are around 0. 