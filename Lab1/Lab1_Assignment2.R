## Assignment 2: Assume that you have asked 10 randomly selected persons about their monthly 
## income(inthousandsSwedishKrona)andobtainedthefollowingtenobservations: 44, 25, 45, 52, 30, 63, 19, 50, 34 and 67.
## A common model for non-negative continuous variables is the log-normal distribution. The log-normal distribution 
## log(N(my, sigma^2)) has density function ... for y > 0, my > 0 and sigma > 0. The log-normal distribution is
## related to the normal distribution as follows: if y ~ log N(my, sigma^2) then log y ~ N(my, sigma^2). Let
## y1,...,yn given my and simga^2 ~ log N(my, sigma^2), where my=3.7 is assumed to be known but sigma^2 is unknown
## with noninformative prior p(sigma^2) is proportional to 1/sigma^2. The posterior for sigma^2 is the Inv - 
## chitwo distribution with X(n, thao^2) distribution, where thao^2=sum((log(yi)-my)^2)/n.

## a) Simulate 10 000 draws from the posterior of sigma^2 (assuming my=3.7) and compare with the theoretical 
## with the theoretical Inv - chitwo distribution with X(n, thao^2) posterior distribution.

x=c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
n=length(x)
my=3.7

#Function for calculating thao^2
calcThao = function(data, my, n) {
  return(sum((log(data)-my)^2)/n)
}

thaosq=calcThao(x, my, n)
thetaVector=c()
for(i in 1:10000) {
  drawX=rchisq(1, n-1)
  sigmasq=(n-1)*thaosq/drawX
  thetaVector=c(thetaVector, rnorm(1, mean(x), sqrt(sigmasq/n)))
}
plot(density(thetaVector), main="Density function of the simulation of the posterior for sigma^2")



