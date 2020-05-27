## Assignment 2: Assume that you have asked 10 randomly selected persons about their monthly 
## income(inthousandsSwedishKrona)andobtainedthefollowingtenobservations: 44, 25, 45, 52, 30, 63, 19, 50, 34 
## and 67. A common model for non-negative continuous variables is the log-normal distribution. The log-normal
## distribution log(N(my, sigma^2)) has density function ... for y > 0, my > 0 and sigma > 0. The log-normal
## distribution is related to the normal distribution as follows: if y ~ log N(my, sigma^2) then 
## log y ~ N(my, sigma^2). Let y1,...,yn given my and simga^2 ~ log N(my, sigma^2), where my=3.7 is assumed to be
## known but sigma^2 is unknown with noninformative prior p(sigma^2) is proportional to 1/sigma^2. The posterior
## for sigma^2 is the Inv - chitwo distribution with X(n, thao^2) distribution, where thao^2=sum((log(yi)-my)^2)/n

## a) Simulate 10 000 draws from the posterior of sigma^2 (assuming my=3.7) and compare with the theoretical 
## with the theoretical Inv - chitwo distribution with X(n, thao^2) posterior distribution.

library(geoR)
x=c(44, 25, 45, 52, 30, 63, 19, 50, 34, 67)
n=length(x)
my=3.7

#Function for calculating thao^2
calcThao = function(data, my, n) {
  return(sum((log(data)-my)^2)/n)
}

invchisquare <- function(x, df, taosq){
  first = ((taosq*df/2)^(df/2))/gamma(df/2)
  second = (exp((-df*taosq)/(2*x)))/(x^(1+df/2))
  return(first*second)
}

thaosq=calcThao(x, my, n)
set.seed(12345)
drawX=rchisq(10000, n)
sigmasq=(n)*thaosq/drawX
xvals=seq(0.001, 3, 0.001)
plot(density(sigmasq), main="Density of simulated sigma^2, black = simulated distrib., red = actual distrib.")
lines(xvals,invchisquare(xvals, n, thaosq), col="red")

## As seen in the plot the theoretical distribution (red line) follows the simulated one with good precision. This
## indicates that the simulation has been made correctly.

## b) The most common measure of income inequality is the Gini coefficient, G, where 0<=G<=1. G=0 means a 
## completely equal income distribution, whereas G=1 means complete income inequality. See Wikipedia for more
## information. It can be shown that G=2*CDF-normal(sigma/sqrt(2))-1 when income follow a log N(my, sigma^2)
## distribution.  Use the posterior draws in a) to compute the posterior distribution of the Gini coefficient G
## for the current data set.

G=2*pnorm(sqrt(sigmasq/2), mean=0, sd=1)-1
hist(G, breaks=100)
plot(density(G), main="Density function of simulated values of the Gini coefficient")

## As seen in the plot the gini coefficient is centered at around 0.2 which means a rather inequal distribution.

## c) Use the posterior draws from b) to compute a 90% equal tail credible interval for G. A 90% equal tail interval
## (a,b) cuts off 5% percent of the posterior probability mass to the left of a, and 5% to the right of b. Also, 
## do a kernel density estimate of the posterior of G using the density function in R with defaultsettings, 
## and use that kernel density estimate to compute a 90% Highest Posterior Density interval for G. Compare the 
## two intervals.

GSorted=sort(G)[(0.05*length(G)+1):(0.95*length(G))]
# 90 % credible interval for G through the simulated draws
G_CredInterval=c(min(GSorted),max(GSorted))
print(G_CredInterval)
plot(density(G), main="Density function of simulated values of the Gini coefficient with credible intervals")
abline(v = G_CredInterval[1], col="blue")
abline(v = G_CredInterval[2], col="blue")

GDensity=density(G)
GDensity.df=data.frame(x=GDensity$x, y=GDensity$y)
GDensity.df=GDensity.df[order(-GDensity.df[,2]),]
index=dim(GDensity.df)[1]
GDensity.df$y=cumsum(GDensity.df$y)/sum(GDensity.df$y)
GDensity_CredInterval_Vals=GDensity.df[GDensity.df$y<0.90,]
GDensity_CredInterval=c(min(GDensity_CredInterval_Vals$x), max(GDensity_CredInterval_Vals$x))
print(GDensity_CredInterval)
abline(v = GDensity_CredInterval[1], col="red")
abline(v = GDensity_CredInterval[2], col="red")
title(sub="Blue = Simulated credible interval, Red = Kernel estimated credible interval")

## As seen in the plot the credible intervals are quite similar with small deviations. 


