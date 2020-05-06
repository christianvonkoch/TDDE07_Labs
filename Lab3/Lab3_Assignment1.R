## Assignment 1: The data rainfall.dat consist of daily records, from the beginning of 1948 to the end of 1983, of
## precipitation (rain or snow in units of 1/100 inch, and records of zero precipitation are exluded) at Snoqualmie
## Falls  Washington. Analyze the data using the following two models.

## a) Assume the daily precipitation (y1,...,yn) are iid normally distributied, 
## y1,...,yn given mu and sigma^2 ~ N(mu,sigma^2) where both mu and sigma^2 are unknown. Let mu ~ N(mu0, tao0^2)
## independently of sigma^2 ~ Inv chisquare(v0, sigma0^2)
## i)  Implement (code!) a Gibbs sampler that simulates from the joint posterior p(mu, sigma^2 given y1,...,yn).
## The full conditional posteriors are given on the slides from Lecture 7. 

# Read data
Rainfall = read.table("rainfall.dat")

# Setup
mu0=mean(Rainfall[,1])
tao0sq=100
v0=1
sigma0sq=1
# Initial sigma value for Gibbs sampling
sigma=1
n=dim(Rainfall)[1]
nDraws=5000

calcTaoN = function(sigmasq,tao0sq,n){
  return(1/(n/sigmasq+1/tao0sq))
}

calcMuN = function(sigmasq, tao0sq, mu0, mean, n) {
  w=(n/sigmasq)/(n/sigmasq+1/tao0sq)
  return(w*mean+(1-w)*mu0)
}

