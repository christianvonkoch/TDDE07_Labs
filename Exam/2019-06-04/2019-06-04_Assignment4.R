## a) Exp model with Gamma prior. Which prior is more informative, Gamma(2,1) or Gamma(10,10)?

thetaGrid=seq(0,20,0.01)
par(mfrow=c(2,1))
plot(thetaGrid, dgamma(thetaGrid, 2, 1), type="l", lwd=2, main="Prior distrib of Gamma(2,1)", xlab=expression(theta),
     ylab="Density")
plot(thetaGrid, dgamma(thetaGrid, 10, 10), type="l", lwd=2, main="Prior distrib of Gamma(10,10)", xlab=expression(theta),
     ylab="Density")
par(mfrow=c(1,1))

## As seen in the graph the prior for theta in the plot below is more informative since it has a tighter peak than
## the graph above it. This means that it is more probable that theta is a specific value whereas in the above plot
## the probability is more spread over a larger interval of possible theta values. 

## b) Compute marginal likelihood for the two models. 

data=cellphones

margLikelihood = function(data, alpha, beta) {
  n=length(data)
  nominator=beta^alpha*gamma(alpha+n)
  denominator=gamma(alpha)*(beta+sum(data))^(alpha+n)
  return(nominator/denominator)
}

margLike1=margLikelihood(data, 2, 1)
margLike2=margLikelihood(data, 10, 10)
postModel1=margLike1*0.5
postModel2=margLike2*0.5
postModel1_norm=postModel1/sum(c(postModel1, postModel2))
postModel2_norm=postModel2/sum(c(postModel1, postModel2))

## Model 1 is more probable and should be chosen!

## c) Compute 90 % posterior predictive interval of x~ given the cellphones dataset. 

predLikelihood = function(data, alpha, beta, xTilde) {
  n=length(data)
  nominator=(beta+sum(data))^(alpha+n)
  denominator=(beta+xTilde+sum(data))^(alpha+n+1)
  return(nominator*(alpha+n)/denominator)
}

xTilde=seq(0,20, 0.01)
xTilde1=sapply(xTilde, predLikelihood, data=data, alpha=2, beta=1)
xTilde2=sapply(xTilde, predLikelihood, data=data, alpha=10, beta=10)
post_xTilde=postModel1_norm*xTilde1+postModel2_norm*xTilde2  

ndraws = 1000000
xTildeDraws = rep(0,ndraws)
for(i in 1:ndraws){
  M = rbinom(1,1,postModel2_norm) + 1 # Simulate which model to use
  if(M==1){
    theta = rgamma(1,shape=2+length(x),rate=1+sum(x))
  } else {
    theta = rgamma(1,shape=10+length(x),rate=10+sum(x))
  }
  xTildeDraws[i] = rexp(1,theta)
}

print(quantile(xTildeDraws,probs = c(.05,.95)))
hist(xTildeDraws, freq=FALSE, breaks=1000)
lines(xTilde, post_xTilde, type="l", lwd=2, col="red")
plot(xTilde, post_xTilde, type="l", lwd=2)
