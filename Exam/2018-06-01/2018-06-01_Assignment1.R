## a) Draw 1000 samples from prior (Gamma) and 1000 samples from posterior (Gamma). Plot prior and posterior using
## both samples and their analytical expressions. 

n=50
x_mean=10
beta=2
nDraws=1000

# We know that posterior distribution is the Gamma(alpha+sum(data), beta+n). Mean for Gamma distrib is alpha/beta.
## If beta=2 then beta+n for posterior is 52. alpha/2=(alpha+500)/52 which yields 50*alpha=1000 and alpha=20
## Check: 20/2=10, (20+500)/52=10 OK!

alpha=20 # According to motivation above
post_draws=rgamma(nDraws, alpha+n*x_mean, beta+n)
prior_draws=rgamma(nDraws, alpha, beta)
gridWidth=0.01
muGrid_post=seq(7,12, gridWidth) # Range taken with inspiration from histogram
muGrid_prior=seq(4,20,gridWidth)
par(mfrow=c(2,1))
hist(post_draws, breaks=50, main="Posterior", xlab=expression(mu),
     freq=FALSE)
lines(muGrid_post, dgamma(muGrid_post, alpha+n*x_mean, beta+n), lwd=2, xlab=expression(mu))
hist(prior_draws, breaks=50, main="Prior", xlab=expression(mu),
     freq=FALSE)
lines(muGrid_prior, dgamma(muGrid_prior, alpha, beta), lwd=2, xlab=expression(mu))

## As seen in the plots the distributions resemble each other. 

## b) Simulate 1000 draws from predictive distribution of new observation and plot distribution.

par(mfrow=c(1,1))
x_pred=rpois(1000, lambda=post_draws)
hist(x_pred, breaks=50, main="Histogram, approximated posterior predictive distribution", xlab=expression(mu),
     freq=FALSE)

## c) Prob that x51=10 based on posterior predictive distribution

sum(x_pred==10)/nDraws

