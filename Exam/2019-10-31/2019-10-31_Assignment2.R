## a) Consider poisson likelihood model. Use conjugate prior and plot posterior in given interval. 
## Compute posterior probability that theta is smaller than 21. 

# Calculations show that alpha=20, beta=1

data=Traffic$y
alpha=20
beta=1
n=length(data)

# We know that Poisson with gamma prior is gamma distributed with alphaNew=alpha+sum(data), betaNew=beta+n

grid=seq(18,24,0.01)
post_distrib=dgamma(grid, shape=alpha+sum(data), rate=beta+n)
plot(grid, post_distrib, type="l", lwd=2, main="Posterior distrib. of theta", xlab=expression(theta))
post_prob=pgamma(21, shape=alpha+sum(data), rate=beta+n)

## Answer: Probability is 0.0557

## b) Two independent poisson models. 

data_model1=Traffic[which(Traffic[,3]=="yes"),]$y
data_model2=Traffic[which(Traffic[,3]=="no"),]$y

alpha_1=20+sum(data_model1)
alpha_2=20+sum(data_model2)
beta_1=1+length(data_model1)
beta_2=1+length(data_model2)
post_distrib_1=rgamma(5000, shape=alpha_1, rate=beta_1)
post_distrib_2=rgamma(5000, shape=alpha_2, rate=beta_2)
hist(post_distrib_1, breaks=50)
hist(post_distrib_2, breaks=50)
post_diff=post_distrib_2-post_distrib_1
hist(post_diff,
     main="Posterior distribution of difference between no speedlimit and speedlimit", xlab="No. of accidents")
quantile(post_diff, prob=c(0.025, 0.975))
mean(post_diff)

## We can see that the difference between the two distributions is larger than 0 with high probability. In this
## case we can say that the difference in traffic accidents between when no speed limit were applied and 
## when a speed limit were applied is between 2.82 and 5.53 approximately with 95 % posterior probability. 
## The conclusion from this is that yes, a speed limit leads to a lower amount of accidents.

## c) A politician claims that the experiment proves that introducing speed limit decreases the number
## of accidents by at least 15 %. 

mean(0.85*post_distrib_2>post_distrib_1)

## Likely that the decrease yields 15 % but 86 % probable and not 95 % probability which is commonly used
## in statistical experiments. 

