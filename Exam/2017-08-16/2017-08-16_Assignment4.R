## a) Simulate 1000 draws from the posterior distrib of theta using conjugate prior for theta with mean 250
## and std = 50. Poisson likelihood. 

data=c(220,323,174,229)
alpha=25
beta=0.1
n=length(data)

logPriorGamma = function(theta, alpha, beta) {
  return(dgamma(theta, 50, beta, log=TRUE))
}

logLike = function(data, theta) {
  n=length(data)
  first=sum(data)*log(theta)
  second=theta*n
  third=0
  for (i in data) {
    for (j in 1:i) {
      third=third+log(j)
    }
  }
  return(first-second-third)
}

logPosterior = function(data, theta, alpha, beta) {
  prior=logPriorGamma(theta, alpha, beta)
  likelihood=logLike(data, theta)
  return(likelihood + prior)
}

# Conjugate prior for poisson is Gamma(alpha, beta), we know that posterior is Gamma(alpha + sum(data), beta+n)
post_draws=rgamma(1000, alpha+sum(data), beta+n)
hist(post_draws, main="Posterior distribution of theta", xlab=expression(theta))

## b) Simulate 1000 draws from the predictive distrib of next quarter's demand, X5, and plot the draws
## in histogram. 

q5=rpois(1000, post_draws)
hist(q5, breaks=50, main="Predictive distribution of quarter 5", xlab="Qty")
sum(q5<=200)/1000

## c) 


utility <- function(a,X5){
  util = rep(0,length(X5))
  util[X5<=a] = 10*X5[X5<=a]-(a-X5[X5<=a])
  util[X5>a] = 10*a-0.05*(X5[X5>a]-a)^2
  return(util)
}

mean(q5)
a=seq(136,336,1)
results = matrix(0,length(q5),length(a))
count=1
nameVec=rep(0,length(a))
for (i in a) {
  results[,count]=utility(i,q5)
  nameVec[count]=as.character(i)
  count=count+1
}
opt_vector=matrix(0,1,length(a))
for (i in 1:length(a)) {
  opt_vector[i]=mean(results[,i])
}
colnames(opt_vector)=nameVec
opt_decision=as.numeric(opt_vector[,which(opt_vector==max(opt_vector))])
names(opt_vector[,which(opt_vector==max(opt_vector))])
plot(a, opt_vector, type="l", lwd=1, col="red")
abline(v=as.numeric(names(opt_vector[,which(opt_vector==max(opt_vector))])), col="blue")

