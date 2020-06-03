## a) Assume following joint prior p(alpha, beta) 1/(alpha*beta)^2

# Reading the data from file
load(file = 'weibull.RData')
library(mvtnorm)

data=weibull

logPrior=function(alphaBeta) {
  return(-2*(log(alphaBeta[1])+log(alphaBeta[2])))
}

logPosterior = function(alphaBeta, data) {
  alpha=alphaBeta[1]
  beta=alphaBeta[2]
  logPrior=logPrior(alphaBeta)
  logLike=sum(dweibull(data, alpha, beta, log=TRUE))
  logPost=logLike+logPrior
  if (abs(logPost) == Inf || is.na(logPost)) logPost = -20000;
  return(logPost)
}

# Defining initial values to be passed on to the optimizer
set.seed(12345)
initVals = rnorm(2)

optimResult = optim(initVals, logPosterior, data=data, method=c("L-BFGS-B"),
                    lower=c(0.0001,0.0001), upper=c(Inf, Inf), control=list(fnscale=-1), hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = -solve(optimResult$hessian)
names(postMode) = c("alpha", "beta")
approx_PostStd = sqrt(diag(postCov))
names(approx_PostStd) = c("alpha", "beta")
colnames(postCov) = c("alpha", "beta")
rownames(postCov) = c("alpha", "beta")
print("The posterior mode is:")
print(postMode)
print("The covariance matrix is:")
print(postCov)

## b) Simulate from the actual posterior using the metropolis algorithm. Denote theta=t(alpha,beta) and use
## proposal density the multivariate normal density (random walk metropolis). Use alpha=1 and beta=1 as starting
## values, 500 iterations burnin and thereafter 2000 samples from the posterior. Run the algorithm for c: 0.1,
## 4 and 100 and use the draws from the best choice of c. Motivate your choice. Compute posterior mean and variance
## for the two params based on your samples. Proposal distrib can be truncated to avoid proposals below zero.

nIter=2000
burnIn=500
alphaGamma=c(1,1)
c=c(0.1, 4, 100)

# Defining function for sampling through metropolishastings
RVMSampler = function(previousVal, postCov, c, myFunction, ...) {
  proposalVal=rmvnorm(1, mean=previousVal, sigma=c*postCov)
  proposalVal[proposalVal<=0]=1e-6
  alpha=min(1, exp(myFunction(proposalVal,...)-myFunction(previousVal, ...)))
  u=runif(1)
  if(u < alpha) {
    return(list(theta=proposalVal, acceptProb=alpha))
  } else {
    return(list(theta=previousVal, acceptProb=alpha))
  }
}

post_matrix = matrix(0, nIter+burnIn, 2*length(c))
colnames=c()
for (i in 1:length(c)) {
  colnames=c(colnames, c(paste("Alpha with c=", c[i], sep=""), paste("Beta with c=", c[i], sep="")))
}
# Setting initial values of beta to same initVals as in the optimizer (taken randomly from normal distrib)
post_matrix[1,]=alphaGamma
accProb=matrix(0, nIter+burnIn, length(c))
colnames(accProb)=c("c=0.1", "c=4", "c=100")
set.seed(12345)

for (j in 1:length(c)) {
  for(i in 1:(nIter+burnIn)) {
    if(i<(nIter+burnIn)) {
      draw=RVMSampler(post_matrix[i,(2*j-1):(2*j)], postCov, c[j], logPosterior, data)
      post_matrix[i+1,(2*j-1):(2*j)]=draw$theta
      accProb[i+1,j]=draw$acceptProb
    }
  }
}
accProb_final=accProb[-(1:burnIn),]
accProb_mean=colMeans(accProb)
accProb_mean
colnames(post_matrix)=colnames

## Should choose c=4 since it yields an acceptance probability of around 0.3 which is a preferable value when
## applying metropolis hastings algorithm. 

post_matrix_final=post_matrix[-(1:burnIn), 3:4]
rownames=seq(501,2500)
rownames(post_matrix_final)=rownames
postMean_alpha=mean(post_matrix_final[,1])
postMean_beta=mean(post_matrix_final[,2])
postVar_alpha=var(post_matrix_final[,1])
postVar_beta=var(post_matrix_final[,2])
postMean_alpha
postMean_beta
postVar_alpha
postVar_beta
