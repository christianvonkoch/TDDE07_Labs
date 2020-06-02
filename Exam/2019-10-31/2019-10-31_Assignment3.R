## c) Make simulations of joins posterior of v and pi using Gibbs sampling.

x=20
lambda=10
alpha=2
beta=2
nIter=2000
burnIn=500

results=matrix(0,burnIn+nIter,2)
initVal=lambda # Since lambda=30
results[1,1]=initVal
results[1,2]=rnorm(1)
for (i in 1:(nIter+burnIn-1)) {
  z=rpois(1, lambda*(1-results[i,2]))
  results[i+1,1]=z+x
  results[i+1,2]=rbeta(1, alpha+x, beta+results[i+1,1]-x)
}

grid=seq(burnIn+1, nIter+burnIn)
barplot(table(results[(burnIn+1):(nIter+burnIn),1]), main="Marginal posterior of nu", xlab=expression(nu))
hist(results[(burnIn+1):(nIter+burnIn),2], breaks=50, main="Marginal posterior of pi", xlab=expression(pi))
plot(grid, results[(burnIn+1):(nIter+burnIn),2],type="l")
plot(grid, results[(burnIn+1):(nIter+burnIn),1], type="l")

## Convergence seems good since markov chain is exploring full posterior and have good mixing. 

