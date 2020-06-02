## a) Model posterior data with prior Gamma and likelihood Poisson, plot the posterior

# We know that posterior mapping with gamma prior and poisson likelihood is gamma distributed

sumBids=sum(bids)
n=length(bids)
alpha=1
beta=1
posterior_theta=dgamma(seq(3,4,0.001), alpha+sumBids, beta+n)
plot(seq(3,4,0.001), posterior_theta, type="l", lwd=2)

# b) Investigate through graphical methods if Poisson model describes data well

xGrid=seq(min(bids), max(bids))
data_norm=bidsCounts/sum(bidsCounts)
nDraws=5000
thetaDraws=rgamma(nDraws, alpha+sumBids, beta+n)
poissonDensity=rep(0, length(xGrid))
for (i in thetaDraws) {
  poissonDensity=poissonDensity+dpois(xGrid, lambda=i)
}

avgPoissonDensity=poissonDensity/nDraws
plot(xGrid, data_norm, xlab="No. of bids", ylab="Density", main="Fitted models", type="o", cex=0.8,
     ylim=c(0,0.25), lwd=2)
lines(xGrid, avgPoissonDensity, col="red", lwd=2, type="o")
legend(x=7, y=0.2, col=c("black", "red"), legend=c("Data", "Poisson mean density"), lty=c(1,1), 
       lwd=c(2,2), pch=c("o", "o"))

## Terrible fit which the plot shows

## c) Use GibbsMixPois.R. Esimate the mixture of Poissons both with K=2 and K=3. nIter=5000.

GibbsMixPois <- function(x, nComp, alpha, alphaGamma, betaGamma, xGrid, nIter){
  
  # Gibbs sampling for a mixture of Poissons
  # Author: Mattias Villani, IDA, Linkoping University. http://mattiasvillani.com
  #
  # INPUTS:
  #   x - vector with data observations (counts)
  #   nComp - Number of mixture components to be fitted
  #   alpha - The prior on the mixture component weights is w ~ Dir(alpha, alpha,..., alpha) 
  #   alphaGamma and betaGamma - 
  #              The prior on the mean (theta) of the Poisson mixture components is 
  #              theta ~ Gamma(alphaGamma, betaGamma) [rate parametrization of the Gamma dist]
  #   xGrid - the grid of data values over which the mixture is evaluated and plotted
  #   nIter - Number of Gibbs iterations
  #
  # OUTPUTS:
  #   results$wSample     - Gibbs sample of mixture component weights. nIter-by-nComp matrix
  #   results$thetaSample - Gibbs sample of mixture component means.   nIter-by-nComp matrix
  #   results$mixDensMean - Posterior mean of the estimated mixture density over xGrid.
  
  
  ####### Defining a function that simulates from a Dirichlet distribution
  rDirichlet <- function(param){
    nCat <- length(param)
    thetaDraws <- matrix(NA,nCat,1)
    for (j in 1:nCat){
      thetaDraws[j] <- rgamma(1,param[j],1)
    }
    thetaDraws = thetaDraws/sum(thetaDraws) # Diving every column of ThetaDraws by the sum of the elements in that column.
    return(thetaDraws)
  }
  
  # Simple function that converts between two different representations of the mixture allocation
  S2alloc <- function(S){
    n <- dim(S)[1]
    alloc <- rep(0,n)
    for (i in 1:n){
      alloc[i] <- which(S[i,] == 1)
    }
    return(alloc)
  }
  
  # Initial values for the Gibbs sampling
  nObs <- length(x)
  S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
  theta <- rep(mean(x), nComp) # Each component is initialized at the mean of the data
  
  # Setting up the grid where the mixture density is evaluated.
  mixDensMean <- rep(0,length(xGrid))
  effIterCount <- 0
  
  # Setting up matrices to store the draws
  wSample <- matrix(0, nIter, nComp)
  thetaSample <- matrix(0, nIter, nComp)
  probObsInComp <- rep(NA, nComp)
  
  # Setting up the priors - the same prior for all components
  alpha <- rep(alpha, nComp) 
  alphaGamma <- rep(alphaGamma, nComp) 
  betaGamma <- rep(betaGamma, nComp) 
  
  # HERE STARTS THE ACTUAL GIBBS SAMPLING
  
  for (k in 1:nIter){
    message(paste('Iteration number:',k))
    alloc <- S2alloc(S) # Function that converts between different representations of the group allocations
    nAlloc <- colSums(S)
    
    # Step 1 - Update components probabilities
    w <- rDirichlet(alpha + nAlloc)
    wSample[k,] <- w
    
    # Step 2 - Update theta's in Poisson components
    for (j in 1:nComp){
      theta[j] <- rgamma(1, shape = alphaGamma + sum(x[alloc == j]), rate = betaGamma + nAlloc[j])
    }
    thetaSample[k,] <- theta
    
    # Step 3 - Update allocation
    for (i in 1:nObs){
      for (j in 1:nComp){
        probObsInComp[j] <- w[j]*dpois(x[i], lambda = theta[j])
      }
      S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
    }
    
    # Computing the mixture density at the current parameters, and averaging that over draws.
    effIterCount <- effIterCount + 1
    mixDens <- rep(0,length(xGrid))
    for (j in 1:nComp){
      compDens <- dpois(xGrid, lambda = theta[j])
      mixDens <- mixDens + w[j]*compDens
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
  }
  return(results = list(wSample = wSample, thetaSample = thetaSample, mixDensMean = mixDensMean))
}

result_comp2=GibbsMixPois(bids, nComp=2, alpha=1, alphaGamma = alpha, betaGamma = beta, 
                    xGrid=xGrid, nIter=500)
result_comp3=GibbsMixPois(bids, nComp=3, alpha=1, alphaGamma = alpha, betaGamma = beta, 
                          xGrid=xGrid, nIter=500)

## c) Use graphical methods to investigate if mixture of poissons fits data well. Is K=2 enough or should we
## use K=3?

plot(xGrid, data_norm, xlab="No. of bids", ylab="Density", main="Fitted models", type="o",
     ylim=c(0,0.25), lwd=2)
lines(xGrid, result_comp2$mixDensMean, col="red", lwd=2, type="o")
lines(xGrid, result_comp3$mixDensMean, col="gray", lwd=2, type="o")
legend(x=7, y=0.2, col=c("black", "red", "gray"), 
       legend=c("Data", "Mixture density with 2 components", "Mixture density with 3 components"), 
       lty=c(1,1,1), lwd=c(2,2, 2), pch=c("o", "o", "o"), cex=1)

## Good enough with 2 components in the mixture density



