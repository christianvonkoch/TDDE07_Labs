## b) Compute the probability that party A gets a majority of the votes (more than 50 %) in the election. 
## Assume that everyone in the population is voting. 

# Assuming uniform Dirichlet prior, i.e. Dirichlet(1,1,1)
nComp=3
alphaPrior=rep(1,nComp)
data=c(184, 67, 149)
alphaPost=alphaPrior+data

simDirich = function(alpha) {
  nComp=length(alpha)
  gammavec=rep(0,nComp)
  for (i in 1:nComp) {
    x=rgamma(1, alpha[i], 1)
    gammavec[i]=x
  }
  z=gammavec/sum(gammavec)
  return(z)
}

nIter=10000
dirichlet=matrix(0, nIter, nComp)
for (i in 1:nIter) {
  dirichlet[i,]=simDirich(alphaPost)
}

prob=mean(dirichlet[,1]>0.5)
prob

## 0.039 % that party A gets a majority of the votes

## Calc probability that A becomes largest party

prob2=mean(dirichlet[,1]>dirichlet[,2] & dirichlet[,1]>dirichlet[,3])
prob2

## 0.971 % chance that party A becomes the biggest party.

## d) Assume that probability in c) was estimated through monte carlo simulation. Compute a 95 % confidence interval
## for estimated probility in c) with respect to the error from the Monte Carlo simulation. 

# If x modeled as Bin(1,p) where p is the probability obtained from c) and x=1 stands for {A becomes largest party}
# and x=0 stands for {A doesn't become largest party}. We get that the expected value of the summation of all
# monte carlo samples divided by the number of samples is p and the variance for same variable is p(1-p)/#Samples

p=prob2
N=10000
stdX=p*(1-p)/N
confInt=c(p-1.96*sqrt(1/N*p*(1-p)), p+1.96*sqrt(1/N*p*(1-p)))
confInt

## How many additional samples would be needed to reduce the width of the interval by half

# To reduce width by half the difference between upper and lower bound of interval divided by 2 needs to be 
# reduced by half. This yields: p+1.96*1/sqrt(N)*sqrt(p*(1-p))-(p-1.96*1/sqrt(N)*sqrt(p*(1-p)))/2=
# 1.96*..., 1.96/2*1/sqrt(N)...=1.96*1/sqrt(N)... So we need to increase N to N~ to obtain the left hand side
# expression. We get that 1/sqrt(N~)=1/(2*sqrt(N)), N~=4*N

NTilde=4*10000
diffN=NTilde-N
diffN

## We need to increase the sample by 30000 to lower the interval by half. 