## b) Simulate predictive draw of max no. of years until next earthquake occurs, 95 % prob. alpha=1, beta=1. 

alpha=1
beta=1
xObs=c(35, 14, 4, 10, 2)
n=length(xObs)i
nIter=5000
predDistrib=rep(0,nIter)
for(i in 1:nIter) {
  posteriorDraw=rbeta(1,alpha+n, beta+sum(xObs))
  predDistrib[i]=rgeom(1,posteriorDraw)
}
predDistrib_maxYear=quantile(predDistrib, probs=0.95)
predDistrib_maxYear