## Bayesian data analysis before an upcoming election. According to model, the posterior of the vote share theta
## that the party will get in the election is Beta(sqrt(c), 20) distributed, where c is the amount (in million SEK)
## that the party spends on the compaign. 
## a) Sample 10000 draws from posterior of thetafor the cases when c=4 and c=16. Based on the samples, plot the 
## posterior of theta for both cases.

postTheta_4=rbeta(10000, sqrt(4), 20)
postTheta_16=rbeta(10000, 4, 20)
hist(postTheta_4, breaks=50, main="Approximated posterior of theta", sub="c=4", xlab=expression(theta), freq=FALSE)
hist(postTheta_16, breaks=50, main="Approximated posterior of theta", sub="c=16", xlab=expression(theta), freq=FALSE)

## Compute probability that party gets at least 10 % of the votes for both cases. 

prob_4=sum(postTheta_4>0.1)/10000
prob_16=sum(postTheta_16>0.1)/10000
prob_4
prob_16

## To evaluate election results, utility function u(theta, c)=100+20log(theta)-c is used. How much money should
## they spend on the campaign? Do at least 10000 draws. Assume that maximum campaign budget is 20 million SEK.
## Consider values of c on the grid seq(4,20,0.5)

cGrid=seq(4,20,0.5)
postC=matrix(0,10000,length(cGrid))
set.seed(12345)
utility = function(c) {
  thetaDraws=rbeta(10000, sqrt(c), 20)
  return(100+20*log(thetaDraws)-c)
}
# Since quadratic loss function, choose posterior mean
postMean=rep(0,length(cGrid))
for(i in 1:length(cGrid)) {
  postC[,i]=utility(cGrid[i])
  postMean[i]=mean(postC[,i])
}

plot(cGrid, postMean, main="Posterior mean of expected utility", xlab="c", ylab="Utility", type="l", lwd=2)
optC=cGrid[which(postMean==max(postMean))]
abline(v=optC, col="blue", lwd=1, lty=2)
legend("topright", legend=c("Utility function", "Optimal c"), col=c("black", "blue"), lty=c(1,2), lwd=c(2,1))
optC

## Optimal decision is to invest 10.5 million SEK into the campaign

