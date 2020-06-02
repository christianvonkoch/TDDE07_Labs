## a) Simulate 1000 draws from predictive distrib of the maximal weight on a given future day, model: y=10*a where 
## y is the weight and a is the build cost. y~N(theta, sigma2). Noninformative prior assumed.

y=c(191, 196, 197, 189)
sigma2=10^2
# Noninformative prior assumed to be constant

yPred_post=rnorm(1000, mean=mean(y), sd=sqrt(sigma2*(1+1/length(y))))

## b) Use simulation to approximate the predictive probability that weight higher than 230

pred_max365=rep(0,1000)
for (i in 1:1000) {
  pred_max365[i]=max(rnorm(365, mean=mean(y), sd=sqrt(sigma2*(1+1/length(y)))))
}
prob_yPred365=sum(pred_max365>230)/1000
print(prob_yPred365)

## The probability is 0.157

## c) The loss function is linear

expectedLoss = function(a, maxWeight) {
  probCollapse=sum(maxWeight>10*a)/1000
  return(a*(1-probCollapse)+probCollapse*(a+100))
}

a=seq(20,30,0.01)
plot(a, sapply(a, expectedLoss, maxWeight=pred_max365), type="l", lwd=2, xlab="a", ylab="EL",
     main="Loss function")
aOpt=a[which(min(sapply(a, expectedLoss, maxWeight=pred_max365)))]
print(aOpt)

## The answer is 23.89


