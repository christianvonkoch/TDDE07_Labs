## Assignment 2:  Consider the logistic regression Pr(y=1 given x)=exp(xT*Beta)/(1+exp(xT*Beta)), where y is the
## binary variable with y = 1 if the woman works and y = 0 if she does not. x is a 8-dimensional vector containing
## the eight features (including a one for the constant term that models the intercept). The goal is to approximate
## the posterior distribution of the 8-dim parameter vector beta with multivariate normal distribution.
## Beta given y and x ~ N(~Beta, JY(~Beta)^(-1)), where ~Beta is the posterior mode and J(~Beta) is the second
## derivative, the observed Hessian evaluated at the posterior mode. It is actually not hard to compute this
## derivative by hand, but don't worry, we will let the computer do it numerically for you. Now, both ~Beta and
## J(~Beta) are computed by the optim function in R. I want you to implement an own version of my example code at
## the website. You can use my code as a template, but I want you to write your own file so that you understand
## every line of your code. Don't just copy my code. Use the prior Beta ~ N(0, thao^2*I) with thao=10. Your report
## should include your code as well as numerical values for ~Beta and JY(~Beta)^(-1) for the WomenWork data. Compute
## an approximate 95% credible interval for the variable NSmallChild. Would you say that this feature is an
## important determinant of the probability that a women works? [Hint: To verify that your results are reasonable,
## you can compare to you get by estimating the parameters using maximum likelihood. 
## glmmodel = glm(Work~0+., data=WomenWork, family=binomial)

# Use of libraries
library(mvtnorm)

# Read data
WomenWork = read.table("WomenWork.dat", header=TRUE)

# User input
tau = 10

# Defining vectors X and Y
X = as.matrix(WomenWork[,2:ncol(WomenWork)])
Y = WomenWork[,1]
nFeatures = dim(X)[2]
covNames=names(WomenWork[,2:ncol(WomenWork)])

# Constructing prior
mu_prior = rep(0,nFeatures)
sigma_prior = tau^2*diag(nFeatures) 

# Defining function for returning the log posterior

logPostLogistic = function(beta, Y, X, mu, sigma) {
  nFeat = length(beta)
  XBeta=X%*%beta
  # Defining loglikelihood
  logLike = sum(Y*XBeta-log(1+exp(XBeta)))
  # Defining prior
  prior = dmvnorm(beta, mean=mu, sigma=sigma, log=TRUE)
  # Adding loglikelihood and logprior together. Since it is log both of them are added instead of multiplied
  return(logLike + prior)
}

# Defining initial values to be passed on to the optimizer
set.seed(12345)
initVals = rnorm(dim(X)[2])

# Finding the optimized betavector
optimResult = optim(initVals, logPostLogistic, Y=Y, X=X, mu=mu_prior, sigma=sigma_prior, method=c("BFGS"),
                    control=list(fnscale=-1), hessian=TRUE)

# Defining the values of interest
postMode = optimResult$par
postCov = -solve(optimResult$hessian)
names(postMode) = covNames
approx_PostStd = sqrt(diag(postCov))
names(approx_PostStd) = covNames
print("The posterior mode is:")
print(postMode)
print("The approximated standard deviations are:")
print(approx_PostStd)

# Compute marginal distribution for nSmallChild
NSmallChild_mode = as.numeric(postMode["NSmallChild"])
NSmallChild_std = as.numeric(approx_PostStd["NSmallChild"])
credInterval_NSmallChild = qnorm(p=c(0.025, 0.975), mean=NSmallChild_mode, sd=NSmallChild_std)
print(paste("The lower bound of the 95 % credible interval for the feature NSmallChild is",
            round(credInterval_NSmallChild[1], 6), "and the upper bound is", 
            round(credInterval_NSmallChild[2], 6)))

# Control that the calculations have been made correctly
glmModel = glm(Work ~ 0+., data=WomenWork, family=binomial)
print(glmModel$coefficients)
print(postMode)

## Since the values for the credible interval for NSmallChild are quite large in the negative direction it is
## reasonable to conclude that the feature NSmallChild affects the response variable farily much towards the 
## response 0 which means that the woman doesn't work. This seems like a reasonable conclusion in terms of how
## it is in reality as well. When checking if the results are reasonable, a comparison was made with an 
## estimation using the maximum likelihood method. The results was very similar which strongly suggests that 
## the results obtained from the code are reasonable. 

## b) Write a function that simulates from the predictive distribution of the response variable in a logistic
## regression. Use your normal approximation from 2(a). Use that function to simulate and plot the predictive
## distribution for the Work variable for a 40 year old woman, with two children (3 and 9 years old), 8 years
## of education, 10 years of experience. and a husband with an income of 10. [Hints: The R package mvtnorm will
## again be handy. Remember my discussion on how Bayesian prediction can be done by simulation.]

sigmoid = function(value) {
  return (exp(value)/(1+exp(value)))
}

makePredLogReg = function(data, mean, sigma, nDraws) {
  betaPred = rmvnorm(nDraws, mean=mean, sigma=sigma)
  linearPred = betaPred %*% data
  logPred = sigmoid(linearPred)
  return(logPred)
}

nDraws=10000
woman=c(1, 10, 8, 10, (10/10)^2, 40, 1, 1)
set.seed(12345)
womanWorkPred=makePredLogReg(woman, postMode, postCov, nDraws)
logistic_distrib=c()
for (i in womanWorkPred) {
  logistic_distrib=c(logistic_distrib, rbinom(1, 1, i))
} 
barplot(table(logistic_distrib), main="Histogram of the predicted probabilities")

## As seen in the plots the calculated probabilities of the woman in question working is fairly low. The highest
## density is seen in the range between 0.2 and 0.3 approximately. This also makes sense if applied to a real
## situation. A woman with a small child is likely to stay at home with the child, i.e. not working. If the 
## classification of the response variable results in "working" if the predicted probabilities is above 0.5 and
## "not working" otherwise, it is clear from the distribution that the classification of a woman working, with
## the parameters inputted, is very unlikely to happen. 

## c) Now, consider 10 women which all have the same features as the woman in 2(b). Rewrite your function and
## plot the predictive distribution for the number of women, out of these 10, that are working.
## [Hint: Which distribution can be described as a sum of Bernoulli random variables?]

makePredLogRegMultiple = function(data, mean, sigma, nDraws, n) {
  multiplePred=c()
  for (i in 1:nDraws) {
    betaDraw = makePredLogReg(data, mean, sigma, 1)
    multiplePred=c(multiplePred, rbinom(1, n, betaDraw))
  }
  barplot(table(multiplePred), main=paste("Distribution for prediction made on", n, "women"), 
       xlab="No. of women")
}

makePredLogRegMultiple(woman, postMode, postCov, 10000, 10)

## As seen in the histogram the binomial case resembles the density of predicted probabilities with the
## highest density found at 2 women. This result seems reasonable since when the number of draws taken from 
## the binomial distribution goes towards infinity the shape of the corresponding distribution will resemble
## the shape of the distribution for the probability p in the Bernoulli case, more and more. 