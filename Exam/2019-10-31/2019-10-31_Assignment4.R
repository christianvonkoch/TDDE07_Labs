## a) Use supplied stan model to do Bayesian inference. Draw 2000 posterior samples and use 500 for burnin. 
## Produce figure with scatter plot, overlay curve for mean of posterior predictive distrib, in range [0,25]. 
## Also overlay curves 90 % equal tail interval for same posterior predictive distrib given values of x in range [0,25]

# Load data
cars = cars

library(rstan)
LinRegModel <- '
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma2;
}
model {
  sigma2 ~ scaled_inv_chi_square(5,10);
  for (n in 1:N)
    y[n] ~ normal(alpha + beta * x[n], sqrt(sigma2));
}
'
x=cars$speed
y=cars$dist
nIter=2000
burnIn=500
N=dim(cars)[1]
data=list(N=N,x=x,y=y)
fit=stan(model_code=LinRegModel, data=data, iter=nIter, warmup = 500, chains=1)
print(fit)
postDraws=extract(fit)
alpha_draws=postDraws$alpha
beta_draws=postDraws$beta
sigma_draws=postDraws$sigma2
xGrid=seq(0,25)
n=length(alpha_draws)
mean_credInt=matrix(0,length(xGrid),3)
count=1
for (i in 1:length(xGrid)) {
  ysim=rep(0,length(nIter-burnIn))
  ysim=alpha_draws+beta_draws*xGrid[i]+rnorm(nIter-burnin, mean=0, sd=sqrt(sigma_draws))
  mean_credInt[count,1]=mean(ysim)
  mean_credInt[count,-1]=quantile(ysim, probs = c(0.05, 0.95))
  count=count+1
}

plot(x,y,xlab="Speed", ylab="Distance", col="blue", main="Plot for model with constant sigma prior")
lines(xGrid, mean_credInt[,1], lwd=2, col="red")
lines(xGrid, mean_credInt[,2], lwd=1, lty=2)
lines(xGrid, mean_credInt[,3], lwd=1, lty=2)
legend("topleft", legend=c("Data", "Posterior mean", "90 % cred interval"), col=c("blue", "red", "grey"), 
       pch=c(1, NaN, NaN), lty=c(NaN, 1, 2), lwd=c(NaN, 2, 1))

## b) Compute 95 % equal tail credible interval for alpha. Give real-world interpret of the interval. 

quantile(alpha_draws, probs=c(0.025, 0.975))

## The interpretation of the credible interval for alpha is that if the car has no speed it travels a negative
## distance between -31 and 4.25 approximately with 95 % posterior probability. This is not realistic. To prevent this
## a prior can be set to alpha with a mean around zero which however would make the linear prediction worse. 
## One can also use the log Normal distribution for y to force it to have a value above zero. 

## c) Reproduce results in b) with heteroscadastic variance. 

LinRegModel_hetero <- '
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma2[N];
  real gamma;
  real phi;
}
model {
  for (n in 1:N)
    sigma2[n] ~ scaled_inv_chi_square(5,exp(gamma+phi*x[n]));
  for (n in 1:N)
    y[n] ~ normal(alpha + beta * x[n], sqrt(sigma2[n]));
}
'
data=list(N=N,x=x,y=y)
fit2=stan(model_code=LinRegModel_hetero, data=data, iter=nIter, warmup = 500, chains=1)
print(fit2)
postDraws2=extract(fit2)
alpha_draws=postDraws2$alpha
beta_draws=postDraws2$beta
sigma_draws=postDraws2$sigma2
xGrid=seq(0,25)
n=length(alpha_draws)
mean_credInt=matrix(0,length(xGrid),3)
count=1
for (i in 1:length(xGrid)) {
  rinv=rchisq(nIter-burnIn, 5)
  sigma_draw=5*exp(postDraws2$gamma + xgrid[i] * postDraws2$phi)^2/rinv
  ysim=rep(0,length(nIter-burnIn))
  ysim=alpha_draws+beta_draws*xGrid[i]+rnorm(nIter-burnin, mean=0, sd=sqrt(sigma_draw))
  mean_credInt[count,1]=mean(ysim)
  mean_credInt[count,-1]=quantile(ysim, probs = c(0.05, 0.95))
  count=count+1
}

plot(x,y,xlab="Speed", ylab="Distance", col="blue", main="Plot of model with heteroscadastic sigma prior")
lines(xGrid, mean_credInt[,1], lwd=2, col="red")
lines(xGrid, mean_credInt[,2], lwd=1, lty=2)
lines(xGrid, mean_credInt[,3], lwd=1, lty=2)
legend("topleft", legend=c("Data", "Posterior mean", "90 % cred interval"), col=c("blue", "red", "grey"), 
       pch=c(1, NaN, NaN), lty=c(NaN, 1, 2), lwd=c(NaN, 2, 1))

## The new model seems to capture the data better than the old one. 


         