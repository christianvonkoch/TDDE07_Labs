## a) Plot posterior density of lognormal likelihood and normal prior

sigma2=0.04
dataLion =lions
library(rstan)

posteriorDens = function(data, mu, sigma2, mu0, sigma2_0) {
  likelihood=sum(dlnorm(data, meanlog=mu, sdlog=sqrt(sigma2), log=TRUE))
  prior=dnorm(mu, mean=mu0, sd=sqrt(sigma2_0), log=TRUE)
  return(likelihood+prior)
}

gridWidth=0.001
muGrid=seq(5,5.5,gridWidth)
postMu=exp(sapply(muGrid, posteriorDens, data=dataLion, sigma2=0.04, mu0=5, sigma2_0=1))
postMu_norm=1/gridWidth*postMu/sum(postMu)
plot(muGrid, postMu_norm, type="l", lwd=2, main="Posterior density of mu", xlab=expression(mu), ylab="Density")

## b) Now assume that also sigma2 is unknown and that sigma2 ~ scaledinvchisq(v0, sigma2_0) a priori independently
## from mu, with v0=5 and sigma2_0=0.04. Implement stan-code that produces at least 2000 samples from the posterior
## of mu and sigma2. Use 500 samples for burnin. Based on samples compute posterior mean and standard deviation
## of mu and sigma2 and plot the joint posterior of mu and sigma2. 

StanModel= '
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma2;
}
model {
  //Priors
  mu ~ normal(5,1);
  sigma2~scaled_inv_chi_square(5,0.2);
  
  //Likelihood
  for (n in 1:N) {
    y[n]~lognormal(mu, sqrt(sigma2));
  }
}
'
n=length(dataLion)
data=list(N=n, y=dataLion)
fit=stan(model_code=StanModel, data=data, warmup=500, iter=2500, chains=1)
print(fit)
post_draws=extract(fit)
mu_post=post_draws$mu
sigma2_post=post_draws$sigma2
mu_postMean=mean(mu_post)
sigma2_postMean=mean(sigma2_post)
mu_postMean
sigma2_postMean
plot(mu_post, sigma2_post, main="Simulated joint posterior of mu and sigma", xlab=expression(mu),
     ylab=expression(sigma))

## c) Compute an estimate of the average weight of male lions. Give an estimate and a 95 % credible interval of 
## the average weight of male lions based on the posterior computed in b. 



estimate=exp(mu_post+1/2*sigma2_post)
mean(estimate)
hist(estimate, breaks=100, main="Predicitve distribution of average weight of lions")
quantile(estimate, probs=c(0.025, 0.975))

## Important here to insert all the draws into the function for the mean and then take the average of that for
## an estimate. 