## a) Consider observations with values above 200. Remaining datapoints assumed to be indep. and follow a 
## truncated normal distribution with density specified. L=200 lower truncation point. Write a function in R
## that computes the (unnormalized) log posterior distribution of mu. Use function to plot the posterior distrib
## of mu for the observations greater than 200 in the data vector sulfur. For the plot, use a grid constructed
## in R with seq(100,400,1)

# Reading the data from file
load(file = 'sulfur.RData')

muGrid=seq(100,400,1)
sigma=100
data=sulfur[sulfur>200]

# Constant prior for mu is assumed

logPost = function(data, mu, sigma, L=200) {
  nominator=dnorm((data-mu)/sigma, mean=0, sd=1, log=TRUE)
  denominator=log(sigma)+log(1-pnorm((L-mu)/sigma))
  return(sum(nominator-denominator+0)) # Assumed constant prior which can be set to 1 which in log scale is 0
}

post_mu=exp(sapply(muGrid, logPost, data=data, sigma=sigma))
post_mu_norm=post_mu/sum(post_mu) # Since gridwidth is 1 we don't have to compensate for it
plot(muGrid, post_mu_norm, type="l", lwd=2, main="Posterior distribution of mu", xlab=expression(mu))

library(rstan)
T = length(sulfur)
T_cens = sum(sulfur <= 200)
censData <- list(T=T, T_cens = T_cens, x=sulfur, L=200)

# Model
censModel <- '
data {
  int<lower=0> T;       // Total number of time points
  int<lower=0> T_cens;  // Number of censored time points
  real x[T];            // Partly censored data
  real<upper=max(x)> L; // Lower truncation point
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<upper=L> x_cens[T_cens]; // Censored values
}

model {
  int t_cens = 0;
  for (t in 1:T){
    if (x[t] > L) 
      x[t] ~ normal(mu,sigma);
    else {
      t_cens += 1;
      x_cens[t_cens] ~ normal(mu,sigma);
    }
  }
}
'

## b) Now condiser all data points. Values below 200 being censored. 

fit=stan(model_code=censModel, data=censData)
print(fit)
post_draws=extract(fit)
grid=seq(1,4000,1)
plot(grid, post_draws$mu, type="l", main="Traceplot of mu", xlab=expression(mu), ylab="Value")
plot(grid, post_draws$sigma, type="l", main="Traceplot of sigma", xlab=expression(sigma), ylab="Value")
par(mfrow=c(4,2))
for (i in 1:8) {
  plot(grid, post_draws$x_cens[,i], type="l", main=paste("Traceplot of ", i, "th obs of obs below 200", sep=""),
       xlab=i, ylab="Value") 
}
par(mfrow=c(1,1))

plot(post_draws$mu, post_draws$sigma, type="p", col="grey", main="Joint posterior of mu and sigma",
     xlab=expression(mu), ylab=expression(sigma))

## c) Instead consider time series model. Assume that observations follow an independent normal distrib
## when conditioned on a latent AR(1) process z, but with values of xi below 200 being censored and set to 200.
## Modify the stan code in order to do inference for this model instead. Also put a normal prior on 
## mu~N(300,100^2) Plot the posterior of phi. Also produce a plot that contains both the data and the posterior 
## mean and 95 % credible intervals for the latent intensity z over time. 

StanModel_AR = '
data {
  int<lower=0> T;       // Total number of time points
  int<lower=0> T_cens;  // Number of censored time points
  real x[T];            // Partly censored data
  real<upper=max(x)> L; // Lower truncation point
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  real<upper=L> x_cens[T_cens]; // Censored values
  vector[T] z;
}

model {
  // Prior
  int t_cens = 0;
  phi ~ uniform(-1,1);
  mu ~ normal(300, 100);
  for (n in 2:T)
    z[n] ~ normal(mu + phi * (z[n-1]-mu), sigma);
    

  // Model/likelihood
  for (t in 1:T){
    if (x[t] > L) 
      x[t] ~ normal(z[t],20);
    else {
      t_cens += 1;
      x_cens[t_cens] ~ normal(z[t],20);
    }
  }
}

generated quantities {
  vector[T] post_mean;
  post_mean = z;
}
'
fitAR=stan(model_code=StanModel_AR, data=censData)
print(fitAR)
post_draws_AR=extract(fitAR)
postPhi=post_draws_AR$phi
postZ=post_draws_AR$post_mean
hist(postPhi, breaks=50, main="Approximated posterior density of phi", xlab=expression(phi), freq=FALSE)
grid=seq(1,31)
plot(grid, sulfur, col="blue", main="Emissions of sulfur dioxide", xlab="Day of month", ylab="mg/Nm^3",
     ylim=c(0,500))
postMean=rep(0,ncol(postZ))
credIntervals=matrix(0,ncol(postZ),2)
for (i in 1:ncol(postZ)) {
  postMean[i]=mean(postZ[,i])
  credIntervals[i,]=quantile(postZ[,i], probs=c(0.025, 0.975))
}
lines(grid, postMean, type="l", col="red", lwd=2)
lines(grid, credIntervals[,1], col="grey", lwd=1, lty=2)
lines(grid, credIntervals[,2], col="grey", lwd=1, lty=2)
legend("topleft", legend=c("Data", "Posterior mean", "95 % cred intervals"), lwd=c(NaN, 2, 1), lty=c(NaN,1,2),
       pch=c(1,NaN, NaN), col=c("blue", "red", "grey"))
