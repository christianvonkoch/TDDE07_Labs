## Assignment 1:
## a) Write a function in R that simulate data from the AR(1)-process: xt=mu+phi(x(t-1)-mu) + epsilon(t), 
## epsilon(t)~N(0,sigma^2), for given values of mu, phi, and sigma^2. Start the process at x1=mu and then simulate
## values for xt for t=2,3,...,T and return the vector x1:T containing all time points. Use mu=10, sigma^2=2 and
## T=200 and look at different realizations (simulation) of x1:T for values of phi between -1 and 1 (this is the
## interval of phi where the AR-process is stable). Include a plot of at least one realization in the report. What
## effect does the value of phi have on x1:t

#install.packages("rstan")
# Setting up values
mu=10
sigma_sq=2
T=200
x_init=mu
phi_vector=seq(-0.9,0.9,0.1)
results_matrix=matrix(0,200,length(phi_vector))
results_matrix[1,]=x_init
counter=1
set.seed(12345)

# Defining function for the AR process
AR_process_function=function(mu, sigma_sq, T, phi) {
  x_init=mu
  result=rep(0,T)
  result[1]=x_init
  for (i in 2:T) {
    epsilon=rnorm(1,0,sqrt(sigma_sq))
    result[i]=mu+phi*(result[i-1]-mu)+epsilon
  }
  return(result)
}

results_matrix=matrix(0,T,length(phi_vector))
counter=1
# Exploring the AR process with different values of phi
for (phi in phi_vector) {
  results_matrix[,counter]=AR_process_function(mu,sigma_sq,T,phi)
  counter=counter+1
}
iter=seq(1,200,1)
counter=1
for (i in 1:length(phi_vector)) {
  if (counter %% 6 == 0) {
    plot(iter, results_matrix[,i], main="Plot of realization of AR-process", 
         sub=paste("Phi =", phi_vector[i]),
         xlab="Iteration", ylab="Value", type="l", col="grey")
  }
  counter=counter+1
}

## Conclusion: With phi-values below zero the process will oscillate faster but with phi-values above zero the process will
## be more correlated. The correlation between the different iterations increases as the phi-value becomes larger. 
## This causes the oscillation to slow down and the process to move more slowly. 

## b) Use your function from a) to simulate two AR(1)-processes, x1:T with phi=0.3 and y1:T with phi=0.95. Now,
## treat the values of mu, phi and sigma^2 as unknown and estimate them using MCMC. Implement Stan-code that
## samples from the posterior of the three parameters, using suitable non-informative priors of your choice. 
## [Hint: Look at the time-series models examples in the Stan user's guide/reference manual, and note the different
## parametizations used here.]
## i) Report the posterior mean, 95% credible intervals and the number of effective posterior samples for the
## three inferred parameters for each of the simulated AR(1)-process. Are you able to estimate the true values? 
## ii)  For each of the two data sets, evaluate the convergence of the samplers and plot the joint posterior of
## mu and phi. Comments?

library(rstan)

x=rep(0,T)
y=rep(0,T)
set.seed(12345)
x=AR_process_function(mu, sigma_sq, T, 0.3)
set.seed(12345)
y=AR_process_function(mu, sigma_sq, T, 0.95)

# Defining Stan model
StanModel= '
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
}
model {
  for (n in 2:N)
    y[n] ~ normal(mu + phi * (y[n-1]-mu), sigma);
}
'

data_x=list(N=T, y=x)
data_y=list(N=T, y=y)
fit_x=stan(model_code=StanModel, data=data_x)
fit_y=stan(model_code=StanModel, data=data_y)
postDraws_x <- extract(fit_x)
postDraws_y <- extract(fit_y)
print(fit_x)
print(fit_y)

# Do traceplots of the first chain
plot(postDraws_x$mu[1000:2000], postDraws_x$phi[1000:2000],ylab="phi", xlab="mu", main="Traceplot")

# Do traceplots of the first chain
plot(postDraws_y$mu[1000:2000],postDraws_y$phi[1000:2000],ylab="mu", xlab="mu",main="Traceplot")

## Conclusion: The posterior mean, number of effective samples as well as 95 % credible interval are shown above for both of the
## simulated AR(1)-processes. It is possible to estimate the true values of the parameters for the sample which
## used a phi=0.3 when obtaining the dataset used in the simulation. However, it is not as obvious to estimate 
## the parameters' true values for the second sample where phi=0.95 were used to obtain the dataset used in this
## particular simulation. The credible intervals for the parameters in this simulation are very wide and it is
## difficult to predict with certainty the true vale of the parameter. This might be due to the higher correlation
## between the lags caused by the higher value of phi. 

## Conclusion: The convergence of the samplers are different. For the first sample which used phi=0.3, the convergence is
## evident whilst for the second sample the posterior distribution is not obvious. This correlates with the fact
## the credible intervals for the parameters on the second sample were very wide. What we can see from the 
## posterior distribution obtained by the second sampler is that for lower values of phi the distribution centers
## around a value between 10 and 20. This is a behaviour similar to what is shown in the posterior for the first
## sampler, where phi was set to 0.3 initially, since this distribution was much tighter around the value of 10
## for mu.

## c) The data campy.dat contain the number of cases of campylobacter infections in the north of the province
## Quebec (Canada) in four week intervals from January 1990 to the end of October 2000. It has 13 observations per
## year and 140 observations in total. Assume that the number of infections ct at each time point follows an 
## independent Poisson distribution when conditioned on a latend AR(1)-process xt, that is
## ct given xt ~ Poisson(exp(xt)), where xt is an AR(1)-process as in a). Implement and estimate the model in Stan,
## using suitable priors of your choice. Produce a plot that contains both the data and the posterior mean and
## 95 % credible intervals for the latent intensity theta_t=exp(xt) over time. 
## [Hint: Should xt be seen as data or parameters]

campy=read.table("campy.dat", header=TRUE)
library(rstan)

# Defining stan model
StanModel_Pois = '
data {
  int<lower=0> T;
  int c[T];
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[T] x;
}

model {
  // Prior
  phi ~ uniform(-1,1);
  for (n in 2:T)
    x[n] ~ normal(mu + phi * (x[n-1]-mu), sigma);
    

  // Model/likelihood
  for (n in 1:T)
    c[n] ~ poisson(exp(x[n]));
    
}

generated quantities {
  vector[T] post_mean;
  post_mean = exp(x);
}
'

data=list(T=dim(campy)[1], c=campy$c)
fit_pois=stan(model_code=StanModel_Pois, data=data)
print(fit_pois)
pois_mean_list=fit_pois@.MISC$summary$msd
post_mean=pois_mean_list[grep("post_mean", rownames(pois_mean_list)),]

plot(campy$c, col="blue", ylab="No. of infected", xlab="Time")
points(post_mean[,1], col="black", type="l")

# Obtaining the quantiles
quantiles=fit_pois@.MISC$summary$quan
# Obtaining the quntiles for the post_mean variable
quantiles_post_mean=quantiles[grep("post_mean", rownames(quantiles)),]
cred_interval_post_mean=matrix(0,dim(quantiles_post_mean)[1], 2)
cred_interval_post_mean[,1]=quantiles_post_mean[,1]
cred_interval_post_mean[,2]=quantiles_post_mean[,ncol(quantiles_post_mean)]

lines(cred_interval_post_mean[,1], col="gray", lty=1)
lines(cred_interval_post_mean[,2], col="gray", lty=1)
title(main="Plot of data vs approximated posterior")
legend("topleft", box.lty= 1, pch=c(1,NaN,NaN), legend=c("Data", "Posterior mean", "95 % cred. interval"),
       col=c("blue", "black", "gray"), lwd=c(NaN,1,1), lty=c(NaN, 1, 1))

## Conclusion: As seen in the plot above the posterior mean follows the data accurately. Almost all of the datapoints are
## inside the credible intervals which aren't that wide which indicates that the approximated posterior
## resembles the reality shown by the data well. 

## d) Now, assume that we have a prior belief that the true underlying intensity theta_t varies more smoothly than
## the data suggests. Change the prior for sigma_sq so that it becomes informative about that the AR(1)-process 
## increments epsilon_t should be small. Re-estimate the model using Stan with the new prior and produce the same
## plot as in c). Has the posterior for theta_t changed?

# Defining Stan Model
StanModel_Pois_Prior = '
data {
  int<lower=0> T;
  int c[T];
}

parameters {
  real mu;
  real phi;
  real<lower=0> sigma;
  vector[T] x;
}

model {
  // Prior
  phi ~ uniform(-1,1);
  sigma ~ scaled_inv_chi_square(140, 0.15);
  for (n in 2:T)
    x[n] ~ normal(mu + phi * (x[n-1]-mu), sigma);
    

  // Model/likelihood
  for (n in 1:T)
    c[n] ~ poisson(exp(x[n]));
    
}

generated quantities {
  vector[T] post_mean;
  post_mean = exp(x);
}
'
fit_pois_prior=stan(model_code=StanModel_Pois_Prior, data=data)
print(fit_pois_prior)
pois_mean_list_prior=fit_pois_prior@.MISC$summary$msd
post_mean_prior=pois_mean_list_prior[grep("post_mean", rownames(pois_mean_list)),]

plot(campy$c, col="blue", ylab="No. of infected", xlab="Time")
points(post_mean_prior[,1], col="black", type="l")

# Obtaining the quantiles
quantiles_prior=fit_pois_prior@.MISC$summary$quan
# Obtaining the quntiles for the post_mean variable
quantiles_post_mean_prior=quantiles_prior[grep("post_mean", rownames(quantiles)),]
cred_interval_post_mean_prior=matrix(0,dim(quantiles_post_mean)[1], 2)
cred_interval_post_mean_prior[,1]=quantiles_post_mean_prior[,1]
cred_interval_post_mean_prior[,2]=quantiles_post_mean_prior[,ncol(quantiles_post_mean)]

lines(cred_interval_post_mean_prior[,1], col="gray", lty=1)
lines(cred_interval_post_mean_prior[,2], col="gray", lty=1)
title(main="Plot of data vs approximated posterior")
legend("topleft", box.lty= 1, pch=c(1,NaN,NaN), legend=c("Data", "Posterior mean", "95 % cred. interval"),
       col=c("blue", "black", "gray"), lwd=c(NaN,1,1), lty=c(NaN, 1, 1))

## Conclusion: Now when we have specified a small prior for sigma it is noteable in the new plot that the posterior mean
## varies less and moves more smoothly. The consequence of this is that more datapoints lie outside of the 
## credible interval which suggests that the approximated posterior does not resemble the reality described by
## the data as accurately as before. However, by doing this one can avoid overfitting when the model is applied
## to a new dataset. 




