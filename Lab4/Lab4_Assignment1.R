## Assignment 1:
## a) Write a function in R that simulate data from the AR(1)-process: xt=mu+phi(x(t-1)-mu) + epsilon(t), 
## epsilon(t)~N(0,sigma^2), for given values of mu, phi, and sigma^2. Start the process at x1=mu and then simulate
## values for xt for t=2,3,...,T and return the vector x1:T containing all time points. Use mu=10, sigma^2=2 and
## T=200 and look at different realizations (simulation) of x1:T for values of phi between -1 and 1 (this is the
## interval of phi where the AR-process is stable). Include a plot of at least one realization in the report. What
## effect does the value of phi have on x1:t

#install.packages("rstan")
mu=10
sigma_sq=2
T=200
x_init=mu
phi_vector=seq(-0.9,0.9,0.1)
results_matrix=matrix(0,200,length(phi_vector))
results_matrix[1,]=x_init
counter=1
set.seed(12345)

AR_process_function=function(mu, sigma_sq, T, phi) {
  x_init=mu
  result=rep(0,T)
  result[1]=x_init
  for (i in 2:T) {
    epsilon=rnorm(1,0,sqrt(sigma_sq))
    result[i]=mu+phi*(result[i-1]-mu)+epsilon
  }
  print(result)
  return(result)
}

results_matrix=matrix(0,T,length(phi_vector))
counter=1
for (phi in phi_vector) {
  results_matrix[,counter]=AR_process_function(mu,sigma_sq,T,phi)
  counter=counter+1
}
iter=seq(1,200,1)
for (i in 1:length(phi_vector)) {
  plot(iter, results_matrix[,i], main="Plot of realization of AR-process", sub=paste("Phi =", phi_vector[i]),
       xlab="Iteration", ylab="Value", type="l", col="grey")
}

## With phi-values below zero the process will oscillate faster but with phi-values above zero the process will
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
par(mfrow = c(1,1))
plot(postDraws_x$mu[1000:2000], postDraws_x$phi[1000:2000],ylab="phi", xlab="mu", main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_x)
# Bivariate posterior plots
pairs(fit_x)

# Do traceplots of the first chain
par(mfrow = c(1,1))
plot(postDraws_y$mu[1000:2000],postDraws_y$phi[1000:2000],ylab="mu", xlab="mu",main="Traceplot")
# Do automatic traceplots of all chains
traceplot(fit_y)
# Bivariate posterior plots
pairs(fit_y)

## The posterior mean, number of effective samples as well as 95 % credible interval are shown above for both of the
## simulated AR(1)-processes. It is possible to estimate the true values of the parameters for the sample which
## used a phi=0.3 when obtaining the dataset used in the simulation. However, it is not as obvious to estimate 
## the parameters' true values for the second sample where phi=0.95 were used to obtain the dataset used in this
## particular simulation. The credible intervals for the parameters in this simulation are very wide and it is
## difficult to predict with certainty the true vale of the parameter. This might be due to the higher correlation
## between the lags caused by the higher value of phi. 

## The convergence of the samplers are different. For the first sample which used phi=0.3, the convergence is
## evident whilst for the second sample the posterior distribution is not obvious. This correlates with the fact
## the credible intervals for the parameters on the second sample were very wide. 


