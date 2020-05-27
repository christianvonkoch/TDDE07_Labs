## Assignment 3: Bayesian inference for the concentration parameter in the von Mises distribution. This exercise is concerned
## with directional data. The point is to show you that the posterior distribution for somewhat weird models can be
## obtained by plotting it over a grid of values. The data points are observed wind directions at a given location on
## ten different days. The data are recorded in degrees: (40, 303, 326, 285, 296, 314, 20, 308, 299, 296) where North
## is located at zero degrees (see Figure 1 on the next page, where the angles are measured clockwise). To fit with 
## Wikipedias description of probability distributions for circular data we convert the data into radians -pi<=y<=pi.
## The 10 observations in radians are (-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02).
## Assume that these data points are independent observations following the von Mises distribution
## p(y given my,k) = exp(k*cos(y-u))/(2*pi*I0(k)), -pi<=y<=pi, where I0(k) is the modified Bessel function of the 
## first kind of order zero (see ?besselI in R). The parameter my (-pi<=my<=pi) is the mean direction and k>0 is
## called the concentration parameter. Large k gives a small variance around my, and vice versa. Assume that my is
## known to be 2.39. Let K ~ Exponential(Lambda=1) a priori, where lambda is the rate parameter of the exponential
## distribution (so that the mean is 1/lambda).

## a) Plot the posterior distribution of k for the wind direction data over a fine grid of k values. 

data_radian=c(-2.44,2.14,2.54,1.83,2.02,2.33,-2.79,2.23,2.07,2.02)
my=2.39
lambda=1

# Function for computing the vonMisesDistrib for a given dataset
vonMisesDistrib = function(kappa, data, my){
  likelihood=1
  for (i in data) {
    likelihood=likelihood*exp(kappa*cos(i-my))/(2*pi*besselI(kappa, 0))
  }
  return(likelihood)
}

# Function for computing the exponential distribution
exponDistrib = function(data, lambda) {
  return(1/lambda*exp(-1/lambda*data))
}

kappa_values=seq(0,10,0.01)

# Function for computing the posterior distribution
posteriorDistrib = function(kappa, lambda, data, my) {
  likelihood=vonMisesDistrib(kappa, data, my)
  prior=exponDistrib(kappa, lambda)
  return(likelihood*prior)
}

posteriorLikelihood=posteriorDistrib(kappa_values, lambda, data_radian, my)
posterior.df=data.frame(kappa=kappa_values, likelihood=posteriorLikelihood)
sumOfPosterior=sum(posterior.df$likelihood)
posterior.df$likelihood=posterior.df$likelihood*(1/sumOfPosterior)
final_sum=sum(posterior.df$likelihood)
plot(kappa_values, posterior.df$likelihood, xlab="Kappa", ylab="Likelihood",
     main="Posterior likelihood for different kappavalues", type="l", col="blue")

## As seen in the plot the likelihood of the posterior peaks between 2 and 4 and then dies off for larger
## kappa-values.

## b) Find the (approximate) posterior mode of k from the information in a).

# Puts likelihood values with corresponding kappa-values to be able to retrieve the kappa-value corresponding to
## the highest likelihood (mode)

posteriorMode=subset(posterior.df, likelihood==max(likelihood), kappa)
print(posteriorMode$kappa)

## The approximated posterior mode is found to be 2.12.
