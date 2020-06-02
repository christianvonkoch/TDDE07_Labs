## c) Do Bayesian model comparison of two models, both geometric likelihood, first with beta prior alpha=0.5, beta=0.5
##  and second is null model assuming theta=0.5. Prior probabilities are p(M1)=0.1 and p(M2)=9/10

alpha=1/2
beta=1/2
data=c(2, 1, 12)
n=length(data)

# We know that marginal likelihood for data is likelihood*prior/posterior which when derived yields the below function

margLikelihood1=gamma(alpha+beta)*gamma(alpha+n)*gamma(beta+sum(data))/(gamma(alpha)*gamma(beta)*
                                                                          gamma(beta+sum(data)+alpha+n))

# Since we have a constant prior the marginal likelihood for model 2 will be the likelihood with theta set to the 
## constant

margLikelihood2=0.5^n*(1-0.5)^sum(data)

prior1=1/10
prior2=9/10

posterior1=margLikelihood1*prior1
posterior2=margLikelihood2*prior2
posterior1_norm=posterior1/sum(c(posterior1, posterior2))
posterior2_norm=posterior2/sum(c(posterior1, posterior2))

## Model 1 has a higher posterior probability than 1 so this model should be chosen. 


