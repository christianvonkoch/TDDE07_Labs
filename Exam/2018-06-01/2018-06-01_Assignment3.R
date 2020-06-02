## c) Choose between three models where two of them use Beta prior and the last one assumes p=0.5. Which model 
## should be chosen?

model1=choose(10,3)*gamma(4)*gamma(8)*gamma(2)/gamma(12)
model2=choose(10,3)*gamma(7)*gamma(11)*gamma(8)/(gamma(4)*gamma(4)*gamma(18))
model3=choose(10,3)*0.5^10
model1_norm=model1/sum(c(model1, model2, model3))
model2_norm=model2/sum(c(model1, model2, model3))
model3_norm=model3/sum(c(model1, model2, model3))
