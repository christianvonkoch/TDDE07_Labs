## Maximizing posterior expected utility.

post_dens = function(x) {
  return(gamma(6+x)/gamma(13+x))
}

barplot(post_dens(seq(0,10,1)), type="l")

# x6=10 seems to yield a low enough probability to be an upper bound for sum

posterior_prob=post_dens(seq(0,10))
posterior_prob=posterior_prob/sum(posterior_prob)
exp_util=c()
for (k in 0:10) {
  exp_util=c(exp_util,(2^k-3))
}

exp_post_dens=sum(posterior_prob*exp_util)
print(exp_post_dens)
