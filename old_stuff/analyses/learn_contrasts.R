rm(list = ls())
library(lsmeans)

trials = 10
x = rep(c(0, 1, 2 , 3), 25) 
y = 1 + 2*x + rnorm(100, 1, 1)
y2 = rbinom(100, trials, 0.1 + x/10)

df = data.frame(y = y, y2 = y2, x = as.factor(x))

sds = aggregate( y ~ x, df, sd)
means = aggregate(y ~ x, df, mean)
mean_prop = aggregate(y2/trials ~ x, df, mean)
means
mean_prop 

sds[2]
se = sds[2]/sqrt(25)
se

m1 = lm(y ~ x, data = df)
summary(m1)
aov(m1)
lsmeans = lsmeans(m1, ~ x)
predicted = as.data.frame(predict(m1, data.frame(x = as.factor(c(0,1,2,3))), se.fit = TRUE))
predicted


m2 = glm(cbind(y2, trials - y2) ~ x, data = df, family = binomial)

summary(m1)
aov(m1)

summary(m2)
aov(m2)

lsmeans(m1, ~ x)

lsmeans(m2,  ~x )




