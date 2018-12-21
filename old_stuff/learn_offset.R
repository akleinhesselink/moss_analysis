#### learn offsets

#library(devtools)
#install_github("lme4", user = "lme4")
library(lme4)
library(MASS)
library(nlme)
library(ggplot2)
library(plyr)
library(lmtest)
library(lsmeans)
library(pbkrtest)
library(car)
library(boot)

##### generate data 
nplants = rpois(100, 20)
fx = rep(c(5, 10), 50)
treat = rep(c('A', 'B'), 50)

ninfls = NA
for (i in 1:length(nplants)){  
  ninfls[i] = rpois(n = 1, nplants[i]*fx[i])
}
fx
hist(ninfls/nplants)


testdata = data.frame(nplants, ninfls, treatment = factor(treat)) 
m1 = glm(ninfls ~ treatment  , offset = nplants, family = 'poisson', data = testdata)
summary(m1)
lsm = lsmeans(m1, pairwise ~ treatment)[[1]]
lsm = lsm$lsmean

coef(m1)[1] + nplants[1]
coef(m1)[1] + coef(m1)[2] + nplants[2]
coef(m1)[1]
coef(m1)[2] 


mypred = coef(m1)[1] + coef(m1)[2]*c(1,0) + nplants
predict(m1) - nplants

plot(predict(m1), mypred)
abline(0,1)

coef(m1)[1] + nplants

lsm[1] - coef(m1)[1]
lsm[2] - coef(m1)[1] - coef(m1)[2]

lsm[1] + mean(nplants)


plot(predict(m1), ninfls)
predict(m1) - (nplants)

plot(predict(m1) - nplants, ninfls/nplants)


log(20.06) 
log(23)


plot(predict(m1), ninfls/nplants)
plot(nplants, ninfls)

round(exp(predict(m1)), 3)

fx2 = seq(1, 100)
fx2
response = NA
for (i in 1:length(fx)){
  response[i] = rpois(1, fx2[i])
}

simpldf = data.frame(response, fx2)
msimp = glm(response ~ fx2, family = 'poisson', data = simpldf)
summary(msimp)
plot( exp(predict(msimp)), response)
abline(0, 1) 

simpldf$response