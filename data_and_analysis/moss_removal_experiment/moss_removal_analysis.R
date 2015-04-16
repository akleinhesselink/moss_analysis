
source('moss_removal_experiment_figures.R')

library(lme4)
library(MASS)
library(nlme)
library(plyr)
library(lmtest)
library(lsmeans)
library(pbkrtest)


#### --- Vulpia final count
vm1 = glmer(cbind(final_count, trials - final_count) ~  stress + treatment + treatment:stress + (1|block), data = v, family = binomial)
vm2 = update(vm1, .~ . - treatment:stress)
vm3 = update(vm2, .~ . - treatment)
vm4 = update(vm3, .~ . - stress)
summary(vm1)
summary(vm2)
summary(vm3)
summary(vm4)
AIC(vm1, vm2, vm3,vm4)
anova(vm1)
anova(vm4, vm3, vm2, vm1)
PBmodcomp(vm1, vm2)
PBmodcomp(vm2, vm3)
PBmodcomp(vm3, vm4)

mctest = glht(vm1, linfct = mcp(treatment = c(0, -1, 1)))
summary(mctest, test = adjusted("single-step"))
v1mct = lsmeans(vm1, pairwise ~ treatment|stress, adjust = 'tukey')

v_success$lsmean = inv.logit(v1mct[[1]]$lsmean)
v_success$lsmeanLCL = inv.logit(v1mct[[1]]$asymp.LCL)
v_success$lsmeanUCL = inv.logit(v1mct[[1]]$asymp.UCL)
v_success

ggplot(data = v, aes (x = stress, y = fitted(vm1), color = treatment)) + geom_point(position = position_dodge(width = 0.2))
ggplot(data = v, aes(x = stress, y = fitted(vm1), color = treatment)) + geom_point(position = position_dodge(width = 0.2))

#### --------- Vulpia mass per plant 
vmass = lmer( l_mass ~ treatment + stress + treatment:stress + (1|block), data = v)
vmass2 = update(vmass, .~. - treatment:stress)
vmass3 = update(vmass2, . ~ . - treatment)
vmass4 = update(vmass3, . ~ . - stress)
summary(vmass)
anova(vmass4, vmass3, vmass2, vmass)
KRmodcomp(vmass, vmass2)
KRmodcomp(vmass2, vmass3)
KRmodcomp(vmass3, vmass4)
v1massmct = lsmeans(vmass, pairwise ~ treatment|stress)

v_size$lsmean = v1massmct[[1]]$lsmean
v_size$lsmeanLCL = v1massmct[[1]]$lower.CL
v_size$lsmeanUCL = v1massmct[[1]]$upper.CL
v_size[, 10:12] = exp(v_size[, 7:9])


hist(v$infls_per_plant)
hist(v$infls)
#### ----------- Vulpia infls. per plant 
vfdata = v[ !v$final_count == 0, ]
vf1 = glmer(infls ~ stress + treatment + treatment:stress + (1|block), data = v, family  = 'poisson', 
            offset = final_count)

summary(vf1)

vftest = glmer(infls ~ offset(final_count) + stress + treatment + treatment:stress + (1|block), 
               data = vfdata, family = 'poisson')

vf2 = update(vf1, . ~ . - treatment:stress)
vf3 = update(vf2, . ~ . - treatment)
vf4 = update(vf3, . ~ . - stress)
summary(vf1)

anova(vf1, vf2, vf3, vf4)
PBmodcomp(vf1, vf2)
PBmodcomp(vf2, vf3)
PBmodcomp(vf3, vf4)

vfmct = lsmeans(vf1, pairwise ~ treatment + stress)
vfmct
plotlsmeans = data.frame(vfmct[[1]])

v_infls$result - exp(v_infls$lsmean)

v_infls$lsmean = vfmct[[1]]$lsmean
v_infls$lsmeanLCL = vfmct[[1]]$asymp.LCL
v_infls$lsmeanUCL = vfmct[[1]]$asymp.UCL

hist(predict(vf1))

vfdata$infls_pred = exp(predict(vf1)) + log(vfdata$final_count)  - 1
vfdata$infls
plot(infls ~ infls_pred, vfdata)
abline(0,1)


#### ---- Bromus Final Count 
bm1 = glmer(cbind(final_count, trials - final_count) ~ treatment + stress + treatment:stress + (1|block), data = b, family = binomial)
bm2 = update(bm1, .~ . - treatment:stress)
bm3 = update(bm2, .~ . - treatment)
bm4 = update(bm2, .~ . - stress)
summary(bm1)
anova(bm1)
anova(bm1, bm2, bm3)
anova(bm1)
PBmodcomp(bm1, bm2)
PBmodcomp(bm2, bm3)

blsmeans = lsmeans(bm1, pairwise ~ treatment|stress, adjust = "tukey")
blsmeans


#### ----------- Bromus mass per plant ---------- 
bmass = lmer( l_mass ~ treatment + stress + treatment:stress + (1|block), data = b)
bmass2 = update(bmass, . ~ . -treatment:stress)
bmass3 = update(bmass2, .~. -treatment)
bmass4 = update(bmass3, .~. -stress)
summary(bmass) 
anova(bmass, bmass2, bmass3, bmass4)
anova(bmass)

KRmodcomp(bmass, bmass2)
KRmodcomp(bmass2, bmass3)
KRmodcomp(bmass3, bmass4)

bmassmct = lsmeans(bmass, pairwise ~ treatment|stress, adjust = "tukey")

b_size$lsmean = bmassmct[[1]]$lsmean
b_size$lsmeanLCL = bmassmct[[1]]$lower.CL
b_size$lsmeanUCL = bmassmct[[1]]$upper.CL
b_size[, 10:12] = exp(b_size[, 7:9])
b_size



#### ----------- Bromus infls. per plant 
bf1 = glmer(infls ~ stress + treatment + treatment:stress + (1|block), data = b, offset=final_count, family  = 'poisson')
bf2 = update(bf1, . ~ . - treatment:stress)
bf3 = update(bf2, . ~ . - treatment)
bf4 = update(bf3, . ~ . - stress)
summary(bf1)
anova(bf1, bf2, bf3, bf4)
anova(bf1)
lrtest(bf1, bf2)
lrtest(bf2, bf3)
lrtest(bf3, bf4)

PBmodcomp(bf1, bf2)
PBmodcomp(bf2, bf3)
PBmodcomp(bf3, bf4)

lsmeans(bf1, pairwise ~ treatment|stress)


