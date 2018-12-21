remove(list = ls())
setwd("~/Documents/moss_manuscript/moss_data")

library(MASS)
library(nlme)
library(lme4)
library(ggplot2)
library(plyr)
library(lmtest)
library(boot)
library(lme4)
library(lsmeans)

v = read.csv("vulpia/vulpia_data.csv")
b = read.csv('bromus/brdi_data.csv')


v$trials = 5
b$trials = 5

v$prop_success = v$final_count/v$trials
b$prop_success = b$final_count/b$trials

vm1 = glmer(cbind(final_count, trials - final_count) ~ treatment + stress + treatment:stress + (1|block), data = v, family = binomial)
vm2 = update(vm1, .~ . - treatment:stress)
vm3 = update(vm2, .~ . - stress)
vm4 = update(vm2, .~ . - treatment)
summary(vm1)
summary(vm2)
summary(vm3)
summary(vm4)
AIC(vm1, vm2, vm3,vm4)

vlsmeans = lsmeans(vm1, pairwise ~ treatment|stress, adjust = "tukey")
vlsmeans
vresults = vlsmeans$'treatment:stress lsmeans'
vresults

vresults$upper_se = vresults[, 3] + vresults[, 4]
vresults$lower_se = vresults[, 3] - vresults[, 4] 
vresults$bt_mean = inv.logit(vresults[, 3])
vresults$bt_upper_se = inv.logit(vresults[ , 8])
vresults$bt_lower_se = inv.logit(vresults[, 9])
vresults

vulp_plot = ggplot(data = vresults, aes( x = stress, y = bt_mean, ymin = bt_lower_se, ymax = bt_upper_se, 
                                         group = treatment, color = treatment))
vulp_plot + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) +  ylim (0, 1) + 
  theme_bw() 

vm2 = glm(cbind(final_count, trials - final_count) ~ stress + treatment + stress:treatment, data = v, family = 'binomial')
df1 = expand.grid(stress = c('high', 'low'), treatment = c('BC', 'MC', 'MR'))
df1 = df1[order(df1$stress), ]

predicted = data.frame(predict(vm2,df1, type = 'response', se.fit= TRUE ))
predicted2 = data.frame(predict(vm2, df1, se.fit = TRUE ))
predicted2
df1 = cbind(df1, predicted)
df1

vulp_plot2 = ggplot(data = df1, aes(x = stress, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, 
                                   group = treatment, fill = treatment, color = treatment))
vulp_plot2 + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) +  ylim (0, 1) + 
  theme_bw() 

bm1 = glmer(cbind(final_count, trials - final_count) ~ treatment + stress + treatment:stress + (stress|block), data = b, family = binomial)
bm2 = update(bm1, .~ . - treatment:stress)
bm3 = update(bm2, .~ . - stress)
summary(bm1)
summary(bm2)
summary(bm3)
AIC(bm1, bm2, bm3)


blsmeans = lsmeans(bm2,  ~ treatment*stress)
blsmeans
blsmeans = blsmeans$'treatment:stress lsmeans'
blsmeans_back = cbind(blsmeans[, c(1,2)], data.frame(lapply(blsmeans[ , c(3:7)], FUN = 'inv.logit')))



bm2 = glm(cbind(final_count, trials-final_count) ~ stress + treatment + stress:treatment, data = b, family = 'binomial')
df2 = expand.grid(stress = c('high', 'low'), treatment = c('BC', 'MC', 'MR'))
bt_predicted = data.frame(predict(bm2,df2, type = 'response', se.fit= TRUE ))
df2 = cbind(bt_predicted, df2)

brom_plot = ggplot(data = df2, aes(x = stress, y = fit, ymin = fit - se.fit, ymax = fit + se.fit, 
                                   group = treatment, fill = treatment, color = treatment))
brom_plot + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar ( position = position_dodge(), width = 0.2) +  ylim(0,1) + 
  theme_bw() 


v$infls_per_plant = v$infls/v$final_count
b$infls_per_plant = b$infls/b$final_count


#### biomass per plant analysis 
hist(v$mass_per_plant_mg)
hist(b$mass_per_plant_mg)
hist(log(v$mass_per_plant_mg))
hist(log(b$mass_per_plant_mg))
hist(v$infls_per_plant)
hist(b$infls_per_plant)
hist()


v_clean = v[ !is.na(v$mass_per_plant_mg), ]
v_clean[, names(v_clean) == 'mass_per_plant_mg']

v_clean = v[ !is.na(v$mass_per_plant_mg), ]
b_clean = b[ !is.na(b$mass_per_plant_mg), ]

v_clean$lmpp = log(v_clean[, "mass_per_plant_mg"])
b_clean$lmpp = log(b_clean[, "mass_per_plant_mg"])


v_bm_mod = lm( )



# Error bars represent standard error of the mean
col = "PuBuGn"
ylabel = c("Final mass per plant (log g)", "Inflorescences per plant", "No. of plants at end of season per replicate (out of five seeds planted)")
title = c("Vulpia", "Bromus")





######## Vulpia plots
pdf(file = 'vulpia_1.pdf', 6, 6)
create_bar_plot(v_lmpp, ylabel = ylabel[1], titlename= title[1])
dev.off()

pdf(file = 'vulpia_2.pdf', 6,6)
create_bar_plot(v_ipp, ylabel = ylabel[2], titlename= title[1])
dev.off()

pdf(file = 'vulpia_3.pdf', 6,6)
create_bar_plot(v_fc, ylabel = ylabel[3], titlename= title[1])
dev.off()

######## Bromus plots

pdf(file = 'bromus_1.pdf', 6, 6)
create_bar_plot(b_lmpp, ylabel = ylabel[1], titlename= title[2])
dev.off()

pdf(file = 'bromus_2.pdf', 6, 6)
create_bar_plot(b_ipp, ylabel = ylabel[2], titlename = title[2]) 
dev.off()

pdf(file = 'bromus_3.pdf', 6, 6)
create_bar_plot(b_fc, ylabel = ylabel[3], titlename = title[2]) 
dev.off()




