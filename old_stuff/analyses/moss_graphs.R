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
v$l_mass = log(v$mass_per_plant_mg)
b$l_mass = log(b$mass_per_plant_mg)
v$infls_per_plant = v$infls/v$final_count 
b$infls_per_plant = b$infls/b$final_count

N = table(v[, c(2,3)])
f_success = formula("prop_success ~ treatment + stress") 
f_size = formula("l_mass  ~ treatment + stress")
f_infls = formula("infls_per_plant ~ treatment + stress")

gen_results = function(df, formu, n, sp_name){ 
  res = aggregate(formu, data = df, 'mean')  
  names(res) [ 3] <- "result" 
  res$stand_dev = aggregate(formu, data = df, 'sd')[, 3]
  res$stand_err = res$stand_dev/sqrt(n)
  res$species = sp_name
  return(res)
} 
v_success = gen_results(v, f_success, 9, "vulpia")
b_success = gen_results(b, f_success, 9, "bromus")    
prop_success = rbind(v_success, b_success)

v_size = gen_results(v, f_size, 9, "vulpia")
b_size = gen_results(b, f_size, 9, "bromus")
size = rbind(v_size, b_size)

v_infls = gen_results(v, f_infls, 9, "vulpia")
b_infls = gen_results(b, f_infls, 9, "bromus")
infls = rbind(v_infls, b_infls)


fig1 = ggplot(data = prop_success, aes( x = stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = treatment, 
                                         color = treatment))
fig1 + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab( "proportion germinated/survived") + 
  facet_grid( . ~ species ) + 
  theme_bw() 

fig2 = ggplot(data = size, aes( x = stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = treatment, 
                                         color = treatment))
fig2 + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab("size per plant (log g)") + 
  facet_grid(. ~ species ) + 
  theme_bw() 


fig3 = ggplot(data = infls, aes( x = stress, y = result, ymin = result - stand_err, 
                                ymax = result + stand_err, group = treatment, 
                                color = treatment))
fig3 + geom_point(position = position_dodge(width = 0.2)) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab('infls per plant') + 
  facet_grid(. ~ species ) + 
  theme_bw() 

