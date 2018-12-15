remove(list = ls())
library(nlme)
library(lme4)
library(MASS)
library(plyr)
library(lmtest)
library(lsmeans)
library(multcomp)
library(boot)
library(ggplot2)
library(grid)
library(pbkrtest)
library(dplyr)
library(tidyr)
library(sjPlot)
library(xtable)

source( 'data_and_analysis/moss_theme.R')

xlab_distance = 'Position towards NW on environmental gradient (m)'

gen_results = function(df, formu, n, sp_name){ 
  res = aggregate(formu, data = df, 'mean')  
  names(res) [ 3] <- "result" 
  res$n = aggregate(formu, data = df, 'length')[,3]
  res$stand_dev = aggregate(formu, data = df, 'sd')[, 3]
  res$stand_err = res$stand_dev/sqrt(n)
  res$species = sp_name
  return(res)
} 

contrast_list <- list( `L:mp - bs` = c(1, -1, 0, 0, 0, 0), `L:mp - mr` = c(1, 0, -1, 0, 0, 0), `L:bs - mr` = c(0, 1, -1, 0, 0, 0),
                       `H:mp - bs` = c(0, 0, 0, 1, -1, 0 ), `H:mp - mr` = c(0, 0, 0, 1, 0, -1), `H:bs - mr` = c(0, 0, 0, 0, 1, -1),
                       `mp:L-H` = c(1, 0, 0, -1, 0, 0), `bs:L-H` = c(0, 1, 0, 0, -1, 0), `mr:L-H` = c(0, 0, 1, 0, 0, -1))

treatment_contrasts <- list(`mp - bs` = c(1, -1, 0, 1, -1, 0)/2, `mp  - mr` = c(1, 0, -1, 1, 0, -1)/2 , `bs - mr` = c(0, 1, -1, 0, 1, -1)/2)


v = read.csv("data_and_analysis/moss_removal_experiment/vulpia/vulpia_data.csv")
b = read.csv('data_and_analysis/moss_removal_experiment/bromus/brdi_data.csv')

v$trials = 5
b$trials = 5

v$stress <- factor(v$stress, levels = c('low', 'high'))
b$stress <- factor (b$stress, levels = c('low', 'high'))
levels(v$stress) = c('Low stress', 'High stress')
levels(b$stress) = c('Low stress', 'High stress')
levels(v$treatment) = c('Bare sand', 'Moss patch', 'Moss removed')
levels(b$treatment) = c('Bare sand', 'Moss patch', 'Moss removed')

v$treatment <- factor(v$treatment, levels= c('Moss patch', 'Bare sand', 'Moss removed'))
b$treatment <- factor(b$treatment, levels= c('Moss patch', 'Bare sand', 'Moss removed'))

v$prop_success = v$final_count/v$trials
b$prop_success = b$final_count/b$trials
v$l_mass = log(v$mass_per_plant_mg)
b$l_mass = log(b$mass_per_plant_mg)
v$infls_per_plant = v$infls/v$final_count 
b$infls_per_plant = b$infls/b$final_count

v$logit_ps = logit(v$prop_success)

N = table(v[, c(2,3)])
f_success = formula("prop_success ~ treatment + stress") 
f_size = formula("l_mass  ~ treatment + stress")
f_infls = formula("infls_per_plant ~ treatment + stress")

v_success = gen_results(v, f_success, 9, "Vulpia")
b_success = gen_results(b, f_success, 9, "Bromus")    
prop_success = rbind(v_success, b_success)
names(prop_success)[1:2] <- c('Treatment', 'Stress')

#

v_size = gen_results(v, f_size, 9, "Vulpia")
b_size = gen_results(b, f_size, 9, "Bromus")
size = rbind(v_size, b_size)
names(size)[1:2] <- c('Treatment', 'Stress')

#

v_infls = gen_results(v, f_infls, 9, "Vulpia")
b_infls = gen_results(b, f_infls, 9, "Bromus")
infls = rbind(v_infls, b_infls)
names(infls)[1:2] <- c('Treatment', 'Stress')

#### --- Vulpia final count -------------------------------------------

vm1 = glmer(cbind(final_count, trials - final_count) ~  stress + treatment + treatment:stress + (1|block), data = v, family = binomial)
vm2 = update(vm1 , . ~ . - treatment:stress)
vm3 = update(vm2, . ~ . - treatment )
vm4 = update(vm2, . ~ . - stress )

summary(vm1)

lrtest(vm1, vm2)
anova(vm1, vm2)

drop1(vm2, test = 'Chisq')
drop1(vm1, test = 'Chisq')

drop1(vm3, test = 'Chisq')
drop1(vm4, test = 'Chisq')

anova(vm1, test = 'Chisq')

summary(vm1)

vm1_lsmeans = lsmeans(vm1, pairwise ~ treatment*stress, adjust = 'sidak')
v1mct.contrasts <-contrast(vm1_lsmeans, treatment_contrasts, options = list(adjust = 'sidak'))
v1mct.contrasts <- as.data.frame( print(v1mct.contrasts))

summary(glm( data = subset(v, stress == 'Low stress'), cbind( final_count, trials - final_count) ~ treatment, family = 'binomial'))
summary(glm( data = subset(v, stress == 'High stress'), cbind( final_count, trials - final_count) ~ treatment, family = 'binomial'))

v_success <- merge( v_success, data.frame( print( vm1_lsmeans$lsmeans ) ), by = c('treatment', 'stress'))

v_success <- v_success %>% 
  mutate( lsmeans_bt = inv.logit(lsmean), LCL_bt = inv.logit(asymp.LCL), UCL_bt = inv.logit(asymp.UCL))

v_success

#### --------- Vulpia mass per plant 
vb1 <- lmer( l_mass ~ stress*treatment + (1|block), data = v, REML = FALSE)
vb2 <- update(vb1, . ~ . - treatment:stress )
vb3 <- update(vb2, . ~ . - treatment )
vb4 <- update(vb2, . ~ . - stress )

lrtest(vb1, vb2)
anova(vb1, vb2)

# pbkrtest::KRmodcomp(vb1, vb2)

drop1(vb1, test = 'Chisq')
drop1(vb2, test = 'Chisq')

vb1 <- update(vb1, REML = TRUE)

vb1_lsmeans <- lsmeans(vb1, pairwise ~ treatment*stress)

contrast(vb1_lsmeans, treatment_contrasts, options = list( adjust = 'sidak'))

vb1_lsmeans <- data.frame( print(vb1_lsmeans$lsmeans) ) 

v_size <- merge( v_size, vb1_lsmeans , by = c('treatment', 'stress')) 

v_size <- v_size %>% 
  mutate( lsmean_bt = exp(lsmean), lcl_bt = exp(lower.CL), ucl_bt = exp(upper.CL))

v_size 

#### ----------- Vulpia infls. per plant 
vfdata = v[ !v$final_count == 0, ]

vfdata <- vfdata %>% select( stress, treatment, block , infls, final_count ) 

sum( vfdata$infls == 0 )

vfdata %>% group_by ( block ) %>% summarise( n() )

vfdata$block <- factor( vfdata$block)

vfdata$log_final_count <- log(vfdata$final_count)

vfdata_summary <-  vfdata %>% 
  group_by ( stress, treatment ) %>% 
  mutate( ipp = infls / final_count ) %>%
  summarise( m_infls = mean(infls), m_final_count = mean(final_count), m_infls_per_plant = mean(ipp))

vf1_glm_poisson <- glm( infls ~ offset( log_final_count) + stress*treatment, data = vfdata, family = 'poisson')
summary(vf1_glm_poisson)

vf1_glm <- glm( infls ~ offset( log_final_count) + stress*treatment, data = vfdata, family = 'quasipoisson')
summary(vf1_glm)

anova(vf1_glm, test = 'F')

drop1(vf1_glm, test = 'Chisq')
vf2_glm <- update( vf1_glm, . ~ . - treatment:stress)
drop1(vf2_glm, test = 'Chisq')
anova(vf1_glm, vf2_glm, test = 'Chisq')

vf1_lsmeans <- lsmeans(vf1_glm, pairwise ~ treatment*stress)

vf1_glm_nb <- glm.nb(infls ~ offset(log_final_count) + stress*treatment, data = vfdata)
summary(vf1_glm_nb)

1 - pchisq(4.9, 2)

vf1glmer <- glmer(infls ~ offset( log_final_count) + stress*treatment + (1|block), data = vfdata, family = 'poisson')
summary(vf1glmer)

vf1_mct <- lsmeans(vf1_glm, pairwise ~ treatment|stress)
vf1_lsmeans <- data.frame( print(vf1_mct$lsmeans) ) 

drop1(vf1_glm_nb, test = 'Chisq')
drop1(update(vf1_glm_nb, . ~ . - stress:treatment), test = 'Chisq')

#

vfdata$pred_glm <- predict(vf1_glm, type = 'response')
vfdata$pred_glm.nb <- predict(vf1_glm_nb, type = 'response')

plot_data <- vfdata %>% 
  gather(key = stat, v = value, c(infls, pred_glm, pred_glm.nb)) 

plot_data <- plot_data %>% 
  group_by( stress, treatment, stat) %>% 
  summarise( val = mean(value), ipp = sum(val)/sum(final_count))

ggplot(data = plot_data, aes( x = treatment, y = val, fill = stat) ) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_wrap( ~ stress)

vf_pred <- vfdata %>% select( treatment, stress) %>% mutate(log_final_count = 0) %>% distinct()

vf_pred <- data.frame(vf_pred, predict( vf1_glm_nb, newdata = vf_pred, se.fit = TRUE) )

vf_pred <- vf_pred %>% mutate( lcl = fit - se.fit*1.96,  ucl = fit + se.fit*1.96, fit_bt = exp(fit), lcl_bt = exp(lcl), ucl_bt = exp(ucl))

vfdata_summary <- merge( v_infls, vf_pred %>% select( - log_final_count), by = c('treatment', 'stress'))

vfdata_summary <- vfdata_summary %>% arrange( stress, treatment )

vfdata_summary <- vfdata_summary %>% select(  -lcl, -ucl)

vfdata_summary

#### ---- Bromus Final Count ------------------------------------------------------------------------------------------------------------------ 
bm1 <- glmer(cbind(final_count, trials - final_count) ~ treatment + stress + treatment:stress + (1|block), data = b, family = binomial)
bm2 <- update(bm1, . ~ . - treatment:stress)
bm3 <- update(bm2, . ~ . - treatment)
bm4 <- update(bm3, . ~ . - stress )

anova(bm1, bm2)
lrtest(bm1, bm2)
drop1(bm1, test = 'Chisq')
drop1(bm2, test = 'Chisq')

summary(bm1)

bm_mct = lsmeans(bm1, pairwise ~ treatment*stress)
bm_mct$lsmeans

contrast(bm_mct, contrast_list, by = NULL, options = list( adjust = 'none'))
contrast(bm_mct, contrast_list, by = NULL, options = list( adjust = 'bonferroni'))
contrast(bm_mct, contrast_list, by = NULL, options = list( adjust = 'sidak'))

summary(contrast(bm_mct, 'pairwise'))[c(1:2, 6, 13:15, 3, 8, 12), ]

bm_lsmeans <- data.frame( print(bm_mct$lsmeans))

b_success <- merge( b_success, bm_lsmeans, by = c('treatment', 'stress'))

b_success <- b_success %>% 
  mutate( lsmeans_bt = inv.logit(lsmean), LCL_bt = inv.logit(asymp.LCL), UCL_bt = inv.logit(asymp.UCL)) %>% 
  arrange( stress, treatment)

#### ----------- Bromus mass per plant ---------- 
bmass1 <- lmer( l_mass ~ stress*treatment + (1|block), data = b, REML = FALSE)
bmass2 <- update(bmass1, . ~ . - treatment:stress )
bmass3 <- update(bmass2, . ~ . - treatment ) 
bmass4 <- update(bmass2, . ~ . - stress )


anova(bmass1, bmass2)
drop1(bmass1, test = 'Chisq')
drop1(bmass2, test = 'Chisq')

#pbkrtest::getLRT(bmass1, bmass2)
#pbkrtest::KRmodcomp(bmass1, bmass2)

bmass1 <- update(bmass1, REML = TRUE)
summary(bmass1)
anova(bmass1)


bmass_lsmeans <- lsmeans(bmass1, pairwise ~ treatment*stress)
bmass_lsmeans
contrast(bmass_lsmeans, treatment_contrasts, by = NULL, options = list(adjust = 'sidak') )

bmass_lsmeans2 <- lsmeans(bmass1, pairwise ~ treatment, adjust = 'sidak')

bmass_lsmeans2

bmass_lsmeans <- data.frame ( print( bmass_lsmeans$lsmeans) ) 

b_size <- merge( b_size, bmass_lsmeans, by = c('treatment', 'stress'))

b_size <- b_size %>% arrange(stress, treatment )

#### ----------- Bromus infls. per plant ----------------------------------------
bfdata <- b[ !b$final_count == 0 , ]

bfdata <- bfdata %>% select( treatment, stress, infls, final_count, block) 

bfdata$block <- factor( bfdata$block)

bfdata$log_final_count <- log(bfdata$final_count)

bf1_glm <- glm( infls ~ stress*treatment + offset( log_final_count), data = bfdata, family = 'quasipoisson')
bf2_glm <- update(bf1_glm, . ~ .  - stress:treatment)
summary(bf2_glm)
summary(bf1_glm)

anova(bf1_glm , test = 'F')
drop1(bf1_glm, test = 'F')

plot( bfdata$infls, predict(bf1_glm, type = 'response'))

bfdata$pred_glm <- predict(bf1_glm, type = 'response')

plot_data <- bfdata %>% gather(key = stat, v = value, c(infls, pred_glm)) 

plot_data <- plot_data %>% group_by( stress, treatment, stat) %>% summarise( val = mean(value) )

ggplot(data = plot_data, aes( x = treatment, y = val, fill = stat) ) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_wrap( ~ stress)


bf1_lsmeans <- lsmeans(bf1_glm, pairwise ~ treatment*stress)

contrast(bf1_lsmeans, contrast_list, options = list(adjust = 'sidak'))

bf_pred <- bfdata %>% select( treatment, stress) %>% mutate(log_final_count = 0) %>% distinct()

bf_pred <- data.frame(bf_pred, predict(bf1_glm, newdata = bf_pred, se.fit = TRUE))

bf_pred <- bf_pred %>% mutate( fit_bt = exp(fit), ucl = fit + se.fit*1.96, lcl = fit -se.fit*1.96, lcl_bt = exp(lcl), ucl_bt = exp(ucl))

bfdata_summary <- merge( b_infls, bf_pred %>% select( -log_final_count), by = c('treatment', 'stress' ) ) %>% arrange( stress, treatment )

# ---------------------- print analysis tables ------------------------------ 
pred_labels <-  c('High stress', 'Bare sand', 'Moss removed', 'High stress:Bare sand', 'High stress:Moss removed')

# vulpia -------------------------------------------------------------------------------------------------------------

vm1_table <- sjt.glmer( vm1, pred.labels = pred_labels, 
                        file = 'data_and_analysis/figures/v_success_glmer_table.html') 


vm1_lrt <- data.frame( rbind(drop1(vm1, test = 'Chisq'),  drop1(vm2, test = 'Chisq')[2:3,] ) )
vm1_lrt$`terms removed` <- row.names(vm1_lrt)
vm1_lrt <- vm1_lrt %>% mutate(  AIC = round(AIC, 0), LRT = round(LRT, 2), Pr.Chi. = round(Pr.Chi., 3))

print( xtable( vm1_lrt), 'html', 'data_and_analysis/figures/vulpia_success_lrt.html')

v_success <- v_success %>% select( treatment, stress, result, n, stand_dev, stand_err, lsmeans_bt, LCL_bt, UCL_bt) %>% 
  rename( `proportion survived` = result, sd = stand_dev, `S.E.` = stand_err, `fitted mean` = lsmeans_bt, 
          `Lower 95% C.L.` = LCL_bt, `Upper 95% C.L.` = UCL_bt)

print( xtable( v_success), 'html', 'data_and_analysis/figures/vulpia_success_stats.html')

#

vm1_table <- sjt.lmer( vb1, pred.labels  = pred_labels, 
                        file = 'data_and_analysis/figures/v_size_lmer_table.html') 

vm1_lrt <- data.frame( rbind(drop1(vb1, test = 'Chisq'),  drop1(vb2, test = 'Chisq')[2:3,] ) )
vm1_lrt$`terms removed` <- row.names(vm1_lrt)
vm1_lrt <- vm1_lrt %>% mutate(  AIC = round(AIC, 1), LRT = round(LRT, 3), Pr.Chi. = round(Pr.Chi., 4))

print( xtable( vm1_lrt), 'html', 'data_and_analysis/figures/vulpia_size_lrt.html')

v_size <- v_size %>% select( treatment, stress, result, n, stand_dev, stand_err, lsmean, lower.CL, upper.CL) %>% 
  rename( `average size (log g)` = result, sd = stand_dev, `S.E.` = stand_err, `fitted mean` = lsmean, 
          `Lower 95% C.L.` = lower.CL, `Upper 95% C.L.` = upper.CL)

print( xtable( v_size), 'html', 'data_and_analysis/figures/vulpia_size_stats.html')

#
print( xtable(vfdata_summary), 'html', file = 'data_and_analysis/figures/vulpia_flowers_summary.html' ) 

vf1_table <- sjt.glmer( vf1_glm, pred.labels = pred_labels, 
           file = 'data_and_analysis/figures/v_infls_glm_table.html') 

# bromus  output ----------------------------------------------------- 
bm1_table <- sjt.glmer( bm1, pred.labels = pred_labels, 
                        file = 'data_and_analysis/figures/bromus_success_glmer_table.html') 

bm1_lrt <- data.frame( rbind(drop1(bm1, test = 'Chisq'),  drop1(bm2, test = 'Chisq')[2:3,] ) )
bm1_lrt$`terms removed` <- row.names(bm1_lrt)
bm1_lrt <- bm1_lrt %>% mutate(  AIC = round(AIC, 1), LRT = round(LRT, 3), Pr.Chi. = round(Pr.Chi., 4))

bm1_lrt

print( xtable( bm1_lrt), 'html', 'data_and_analysis/figures/bromus_success_lrt.html')

print( xtable ( data.frame( print( lsmeans(bm1, pairwise ~ treatment|stress))$contrasts)), 'html', 'data_and_analysis/figures/bromus_succuss_mct.html')

b_success <- b_success %>% select( treatment, stress, result, n, stand_dev, stand_err, lsmeans_bt, LCL_bt, UCL_bt) %>% 
  rename( `proportion survived` = result, sd = stand_dev, `S.E.` = stand_err, `fitted mean` = lsmeans_bt, 
          `Lower 95% C.L.` = LCL_bt, `Upper 95% C.L.` = UCL_bt)

print( xtable( b_success), 'html', 'data_and_analysis/figures/bromus_success_stats.html')

#

bm1_table <- sjt.lmer( bmass1, pred.labels  = pred_labels, 
                       file = 'data_and_analysis/figures/bromus_size_lmer_table.html') 

bm1_lrt <- data.frame( rbind(drop1(bmass1, test = 'Chisq'),  drop1(bmass2, test = 'Chisq')[2:3,] ) )
bm1_lrt$`terms removed` <- row.names(bm1_lrt)
bm1_lrt <- bm1_lrt %>% mutate(  AIC = round(AIC, 0), LRT = round(LRT, 2), Pr.Chi. = round(Pr.Chi., 3))

print( xtable( bm1_lrt), 'html', 'data_and_analysis/figures/bromus_size_lrt.html')


b_size <- b_size %>% select( treatment, stress, result, n, stand_dev, stand_err, lsmean, upper.CL, lower.CL) %>% 
  rename( `average size (log g)` = result, sd = stand_dev, `S.E.` = stand_err, `fitted mean` = lsmean, 
          `Lower 95% C.L.` = lower.CL, `Upper 95% C.L.` = upper.CL)

print( xtable( b_size), 'html', 'data_and_analysis/figures/bromus_size_stats.html')

#
print( xtable(bfdata_summary), 'html', file = 'data_and_analysis/figures/bromus_flowers_summary.html' ) 

bf1_table <- sjt.glmer( bf1_glm, pred.labels = pred_labels, 
                        file = 'data_and_analysis/figures/b_infls_glm_table.html') 
bf1_table
summary(bf1_glm)
