remove(list = ls())

library(tidyverse)
library(lme4)
library(MASS)
library(lmtest)
library(emmeans)
library(sjPlot)
library(xtable)

pred_labels <-
  c(
    'Intercept', 
    'High stress',
    'Bare sand',
    'Moss removed',
    'High stress:Bare sand',
    'High stress:Moss removed'
  )


gen_results = function(df, formu, n, sp_name) {
  res = aggregate(formu, data = df, 'mean')
  names(res) [3] <- "result"
  res$n = aggregate(formu, data = df, 'length')[, 3]
  res$stand_dev = aggregate(formu, data = df, 'sd')[, 3]
  res$stand_err = res$stand_dev / sqrt(n)
  res$species = sp_name
  return(res)
}

logit <- function( p ) { 
  
  log ( p / ( 1 - p ) )
  
}

inv.logit <- function( x ) {
  
  exp(x)/( 1 + exp(x))

}

write_lrt_table <- function( mod1, mod2, outfile ) { 
  lrt <-
    data.frame(rbind(drop1(mod1, test = 'Chisq'),  
                     drop1(mod2, test = 'Chisq')[2:3, ]))
  
  lrt$`terms removed` <- row.names(lrt)
  
  lrt <-
    lrt %>% mutate(
      AIC = round(AIC, 0),
      LRT = round(LRT, 2),
      Pr.Chi. = round(Pr.Chi., 3)
    )
  
  print(xtable(lrt), 'html', outfile )
}


contrast_list <-
  list(
    `L:mp - bs` = c(1, -1, 0, 0, 0, 0),
    `L:mp - mr` = c(1, 0, -1, 0, 0, 0),
    `L:bs - mr` = c(0, 1, -1, 0, 0, 0),
    `H:mp - bs` = c(0, 0, 0, 1, -1, 0),
    `H:mp - mr` = c(0, 0, 0, 1, 0, -1),
    `H:bs - mr` = c(0, 0, 0, 0, 1, -1),
    `mp:L-H` = c(1, 0, 0, -1, 0, 0),
    `bs:L-H` = c(0, 1, 0, 0, -1, 0),
    `mr:L-H` = c(0, 0, 1, 0, 0, -1)
  )

treatment_contrasts <-
  list(
    `mp - bs` = c(1, -1, 0, 1, -1, 0) / 2,
    `mp  - mr` = c(1, 0, -1, 1, 0, -1) / 2 ,
    `bs - mr` = c(0, 1, -1, 0, 1, -1) / 2
  )


v <- read.csv("data/vulpia_data.csv")
b <- read.csv('data/brdi_data.csv')

v$trials <- 5
b$trials <- 5

v$stress <- factor(v$stress, levels = c('low', 'high'))
b$stress <- factor (b$stress, levels = c('low', 'high'))
levels(v$stress) <- c('Low stress', 'High stress')
levels(b$stress) <- c('Low stress', 'High stress')
levels(v$treatment) <- c('Bare sand', 'Moss patch', 'Moss removed')
levels(b$treatment) <- c('Bare sand', 'Moss patch', 'Moss removed')

v$treatment <-
  factor(v$treatment, levels = c('Moss patch', 'Bare sand', 'Moss removed'))
b$treatment <-
  factor(b$treatment, levels = c('Moss patch', 'Bare sand', 'Moss removed'))

v$prop_success <- v$final_count / v$trials
b$prop_success <- b$final_count / b$trials
v$l_mass <- log(v$mass_per_plant_mg)
b$l_mass <- log(b$mass_per_plant_mg)
v$infls_per_plant <- v$infls / v$final_count
b$infls_per_plant <- b$infls / b$final_count

v$logit_ps <- logit(v$prop_success)

N <- table(v[, c(2, 3)])
f_success <- formula("prop_success ~ treatment + stress")
f_size <- formula("l_mass  ~ treatment + stress")
f_infls <- formula("infls_per_plant ~ treatment + stress")

v_success <- gen_results(v, f_success, 9, "Vulpia")
b_success <- gen_results(b, f_success, 9, "Bromus")
prop_success <- rbind(v_success, b_success)
names(prop_success)[1:2] <- c('Treatment', 'Stress')

# Final biomass data ------------------------ # 

v_size <- gen_results(v, f_size, 9, "Vulpia")
b_size <- gen_results(b, f_size, 9, "Bromus")
size <- rbind(v_size, b_size)
names(size)[1:2] <- c('Treatment', 'Stress')

# Final inflorescence numbers ----------------- # 
vfdata <- 
  v[!v$final_count == 0,] %>% 
  dplyr::select(stress, treatment, block , infls, final_count) %>% 
  mutate( block = factor(block), 
          log_final_count = log(final_count), 
          infls_per_plant = infls/final_count)


bfdata <- 
  b[!b$final_count == 0,] %>% 
  dplyr::select(stress, treatment, block , infls, final_count) %>% 
  mutate( block = factor(block), 
          log_final_count = log(final_count) ,
          infls_per_plant = infls/final_count)


v_infls <- gen_results(vfdata, f_infls, 9, "Vulpia")
b_infls <- gen_results(bfdata, f_infls, 9, "Bromus")
infls <- rbind(v_infls, b_infls)
names(infls)[1:2] <- c('Treatment', 'Stress')

#### --- Vulpia final count -------------------------------------------

vm1 <-
  glmer(
    cbind(final_count, trials - final_count) ~  stress + treatment + treatment:stress + (1 | block),
    data = v,
    family = binomial
  )

vm2 <- update(vm1 , . ~ . - treatment:stress)

vm1_emmeans <-
  emmeans(vm1, pairwise ~ treatment | stress, 
          type = 'response')

v_success <-
  merge(v_success, 
        data.frame(vm1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# success/survival output ------------------------------- # 

tab_model(vm1, pred.labels = pred_labels, file = 'output/v_success_glmer_table.html')

write_lrt_table(mod1 = vm1, vm2, outfile = 'output/v_success_lrt.html')

v_success %>% 
  rename(
    `proportion survived` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  prob,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print('html',
      'output/v_success_stats.html')


#### --------- Vulpia mass per plant ---------- # 
vb1 <-
  lmer(l_mass ~ stress * treatment + (1 | block),
       data = v,
       REML = FALSE)

vb2 <- update(vb1, . ~ . - treatment:stress)
vb1 <- update(vb1, REML = TRUE)

vb1_emmeans <- emmeans(vb1, pairwise ~ treatment | stress , 
                       type = 'response')

v_size <-
  merge(v_size, 
        data.frame(vb1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# size output -------- # 

vm1_table <- tab_model(vb1, pred.labels  = pred_labels,
                       file = 'output/v_size_lmer_table.html')


write_lrt_table(vb1, vb2, 'output/v_size_lrt.html')

v_size %>% 
  rename(
    `avg. size` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  emmean,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print(  'html', 'output/v_size_stats.html' ) 

#### ----------- Vulpia infls. per plant

# Doesn't converge -------------------------------------------------------- 
# vf1glmer <-
#   glmer(
#     infls ~ offset(log_final_count) + stress * treatment + (1 |
#                                                               block),
#     data = vfdata,
#     family = 'poisson'
#   )
# 
# summary(vf1glmer)

vf1_glm_pois <-
  glm(
    infls ~ offset(log_final_count) + stress * treatment,
    data = vfdata,
    family = 'poisson'
  )

summary(vf1_glm_pois)

vf1_glm_qpois <-
  glm(
    infls ~ offset(log_final_count) + stress * treatment,
    data = vfdata,
    family = 'quasipoisson'
  )

summary(vf1_glm_qpois)

vf2_glm_qpois <- update(vf1_glm_qpois, . ~ . - treatment:stress)

vf1_emmeans <- emmeans(vf1_glm_qpois, 
                   pairwise ~ treatment | stress, 
                   offset = 0 , 
                   type = 'response')

v_infls <-
  merge(v_infls, 
        data.frame(vf1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# write infls output --------------- # 

vf1_table <- tab_model(vf1_glm_qpois, pred.labels = pred_labels,
                       file = 'output/v_infls_glm_table.html')

data.frame( rbind( drop1( vf1_glm_qpois, test = 'Chisq'), 
                   drop1(vf2_glm_qpois, test = 'Chisq')[2:3, ])) %>% 
  mutate( `terms removed` = row.names(. ) , 
          Pr.Chi. = round(Pr..Chi., 3)) %>% 
  xtable %>% 
  print( 'html', 'output/v_infls_lrt.html' )

v_infls %>% 
  rename(
    `avg. infls. per plant` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  rate,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print(  'html', 'output/v_flowers_summary.html' ) 


#### ---- Bromus success/survival ------------------------------ # 

bm1 <-
  glmer(
    cbind(final_count, trials - final_count) ~  stress + treatment + treatment:stress + (1 | block),
    data = b,
    family = binomial
  )

bm2 <- update(bm1 , . ~ . - treatment:stress)

bm1_emmeans <-
  emmeans(bm1, pairwise ~ treatment | stress, 
          type = 'response')

b_success <-
  merge(b_success, 
        data.frame(bm1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# success/survival output ------------------------------- # 

tab_model(bm1, pred.labels = pred_labels, file = 'output/b_success_glmer_table.html')

write_lrt_table(mod1 = bm1, bm2, outfile = 'output/b_success_lrt.html')

b_success %>% 
  rename(
    `proportion survived` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  prob,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print('html',
        'output/b_success_stats.html')


#### --------- Bromus mass per plant ---------- # 
bb1 <-
  lmer(l_mass ~ stress * treatment + (1 | block),
       data = b,
       REML = FALSE)

bb2 <- update(bb1, . ~ . - treatment:stress)
bb1 <- update(bb1, REML = TRUE)

bb1_emmeans <- emmeans(bb1, pairwise ~ treatment | stress , 
                       type = 'response')

b_size <-
  merge(b_size, 
        data.frame(bb1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# size output -------- # 

bm1_table <- tab_model(bb1, pred.labels  = pred_labels,
                       file = 'output/b_size_lmer_table.html')


write_lrt_table(bb1, bb2, 'output/b_size_lrt.html')

b_size %>% 
  rename(
    `avg. size` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  emmean,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print(  'html', 'output/b_size_stats.html' ) 

#### ----------- bromus infls. per plant

# Doesn't converge -------------------------------------------------------- 
# bf1glmer <-
#   glmer(
#     infls ~ offset(log_final_count) + stress * treatment + (1 | block),
#     data = bfdata,
#     family = 'poisson'
#   )
# 
# summary(bf1glmer)

bf1_glm_pois <-
  glm(
    infls ~ offset(log_final_count) + stress * treatment,
    data = bfdata,
    family = 'poisson'
  )

summary(bf1_glm_pois)

bf1_glm_qpois <-
  glm(
    infls ~ offset(log_final_count) + stress * treatment,
    data = bfdata,
    family = 'quasipoisson'
  )

summary(bf1_glm_qpois)

bf2_glm_qpois <- update(bf1_glm_qpois, . ~ . - treatment:stress)

bf1_emmeans <- emmeans(bf1_glm_qpois, 
                       pairwise ~ treatment | stress, 
                       offset = 0 , 
                       type = 'response')

b_infls <-
  merge(b_infls, 
        data.frame(bf1_emmeans$emmeans), 
        by = c('treatment', 'stress'))

# write infls output --------------- # 

bf1_table <- tab_model(bf1_glm_qpois, pred.labels = pred_labels,
                       file = 'output/b_infls_glm_table.html')

data.frame( rbind( drop1( bf1_glm_qpois, test = 'Chisq'), 
                   drop1( bf2_glm_qpois, test = 'Chisq')[2:3, ])) %>% 
  mutate( `terms removed` = row.names(. ) , 
          Pr.Chi. = round(Pr..Chi., 3)) %>% 
  xtable %>% 
  print( 'html', 'output/b_infls_lrt.html' )

b_infls %>% 
  rename(
    `avg. infls. per plant` = result,
    sd = stand_dev,
    `S.E.` = stand_err,
    `fitted mean` =  rate,
    `Lower 95% C.L.` = asymp.LCL,
    `Upper 95% C.L.` = asymp.UCL
  ) %>% 
  arrange( stress, treatment ) %>% 
  xtable %>% 
  print(  'html', 'output/b_flowers_summary.html' ) 


