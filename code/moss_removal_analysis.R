remove(list = ls())

library(tidyverse)
library(lme4)
library(emmeans)

pred_labels <-
  c(
    'Intercept', 
    'NW',
    'Bare sand',
    'Moss removed',
    'NW:Bare sand',
    'NW:Moss removed'
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


v <- read.csv("data/vulpia_data.csv")
b <- read.csv('data/brdi_data.csv')

v$trials <- 5
b$trials <- 5

v$position <- factor(v$position, levels = c('SE', 'NW'))
b$position <- factor(b$position, levels = c('SE', 'NW'))

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
f_success <- formula("prop_success ~ treatment + position")
f_size <- formula("l_mass  ~ treatment + position")
f_infls <- formula("infls_per_plant ~ treatment + position")

v_success <- gen_results(v, f_success, 9, "Vulpia")

b_success <- gen_results(b, f_success, 9, "Bromus")
prop_success <- rbind(v_success, b_success)
names(prop_success)[1:2] <- c('Treatment', 'Position')

# Final biomass data ------------------------ # 

v_size <- gen_results(v, f_size, 9, "Vulpia")
b_size <- gen_results(b, f_size, 9, "Bromus")
size <- rbind(v_size, b_size)
names(size)[1:2] <- c('Treatment', 'Position')

# Final inflorescence numbers ----------------- # 
vfdata <- 
  v[!v$final_count == 0,] %>% 
  dplyr::select(position, treatment, block , infls, final_count) %>% 
  mutate( block = factor(block), 
          log_final_count = log(final_count), 
          infls_per_plant = infls/final_count)


bfdata <- 
  b[!b$final_count == 0,] %>% 
  dplyr::select(position, treatment, block , infls, final_count) %>% 
  mutate( block = factor(block), 
          log_final_count = log(final_count) ,
          infls_per_plant = infls/final_count)


v_infls <- gen_results(vfdata, f_infls, 9, "Vulpia")
b_infls <- gen_results(bfdata, f_infls, 9, "Bromus")
infls <- rbind(v_infls, b_infls)
names(infls)[1:2] <- c('Treatment', 'Position')

#### --- Vulpia final count -------------------------------------------

vm1 <-
  glmer(
    cbind(final_count, trials - final_count) ~  position + treatment + treatment:position + (1 | block),
    data = v,
    family = binomial
  )

vm2 <- update(vm1 , . ~ . - treatment:position)

vm1_emmeans <-
  emmeans(vm1, pairwise ~ treatment | position, 
          type = 'response')

v_success <-
  merge(v_success, 
        data.frame(vm1_emmeans$emmeans), 
        by = c('treatment', 'position'))

position_effect <- glm(
  cbind(final_count, trials - final_count ) ~ position, data = v[v$treatment == 'Bare sand', ], 
  family = binomial
  )

summary( position_effect )
drop1(position_effect, test = 'Chisq')

#### --------- Vulpia mass per plant ---------- # 
vb1 <-
  lmer(l_mass ~ position * treatment + (1 | block),
       data = v,
       REML = FALSE)

vb2 <- update(vb1, . ~ . - treatment:position)
vb1 <- update(vb1, REML = TRUE)

vb1_emmeans <- emmeans(vb1, pairwise ~ treatment | position , 
                       type = 'response')

v_size <-
  merge(v_size, 
        data.frame(vb1_emmeans$emmeans), 
        by = c('treatment', 'position'))


position_effect <- lm(
  l_mass ~ position, data = v[v$treatment == 'Bare sand', ]
)

summary( position_effect )
drop1(position_effect, test = 'Chisq')

#### ----------- Vulpia infls. per plant

# Doesn't converge -------------------------------------------------------- 
# vf1glmer <-
#   glmer(
#     infls ~ offset(log_final_count) + position * treatment + (1 |
#                                                               block),
#     data = vfdata,
#     family = 'poisson'
#   )
# 
# summary(vf1glmer)

vf1_glm_pois <-
  glm(
    infls ~ offset(log_final_count) + position * treatment,
    data = vfdata,
    family = 'poisson'
  )

summary(vf1_glm_pois)

vf1_glm_qpois <-
  glm(
    infls ~ offset(log_final_count) + position * treatment,
    data = vfdata,
    family = 'quasipoisson'
  )

summary(vf1_glm_qpois)

vf2_glm_qpois <- update(vf1_glm_qpois, . ~ . - treatment:position)

vf1_emmeans <- emmeans(vf1_glm_qpois, 
                   pairwise ~ treatment | position, 
                   offset = 0 , 
                   type = 'response')

v_infls <-
  merge(v_infls, 
        data.frame(vf1_emmeans$emmeans), 
        by = c('treatment', 'position'))


position_effect <- glm(
  infls ~ offset(log_final_count) + position, data = vfdata[vfdata$treatment == 'Bare sand', ], family = 'quasipoisson'
)
summary( position_effect )
drop1(position_effect, test = 'Chisq')

#### ---- Bromus success/survival ------------------------------ # 

bm1 <-
  glmer(
    cbind(final_count, trials - final_count) ~  position + treatment + treatment:position + (1 | block),
    data = b,
    family = binomial
  )


bm2 <- update(bm1 , . ~ . - treatment:position)

summary(bm1 )
drop1(bm1, test = 'Chisq')

bm1_emmeans <-
  emmeans(bm1, pairwise ~ treatment | position, 
          type = 'response')



b_success <-
  merge(b_success, 
        data.frame(bm1_emmeans$emmeans), 
        by = c('treatment', 'position'))


position_effect <- glm(
  cbind(final_count, trials - final_count ) ~ position, data = b[b$treatment == 'Bare sand', ], 
  family = binomial
)

summary( position_effect )

drop1(position_effect, test = 'Chisq')

#### --------- Bromus mass per plant ---------- # 
bb1 <-
  lmer(l_mass ~ position * treatment + (1 | block),
       data = b,
       REML = FALSE)

bb2 <- update(bb1, . ~ . - treatment:position)
bb1 <- update(bb1, REML = TRUE)

bb1_emmeans <- emmeans(bb1, pairwise ~ treatment | position , 
                       type = 'response')

emmeans(bb1, pairwise ~ position | treatment , 
                       type = 'response')

bb2_emmeans <- emmeans(bb2, pairwise  ~ treatment, type = 'response')
bb2_emmeans


b_size <-
  merge(b_size, 
        data.frame(bb1_emmeans$emmeans), 
        by = c('treatment', 'position'))


position_effect <- lm(
  l_mass ~ position, data = b[b$treatment == 'Bare sand', ]
)

summary( position_effect )
drop1( position_effect , test = 'Chisq')

#### ----------- Bromus infls. per plant

# Doesn't converge -------------------------------------------------------- 
# bf1glmer <-
#   glmer(
#     infls ~ offset(log_final_count) + position * treatment + (1 | block),
#     data = bfdata,
#     family = 'poisson'
#   )
# 
# summary(bf1glmer)

bf1_glm_pois <-
  glm(
    infls ~ offset(log_final_count) + position * treatment,
    data = bfdata,
    family = 'poisson'
  )

summary(bf1_glm_pois)

bf1_glm_qpois <-
  glm(
    infls ~ offset(log_final_count) + position * treatment,
    data = bfdata,
    family = 'quasipoisson'
  )

summary(bf1_glm_qpois)

bf2_glm_qpois <- update(bf1_glm_qpois, . ~ . - treatment:position)

bf1_emmeans <- emmeans(bf1_glm_qpois, 
                       pairwise ~ treatment | position, 
                       offset = 0 , 
                       type = 'response')

b_infls <-
  merge(b_infls, 
        data.frame(bf1_emmeans$emmeans), 
        by = c('treatment', 'position'))

position_effect <- glm(
  infls ~ offset( log_final_count) + position, data = bfdata[bfdata$treatment == 'Bare sand', ]
)

summary( position_effect )

drop1( position_effect, test = 'Chisq')

# Save processed data for figures -------------------------- #

save(prop_success, size, infls, 
     file = 'output/processed_experiment_data.rda')


# Save models for generation of statistical tables --------- # 

save(vm1, 
     vb1, 
     vf1_glm_qpois, 
     bm1, 
     bb1, 
     bf1_glm_qpois, 
     v, 
     b,
     v_infls,
     b_infls,
     vfdata, 
     bfdata, 
     file = 'output/experiment_models.rda')

