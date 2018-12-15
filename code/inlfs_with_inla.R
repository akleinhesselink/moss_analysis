# check moss removal infls effects with inla 
remove(list = ls())
library(dplyr)
library(tidyr)
library(INLA)

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

b$infls

# --------------------------------------------------------------------------------------------------------------------- 
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

no_block <- 
  vfdata %>%
  mutate(block = NA, log_final_count = 0, final_count = 1, infls = NA)  %>% 
  distinct()

all_blocks <- 
  vfdata %>% 
  mutate( log_final_count = 0, final_count = 1, infls = NA) %>% 
  distinct()

vfdata_inla <- rbind(vfdata, no_block)
formula <- infls ~ offset( log_final_count) + treatment + stress + treatment:stress + f(block, model = "iid")

link <- rep(NA, length(vfdata_inla$infls))
link[which(is.na(vfdata$infls))] <- 1

result_fixed <- inla( formula , 
                      family = 'nbinomial', 
                      data = vfdata_inla , 
                      control.family = list( initial = 6, fixed = TRUE ),
                      control.predictor = list(compute = TRUE, link = 1))


result_fixed$summary.fixed
result_fixed$summary.linear.predictor
result_fixed$summary.fitted.values

linear_predictors <- vfdata_inla[is.na(vfdata_inla$infls), ] 
linear_predictors <- cbind(linear_predictors, result_fixed$summary.linear.predictor[which(is.na(vfdata_inla$infls)), ] )

ggplot( linear_predictors, aes( x = stress, y = exp(`0.5quant`) , color = treatment, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`) )) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar( position = position_dodge( width = 0.5))


result <- inla( formula , 
                family = 'nbinomial', 
                data = vfdata_inla , 
                control.predictor = list(compute = TRUE, link = 1))

result_fixed$summary.fixed
result$summary.fixed
result$summary.linear.predictor
result$summary.fitted.values

linear_predictors <- vfdata_inla[is.na(vfdata_inla$infls), ] 
linear_predictors <- cbind(linear_predictors, result$summary.linear.predictor[which(is.na(vfdata_inla$infls)), ] )

ggplot( linear_predictors, aes( x = stress, y = exp(`0.5quant`) , color = treatment, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`) )) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar( position = position_dodge( width = 0.5))


# ------------------------------------------------------------------------------------------------------------------------------
bfdata = b[ !b$final_count == 0, ]

bfdata <- bfdata %>% select( stress, treatment, block , infls, final_count ) 

sum( bfdata$infls == 0 )

bfdata %>% group_by ( block ) %>% summarise( n() )

bfdata$block <- factor( bfdata$block)

bfdata$log_final_count <- log(bfdata$final_count)

bfdata_summary <-  bfdata %>% 
  group_by ( stress, treatment ) %>% 
  mutate( ipp = infls / final_count ) %>%
  summarise( m_infls = mean(infls), m_final_count = mean(final_count), m_infls_per_plant = mean(ipp))

no_block <- 
  bfdata %>%
  mutate(block = NA, log_final_count = 0, final_count = 1, infls = NA)  %>% 
  distinct()

all_blocks <- 
  bfdata %>% 
  mutate( log_final_count = 0, final_count = 1, infls = NA) %>% 
  distinct()

bfdata_inla <- rbind(bfdata, no_block)
formula <- infls ~ offset( log_final_count) + treatment + stress + treatment:stress + f(block, model = "iid")

link <- rep(NA, length(bfdata_inla$infls))
link[which(is.na(bfdata$infls))] <- 1

result_fixed <- inla( formula , 
                      family = 'nbinomial', 
                      data = bfdata_inla , 
                      control.family = list( initial = 6, fixed = TRUE ),
                      control.predictor = list(compute = TRUE, link = 1))

result_fixed$summary.fixed
result_fixed$summary.linear.predictor
result_fixed$summary.fitted.values

linear_predictors <- bfdata_inla[is.na(bfdata_inla$infls), ] 
linear_predictors <- cbind(linear_predictors, result_fixed$summary.linear.predictor[which(is.na(bfdata_inla$infls)), ] )

ggplot( linear_predictors, aes( x = stress, y = exp(`0.5quant`) , color = treatment, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`) )) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar( position = position_dodge( width = 0.5))


result <- inla( formula , 
                family = 'nbinomial', 
                data = bfdata_inla , 
                control.predictor = list(compute = TRUE, link = 1))

result_fixed$summary.fixed
result$summary.fixed
result$summary.linear.predictor
result$summary.fitted.values

linear_predictors <- bfdata_inla[is.na(bfdata_inla$infls), ] 
linear_predictors <- cbind(linear_predictors, result$summary.linear.predictor[which(is.na(bfdata_inla$infls)), ] )

ggplot( linear_predictors, aes( x = stress, y = exp(`0.5quant`) , color = treatment, ymin = exp(`0.025quant`), ymax = exp(`0.975quant`) )) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar( position = position_dodge( width = 0.5))

