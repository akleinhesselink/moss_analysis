#### vulpia and bromus final biomass and fecundity
#### 
remove(list = ls())
setwd("~/Documents/moss_manuscript/moss_data")
library(MASS)
library(nlme)
library(multcomp)
library(lsmeans)


source('moss_analyses/moss_analysis_tools.R')
v = read.csv("vulpia/vulpia_data.csv")
b = read.csv('bromus/brdi_data.csv')

v$infls_per_plant = v$infls/v$final_count
b$infls_per_plant = b$infls/b$final_count

#### biomass per plant analysis 
hist(v$mass_per_plant_mg)
hist(b$mass_per_plant_mg)
hist(log(v$mass_per_plant_mg))
hist(log(b$mass_per_plant_mg))

v_clean = v[ !is.na(v$mass_per_plant_mg), ]
b_clean = b[ !is.na(b$mass_per_plant_mg), ]

v_clean[, names(v_clean) == 'mass_per_plant_mg']


treatment_by_stress_plots(df = v_clean, y_var='mass_per_plant_mg', lg = T)

vmm1 = lme(log(mass_per_plant_mg) ~ stress*treatment, random = ~1|block,  data = v_clean)
summary(vmm1)
anova(vmm1)
vlsmeans = lsmeans(vmm1, ~ stress*treatment)

bmm1 = lme(log(mass_per_plant_mg) ~ stress*treatment, random = ~1|block,  data = b_clean)
summary(bmm1)
anova(bmm1)
lsmeans(bmm1, ~ stress*treatment)

contrast.matrix <- rbind("stresshigh:BC vs. stresshigh:MC" =  c(0,  0,  1,  0,  0,  0),
                         "stresshigh:BC vs. stresshigh:MR" =  c(0,  0,  0,  1,  0,  0),
                         "stresshigh:MC vs. stresshigh:MR" =  c(0,  0, -1,  1,  0,  0),
                                         
                         "stresslow:BC vs. stresslow:MC" =    c(0, -1,  0,  0,  1,  0), 
                         "stresslow:BC vs. stresslow:MR" =    c(0, -1,  0,  0,  0,  1), 
                         "stresslow:MC vs. stresslow:MR" =    c(0,  0,  0,  0, -1,  1))

summary(glht(vmm1, linfct =  mcp(treatment ="Tukey")))
summary(glht(vmm1, contrast.matrix))

summary(glht(bmm1, linfct = mcp(treatment = "Tukey")))
summary(glht(bmm1, contrast.matrix))

#### reproductive output analysis 
treatment_by_stress_plots(df = v_clean, y_var = 'infls_per_plant', lg = FALSE)
hist(v_clean$infls_per_plant)

vipp2 = lme(infls_per_plant ~ stress*treatment, random = ~ 1|block, data = v_clean)
summary(vipp2)
plot(vipp2)

bipp2 = lme(infls_per_plant ~ stress*treatment, random = ~ 1|block, data = b_clean)
summary(bipp2)
plot(bipp2)


#### number successful
treatment_by_stress_plots(df = v, y_var = 'final_count', lg = FALSE)
treatment_by_stress_plots(df = b, y_var = 'final_count', lg = FALSE)

v_prop_success = cbind(v$final_count, 5-v$final_count)

v_prop_success

vfc = glmmPQL(v_prop_success ~ stress*treatment, random = ~1|block, data = v, family = 'binomial')
summary(vfc)
plot(vfc)

b_prop_success = cbind(b$final_count, 5-b$final_count)
bfc = glmmPQL(b_prop_success ~ stress*treatment, random = ~1|block, family = 'binomial', data = b)
summary(bfc)
plot(bfc)

