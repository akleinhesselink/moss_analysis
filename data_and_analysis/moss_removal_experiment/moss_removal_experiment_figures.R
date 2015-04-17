remove(list = ls())

library(car)
library(boot)
library(ggplot2)
library(grid)

source( '../moss_theme.R')

xlab_distance = 'Position towards NW on environmental gradient (m)'

gen_results = function(df, formu, n, sp_name){ 
  res = aggregate(formu, data = df, 'mean')  
  names(res) [ 3] <- "result" 
  res$stand_dev = aggregate(formu, data = df, 'sd')[, 3]
  res$stand_err = res$stand_dev/sqrt(n)
  res$species = sp_name
  return(res)
} 

v = read.csv("vulpia/vulpia_data.csv")
b = read.csv('bromus/brdi_data.csv')

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
#v.aov = aov(logit_ps ~ stress*treatment + Error(block), data = v)
#summary(v.aov)
#1- pf(8.566, 2, 47)

N = table(v[, c(2,3)])
f_success = formula("prop_success ~ treatment + stress") 
f_size = formula("l_mass  ~ treatment + stress")
f_infls = formula("infls_per_plant ~ treatment + stress")

v_success = gen_results(v, f_success, 9, "Vulpia")
b_success = gen_results(b, f_success, 9, "Bromus")    
prop_success = rbind(v_success, b_success)
names(prop_success)[1:2] <- c('Treatment', 'Stress')

#write.table(prop_success, "prop_success.csv", sep = ',', row.names = FALSE)

v_size = gen_results(v, f_size, 9, "Vulpia")
b_size = gen_results(b, f_size, 9, "Bromus")
size = rbind(v_size, b_size)
names(size)[1:2] <- c('Treatment', 'Stress')

#write.table(size, "plant_size.csv", sep = ',', row.names = FALSE)

v_infls = gen_results(v, f_infls, 9, "Vulpia")
b_infls = gen_results(b, f_infls, 9, "Bromus")
infls = rbind(v_infls, b_infls)
names(infls)[1:2] <- c('Treatment', 'Stress')

#write.table(infls, "plant_infls.csv", sep = ',', row.names = FALSE)

#### plot parameters
xTitle = "Position on environmental gradient"
yTitle1 = 'Probability of seed germinating\nand surviving to mature plant'
yTitle2 = 'Mean plant size (log g)'
yTitle3 = 'Inflorescences per plant (no.)'

theme_update( strip.text = element_text(face = 'italic', size = 15))

fig1 = ggplot(data = prop_success, aes( x = Stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = Treatment, 
                                         color = Treatment, shape = Treatment))

fig1 = fig1 + geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab(yTitle1) + 
  facet_grid( . ~ species) +   
  scale_color_grey() +
  xlab(xTitle) 

fig1

fig2 = ggplot(data = size, aes( x = Stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = Treatment, 
                                         color = Treatment, shape = Treatment))

fig2 = fig2 + geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab(yTitle2) + 
  facet_grid(. ~ species ) + 
  scale_color_grey()  + 
  xlab(xTitle) 

fig2

fig3 = ggplot(data = infls, aes( x = Stress, y = result, ymin = result - stand_err, 
                                ymax = result + stand_err, group = Treatment, 
                                color = Treatment, shape = Treatment))
fig3 = fig3 + geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab(yTitle3) + 
  facet_grid(. ~ species )  + 
  scale_color_grey() +
  xlab(xTitle) 

fig3

dir.create(path = '../figures')
ggsave(path= '../figures', filename='survival.png', plot=fig1, width= 8, height = 5, units = 'in', dpi= 300)
ggsave(path = '../figures', filename = 'biomass.png', plot= fig2, width = 8, height = 5, units = 'in', dpi = 300)
ggsave(path = '../figures', filename = "infls.png", plot = fig3, width= 8, height = 5, units= "in", dpi = 300)
