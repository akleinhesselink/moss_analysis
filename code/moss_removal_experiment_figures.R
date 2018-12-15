source('data_and_analysis/moss_removal_experiment/moss_removal_analysis.R')

library(scales)
library(grid)
library(gridExtra)


#### plot parameters
xTitle = "Position on environmental gradient"
yTitle1 = 'Percent success\n(seeds surviving to mature plants)'
yTitle2 = 'Mean plant size (log g)\n'
yTitle3 = 'Inflorescences per plant (no.)\n'

theme_update( strip.text = element_text(face = 'italic', size = 15))

label_df <- data.frame(label = c('a)', 'b)', 'c)', 'd)', 'e)', 'f)'), species = c('Bromus', 'Vulpia'), x_pos = 0.6, y_pos = 1, Treatment = 'Moss patch')

fig1 = 
  ggplot(data = prop_success, aes( x = Stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = Treatment, 
                                         color = Treatment, shape = Treatment)) + 
  theme(axis.ticks.length = unit( 5, unit = 'pt'))
        
fig1 <-  
  fig1 + geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) +
  ylab(yTitle1) + 
  facet_grid( . ~ species) +   
  scale_color_grey() +
  xlab(xTitle) + 
  scale_y_continuous(labels = percent_format()) + 
  geom_text(data = label_df[1:2,], aes( x = x_pos, y = y_pos, label = label, ymin = 0, ymax = 0), size = 5, show.legend = FALSE)  


fig2 = ggplot(data = size, aes( x = Stress, y = result, ymin = result - stand_err, 
                                         ymax = result + stand_err, group = Treatment, 
                                         color = Treatment, shape = Treatment)) + 
  theme(axis.ticks.length = unit(5, unit = 'pt'))

fig2 = 
  fig2 + 
  geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab(yTitle2) + 
  facet_grid(. ~ species ) + 
  scale_color_grey()  + 
  xlab(xTitle) + 
  geom_text(data = label_df[3:4,], aes( x = x_pos, y = y_pos*max(size$result) + 0.3, label = label, ymin = 0, ymax = 0), size = 5, vjust = 0.1, show.legend = FALSE) 

fig3 = ggplot(data = infls, aes( x = Stress, y = result, ymin = result - stand_err, 
                                ymax = result + stand_err, group = Treatment, 
                                color = Treatment, shape = Treatment)) + 
  theme(axis.ticks.length = unit(5, unit = 'pt'))


fig3 = 
  fig3 + geom_point(position = position_dodge(width = 0.2), size = 4) + 
  geom_line(stat = 'identity', position = position_dodge( width = 0.2)) + 
  geom_errorbar (position = position_dodge(), width = 0.2) + 
  ylab(yTitle3) + 
  facet_grid(. ~ species )  + 
  scale_color_grey() +
  xlab(xTitle) + 
  geom_text(data = label_df[5:6,], aes( x = x_pos, y = y_pos*max(infls$result) + 0.2, label = label, ymin = 0, ymax = 0), size = 5, vjust = 0.1, show.legend = FALSE) 


fig3

fmt_dcimals <- function(decimals=0){
  function(x) format(x,nsmall = decimals,scientific = FALSE)
}

top_panel <- 
  fig1 + 
  scale_color_grey(guide = 'none') + 
  scale_shape_discrete(guide = 'none') + 
  theme( 
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    plot.margin = margin( c(5,150,1,20))) 
  
mid_panel <- 
  fig2 + 
  scale_y_continuous(labels = fmt_dcimals(1)) + 
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(), 
    strip.text = element_blank(), 
    plot.margin = margin(c(0, 0, 1, 34)))

bottom_panel <- 
  fig3 + 
  scale_shape_discrete(guide = 'none') + 
  scale_color_grey( guide = 'none')  + 
  theme(strip.text = element_blank(), 
        plot.margin = margin( c(0, 150, 15, 34)), 
        axis.title.x = element_text(margin=margin(15,0,0,0)) )  


plot_all <- 
  grid.arrange(arrangeGrob(top_panel, mid_panel, bottom_panel, ncol = 1, nrow = 3, heights = c(1, 0.9, 1.1) ))


ggsave(path = 'figures', filename = 'all.plots.png', plot = plot_all, width = 8, height = 10, units = 'in', dpi = 300)

# ggsave(path= 'figures', filename='survival.png', plot=fig1, width= 8, height = 5, units = 'in', dpi= 300)
# ggsave(path = 'figures', filename = 'biomass.png', plot= fig2, width = 8, height = 5, units = 'in', dpi = 300)
# ggsave(path = 'figures', filename = "infls.png", plot = fig3, width= 8, height = 5, units= "in", dpi = 300)


# # 
# 
# library("lme4")
# library(scales)
# library(dplyr)
# library(tidyr)
# library("ggplot2") # Plotting
# 
# newdat <- unique( cbind(trials = 5, final_count = 0, vm1@frame[, c('stress', 'treatment')]))
# 
# mm <- model.matrix(terms(vm1),newdat)
# 
# newdat$surv <- predict(vm1,newdat,re.form=NA, type = 'link')
# 
# pvar1 <- diag(mm %*% tcrossprod(vcov(vm1),mm))
# tvar1 <- pvar1+VarCorr(vm1)$block[1]  ## must be adapted for more complex models
# cmult <- 1.96 ## could use 1.96
# 
# newdat <- data.frame(
#   newdat
#   , ploSE = newdat$surv - sqrt(pvar1)
#   , phiSE = newdat$surv + sqrt(pvar1)
#   , plo = newdat$surv - cmult*sqrt(pvar1)
#   , phi = newdat$surv + cmult*sqrt(pvar1)
#   , tlo = newdat$surv - cmult*sqrt(tvar1)
#   , thi = newdat$surv + cmult*sqrt(tvar1)
# )
# 
# #plot confidence
# 
# head(newdat)
# newdat
# a
# 
# 
# plot_dat <- newdat %>% 
#   gather( type, val , surv, ploSE:thi) %>% 
#   mutate( pSurv = inv.logit( val ) ) %>%
#   select( stress, treatment, type,  pSurv) %>% 
#   spread( type, pSurv )
# 
# g0 <- ggplot(plot_dat, aes(x=stress, y=surv, group = treatment, colour=treatment, ymin = plo, ymax = phi))+
#   geom_point( position = position_dodge( width = 0.5)) + 
#   geom_errorbar( position = position_dodge( width = 0.5), color = 'gray', width = 0.2)+
#   geom_errorbar( aes( ymin  = ploSE, ymax = phiSE), position = position_dodge( width = 0.5), width = 0.1) +
#   geom_line( aes( group = treatment), position = position_dodge( width = 0.5), linetype = 'dashed') + 
#   labs(title="CI based on fixed-effects uncertainty ONLY") + 
#   ylab ( "Percent suvival") + 
#   scale_y_continuous(labels = percent_format())
# 
# 
# g0
# 
# 
