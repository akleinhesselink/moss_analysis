#### Moss Point Intercept Analysis

rm(list = ls())

library(tidyverse)
library(scales)
library(grid)
library(cowplot)

source('code/moss_theme.R')

expit = function(v){
  p = exp(v)/(exp(v) + 1)
  return (p)
}

gen_agg_hits_df <- function( my_mod ) { 
  # add NA points to fill out prediction line in figure 
  
  agg_data <- 
    my_mod$data %>% 
    group_by( transect, cover_category) %>% 
    summarise( plant_points = sum(hit), total_points = n()) 
  
  agg_data <- 
    expand.grid( transect = unique( my_mod$data$transect),  
                 cover_category = c('bare', 'moss')) %>% 
    left_join(agg_data) %>% 
    mutate( total_points = ifelse( is.na( total_points), 0, total_points ) )
  
  predicted <- data.frame(predict(my_mod, newdata = agg_data, se.fit = TRUE))
  
  predicted$upperSE <- predicted[, 1] + predicted[, 2]
  predicted$lowerSE <- predicted[, 1] - predicted[, 2]
  
  hit_data <- cbind(agg_data, expit(predicted[,c(1,4,5)]))
  
  hit_data <- 
    hit_data %>% 
    mutate( hit = plant_points/total_points) %>% 
    mutate( Category = factor( cover_category , levels = c('moss', 'bare'), ordered = T)) %>% 
    mutate( Category = factor( Category, labels = c('Moss patch', 'Bare sand')))
  
  return( hit_data ) 
}


pid <- read.csv('data/moss_cover.csv')

pid$cover_category[pid$cover_category == 'ercameria'] <- 'ericameria'
cover_counts <- t(table(pid$cover_category, pid$transect))
cover_counts <- cover_counts[ , -2]
cover_counts <- as.data.frame( cover_counts ) 
cover_counts <- 
  cover_counts %>% 
  spread( key = Var2, Freq) %>% 
  rename( distance = Var1) %>% 
  mutate( shrubs = ericameria + lupine ) 

prop_cover <- 
  cover_counts %>% 
  dplyr::select( bare, ericameria, lupine, moss, shrubs )/25 

prop_cover$distance <- cover_counts$distance

#### plot cover over gradient ------------------------------ 
coverLong <- 
  prop_cover %>% 
  gather( Category, percentCover, -distance ) %>% 
  mutate( distance = as.numeric(levels(distance)[distance])) %>% 
  mutate( Category = factor( Category, 
                             labels = c("Bare sand", "ericameria", "lupine", "Moss patch", "Shrub patch")))

xTitle <- xlab_distance
yTitle1 <- "Proportion cover"
yTitle2 <- "Probability of plant rooted at point"

p1 <- ggplot(subset(coverLong, Category %in% c('Moss patch', 'Bare sand', 'Shrub patch')), 
             aes(x = distance, 
                 y = percentCover, 
                 color = Category, 
                 group = Category, 
                 linetype = Category, 
                 shape= Category)) + 
  geom_point(size = 3) +  
  geom_smooth( method = 'lm',formula = y ~ poly(x, 2), se = FALSE, size = 0.8) + 
  xlab(xTitle) + 
  ylab(yTitle1) + 
  ylim( 0, 1) + 
  scale_color_manual(values = c('black', 'black', 'black')) + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) + 
  moss_theme + 
  scale_shape_manual(values = c(19, 2, 0)) + 
  theme(legend.position = 'right',
        legend.justification = 0, 
        legend.key.width = unit(3.5, 'line'),
        plot.margin = margin(c(10, 10, 10, 20)), 
        axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))), 
        axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 

p1_mod <- 
  ggdraw(p1 + ylim( 0, 1))
# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.05, 0.75), 
#             y = 0.15, size = 15)

coverWide <- 
  coverLong %>% 
  spread(Category, percentCover)

moss_cover_lm <- lm(data = coverWide, `Moss patch` ~ poly(distance, 2) )

# model hits per patch type across the gradient ------- 

pid$hit[pid$species != 0 ] = 1
pid$hit[pid$species == 0 ] = 0

no_shrub <- subset(pid, cover_category %in% c('bare', 'moss'))

no_shrub$cover_category<-factor(no_shrub$cover_category)
no_shrub <- no_shrub[, -c(2, 4)]

hit_mod.binomial <- glm(hit ~ transect*cover_category, 
                        data = no_shrub, family = 'binomial')

summary(hit_mod.binomial)

all_hits_glm <- glm(hit ~ transect*cover_category, 
                   data = no_shrub, family = 'quasibinomial')

anova(all_hits_glm, test = 'F')
summary(all_hits_glm)


hit_data <- gen_agg_hits_df(all_hits_glm)


hit_data %>% 
  ungroup() %>%
  group_by( Category ) %>% 
  summarise( sum( plant_points, na.rm= T), sum(total_points))

View(hit_data %>% 
  ungroup() %>% 
  group_by(transect) %>% 
  summarise( sum(total_points)))

#### plot probability of plant being rooted in moss vs. bare ground 

p2 <-  
  hit_data %>%
  ggplot( aes( x = transect, y = hit )) + 
  geom_ribbon(aes(ymin = lowerSE, 
                  ymax = upperSE, 
                  fill = Category, 
                  group = Category), alpha = 0.5) + 
  geom_line(aes(y = fit, 
                linetype = Category), color = 'black', size = 0.8)   + 
  geom_point( aes(shape = Category), size = 5) + 
  scale_shape_manual(values = c(19, 2)) + 
  scale_linetype_manual(values= c(2, 1) ) + 
  scale_fill_grey() +
  moss_theme + 
  theme(
    legend.key = element_rect(color = 'white'), 
    axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))),
    axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 

  
p2_legend <- get_legend(p2)

p3 <- 
  hit_data %>% 
  ggplot( aes( x = transect, y = hit )) + 
  geom_ribbon(aes(ymin = lowerSE, 
                  ymax = upperSE, 
                  fill = Category, 
                  group = Category), alpha = 0.5) + 
  geom_line(aes(y = fit, 
                linetype = Category), color = 'black', size = 0.8)   + 
  geom_point( aes(shape = Category, size = total_points)) + 
  scale_shape_manual(values = c(19, 2)) + 
  scale_linetype_manual(values= c(1, 2) ) + 
  scale_fill_grey() +
  scale_size(range = c(1,7)) + 
  scale_x_continuous(name = xlab_distance) + 
  scale_y_continuous(name = 'Proportion of points with plants' ) + 
  moss_theme + 
  theme(
    legend.position = "none", 
    axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))),
    axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 


p3_comb <- cowplot::plot_grid(p3, 
                              p2_legend, ncol = 2, rel_widths = c(1, 0.3))

p3_mod <- 
  ggdraw(p3_comb) 
# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.06, 0.82), 
#             y = 0.17, size = 15)

p3_mod

### Print plots out ----------------------------- # 

ggsave(filename='figures/cover_on_transect.png', 
       plot= p1_mod, 
       height= 4.5, 
       width = 7, 
       units= 'in', 
       dpi = 300)

ggsave(filename= 'figures/rooted_in_moss.png', 
       plot = p3_mod, 
       height = 4.5, 
       width = 7, 
       units = 'in', 
       dpi = 300)


######
# Additional species specific analyses of co-occurence with moss
######

points <- aggregate(hit ~ transect + cover_category, data = no_shrub, length)
names(points)[3] <- "points"
sum_hits <- aggregate(hit ~ transect + cover_category, data = no_shrub, sum)

prop_hit <- merge(points, sum_hits)
prop_hit$success_prob <- prop_hit$hit/prop_hit$points
prop_hit$no_hit <- prop_hit$points - prop_hit$hit 

prop_hit
######## species specific analysis          #####################
no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

table(no_shrub$cover_category)
speciesTable <- table(no_shrub$species)

#### total non-shrub points: 
sum(table(no_shrub$cover_category))
#### total non-shrub points with plants: 
sum(speciesTable[ -1])

spTableByCat <- table(no_shrub$species, no_shrub$cover_category)

spTableByCat <- data.frame(spTableByCat)

spTableByCat <- spTableByCat[ spTableByCat$Var2 %in% c('moss', 'bare'),  ]

ggplot(spTableByCat, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar( stat = 'identity', position = 'dodge') + coord_flip()

spTableByCat %>% filter( Var1 %in% c('vubr', 'brdi'))

######## vulpia hits in moss and bare sand sub-analysis #########

vdata <- 
  pid %>% 
  filter( cover_category %in% c('moss', 'bare')) %>% 
  mutate( hit = ifelse( species == 'vubr', 1, 0) ) %>% 
  dplyr::select(-species )

vm1 <- glm(hit ~ transect*cover_category, vdata, family = 'quasibinomial')

summary(vm1)
anova(vm1, test = 'F')

hit_data <- gen_agg_hits_df(vm1)

hit_plot <- p3 %+% hit_data + ggtitle('Vulpia')

hit_plot <- cowplot::plot_grid(hit_plot, 
                               p2_legend, ncol = 2, rel_widths = c(1, 0.3))

hit_plot <- 
  ggdraw(hit_plot) 
# + 
  # draw_text(c('Low\nStress', 'High\nStress'), 
  #           x = c(0.06, 0.82), 
  #           y = 0.17, size = 15)

ggsave(filename= 'figures/vulpHits.png',  
       plot = hit_plot, height= 5, width = 8, units= 'in', dpi = 300 )

###################### Bromus hits in moss analysis ###############

bdata <- 
  pid %>% 
  filter( cover_category %in% c('moss', 'bare')) %>% 
  mutate( hit = ifelse( species == 'brdi', 1, 0) ) %>% 
  dplyr::select(-species )

bm1 <- glm(hit ~ transect*cover_category, bdata, family = 'quasibinomial')
summary(bm1)
anova(bm1, test = 'F')

hit_data <- gen_agg_hits_df(bm1)

hit_plot <- p3 %+% hit_data + ggtitle('Bromus')

hit_plot <- cowplot::plot_grid(hit_plot, p2_legend , ncol = 2, rel_widths = c(1, 0.3))

hit_plot <- 
  ggdraw(hit_plot) 

# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.06, 0.82), 
#             y = 0.17, size = 15)

ggsave(filename= 'figures/brdiHits.png',  
       plot = hit_plot, height= 5, width = 8, units= 'in', dpi = 300 )

###################### Annual grass hits in moss analysis ###############
agdata <- 
  pid %>% 
  filter( cover_category %in% c('moss', 'bare')) %>% 
  mutate( hit = ifelse( species %in% c('brdi', 'vubr'), 1, 0) ) %>% 
  dplyr::select(-species )

agm1 <- glm(hit ~ transect*cover_category, agdata, family = 'quasibinomial')
summary(agm1)
anova(agm1, test = 'F')

hit_data <- gen_agg_hits_df(agm1)

hit_plot <- p3 %+% hit_data + ggtitle('Annual Grass')

hit_plot <- cowplot::plot_grid(hit_plot, 
                               p2_legend, ncol = 2, rel_widths = c(1, 0.3))

hit_plot <- 
  ggdraw(hit_plot) 
# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.06, 0.82), 
#             y = 0.17, size = 15)

ggsave(filename= 'figures/agrassHits.png',  
       plot = hit_plot, height= 5, width = 8, units= 'in', dpi = 300 )

# hits by origin ---------------------------------------------------------------------------------------------------- 
origin_table <-  read.csv('data/species_info.csv')
origin_table$origin <- as.character( origin_table$origin )

pid2 <- left_join(pid, origin_table, by = 'species')

# exotic ------------------------------------------------------------------------------------------------------------
edata <- 
  pid2 %>% 
  filter( cover_category %in% c('moss', 'bare')) %>%
  mutate( origin =  ifelse(is.na(origin), 0, origin)) %>% 
  mutate( hit = ifelse( origin == 'exotic', 1, 0) ) %>% 
  dplyr::select(-species, -origin ) 

em1 <- glm(hit ~ transect*cover_category, edata, family = 'quasibinomial')
summary(em1)
anova(em1, test = 'F')

hit_data <- gen_agg_hits_df(em1)

hit_plot <- 
  p3 %+% hit_data + 
  ggtitle('Exotic Species')

hit_plot <- cowplot::plot_grid(hit_plot, 
                               p2_legend, ncol = 2, rel_widths = c(1, 0.3))

hit_plot <- 
  ggdraw(hit_plot) 
# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.06, 0.82), 
#             y = 0.17, size = 15)


ggsave(filename= 'figures/eHits.png',  
       plot = hit_plot, height= 5, width = 8, units= 'in', dpi = 300 )

# native ------------------------------------------------------------------------------------------------------------
ndata <- 
  pid2 %>% 
  filter( cover_category %in% c('moss', 'bare')) %>%
  mutate( origin =  ifelse(is.na(origin), 0, origin)) %>% 
  mutate( hit = ifelse( origin == 'native', 1, 0) ) %>% 
  dplyr::select(-species, -origin ) 

nm1 <- glm(hit ~ transect*cover_category, ndata, family = 'quasibinomial')
summary(nm1)
anova(nm1, test = 'F')

hit_data <- gen_agg_hits_df(nm1)

hit_plot <- p3 %+% 
  hit_data + 
  ggtitle('Native Species')

hit_plot <- cowplot::plot_grid(hit_plot, 
                               p2_legend, ncol = 2, rel_widths = c(1, 0.3))

hit_plot <- 
  ggdraw(hit_plot) 
# + 
#   draw_text(c('Low\nStress', 'High\nStress'), 
#             x = c(0.06, 0.82), 
#             y = 0.17, size = 15)

ggsave(filename= 'figures/nHits.png',  
       plot = hit_plot, height= 5, width = 8, units= 'in', dpi = 300 )


save(moss_cover_lm, all_hits_glm, vm1, bm1, agm1, em1, nm1, 
     file = 'output/point_intercept_models.rda')
