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

pid <- read.csv('data/moss_cover.csv')

pid$cover_category[pid$cover_category == 'ercameria'] = 'ericameria'

cover_counts = t(table(pid$cover_category, pid$transect))

cover_counts <- cover_counts[ , -2]


cover_counts <- as.data.frame( cover_counts ) 
cover_counts <- cover_counts %>% spread( key = Var2, Freq) %>% rename( distance = Var1)

cover_counts <- cover_counts %>% mutate( shrubs = ericameria + lupine ) 

prop_cover <- cover_counts %>% dplyr::select( bare, ericameria, lupine, moss, shrubs )/25 

prop_cover$distance <- cover_counts$distance

coverLong <- cbind(distance = as.numeric( levels(prop_cover$distance)[ prop_cover$distance ]), stack(prop_cover[, 1:5]))
names(coverLong) <- c("distance", "percentCover", "Category")

pid$hit[pid$species != 0 ] = 1
pid$hit[pid$species == 0 ] = 0

no_shrub <- subset(pid, cover_category %in% c('bare', 'moss'))

no_shrub$cover_category<-factor(no_shrub$cover_category)
no_shrub$cover_category

no_shrub <- no_shrub[, -c(2, 4)]
head(no_shrub)
summary(no_shrub)

testTable <- table( no_shrub$cover_category, no_shrub$hit) [ -c(2:4), ]
testTable
chisq.test(testTable)

no_shrub$transect2 <- no_shrub$transect - median(range(no_shrub$transect)) ###### center distances at midpoint of range 

# model hits per patch type across the gradient ------- 
hit_mod.binomial <- glm(hit ~ transect*cover_category, data = no_shrub, family = 'binomial')
summary(hit_mod.binomial)
anova(hit_mod.binomial)

hit_mod1 <- glm(hit ~ transect*cover_category, data = no_shrub, family = 'quasibinomial')
summary(hit_mod1)
anova(hit_mod1, test= 'F')

drop1(hit_mod1, test = 'LRT') # likelihood ratio test 

2*pt(0.825, df = 397, lower.tail= FALSE) 

# add NA points to fill out prediction line in figure 
pred_frame <- rbind( hit_mod1$data, data.frame( transect = c(9, 23, 47, 66),   
                                                cover_category = 'moss', 
                                                hit = NA, 
                                                transect2 = NA))

predicted <- data.frame(predict(hit_mod1, newdata = pred_frame, se.fit = TRUE))
predicted$upperSE <- predicted[, 1] + predicted[, 2]
predicted$lowerSE <- predicted[, 1] - predicted[, 2]
predicted

hit_data <- cbind(pred_frame, predicted = expit(predicted[,c(1,4,5)]))
hit_data2 <- aggregate(hit ~ transect*cover_category, data = hit_data , FUN= 'mean', na.action = 'na.pass')

hit_data2$predicted  <- aggregate(predicted.fit ~ transect*cover_category, data = hit_data, FUN = 'mean')[,3]
hit_data2$upperSE <- aggregate(predicted.upperSE ~ transect*cover_category, data = hit_data, FUN = 'mean')[, 3]
hit_data2$lowerSE <- aggregate(predicted.lowerSE ~ transect*cover_category, data = hit_data, FUN = 'mean')[, 3]
head(hit_data2)

##### add counts of moss and bare to show the points weight 
cover_counts$transect <- row.names(cover_counts)
hit_data2$counts <- NA
tail(hit_data2)

hit_data2[ hit_data2$cover_category == 'bare', 'counts'] <- cover_counts$bare
hit_data2[ hit_data2$cover_category == 'moss', 'counts'] <- cover_counts$moss

levels(coverLong$Category) <- c("Bare sand", "ericameria", "lupine", "Moss patch", "Shrub patch")
coverLong$Category

#### plot cover over gradient ------------------------------ 
xTitle <- xlab_distance
yTitle1 <- "Proportion cover"
yTitle2 <- "Probability of plant rooted at point"

coverLong$Category <- factor( coverLong$Category,levels=c('Moss patch', 'Bare sand', 'Shrub patch', 'ericameria', 'lupine'))

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
  scale_color_manual(values = c('black', 'black', 'darkgray')) + 
  scale_y_continuous(labels = percent_format(accuracy = 1)) + 
  moss_theme + 
  scale_shape_manual(values = c(19, 2, 15)) + 
  theme(legend.position = 'right',
        legend.justification = 0, 
        legend.key.width = unit(3.5, 'line'),
        plot.margin = margin(c(10, 10, 10, 20)), 
        axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))), 
        axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 



p1_mod <- 
  ggdraw(p1 + ylim( 0, 1)) + 
  draw_text(c('Low\nStress', 'High\nStress'), 
            x = c(0.05, 0.75), 
            y = 0.15, size = 15)

coverWide <- coverLong %>% spread(Category, percentCover)
moss_cover <- lm(data =  coverWide , `Moss patch` ~ poly(distance, 2) )

sand_cover <- lm(data =  coverWide , `Bare sand` ~ poly(distance, 2) )

shrub_cover <- lm(data =  coverWide , `Shrub patch` ~ poly(distance, 2) )

summary( moss_cover ) 
summary(sand_cover)
summary(shrub_cover)

#### plot probability of plant being rooted in moss vs. bare ground 


hit_data2$cover_category
hit_data2$Category <- factor( hit_data2$cover_category , levels = c('moss', 'bare'), ordered = T)
hit_data2$cover_category
hit_data2$Category <- factor(hit_data2$Category, labels = c('Moss patch', 'Bare sand'))
hit_data2$Category


p2 <- ggplot( hit_data2, aes( x = transect, y = hit )) + 
  geom_ribbon(aes(ymin = lowerSE, ymax = upperSE, fill = Category, group = Category), alpha = 0.5) + 
  geom_line(aes(y = predicted, group = Category, linetype = Category), color = 'black', size = 0.8) + 
  geom_point( aes(shape = Category), size = 6) + 
  scale_shape_manual(values = c(19, 2)) + 
  scale_linetype_manual(values= c(1, 2) ) + 
  scale_fill_grey() +
  moss_theme + 
  theme(
    legend.key = element_rect(color = 'white'), 
    axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))),
    axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 

p2_legend <- get_legend(p2)


p3 <- 
  ggplot( hit_data2, aes( x = transect, y = hit )) + 
  geom_ribbon(aes(ymin = lowerSE, ymax = upperSE, fill = Category, group = Category), alpha = 0.5) + 
  geom_line(aes(y = predicted, group = Category, linetype = Category), color = 'black', size = 0.8) + 
  geom_point( aes(shape = Category, size = counts)) + 
  scale_shape_manual(values = c(19, 2)) + 
  scale_linetype_manual(values= c(1, 2) ) + 
  scale_fill_grey() +
  scale_size(range = c(1,7), guide = F) + 
  scale_x_continuous(name = xlab_distance) + 
  scale_y_continuous(name = 'Proportion of points with plants') + 
  moss_theme + 
  theme(
    legend.position = "none", 
    legend.key = element_rect(color = 'white'), 
    axis.title.y = element_text(margin = margin(c(0, 15, 0, 0))),
    axis.title.x = element_text(margin = margin(c(10, 0, 0, 0)))) 


p3 <- cowplot::plot_grid(p3 , p2_legend , ncol = 2, rel_widths = c(1, 0.3))


p3 #### for plot 


p3_mod <- 
  ggdraw(p3) + 
  draw_text(c('Low\nStress', 'High\nStress'), 
            x = c(0.06, 0.82), 
            y = 0.17, size = 15)
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

######## species specific analysis          #####################
no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

table(no_shrub$cover_category)
speciesTable <- table(no_shrub$species)

length(speciesTable)
#### total non-shrub points: 
sum(table(no_shrub$cover_category))
#### total non-shrub points with plants: 
sum(speciesTable[ -1])

barplot(speciesTable, horiz= TRUE, las = 2)

spTableByCat <- table(no_shrub$species, no_shrub$cover_category)

spTableByCat <- data.frame(spTableByCat)

spTableByCat <- spTableByCat[ spTableByCat$Var2 %in% c('moss', 'bare'),  ]

ggplot(spTableByCat, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar( stat = 'identity', position = 'dodge') + coord_flip()

######## vulpia hits in moss and bare sand sub-analysis #########
pid$vHit[ pid$species == 'vubr' ] <- 1
pid$vHit[ is.na(pid$vHit)  ] <- 0

no_shrub = subset(pid, cover_category %in% c('moss', 'bare'))

testTable <- table(no_shrub$cover_category, no_shrub$vHit)[ -c(2:4), ]
chisq.test(testTable)

vm1 <- glm(vHit ~ transect*cover_category, no_shrub, family = 'binomial')
summary(vm1)
vulpPredDF <- cbind( no_shrub, predictedVulpia = predict(vm1), se = predict(vm1, se.fit = TRUE)$se.fit )

vulpPredDF$upperSE <- expit(vulpPredDF$predictedVulpia + vulpPredDF$se)
vulpPredDF$lowerSE <- expit(vulpPredDF$predictedVulpia - vulpPredDF$se)
vulpPredDF$pred <- expit(vulpPredDF$predictedVulpia)

vulpPredDF$cover_category <- factor(vulpPredDF$cover_category, levels = c('moss', 'bare'))
levels( vulpPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

names(vulpPredDF )[ 3] <- 'Category'

vHitsAgg <- aggregate( vHit ~ transect*Category, vulpPredDF , FUN = 'mean')

vHitsAgg[ vHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
vHitsAgg[ vHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

vulpPlot <- ggplot(vulpPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of Vulpia rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = vHitsAgg, aes(y = vHit, size = counts))
             
ggsave(filename= 'figures/vulpHits.png',  plot = vulpPlot, height= 5, width = 8, units= 'in', dpi = 300 )

###################### Bromus hits in moss analysis ###############
pid$bHit[ pid$species == 'brdi' ] <- 1
pid$bHit[ is.na(pid$bHit)  ] <- 0

no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

vm1 <- glm(bHit ~ transect*cover_category, no_shrub, family = 'binomial')
summary(vm1)
brdiPredDF = cbind( no_shrub, predictedHits = predict(vm1), se = predict(vm1, se.fit = TRUE)$se.fit )

brdiPredDF$upperSE <- expit(brdiPredDF$predictedHits + brdiPredDF$se)
brdiPredDF$lowerSE <- expit(brdiPredDF$predictedHits - brdiPredDF$se)
brdiPredDF$pred <- expit(brdiPredDF$predictedHits)

brdiPredDF$cover_category <- factor(brdiPredDF$cover_category, levels = c('moss', 'bare'))
levels( brdiPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

names(brdiPredDF )[ 3] <- 'Category'

bHitsAgg = aggregate( bHit ~ transect*Category, brdiPredDF , FUN = 'mean')

bHitsAgg[ bHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
bHitsAgg[ bHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

brdiPlot = ggplot(brdiPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of Vulpia rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = bHitsAgg, aes(y = bHit, size = counts))

ggsave(filename= 'figures/brdiHits.png',  plot = brdiPlot, height= 5, width = 8, units= 'in', dpi = 300 )

###################### Annual grass hits in moss analysis ###############
pid$aHit[ pid$species %in% c('brdi', 'vubr') ] <- 1
pid$aHit[ is.na(pid$aHit)  ] <- 0

no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

testMat <- table( no_shrub$cover_category, no_shrub$aHit) [ -c(2:4), ]
class(testMat)
chisq.test( x= testMat)

head(no_shrub)

vm1 <- glm(aHit ~ transect*cover_category, no_shrub, family = 'binomial')
summary(vm1)
anova(vm1, test = 'Chisq')

agrassPredDF <- cbind( no_shrub, predictedHits = predict(vm1), se = predict(vm1, se.fit = TRUE)$se.fit )

agrassPredDF$upperSE <- expit(agrassPredDF$predictedHits + agrassPredDF$se)
agrassPredDF$lowerSE <- expit(agrassPredDF$predictedHits - agrassPredDF$se)
agrassPredDF$pred <- expit(agrassPredDF$predictedHits)

agrassPredDF$cover_category <- factor(agrassPredDF$cover_category, levels = c('moss', 'bare'))
levels( agrassPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

names(agrassPredDF )[ 3] <- 'Category'

aHitsAgg <- aggregate( aHit ~ transect*Category, agrassPredDF , FUN = 'mean')

aHitsAgg[ aHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
aHitsAgg[ aHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

agrassPlot <- ggplot(agrassPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of Vulpia rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = aHitsAgg, aes(y = aHit, size = counts))

ggsave(filename= 'figures/agrassHits.png',  plot = agrassPlot, height= 5, width = 8, units= 'in', dpi = 300 )

# hits by origin ---------------------------------------------------------------------------------------------------- 
d <- as.matrix( head( read.csv('data/moss_association_data.csv'), 1 ) )
origin <- colnames(d)
species <-  d[1, ]

origin_table <- data.frame(  origin, species ) %>% mutate( origin = str_extract(pattern = '[a-z]+', origin ))

pid <- merge( pid, origin_table , by = 'species', all.x = TRUE)

pid$eHit <- 0
pid$eHit[ pid$origin == 'exotic' ] <- 1

pid$nHit <- 0
pid$nHit[ pid$origin == 'native' ] <- 1

# exotic ------------------------------------------------------------------------------------------------------------

no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

testMat <- table( no_shrub$cover_category, no_shrub$eHit)[ -c(2:4), ]
class(testMat)
chisq.test( x= testMat)

em1 <- glm(eHit ~ transect*cover_category, no_shrub, family = 'quasibinomial')
summary(em1)
anova(em1, test = 'F')

ePredDF <- cbind( no_shrub, predictedHits = predict(em1), se = predict(em1, se.fit = TRUE)$se.fit )

ePredDF$upperSE <- expit(ePredDF$predictedHits + ePredDF$se)
ePredDF$lowerSE <- expit(ePredDF$predictedHits - ePredDF$se)
ePredDF$pred <- expit(ePredDF$predictedHits)

ePredDF$cover_category <- factor(ePredDF$cover_category, levels = c('moss', 'bare'))
levels(ePredDF$cover_category ) <- c('Moss patch', 'Bare sand')

ePredDF$Category <- ePredDF$cover_category

eHitsAgg <- aggregate( eHit ~ transect*Category, ePredDF , FUN = 'mean')

eHitsAgg[ eHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
eHitsAgg[ eHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

ePlot <- ggplot(ePredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of exotic species rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = eHitsAgg, aes(y = eHit, size = counts))

ggsave(filename= 'figures/eHits.png',  plot = ePlot, height= 5, width = 8, units= 'in', dpi = 300 )

# native ------------------------------------------------------------------------------------------------------------

no_shrub <- subset(pid, cover_category %in% c('moss', 'bare'))

testMat <- table( no_shrub$cover_category, no_shrub$nHit)[ -c(2:4), ]
class(testMat)
chisq.test( x= testMat)

nm1 <- glm(nHit ~ transect*cover_category, no_shrub, family = 'quasibinomial')
summary(nm1)
anova(nm1, test = 'Chisq')

nPredDF <- cbind( no_shrub, predictedHits = predict(nm1), se = predict(nm1, se.fit = TRUE)$se.fit )

nPredDF$upperSE = expit(nPredDF$predictedHits + nPredDF$se)
nPredDF$lowerSE = expit(nPredDF$predictedHits - nPredDF$se)
nPredDF$pred = expit(nPredDF$predictedHits)

nPredDF$cover_category <- factor(nPredDF$cover_category, levels = c('moss', 'bare'))
levels(nPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

nPredDF$Category <- nPredDF$cover_category

nHitsAgg <- aggregate( nHit ~ transect*Category, nPredDF , FUN = 'mean')

nHitsAgg[ nHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
nHitsAgg[ nHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

nPlot <- ggplot(nPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of exotic species rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = nHitsAgg, aes(y = nHit, size = counts))


ggsave(filename= 'figures/nHits.png',  plot = nPlot, height= 5, width = 8, units= 'in', dpi = 300 )


summary(hit_mod1)
anova(hit_mod1, test= 'F')

