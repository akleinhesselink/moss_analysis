#### Moss Point Intercept Data

rm(list = ls())

library(coin)
library(mgcv)
library(sm)
library(ggplot2)
library(grid)

source('../moss_theme.R')

expit = function(v){
  p = exp(v)/(exp(v) + 1)
  return (p)
}

pid = read.csv('moss_cover.csv')

pid$cover_category[pid$cover_category == 'ercameria'] = 'ericameria'

cover_counts = t(table(pid$cover_category, pid$transect))

cover_counts = as.data.frame(cover_counts[ , -2])

cover_counts$shrubs = rowSums(cover_counts[, c(2:3)])

prop_cover = cover_counts/25
prop_cover$distance = row.names(cover_counts)

coverLong = cbind(distance = as.numeric(prop_cover$distance), stack(prop_cover[, 1:5]))
names(coverLong)  = c("distance", "percentCover", "Category")

#### plot cover over gradient 
p = ggplot(subset(coverLong, Category %in% c('moss', 'bare', 'shrubs')), 
           aes(x = distance, y = percentCover, color = Category, group = Category, shape = Category, linetype = Category)) 

p + geom_point(size = 3.5) + 
  geom_smooth( method = 'loess', se = FALSE, size = 1.25) + 
  xlab(xlab_distance) + 
  ylab("Proportion of points per transect")

barplot(colSums(cover_counts)  , ylab = "total hits")

pid$hit[pid$species != 0 ] = 1
pid$hit[pid$species == 0 ] = 0

no_shrub = subset(pid, cover_category %in% c('bare', 'moss'))

no_shrub$cover_category = factor(no_shrub$cover_category)
no_shrub$cover_category

no_shrub = no_shrub[, -c(2, 4)]
head(no_shrub)
summary(no_shrub)

###### center distances at midpoint of range 
no_shrub$transect2 = no_shrub$transect - median(range(no_shrub$transect)) 

hit_mod1 = glm(hit ~ transect2*cover_category, data = no_shrub, family = 'quasibinomial')
summary(hit_mod1)
anova(hit_mod1)

predicted = data.frame(predict(hit_mod1, se.fit = TRUE))
predicted$upperSE = predicted[, 1] + predicted[, 2]
predicted$lowerSE = predicted[, 1] - predicted[, 2]
predicted

hit_data = cbind(no_shrub, predicted = expit(predicted[,c(1,4,5)]))
head(hit_data)

hit_data2 = aggregate(hit ~ transect*cover_category, data = hit_data , FUN= 'mean')
hit_data2$predicted  = aggregate(predicted.fit ~ transect*cover_category, data = hit_data, FUN = 'mean')[,3]
hit_data2$upperSE = aggregate(predicted.upperSE ~ transect*cover_category, data = hit_data, FUN = 'mean')[, 3]
hit_data2$lowerSE = aggregate(predicted.lowerSE ~ transect*cover_category, data = hit_data, FUN = 'mean')[, 3]
head(hit_data2)

##### add counts of moss and bare to show the points weight 
cover_counts$transect = row.names(cover_counts)
hit_data2$counts = NA
tail(hit_data2)

hit_data2[ hit_data2$cover_category == 'bare', 'counts'] <- cover_counts$bare
hit_data2[ hit_data2$cover_category == 'moss', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

p = ggplot(hit_data2, aes(x = transect, y = hit, fill = cover_category, color = cover_category))
p + geom_point( aes(size = counts))

p + geom_point(aes(size = counts)) + 
  geom_line(aes(y = predicted)) + 
  geom_ribbon(aes(ymin=lowerSE,ymax=upperSE),alpha=0.3) + 
  xlab(xlab_distance) + ylab('Probability of plant')

levels(coverLong$Category) = c("Bare sand", "ericameria", "lupine", "Moss patch", "Shrub patch")
coverLong$Category

#### plot cover over gradient 
xTitle = xlab_distance
yTitle1 = "Proportion cover"
yTitle2 = "Probability of plant rooted at point"

coverLong$Category = factor( coverLong$Category,levels=c('Moss patch', 'Bare sand', 'Shrub patch', 'ericameria', 'lupine'))

p1 = ggplot(subset(coverLong, Category %in% c('Moss patch', 'Bare sand', 'Shrub patch')), 
            aes(x = distance, y = percentCover, color = Category, 
                group = Category, linetype = Category, shape= Category)) 

p1 = p1 + 
  geom_point(size = 3) +  
  geom_smooth( method = 'loess', se = FALSE, size = 1) + 
  xlab(xTitle) + 
  ylab(yTitle1) + 
  scale_color_grey()

p1

#### plot probability of plant being rooted in moss vs. bare ground 
names(hit_data2)[2] <- "Category"
levels(hit_data2$Category) <- c("Bare sand", "Moss patch")
hit_data2$Category <- factor(hit_data2$Category, levels = c('Moss patch', 'Bare sand'))

p2.base = ggplot(hit_data2, aes(x = transect, 
                           y = hit, 
                           color = Category, 
                           fill = Category, 
                           shape = Category, 
                           linetype = Category))

p2 = p2.base +  
  geom_line(aes(y = predicted), color = 'black', size = 0.8) + 
  ylab(yTitle2) + 
  xlab(xTitle) + 
  scale_color_grey(start= 0, end = 0.5) + 
  scale_linetype_manual(values= c(1, 2))  

p3 = p2 + geom_point(size = 3) + 
  geom_ribbon(aes(ymin = lowerSE, ymax = upperSE, color = NULL), alpha = 0.4) + 
  scale_fill_grey()  + 
  guides(fill = guide_legend(override.aes = list(colour = NULL))) 

p3

p2.withSize = p2 + 
  geom_ribbon(aes(ymin =lowerSE,
                  ymax = upperSE, color = NULL), alpha = 0.4) + 
  geom_point(aes(size = counts)) +   
  scale_size(range = c(2, 7), guide = 'none') +   
  scale_fill_grey() + 
  guides(fill = guide_legend(override.aes = list(colour = NULL)))

p2.withSize #### for plot 

dir.create('../figures' )
ggsave(filename='../figures/cover_on_transect.png', plot= p1, height= 5, width = 8, units= 'in', dpi = 300)
ggsave(filename= '../figures/rooted_in_moss.png', plot = p2.withSize, height = 5, width = 8, units = 'in', dpi = 300)
ggsave(filename= '../figures/rooted_in_moss_legend.png', plot = p3, height = 5, width = 8, units = 'in', dpi = 300 )

######
######


1- pchisq(54.3, df = 1) 

1- pchisq(17.8, df = 1)

1- pchisq(0.698, df =1 )

p_moss = expit(-1.648462 + 1.021037)
p_bare = expit(-1.648462)
p_moss
p_bare
no_shrub

points = aggregate(hit ~ transect + cover_category, data = no_shrub, length)
names(points)[3] = "points"
sum_hits = aggregate(hit ~ transect + cover_category, data = no_shrub, sum)

prop_hit = merge(points, sum_hits)
prop_hit$success_prob = prop_hit$hit/prop_hit$points
prop_hit$no_hit = prop_hit$points - prop_hit$hit 

p2 = ggplot(prop_hit, aes(x = transect, y = success_prob, color = cover_category)) 
p2 + geom_point() + geom_smooth(method = 'lm', se = FALSE) + theme_bw() +
  ylab(yTitle2) + xlab(xTitle) 

gam1 = with (prop_hit [ prop_hit$cover_category == 'bare', ], gam(cbind(hit, no_hit) ~  s(transect) ,family=binomial))
gam2 = with (prop_hit [prop_hit$cover_category == 'moss', ], gam(cbind(hit, no_hit) ~ s(transect), family = binomial))
summary(gam1)
summary(gam2)
anova(gam1)

plot(gam1)
plot(gam2)

response1 <- predict(gam1, type="response", se.fit=T)
response0 <- predict(gam2, type="response", se.fit=T)

par(mfcol=c(1,1))
plot(0, type="n", bty="n", main="Fancy GAM plot", xlab="MyCovariate", ylab="MyResponse", lwd=3,ylim=c(0,60), xlim=c(0,200))
legend("bottomright", bty="n", lwd=5, col=c("green","red"), legend=c("Strata = 0", "Strata = 1"))

#lines(spline(gam1$model$Covariate , response1$fit) , lwd = 3 , col = "red")
#lines(sm.spline(gam1$model$Covariate , response1$fit+1.96*response1$se) , lty = 3 , lwd = 2 , col = "red")
#lines(sm.spline(gam1$model$Covariate , response1$fit-1.96*response1$se) , lty = 3 , lwd = 2 , col = "red")

#lines(sm.spline(gam2$model$Covariate , response0$fit) , lwd = 3 , col = "green")
#lines(sm.spline(gam2$model$Covariate, response0$fit + 1.96 * response0$se) , lty = 3 , lwd = 2, col = "green")
#lines(sm.spline(gam2$model$Covariate, response0$fit - 1.96 * response0$se) , lty = 3 , lwd = 2 , col = "green")

######## species specific analysis          #####################
no_shrub = subset(pid, cover_category %in% c('moss', 'bare'))

table(no_shrub$cover_category)
speciesTable = table(no_shrub$species)
barplot(speciesTable, horiz= TRUE, las = 2)

spTableByCat = table(no_shrub$species, no_shrub$cover_category)

spTableByCat = data.frame(spTableByCat)

spTableByCat = spTableByCat[ spTableByCat$Var2 %in% c('moss', 'bare'),  ]

ggplot(spTableByCat, aes(x = Var1, y = Freq, fill = Var2)) + 
  geom_bar( stat = 'identity', position = 'dodge') + coord_flip()

######## vulpia hits in moss and bare sand sub-analysis #########
pid$vHit[ pid$species == 'vubr' ] <- 1
pid$vHit[ is.na(pid$vHit)  ] <- 0

no_shrub = subset(pid, cover_category %in% c('moss', 'bare'))

vm1 = glm(vHit ~ transect*cover_category, no_shrub, family = 'binomial')
summary(vm1)
vulpPredDF = cbind( no_shrub, predictedVulpia = predict(vm1), se = predict(vm1, se.fit = TRUE)$se.fit )

vulpPredDF$upperSE = expit(vulpPredDF$predictedVulpia + vulpPredDF$se)
vulpPredDF$lowerSE = expit(vulpPredDF$predictedVulpia - vulpPredDF$se)
vulpPredDF$pred = expit(vulpPredDF$predictedVulpia)

vulpPredDF$cover_category <- factor(vulpPredDF$cover_category, levels = c('moss', 'bare'))
levels( vulpPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

names(vulpPredDF )[ 3] <- 'Category'

vHitsAgg = aggregate( vHit ~ transect*Category, vulpPredDF , FUN = 'mean')

vHitsAgg[ vHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
vHitsAgg[ vHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

ggplot(vulpPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of Vulpia rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = vHitsAgg, aes(y = vHit, size = counts))
             

###################### Bromus hits in moss analysis ###############
pid$bHit[ pid$species == 'brdi' ] <- 1
pid$bHit[ is.na(pid$bHit)  ] <- 0

no_shrub = subset(pid, cover_category %in% c('moss', 'bare'))

vm1 = glm(bHit ~ transect*cover_category, no_shrub, family = 'binomial')
summary(vm1)
brdiPredDF = cbind( no_shrub, predictedVulpia = predict(vm1), se = predict(vm1, se.fit = TRUE)$se.fit )

brdiPredDF$upperSE = expit(brdiPredDF$predictedVulpia + brdiPredDF$se)
brdiPredDF$lowerSE = expit(brdiPredDF$predictedVulpia - brdiPredDF$se)
brdiPredDF$pred = expit(brdiPredDF$predictedVulpia)

brdiPredDF$cover_category <- factor(brdiPredDF$cover_category, levels = c('moss', 'bare'))
levels( brdiPredDF$cover_category ) <- c('Moss patch', 'Bare sand')

names(brdiPredDF )[ 3] <- 'Category'

bHitsAgg = aggregate( bHit ~ transect*Category, brdiPredDF , FUN = 'mean')

bHitsAgg[ bHitsAgg$Category == 'Bare sand', 'counts'] <- cover_counts$bare
bHitsAgg[ bHitsAgg$Category == 'Moss patch', 'counts'] <- cover_counts$moss[cover_counts$moss > 0]

ggplot(brdiPredDF , aes(x = transect, shape = Category, fill = Category, color = Category, y = pred) ) + 
  geom_line( aes(linetype = Category)) + 
  geom_ribbon (aes( ymin = upperSE, ymax = lowerSE), alpha = 0.2, color = NA) + 
  scale_color_grey() +
  scale_fill_grey() + 
  ylab('Probability of Vulpia rooted at point') + 
  xlab( xlab_distance) +   
  geom_point(data = bHitsAgg, aes(y = bHit, size = counts))


