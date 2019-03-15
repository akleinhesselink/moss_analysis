rm(list = ls())

library(tidyverse)

sp_info <- read_csv('data/species_info.csv')
biomass <- readxl::read_xls('data/Bodega biomass 2010 v3.xls')

biomass <- 
  biomass %>% 
  gather( species, biomass , c(`total biomass`:aica)) 

tbio <- 
  biomass %>% 
  filter( species == 'total biomass')

tbio$Excluded <- factor(tbio$`Treat Combo`, labels = c('None', 'Deer', 'Rabbits', 'Deer & Rabbits'))

gg_biomass <- 
  tbio %>% 
  ggplot(aes( x = Distance, y = biomass , color = Excluded, shape = Excluded)) + 
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  ylab( 'Total biomass (g)')

ggsave(gg_biomass, filename = 'figures/total_biomass_on_gradient.png', width = 5, height = 3 )

mod1 <- lm(data = tbio, log( biomass) ~ Distance)
summary(mod1)  

drop1(mod1, test = 'Chisq')

mod2 <- lm(data = tbio, log(biomass) ~ Distance + Jackrabbits + Deer + Jackrabbits*Deer)
summary(mod2)

drop1(mod2, test = 'Chisq')


unique( biomass$species )
biomass %>% 
  group_by( Distance, `Treat Combo`) %>% 
  summarise( total = biomass[species == 'total biomass'] , 
             total2 = sum(biomass[species != 'total biomass'], na.rm = T), 
             mix =  sum(biomass[species == 'mixture'], na.rm = T)) %>% 
  mutate( nmix = total - mix)

biomass %>% 
  filter( species == 'mixture' ) %>% 
  ggplot(aes( x = Distance, y = biomass , color = `Treat Combo`)) + 
  geom_point() + 
  geom_line() + 
  geom_smooth(aes( group = 1), method = 'lm', se = F, col = 'gray', linetype = 2)
 


