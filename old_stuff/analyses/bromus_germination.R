#### Bromus Germination plotting

#### 
remove(list = ls())
setwd("~/Documents/moss_manuscript/moss_data")
library(MASS)

source('moss_analyses/moss_analysis_tools.R')
b = read.csv('bromus/bromus_longitudinal.csv')
b_wide = read.csv('bromus/bromus_wide_corrected.csv')


head(b)

b = b[ order(b$days, b$block, b$treatment, b$days), ]
head(b)

b$plot = 1:54
head(b)
b_plot_wide = reshape(b, v.names = 'count', idvar = 'plot', timevar = 'days', direction = 'wide')
head(b_plot_wide)


for (i in 1:nrow(b_plot_wide)){
  print(paste(i, counter(b_plot_wide[i, 5:14])))
}

#write.table(b_plot_wide, "bromus/b_wide.csv", sep = ',', row.names = FALSE)
b_wide_corrected = read.csv('bromus/bromus_wide_corrected.csv')


for (i in 1:nrow(b_wide_corrected)){
  print(paste(i, counter(b_wide_corrected[i, 5:14])))
}

#### Bromus plots

b_wide_corrected

days = c(0,7,14, 21, 28, 36,43,60,77,130)

bromus_long = reshape(b_wide_corrected, idvar = 'plot', varying = list(5:14), v.names = 'count', times = days, direction = 'long')
head(bromus_long)
max(bromus_long$time)


b_ag = aggregate(count ~ treatment + time + stress, data = bromus_long, FUN = 'mean')
head(b_ag)
b_low = reshape(subset(b_ag, stress == 'low'), idvar =  'time', timevar = 'treatment', v.names = 'count',  direction = 'wide')
b_high = reshape(subset(b_ag, stress == 'high'), idvar =  'time', timevar = 'treatment', v.names = 'count',  direction = 'wide')

b_wide = rbind(b_low, b_high)

head(b_wide)

emergence_over_time(b_wide, b_low, b_high, rows = c(3:5))

bromus_mass = read.csv('bromus/brdi_data.csv')

head(bromus_mass)

#### Check that number of plants in the final count are the same in the biomass data
b_check = merge(bromus_long[bromus_long$time == 130, ], bromus_mass, INDEX = c('block', 'treatment') )
head(b_check)
b_check$final_count == b_check$count