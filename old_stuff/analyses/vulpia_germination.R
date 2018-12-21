#### Vulpia Germination Analysis 

#### 
remove(list = ls())
setwd("~/Documents/moss_manuscript/moss_data")
library(MASS)

source('moss_analyses/moss_analysis_tools.R')
v_wide = read.csv('vulpia/vulpia_wide.csv')
v_wide_corrected = read.csv('vulpia/vulpia_wide_corrected.csv')
head(v_wide)

for (i in 1:nrow(v_wide)){
  print(paste(i, counter(v_wide[i, 4:14])))
}

for (i in 1:nrow(v_wide_corrected)){
  print(paste(i, counter(v_wide_corrected[i, 4:14])))
}


v_wide_corrected$plot = 1:54

days = c(16,23, 30, 37,44, 51,59,66,83,100, 159)

vulpia_long = reshape(v_wide_corrected, idvar = 'plot', varying = list(4:14), v.names = 'count', times = days, direction = 'long')
vulpia_long

v = aggregate(count ~ time + treatment + stress + block, data = vulpia_long, FUN = 'sum')

v_ag = aggregate(count ~ treatment + time + stress, data = v, FUN = 'mean')

v_low = reshape(subset(v_ag, stress == 'Low'), idvar =  'time', timevar = 'treatment', v.names = 'count',  direction = 'wide')
v_high = reshape(subset(v_ag, stress == 'High'), idvar =  'time', timevar = 'treatment', v.names = 'count',  direction = 'wide')

v_wide = rbind(v_low, v_high)
head(v_wide)

emergence_over_time(v_wide, v_low, v_high)

#### Check that number of plants in the final count are the same in the biomass data

v_mass = read.csv("vulpia/vulpia_data.csv")

v_final = vulpia_long[ vulpia_long[,5] == 159, c(1,3,4,5,6)]

v_check = merge(v_mass, v_final, INDEX = c('block', 'treatment'))

v_check$final_count == v_check$count