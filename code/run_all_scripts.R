# All analyses for "Effects of native bryophytes on exotic grass invasion: 
# a test of the stress gradient hypothesis" 
# 
# Instructions: To recreate analyses and figures, open the "moss_analysis.Rproj" file 
# in Rstudio and then run this script. Alternatively set the working directory to the 
# "moss_analysis" folder and then run this script. All R scripts must be in the code 
# folder and all the data files must be in the data folder. 

rm(list = ls())

library(tidyverse)
library(scales)
library(grid)
library(gridExtra)
library(cowplot)
library(lme4)
library(emmeans)

source( 'code/point_intercept_analysis.R' )

source( 'code/moss_removal_analysis.R' ) 

source( 'code/moss_removal_experiment_figures.R')

rmarkdown::render('code/print_association_model_tables.Rmd')

rmarkdown::render('code/print_experiment_tables.Rmd')