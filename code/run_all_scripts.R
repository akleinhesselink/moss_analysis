# All analyses for "Effects of native bryophytes on exotic grass 
# invasion across an environmental gradient" 
# 
# Note:  To run all the scripts for the analysis set the working directory to 
# the 'moss_manuscript' folder then run this script. All scripts and data files
# must be present. 

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

rmarkdown::render('code/point_intercept_tables.Rmd')

rmarkdown::render('code/experiment_tables.Rmd')