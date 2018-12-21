# All analyses for "Moss controls exotic annual grass invasion: 
# a test of the stress gradient hypothesis" 

rm(list = ls())

source( 'code/point_intercept_analysis.R' )

source( 'code/moss_removal_analysis.R' ) 

source( 'code/moss_removal_experiment_figures.R')

rmarkdown::render('code/print_association_model_tables.Rmd')

rmarkdown::render('code/print_experiment_tables.Rmd')