## Computer Code for "Effects of native bryophytes on exotic grass invasion across an environmental gradient"

This directory contains all R scripts necessary to run the simulations and analyses and generate the figures presented in the manuscript. 

### Reproducing analyses 

The 'code' directory contains all the R scripts. The best way to recreate the analyses in the manuscript is to open the moss_analysis.Rproj file in Rstudio.  This will automatically set the working directory appropriately.  Then run the "code/run_all_scripts.R" file. Alternatively, one can set the working directory manually to the "moss_analysis" folder containing the code, data, figures, and output folders.  Then run the "run_all_scripts.R" file. 

### File details 

1. run_all_scripts.R

  Runs all the data analysis and plotting scripts in the correct order. 

2. moss_theme.R

  plotting parameters for ggplot2 plots 
  
3. point_intercept_analysis.R

  Perform analyses and plot observational point intercept data.  

4. moss_removal_analysis.R

  Perform analysis of experimental moss removal data. 

5. moss_removal_experiment_figures.R
 
  Generate figure showing response of annual grasses to moss removal across stress gradient. 
  
6. print_association_model_tables.Rmd

  Generate summary statistics tables and summarize the statistcal tests for observational data. Saves output as a word document.  Requires that the previous scripts have been run. 

7. print_experiment_tables.Rmd
  
  Generate summary statistics tables and summarize the statistcal tests for experimental data. Saves output as a word document.  Requires that the previous scripts have been run. 

8. species_info.csv

  Comma seperated value file containing information on the vascular plant species in the observational data. 

9. vulpia_data.csv

   Comma seperated value file containing demographic data for Vulpia from moss removal experiment. 
  
10. brid_data.csv

  Comma seperated value file  containing demographic data for Bromus from moss removal experiment. 
  
11. moss_cover.csv

  Comma seperated value file containing observational point intercept data on moss and plant cooccurence across the gradient. 

### Built With 

platform       x86_64-apple-darwin15.6.0   
arch           x86_64                      
os             darwin15.6.0                
system         x86_64, darwin15.6.0        
status                                     
major          3                           
minor          5.3                         
year           2019                        
month          03                          
day            11                          
svn rev        76217                       
language       R                           
version.string R version 3.5.3 (2019-03-11)
nickname       Great Truth 

### Required R packages 

1. tidyverse_1.2.1
2. gridExtra_2.3   
3. lme4_1.1-19  
4. cowplot_0.9.3
5. emmeans_1.3.0
6. scales_1.0.0
7. stringr_1.3.1
8. ggplot2_3.1.0
9. rmarkdown_1.10

### Authors 

Andrew Kleinhesselink 
arklein@ucla.edu 





