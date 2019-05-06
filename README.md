## Computer Code for "Effects of native bryophytes on exotic grass invasion across an environmental gradient" a manuscript accepted for publication in *Ecosphere*. 

### Reproducing the analyses: 

Open the moss_analysis.Rproj file in the program Rstudio.  This will automatically set the working directory appropriately. Then run "run_all_scripts.R" from within Rstudio. Alternatively, one can set the working directory manually to the "moss_analysis" folder containing the code, data, figures, and output folders.  Then run the "run_all_scripts.R" file. 

When using this data please cite Kleinhesselink, Andrew, R. and J. Hall Cushman. 2019. "Effects of native bryophytes on exotic grass invasion across an environmental gradient". Ecosphere. 

### Files 


#### 1. data/species_info.csv

> Comma seperated value file containing information on the vascular plant species in the observational data. Species names and origin information taken from the Jepson Flora: http://ucjeps.berkeley.edu/eflora/

> #### Data Fields: 
- origin:  native or exotic.  Classify each species as either being native to the Bodega Bay dune community or non-native (exotic). 
- species: four letter abbreviation used for identifying species.  The four letters correspond to the first two letters of the genus and the first two letters of the species.  
- full_name:  latin binomial taxon name. 

#### 2. data/vulpia_data.csv

>Comma seperated value file containing demographic data for Vulpia from moss removal experiment. Each row represents data from one experimental patch. 

> #### Data Fields: 
- block: block group 1 - 18
- position: position on the environmental gradient, "SE" or "NW"
- treatment: moss patch treatment, "MC" = Moss Covered, "MR" = Moss Removed, "BC" = Bare Sand
- final_count: final number of surviving plants in each experimental patch
- infls: total number of inflorescences produced by all surviving plants in each patch. 
- total_mass_mg: total aboveground oven dried biomass (mg) of all surviving plants in each experimental patch. 
- mass_per_plant_mg: average mass per surviving plant, (total_mass_mg/final_count) 
  
#### 3. data/brid_data.csv

>Comma seperated value file  containing demographic data for Bromus from moss removal experiment. Each row represents data from one experimental patch. 

> #### Data Fields: 
- block: block group 1 - 18
- position: position on the environmental gradient, "SE" or "NW"
- treatment: moss patch treatment, "MC" = Moss Covered, "MR" = Moss Removed, "BC" = Bare Sand
- final_count: final number of surviving plants in each experimental patch
- infls: total number of inflorescences produced by all surviving plants in each patch. 
- total_mass_mg: total aboveground oven dried biomass (mg) of all surviving plants in each experimental patch. 
- mass_per_plant_mg: average mass per surviving plant, (total_mass_mg/final_count) 
  
#### 4. data/moss_cover.csv

>Comma seperated value file containing observational point intercept data on moss and plant cooccurence across the gradient. 

> #### Data Fields: 
- transect: position in (m) towards the NW on the environmental gradient
- point: sampling point along transect, 1 - 25
- cover_category: category for cover at sampled point: "bare", "moss", "ericameria" or "lupine".  Points marked lupine and Ericameria fell within these dominant shrub species at this site. 
- species: four letter species code for vascular plant intersected at sampled point. "0" if no plant was sampled at that point. 

#### 5. run_all_scripts.R

> Runs all the data analysis and plotting scripts in the correct order.

#### 6. code/moss_theme.R

> Plotting parameters for ggplot2 plots 
  
#### 7. code/point_intercept_analysis.R

> Perform analyses and plot observational point intercept data.  

#### 8. code/moss_removal_analysis.R

> Perform analysis of experimental moss removal data. 

#### 9. code/moss_removal_experiment_figures.R
 
> Generate figure showing response of annual grasses to moss removal across stress gradient. 
  
#### 10. code/point_intercept_tables.Rmd

> Generate summary statistics tables and summarize the statistcal tests for observational data. Saves output as a PDF document.  Requires that the previous scripts have been run. 

#### 11. code/experiment_tables.Rmd
  
> Generate summary statistics tables and summarize the statistcal tests for experimental data. Saves output as a PDF document.  Requires that the previous scripts have been run. 


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
9. rmarkdown_1.10
10. pander_0.6.3
11. broom_0.5.2

### Authors 

Andrew Kleinhesselink 
arklein@ucla.edu 





