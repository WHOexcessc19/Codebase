

# WHO COVID-19 Excess Mortality Estimates Code and Data 

This repository contains code for fitting of the excess mortality models developed for the 2022 WHO COVID-19 excess mortality estimates, and scripts for data pre-processing, intermediate model fitting, and analysis for the manuscripts, "Estimating country-specific excess mortality during the COVID-19 pandemic" and "Estimates of the excess mortality associated with the COVID-19 pandemic from the World Health Organization."

## Directories

* `Excess/`  contains R scripts for fitting the excess death models for pandemic times. It contains the individual models for all cause mortality (ACM) for the special subnational data and annual data countries, as well as Final_Sampling.R which combines everything in a unified sampling scheme. 
* `Expected/` contains R scripts for fitting the expected death models for pandemic times. 
* `Generated_Data/` contains data objects generated from the Excess and Expected models. 
* `Imported_Data/` contains imported data objects used in the models. 
* `PreData_Structures/` contains preprocessed data files and the R scripts for computing measures and formatting the data. 

## Dependencies

* R version 3.6.3 and the following R packages: "tidyverse", "posterior", "INLA", "zoo", "readxl", "cmdstanr", "gridExtra", "mgcv", "matrixStats".
  * IMPORTANT: The INLA package should be installed as the latest testing (not stable) branch, with instructions on the https://www.r-inla.org/download-install website
	* To install any missing packages, run `sapply(setdiff(c("tidyverse", "posterior", "zoo", "readxl", "cmdstanr", "gridExtra", "mgcv", "matrixStats"), installed.packages()), install.packages)`

## Workflow
Note: If only interested in the excess model fitting and not the data processing or expected death model, skip to the Excess section/code and the data for the pre-processing and expected components are pre-computed in the Generated_Data folder.

### Data Pre-Processing

### Expected
 
### Excess
* If not interested in running the annual data ACM models or subnational data ACM models, skip to the Final_Sampling.R file and the data for the the annual/subnational model components are pre-computed in the Generated_Data folder.
* Argentina_Estimates.R: Data processing and fitting of the subnational model for 2021 ACM in Argentina, using ACM data from the Cordoba region.
* China_Estimates.R: Data processing and fitting of the annual data model for monthly ACM in China, using annual ACM data for 2020 and 2021.
* India_Estimates.R: Data processing and fitting of the subnational model for ACM in India, using ACM data from the 17 states in India.
* Indonesia_Estimates.R: Data processing and fitting of the subnational model for ACM in Indonesia, using ACM data from Jakarta.
* Turkey_Estimates.R: Data processing and fitting of the subnational model for ACM in Turkey, using ACM data from a sample of Turkey.
* AnnualData_Estimates.R: Data processing and fitting of the annual data model for 2020 monthly ACM in Viet Nam, Grenada, Saint Kitts and Nevis, and Saint Vincent and the Grenadines, using annual ACM data for 2020.
* Final_Sampling.R: The final combined sampling of the expecteds and ACM from the observed data country approach, subnational data and annual data country approaches, and the no ACM data covariate model approach.
	* Warning: This may take a while. 

## Contact

Maintainer contact information: William Msemburi (msemburiw@who.int); Victoria Knutson (vknuts@uw.edu); Serge Aleshin-Guendel (aleshing@uw.edu)

## Package citations

  R Core Team (2020). R: A language and environment for statistical
  computing. R Foundation for Statistical Computing, Vienna, Austria.
  URL https://www.R-project.org/.

  H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
  Springer-Verlag New York, 2016.

  Baptiste Auguie (2017). gridExtra: Miscellaneous Functions for "Grid"
  Graphics. R package version 2.3.
  https://CRAN.R-project.org/package=gridExtra
  
  H. Rue, S. Martino, and N. Chopin. Approximate Bayesian inference for latent Gaussian models
  using integrated nested Laplace approximations (with discussion). Journal of the Royal      
  Statistical Society, Series B, 71(2):319{392, 2009. 
  
  Stan Development Team (2022). “RStan: the R interface to Stan.” R package version 2.21. 5,
  https://mc-stan.org/.

