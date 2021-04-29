# SDMcollembola
Repository with the data and the code used to run SDMs of invasive Antarctic collembola
## Purpose of this repository
The purpose of this repository is to share the code and files that have been used to carry out the analyses presented in: _Ensemble forecasting of the invasive risk in four exotic springtail species in Antarctica_ under review in _Polar Biology_ (POBI-D-20-00115)
## Requirements
All the analyses have been carried out using R.
### libraries used

#### modelling
biomod2
adehabitatHS
adehabitatMA

#### handling of spatial data
rgdal
raster
sp

#### data wrangling
scales
plyr
reshape2

#### plotting
ggplot2
gridExtra
grid
RColorBrewer
## Contained files
### scripts
`script_ENFA_and figures.R` contains the ENFA analysis and the code to create the figures.
`Biomod2_script.R` R code with the code ran to model the species distribution using two sets of datapoints: (1) only with the native observed distribution and (2) with the native and non-native distribution.
### BIOMOD2_inputs Datapoints
Species observations used in the modelling. 
