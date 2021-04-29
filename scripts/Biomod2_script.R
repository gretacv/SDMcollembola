library(biomod2)
library(rgdal)
require(raster)

# setwd("~/Desktop/biomod2/output_28052016")
setwd("~/Desktop/biomod2/output_29052016_native")

#bios_names<-c('bio1_Vmin_buff.tif', 'bio10_Vmin_buff.tif', 'bio11_Vmin_buff.tif', 'bio12_Vmin_buff.tif', 'bio13_Vmin_buff.tif', 'bio14_Vmin_buff.tif', 'bio15_Vmin_buff.tif', 'bio16_Vmin_buff.tif', 'bio17_Vmin_buff.tif', 'bio18_Vmin_buff.tif', 'bio19_Vmin_buff.tif', 'bio2_Vmin_buff.tif', 'bio3_Vmin_buff.tif', 'bio4_Vmin_buff.tif', 'bio5_Vmin_buff.tif', 'bio6_Vmin_buff.tif', 'bio7_Vmin_buff.tif', 'bio8_Vmin_buff.tif', 'bio9_Vmin_buff.tif')

vars <- c('bio1_Vmin_buff.tif', 'bio10_Vmin_buff.tif', 'bio11_Vmin_buff.tif', 'bio12_Vmin_buff.tif', 'bio13_Vmin_buff.tif', 'bio14_Vmin_buff.tif', 'bio15_Vmin_buff.tif', 'bio16_Vmin_buff.tif', 'bio17_Vmin_buff.tif', 'bio18_Vmin_buff.tif', 'bio19_Vmin_buff.tif', 'bio2_Vmin_buff.tif', 'bio3_Vmin_buff.tif', 'bio4_Vmin_buff.tif', 'bio5_Vmin_buff.tif', 'bio6_Vmin_buff.tif', 'bio7_Vmin_buff.tif', 'bio8_Vmin_buff.tif', 'bio9_Vmin_buff.tif') 


myExpl <- stack(paste("~/Desktop/biomod2/environmentalBIOVmin_buff_clip", vars, sep="/"))


# define the species of interest

sp.data.hv<- read.csv("~/Desktop/biomod2/datapoints/hypogastrura_viatica_01_less0_native.csv", header=TRUE)
sp.data.fc<- read.csv("~/Desktop/biomod2/datapoints/folsomia_candida_01_less0_native.csv", header=TRUE)
sp.data.mm<- read.csv("~/Desktop/biomod2/datapoints/mesaphorura_macrochaeta_01_less0_native.csv", header=TRUE)
sp.data.pm<- read.csv("~/Desktop/biomod2/datapoints/proisotoma_minuta_01_less0_native.csv", header=TRUE)

# Hypogastrura viatica bio2, bio4, bio5, bio15

  myRespName = "Hypogastrura_viatica"
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  
  myResp <- as.numeric(sp.data.hv[,3])
  
  myRespCoord = sp.data.hv[c("POINT_X","POINT_Y")]  
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = subset(myExpl,c(7,12,14,15)),   #we choose the bios we want to use, we have to be careful because the bios are not in order
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
									   PA.nb.rep=0, 
									   PA.nb.absences=NULL)
  
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  #Defining the quadratic interaction for the four variables in the GLM
  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list( myFormula = formula("Hypogastrura.viatica ~ bio15_Vmin_buff + bio15_Vmin_buff^2 + bio2_Vmin_buff + bio2_Vmin_buff^2 + bio4_Vmin_buff + bio4_Vmin_buff^2 + bio5_Vmin_buff + bio5_Vmin_buff^2")))

  myBiomodOption
  
  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE', 'CTA', 'RF','GLM', 'GAM','ANN','GBM','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=NULL, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC','KAPPA'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = TRUE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
   
   cat('\n',myRespName,'modeling finished')  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
  
  
 
  
  ### Make projections on current variable
 # myBiomodProjTSS <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(7,12,14,15)),   #variables have to be changed
 #   proj.name = 'current_TSS',
 #   selected.models = 'all',
 #   binary.meth = 'TSS',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	### Make projections on current variable
 # myBiomodProjKAPPA <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(7,12,14,15)), #variables have to be changed
 #   proj.name = 'current_KAPPA',
 #   selected.models = 'all',
 #   binary.meth = 'KAPPA',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	## Make projections on current variable
  myBiomodProjROC <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = subset(myExpl,c(7,12,14,15)),#variables have to be changed
    proj.name = 'current_ROC',
    selected.models = 'all',
    binary.meth = 'ROC',
    compress = 'xz',
    clamping.mask = TRUE,
    output.format = '.grd')
  

# Folsomia candida, bio2, bio4, bio5, bio17

  myRespName = "Folsomia_candida"
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  
  myResp <- as.numeric(sp.data.fc[,3])
  
  myRespCoord = sp.data.fc[c("POINT_X","POINT_Y")]  
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = subset(myExpl,c(9,12,14,15)),#variables have to be changed
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
									   PA.nb.rep=0, 
									   PA.nb.absences=NULL)
  
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
   #Defining the quadratic interaction for the four variables in the GLM
  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list( myFormula = formula("Folsomia.candida ~ bio17_Vmin_buff + bio17_Vmin_buff^2 + bio2_Vmin_buff + bio2_Vmin_buff^2 + bio4_Vmin_buff + bio4_Vmin_buff^2 + bio5_Vmin_buff + bio5_Vmin_buff^2")))

myBiomodOption

  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE', 'CTA', 'RF','GLM', 'GAM','ANN','GBM','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=NULL, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC','KAPPA'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = TRUE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
   
   cat('\n',myRespName,'modeling finished')  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
  
  
  
  ### Make projections on current variable
 # myBiomodProjTSS <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),#variables have to be changed
 #   proj.name = 'current_TSS',
 #   selected.models = 'all',
 #   binary.meth = 'TSS',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	### Make projections on current variable
 # myBiomodProjKAPPA <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),#variables have to be changed
 #   proj.name = 'current_KAPPA',
 #   selected.models = 'all',
 #   binary.meth = 'KAPPA',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	## Make projections on current variable
  myBiomodProjROC <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = subset(myExpl,c(9,12,14,15)),#variables have to be changed
    proj.name = 'current_ROC',
    selected.models = 'all',
    binary.meth = 'ROC',
    compress = 'xz',
    clamping.mask = TRUE,
    output.format = '.grd')
  

	
	 
# Mesaphorura macrochaeta, bio2, bio4, bio5, bio17

  myRespName = "Mesaphorura_macrochaeta"
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  
  myResp <- as.numeric(sp.data.mm[,3])
  
  myRespCoord = sp.data.mm[c("POINT_X","POINT_Y")]  
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = subset(myExpl,c(9,12,14,15)),
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
									   PA.nb.rep=0, 
									   PA.nb.absences=NULL)
  
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  #Defining the quadratic interaction for the four variables in the GLM
  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list( myFormula = formula("Mesaphorura.macrochaeta ~ bio17_Vmin_buff + bio17_Vmin_buff^2 + bio2_Vmin_buff + bio2_Vmin_buff^2 + bio4_Vmin_buff + bio4_Vmin_buff^2 + bio5_Vmin_buff + bio5_Vmin_buff^2")))

myBiomodOption

  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE', 'CTA', 'RF','GLM', 'GAM','ANN','GBM','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=NULL, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC','KAPPA'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = TRUE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
   
   cat('\n',myRespName,'modeling finished')  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
  
  
  
  
  ### Make projections on current variable
 # myBiomodProjTSS <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),
 #   proj.name = 'current_TSS',
 #   selected.models = 'all',
 #   binary.meth = 'TSS',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	### Make projections on current variable
 # myBiomodProjKAPPA <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),
 #   proj.name = 'current_KAPPA',
 #   selected.models = 'all',
 #   binary.meth = 'KAPPA',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	## Make projections on current variable
  myBiomodProjROC <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = subset(myExpl,c(9,12,14,15)),
    proj.name = 'current_ROC',
    selected.models = 'all',
    binary.meth = 'ROC',
    compress = 'xz',
    clamping.mask = TRUE,
    output.format = '.grd')
  
 
# Proisotoma minuta, bio2, bio4, bio5, bio17

  myRespName = "Proisotoma_minuta"
  
  cat('\n',myRespName,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  
  myResp <- as.numeric(sp.data.pm[,3])
  
  myRespCoord = sp.data.pm[c("POINT_X","POINT_Y")]  
  
  
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = subset(myExpl,c(9,12,14,15)),
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
									   PA.nb.rep=0, 
									   PA.nb.absences=NULL)
  
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  #Defining the quadratic interaction for the four variables in the GLM
  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list( myFormula = formula("Proisotoma.minuta ~ bio17_Vmin_buff + bio17_Vmin_buff^2 + bio2_Vmin_buff + bio2_Vmin_buff^2 + bio4_Vmin_buff + bio4_Vmin_buff^2 + bio5_Vmin_buff + bio5_Vmin_buff^2")))

myBiomodOption


  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
    myBiomodData, 
    models = c('SRE', 'CTA', 'RF','GLM', 'GAM','ANN','GBM','MARS','FDA'), 
    models.options = myBiomodOption, 
    NbRunEval=3, 
    DataSplit=80, 
    Prevalence=NULL, 
    VarImport=3,
    models.eval.meth = c('TSS','ROC','KAPPA'),
    SaveObj = TRUE,
    rescal.all.models = TRUE,
    do.full.models = TRUE,
    modeling.id = paste(myRespName,"FirstModeling",sep=""))
   
   cat('\n',myRespName,'modeling finished')  
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_evaluation.txt", sep="")))
  
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(gsub("_",".",myRespName), 
                                paste(myRespName,"_formal_models_variables_importance.txt", sep="")))               
  
 
  
  ### Make projections on current variable
 # myBiomodProjTSS <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),
 #   proj.name = 'current_TSS',
 #   selected.models = 'all',
 #   binary.meth = 'TSS',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	### Make projections on current variable
 # myBiomodProjKAPPA <- BIOMOD_Projection(
 #   modeling.output = myBiomodModelOut,
 #   new.env = subset(myExpl,c(9,12,14,15)),
 #   proj.name = 'current_KAPPA',
 #   selected.models = 'all',
 #   binary.meth = 'KAPPA',
 #   compress = 'xz',
 #   clamping.mask = TRUE,
 #   output.format = '.grd')
	
	## Make projections on current variable
  myBiomodProjROC <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = subset(myExpl,c(9,12,14,15)),
    proj.name = 'current_ROC',
    selected.models = 'all',
    binary.meth = 'ROC',
    compress = 'xz',
    clamping.mask = TRUE,
    output.format = '.grd')
	 