library(scales)
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

###going back to biomod methods and analyses

setwd("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2_inputs/datapoints")
spfiles=list.files(pattern="01_less")
spdata=c(
  fc=list(non=read.csv(spfiles[1]),nat=read.csv(spfiles[2])),
  hv=list(non=read.csv(spfiles[3]),nat=read.csv(spfiles[4])),
  mm=list(non=read.csv(spfiles[5]),nat=read.csv(spfiles[6])),
  pm=list(non=read.csv(spfiles[7]),nat=read.csv(spfiles[8]))
)

str(spdata)

lapply(spdata, function(X)summary(as.factor(X[,3])))

##en los invasores se incluye como invasor los puntos fuera de Europa, tanto los de las SOI+Antartida como los américa del norte



spdata_2 <- lapply(names(spdata), function(x) setNames(spdata[[x]], c("x","y","pa")) )
names(spdata_2)=names(spdata)
str(spdata_2)
library(plyr)
library(reshape2)
library(ggplot2)


d=ldply(spdata_2)
head(d)
summary(d)
names(d)[1]="sp.status"
d$sp.status=as.factor(d$sp.status)


ggplot(data=d, aes(x=x,y=y, color=as.factor(pa)))+
  geom_point()+
  facet_grid(~sp.status)
lvl=unique(d$sp.status)
d2=subset(d, sp.status%in%lvl[c(1,3,5,7)])
summary(d2)
d2$pa=as.factor(d2$pa)
d2$status=ifelse(d2$y>0&d2$x>=-14&d2$x<=73,"eur","inv") 
d2$status=as.factor(d2$status)
summary(d2$status)

d2$comb=paste(d2$sp.status, d2$status, sep="_")
d2$comb=as.factor(d2$comb)

ggplot(d2,aes(x=x,y=y, color=pa))+
  geom_point()+
  facet_wrap(status~sp.status, scales="free_y", ncol=4)

ggplot(d2,aes(x=x,y=y, color=pa))+
  geom_point()+
  facet_wrap(~comb, scales="free_y", ncol=4)

#hay que ordenar los puntos para que aparezcan las presencias por encima
d2_order=d2[order(d2$pa,decreasing = FALSE),]
levels(d2_order$comb)[c(1,3,5,7)]=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta","Proisotoma minuta")
library(rgdal)
europe_outline = readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS", layer="europe_dis_all2")
eu_fort=fortify(europe_outline)


europePlot=ggplot(subset(d2_order, status=="eur"), aes(x=x,y=y, color=pa))+
  geom_polygon(data=eu_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(size=0.5)+
  scale_color_manual(values = c("#3182bd","#de2d26"), name=" ", labels=c("pseudo-absence","presence"))+
  #scale_y_continuous(limits =c(30.50000, 83.71375) ) +
  #scale_x_continuous(limits = c(-10,  65.30000))+
  coord_map("conic", lat0 = 30)+
  #coord_fixed(ratio=3)+
  facet_wrap(~comb, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_line(size=0.25), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(.6)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom")+
  labs(x="Longitude",y="Latitude")+
  guides(color = guide_legend(override.aes = list(size=3)))


ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/europe_small.png", dpi=300, height=234, width=84, units="mm")





##ENFA, histniches. 
library(raster)
library(rgdal)
library(adehabitatHS)
library(adehabitatMA)
library(sp)

setwd("C:/Users/greta.vega/Documents/collembola/environmentalBIOVmin_buff_clip") 

bios_names<-c('bio1_Vmin_buff.tif', 'bio2_Vmin_buff.tif', 'bio3_Vmin_buff.tif', 'bio4_Vmin_buff.tif', 'bio5_Vmin_buff.tif', 'bio6_Vmin_buff.tif', 'bio7_Vmin_buff.tif', 'bio8_Vmin_buff.tif', 'bio9_Vmin_buff.tif','bio10_Vmin_buff.tif', 'bio11_Vmin_buff.tif', 'bio12_Vmin_buff.tif', 'bio13_Vmin_buff.tif', 'bio14_Vmin_buff.tif', 'bio15_Vmin_buff.tif', 'bio16_Vmin_buff.tif', 'bio17_Vmin_buff.tif', 'bio18_Vmin_buff.tif', 'bio19_Vmin_buff.tif')

bios<-stack(bios_names)




#To use adehabitatHS functions the raster files have to be SpatialPixelsDataFrame
biosAA<- as(bios, "SpatialPixelsDataFrame")
names(biosAA)<-c("bio1","bio2","bio3","bio4", "bio5","bio6","bio7","bio8","bio9","bio10","bio11","bio12","bio13","bio14","bio15","bio16","bio17","bio18","bio19")
class(biosAA)



#after looking at the correlation of the variables we choose 5
biosAA<-biosAA[,c(2,  4, 5,  15, 17)]
class(biosAA)



#using all the data points from GBIF
fc<-read.csv("C:\\Users\\greta.vega\\Documents\\collembola\\Point_Data\\folsomia_candida_merge_1.csv")
speciesPointsFc<-SpatialPoints(fc[,1:2])
speciesPointsFc@proj4string<- biosAA@proj4string


hv<-read.csv("C:\\Users\\greta.vega\\Documents\\collembola\\Point_Data\\hypogastrura_viatica_merge_1.csv")
speciesPointsHv<-SpatialPoints(hv[,1:2])
speciesPointsHv@proj4string<- biosAA@proj4string

mm<-read.csv("C:\\Users\\greta.vega\\Documents\\collembola\\Point_Data\\mesaphorura_macrochaeta_merge_1.csv")
speciesPointsMm<-SpatialPoints(mm[,1:2])
speciesPointsMm@proj4string<- biosAA@proj4string

pm<-read.csv("C:\\Users\\greta.vega\\Documents\\collembola\\Point_Data\\proisotoma_minuta_merge_1.csv")
speciesPointsPm<-SpatialPoints(pm[,1:2])
speciesPointsPm@proj4string<- biosAA@proj4string

#using native range
fc_nat<-read.csv("C:/Users/greta.vega/Documents/BIOMOD2_inputs/datapoints/folsomia_candida_01_less0_native.csv")
speciesPointsFc_nat<-SpatialPoints(fc_nat[fc_nat$fc_pres==1,1:2])
speciesPointsFc_nat@proj4string<- biosAA@proj4string

hv_nat<-read.csv("C:/Users/greta.vega/Documents/BIOMOD2_inputs/datapoints/hypogastrura_viatica_01_less0_native.csv")
speciesPointsHv_nat<-SpatialPoints(hv_nat[hv_nat$hv_pres==1,1:2])
speciesPointsHv_nat@proj4string<- biosAA@proj4string

mm_nat<-read.csv("C:/Users/greta.vega/Documents/BIOMOD2_inputs/datapoints/mesaphorura_macrochaeta_01_less0_native.csv")
speciesPointsMm_nat<-SpatialPoints(mm_nat[mm_nat$mm_pres==1,1:2])
speciesPointsMm_nat@proj4string<- biosAA@proj4string

pm_nat<-read.csv("C:/Users/greta.vega/Documents/BIOMOD2_inputs/datapoints/proisotoma_minuta_01_less0_native.csv")
speciesPointsPm_nat<-SpatialPoints(pm_nat[pm_nat$pm_pres==1,1:2])
speciesPointsPm_nat@proj4string<- biosAA@proj4string



tab<-slot(biosAA, "data")


#all the points
cpfc <- count.points(speciesPointsFc,biosAA)
cphv <- count.points(speciesPointsHv,biosAA)
cpmm <- count.points(speciesPointsMm,biosAA)
cppm <- count.points(speciesPointsPm,biosAA)

prfc<-slot(cpfc,"data")[,1]
prhv<-slot(cphv,"data")[,1]
prmm<-slot(cpmm,"data")[,1]
prpm<-slot(cppm,"data")[,1]

#native points
cpfc_nat <- count.points(speciesPointsFc_nat,biosAA)
cphv_nat <- count.points(speciesPointsHv_nat,biosAA)
cpmm_nat <- count.points(speciesPointsMm_nat,biosAA)
cppm_nat <- count.points(speciesPointsPm_nat,biosAA)

prfc_nat<-slot(cpfc_nat,"data")[,1]
prhv_nat<-slot(cphv_nat,"data")[,1]
prmm_nat<-slot(cpmm_nat,"data")[,1]
prpm_nat<-slot(cppm_nat,"data")[,1]

#Sub antarctic and antarctic points
ant_fc<-tail(speciesPointsFc)
cpfc_ant <- count.points(ant_fc,biosAA)
prfc_ant<-slot(cpfc_ant,"data")[,1]

ant_hv<-SpatialPoints(hv[hv$POINT_Y<0,1:2])
ant_hv@proj4string<- biosAA@proj4string
cphv_ant <- count.points(ant_hv,biosAA)
prhv_ant<-slot(cphv_ant,"data")[,1]

ant_mm<-SpatialPoints(mm[mm$POINT_Y<0,1:2])
ant_mm@proj4string<- biosAA@proj4string
cpmm_ant <- count.points(ant_mm,biosAA)
prmm_ant<-slot(cpmm_ant,"data")[,1]

ant_pm<-SpatialPoints(pm[pm$POINT_Y<0,1:2])
ant_pm@proj4string<- biosAA@proj4string
cppm_ant <- count.points(ant_pm,biosAA)
prpm_ant<-slot(cppm_ant,"data")[,1]

#putting all the data together
ndf<-cbind(tab, prfc, prfc_nat, prfc_ant, prhv, prhv_nat, prhv_ant, prmm, prmm_nat,prmm_ant, prpm, prpm_nat, prpm_ant)



all_prfc<-ndf[ndf$prfc>0,]
alien_prfc<-all_prfc[all_prfc$prfc_nat==0,]
ant_prfc<-ndf[ndf$prfc_ant>0,]

all_prhv<-ndf[ndf$prhv>0,]
alien_prhv<-all_prhv[all_prhv$prhv_nat==0,]
ant_prhv<-ndf[ndf$prhv_ant>0,]

all_prmm<-ndf[ndf$prmm>0,]
alien_prmm<-all_prmm[all_prmm$prmm_nat==0,]
ant_prmm<-ndf[ndf$prmm_ant>0,]

all_prpm<-ndf[ndf$prpm>0,] #all presences
alien_prpm<-all_prpm[all_prpm$prpm_nat==0,] #all alien
ant_prpm<-ndf[ndf$prpm_ant>0,] #all antarctic


all_prfc_long=melt(subset(all_prfc, select=c(1:5, 6:8)), id.vars = c("bio2", "bio4", "bio5", "bio15", "bio17"))

all_prhv_long=melt(subset(all_prhv, select=c(1:5, 9:11)), id.vars = c("bio2", "bio4", "bio5", "bio15", "bio17"))

all_prmm_long=melt(subset(all_prmm, select=c(1:5, 12:14)), id.vars = c("bio2", "bio4", "bio5", "bio15", "bio17"))

all_prpm_long=melt(subset(all_prpm, select=c(1: 5, 15:17)), id.vars = c("bio2", "bio4", "bio5", "bio15", "bio17"))

head(all_prpm_long)
bin_number=30


##marginality factor plots
fc_bio2=ggplot(data=subset(all_prfc_long, value>0&variable!="prfc"))+
  geom_density(data=ndf, aes(x=bio2, y=..density..*(dim(subset(all_prfc_long, value>0&variable!="prfc"))[1]*bin_number)), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio2, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

hv_bio2=ggplot(data=subset(all_prhv_long, value>0&variable!="prhv"))+
  geom_density(data=ndf, aes(x=bio2, y=..density..*(dim(subset(all_prhv_long, value>0&variable!="prfc"))[1]*bin_number)), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio2, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

mm_bio2=ggplot(data=subset(all_prmm_long, value>0&variable!="prmm"))+
  geom_density(data=ndf, aes(x=bio2, y=..density..*(dim(subset(all_prhv_long, value>0&variable!="prmm"))[1]*bin_number)), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio2, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pm_bio2=ggplot(data=subset(all_prpm_long, value>0&variable!="prpm"))+
  geom_density(data=ndf, aes(x=bio2, y=..density..*(dim(subset(all_prpm_long, value>0&variable!="prpm"))[1]*bin_number)), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio2, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
library(gridExtra)
library(grid)  
grid.arrange(fc_bio2, hv_bio2, mm_bio2, pm_bio2, ncol=1)
ggsave(grid.arrange(fc_bio2, hv_bio2, mm_bio2, pm_bio2, ncol=1), filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/marginalityPlots.png", dpi=300, height=234, width=58.5, units="mm")


##specialisation factor plots
r=70000
bin_number_spe=30
fc_bio17=ggplot(data=subset(all_prfc_long, value>0&variable!="prfc"), aes(x=bio17))+
  geom_density(data=ndf, aes(x=bio17, y=..density..*(dim(subset(all_prfc_long, value>0&variable!="prfc"))[1]*bin_number_spe)/r), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio17, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number_spe)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

hv_bio15=ggplot(data=subset(all_prhv_long, value>0&variable!="prhv"))+
  geom_density(data=ndf, aes(x=bio15, y=..density..*(dim(subset(all_prhv_long, value>0&variable!="prfc"))[1]*bin_number_spe)/800), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio15, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number_spe)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

mm_bio17=ggplot(data=subset(all_prmm_long, value>0&variable!="prmm"))+
  geom_density(data=ndf, aes(x=bio17, y=..density..*(dim(subset(all_prhv_long, value>0&variable!="prmm"))[1]*bin_number_spe)/r), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio17, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number_spe)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


pm_bio17=ggplot(data=subset(all_prpm_long, value>0&variable!="prpm"))+
  geom_density(data=ndf, aes(x=bio17, y=..density..*(dim(subset(all_prpm_long, value>0&variable!="prpm"))[1]*bin_number_spe)/r), fill = "antiquewhite3", alpha=.2, colour="darkgray")+
  geom_histogram(aes(x=bio17, fill=variable, y=..count..), alpha=.7, 
                 #colour="black", 
                 bins=bin_number_spe)+
  scale_fill_manual(values = c("#3182bd","#de2d26"),name="Range", labels=c("Native","Antarctic"), guide=FALSE)+
  scale_y_continuous(trans="S_sqrt",name="Count")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

grid.arrange(fc_bio17, hv_bio15, mm_bio17, pm_bio17, ncol=1)


ggsave(grid.arrange(fc_bio17, hv_bio15, mm_bio17, pm_bio17, ncol=1), filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/specialisationPlots.png", dpi=300, height=234, width=58.5, units="mm")



grid.arrange(fc_bio2, hv_bio2, mm_bio2, pm_bio2,fc_bio17, hv_bio15, mm_bio17, pm_bio17, ncol=2, as.table=FALSE)

ggsave(grid.arrange(fc_bio2, hv_bio2, mm_bio2, pm_bio2,fc_bio17, hv_bio15, mm_bio17, pm_bio17, ncol=2, as.table=FALSE), filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/NichePlots.png", dpi=300, height=234, width=174, units="mm")
lay=rbind(c(1, 6,6,6,11,11,11),
          c(2, 7,7,7,12,12,12),
          c(3, 8,8,8,13,13,13),
          c(4, 9,9,9,14,14,14),
          c(5, 10,10,10,15,15,15))
t_empty=textGrob(" ")
t_fc <- textGrob("Folsomia Candida", rot=90, just="centre", gp=gpar(fontface="italic"))
t_hv <- textGrob("Hypogastrura viatica", rot=90, just="centre", gp=gpar(fontface="italic"))
t_mm <- textGrob("Mesaphorura macrochaeta", rot=90, just="centre", gp=gpar(fontface="italic"))
t_pm <- textGrob("Proisotoma minuta", rot=90, just="centre", gp=gpar(fontface="italic"))
t_mar=textGrob("Marginality factor", just="bottom")
t_spe=textGrob("Specialisation factor", just="bottom")

ggsave(grid.arrange(t_empty, t_fc, t_hv, t_mm, t_pm, t_mar, fc_bio2, hv_bio2, mm_bio2, pm_bio2,t_spe, fc_bio17, hv_bio15, mm_bio17, pm_bio17, layout_matrix = lay, as.table=FALSE), filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/NichePlots_labs.png", dpi=300, height=234, width=174, units="mm")


##extracting CA and CV from biomod to create the histograms
###using native and invaded range

fc_ca=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida_Extracted/Folsomia.candida_Ensembling/prob/Folsomia.candida_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
fc_cv=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida_Extracted/Folsomia.candida_Ensembling/prob/Folsomia.candida_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


hv_ca=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_Ensembling/prob/Hypogastrura.viatica_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
hv_cv=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_Ensembling/prob/Hypogastrura.viatica_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")

mm_ca=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_Ensembling/prob/Mesaphorura.macrochaeta_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
mm_cv=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_Ensembling/prob/Mesaphorura.macrochaeta_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


pm_ca=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta_Extracted/Proisotoma.minuta_Ensembling/prob/Proisotoma.minuta_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
pm_cv=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta_Extracted/Proisotoma.minuta_Ensembling/prob/Proisotoma.minuta_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


extent(europe_outline)
sp_ca_cv=stack(fc_ca, fc_cv, hv_ca, hv_cv, mm_ca, mm_cv, pm_ca, pm_cv)

sp_ca_cv_crop=crop(x = sp_ca_cv, y = europe_outline)
sp_ca_cv_cMask=mask(x = sp_ca_cv_crop, mask = europe_outline)

clamping_fc=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida_Extracted/Folsomia.candida_modelling/Folsomia.candida_ROC_ClampingMask.tif")

clamping_hv=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_modelling/Hypogastrura.viatica_ROC_ClampingMask.tif")

clamping_mm=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_modelling/Mesaphorura.macrochaeta_ROC_ClampingMask.tif")

clamping_pm=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta_Extracted/Proisotoma.minuta_modelling/Proisotoma.minuta_ROC_ClampingMask.tif")
clamping=stack(clamping_fc, clamping_hv, clamping_mm,clamping_pm)

clamping_crop=crop(x = clamping, y = europe_outline)
clamping_cMask=mask(x = clamping_crop, mask = europe_outline)

sp_ca_cv_clamping=stack(sp_ca_cv_cMask, clamping_cMask)

names(sp_ca_cv_clamping)=c("fc_ca", "fc_cv", "hv_ca", "hv_cv", "mm_ca", "mm_cv", "pm_ca", "pm_cv", "clamping_fc", "clamping_hv", "clamping_mm","clamping_pm")
ens_data=as.data.frame(rasterToPoints(sp_ca_cv_clamping))
dim(ens_data)
names(ens_data)
head(ens_data)

ens_data_fc=ens_data[, c(1:2,grep("fc", names(ens_data)))]
ens_data_hv=ens_data[, c(1:2,grep("hv", names(ens_data)))]
ens_data_mm=ens_data[, c(1:2,grep("mm", names(ens_data)))]
ens_data_pm=ens_data[, c(1:2,grep("pm", names(ens_data)))]

ens_data_fc$ca_fct=cut(ens_data_fc$fc_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)
table(ens_data_fc$clamping_fc)
ens_data_hv$ca_fct=cut(ens_data_hv$hv_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_mm$ca_fct=cut(ens_data_mm$mm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_pm$ca_fct=cut(ens_data_pm$pm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_list=list(ens_data_fc, ens_data_hv, ens_data_mm,ens_data_pm )
names(ens_data_list)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")
ens_data_list=lapply(ens_data_list, setNames, c("x","y","ca","cv", "clamping","ca_fct"))
str(ens_data_list)
ens_data_all=ldply(ens_data_list)
summary(ens_data_all)
names(ens_data_all)[1]="species"
ens_data_all$species=as.factor(ens_data_all$species)
library(viridis)
ggplot(data=na.omit(ens_data_all))+
  geom_histogram(aes(x=cv, fill=ca_fct, alpha=as.factor(clamping)), bins = 20, color="black")+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_europe20_bis.png", dpi=300, height=234*3/4, width=174, units="mm")



ggplot(data=na.omit(ens_data_all)[1:100000,])+
  geom_histogram(aes(x=cv, fill=ca_fct, alpha=as.factor(clamping)), bins = 20, color="black")+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Variables extrapolated")+
  scale_fill_viridis(discrete = TRUE,  labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage", direction = -1, option = "plasma")+
  scale_y_continuous(trans="S_sqrt")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_europe20_magma.png", dpi=300, height=234*3/4, width=174, units="mm")


summary(ens_data_all)
ens_data_all$clamping=as.factor(ens_data_all$clamping)
cv.mean=acast(ens_data_all, formula=clamping~ca_fct~species, fun.aggregate = mean, value.var="cv", margins=TRUE)
cv.mean

cv.count=acast(ens_data_all, formula=clamping~ca_fct~species, fun.aggregate = length, value.var="cv", margins=TRUE)
cv.count
str(cv.count)
for (i in 1:4){
  print(dimnames(cv.count)[[3]][i])
  print(round(cv.count[,,i]*100/cv.count[,7,i],2))}
for (i in 1:4){
  print(dimnames(cv.count)[[3]][i])
  print(round(t(t(cv.count[,1:5,i])*100/cv.count[5,1:5,i]),2))}

t(t(cv.count[,1:5,1])*100/cv.count[5,1:5,1])
#####Antarctic maps of proyection. 

ant_extent=extent(-80,-40,-75,-60)
sp_ca_cv
clamping
world_all=stack(sp_ca_cv, clamping)
ant_all=crop(world_all, ant_extent)
names(ant_all)=c("fc_ca", "fc_cv", "hv_ca", "hv_cv", "mm_ca", "mm_cv", "pm_ca", "pm_cv", "clamping_fc", "clamping_hv", "clamping_mm","clamping_pm")

ant_all_df=as.data.frame(rasterToPoints(ant_all))

ant_fc=subset(ant_all_df, select=c(1:2,grep("fc", names(ant_all_df))))

ant_hv=subset(ant_all_df, select=c(1:2,grep("hv", names(ant_all_df))))

ant_mm=subset(ant_all_df, select=c(1:2,grep("mm", names(ant_all_df))))

ant_pm=subset(ant_all_df, select=c(1:2,grep("pm", names(ant_all_df))))

ant_fc$ca_fct=cut(ant_fc$fc_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_hv$ca_fct=cut(ant_hv$hv_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_mm$ca_fct=cut(ant_mm$mm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_pm$ca_fct=cut(ant_pm$pm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)




ant_data_list=list(ant_fc, ant_hv, ant_mm,ant_pm )
names(ant_data_list)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")
ant_data_list=lapply(ant_data_list, setNames, c("x","y","ca","cv", "clamping", "ca_fct"))
str(ant_data_list)
ant_data_all=ldply(ant_data_list)
summary(ant_data_all)
names(ant_data_all)[1]="species"
ant_data_all$species=as.factor(ant_data_all$species)
ant_data_all=ant_data_all[order(ant_data_all$ca),]

library(RColorBrewer)
getPaletteCA=colorRampPalette(brewer.pal(9, "RdYlBu"))

ggplot(data=ant_data_all)+
  geom_point(aes(x=x,y=y, colour=ca))+
  scale_colour_gradientn(colours=rev(getPaletteCA(12)))+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_line(size=0.25), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="right")+
  labs(x="Longitude",y="Latitude")
  

ggplot(data=na.omit(ant_data_all))+
  geom_histogram(aes(x=cv, fill=ca_fct, alpha=as.factor(clamping)), bins = 20)+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20.png", dpi=300, height=234*3/4, width=174, units="mm")

summary(ant_data_all)
ant_data_all$clamping=as.factor(ant_data_all$clamping)
cv.mean_antAll=acast(ant_data_all, formula=clamping~ca_fct~species, fun.aggregate = mean, value.var="cv", margins=TRUE)
cv.mean_antAll

cv.count_antAll=acast(ant_data_all, formula=clamping~ca_fct~species, fun.aggregate = length, value.var="cv", margins=TRUE)
cv.count_antAll
str(cv.count_antAll)
for (i in 1:4){
  print(dimnames(cv.count_antAll)[[3]][i])
  print(round(cv.count_antAll[,,i]*100/cv.count_antAll[,6,i],2))}

for (i in 1:4){
  print(dimnames(cv.count_antAll)[[3]][i])
  print(round(t(t(cv.count_antAll[,1:5,i])*100/cv.count_antAll[5,1:5,i]),2))}



###loading results with only the native range

fc_ca_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida_Extracted/Folsomia.candida_Ensembling/prob/Folsomia.candida_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
fc_cv_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida_Extracted/Folsomia.candida_Ensembling/prob/Folsomia.candida_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


hv_ca_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_Ensembling/prob/Hypogastrura.viatica_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
hv_cv_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_Ensembling/prob/Hypogastrura.viatica_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")

mm_ca_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_Ensembling/prob/Mesaphorura.macrochaeta_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
mm_cv_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_Ensembling/prob/Mesaphorura.macrochaeta_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


pm_ca_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta_Extracted/Proisotoma.minuta_Ensembling/prob/Proisotoma.minuta_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")
pm_cv_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta_Extracted/Proisotoma.minuta_Ensembling/prob/Proisotoma.minuta_EMcvByROC_mergedAlgo_mergedRun_mergedData_ROC.tif")


sp_ca_cv_nat=stack(fc_ca_nat, fc_cv_nat, hv_ca_nat, hv_cv_nat, mm_ca_nat, mm_cv_nat, pm_ca_nat, pm_cv_nat)

sp_ca_cv_crop_nat=crop(x = sp_ca_cv_nat, y = europe_outline)
sp_ca_cv_nat_cMask=mask(x = sp_ca_cv_crop_nat, mask = europe_outline)

clamping_fc_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida_Extracted/Folsomia.candida_modelling/Folsomia.candida_ROC_ClampingMask.tif")

clamping_hv_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica_Extracted/Hypogastrura.viatica_modelling/Hypogastrura.viatica_ROC_ClampingMask.tif")

clamping_mm_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta_Extracted/Mesaphorura.macrochaeta_modelling/Mesaphorura.macrochaeta_ROC_ClampingMask.tif")

clamping_pm_nat=raster("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta_Extracted/Proisotoma.minuta_modelling/Proisotoma.minuta_ROC_ClampingMask.tif")
clamping_nat=stack(clamping_fc_nat, clamping_hv_nat, clamping_mm_nat,clamping_pm_nat)

clamping_crop_nat=crop(x = clamping_nat, y = europe_outline)
clamping_nat_cMask=mask(x = clamping_crop_nat, mask = europe_outline)

sp_ca_cv_clamping_nat=stack(sp_ca_cv_nat_cMask, clamping_nat_cMask)

names(sp_ca_cv_clamping_nat)=c("fc_ca", "fc_cv", "hv_ca", "hv_cv", "mm_ca", "mm_cv", "pm_ca", "pm_cv", "clamping_fc", "clamping_hv", "clamping_mm","clamping_pm")
ens_nat_data=as.data.frame(rasterToPoints(sp_ca_cv_clamping_nat))


ens_data_fc_nat=ens_nat_data[, c(1:2,grep("fc", names(ens_nat_data)))]
ens_data_hv_nat=ens_nat_data[, c(1:2,grep("hv", names(ens_nat_data)))]
ens_data_mm_nat=ens_nat_data[, c(1:2,grep("mm", names(ens_nat_data)))]
ens_data_pm_nat=ens_nat_data[, c(1:2,grep("pm", names(ens_nat_data)))]

ens_data_fc_nat$ca_fct=cut(ens_data_fc_nat$fc_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_hv_nat$ca_fct=cut(ens_data_hv_nat$hv_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_mm_nat$ca_fct=cut(ens_data_mm_nat$mm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_pm_nat$ca_fct=cut(ens_data_pm_nat$pm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ens_data_nat_list=list(ens_data_fc_nat, ens_data_hv_nat, ens_data_mm_nat, ens_data_pm_nat )
names(ens_data_nat_list)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")
ens_data_nat_list=lapply(ens_data_nat_list, setNames, c("x","y","ca","cv", "clamping","ca_fct"))
str(ens_data_nat_list)
ens_data_all_nat=ldply(ens_data_nat_list)
summary(ens_data_all_nat)
names(ens_data_all_nat)[1]="species"
ens_data_all_nat$species=as.factor(ens_data_all_nat$species)

ggplot(data=na.omit(ens_data_all_nat))+
  geom_histogram(aes(x=cv, fill=ca_fct, alpha=as.factor(clamping)), bins = 20)+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_europe20_nat.png", dpi=300, height=234*3/4, width=174, units="mm")


ant_extent=extent(-80,-40,-75,-60)

world_all_nat=stack(sp_ca_cv_nat, clamping_nat)
ant_all_nat=crop(world_all_nat, ant_extent)
names(ant_all_nat)=c("fc_ca", "fc_cv", "hv_ca", "hv_cv", "mm_ca", "mm_cv", "pm_ca", "pm_cv", "clamping_fc", "clamping_hv", "clamping_mm","clamping_pm")

ant_all_df_nat=as.data.frame(rasterToPoints(ant_all_nat))

ant_fc_nat=subset(ant_all_df_nat, select=c(1:2,grep("fc", names(ant_all_df_nat))))

ant_hv_nat=subset(ant_all_df_nat, select=c(1:2,grep("hv", names(ant_all_df_nat))))

ant_mm_nat=subset(ant_all_df_nat, select=c(1:2,grep("mm", names(ant_all_df_nat))))

ant_pm_nat=subset(ant_all_df_nat, select=c(1:2,grep("pm", names(ant_all_df_nat))))

ant_fc_nat$ca_fct=cut(ant_fc_nat$fc_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_hv_nat$ca_fct=cut(ant_hv_nat$hv_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_mm_nat$ca_fct=cut(ant_mm_nat$mm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

ant_pm_nat$ca_fct=cut(ant_pm_nat$pm_ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)




ant_nat_data_list=list(ant_fc_nat, ant_hv_nat, ant_mm_nat,ant_pm_nat )
names(ant_nat_data_list)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")
ant_nat_data_list=lapply(ant_nat_data_list, setNames, c("x","y","ca","cv", "clamping", "ca_fct"))
str(ant_nat_data_list)
ant_data_all_nat=ldply(ant_nat_data_list)
summary(ant_data_all_nat)
names(ant_data_all_nat)[1]="species"
ant_data_all_nat$species=as.factor(ant_data_all_nat$species)
ant_data_all_nat=ant_data_all_nat[order(ant_data_all_nat$ca),]

library(RColorBrewer)
getPaletteCA=colorRampPalette(rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")))

ggplot(data=ant_data_all_nat)+
  geom_point(aes(x=x,y=y, colour=ca))+
  scale_colour_gradientn(colours=rev(getPaletteCA(12)))+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_line(size=0.25), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="right")+
  labs(x="Longitude",y="Latitude")


ggplot(data=na.omit(ant_data_all_nat))+
  geom_histogram(aes(x=cv, fill=ca_fct, 
                     alpha=as.factor(clamping)
                     ), 
                 bins = 20)+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_nat.png", dpi=300, height=234*3/4, width=174, units="mm")


##use disaggregate to clip by ice free ACBRs
acbrs_antpen=readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS", layer="ACBRs_clipAntPen")

summary(acbrs_antpen)

ant_all
ant_all_nat

acbrs_antpen_1984=spTransform(acbrs_antpen, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


ant_all_dis=disaggregate(ant_all, fact=200)
ant_all_iceFree=mask(ant_all_dis, acbrs_antpen_1984)
ant_all_nat_dis=disaggregate(ant_all_nat, fact=200)
ant_all_nat_iceFree=mask(ant_all_nat_dis, acbrs_antpen_1984)


acbrs_raster=rasterize(acbrs_antpen_1984, ant_all_nat_iceFree, field=acbrs_antpen_1984$ACBR_ID)
ant_all_iceFree_acbr=stack(ant_all_iceFree, acbrs_raster)
ant_all_iceFree_nat_acbr=stack(ant_all_nat_iceFree, acbrs_raster)

df_ant_all_iceFree=as.data.frame(rasterToPoints(ant_all_iceFree_acbr))
df_ant_all_iceFree_nat=as.data.frame(rasterToPoints(ant_all_iceFree_nat_acbr))

###check the acbr numbers, make as factor and use map values
names(df_ant_all_iceFree_nat)[15]="acbr"
names(df_ant_all_iceFree)[15]="acbr"
df_ant_all_iceFree$acbr=as.factor(df_ant_all_iceFree$acbr)
df_ant_all_iceFree_nat$acbr=as.factor(df_ant_all_iceFree_nat$acbr)

ggplot(data=df_ant_all_iceFree_nat)+
  geom_point(aes(x=x, y=y, colour=acbr))
library(plyr)
df_ant_all_iceFree_nat$acbr=mapvalues(df_ant_all_iceFree_nat$acbr, 
          from =c("1", "2","3","4","5"),
          to=c("1","15", "2", "3", "4"))
df_ant_all_iceFree$acbr=mapvalues(df_ant_all_iceFree$acbr, 
          from =c("1", "2","3","4","5"),
          to=c("1","15", "2", "3", "4"))

#list and ldply
list_df_ant_all_iceFree=list(
  df_ant_all_iceFree[, c(1:2, 15, grep("fc", names(df_ant_all_iceFree)))],
  df_ant_all_iceFree[, c(1:2, 15, grep("hv", names(df_ant_all_iceFree)))],
  df_ant_all_iceFree[, c(1:2, 15, grep("mm", names(df_ant_all_iceFree)))],
  df_ant_all_iceFree[, c(1:2, 15, grep("pm", names(df_ant_all_iceFree)))])
names(list_df_ant_all_iceFree)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")

list_df_ant_all_nat_iceFree=list(
  df_ant_all_iceFree_nat[, c(1:2, 15, grep("fc", names(df_ant_all_iceFree_nat)))],
  df_ant_all_iceFree_nat[, c(1:2, 15, grep("hv", names(df_ant_all_iceFree_nat)))],
  df_ant_all_iceFree_nat[, c(1:2, 15, grep("mm", names(df_ant_all_iceFree_nat)))],
  df_ant_all_iceFree_nat[, c(1:2, 15, grep("pm", names(df_ant_all_iceFree_nat)))])
names(list_df_ant_all_nat_iceFree)=c("Folsomia candida", "Hypogastrura viatica", "Mesaphorura macrochaeta", "Proisotoma minuta")

list_df_ant_all_iceFree=lapply(list_df_ant_all_iceFree, setNames, c("x","y","acbr", "ca","cv", "clamping"))
list_df_ant_all_nat_iceFree=lapply(list_df_ant_all_nat_iceFree, setNames, c("x","y","acbr", "ca","cv", "clamping"))

short_df_ant_all_iceFree=ldply(list_df_ant_all_iceFree)
short_df_ant_all_nat_iceFree=ldply(list_df_ant_all_nat_iceFree)

names(short_df_ant_all_iceFree)[1]="species"
names(short_df_ant_all_nat_iceFree)[1]="species"

short_df_ant_all_iceFree$species=as.factor(short_df_ant_all_iceFree$species)
short_df_ant_all_nat_iceFree$species=as.factor(short_df_ant_all_nat_iceFree$species)

short_df_ant_all_iceFree$ca_fct=cut(short_df_ant_all_iceFree$ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)
short_df_ant_all_nat_iceFree$ca_fct=cut(short_df_ant_all_nat_iceFree$ca, breaks=c(0,200,400,600,800,1000), include.lowest=TRUE)

#redo histograms, only with icefree data.
nat_acbr_hist=ggplot(data=short_df_ant_all_nat_iceFree)+
  geom_histogram(aes(x=cv, fill=ca_fct, 
                     alpha=as.factor(clamping)
  ), 
  bins = 20,
  color="black",
  size=.5)+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_grid(acbr~species)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")


all_acbr_hist=ggplot(data=short_df_ant_all_iceFree)+
  geom_histogram(aes(x=cv, fill=ca_fct, 
                     alpha=as.factor(clamping)
  ), 
  bins = 20,
  color="black",
  size=.5)+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Clamping")+
  scale_fill_manual(values = rev(c("#d7191c", "#fdae61","#ffffbf","#abd9e9","#2c7bb6")), labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage")+
  theme_bw()+
  facet_grid(acbr~species)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(plot = nat_acbr_hist, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_nat_iceFree3.png", dpi=300, height=234*3/4, width=174, units="mm")
ggsave(plot = all_acbr_hist, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree3.png", dpi=300, height=234*3/4, width=174, units="mm")

summary(short_df_ant_all_nat_iceFree)
short_df_ant_all_nat_iceFree$clamping=as.factor(short_df_ant_all_nat_iceFree$clamping)
cv.mean_antNat=acast(short_df_ant_all_nat_iceFree, formula=clamping~ca_fct~species, fun.aggregate = mean, value.var="cv", margins=TRUE)
cv.mean_antNat

cv.count_antNat=acast(short_df_ant_all_nat_iceFree, formula=clamping~ca_fct~species, fun.aggregate = length, value.var="cv", margins=TRUE)
cv.count_antNat


for (i in 1:4){
  print(dimnames(cv.count_antNat)[[3]][i])
  print(round(cv.count_antNat[,,i]*100/cv.count_antNat[,5,i],2))}

for (i in 1:4){
  print(dimnames(cv.count_antNat)[[3]][i])
  print(round(t(t(cv.count_antNat[1:3,,i])*100/cv.count_antNat[4,1:5,i]),2))}


summary(short_df_ant_all_iceFree)
short_df_ant_all_iceFree$clamping=as.factor(short_df_ant_all_iceFree$clamping)
cv.mean_antAll=acast(short_df_ant_all_iceFree, formula=clamping~ca_fct~species, fun.aggregate = mean, value.var="cv", margins=TRUE)
cv.mean_antAll

cv.count_antAll=acast(short_df_ant_all_iceFree, formula=clamping~ca_fct~species, fun.aggregate = length, value.var="cv", margins=TRUE)
cv.count_antAll


for (i in 1:4){
  print(dimnames(cv.count_antAll)[[3]][i])
  print(round(cv.count_antAll[,,i]*100/cv.count_antAll[,6,i],2))}

for (i in 1:4){
  print(dimnames(cv.count_antAll)[[3]][i])
  print(round(t(t(cv.count_antAll[1:4,,i])*100/cv.count_antAll[5,1:6,i]),2))}



###redo map
pen_outline = readOGR(dsn="E:/BAS 2018/ArcGIS", layer="Pen_seamask2")
#spTransform
pen_fort=fortify(pen_outline)

short_df_ant_all_nat_iceFree=short_df_ant_all_nat_iceFree[order(short_df_ant_all_nat_iceFree$ca),]
short_df_ant_all_iceFree=short_df_ant_all_iceFree[order(short_df_ant_all_iceFree$ca),]

nat_acbr_map=ggplot(data=short_df_ant_all_nat_iceFree)+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  scale_colour_gradientn(colours=getPaletteCA(10)[1:8], breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/30,"mm"))+
  labs(x="Longitude",y="Latitude")


all_acbr_map=ggplot(data=short_df_ant_all_iceFree)+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  scale_colour_gradientn(colours=getPaletteCA(10), breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")


ggsave(plot = nat_acbr_map, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_nat_iceFree_Map3.png", dpi=300, height=234*3/4, width=174, units="mm")
ggsave(plot = all_acbr_map, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree_Map3.png", dpi=300, height=234*3/4, width=174, units="mm")




####figures with viridis and square root axis
cvcaEuropeplot=ggplot(data=na.omit(ens_data_all))+
  geom_histogram(aes(x=cv, fill=ca_fct, alpha=as.factor(clamping)), bins = 20, color="black")+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Variables extrapolated")+
  scale_fill_viridis(discrete = TRUE,  labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage", direction = -1, option = "plasma")+
  scale_y_continuous(trans="S_sqrt")+
  theme_bw()+
  facet_wrap(~species, ncol=2)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")

ggsave(plot=cvcaEuropeplot, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_europe20_magma.png", dpi=300, height=234*3/4, width=174, units="mm")


nat_acbr_hist=ggplot(data=short_df_ant_all_nat_iceFree)+
  geom_histogram(aes(x=cv, fill=ca_fct, 
                     alpha=as.factor(clamping)
  ), 
  bins = 20,
  color="black",
  size=.5)+
  scale_fill_viridis(discrete = TRUE, begin=0.2, labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage", direction = -1, option = "plasma")+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25)[2:4], name="Variables extrapolated")+
  scale_y_continuous(trans="S_sqrt")+
  theme_bw()+
  facet_grid(acbr~species)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")
ggsave(plot = nat_acbr_hist, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_nat_iceFree3_plasma.png", dpi=300, height=234*3/4, width=174, units="mm")

all_acbr_hist=ggplot(data=short_df_ant_all_iceFree)+
  geom_histogram(aes(x=cv, fill=ca_fct, 
                     alpha=as.factor(clamping)
  ), 
  bins = 20,
  color="black",
  size=.5)+
  scale_fill_viridis(discrete = TRUE, labels=c("0-200","200-400","400-600","600-800", "800-1000"),name="Committee \naverage", direction = -1, option = "plasma")+
  scale_alpha_manual(values=c(1,0.75,0.5,0.25), name="Variables extrapolated")+
  scale_y_continuous(trans="S_sqrt")+
  theme_bw()+
  facet_grid(acbr~species)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "vertical",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom")+
  labs(x="Coefficient of Variation",y="Count")


ggsave(plot = all_acbr_hist, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree3_plasma.png", dpi=300, height=234*3/4, width=174, units="mm")



all_acbr_map_vir=ggplot(data=subset(short_df_ant_all_iceFree, clamping==0))+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  #scale_colour_viridis(colours=getPaletteCA(10), breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  scale_colour_viridis(discrete = FALSE,name="Committee \naverage", direction = -1, option = "plasma")+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")


ggsave(plot = all_acbr_map_vir, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree_Map_viridis_0clamped.png", dpi=300, height=234*3/4, width=174, units="mm")

short_df_ant_all_iceFreeNAs=short_df_ant_all_iceFree
short_df_ant_all_iceFreeNAs$ca[short_df_ant_all_iceFreeNAs$clamping!=0]=NA
short_df_ant_all_iceFreeNAs=short_df_ant_all_iceFreeNAs[order(short_df_ant_all_iceFreeNAs$ca, na.last=FALSE, decreasing=FALSE),]

###obtain longitude, latitude of the non extraolated areas for the three species except Hv. 
ss=subset(short_df_ant_all_iceFreeNAs, species=="Folsomia candida"&clamping==0)
dim(ss)
ggplot(data=ss)+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  coord_cartesian(xlim=c(
    min(ss$x), 
    max(ss$x)), ylim=c(min(ss$y), max(ss$y)),clip = "on")


all_acbr_map_vir_clamping0=ggplot(data=short_df_ant_all_iceFreeNAs)+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  #scale_colour_viridis(colours=getPaletteCA(10), breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  scale_colour_viridis(discrete = FALSE,name="Committee \naverage", direction = -1, option = "plasma", 
                       end=1-(min(subset(short_df_ant_all_iceFree, clamping==0, select="ca"))/1000))+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")


ggsave(plot = all_acbr_map_vir_clamping0, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree_Map_viridis_0clamped.png", dpi=300, height=234*3/4, width=174, units="mm")


pen_outline_in = readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS", layer="penOutline")
pen_out_fort=fortify(pen_outline_in)
xx=na.omit(short_df_ant_all_iceFreeNAs)

all_acbr_map_vir_clamping0_stereo=ggplot(data=short_df_ant_all_iceFreeNAs)+
  geom_polygon(data=pen_out_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  #scale_colour_viridis(colours=getPaletteCA(10), breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  scale_colour_viridis(discrete = FALSE,name="Committee \naverage", direction = -1, option = "plasma", 
                       end=1-(min(subset(short_df_ant_all_iceFree, clamping==0, select="ca"))/1000),
                       na.value="darkgray")+
  scale_y_continuous(breaks = c(-60, -65, -70, -75)) +
  scale_x_continuous(breaks = c(-80,-70,-60,-50,-40))+
  coord_map("ortho", orientation = c(-90, 0, 65),ylim = c(-60,-75),clip = "on", xlim=c(-80,-40))+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_line(color="gray"), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")

#all_acbr_map_vir_clamping0_stereo

ggsave(plot = all_acbr_map_vir_clamping0_stereo, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree_Map_viridis_0clamped_ortho2.png", dpi=300, height=234, width=174, units="mm")





 all_acbr_map_vir_noclamping=ggplot(data=short_df_ant_all_iceFree)+
  geom_polygon(data=pen_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_point(aes(x=x,y=y, colour=ca), size=1)+
  #scale_colour_viridis(colours=getPaletteCA(10), breaks=c(0,200,400,600,800,1000), name="Committee \naverage")+
  scale_colour_viridis(discrete = FALSE,name="Committee \naverage", direction = -1, option = "plasma")+
  facet_wrap(~species, ncol=2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")


ggsave(plot = all_acbr_map_vir_noclamping, filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/cv_ca_AntPen20_all_iceFree_Map_viridis_noclamping.png", dpi=300, height=234*3/4, width=174, units="mm")




#finding out what happend to the alien points during the runs
str(spdata)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/.BIOMOD_DATA/Folsomia_candidaFirstModeling/calib.lines")
str(calib.lines)
head(calib.lines)
calib.lines[1:10,,]
fc.calib.lines=calib.lines

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica/.BIOMOD_DATA/Hypogastrura_viaticaFirstModeling/calib.lines")
hv.calib.lines=calib.lines

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta/.BIOMOD_DATA/Mesaphorura_macrochaetaFirstModeling/calib.lines")
mm.calib.lines=calib.lines
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta/.BIOMOD_DATA/Proisotoma_minutaFirstModeling/calib.lines")
pm.calib.lines=calib.lines

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/.BIOMOD_DATA/Folsomia_candidaFirstModeling/formated.input.data")

fc.data.coord=data@coord

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica/.BIOMOD_DATA/Hypogastrura_viaticaFirstModeling/formated.input.data")

hv.data.coord=data@coord

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta/.BIOMOD_DATA/Mesaphorura_macrochaetaFirstModeling/formated.input.data")

mm.data.coord=data@coord

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta/.BIOMOD_DATA/Proisotoma_minutaFirstModeling/formated.input.data")

pm.data.coord=data@coord


list_calib.lines=list(
  Folsomia_candida=data.frame(fc.data.coord, fc.calib.lines),
  Hypogastrura_viatica=data.frame(hv.data.coord, hv.calib.lines),
  Mesaphorura_macrochaeta=data.frame(mm.data.coord, mm.calib.lines),
  Proisotoma_minuta=data.frame(pm.data.coord, pm.calib.lines)
)
str(list_calib.lines)

df_calib.lines=ldply(list_calib.lines)
head(df_calib.lines)
names(df_calib.lines)[1]="species"
df_calib.lines$species=as.factor(df_calib.lines$species)

df_calib.linesAnt=subset(df_calib.lines, POINT_Y<0)
dim(df_calib.linesAnt)
head(df_calib.linesAnt)
df_calib.linesAntLong=melt(df_calib.linesAnt, id.vars = c("species", "POINT_X", "POINT_Y"))
df_calib.linesAntLong=df_calib.linesAntLong[df_calib.linesAntLong$variable!="X_Full._AllData",]
head(df_calib.linesAntLong)
summary(df_calib.linesAntLong)

ggplot(data=subset(df_calib.linesAntLong, species%in%c("Mesaphorura_macrochaeta", "Proisotoma_minuta")))+
  geom_jitter(aes(x=POINT_X, y=POINT_Y, fill=value, shape=variable, size=value),width = .0, height = .0, alpha=.5)+
  scale_size_manual(values = c(4,1))+
  scale_shape_manual(values = c(21, 22, 24))+
  facet_wrap(~species, nrow=2)+
  theme_bw()
  
write.csv(df_calib.linesAntLong, "C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS/calibLines.csv")


###check was has happened with the GLMs
grep("GLM",list.files(path="C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/models/Folsomia_candidaFirstModeling"),value = TRUE)

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN1_GLM")
get_formal_model(Folsomia.candida_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN2_GLM")
get_formal_model(Folsomia.candida_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN3_GLM")
get_formal_model(Folsomia.candida_AllData_RUN3_GLM)

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN1_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN2_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN3_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN3_GLM)

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN1_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN2_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN3_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN3_GLM)


load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN1_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN2_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_28052016/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN3_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN3_GLM)







load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN1_GLM")
get_formal_model(Folsomia.candida_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN2_GLM")
get_formal_model(Folsomia.candida_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Folsomia.candida/models/Folsomia_candidaFirstModeling/Folsomia.candida_AllData_RUN3_GLM")
get_formal_model(Folsomia.candida_AllData_RUN3_GLM)

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN1_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN2_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Hypogastrura.viatica/models/Hypogastrura_viaticaFirstModeling/Hypogastrura.viatica_AllData_RUN3_GLM")
get_formal_model(Hypogastrura.viatica_AllData_RUN3_GLM)

load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN1_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN2_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Mesaphorura.macrochaeta/models/Mesaphorura_macrochaetaFirstModeling/Mesaphorura.macrochaeta_AllData_RUN3_GLM")
get_formal_model(Mesaphorura.macrochaeta_AllData_RUN3_GLM)


load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN1_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN1_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN2_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN2_GLM)
load("C:/Users/greta.vega/Dropbox/Capitulo siguiente/BIOMOD2/BIOMOD2/output_29052016_native/Proisotoma.minuta/models/Proisotoma_minutaFirstModeling/Proisotoma.minuta_AllData_RUN3_GLM")
get_formal_model(Proisotoma.minuta_AllData_RUN3_GLM)




##making WSS islands map
pen_outline_in = readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS", layer="penOutline")
pen_out_fort=fortify(pen_outline_in)

acbrs_antpen_fort=fortify(acbrs_antpen_1984, region="ACBR_ID")

summary(acbrs_antpen_fort)
acbrs_antpen_fort$id_fct=factor(acbrs_antpen_fort$id, levels=c("1","2","3","4","15"))
levels(acbrs_antpen_fort$id_fct)
acbrs_antpen_fort_order=acbrs_antpen_fort[order(acbrs_antpen_fort$id_fct, decreasing=TRUE),]

rectShape=data.frame(xmin=-63, xmax=-60.4, ymax=-63.6, ymin=-62.4)


acbr_map_stereo=ggplot()+
  geom_polygon(data=pen_out_fort, aes(x=long, y=lat, group=group), fill="lightgray", colour="transparent")+
  geom_polygon(data=acbrs_antpen_fort_order, aes(x=long, y=lat, group=group, fill=id_fct, colour=id_fct), size=2)+
  scale_colour_manual(values = c(rgb(115,0,0, maxColorValue = 255), rgb(255,127,127, maxColorValue = 255), rgb(255,0,0, maxColorValue = 255), rgb(255,255,0, maxColorValue = 255), rgb(245,122,182, maxColorValue = 255)), name="ACBR")+
  scale_fill_manual(values = c(rgb(115,0,0, maxColorValue = 255), rgb(255,127,127, maxColorValue = 255), rgb(255,0,0, maxColorValue = 255), rgb(255,255,0, maxColorValue = 255), rgb(245,122,182, maxColorValue = 255)), name="ACBR")+
  geom_polygon(data=subset(acbrs_antpen_fort_order, id_fct==1), aes(x=long, y=lat, group=group), fill=rgb(115,0,0, maxColorValue = 255), colour=rgb(115,0,0, maxColorValue = 255), size=2)+
  geom_rect(data=rectShape, aes(xmin=xmin, xmax=xmax, ymax=ymax, ymin=ymin), fill="transparent", colour="black", size=1)+
  scale_y_continuous(breaks = c(-60, -65, -70, -75)) +
  scale_x_continuous(breaks = c(-80,-70,-60,-50,-40))+
  coord_map("ortho", orientation = c(-90, 0, 65),ylim = c(-60,-75),clip = "on", xlim=c(-80,-40))+
  
  theme_bw()+
  theme(panel.grid.major = element_line(color="gray"), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="right",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")
acbr_map_stereo


###read data for W SSI
WSScoast=readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS/WSSIs", layer="coast_WSSIslands2")
WSScoast_fort=fortify(WSScoast)


IBApoints=readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS/WSSIs", layer="IBA_WSSIslands")
IBApoints_fort=fortify(IBApoints)

ice_freeWSS=readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS/WSSIs", layer="iceFree_WSSIslands2")
ice_freeWSS_1984=spTransform(ice_freeWSS, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ice_freeWSS_fort=fortify(ice_freeWSS_1984)

ASPAS=readOGR(dsn="C:/Users/greta.vega/Dropbox/Capitulo siguiente/arcGIS/WSSIs", layer="ASPAs_WSSIslands")
ASPAS_1984=spTransform(ASPAS, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

WSS_map_names=data.frame(x=c(-62.386,-61.018,-60.76, -60.644,-60.515, -61.233, -61.494, -61.022,-62.299, -61.867),
                         y=c(-63.377,-62.574,-62.42, -62.852,-63.046, -62.475, -62.863, -63.041, -63.059, -63.221),
                         lab=c("152","126","149","145","140", "Livingston I.", "Snow I.", "Deception I.", "Smith I.", "Low I."), type=c("ASPA","ASPA","ASPA","ASPA","ASPA", "Island", "Island", "Island", "Island", "Island"))

library(ggsn)

WSS_map=ggplot()+
  geom_polygon(data=WSScoast_fort, aes(x=long, y=lat, group=group), fill="lightgray")+
  geom_polygon(data=ice_freeWSS_fort, aes(x=long, y=lat, group=group), fill=rgb(255,0,0, maxColorValue = 255), colour=rgb(255,0,0,maxColorValue = 255))+
  geom_polygon(data=ASPAS_fort, aes(x=long, y=lat, group=group), fill="transparent", colour="blue")+
  geom_point(data=IBApoints@data, aes(x=Longitude, y=Latitude), shape=21, fill="white", alpha=.5, size=3)+
  scale_y_continuous(breaks = c(-63)) +
  scale_x_continuous(breaks = c(-63,-62,-61))+
  geom_text(data=WSS_map_names, aes(x=x, y=y, label=lab, colour=type))+
  scale_colour_manual(values=c("blue","black"), guide=FALSE)+
  coord_map("ortho", orientation = c(-90, 0, 61.5),ylim = c(-62.4,-63.6),clip = "on", xlim=c(-63,-60.4))+
  theme_bw()+
  theme(panel.grid.major = element_line(color="gray"), 
        panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        strip.text=element_text(face="italic", size = rel(1)),
        axis.ticks.length=unit(-0.125, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.25,0.25,0.25,0.25), "cm")),
        legend.position="bottom",
        legend.key.width=unit(174/15,"mm"))+
  labs(x="Longitude",y="Latitude")
WSS_map
ggsave(grid.arrange( WSS_map,acbr_map_stereo), filename = "C:/Users/greta.vega/Dropbox/Capitulo siguiente/_writing/figures/WSS_ggplot.png", dpi=300, height=234, width=174, units="mm")
