library(ggplot2)
library(grid)
library(RColorBrewer)
library(viridis)
setwd("C:/Users/greta.vega/Dropbox/Capitulo siguiente/evaluations")
setwd("/Users/gretom/Desktop/Greta/polar biology/r_stuff/")
eval_data=read.csv("allEvals07122018.csv", sep=";")
summary(eval_data)

eval_data$TSS_fct=cut(eval_data$TSS, breaks = c(0,0.05,0.2,0.4,0.55,0.7,0.8,0.99,1),
                      labels = c("No agreement", "Very Poor", "Poor", "Fair", "Good", "Very Good", "Excellent", "Perfect"),
                      include.lowest = TRUE)



eval_data$TSS_fct2=cut(eval_data$TSS, breaks = c(0,0.05,0.2,0.4,0.55,0.7,0.85,1),
                       labels = c("No agreement", "Very Poor", "Poor", "Fair", "Good", "Very Good", "Excellent"),
                       include.lowest = TRUE)

n2=factor(c("No agreement", "Very Poor", "Poor", "Fair", "Good", "Very Good", "Excellent", "Perfect"), 
          levels = c("No agreement","Very Poor", "Poor", "Fair", "Good", "Very Good", "Excellent", "Perfect"))
n3=c(0.015,0.1,0.3,0.475,0.625,0.775,0.9,1.01)
coln= palG(length(n3))
let= c(0, rep(100,6), 0)
N2=data.frame(n2,n3,coln, let)




p2=ggplot( )+
  geom_jitter(data=subset(eval_data, run!="Full"),aes(y=TSS, x=algo, color=Sensitivity, shape=status), alpha=1, width = 0.25, height = 0, size = 2)+
  #scale_fill_gradientn(colours = c(palG(6)))+
  #scale_fill_viridis(direction = -1, option = "viridis")+
  scale_color_gradient(low="white", high="black")+
  scale_shape_manual(values = c(16,17),name="Range", label=c("Alien and Native", "Only Native"))+
  geom_label(data=N2, aes(y=n3,label=n2, fill = n2, color = let),
             x = 9, # Set the position of the text 
            hjust = 0,
            #vjust = 1,
            size = rel(3)) +
  facet_wrap(~species)+
  theme_classic()+
  scale_y_continuous(name = "TSS value",
                     breaks = c(0.05,0.2,0.4,0.55,0.7,0.85,1),
                     limits = c(0,1.05))+
  scale_x_discrete(name="Algorithms", labels=c("ANN", "CTA", "FDA", "GAM", "BRT", "GLM", "MARS", "RF"))+
  scale_fill_manual(name= "category", values = N2$coln, guide = "none")+
  theme(
        panel.grid.major.y = element_line(color="darkgray",size = .25,linetype = 1), 
        #panel.grid.minor = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.box = "horizontal",
        strip.text=element_text(face="italic", size = rel(1)),
        legend.position="bottom", 
        plot.margin = unit(c(1,6,1,1), "lines"))+ # This widens the right margin
    coord_cartesian(xlim = c(1,8), # This focuses the x-axis on the range of interest
                 # ylim=c(0,1),
                  clip = 'off')    # This keeps the labels from disappearing
p2

ggsave(plot = p2, filename = "evalPlot_bw_small.png",width = 100,units = "mm",dpi = 300)






library(reshape2)

eval_TSS=subset(eval_data,run!="Full", select=c("species", "algo", "TSS", "status"))
TSS.mean=acast(data=eval_TSS,formula=algo~species~status,fun.aggregate=mean, value.var = "TSS",margins = TRUE, na.rm=TRUE)
TSS.mean[,,1]-TSS.mean[,,2]


eval_Sens=subset(eval_data,run!="Full", select=c("species", "algo", "Sensitivity", "status"))

Sensitivity.mean=acast(data=eval_Sens,formula=algo~species~status,fun.aggregate=mean, value.var = "Sensitivity",margins = TRUE, na.rm=TRUE)
Sensitivity.mean[,,1]-Sensitivity.mean[,,2]



acc_class=subset(eval_data,run!="Full", select=c("species", "algo", "TSS_fct", "status"))
acc_class_count=acast(data=acc_class,formula=TSS_fct~species~status,fun.aggregate=length, value.var = "algo", margins = TRUE)

acc_class_count[,,1]-acc_class_count[,,2]


##checking what is happening for Pm and Mm as Mm has only one run with the alien presence. 
evalTSS_pm_mm=subset(eval_data,run!="Full"&species%in%c("Mesaphorura macrochaeta", "Proisotoma minuta"), select=c("species", "algo", "TSS", "status", "run") )
summary(evalTSS_pm_mm)
acast(data=evalTSS_pm_mm,formula=run~species~status,fun.aggregate=mean, value.var = "TSS",margins = TRUE, na.rm=TRUE)



