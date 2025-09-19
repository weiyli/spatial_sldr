# rm(list = ls())

# Data from sldr_fit.R

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Load packages----------#
library(ggsignif)   # gemo_signif()
library(ggtext)     # element_markdown()


#----------Load msa----------#
Dname.msa<-c('Atlanta',
             'Boston',
             'Chicago',
             'Dallas', 
             'Houston',
             'Los Angeles',
             'Miami',
             'New York',
             'Philadelphia',
             'San Francisco',
             'Seattle',
             'Washington, D.C.')
Dname.dis<-c('Houston',
             'Jacksonville',
             'Houston',
             'Santa Rosa-Petaluma')
Dname.set <- c(Dname.msa,Dname.dis)
d<<-1
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
Flow.name<-"Intra_Flow"


Dname<<-Dname.set[1]
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif<-date.sep(Datevalue,durdate)

#----------Data----------#
for(yindex in 1:Ynum){
  
  SMAPE<-NULL
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    date.gif<-date.sep(Datevalue,durdate)
    
    #----------Data from sldr_fit.R: SLDR_fit----------# 
    sldr.fit <- read.csv(paste(datapath,"/msa/fit/",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
    sldr.fit <- melt(sldr.fit, id.vars = c("day", "empirical"), 
                     measure.vars = c("homo", "exp", "power"),
                     variable.name = "model", value.name = "predicted")
    error.msa <- sldr.fit[, .(error = error_test(empirical, predicted)), by = .(model,day)]
    error.msa$msa <- Dname
    SMAPE <- plyr::rbind.fill(SMAPE, error.msa)
  } # msa
  write.csv(SMAPE, file=paste(datapath,"/msa/SMAPE_",Yname[yindex],".csv",sep=""), row.names = FALSE)
} # yindex


#----------Plot----------# 
A.msa <- B.msa <- list()
model.info<-data.frame(breaks=c("homo","power","exp"),
                       labels=c("No decay","Power law","Exponential"),
                       col=col.all[1:3],
                       shape=c(21,22,24))
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")
id.week.during <- paste(format(Datevalue[Day/2+(1:(Day/2))],"%m-%d"),format(Datevalue[Day/2+(1:(Day/2))], "%a"),sep=" ")
for(yindex in 1:Ynum){
  
  SMAPE <- read.csv(paste(datapath,"/msa/SMAPE_",Yname[yindex],".csv",sep=""), header=TRUE)
  SMAPE$day <- as.Date(SMAPE$day)
  SMAPE <- left_join(SMAPE,date.gif,by="day")
  SMAPE.before <- subset(SMAPE,period==pervalue[1])
  SMAPE.before <- left_join(SMAPE.before,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')
  
  SMAPE.during <- subset(SMAPE,period==pervalue[2])
  SMAPE.during <- left_join(SMAPE.during,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')
  
  #----------Before----------# 
  A.msa[[yindex]] <- ggplot(SMAPE.before, aes(x = day.id, y = error)) +
    geom_line(aes(color = model), linewidth = 1) +
    geom_point(aes(color = model, fill = model, shape=model), size=2, stroke = 1) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    scale_x_continuous(breaks = c(1:(Day/2)),labels = id.week)+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = "Date", y = "SMAPE", title = NULL) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain",size=15),
          axis.text.y = element_text(hjust = 0),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          legend.justification = c(0,0.5)) 
  ggsave(A.msa[[yindex]], filename = paste(figpath,"/msa/SI_fit_model_before_",Yname[yindex],".pdf",sep=""), width = 4*4, height = 4*3)
  
  #----------During----------# 
  B.msa[[yindex]] <- ggplot(SMAPE.during, aes(x = day.id, y = error)) +
    geom_line(aes(color = model), linewidth = 1) +
    geom_point(aes(color = model, fill = model, shape=model), size=2, stroke = 1) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    scale_x_continuous(breaks = c(Day/2+(1:(Day/2))),labels = id.week.during)+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = "Date", y = "SMAPE", title = NULL) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain",size=15),
          axis.text.y = element_text(hjust = 0),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          legend.justification = c(0,0.5)) 
  ggsave(B.msa[[yindex]], filename = paste(figpath,"/msa/SI_fit_model_during_",Yname[yindex],".pdf",sep=""), width = 4*4, height = 4*3)
  
} # Yname



