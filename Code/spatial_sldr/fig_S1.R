# rm(list = ls())

# Data from sldr_fit_msa.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


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


#----------population data----------#
# gini of population: urban structure, compact and dispersedCompact and dispersal Gini of population
# method ref: https://doi.org/10.1038/s43588-023-00484-5
# data ref (acs_7_categories_geodata.R): ACS 2015-2019 5 years data, with 1055 cbgs with 0 pop and NA income, the total NA income is 8065
# demo <-  read.csv("/work/rwuirlab/Xinhua/WeiyuLi/Data/Geo/census/acs_2019_5years_cbg.csv", header=TRUE)%>%setDT

demo <-  read.csv("D:/ood/Data/Geo/census/acs_2019_5years_cbg.csv", header=TRUE)%>%setDT
dis.demo <- demo[, .(
  CensusBlockGroup = cbg_fips,
  pop = total_population
)]

pop<-NULL
for(d in 1:(Nmsa+Ndis)){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  
  #----------BlockID, NBlock----------#
  block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  BlockID <- unique(block.msa$CensusBlockGroup)
  NBlock <- length(BlockID)
  
  #----------pop----------#
  demo.msa <- subset(dis.demo,CensusBlockGroup%in%BlockID)%>%setDT
  demo.msa <- left_join(data.frame(CensusBlockGroup=BlockID, msa=Dname), demo.msa, by="CensusBlockGroup")
  pop <- plyr::rbind.fill(pop,demo.msa)
}

#----------Plot the max lag h with sum(W[[h]])>0----------#
pop<- pop%>%as.data.table
for(s in 2:Nregion){
  
  queen.lag <- read.csv(paste(datapath,"/",region[s],"/queen_lag.csv",sep=""), header=TRUE)%>%setDT
  queen.max <- queen.lag[, .SD[which.max(nozero)], by = msa, .SDcols = c("lag", "nozero", "pop")]
  
  #----------demographic data: pop----------#
  pop.msa <- pop[, .(pop.msa = sum(pop)), by = msa]
  queen.max <- merge(queen.max,pop,by='msa')
  fit <- queen.lag[, .(fit0 = coef(lm(pop.adj ~ nozero + 0)),
                       fit = coef(lm(pop.adj ~ nozero))[2]), by = .(msa)]
  queen.max <- merge(queen.max,fit,by='msa')
  
  #----------lag and neighbors' number----------#
  fig.lag.num <- ggplot(data=queen.lag, aes(lag, nozero/10^6))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$lag))+
    scale_y_continuous(limits=range(queen.lag$nozero/10^6))+
    labs(x = "Lag", y = TeX('The number of neighbors (millions)'))+ # expression("The number of neighbors (" %*% 10^6 * ")")
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------lag and neighbors' population----------#
  fig.lag.pop <- ggplot(data=queen.lag, aes(lag, pop.adj/10^6))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$lag))+
    scale_y_continuous(limits=range(queen.lag$pop.adj/10^6))+
    labs(x = "Lag", y = TeX('Population (millions)'))+  # expression("The population of neighbors (" %*% 10^6 * ")")
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------lag and neighbors' area----------#
  fig.lag.area <- ggplot(data=queen.lag, aes(lag, area.adj/10^6))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$lag))+
    scale_y_continuous(limits=range(queen.lag$area.adj/10^6))+
    labs(x = "Lag", y = TeX('Area (km$^2$)'))+  # expression("The area of neighbors (" %*% 10^6 * ")")
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "right",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------lag and neighbors' population/neighbors----------#
  fig.mean.pop <- ggplot(data=queen.lag, aes(lag, pop.adj/nozero))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$lag))+
    scale_y_continuous(limits=range(queen.lag$pop.adj/queen.lag$nozero))+
    labs(x = "Lag", y = "The mean population of neighbors")+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------neighbors and neighbors' population----------#
  fig.num.pop <- ggplot(data=queen.lag, aes(nozero/10^6, pop.adj/10^6))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$nozero/10^6))+
    scale_y_continuous(limits=range(queen.lag$pop.adj/10^6))+
    labs(x = expression("The number of neighbors (" %*% 10^6 * ")"), 
         y = expression("The population of neighbors (" %*% 10^6 * ")"))+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------neighbors and neighbors' area----------#
  fig.num.area <- ggplot(data=queen.lag, aes(nozero/10^6, area.adj/10^6))+
    geom_line(aes(color = msa, group=msa), na.rm = TRUE) +
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.lag$nozero/10^6))+
    scale_y_continuous(limits=range(queen.lag$area.adj/10^6))+
    labs(x = expression("The number of neighbors (" %*% 10^6 * ")"), 
         y = expression("The area of neighbors (" %*% 10^6 * ")"))+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "right",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------lag and max neighbors----------#
  fig.max.num <- ggplot(data=queen.max, aes(lag, nozero/10^6))+
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.max$lag))+
    scale_y_continuous(limits=range(queen.max$nozero/10^6))+
    labs(x = "Lag", y = expression("Maximum neighbors (" %*% 10^6 * ")"))+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          # legend.position = c(0.55, 0.85),  # new
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------pop and max neighbors----------#
  fig.max.pop <- ggplot(data=queen.max, aes(pop/10^6, nozero/10^6))+
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.max$pop/10^6))+
    scale_y_continuous(limits=range(queen.max$nozero/10^6))+
    labs(x = expression("Population size (" %*% 10^6 * ")"), 
         y = expression("Maximum neighbors (" %*% 10^6 * ")"))+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  #----------pop and max neighbors----------#
  fig.pop.linear <- ggplot(data=queen.max, aes(pop/10^6, fit0))+
    geom_point(aes(fill= msa, color=msa, shape=msa),
               size=3, stroke = 0.5) +
    scale_fill_manual(name="MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    guides(fill = guide_legend(override.aes = list(size=3)))+
    scale_x_continuous(limits=range(queen.max$pop/10^6))+
    scale_y_continuous(limits=range(queen.max$fit0))+
    labs(x = expression("Population size (" %*% 10^6 * ")"), 
         y = expression("pop~num"))+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = "right",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  
  fig.queen <- (fig.lag.num|fig.lag.pop|fig.lag.area) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.queen,filename = paste(figpath,"/",region[s],"/SI_queen_lag.pdf",sep=""), width =5*3, height = 4*1)
  
} # region



