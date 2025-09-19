# rm(list = ls())

# Data from sldr_fit.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Load packages----------#
library(ggsignif)   # gemo_signif()
library(ggtext)     # element_markdown()
library(ggrepel)
library(sf)


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


#----------Part1: COVID-19----------#
datapath <- '/home/weiy.li/Data/spatial_sldr'
figpath <- '/home/weiy.li/Figure/spatial_sldr'


#----------Plot the fitting results: spearman rank correlation and map----------#
A.msa <- B.msa <- C.msa <- D.msa <- list()
id.week<-paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")
Dname<<-Dname.set[1]
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif<-date.sep(Datevalue,durdate)

# for(yindex in 1:Ynum){
for(yindex in 1:2){
  
  #----------data for plotting Pearson's correlation: data form sldr_fit.R----------#
  dis.cor <- read.csv(paste(datapath,"/msa/SLDR_rank_cor_",Yname[yindex],".csv",sep=""), header=TRUE)
  dis.cor$day <- as.Date(dis.cor$day)
  dis.cor <- left_join(dis.cor,date.gif,by='day')
  dis.cor.before <- subset(dis.cor,period==pervalue[1])
  dis.cor.before <- left_join(dis.cor.before,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')
  
  #----------data for plotting the map with max freq R2----------#
  dis.r2 <- read.csv(paste(datapath,"/msa/SLDR_r2_",Yname[yindex],".csv",sep=""), header=TRUE)
  dis.r2$day <- as.Date(dis.r2$day)
  dis.r2 <- left_join(dis.r2,date.gif,by='day')
  max.r2.day <- setDT(dis.r2)[, .SD[which.max(r2_exp)], by = .(msa, period), .SDcols = c("day", "r2_exp")]
  max.r2 <- max.r2.day[, .N, by = .(period, day)][order(-N)][, .SD[1], by = period]
  max.date <- c(as.Date(max.r2[period==pervalue[1],day]),as.Date(max.r2[period==pervalue[2],day]))
  print(max.date[1])
  max.day <- length(max.date)
  
  #----------plot Pearson's correlation: data form sldr_fit.R----------#
  A.msa[[yindex]] <- ggplot(dis.cor.before, aes(day.id, pearson_exp))+  
    geom_line(aes(color = msa, group = msa), linetype = "dotted", na.rm = TRUE) +
    geom_point(aes(fill= msa, color = msa, shape = msa), size=3, stroke = 0.5) +
    # geom_label_repel(data = subset(dis.cor.before, day == max.date[1] & msa%in%c("Atlanta","Boston")),
    #                  aes(label = msa,color=msa), max.overlaps = Inf,
    #                  box.padding = 1, size = 3, fill = 'white') +
    # geom_label_repel(data = subset(dis.cor.before, day == max.date[1] & msa %in% c("Atlanta", "Boston")),
    #                  aes(label = paste0(msa, "\nr = ", round(pearson_exp, 2),
    #                                     ifelse(msa=="Atlanta", ", corr(a,b)", ", corr(d,e)")), color = msa),
    #   max.overlaps = Inf, box.padding = 1, size = 3, fill = 'white') +
    geom_label_repel(data = subset(dis.cor.before, day == max.date[1] & msa == "Atlanta"),
                     aes(label = paste0(msa, "\nr = ", round(pearson_exp, 2),", corr(a,b)")),
                     max.overlaps = Inf, box.padding = 1, size = 4, fill = 'white',
                     color=msa.info$col[which(msa.info$dname=="Atlanta")]) +
    geom_label_repel(data = subset(dis.cor.before, day == max.date[1] & msa == "Boston"),
                     aes(label = paste0(msa, "\nr = ", round(pearson_exp, 2),", corr(d,e)")),
                     max.overlaps = Inf, box.padding = 1, size = 4, fill = 'white',
                     color=msa.info$col[which(msa.info$dname=="Boston")]) +
    scale_fill_manual(name = 'MSAs',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = 'MSAs',
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col) +
    scale_shape_manual(name = 'MSAs',
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    scale_x_continuous(breaks = c(1:(Day/2)), labels = id.week)+
    labs(x = "Date", y = "Pearson's correlation coefficient")+
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          # legend.title = element_blank(), 
          legend.key.size = unit(0.5, "cm"),
          legend.position = 'right',
          legend.text = element_text(hjust=0),
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  
  
  #----------Plot the map with max freq R2----------#
  for(d in 1:Nmsa){
    # for(d in c(1,2)){
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    date.gif<-date.sep(Datevalue,durdate)
    
    #----------MSA BlockID, NBlock----------#
    block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
    BlockID<-unique(block.msa$CensusBlockGroup)
    NBlock<-length(BlockID)
    #----------Center point of each block----------#
    block.msa$centroid <- st_centroid(block.msa$geometry,of_largest_polygon = TRUE)
    s1 <- st_coordinates(block.msa$centroid)
    block.msa$lon<-s1[,1]
    block.msa$lat<-s1[,2]
    
    #----------Generate edges of network and W: 1 lag queen adjacency matrices----------#
    W <- spdep::poly2nb(block.msa, queen = T) 
    edges<-NULL
    for(i in 1:length(W)){
      edges<-plyr::rbind.fill(edges,data.frame(from=i,to=W[[i]]))
    }
    dis.edges <- distinct(edges)
    colnames(dis.edges) <- c("from", "to")
    s1<-data.frame(from=1:NBlock,to=1:NBlock,lon=block.msa$lon,lat=block.msa$lat)
    s2<-left_join(dis.edges,s1,by="from")
    s3<-left_join(dis.edges,s1,by="to")
    dis.edges$Olon<-s2$lon
    dis.edges$Olat<-s2$lat
    dis.edges$Dlon<-s3$lon
    dis.edges$Dlat<-s3$lat
    
    #----------Data from sldr_fit.R: SLDR_fit----------# 
    sldr.fit <- read.csv(paste(datapath,"/msa/fit/",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
    sldr.fit$day<-as.Date(sldr.fit$day)
    sldr.fit<-left_join(sldr.fit,date.gif,by="day")
    sldr.fit$period<-left_join(sldr.fit,data.frame(period=pervalue,ID=ID),by="period")$ID
    sldr.fit$res<-sldr.fit$empirical-sldr.fit$exp
    sldr.fit <- left_join(sldr.fit, data.frame(day=Datevalue,day.id=factor(c(1:Day))),by='day')
    
    #----------Fig.2: boxplot for residual----------#
    lower_bound<--500
    upper_bound<-500
    sldr.fit.box <- subset(sldr.fit,period==1) %>%filter(res >= lower_bound & res <= upper_bound)
    B.msa[[d]] <- ggplot(sldr.fit.box, aes(x=day.id, y=res)) +
      geom_violin(alpha = 0.5, linewidth = 0.5, color = "gray2", fill=msa.info$col[d]) +
      geom_boxplot(alpha = 0.5, width=0.05, color="gray2", fill=msa.info$col[d],
                   notch=F, outlier.shape=NA)+
      geom_signif(comparisons = list(c("1", "2"), c("2", "3"), c("3", "4"),
                                     c("4", "5"), c("5", "6"), c("6", "7")), 
                  y_position = 450,  tip_length = 0, vjust = 2,
                  map_signif_level = TRUE, test = "t.test") +
      labs(title = Dname, x="Date", y = "<span style = 'color: #000000;'>Residual R (</span><span style = 'color: #eaac8b;'>Empirical</span>-<span style = 'color: #355070;'>Model<span style = 'color: #000000;'>)</span>") +
      # TeX('Residual $\\textit{R}$') +
      scale_x_discrete(breaks = factor(c(1:(Day/2))), labels = id.week)+
      # scale_y_continuous(limits=c(lower_bound,upper_bound),expand=c(0.2,0.2))+
      theme_wy() +
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
            plot.title = element_text(color = msa.info$col[d]),
            axis.title.y = ggtext::element_markdown(),
            # axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")
    
    #----------Fig2: map of flows----------# 
    land.fit<-subset(sldr.fit,day==max.date[1])
    land.fit$lon<-block.msa$lon
    land.fit$lat<-block.msa$lat
    #----------Empirical flows----------#
    df.nodes<-data.frame(lon=land.fit$lon,lat=land.fit$lat,mob=land.fit$empirical)
    df.nodes$mob<-df.nodes$mob/max(df.nodes$mob)
    C.msa[[d]] <- ggplot() +
      geom_segment(data = na.omit(dis.edges), aes(x = Olon, y = Olat, xend = Dlon, yend = Dlat),
                   arrow = arrow(length = unit(0, "inches")), color = "gray")+
      geom_point(data = df.nodes, aes(x = lon, y = lat, color=mob, size= mob), alpha=0.8) +
      # scale_colour_viridis_c(name='Flows',option = 'inferno') +
      scale_colour_gradientn(name = Dname, colors = col.fit) +
      scale_size_continuous(name='Flows', range = c(0,5))+
      labs(title=DM[1])+
      theme_wy()+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
            plot.title = element_text(color = tail(col.fit,1)), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none")+
      guides(color = guide_colorbar(order = 1),
             size = guide_legend(order = 2),
             override.aes = list(size=3))
    #----------Model flows----------#
    df.nodes<-data.frame(lon=land.fit$lon,lat=land.fit$lat,mob=land.fit$exp)
    df.nodes$mob<-df.nodes$mob/max(df.nodes$mob)
    D.msa[[d]] <- ggplot() +
      geom_segment(data = na.omit(dis.edges), aes(x = Olon, y = Olat, xend = Dlon, yend = Dlat),
                   arrow = arrow(length = unit(0, "inches")), color = "gray")+
      geom_point(data = df.nodes, aes(x = lon, y = lat, color=mob, size= mob), alpha=0.8) +
      scale_colour_gradientn(name = Dname, colors = col.fit) +
      scale_size_continuous(name='Flows', range = c(0,5))+
      labs(title=DM[2])+
      theme_wy()+
      theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
            plot.title = element_text(color = col.fit[1]), 
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_text(color = msa.info$col[d]),
            legend.position = "right",
            legend.justification = c(0,0.5))+
      guides(color = guide_colorbar(order = 1),
             size = guide_legend(order = 2),
             override.aes = list(size=3))
  } # msa
  
  #----------Model performance for Atlanta and Boston----------#
  fig.fit <- ((((C.msa[[1]]|D.msa[[1]])/B.msa[[1]])|((C.msa[[2]]|D.msa[[2]])/B.msa[[2]]))/A.msa[[yindex]]) + plot_layout(heights = c(2, 1)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.fit, filename = paste(figpath,"/msa/fit_",Yname[yindex],".pdf",sep=""), width = 5*4, height = 4.8*3)
  
  #----------SI fitting: Model performance for additional MSAs----------#
  fig.fit.si <- (C.msa[[3]]|D.msa[[3]]|C.msa[[4]]|D.msa[[4]])/
    (C.msa[[5]]|D.msa[[5]]|C.msa[[6]]|D.msa[[6]])/
    (C.msa[[7]]|D.msa[[7]]|C.msa[[8]]|D.msa[[8]])/
    (C.msa[[9]]|D.msa[[9]]|C.msa[[10]]|D.msa[[10]])/
    (C.msa[[11]]|D.msa[[11]]|C.msa[[12]]|D.msa[[12]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.fit.si, filename = paste(figpath,"/msa/SI_fit_",Yname[yindex],".pdf",sep=""), width =5*4, height = 4.5*5)
  #----------SI residuals: Model performance for additional MSAs----------#
  fig.res.si <- (B.msa[[3]]|B.msa[[4]])/(B.msa[[5]]|B.msa[[6]])/
    (B.msa[[7]]|B.msa[[8]])/(B.msa[[9]]|B.msa[[10]])/
    (B.msa[[11]]|B.msa[[12]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.res.si, filename = paste(figpath,"/msa/SI_fit_res_",Yname[yindex],".pdf",sep=""), width = 5*4, height = 5*5)
  
} # Yname



