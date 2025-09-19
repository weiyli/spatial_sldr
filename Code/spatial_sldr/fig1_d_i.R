# rm(list = ls())

# county level data from sldr_fit_county.R; msa level data from sldr_fit_msa.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(sf)          # read_sf() 
library(spdep)       # poly2nb() 
library(gridExtra)
library(ggsignif)   # gemo_signif()
library(ggtext)     # element_markdown()
library(cowplot)

#----------Error-measure: R2----------#
R2_test <- function(y_test,y_predict){
  SS.tot = sum((y_test-mean(y_test))^2)
  SS.err = sum((y_predict-y_test)^2)
  result <- 1-SS.err/SS.tot
  return(result)
}

get.significance <- function(p) {
  if (p < 0.001) {
    return("p<0.001")
  } else if (p < 0.01) {
    return("p<0.01")
  } else if (p < 0.05) {
    return("p<0.05")
  } else {
    return("")
  }
}

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


#----------rho and msa or county----------#
A.msa <- B.msa <- C.msa <- D.msa <- list()

#----------pop and area for county and msa----------#
queen.lag.county <- read.csv(paste(datapath,"/county/queen_lag.csv",sep=""), header=TRUE)%>%setDT
queen.lag.county$pop <- queen.lag.county$pop/10^6  # million
queen.lag.county$area <- queen.lag.county$area/10^6 # km2
queen.max.county <- queen.lag.county[, .(msa, county, nblock, lag.max, pop, area)]%>%distinct

queen.lag.msa <- read.csv(paste(datapath,"/msa/queen_lag.csv",sep=""), header=TRUE)%>%setDT
queen.lag.msa$pop <- queen.lag.msa$pop/10^6  # million
queen.lag.msa$area <- queen.lag.msa$area/10^6 # km2
queen.max.msa <- queen.lag.msa[, .(msa, nblock, lag.max, pop, area)]%>%distinct

for(yindex in 2:2){
  
  #----------data from sldr_fit.R----------#
  dis.rho <- read.csv(paste(datapath,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  
  #----------rho and msa sorted by median_rho----------#
  rho.msa <- dis.rho[, .(
    min_rho = min(rho),
    median_rho = median(rho),
    max_rho = max(rho)
  ), by = .(msa, index, period)]
  rho.msa <- rho.msa[index == "exp"&period==pervalue[1]]
  rho.rank <- rho.msa[order(median_rho)]
  rho.rank[, msa_rank := .I]
  
  all.rho <- c(rho.rank$min_rho, rho.rank$median_rho, rho.rank$max_rho)
  mean(all.rho)
  sd(all.rho)
  
  
  A.msa[[yindex]] <- ggplot(rho.rank, aes(x = msa_rank, y = median_rho, color = msa, shape = msa)) +
    geom_errorbar(aes(ymin = min_rho, ymax = max_rho), width = 0.3, linewidth = 0.5) +  # Vertical error bars
    geom_point(aes(fill = msa, color = msa, shape = msa),
               size = 3, stroke = 0.5) +
    scale_fill_manual(name = "MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col, alpha = 0.8)) +
    scale_colour_manual(name = "MSA",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col) +
    scale_shape_manual(name = "MSA",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    scale_x_continuous(breaks = rho.rank$msa_rank, labels = rho.rank$msa) +  
    scale_y_continuous(limits = c(0.35, 0.5)) +  
    # labs(x="MSA", y = TeX('Spatial range exponent $\\rho$')) +  
    labs(x=NULL, y = TeX('Spatial range exponent $\\rho$')) +  
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.justification = c(0,0.5),
          legend.position = "right",
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  
  #----------power-law scaling model: log(rho)=rho0+alpha*log(pop)----------#
  #----------msa level: data from sldr_fit.R----------#
  rho.msa <- read.csv(paste(datapath,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  sldr.msa <- rho.msa[, .SD[rho == median(rho)], by = .(msa,period,index)]
  sldr.msa <- sldr.msa[index == "exp"&period==pervalue[1]]
  sldr.msa <- sldr.msa[, .(msa, day, rho)]
  sldr.msa[, county := region[2]]
  
  rho.msa <- rho.msa[, .(
    min_rho = quantile(rho, probs = 0),
    Q1_rho = quantile(rho, probs = 0.25),
    median_rho = quantile(rho, probs = 0.5),
    Q3_rho = quantile(rho, probs = 0.75),
    max_rho = quantile(rho, probs = 1),
    median_day = day[which(rho==quantile(rho, probs = 0.5))]
  ), by = .(msa, index, period)]
  rho.msa <- rho.msa[index == "exp"&period==pervalue[1]]
  rho.msa <- merge(rho.msa, queen.max.msa, by=c('msa'))
  var.msa <- rho.msa[, .(msa, median_day, median_rho, pop, area)]
  var.msa[, county := "msa"]
  
  #----------conuty level: data from urban_sldr_fit.R----------#
  rho.county <- read.csv(paste(datapath,"/county/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  sldr.county <- rho.county[, .SD[rho == median(rho)], by = .(msa,county,period,index)]
  sldr.county <- sldr.county[index == "exp"&period==pervalue[1]]
  sldr.county <- sldr.county[, .(msa, day, rho, county)]
  dis.sldr <- plyr::rbind.fill(sldr.msa,sldr.county)
  
  rho.county <- rho.county[, .(
    min_rho = quantile(rho, probs = 0),
    Q1_rho = quantile(rho, probs = 0.25),
    median_rho = quantile(rho, probs = 0.5),
    Q3_rho = quantile(rho, probs = 0.75),
    max_rho = quantile(rho, probs = 1),
    median_day = day[which(rho==quantile(rho, probs = 0.5))]
  ), by = .(county, msa, index, period)]
  rho.county <- rho.county[index == "exp"&period==pervalue[1]]
  rho.county <- merge(rho.county, queen.max.county, by=c('msa','county'))
  var.county <- rho.county[, .(msa, median_day, median_rho, pop, area,county)]
  
  dis.var <- plyr::rbind.fill(var.msa, var.county)
  
  #----------linear fitting rho and pop for county and msa----------#
  # rho.pop.model <- lm(log(dis.var$median_rho) ~ log(dis.var$pop))
  dis.var <- subset(dis.var,pop>0.1)
  all.rho <- dis.var$median_rho
  mean(all.rho)
  sd(all.rho)
  
  rho.pop.model <- lm(log(dis.var$median_rho) ~ dis.var$pop)
  rho.pop.coef <- coef(rho.pop.model)
  rho.pop.coef[1] <- exp(rho.pop.coef[1])
  rho0 <- round(rho.pop.coef[1],2)
  alpha <- round(rho.pop.coef[2],2)
  R2 <- summary(rho.pop.model)$r.squared%>%round(2)
  
  rho.pop.pval <- coef(summary(rho.pop.model))[, "Pr(>|t|)"]
  intercept.signif <- get.significance(rho.pop.pval[1])
  pop.signif <- get.significance(rho.pop.pval[2])
  
  # #----------plot linear fitting of rho and pop----------#
  # B.msa[[yindex]] <- ggplot(dis.var, aes(x = pop, y = median_rho)) +
  #   geom_point(data=subset(dis.var,county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 2, stroke = 0.5) +
  #   geom_point(data=subset(dis.var,county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 4, stroke = 1) +
  #   geom_smooth(method = 'lm', color = "#69b3a2") +
  #   scale_fill_manual(name="MSA: counties",
  #                     breaks = msa.info$dname,
  #                     labels = msa.info$dname,
  #                     values = scales::alpha(msa.info$col,alpha=0.8))+
  #   scale_colour_manual(name = "MSA: counties",
  #                       breaks = msa.info$dname,
  #                       labels = msa.info$dname,
  #                       values = msa.info$col)+
  #   scale_shape_manual(name="MSA: counties",
  #                      breaks = msa.info$dname,
  #                      labels = msa.info$dname,
  #                      values = msa.info$shape) +
  #   scale_x_log10(limits=range(dis.var$pop)) +
  #   scale_y_log10(limits=range(dis.var$median_rho)) +
  #   # scale_x_continuous(
  #   #   trans = log10_trans(),
  #   #   breaks = trans_breaks("log10", function(x) 10^x),
  #   #   labels = trans_format("log10", math_format(10^.x)))+
  #   # scale_y_continuous(
  #   #   trans = log10_trans(),
  #   #   breaks = trans_breaks("log10", function(x) 10^x),
  #   #   labels = trans_format("log10", math_format(10^.x)))+
  #   labs(x = TeX('Population (millions), $\\log P$'),
  #        y = TeX('Spatial range exponent, $\\log \\rho$')) +
  #   theme_wy() +
  #   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
  #         legend.title = element_blank(),
  #         legend.position = "none") +
  #   annotation_custom(grob = textGrob(TeX(sprintf("$\\rho = %.2f P^{%.2f}$", rho0, alpha)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = 0.9, xmax = 0.9, ymin = -0.29, ymax = -0.29) +
  #   annotation_custom(grob = textGrob(TeX(sprintf("$R^2 = %.2f$", R2)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = 0.9, xmax = 0.9, ymin = -0.3, ymax = -0.3)
  
  #----------plot linear fitting of rho and pop----------#
  B.msa[[yindex]] <- ggplot(dis.var, aes(x = pop, y = median_rho)) +
    geom_point(data=subset(dis.var,county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 1.5, stroke = 0.5) +
    geom_point(data=subset(dis.var,county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
    scale_fill_manual(name = "MSA (large) and \ncounty (small)",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA (large) and \ncounty (small)",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name = "MSA (large) and \ncounty (small)",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    scale_x_log10(limits=range(dis.var$pop)) +
    labs(x = TeX('Population (millions), $\\log P$'),
         y = TeX('Spatial range exponent $\\rho$')) +
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.justification = c(0,0.5),
          legend.position = "right") +
    guides(alpha = "none")
  
  
  # B.msa[[yindex]] <- ggplot(dis.var, aes(x = pop, y = median_rho)) +
  #   geom_point(data=subset(dis.var,county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 2, stroke = 0.5) +
  #   geom_point(data=subset(dis.var,county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 4, stroke = 1) +
  #   scale_fill_manual(name="MSAs",
  #                     breaks = msa.info$dname,
  #                     labels = msa.info$dname,
  #                     values = scales::alpha(msa.info$col,alpha=0.8))+
  #   scale_colour_manual(name = "MSAs",
  #                       breaks = msa.info$dname,
  #                       labels = msa.info$dname,
  #                       values = msa.info$col)+
  #   scale_shape_manual(name="MSAs",
  #                      breaks = msa.info$dname,
  #                      labels = msa.info$dname,
  #                      values = msa.info$shape) +
  #   # scale_x_log10(limits=range(dis.var$pop)) +
  #   labs(x = TeX('Population (millions) $\\P$'),
  #        y = TeX('Spatial range exponent $\\rho$')) +
  #   theme_wy() +
  #   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
  #         legend.title = element_blank(),
  #         legend.position = "none")
} # yindex
fig.rho.fit <- (wrap_elements(A.msa[[2]])/wrap_elements(B.msa[[2]])) + plot_annotation(tag_levels = list(letters[4:5])) & theme(plot.tag = element_text(size = 20))&
  plot_layout(heights = c(1, 0.95))
ggsave(fig.rho.fit, filename = paste(figpath,"/msa/brief_rho_fit_",Yname[yindex],"_de.pdf",sep=""), width = 4*2, height = 4*2)


#----------Plot the fitting results: spearman rank correlation and map----------#
id.week<-paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")
Dname<<-Dname.set[1]
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif <- date.sep(Datevalue,durdate)
for(yindex in 2:2){
  
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
    
    #----------Fig1: map of flows----------# 
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
  
} # Yname
fig.rho.fit <- ((A.msa[[2]]/B.msa[[2]])|((C.msa[[1]]|D.msa[[1]])/(C.msa[[2]]|D.msa[[2]]))) + plot_layout(widths = c(3, 4)) + plot_annotation(tag_levels = list(letters[4:9])) & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.fit, filename = paste(figpath,"/msa/brief_rho_fit_",Yname[yindex],".pdf",sep=""), width = 4*4, height = 4*2)


