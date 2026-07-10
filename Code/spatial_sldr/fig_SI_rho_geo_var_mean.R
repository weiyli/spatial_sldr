# rm(list = ls())
# Mean-CBG-area version

# [ ρ vs travel distance ]   [ ρ vs travel mode ]   [ ρ vs area ]
# county level data from urban_sldr_fit.R; msa level data from sldr_fit.R; variable from rho_geo_var.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(sf)          # read_sf() 
library(spdep)       # poly2nb() https://blog.csdn.net/weixin_54000907/article/details/116247097
library(gridExtra)
library(ggtext)

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
Dname.set<-c(Dname.msa,Dname.dis)
d<<-1
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
Flow.name<-"Intra_Flow"



#----------rho and msa or county----------#
A.msa <- B.msa <- C.msa <- D.msa <- E.msa <- F.msa <- G.msa <- H.msa <- list()
#----------rho and msa----------#
id.week<-paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")

#----------Load msa-county / msa-tract relationship tables----------#
county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

#----------rho and var----------#
geo_var <- fread(paste0(datapath,"/msa/geo_var.csv",sep=""))
geo_var_county <- fread(paste0(datapath,"/county/geo_var.csv",sep=""))

#----------pop and area for county and msa----------#
queen.lag.msa <- fread(file.path(datapath, "msa", "queen_lag.csv"))
queen.lag.msa[, pop  := pop  / 1e6]  # million
queen.lag.msa[, area := area / 1e6]  # km^2
queen.max.msa <- unique(queen.lag.msa[, .(msa, nblock, lag.max, pop, area)])

#----------pop and area for county and msa----------#
queen.lag.county <- fread(file.path(datapath, "county", "queen_lag.csv"))
queen.lag.county[, pop  := pop  / 1e6]  # million
queen.lag.county[, area := area / 1e6]  # km^2
queen.max.county <- unique(queen.lag.county[, .(msa, nblock, lag.max, pop, area,county)])

#----------power-law scaling model: log(rho)=rho0+alpha*log(pop)----------#
for(yindex in 2:2){
  
  #----------msa level: data from sldr_fit.R and od_dist.R----------#
  rho.msa <- fread(paste0(datapath,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""))
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
  # distance and destination
  dist.msa.all<-NULL
  for(d in 1:Nmsa){
    Dname<<-Dname.set[d]
    dist <- fread(paste0(flowpath,"/msa/",Dname, "/dist_", var.msa[msa==Dname]$median_day, ".csv", sep=""))
    dist$msa <- Dname
    dist.msa.all <- plyr::rbind.fill(dist.msa.all, dist)
  }
  dist.msa.all$weighted_haversine <- dist.msa.all$weighted_haversine/1000  # m->km
  dist.msa <- as.data.table(dist.msa.all)[, .(
    mean_haversine     = mean(weighted_haversine, na.rm = TRUE),
    mean_destinations  = mean(n_destinations, na.rm = TRUE)
  ), by = msa]
  
  rho.dist <- left_join(var.msa, dist.msa, by="msa")
  
  #----------Connecting rho with urban indicators----------#
  rho_var <- left_join(rho.dist, geo_var, by="msa")
  rho_var[, c("pop_total", "area_km2") := NULL]
  cols_keep <- c(
    "msa",
    "county",
    "median_rho",
    "pop",
    "area",
    "avg_density",
    "mean_haversine",
    # "mean_destinations",
    "radial_slope",
    # "radial_r2",
    "compactness",
    "poly_gini",
    "rank_size_slope"
  )
  rho_sub <- rho_var[, ..cols_keep]
  rho_long <- melt(
    rho_sub,
    id.vars = c("msa","county", "median_rho"),
    variable.name = "indicator",
    value.name = "x"
  )
  rho_dt_msa <- as.data.table(rho_long)
  
  
  #----------county level: data from sldr_fit_county.R and od_dist.R----------#
  rho.county <- fread(paste0(datapath,"/county/SLDR_params_",Yname[yindex],".csv",sep=""))
  sldr.county <- rho.county[, .SD[rho == median(rho)], by = .(msa,period,index,county)]
  sldr.county <- sldr.county[index == "exp"&period==pervalue[1]]
  sldr.county <- sldr.county[, .(msa, day, rho, county)]

  rho.county <- rho.county[, .(
    min_rho = quantile(rho, probs = 0),
    Q1_rho = quantile(rho, probs = 0.25),
    median_rho = quantile(rho, probs = 0.5),
    Q3_rho = quantile(rho, probs = 0.75),
    max_rho = quantile(rho, probs = 1),
    median_day = day[which(rho==quantile(rho, probs = 0.5))]
  ), by = .(msa, index, period,county)]
  rho.county <- rho.county[index == "exp"&period==pervalue[1]]
  rho.county <- merge(rho.county, queen.max.county, by=c('msa','county'))
  var.county <- rho.county[, .(msa, median_day, median_rho, pop, area, county)]

  # distance and destination
  dist.county.all<-NULL
  for(d in 1:Nmsa){
    Dname<<-Dname.set[d]
    county.msa <- subset(county.info,msa_id==d)
    Ncounty <- nrow(county.msa)
    for(f in 1:Ncounty){
      county.id <- county.msa$county_fips[f]%>%as.numeric
      dist <- fread(file = paste0(flowpath,"/county/",Dname, "/dist_", var.county[msa==Dname&county==county.id]$median_day, ".csv", sep=""))
      dist$msa <- Dname
      dist.county.all <- plyr::rbind.fill(dist.county.all, dist)
    }
  }

  dist.county.all$weighted_haversine <- dist.county.all$weighted_haversine/1000  # m->km
  dist.county.all <- as.data.table(dist.county.all)
  setnames(dist.county.all, "county_fips", "county")
  dist.county <- dist.county.all[, .(
    mean_haversine    = mean(weighted_haversine, na.rm = TRUE),
    mean_destinations = mean(n_destinations, na.rm = TRUE)
  ), by = .(msa, county)]
  rho.dist.county <- left_join(var.county, dist.county, by=c("msa","county"))
  
  #----------Connecting rho with urban indicators----------#
  rho_var_county <- left_join(rho.dist.county, geo_var_county, by=c("msa","county"))
  rho_var_county[, c("pop_total", "area_km2") := NULL]
  cols_keep <- c(
    "msa",
    "county",
    "median_rho",
    "pop",
    "area",
    "avg_density",
    "mean_haversine",
    "compactness"
  )
  rho_sub <- rho_var_county[, ..cols_keep]
  rho_long <- melt(
    rho_sub,
    id.vars = c("msa", "county","median_rho"),
    variable.name = "indicator",
    value.name = "x"
  )
  rho_dt_county <- as.data.table(rho_long)
  rho_dt_county[, county := as.character(county)]
  rho_dt_county[, county := "county"]
  
  #----------row bind msa and county----------#
  rho_dt <- rbind(rho_dt_msa, rho_dt_county, fill = TRUE)
  
  stats_r2 <- rho_dt[is.finite(x) & is.finite(median_rho),
                     {
                       fit <- lm(median_rho ~ x)
                       s <- summary(fit)
                       r2 <- s$r.squared
                       
                       # p-value of slope term
                       p  <- coef(s)[2, "Pr(>|t|)"]
                       
                       # Pearson correlation
                       corr  <- suppressWarnings(cor(x, median_rho, method = "pearson", use = "complete.obs"))
                       
                       # choose an annotation position inside each facet
                       # (top-left corner)
                       x_pos <- min(x, na.rm = TRUE)
                       y_pos <- max(median_rho, na.rm = TRUE)
                       
                       list(
                         R2 = r2,
                         pval = p,
                         corr = corr,
                         x_pos = x_pos,
                         y_pos = y_pos
                       )
                     },
                     by = indicator]
  
  # format label (R2 + stars)
  stats_r2[, stars := fifelse(pval < 0.001, "***", fifelse(pval < 0.01,  "**", fifelse(pval < 0.05,  "*", "")))]
  # R2+ p-value：
  # stats_r2[, label := sprintf("R^2 = %.2f%s", R2, stars)]
  stats_r2[, label := sprintf(
    "R^2 == %.2f*'%s,'~'\n'~italic(r) == %.2g",
    R2, stars, corr
  )]
  
  # small nudges so text isn't exactly on the border
  stats_r2[, x_pos := x_pos + 0.03 * (max(rho_dt[indicator==.BY$indicator]$x, na.rm=TRUE) - 
                                        min(rho_dt[indicator==.BY$indicator]$x, na.rm=TRUE)),
           by = indicator]
  stats_r2[, y_pos := y_pos - 0.05 * (max(rho_dt$median_rho, na.rm=TRUE) - min(rho_dt$median_rho, na.rm=TRUE))]
  
  #---- 2) Plot with R^2 annotation ----#
  indicator_labels <- c(
    pop             = "Population~(million)",
    area            = "Urban~area~(km^2)",
    avg_density     = "Population~density",
    mean_haversine  = "Average~travel~distance~(km)",
    radial_slope    = "Radial~density~gradient",
    compactness     = "Compactness",
    poly_gini       = "Polycentricity~(Gini)",
    rank_size_slope = "Rank-size~slope"
  )
  
  A.msa[[yindex]] <- ggplot(rho_dt, aes(x = x, y = median_rho)) +
    # geom_point(aes(fill = msa, color = msa, shape = msa), size = 3, stroke = 0.5) +
    geom_point(data=subset(rho_dt, county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 1.5, stroke = 0.5) +
    geom_point(data=subset(rho_dt, county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
    geom_smooth(method = "lm", color = "#69b3a2", se = FALSE) +
    facet_wrap(~ indicator, scales = "free_x", nrow = 2,
               labeller = labeller(indicator = as_labeller(indicator_labels, label_parsed))) +
    geom_text(
      data = stats_r2,
      aes(x = x_pos, y = 0.8, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 5,
      parse = TRUE
    ) +
    scale_fill_manual(
      name = "MSA (large) and \ncounty (small)",
      breaks = msa.info$dname,
      labels = msa.info$dname,
      values = scales::alpha(msa.info$col, alpha = 0.8)
    ) +
    scale_colour_manual(
      name = "MSA (large) and \ncounty (small)",
      breaks = msa.info$dname,
      labels = msa.info$dname,
      values = msa.info$col
    ) +
    scale_shape_manual(
      name = "MSA (large) and \ncounty (small)",
      breaks = msa.info$dname,
      labels = msa.info$dname,
      values = msa.info$shape
    ) +
    labs(x = "Geographic/structural indicators",
         y = TeX("Spatial range exponent $\\rho$")) +
    theme_wy() +
    theme(
      panel.background = element_blank(),
      panel.spacing = unit(0.5, "lines"),
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
      strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
      strip.text = element_text(face = "plain", size = 15),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"
    ) +
    guides(alpha = "none")
} # yindex






#----------pop and area for county and msa----------#
queen.lag.county <- fread(file.path(datapath, "county", "queen_lag.csv"))
queen.lag.county[, pop  := pop  / 1e6]  # million
queen.lag.county[, area := area / 1e6]  # km^2
queen.max.county <- unique(queen.lag.county[, .(msa, county, nblock, lag.max, pop, area)])

queen.lag.msa <- fread(file.path(datapath, "msa", "queen_lag.csv"))
queen.lag.msa[, pop  := pop  / 1e6]  # million
queen.lag.msa[, area := area / 1e6]  # km^2
queen.max.msa <- unique(queen.lag.msa[, .(msa, nblock, lag.max, pop, area)])

#----------power-law scaling model: log(rho)=rho0+alpha*log(pop)----------#
for(yindex in 2:2){
  
  #----------msa level: data from sldr_fit.R and od_dist.R----------#
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
  # distance and destination
  dist.msa.all<-NULL
  for(d in 1:Nmsa){
    Dname<<-Dname.set[d]
    dist <- read.csv(paste(flowpath,"/msa/",Dname, "/dist_", var.msa[msa==Dname]$median_day, ".csv", sep=""), header=TRUE)
    dist$msa <- Dname
    dist.msa.all <- plyr::rbind.fill(dist.msa.all, dist)
  }
  dist.msa.all$weighted_haversine <- dist.msa.all$weighted_haversine/1000  # m->km
  dist.msa <- as.data.table(dist.msa.all)[, .(
    mean_haversine     = mean(weighted_haversine, na.rm = TRUE),
    mean_destinations  = mean(n_destinations, na.rm = TRUE)
  ), by = msa]
  
  rho.dist <- left_join(var.msa,dist.msa,by="msa")
  
  #----------SI msa level: density plot of travel distance and number of destinations----------#
  dist.mean <- dist.msa
  # dist.mean[, label_haversine := paste0("mean = ", round(mean_haversine, 1), " km")]
  dist.mean[, `:=`(
    label_haversine = paste0("mean = ", round(mean_haversine, 1), " km"),
    label_destinations = paste0("mean = ", round(mean_destinations))
  )]
  fig.dist.si <- ggplot(dist.msa.all, aes(x = weighted_haversine)) +
    geom_density(aes(fill = msa), alpha = 0.6, color = NA) +
    # geom_vline(data = dist.msa, aes(xintercept = mean_haversine,color = msa), linetype = "dashed") +
    geom_text(data = dist.mean, aes(x = mean_haversine, y = 0.2, label = label_haversine, color=msa),
              size = 4, hjust = -0.1) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_color_manual(name= "MSAs",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$col) +
    scale_fill_manual(name= "MSAs",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
    xlim(0, 50) +
    labs(# title = "Density of weighted OD travel distance across MSAs",
      x = "Weighted average OD travel distance (km)",
      y = "Density") +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          strip.text = element_text(face = "plain"),     
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          legend.position = "none")
  fig.ndest.si <- ggplot(dist.msa.all, aes(x = n_destinations)) +
    geom_density(aes(fill = msa), alpha = 0.6, color = NA) +
    geom_text(data = dist.mean, aes(x = mean_destinations, y = 0.015, label = label_destinations, color=msa),
              size = 4, hjust = -0.11) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_color_manual(name= "MSAs",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$col) +
    scale_fill_manual(name= "MSAs",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
    xlim(0, 250) +
    labs(x = "Number of destinations", y = "Density") +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          strip.text = element_text(face = "plain"),     
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          legend.position = "none")
  fig.dist.ndest.si <- (fig.dist.si/fig.ndest.si)+ plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.dist.ndest.si, filename = paste(figpath,"/msa/SI_dist_ndest_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*4)
  ggsave(fig.dist.si, filename = paste(figpath,"/msa/SI_dist_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*2)
  
  
  #----------rho and distance and destinations----------#
  D.msa[[yindex]] <- ggplot(rho.dist, aes(x = mean_haversine, y = median_rho)) +
    geom_point(aes(fill= msa, color = msa, shape = msa), size=3,stroke = 0.5) +
    geom_smooth(method = "lm", color = "#69b3a2") + 
    scale_fill_manual(name = 'MSA',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = 'MSA',
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col) +
    scale_shape_manual(name = 'MSA',
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    labs(x = "Average travel distance (km)", 
         y = TeX('Spatial range exponent $\\rho$')) + 
    theme_wy() +
    theme(panel.border = element_blank(),
          legend.position = "none",
          axis.line.x = axis.arrow,    
          axis.line.y = axis.arrow)
  
  E.msa[[yindex]] <- ggplot(rho.dist, aes(x = mean_destinations, y = median_rho)) +
    geom_point(aes(fill= msa, color = msa, shape = msa), size=3,stroke = 0.5) +
    geom_smooth(method = "lm", color = "#69b3a2") + 
    scale_fill_manual(name = 'MSA',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = 'MSA',
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col) +
    scale_shape_manual(name = 'MSA',
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    labs(x = "Average number of destinations", 
         y = TeX('Spatial range exponent $\\rho$')) + 
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          axis.line.x = axis.arrow,    
          axis.line.y = axis.arrow)
  
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
  rho.pop.model <- lm(log(dis.var$median_rho) ~ dis.var$pop)
  rho.pop.coef <- coef(rho.pop.model)
  rho.pop.coef[1] <- exp(rho.pop.coef[1])
  rho0 <- round(rho.pop.coef[1],2)
  alpha <- round(rho.pop.coef[2],2)
  R2 <- summary(rho.pop.model)$r.squared%>%round(2)
  
  rho.pop.pval <- coef(summary(rho.pop.model))[, "Pr(>|t|)"]
  intercept.signif <- get.significance(rho.pop.pval[1])
  pop.signif <- get.significance(rho.pop.pval[2])
  
  #----------fitted rho for msa and county----------#
  sldr.rho <- left_join(dis.sldr,dis.var,by=c("msa","county"))
  # write.csv(sldr.rho, file=paste(datapath,"/msa/rho_pop_",Yname[yindex],".csv",sep=""), row.names = FALSE)
  
  #----------plot linear fitting of rho and pop----------#
  C.msa[[yindex]] <- ggplot(dis.var, aes(x = pop, y = median_rho)) +
    geom_point(data=subset(dis.var,county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 2, stroke = 0.5) +
    geom_point(data=subset(dis.var,county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 4, stroke = 1) +
    geom_smooth(method = 'lm', color = "#69b3a2") +
    scale_fill_manual(name="MSA: counties",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA: counties",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA: counties",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    # scale_x_log10(limits=range(dis.var$pop)) +
    # scale_y_log10(limits=range(dis.var$median_rho)) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = TeX('Population (millions), $\\log P$'),
         y = TeX('Spatial range exponent, $\\log \\rho$')) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(),
          legend.position = "none") +
    annotation_custom(grob = textGrob(TeX(sprintf("$\\rho = %.2f P^{%.2f}$", rho0, alpha)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = 0.6, xmax = 0.6, ymin = 0.8, ymax = 0.8) +
    annotation_custom(grob = textGrob(TeX(sprintf("$R^2 = %.2f,$", R2)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = 0.6, xmax = 0.6, ymin = 0.7, ymax = 0.7) +
    annotation_custom(grob = textGrob(TeX(sprintf("$\\italic(p) < %.3f$", 0.001)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = 2.5, xmax = 2.5, ymin = 0.7, ymax = 0.7)
  
  
  
  #----------nonlinear fitting rho and pop for county and msa----------#
  nls.model <- nls(median_rho ~ a * exp(-b * pop) + c, data = dis.var, start = list(a = 0.5, b = 0.1, c = 0.01))
  dis.var$model_rho <- predict(nls.model)
  R2 <- R2_test(y_test=dis.var$median_rho,y_predict=dis.var$model_rho)
  nls.coef <- coef(nls.model)
  a <- round(nls.coef[1], 2)
  b <- round(nls.coef[2], 2)
  c <- round(nls.coef[3], 2)
  
  H.msa[[yindex]] <- ggplot(dis.var, aes(x = pop, y = median_rho)) +
    geom_point(data=subset(dis.var,county!="msa"),aes(fill = msa, color = msa, shape = msa, alpha=0.1), size= 2, stroke = 0.5) +
    geom_point(data=subset(dis.var,county=="msa"),aes(fill = msa, color = msa, shape = msa), size= 4, stroke = 1) +
    geom_line(data = dis.var, aes(x = pop, y= model_rho), color = "#69b3a2", linewidth = 1.2) +
    scale_fill_manual(name="MSA: counties",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "MSA: counties",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="MSA: counties",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = TeX('Population (millions), $\\log P$'),
         y = TeX('Spatial range exponent, $\\log \\rho$')) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.title = element_blank(),
          legend.position = "none") +
    annotation_custom(grob = textGrob(TeX(sprintf("$\\rho = %.2f e^{-%.2f P} + %.2f$", a, b, c)), 
                                      gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = -0.5, xmax = -0.5, ymin = -0.1, ymax = -0.1) +
    annotation_custom(grob = textGrob(TeX(sprintf("$R^2 = %.2f$", R2)), gp = gpar(fontsize = 12, col = "#69b3a2"), hjust = 0), xmin = -0.5, xmax = -0.5, ymin = -0.13, ymax = -0.13)
  
  
} # yindex




#----------rho and travel mode for msa----------#
travel.mode <- read.csv(paste(geopath,"/census/acs_2019_5years_travel_mode_msa.csv",sep=""), header=TRUE)
travel.mode<-left_join(travel.mode, data.frame(msa_fips=msa.info$did, msa=msa.info$dname), by='msa_fips')
travel.mode <- left_join(queen.max.msa,travel.mode,by='msa')
#----------plot----------#
for(yindex in 2:2){
  
  dis.rho <- read.csv(paste(datapath,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  dis.rho <- dis.rho[, .(
    min_rho = min(rho),
    median_rho = median(rho),
    max_rho = max(rho)
  ), by = .(msa, index, period)]
  dis.rho <- dis.rho[index == "exp"&period==pervalue[1]]
  
  #----------travel mode----------#
  rho.mode <- left_join(data.frame(msa=dis.rho$msa,median_rho=dis.rho$median_rho),travel.mode,by=c("msa"))
  index <- colnames(rho.mode)%in%c("msa","pop","median_rho","drive_alone","carpool","public_transit","bicycle","walk","work_at_home")
  rho.mode<-rho.mode[,index]
  for(i in 4:ncol(rho.mode)){
    rho.mode[,i] <- rho.mode[,i]/rho.mode$pop
    # rho.mode[,i] <- rho.mode[,i]/rho.mode$total_commute
  }
  
  #----------correlation between rho and travel mode,time----------#
  var.name<-c("drive_alone","carpool","public_transit","bicycle","walk","work_at_home")
  var.name <- stringr::str_to_sentence(gsub("_", " ", var.name))
  var.group <- rep('Travel mode',8)
  var.corr<-NULL
  for(i in 1:length(var.name)){
    pearson_corr <- cor(rho.mode$median_rho, rho.mode[,3+i], method = "pearson", use = "complete.obs")
    spearman_corr <- cor(rho.mode$median_rho, rho.mode[,3+i], method = "spearman", use = "complete.obs")
    var.corr <- plyr::rbind.fill(var.corr, data.frame(corr_value=c(pearson_corr,spearman_corr),
                                                      corr_method=c("Pearson","Spearman"),
                                                      var_name=var.name[i],
                                                      var_group=var.group[i]))
  }
  #----------plot travel mode for pearson----------# 
  mode.corr <- subset(var.corr,corr_method=="Pearson"&var_group=="Travel mode")
  max_corr <- max(abs(mode.corr$corr_value), na.rm = TRUE)
  F.msa[[yindex]] <- ggplot(mode.corr, aes(x = corr_value, y = reorder(var_name, corr_value), 
                                           fill = corr_value)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_x_continuous(limits = c(-max_corr, max_corr), breaks = seq(-1, 1, by = 0.2)) +
    scale_fill_gradient2(low = "#698CC7", mid = "white", high = "#C7698C", midpoint = 0) +
    # facet_grid( ~ corr_method, scales = "free", space = "free") +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    labs(
      # title = 'Spatial range exponent vs. Population share by travel mode ',
      # TeX('Population share by travel mode vs. spatial range exponent $\\rho_{mob}$'), 
      x = "Pearson's correlation coefficient", y = "Travel mode", fill = "Correlation") +
    theme_wy() +
    theme(panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0.5), 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "gray", linetype = "dashed"),
          legend.position = "none",
          axis.ticks = element_blank())
  
  #----------plot travel mode for spearman----------# 
  mode.corr <- subset(var.corr,corr_method=="Spearman"&var_group=="Travel mode")
  max_corr <- max(abs(mode.corr$corr_value), na.rm = TRUE)
  G.msa[[yindex]] <- ggplot(mode.corr, aes(x = corr_value, y = reorder(var_name, corr_value), 
                                           fill = corr_value)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_x_continuous(limits = c(-max_corr, max_corr), breaks = seq(-1, 1, by = 0.2)) +
    scale_fill_gradient2(low = "#698CC7", mid = "white", high = "#C7698C", midpoint = 0) +
    # facet_grid(var_group ~ ., scales = "free", space = "free") +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    labs(title = 'Spatial range exponent vs. Population share by travel mode ',
         # TeX('Population share by travel mode vs. spatial range exponent $\\rho_{mob}$'), 
         x = "Spearman's rank correlation coefficient", y = "Travel mode", fill = "Correlation") +
    theme_wy() +
    theme(panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0.5), 
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line(color = "gray", linetype = "dashed"),
          axis.ticks = element_blank())
  ggsave(G.msa[[yindex]], filename = paste(figpath,"/msa/SI_rho_travel_mode_",Yname[yindex],".pdf",sep=""), width = 4*2, height = 3*2)
  
}

#----------Model performance for Atlanta and Boston----------#
yindex<-2
fig.rho.var <- ((A.msa[[yindex]]) / (C.msa[[yindex]] | F.msa[[yindex]])) + plot_layout(heights = c(1.6, 1), widths = c(1, 1)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.var,filename = paste(figpath,"/msa/SI_rho_var_",Yname[yindex],".pdf",sep=""), width =4*4, height = 4.5*3)


#----------MSA-only bivariate relationships for TO_in_OO----------#
yindex <- 2
rho.msa <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], ".csv"))
rho.msa <- rho.msa[, .(
  min_rho = quantile(rho, probs = 0),
  Q1_rho = quantile(rho, probs = 0.25),
  median_rho = quantile(rho, probs = 0.5),
  Q3_rho = quantile(rho, probs = 0.75),
  max_rho = quantile(rho, probs = 1),
  median_day = day[which(rho == quantile(rho, probs = 0.5))[1]]
), by = .(msa, index, period)]
rho.msa <- rho.msa[index == "exp" & period == pervalue[1]]
rho.msa <- merge(rho.msa, queen.max.msa, by = c("msa"))
var.msa <- rho.msa[, .(msa, day = median_day, median_rho, pop, area)]
var.msa[, county := "msa"]

#----------distance and destination (daily, msa-level)----------#
dist.msa.all <- NULL
for (d in 1:Nmsa) {
  Dname <- Dname.set[d]
  days_msa <- unique(var.msa[msa == Dname]$day)
  
  for (dd in 1:length(days_msa)) {
    fpath <- paste0(flowpath, "/msa/", Dname, "/dist_", days_msa[dd], ".csv")
    if (!file.exists(fpath)) next
    
    dist <- fread(fpath)
    dist[, msa := Dname]
    dist[, day := days_msa[dd]]
    
    dist.msa.all <- plyr::rbind.fill(dist.msa.all, dist)
  }
}
dist.msa.all <- as.data.table(dist.msa.all)
dist.msa.all[, weighted_haversine := weighted_haversine / 1000]
dist.msa <- dist.msa.all[, .(
  mean_haversine = mean(weighted_haversine, na.rm = TRUE),
  mean_destinations = mean(n_destinations, na.rm = TRUE)
), by = .(msa, day)]
rho.dist <- left_join(var.msa, dist.msa, by = c("msa", "day"))

rho_var_msa <- left_join(rho.dist, geo_var, by = "msa")
travel.mode.msa <- read.csv(paste(geopath, "/census/acs_2019_5years_travel_mode_msa.csv", sep = ""), header = TRUE)
travel.mode.msa <- left_join(
  travel.mode.msa,
  data.frame(msa_fips = msa.info$did, msa = msa.info$dname),
  by = "msa_fips"
)
rho_var_msa <- left_join(rho_var_msa, travel.mode.msa, by = "msa")
vars_norm_msa <- c(
  "total_commute",
  "drive_alone",
  "carpool",
  "public_transit",
  "bicycle",
  "walk",
  "other_means",
  "work_at_home"
)
vars_norm_msa <- vars_norm_msa[vars_norm_msa %in% names(rho_var_msa)]
rho_var_msa[, (vars_norm_msa) := lapply(.SD, function(x) x / pop_total), .SDcols = vars_norm_msa]
rho_var_msa[, mean_cbg_area := mean_cbg_area_km2]
rho_var_msa[, c("pop_total", "area_km2") := NULL]
cols_keep_msa <- c(
  "msa",
  "county",
  "median_rho",
  "pop",
  "area",
  "mean_cbg_area",
  "mean_haversine",
  "drive_alone",
  "public_transit",
  "radial_slope",
  "poly_gini",
  "rank_size_slope"
)
rho_sub_msa <- rho_var_msa[, ..cols_keep_msa]
rho_dt_msa_only <- melt(
  rho_sub_msa,
  id.vars = c("msa", "county", "median_rho"),
  variable.name = "indicator",
  value.name = "x"
)
rho_dt_msa_only <- as.data.table(rho_dt_msa_only)

indicator_labels_msa <- c(
  pop = "Population~(million)",
  area = "Urban~area~(km^2)",
  mean_cbg_area = "Mean~CBG~area~(km^2)",
  mean_haversine = "Average~travel~distance~(km)",
  drive_alone = "Drive~alone",
  public_transit = "Public~transit",
  radial_slope = "Radial~density~gradient",
  poly_gini = "Polycentricity~(Gini)",
  rank_size_slope = "Rank-size~slope"
)
indicator_order_msa <- names(indicator_labels_msa)
rho_dt_msa_only <- rho_dt_msa_only[indicator %in% indicator_order_msa]
rho_dt_msa_only[, indicator := factor(indicator, levels = indicator_order_msa)]

stats_r2_msa <- rho_dt_msa_only[is.finite(x) & is.finite(median_rho), {
  fit <- lm(median_rho ~ x)
  s <- summary(fit)
  pval <- coef(s)[2, 4]
  corr <- suppressWarnings(cor(x, median_rho, method = "pearson", use = "complete.obs"))
  .(
    R2 = s$r.squared,
    pval = pval,
    corr = corr,
    x_pos = min(x, na.rm = TRUE),
    y_pos = max(median_rho, na.rm = TRUE)
  )
}, by = indicator]
stats_r2_msa[, stars := fifelse(pval < 0.001, "***",
                                fifelse(pval < 0.01, "**",
                                        fifelse(pval < 0.05, "*", "")))]
stats_r2_msa[, label := sprintf(
  "R^2 == %.2f*'%s,'~italic(r) == %.2g",
  R2, stars, corr
)]
stats_r2_msa[, x_pos := x_pos + 0.03 * (
  max(rho_dt_msa_only[indicator == .BY$indicator]$x, na.rm = TRUE) -
    min(rho_dt_msa_only[indicator == .BY$indicator]$x, na.rm = TRUE)
), by = indicator]
stats_r2_msa[, y_pos := y_pos - 0.05 * (
  max(rho_dt_msa_only$median_rho, na.rm = TRUE) -
    min(rho_dt_msa_only$median_rho, na.rm = TRUE)
)]

fig.rho.var.msa <- ggplot(rho_dt_msa_only, aes(x = x, y = median_rho)) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 3, stroke = 0.5) +
  geom_smooth(method = "lm", color = "#69b3a2", se = FALSE) +
  facet_wrap(
    ~ indicator,
    scales = "free_x",
    nrow = 3,
    labeller = labeller(indicator = as_labeller(indicator_labels_msa, label_parsed))
  ) +
  geom_text(
    data = stats_r2_msa,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 5,
    parse = TRUE
  ) +
  scale_fill_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col, alpha = 0.8)
  ) +
  scale_colour_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$col
  ) +
  scale_shape_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$shape
  ) +
  labs(
    x = "Geographic/structural indicators",
    y = TeX("Spatial range exponent $\\rho$")
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

fig.rho.var.msa <- fig.rho.var.msa +
  plot_annotation(tag_levels = list(letters[3])) &
  theme(plot.tag = element_text(size = 20))

ggsave(
  fig.rho.var.msa,
  filename = paste(figpath, "/msa/SI_rho_var_msa_", Yname[yindex], "_mean.pdf", sep = ""),
  width = 4 * 4,
  height = 4.5 * 3
)

