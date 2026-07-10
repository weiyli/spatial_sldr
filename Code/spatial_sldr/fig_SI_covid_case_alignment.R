# rm(list = ls())


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath.msa <- 'D:/ood/Data/spatial_sldr'
datapath.dis <- 'D:/ood/Data/spatial_sldr/disaster'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Load packages----------#
library(ggrepel)     # geom_text_repel()
library(stringr)

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
Dname<<-Dname.set[d]
date.gif<-date.sep(Datevalue,durdate)


#----------confirmed cases for COVID-19----------#
con.us <- read.csv(paste(datapath.dis,"/us-counties-daily-cumulative-case.csv",sep=""), header=TRUE)
con.us$fips <- str_pad(con.us[,1], width = 5, pad = "0")
con.date <- as.Date(sub("X", "", colnames(con.us[,-1])), format = "%Y%m%d")
cases <- NULL
for(d in 1:Nmsa){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  date.gif<-date.sep(Datevalue,durdate)
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  date.gif <- date.sep(Datevalue,durdate)
  
  #----------BlockID, NBlock----------#
  block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.msa$CountyID <- str_pad(block.msa$CensusBlockGroup, width = 12, pad = "0")
  block.msa$CountyID <- substring(block.msa$CountyID,first=1,last=5)
  CountyID <- unique(block.msa$CountyID)
  con.county <- subset(con.us,fips%in%CountyID)
  con <- data.frame(day=con.date, case=colSums(con.county[,-1]))
  rownames(con) <- NULL
  con <- con %>% arrange(day) %>% mutate(new_case = case - lag(case, default = first(case)))
  con <- subset(con, day%in%Datevalue)
  con$msa <- Dname
  cases <- plyr::rbind.fill(cases, con)
}
write.csv(cases,file=paste(datapath.dis,"/us-msas-daily-cumulative-case.csv",sep=""),row.names = FALSE)


#----------percent change of median rho and confirmed cases for COVID-19----------#
con.msa <- read.csv(paste(datapath.dis,"/us-msas-daily-cumulative-case.csv",sep=""), header=TRUE)
queen.lag.msa <- read.csv(paste(datapath.msa,"/msa/queen_lag.csv",sep=""), header=TRUE)
pop.msa <- queen.lag.msa[,c("msa","pop")]%>%distinct
con.msa <- left_join(con.msa,pop.msa,by="msa")
# plot rho of pop and mob for different disasters
A.msa <- B.msa <- C.msa <- D.msa <- list()
for(yindex in 2:2){
  
  #----------daily percent change of all rho----------#
  # daily rho
  rho.mob.msa <- read.csv(paste(datapath.msa,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  id.week<-data.frame(day=as.character(Datevalue),week.id=factor(rep(1:(Day/2),2)),week.name=format(Datevalue, "%a"))
  daily.rho<-left_join(rho.mob.msa[index == "exp"],id.week,by="day")
  daily.rho.diff <- daily.rho[, .(det_rho = rho[period == "during"]-rho[period == "before"],
                                  percent_det_rho = (rho[period == "during"] - rho[period == "before"])/rho[period == "during"]), by = .(msa, week.id, week.name)]
  # case and rho
  con.msa$day<-as.Date(con.msa$day)
  case.msa<-left_join(con.msa,date.gif,by="day")
  case.msa<-subset(case.msa,period=="during")
  case.msa$day<-as.character(case.msa$day)
  case.msa<-left_join(case.msa,id.week,by="day")
  rho.case.msa <- left_join(daily.rho.diff,case.msa,by=c("msa","week.id","week.name" ))
  rho.case.msa$percent_case <- rho.case.msa$case/rho.case.msa$pop
  rho.case.msa$percent_new_case <- rho.case.msa$new_case/rho.case.msa$pop
  rho.case.median <-  setDT(rho.case.msa)[, .SD[percent_det_rho == median(percent_det_rho)], by = .(msa)]
  cor(rho.case.msa[,c("percent_det_rho","case","new_case","percent_case","percent_new_case")])
  pearson.cor <- cor(rho.case.msa$percent_det_rho,log(rho.case.msa$case/rho.case.msa$pop,10), method="pearson")
  spearman.cor <- cor(rho.case.msa$percent_det_rho,log(rho.case.msa$case/rho.case.msa$pop,10), method = "spearman")
  
  #----------Bin x (severity) and average y within bins, then compute correlation----------#
  # equal-width bins in log10-space
  rho.case.msa$log_percent_case <- log(rho.case.msa$percent_case,10)
  nbin <- 15
  bin_breaks <- seq(min(rho.case.msa$log_percent_case, na.rm = TRUE), max(rho.case.msa$log_percent_case, na.rm = TRUE), length.out = nbin + 1)
  rho.case.bin <- as.data.table(rho.case.msa)
  rho.case.bin[, xbin := cut(log_percent_case, breaks = bin_breaks, include.lowest = TRUE, right = FALSE)]
  # mean(y) within each x-bin (and keep a representative x value for correlation)
  rho.case.bin_sum <- rho.case.bin[!is.na(xbin), .(
    x_log_mean = mean(log_percent_case, na.rm = TRUE),   # representative x (log10)
    x_log_mid  = (min(log_percent_case, na.rm = TRUE) + max(log_percent_case, na.rm = TRUE)) / 2,
    y_mean     = mean(percent_det_rho, na.rm = TRUE),
    n          = .N
  ), by = xbin][n >= 3]   # drop sparse bins
  # Pearson correlation on binned means
  pearson.cor.bin <- cor(rho.case.bin_sum$y_mean, rho.case.bin_sum$x_log_mean, method = "pearson")
  spearman.cor.bin <- cor(rho.case.bin_sum$y_mean, rho.case.bin_sum$x_log_mean, method = "spearman")
  rho.case.bin_sum[, x_mean := 10^x_log_mean]   # If you prefer to report/plot x back on original scale:
  
  
  #--------------------------------------#
  # Correlation (raw vs binned) for annotation
  #--------------------------------------#
  # raw scatter: x is log10(percent_case), y is percent_det_rho
  dt_raw <- as.data.table(rho.case.msa)
  dt_raw[, x_raw := log10(percent_case)]
  dt_raw <- dt_raw[is.finite(x_raw) & is.finite(percent_det_rho)]
  
  pearson_raw  <- suppressWarnings(cor(dt_raw$percent_det_rho, dt_raw$x_raw, method = "pearson"))
  spearman_raw <- suppressWarnings(cor(dt_raw$percent_det_rho, dt_raw$x_raw, method = "spearman"))
  
  # binned means: x is bin mean in log10-space, y is bin mean
  dt_bin <- as.data.table(rho.case.bin_sum)
  dt_bin <- dt_bin[is.finite(x_log_mean) & is.finite(y_mean)]
  
  pearson_bin  <- suppressWarnings(cor(dt_bin$y_mean, dt_bin$x_log_mean, method = "pearson"))
  spearman_bin <- suppressWarnings(cor(dt_bin$y_mean, dt_bin$x_log_mean, method = "spearman"))
  
  # anchor position: top-right inside panel
  x_pos <- max(dt_raw$percent_case, na.rm = TRUE)*0.155
  y_pos <- max(dt_raw$percent_det_rho, na.rm = TRUE)
  
  # vertical spacing
  y_range <- range(dt_raw$percent_det_rho, na.rm = TRUE)
  dy <- 0.06 * diff(y_range)

  ann_cor <- data.table::rbindlist(list(
    data.table(line = 1, x = x_pos, y = y_pos,
               label = sprintf("Raw~~Pearson~~italic(r)==%.2f",  pearson_raw)),
    data.table(line = 2, x = x_pos, y = y_pos - dy,
               label = sprintf("Raw~~Spearman~~italic(r)[S]==%.2f", spearman_raw)),
    data.table(line = 3, x = x_pos, y = y_pos - 2*dy,
               label = sprintf("Binned~~Pearson~~italic(r)==%.2f", pearson_bin)),
    data.table(line = 4, x = x_pos, y = y_pos - 3*dy,
               label = sprintf("Binned~~Spearman~~italic(r)[S]==%.2f", spearman_bin))
  ), use.names = TRUE)
  
  
  #----------plot the percent change of rho and cum confirmed cases, then add binned-mean points on top of the original scatter----------#
  A.msa[[yindex]] <- ggplot(rho.case.msa, aes(x = percent_case, y = percent_det_rho)) +
    geom_point(aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
    geom_point(data = rho.case.bin_sum, aes(x = x_mean, y = y_mean),
               inherit.aes = FALSE, shape = 24, size = 5, stroke = 1, fill = "#69b3a2", color = "#69b3a2", alpha = 1) +
    geom_text(data = ann_cor, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, 
              size = 4, color = "#69b3a2", parse = TRUE) + 
    scale_fill_manual(name="COVID-19",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha = 0.5))+
    scale_colour_manual(name = "COVID-19",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = scales::alpha(msa.info$col,alpha = 0.5))+
    scale_shape_manual(name="COVID-19",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    labs(x = TeX('Cumulative confirmed cases/population ($\\log$)'),
         y = TeX('The percentage change of $\\rho$')) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
    theme_wy() +
    theme(panel.border = element_rect(fill=NA, color="black", linewidth=0.5, linetype="solid"),
          panel.background = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "right",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())

  ggsave(A.msa[[yindex]], filename = paste(figpath,"/msa/SI_rho_case_",Yname[yindex],".pdf",sep=""), width = 7.2, height = 5)
  
} # yindex



#----------plot SIR results----------#
E.msa <- F.msa <- G.msa <- list()
yindex <- 2
#----------data from the spreading_model.R----------#
sir.model <- read.csv(paste(datapath.msa,"/msa/SIR_model_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.model$DM <- "Model"
sir.model[which(sir.model$det.rho==0),]$DM<-"Empirical"
sir.model.main <- subset(sir.model,beta==0.3&mu==0.1&msa%in%c("Atlanta","San Francisco"))
#----------plot infection peak and rho----------#
sir.model.peak <- as.data.table(sir.model.main)[, .SD[which.max(infected)], by = .(msa, day, rho, det.rho)]
baseline <- sir.model.peak[det.rho == 0, .(msa, day, time.base = time, infected.base = infected)]
sir.model.peak <- merge(sir.model.peak, baseline, by = c("msa", "day"), all.x = TRUE)
sir.model.peak[, ':='(
  time.pct =  (time - time.base) / time.base,
  infected.pct = (infected - infected.base) / infected.base
)]
E.msa[[yindex]] <- ggplot(sir.model.peak, aes(det.rho,round(infected.pct,2))) +  
  geom_line(aes(color = msa, group = msa), linetype = "dotted", na.rm = TRUE) +
  geom_point(aes(fill= msa, color = msa, shape = msa),size=3,stroke = 0.5) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
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
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.2),
                     labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.2),
                     labels = scales::percent_format(scale = 100)) +
  labs(x=TeX('Percentage change in $\\rho_{mob}$'), y='Percentage change in infection peak')+
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.title = element_blank(), 
        legend.key.size = unit(0.5, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.05,0.85),
        legend.text = element_text(hjust=0),
        legend.justification = c(0,0.5),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())
#----------plot infected pop and time----------#
sir.model.time <- subset(sir.model.main,det.rho%in%c(-0.5,0,0.5))
det.rho <- sort(unique(sir.model.time$det.rho),decreasing=TRUE)
rho.breaks <- rho.labels <- paste0(round(det.rho,4)*100,"%")
rho.labels[which(rho.breaks=="0%")] <- "0% empirical"
rho.info <- data.frame(det.rho = det.rho, 
                       breaks = rho.breaks,
                       labels = rho.labels,
                       color = col.fit[c(5,3,1)],
                       line = c("dotted","solid","dotdash"))
# color = col.all[1:length(det.rho)])
sir.model.time <- left_join(sir.model.time, rho.info, by="det.rho")
emp.rho <- as.data.table(sir.model.time)[DM == "Empirical", .SD[which.max(infected)], by = msa]
#----------plot infected rate for two msas with (beta=0.3,mu=0.1)----------#
F.msa[[yindex]] <- ggplot(data = sir.model.time, aes(x = time, y = infected/pop)) +
  geom_line(aes(color = labels, linetype = labels, group = rho), linewidth = 0.8) +
  geom_label_repel(data = emp.rho,
                   aes(label = paste("empirical",round(rho, 4),sep=" ")), 
                   # aes(label = sprintf("rho == %.4f", round(rho, 4))), 
                   max.overlaps = Inf,
                   box.padding = 1, size = 3, fill = 'white',color="black") +
  facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  scale_colour_manual(name = TeX('% change in $\\rho_{mob}$'),
                      breaks = rho.info$labels,
                      labels = rho.info$labels,
                      values = rho.info$col) +
  scale_linetype_manual(name = TeX('% change in $\\rho_{mob}$'), 
                        breaks = rho.info$labels,
                        labels = rho.info$labels,
                        values = rho.info$line) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(face = "plain"),     
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
        legend.position = "right",
        legend.justification = c(0,0.5))
fig.rho.sir <- ((A.msa[[yindex]]|B.msa[[yindex]])/C.msa[[yindex]]/(E.msa[[yindex]]|F.msa[[yindex]]))+ 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/rho_pop_mob_sir_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*3)
#----------plot infected rate for each msa with all (beta,mu)----------#
dset <- unique(sir.model$msa)%>%sort
dset <- c("Atlanta","San Francisco")
sir.model$beta_lab <- paste0("beta==", sir.model$beta)
sir.model$mu_lab   <- paste0("mu==", sir.model$mu)
for(d in 1:length(dset)){
  G.msa[[d]] <- ggplot(data = subset(sir.model,msa==dset[d]), aes(x = time, y = infected/pop)) +
    geom_line(aes(color = rho, group = rho), linewidth = 0.8) +
    facet_grid(mu_lab ~ beta_lab, scales = "fixed", labeller = label_parsed) +
    scale_colour_gradientn(name = dset[d], colors = col.fit) +
    labs(x = "Time", y = "Infected fraction", title = NULL) +
    xlim(0, 400) +  
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain", size = 15),
          axis.text.y = element_text(hjust = 0),
          legend.position = "right",
          legend.justification = c(0, 0.5))
}
fig.rho.sir <- (G.msa[[1]]/G.msa[[2]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/SI_sir_beta_mu_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*4)












#----------plot SIR results----------#
E.msa <- F.msa <- list()
yindex <- 2
#----------data from the spreading_data.R----------#
sir.data <- read.csv(paste(datapath.msa,"/msa/SIR_data_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.data.main <- subset(sir.data, beta==0.3&mu==0.1&msa%in%c("Atlanta","San Francisco"))
sir.data.main$DM <-"Empirical"
#----------data from the spreading_model.R----------#
sir.model <- read.csv(paste(datapath.msa,"/msa/SIR_model_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.model <- subset(sir.model,det.rho<=0)
# sir.model <- sir.model %>% group_by(msa) %>% filter(rho == max(rho) | rho == min(rho))
sir.model.main <- subset(sir.model, beta==0.3&mu==0.1&msa%in%c("Atlanta","San Francisco"))
sir.model.main$DM <-"Model"
sir.all <- plyr::rbind.fill(sir.model.main, sir.data.main)
sir.all$day <- as.Date(sir.all$day)
sir.all <- left_join(sir.all,date.gif,by="day")

sir1 <- as.data.table(subset(sir.all,det.rho == 0))[, .(msa, rho, det.rho, DM, period)]%>%unique
sir2 <- sir1[, .(
  det.rho.new = (rho[period == "during"]-rho[period == "before"])/rho[period == "before"] 
), by = .(msa, det.rho, DM)]
sir3 <- left_join(sir.all,sir2,by=c("msa","det.rho","DM"))%>%as.data.table
sir3[det.rho == 0 & period == 'during', det.rho := det.rho.new]
sir.all$det.rho <- sir3$det.rho

det.rho <- sort(unique(sir.all$det.rho),decreasing=TRUE)
rho.info <- data.frame(det.rho = det.rho, 
                       labels = paste0(round(det.rho,4)*100,"%"),
                       color= col.all[1:length(det.rho)])
sir.all<-left_join(sir.all,rho.info,by="det.rho")

max.rho <- sir.all %>% group_by(msa) %>%  filter(rho == max(rho, na.rm = TRUE)) %>% slice(1) 
#----------plot infected rate for two msas with one (beta,mu)----------#
E.msa[[yindex]] <- ggplot(data = sir.all, aes(x = time, y = infected/pop)) +
  geom_line(aes(color = labels, linetype=DM, group = interaction(rho, DM)), linewidth = 0.8) +
  # geom_label_repel(data = max.rho, aes(label = round(rho, 4), x = time, y = infected / pop),
  #   color = "black",fill = "white", size = 4, fontface = "bold", box.padding = 0.3, 
  #   point.padding = 0.3, segment.color = "grey50"
  # ) +
  geom_label_repel(data = max.rho,
                   aes(label = round(rho, 4)), max.overlaps = Inf,
                   box.padding = 1, size = 3, fill = 'white',color="black") +
  facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  scale_colour_manual(name= "% change",
                      breaks = rho.info$labels,
                      labels = rho.info$labels,
                      values = rho.info$col) +
  scale_linetype_manual(name=TeX('$\\rho$'),
                        breaks = DM.info$breaks,
                        labels = DM.info$labels,
                        values = DM.info$line) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  # xlim(0, 400) +
  # theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        panel.spacing = unit(0.1, "lines"),
        text = element_text(size = 15, family = "Helvetica Neue"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
        strip.text = element_text(face = "plain",size=15),
        axis.text.y = element_text(hjust = 0),
        legend.position = "right",
        legend.justification = c(0,0.5))
fig.rho.sir <- ((A.msa[[yindex]]|B.msa[[yindex]])/C.msa[[yindex]]/E.msa[[yindex]])+ 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/rho_pop_mob_sir_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*3)
#----------plot infected rate for each msa with all (beta,mu)----------#
dset<-unique(sir.all$msa)
dset<-c("Atlanta","San Francisco")
for(d in 1:length(dset)){
  F.msa[[d]] <- ggplot(data = subset(sir.all,msa==dset[d]), aes(x = time, y = infected/pop)) +
    geom_line(aes(color = rho, group = rho), linewidth = 0.8) +
    facet_grid(mu ~ beta, scales = "fixed", labeller = labeller(
      beta = function(x) paste0(expression(beta), " = ", x),
      mu = function(x) paste0(expression(mu), " = ", x)
    )) + 
    scale_colour_gradientn(name = dset[d], colors = col.fit) +
    labs(x = "Time", y = "Infected fraction", title = NULL) +
    xlim(0, 500) +  # Limit time to 500
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain", size = 15),
          axis.text.y = element_text(hjust = 0),
          legend.position = "right",
          legend.justification = c(0, 0.5))
}
fig.rho.sir <- (F.msa[[1]]/F.msa[[2]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/SI_sir_beta_mu_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*4)





#----------plot SIR results----------#
#----------data from the spreading.R----------#
yindex<-2
sir.all <- read.csv(paste(datapath.msa,"/msa/SIR_",Yname[yindex],".csv",sep=""), header=TRUE)
E.msa <- F.msa <- list()
#----------plot infected rate for two msas with one (beta,mu)----------#
sir.all.main <- subset(sir.all,beta==0.3&mu==0.1&msa%in%c("Atlanta","San Francisco"))
E.msa[[yindex]] <- ggplot(data = sir.all.main, aes(x = time, y = infected/pop)) +
  geom_line(aes(color = rho, group = rho), linewidth = 0.8) +
  # geom_point(aes(color = rho, fill = rho), size=0.5, stroke = 0.1) +
  facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  scale_colour_gradientn(name = TeX('$\\rho$'), colors = col.fit) +
  scale_fill_gradientn(name = TeX('$\\rho$'), colors = col.fit) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        panel.spacing = unit(0.1, "lines"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
        strip.text = element_text(face = "plain",size=15),
        axis.text.y = element_text(hjust = 0),
        legend.position = "right",
        legend.justification = c(0,0.5))
fig.rho.sir <- ((A.msa[[yindex]]|B.msa[[yindex]])/C.msa[[yindex]]/E.msa[[yindex]])+ 
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/rho_pop_mob_sir_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*3)
#----------plot infected rate for each msa with all (beta,mu)----------#
dset<-unique(sir.all$msa)
dset<-c("Atlanta","San Francisco")
for(d in 1:length(dset)){
  F.msa[[d]] <- ggplot(data = subset(sir.all,msa==dset[d]), aes(x = time, y = infected/pop)) +
    geom_line(aes(color = rho, group = rho), linewidth = 0.8) +
    facet_grid(mu ~ beta, scales = "fixed", labeller = labeller(
      beta = function(x) paste0(expression(beta), " = ", x),
      mu = function(x) paste0(expression(mu), " = ", x)
    )) + 
    scale_colour_gradientn(name = dset[d], colors = col.fit) +
    labs(x = "Time", y = "Infected fraction", title = NULL) +
    xlim(0, 500) +  # Limit time to 500
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain", size = 15),
          axis.text.y = element_text(hjust = 0),
          legend.position = "right",
          legend.justification = c(0, 0.5))
}
fig.rho.sir <- (F.msa[[1]]/F.msa[[2]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/SI_sir_beta_mu_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*4)











#----------mobility map for each county----------#
col.mob<-colorRampPalette(rev(c('#22223b', '#4a4e69', '#9a8c98', '#c9ada7', '#f2e9e4')))(10)
# col.mob<-rev(c('#22223b', '#4a4e69', '#9a8c98', '#c9ada7', '#f2e9e4'))
E.msa <- F.msa<- G.msa <- list()
dset<-c(2,7,8,10)
for(yindex in 2:Ynum){
  
  #----------rho and msa sorted by median_rho----------#
  rho.mob.msa <- read.csv(paste(datapath.msa,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  rho.mob.msa <- rho.mob.msa[index == "exp"]
  rho.mob.msa$week.name<-format(as.Date(rho.mob.msa$day), "%a")
  
  rho.diff.msa <- rho.mob.msa[, .(before.rho = rho[period == "before"], 
                                  during.rho = rho[period == "during"]), 
                              by = .(msa, event, week.name)]
  rho.diff.msa[, `:=`(det_rho = during.rho - before.rho,
                      percent_det_rho = (during.rho - before.rho) / before.rho)]
  rho.select <- rho.diff.msa[, .(
    min_week = week.name[which.min(percent_det_rho)],
    max_week = week.name[which.max(percent_det_rho)]
  ), by = .(msa, event)]
  
  
  for(d in 1:4){
    
    Dname<<-Dname.set[dset[d]]
    # source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    s1<-subset(rho.select,msa==Dname)
    rho.min <- rho.mob.msa[msa == Dname & event=="COVID-19" & week.name==s1$min_week]
    rho.max <- rho.mob.msa[msa == Dname & event=="COVID-19" & week.name==s1$max_week]
    rho.min.max<-plyr::rbind.fill(rho.min,rho.max)
    day.min.max<-c(rho.min$day,rho.max$day)%>%as.Date
    
    #----------BlockID, NBlock----------#
    block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
    BlockID<-unique(block.msa$CensusBlockGroup)
    NBlock<-length(BlockID)
    
    #----------mob----------#
    mob <- NULL
    for(i in 1:Day){
      s1 <- read.csv(paste(flowpath,"/msa/", Dname, "/Intra_Flow_", Datevalue[i], ".csv", sep = ""), header=TRUE)
      s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                      data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                      by="CensusBlockGroup")
      s2$mob[is.na(s2$mob)] <- 0
      s2$day <- Datevalue[i]
      mob <- plyr::rbind.fill(mob,s2)
    } # day
    min.mob <- min(mob$mob, na.rm = TRUE)
    max.mob <- max(mob$mob, na.rm = TRUE)
    
    #----------plot mob map----------#
    for(i in 1:2){
      
      SSblock<-left_join(block.msa, subset(mob,day==day.min.max[2*i-1]), by="CensusBlockGroup")
      E.msa[[i]] <- ggplot() +
        geom_sf(data = SSblock,aes(fill= mob), size=0.01, lwd = 0.05,
                color="#7B8294", 
                show.legend = TRUE) +
        scale_fill_gradientn(name = paste(
          paste("MSA:", Dname),
          paste("Date:", format(day.min.max[2*i-1], "%m-%d")), 
          paste("Decay rate: ",round(rho.min.max$rho[2*i-1], 2)),     
          sep = '\n'                                 
        ),
        # limits = c(min.mob, max.mob),
        # values = scales::rescale(seq(min.mob, max.mob, length.out = length(col.map) + 1)), 
        values = scales::rescale(seq(0, 1, length.out = length(col.mob)+1)),
        colors = col.mob) +
        guides(fill = guide_legend(title.position = "top", title.hjust = 0))+
        theme_wy() +
        theme(legend.title = element_text(size=12),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())
      
      
      SSblock<-left_join(block.msa, subset(mob,day==day.min.max[2*i]), by="CensusBlockGroup")
      F.msa[[i]] <- ggplot() +
        geom_sf(data = SSblock,aes(fill= mob), size=0.01, lwd = 0.05,
                color="#7B8294", 
                show.legend = TRUE) +
        scale_fill_gradientn(name = paste(
          paste("MSA:", Dname),
          paste("Date:", format(day.min.max[2*i], "%m-%d")), 
          paste("Decay rate: ",round(rho.min.max$rho[2*i], 2)),  
          sep = '\n'                                 
        ),
        # limits = c(min.mob, max.mob),
        # values = scales::rescale(seq(min.mob, max.mob, length.out = length(col.map) + 1)), 
        values = scales::rescale(seq(0, 1, length.out = length(col.mob)+1)),
        colors = col.mob) +
        guides(fill = guide_legend(title.position = "top", title.hjust = 0))+
        theme_wy() +
        theme(legend.title = element_text(size=12),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank())
      
    } # day
    
    G.msa[[d]] <- (E.msa[[1]]|F.msa[1]|E.msa[[2]]|F.msa[2])
    
  } # msa
  
  fig.mob.map <- (G.msa[[1]]/G.msa[[2]]/G.msa[[3]]/G.msa[[4]])+ 
    plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.mob.map,filename = paste(figpath,"/mob_map_",Yname[yindex],".pdf",sep=""), width =5*4, height = 5*4)
} # yindex


