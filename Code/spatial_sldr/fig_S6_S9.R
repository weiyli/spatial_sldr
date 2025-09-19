# rm(list = ls())


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
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
  
  #----------BlockID, NBlock----------#
  block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.msa$CountyID <- str_pad(block.msa$CensusBlockGroup, width = 12, pad = "0")
  block.msa$CountyID <- substring(block.msa$CountyID,first=1,last=5)
  CountyID<-unique(block.msa$CountyID)
  con.county<-subset(con.us,fips%in%CountyID)
  con <-data.frame(day=con.date, case=colSums(con.county[,-1]))
  rownames(con)<-NULL
  con <- con %>% arrange(day) %>% mutate(new_case = case - lag(case, default = first(case)))
  con <- subset(con, day%in%Datevalue)
  con$msa<-Dname
  cases <- plyr::rbind.fill(cases, con)
}
write.csv(cases,file=paste(datapath.dis,"/us-msas-daily-cumulative-case.csv",sep=""),row.names = FALSE)


#----------percent change of median rho and confirmed cases for COVID-19----------#
con.msa <- read.csv(paste(datapath.dis,"/us-msas-daily-cumulative-case.csv",sep=""), header=TRUE)
queen.lag.msa <- read.csv(paste(datapath.msa,"/msa/queen_lag.csv",sep=""), header=TRUE)
pop.msa <- queen.lag.msa[,c("msa","pop")]%>%distinct
con.msa <- left_join(con.msa,pop.msa,by="msa")

#----------SLDR fitting for pop from sldr_fit_pop.R----------#
# rho_pop_msa
rho.pop.msa <- read.csv(paste(datapath.msa, "/msa/SLDR_params_pop.csv",sep=""), header=TRUE)%>%setDT
rho.pop.msa <- rho.pop.msa[index == "exp"]
rho.pop.msa[, c("error", "index") := NULL]

# rho_pop_dis
rho.pop.dis <- read.csv(paste(datapath.dis,"/msa/SLDR_params_pop.csv",sep=""), header=TRUE)%>%setDT
rho.pop.dis <- rho.pop.dis[index == "exp"]
rho.pop.dis[, c("error", "index") := NULL]
rho.pop.dis <- rho.pop.dis %>% distinct()  # Removes exact duplicates

# plot rho of pop and mob for different disasters
A.msa <- B.msa <- C.msa <- D.msa <- list()
for(yindex in 1:Ynum){
  
  # rho_mob_msa
  rho.mob.msa <- read.csv(paste(datapath.msa,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  rho.mob.msa <- rho.mob.msa[, .(
    min_rho = min(rho),
    median_rho = median(rho),
    max_rho = max(rho)
  ), by = .(msa, event, index, period)]
  # rho_pop_msa vs rho_mob_msa
  rho.pop.mob.msa <- left_join(rho.mob.msa[index == "exp"], rho.pop.msa, by=c('msa','event'))
  
  # rho_mob_dis
  rho.mob.dis <- read.csv(paste(datapath.dis,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  rho.mob.dis <- rho.mob.dis[period!="after"][, .(
    min_rho = min(rho),
    median_rho = median(rho),
    max_rho = max(rho)
  ), by = .(msa, event, index, period)]
  # rho_pop_dis vs rho_mob_dis
  rho.pop.mob.dis <- left_join(rho.mob.dis[index == "exp"], rho.pop.dis, by=c('msa','event'))
  
  rho.pop.mob <- rbind(rho.pop.mob.msa, rho.pop.mob.dis)
  xy.range <- range(c(rho.pop.mob$rho,rho.pop.mob$median_rho))
  
  rho.diff <- rho.pop.mob[, .(det_rho = median_rho[period == "during"] - median_rho[period == "before"],
                              percent_det_rho = (median_rho[period == "during"] - median_rho[period == "before"])/median_rho[period == "during"]), by = .(msa, event)]
  rho.rank <- rho.diff[order(percent_det_rho)]
  rho.rank[, msa_rank := .I]
  rho.rank<-left_join(rho.rank,data.frame(event=event.info$breaks,event_name=event.info$labels)%>%distinct,by="event")
  
  #----------rho_mob and rho_pop----------#
  A.msa[[yindex]] <- ggplot(rho.pop.mob, aes(x = rho, y = median_rho)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
    geom_label_repel(data = subset(rho.pop.mob, period == "before" & msa%in%c("Atlanta","New York")),
                     aes(label = msa), max.overlaps = Inf,
                     box.padding = 1, size = 3, fill = 'white',color=event.info$col[1]) +
    geom_line(aes(color=event, group = interaction(msa, event)), linewidth = 0.5, linetype = "dashed") +
    geom_point(aes(fill = event, color = event, shape = period), size = 3, stroke = 1) +
    scale_fill_manual(name =  'Disaster',
                      breaks = event.info$breaks,
                      labels = event.info$labels,
                      values = event.info$col)+
    scale_color_manual(name =  'Disaster',
                       breaks = event.info$breaks,
                       labels = event.info$labels,
                       values = event.info$col) +
    scale_shape_manual(name = "Period",
                       breaks = period.info$breaks,
                       labels = period.info$labels.new,
                       values = c(21,25)) +
    scale_x_continuous(limits = xy.range) +
    scale_y_continuous(limits = xy.range) +
    labs(x = TeX('Spatial range exponent $\\rho_{pop}$'), y = TeX('Spatial range exponent $\\rho_{mob}$')) +
    guides(fill = guide_legend(order = 1, override.aes = list(size=3)),  
           color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          legend.position ="right",
          legend.text = element_text(hjust=0),
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank())
  
  
  
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
  pearson.cor<-cor(rho.case.msa$percent_det_rho,log(rho.case.msa$case/rho.case.msa$pop,10),method="pearson")
  
  #----------plot the percent change of rho and cum confirmed cases----------#
  xrange<-range(rho.case.msa$percent_case)
  B.msa[[yindex]] <- ggplot(rho.case.msa, aes(x = percent_case, y = percent_det_rho)) +
    geom_abline(intercept = 0, slope = -1, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
    geom_point(aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
    scale_fill_manual(name="COVID-19",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = scales::alpha(msa.info$col,alpha=0.8))+
    scale_colour_manual(name = "COVID-19",
                        breaks = msa.info$dname,
                        labels = msa.info$dname,
                        values = msa.info$col)+
    scale_shape_manual(name="COVID-19",
                       breaks = msa.info$dname,
                       labels = msa.info$dname,
                       values = msa.info$shape) +
    labs(x = TeX('Cumulative confirmed cases/population ($\\log$)'),
         y = TeX('The percentage change of $\\rho_{mob}$')) +
    # scale_x_log10() +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          # legend.title = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "right",
          legend.justification = c(0,0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank()) +
    annotation_custom(grob = textGrob(TeX(sprintf("$Pearson's\\,r = %.2f$", round(pearson.cor, 2))), gp = gpar(fontsize = 12), hjust = 0), xmin = log(xrange[2]/7,10), xmax = log(xrange[2]/7,10), ymin = 0.018, ymax = 0.02)
  
  
  #----------rho and msa sorted by median_rho----------#
  rho.mob.msa <- read.csv(paste(datapath.msa,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  s1 <- rho.mob.msa[index == "exp"]
  s1$week.name<-format(as.Date(s1$day), "%a")
  
  rho.diff.msa <- s1[, .(before.rho = rho[period == "before"], 
                         during.rho = rho[period == "during"]), 
                     by = .(msa, event, week.name)]
  rho.diff.msa[, `:=`(det_rho = during.rho - before.rho,
                      percent_det_rho = (during.rho - before.rho) / before.rho)]
  rho.diff.msa<-rho.diff.msa[, .(msa, event, det_rho, percent_det_rho)]
  # rho.diff.msa <- s1[, .(det_rho = rho[period == "during"] - rho[period == "before"],
  #                        percent_det_rho = (rho[period == "during"] - rho[period == "before"])/rho[period == "during"]), by = .(msa, event)]
  
  # disaster
  rho.mob.dis <- read.csv(paste(datapath.dis,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
  rho.mob.dis.exp <- rho.mob.dis[index == "exp"&period!="after"]
  rho.mob.dis.exp$week.name<-format(as.Date(rho.mob.dis.exp$day), "%a")
  
  rho.mob.dis.exp <- rho.mob.dis.exp[,.(rho=mean(rho)),by=.(msa,event,period,week.name)]
  rho.mob.dis.exp <- rho.mob.dis.exp[, 
                                     .SD[week.name %in% intersect(week.name[period == "before"], week.name[period == "during"])], by = .(msa, event)]
  rho.diff.dis <- rho.mob.dis.exp[, .(before.rho = rho[period == "before"], 
                                      during.rho = rho[period == "during"]), 
                                  by = .(msa, event, week.name)]
  rho.diff.dis[, `:=`(det_rho = during.rho - before.rho,
                      percent_det_rho = (during.rho - before.rho) / before.rho)]
  rho.diff.dis<-rho.diff.dis[, .(msa, event, det_rho, percent_det_rho)]
  # s1 <- rho.mob.dis.exp[,.N,by=.(msa,event,period)]
  # s2 <- s1[, .SD[which.max(N)], by = .(msa, event)]
  # s3 <- s1[, .(N_diff = abs(diff(N))), by = .(msa, event)]
  # s4 <- merge(s2, s3, by = c("msa", "event"))
  # for (i in 1:nrow(s4)) {
  #   msa_curr <- s4[i, msa]
  #   event_curr <- s4[i, event]
  #   period_curr <- s4[i, period]
  #   N_diff_curr <- s4[i, N_diff]
  #   df_subset <- rho.mob.dis.exp[msa == msa_curr & event == event_curr & period == period_curr]
  #   if(N_diff_curr==0){
  #     rho.mob.dis.exp <- rho.mob.dis.exp
  #   }else{
  # 
  #     if (period_curr == "before") {
  #       days_to_remove <- df_subset[order(day)][1:N_diff_curr, day]
  #     } else {
  #       days_to_remove <- df_subset[order(day)][(.N - N_diff_curr + 1):.N, day]
  #     }
  #     rho.mob.dis.exp <- rho.mob.dis.exp[!(msa == msa_curr & event == event_curr & period == period_curr & day %in% days_to_remove)]
  #   }
  # }
  # rho.diff.dis <- rho.mob.dis.exp[, .(det_rho = rho[period == "during"] - rho[period == "before"],
  #                             percent_det_rho = (rho[period == "during"] - rho[period == "before"])/rho[period == "during"]), by = .(msa, event)]
  rho.diff <- plyr::rbind.fill(rho.diff.msa, rho.diff.dis)%>%setDT
  rho.diff.bar <- rho.diff[, .(
    min_rho = min(percent_det_rho),
    q25_rho = quantile(percent_det_rho, probs = 0.25),
    median_rho = median(percent_det_rho),
    q75_rho = quantile(percent_det_rho, probs = 0.75),
    max_rho = max(percent_det_rho)
  ), by = .(msa,event)]
  rho.rank <- rho.diff.bar[order(median_rho)]
  rho.rank[, msa_rank := .I]
  rho.rank<-left_join(rho.rank,data.frame(event=event.info$breaks,event_name=event.info$labels)%>%distinct,by="event")
  
  C.msa[[yindex]] <- ggplot(rho.rank, aes(x = msa_rank, y = median_rho, color = event)) +
    # geom_errorbar(aes(ymin = min_rho, ymax = max_rho), width = 0.2, linewidth = 0.5) +  
    geom_errorbar(aes(ymin = q25_rho, ymax = q75_rho), width = 0.2, linewidth = 0.5) +  
    geom_point(aes(fill = event, color = event), size = 3, stroke = 0.5) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
    geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Harvey","Storm_Texas")|(event=="COVID-19"&msa=="Houston")), aes(label = event_name, color = event), 
                     nudge_x = -0.5, nudge_y = 0.008, box.padding = 0.3, 
                     point.padding = 1, size = 3, fill = 'white') +
    geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Dorian","Fire_Kincade")), 
                     aes(label = event_name, color=event), nudge_x = 0.5, nudge_y = -0.008, 
                     box.padding = 0.3, point.padding = 1, size = 3, fill = 'white') +
    scale_fill_manual(name =  'Disaster',
                      breaks = event.info$breaks,
                      labels = event.info$labels,
                      values = scales::alpha(event.info$col,alpha=0.5))+
    scale_color_manual(name =  'Disaster',
                       breaks = event.info$breaks,
                       labels = event.info$labels,
                       values = event.info$col) +
    scale_x_continuous(breaks = rho.rank$msa_rank, labels = rho.rank$msa) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
    labs(x="MSA", y = TeX('Percentage change in $\\rho_{mob}$')) +
    theme_wy() +
    theme(panel.border = element_blank(),
          legend.key.size = unit(0.5, "cm"),
          legend.position = "none",
          legend.justification = c(0, 0.5),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.line.x = axis.arrow,
          axis.line.y = axis.arrow)
  
  #----------percent change of median rho----------#
  D.msa[[yindex]] <- ggplot(rho.rank, aes(x = msa_rank, y = median_rho)) +
    geom_point(aes(fill = event, color = event), shape=21, size = 3, stroke = 2) +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
    geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Harvey","Storm_Texas")|(event=="COVID-19"&msa=="Houston")), aes(label = event_name, color=event), nudge_x = -0.5, nudge_y = 0.004, 
                     box.padding = 0.3, point.padding = 1, size = 3, fill = 'white') +
    geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Dorian","Fire_Kincade")), 
                     aes(label = event_name, color=event), nudge_x = 0.5, nudge_y = -0.004, 
                     box.padding = 0.3, point.padding = 1, size = 3, fill = 'white') +
    scale_fill_manual(name =  'Disaster',
                      breaks = event.info$breaks,
                      labels = event.info$labels,
                      values = scales::alpha(event.info$col,alpha=0.5))+
    scale_color_manual(name =  'Disaster',
                       breaks = event.info$breaks,
                       labels = event.info$labels,
                       values = event.info$col) +
    scale_x_continuous(breaks = rho.rank$msa_rank, labels = rho.rank$msa) +  
    labs(x="MSA",y = TeX('Change in spatial range exponent $\\rho_{mob}$')) +
    theme_wy() +
    theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          # axis.title.x = element_blank(),
          legend.position ="none")
  
  
  fig.rho <- ((A.msa[[yindex]]|B.msa[[yindex]])/C.msa[[yindex]]/D.msa[[yindex]])+ 
    plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig.rho,filename = paste(figpath,"/msa/rho_pop_mob_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5*3)
  
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
# theme(plot.title = element_text(size = 15),
#       axis.title = element_text(size = 12),
#       axis.text.x = element_text(size = 12),
#       axis.text.y = element_text(hjust = 0,size = 12),
#       strip.text = element_text(size = 12),        
#       legend.title = element_text(size = 12),
#       legend.text = element_text(size = 12),
#       panel.background = element_blank(),
#       panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
#       panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
#       panel.spacing = unit(0.1, "lines"),
#       axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
#       legend.position = "right",
#       legend.justification = c(0,0.5))
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


