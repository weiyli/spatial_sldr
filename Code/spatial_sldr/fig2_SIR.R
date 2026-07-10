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
Flow.name <- "Intra_Flow"
Dname<<-Dname.set[d]
date.gif <- date.sep(Datevalue,durdate)


#----------confirmed cases for COVID-19----------#
con.us <- read.csv(paste(datapath.dis,"/us-counties-daily-cumulative-case.csv",sep=""), header=TRUE)
con.us$fips <- stringr::str_pad(con.us[,1], width = 5, pad = "0")
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
C.msa <- list()
for(yindex in 2:2){
 
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
    labs(x="MSA", y = TeX('Percentage change in $\\rho$')) +
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
  
} # yindex



#----------plot SIR results----------#
E.msa <- F.msa <- G.msa <- H.msa <- list()
yindex <- 2
#----------data from the spreading_model.R----------#
sir.model <- read.csv(paste(datapath.msa,"/msa/SIR_model_as_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.model$DM <- "Model"
sir.model[which(sir.model$det.rho==0),]$DM<-"Empirical"
sir.model.main <- subset(sir.model,beta==0.3 & mu==0.1 & msa%in%c("Atlanta","San Francisco"))

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
  scale_x_continuous(limits = c(-0.1, 0.1),
                     breaks = seq(-0.1, 0.1, 0.02),
                     labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(limits = c(-0.15, 0.12),
                     breaks = seq(-0.14, 0.12, 0.02),
                     labels = scales::percent_format(scale = 100)) +
  labs(x=TeX('Percentage change in $\\rho$'), y='Percentage change in infection peak')+
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
sir.model.time <- subset(sir.model.main, det.rho%in%c(-0.03,0,0.03))
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
                   max.overlaps = Inf,  nudge_x = 40,
                   box.padding = 1, size = 3, fill = 'white',color="black") +
  facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  scale_colour_manual(name = TeX('% change in $\\rho$'),
                      breaks = rho.info$labels,
                      labels = rho.info$labels,
                      values = rho.info$col) +
  scale_linetype_manual(name = TeX('% change in $\\rho$'), 
                        breaks = rho.info$labels,
                        labels = rho.info$labels,
                        values = rho.info$line) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, by = 30)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        strip.background = element_rect(fill = "gray90", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain", size=15, color = event.info$col[d]),
        legend.position = "top",
        legend.justification = c(0.5,0.5))

#----------plot infection peak and rho----------#
sir.model.peak <- as.data.table(sir.model.main)[, .SD[which.max(infected)], by = .(msa, day, rho, det.rho)]
baseline <- sir.model.peak[det.rho == 0, .(msa, day, time.base = time, infected.base = infected)]
sir.model.peak <- merge(sir.model.peak, baseline, by = c("msa", "day"), all.x = TRUE)
sir.model.peak[, ':='(
  time.change =  (time - time.base),
  infected.change = (infected - infected.base)
)]

sir.peak.change <- melt(
  sir.model.peak,
  id.vars = c("msa", "day", "det.rho"),
  measure.vars = c("time.change", "infected.change"),
  variable.name = "metric",
  value.name = "change"
)
# make metric labels prettier (optional)
sir.peak.change[metric == "time.change", metric := "Peak time change"][metric == "infected.change",  metric := "Peak infected change"]
sir.peak.change <- sir.peak.change[det.rho%in%c(-0.09,-0.06,-0.03,0.03,0.06,0.09)]
G.msa[[yindex]] <- ggplot(sir.peak.change, aes(x = factor(det.rho), y = change, fill = msa)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "white") +
  geom_abline(intercept = 0, slope = 0, linetype = "solid", linewidth = 0.5, color = "gray") +
  geom_text(aes(label = round(change), color = msa, vjust = ifelse(change < 0, 0.3, -0.3)),  
            hjust = 0.5, size = 4,
            position = position_dodge(width = 0.6),
            alpha = 0.8, 
            show.legend = FALSE) +
  scale_x_discrete(labels = function(x) paste0(round(as.numeric(as.character(x)) * 100), "%")) + 
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(name = 'MSA',
                    breaks = msa.info$dname,
                    labels = msa.info$dname,
                    values = scales::alpha(msa.info$col,alpha=0.8))+
  scale_colour_manual(name = 'MSA',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
  labs(x = TeX('Percentage change in $\\rho$'), y='Change in time or infected population', fill = "MSA") +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "inside",
        legend.position.inside = c(0.7, 0.8),
        legend.justification = c(0, 0))
H.msa[[yindex]] <- ggplot(sir.peak.change[metric=="Peak infected change"], aes(x = factor(det.rho), y = change, fill = msa)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "white") +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_text(aes(label = round(change), color = msa, 
                vjust = ifelse(change < 0,  1.5 * ifelse(msa == "Atlanta", 1.5, 1), -0.3 * ifelse(msa == "Atlanta", 1.5, 1))),  
            hjust = 0.5, size = 3, position = position_dodge(width = 0.6), alpha = 1, show.legend = FALSE) +
  scale_x_discrete(labels = function(x) paste0(round(as.numeric(as.character(x)) * 100), "%")) + 
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08)), breaks = scales::pretty_breaks(n = 10)) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(name = 'MSA',
                    breaks = msa.info$dname,
                    labels = msa.info$dname,
                    values = scales::alpha(msa.info$col, alpha = 0.7))+
  scale_colour_manual(name = 'MSA',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
  labs(x = TeX('Percentage change in $\\rho$'), y='Change in infected population', fill = "MSA") +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "top",
        legend.justification = c(0.5,0.5))
# fig.rho.sir <- (C.msa[[yindex]]/(E.msa[[yindex]]|F.msa[[yindex]])) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
fig.rho.sir <- (C.msa[[yindex]]/(F.msa[[yindex]]|H.msa[[yindex]])) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/brief_rho_pop_mob_sir_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 6.5*2)



