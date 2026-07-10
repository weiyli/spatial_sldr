# rm(list = ls())

# depend on R code from fig_SI_rho_comp_covid.R and fig_SI_covid_case_alignment.R


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
library(ggstar)


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
yindex <- 2   
date.gif <- date.sep(Datevalue,durdate)
period_norm  <- pervalue[1]
period_covid <- pervalue[2]


#----------Load data: rho----------#
rho_covid <- fread(paste0(datapath.msa, "/msa/SLDR_params_", Yname[yindex], ".csv"))
rho_dis <- fread(paste0(datapath.dis, "/msa/SLDR_params_", Yname[yindex], ".csv"))


#-------------------------
# Part 1: |螖rho_COVID|/inter-city variation in rho, |螖rho_COVID|/normal daily variation in rho
#-------------------------
#-------------------------
# 1) City-level mean rho by period and one-week COVID-induced changes (delta rho)
#-------------------------
rho_all <- rbind(rho_covid, rho_dis)
rho_all <- as.data.table(rho_all)[, week := NULL]
rho_all <- rho_all[index == "exp" & period != "after"]
rho_all[, day := as.Date(day)]
Sys.setlocale("LC_TIME", "C") 
rho_all[, wday := weekdays(day)]
rho_all[, .N, by = .(msa, event, wday, period)][N > 1]
rho_all <- rho_all[,.(rho = mean(rho)), by = .(msa, event, wday, period)]

#-------------------------
# Analysis 1: |螖rho_COVID| vs inter-city variation in rho
#-------------------------
# |螖rho| per city |COVID - normal|
rho_wide <- dcast(rho_all, msa + event + wday ~ period, value.var = "rho")%>%na.omit()
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
# median(|螖rho|)
rho_delta_city <- rho_wide[,  .(abs_delta_rho_covid = median(abs_delta_rho)), by = .(msa, event)]
# Inter-city SD of median rho in the normal period 
rho_city_norm <- rho_wide[, .(rho_city = median(before, na.rm = TRUE)), by = .(msa, event)]
sd_across_cities <- sd(rho_city_norm$rho_city, na.rm = TRUE)
iqr_across_cities <- IQR(rho_city_norm$rho_city, na.rm = TRUE)
# Ratio per city (using SD across cities in normal period, as in your note)
rho_delta_city[, ratio_sd := abs_delta_rho_covid / sd_across_cities]
rho_delta_city[, ratio_iqr := abs_delta_rho_covid / iqr_across_cities]
mean(rho_delta_city$ratio_sd)
# Plot: a compact bar/point plot for ratio (one number per city)
fig_city <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = event, color = event), size = 4, stroke = 1) +
  scale_fill_manual(name =  'Disaster',
                    breaks = event.info$breaks,
                    labels = event.info$labels,
                    values = scales::alpha(event.info$col,alpha=0.5))+
  scale_color_manual(name =  'Disaster',
                     breaks = event.info$breaks,
                     labels = event.info$labels,
                     values = event.info$col) +
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[disaster] * "|", SD(rho[inter-city])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        # panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        # panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

#-------------------------
# Analysis 2: |螖rho_COVID| vs null distribution of |螖rho| in normal period
#   螖rho_null(t)=rho_{t+1}-rho_t  (within each city, normal period)
#   螖rho_COVID(t)=rho_{t+1}-rho_t (within each city, COVID period)
#-------------------------
# |螖rho| per city |COVID - normal|
rho_wide <- dcast(rho_all, msa + event + wday ~ period, value.var = "rho")%>%na.omit()
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
# median(|螖rho|)
rho_delta_city <- rho_wide[,  .(abs_delta_rho_covid = median(abs_delta_rho)), by = .(msa, event)]
# Within-city SD of rho during normal period 
rho_city_norm_sd <- rho_all[period == period_norm,
                            .(sd_rho_normal = sd(rho, na.rm = TRUE)),
                            by = .(msa, event)]
# Merge and compute E = |螖rho_COVID| / SD(rho_normal)
rho_delta_city <- merge(rho_delta_city, rho_city_norm_sd, by = c("msa", "event"), all.x = TRUE)
rho_delta_city[, ratio_sd := abs_delta_rho_covid / sd_rho_normal]
# plot: a compact bar/point plot for ratio (one number per city)
fig_normal <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = event, color = event),
             size = 4, stroke = 1) +
  scale_fill_manual(name =  'Disaster',
                    breaks = event.info$breaks,
                    labels = event.info$labels,
                    values = scales::alpha(event.info$col,alpha=0.5)) +
  scale_color_manual(name =  'Disaster',
                     breaks = event.info$breaks,
                     labels = event.info$labels,
                     values = event.info$col) +
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[disaster] * "|", SD(rho[normal])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        # panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        # panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")
# fig_comp <- (fig_city|fig_normal) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
# ggsave(fig_comp, filename = paste(figpath,"/msa/SI_rho_comp_",Yname[yindex],".pdf",sep=""), width = 7*2, height = 7)



#-------------------------
# Part 2: change in rho and confirmed cases for COVID-19
#-------------------------
# #----------data: confirmed cases for COVID-19----------#
# con.us <- fread(paste0(datapath.dis, "/us-counties-daily-cumulative-case.csv"))
# con.us$fips <- str_pad(con.us[,1], width = 5, pad = "0")
# con.date <- as.Date(sub("X", "", colnames(con.us[,-1])), format = "%Y%m%d")
# cases <- NULL
# for(d in 1:Nmsa){
#   
#   Dname<<-Dname.set[d]
#   source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
#   date.gif<-date.sep(Datevalue,durdate)
#   
#   Dname<<-Dname.set[d]
#   source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
#   date.gif <- date.sep(Datevalue,durdate)
#   
#   #----------BlockID, NBlock----------#
#   block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
#   block.msa$CountyID <- str_pad(block.msa$CensusBlockGroup, width = 12, pad = "0")
#   block.msa$CountyID <- substring(block.msa$CountyID,first=1,last=5)
#   CountyID <- unique(block.msa$CountyID)
#   con.county <- subset(con.us,fips%in%CountyID)
#   con <- data.frame(day=con.date, case=colSums(con.county[,-1]))
#   rownames(con) <- NULL
#   con <- con %>% arrange(day) %>% mutate(new_case = case - lag(case, default = first(case)))
#   con <- subset(con, day%in%Datevalue)
#   con$msa <- Dname
#   cases <- plyr::rbind.fill(cases, con)
# }
# write.csv(cases,file=paste(datapath.dis,"/us-msas-daily-cumulative-case.csv",sep=""),row.names = FALSE)

#----------percent change of median rho and confirmed cases for COVID-19----------#
con.msa <- fread(paste0(datapath.dis,"/us-msas-daily-cumulative-case.csv"))
queen.lag.msa <- fread(paste0(datapath.msa,"/msa/queen_lag.csv"))
pop.msa <- queen.lag.msa[,c("msa","pop")]%>%distinct
con.msa <- left_join(con.msa,pop.msa,by="msa")
#----------daily percent change of all rho----------#
# daily rho
rho.mob.msa <- rho_covid
rho.mob.msa[, day := as.Date(day)]
id.week<-data.frame(day=Datevalue,week.id=factor(rep(1:(Day/2),2)),week.name=format(Datevalue, "%a"))
daily.rho<-left_join(rho.mob.msa[index == "exp"],id.week,by="day")
daily.rho.diff <- daily.rho[, .(det_rho = rho[period == "during"]-rho[period == "before"],
                                percent_det_rho = (rho[period == "during"] - rho[period == "before"])/rho[period == "during"]), by = .(msa, week.id, week.name)]
# case and rho
con.msa$day<-as.Date(con.msa$day)
case.msa<-left_join(con.msa,date.gif,by="day")
case.msa<-subset(case.msa,period=="during")
case.msa <- left_join(case.msa,id.week,by="day")
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


#----------plot the percent change of rho and cum confirmed cases, then add binned-mean points----------#
fig_rho_case <- ggplot(rho.case.msa, aes(x = percent_case, y = percent_det_rho)) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
  geom_point(data = rho.case.bin_sum, aes(x = x_mean, y = y_mean, fill = "Bin-averaged", color = "Bin-averaged", shape = "Bin-averaged"), inherit.aes = FALSE, size = 5, stroke = 1, alpha = 1) +
  geom_text(data = ann_cor, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, 
            size = 4, color = "#009E73", parse = TRUE) + 
  scale_fill_manual(
    name = "COVID-19",
    breaks = c(msa.info$dname, "Bin-averaged"),
    labels = c(msa.info$dname, "Bin-averaged"),
    values = c(setNames(scales::alpha(msa.info$col, alpha = 0.5), msa.info$dname), "Bin-averaged" = "#69b3a2")
  ) +
  scale_color_manual(
    name = "COVID-19",
    breaks = c(msa.info$dname, "Bin-averaged"),
    labels = c(msa.info$dname, "Bin-averaged"),
    values = c(setNames(scales::alpha(msa.info$col, alpha = 0.5), msa.info$dname), "Bin-averaged" = "#69b3a2")
  ) +
  scale_shape_manual(
    name = "COVID-19",
    breaks = c(msa.info$dname, "Bin-averaged"),
    labels = c(msa.info$dname, "Bin-averaged"),
    values = c(setNames(msa.info$shape, msa.info$dname), "Bin-averaged" = 24)
  ) +
  labs(x = TeX('Cumulative confirmed cases/population ($\\log$)'),
       y = TeX('Percentage change in $\\rho$')) +
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
ggsave(fig_rho_case, filename = paste(figpath,"/msa/SI_rho_case_",Yname[yindex],".pdf",sep=""), width = 7.2, height = 5)
  

  

#-------------------------
# Part 3: percentage change in rho for all disasters
#-------------------------
#----------rho and msa sorted by median_rho----------#
s1 <- rho_covid[index == "exp"]
s1$week.name <- format(as.Date(s1$day), "%a")
rho.diff.msa <- s1[, .(before.rho = rho[period == "before"], 
                       during.rho = rho[period == "during"]), 
                   by = .(msa, event, week.name)]
rho.diff.msa[, `:=`(det_rho = during.rho - before.rho,
                    percent_det_rho = (during.rho - before.rho) / before.rho)]
rho.diff.msa<-rho.diff.msa[, .(msa, event, det_rho, percent_det_rho)]
# disaster
rho.mob.dis.exp <- rho_dis[index == "exp"&period!="after"]
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
rho.rank <- left_join(rho.rank,data.frame(event=event.info$breaks,event_name=event.info$labels)%>%distinct,by="event")
fig_rho_change <- ggplot(rho.rank, aes(x = msa_rank, y = median_rho, color = event)) +
  # geom_errorbar(aes(ymin = min_rho, ymax = max_rho), width = 0.2, linewidth = 0.5) +  
  geom_errorbar(aes(ymin = q25_rho, ymax = q75_rho), width = 0.2, linewidth = 0.5) +  
  geom_point(aes(fill = event, color = event), shape = 21, size = 4, stroke = 1) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Harvey","Storm_Texas")|(event=="COVID-19"&msa=="Houston")), aes(label = event_name, color = event), 
                   nudge_x = -0.5, nudge_y = 0.008, box.padding = 0.3, 
                   point.padding = 1, size = 4, fill = 'white') +
  geom_label_repel(data = subset(rho.rank, event %in% c("Hurricane_Dorian","Fire_Kincade")), 
                   aes(label = event_name, color=event), nudge_x = 0.5, nudge_y = -0.008, 
                   box.padding = 0.3, point.padding = 1, size = 3, fill = 'white') +
  scale_fill_manual(name =  'Disaster',
                    breaks = event.info$breaks,
                    labels = event.info$labels,
                    values = scales::alpha(event.info$col,alpha=1))+
  scale_color_manual(name =  'Disaster',
                     breaks = event.info$breaks,
                     labels = event.info$labels,
                     values = event.info$col) +
  scale_x_continuous(breaks = rho.rank$msa_rank, labels = rho.rank$msa) +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = c(0.0001,0.01)) + 
  labs(x="MSA", y = TeX('Percentage change in $\\rho$')) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, stroke = 0, colour = NA)), color = "none") +
  theme_wy() +
  theme(panel.border = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "right",
        legend.justification = c(0, 0.5),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line.x = axis.arrow,
        axis.line.y = axis.arrow)

fig2 <- (fig_rho_change/((fig_city|fig_normal|fig_rho_case) + plot_layout(widths = c(1, 1, 1.5)))) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig2, filename = paste(figpath,"/msa/fig2.pdf",sep=""), width = 18, height = 6*2)







