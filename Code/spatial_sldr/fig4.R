# rm(list = ls())

# depend on figure from fig_SI_covid_case_alignment.R, fig2_SIR.R
# depend on data from spreading_model.R

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath.msa <- 'D:/ood/Data/spatial_sldr'
datapath.dis <- 'D:/ood/Data/spatial_sldr/disaster'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(ggh4x)   # facet_wrap2()
library(ggrepel)

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

#----------Load data: rho----------#
rho_covid <- fread(paste0(datapath.msa, "/msa/SLDR_params_", Yname[yindex], ".csv"))

#-------------------------
# Part 1: change in rho and confirmed cases for COVID-19
#-------------------------
#----------data: confirmed cases for COVID-19----------#
con.us <- fread(paste0(datapath.dis, "/us-counties-daily-cumulative-case.csv"))
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
con.msa <- fread(paste0(datapath.dis,"/us-msas-daily-cumulative-case.csv"))
queen.lag.msa <- fread(paste0(datapath.msa,"/msa/queen_lag.csv"))
pop.msa <- queen.lag.msa[,c("msa","pop")]%>%distinct
con.msa <- left_join(con.msa,pop.msa,by="msa")
#----------daily percent change of all rho----------#
# daily rho
rho.mob.msa <- rho_covid
rho.mob.msa[, day := as.Date(day)]
id.week<-data.frame(day=Datevalue,week.id=factor(rep(1:(Day/2),2)),week.name=format(Datevalue, "%a"))
daily.rho <-left_join(rho.mob.msa[index == "exp"],id.week,by="day")
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
x_pos <- max(dt_raw$percent_case, na.rm = TRUE)*0.13
y_pos <- max(dt_raw$percent_det_rho, na.rm = TRUE)

# vertical spacing
y_range <- range(dt_raw$percent_det_rho, na.rm = TRUE)
dy <- 0.06 * diff(y_range)

ann_cor <- data.table::rbindlist(list(
  data.table(line = 1, x = x_pos, y = y_pos,
             label = sprintf("Raw~~Pearson~~italic(r)==%.2f",  pearson_raw)),
  data.table(line = 2, x = x_pos, y = y_pos - dy,
             label = sprintf("Raw~~Spearman~~italic(r)[s]==%.2f", spearman_raw)),
  data.table(line = 3, x = x_pos, y = y_pos - 2*dy,
             label = sprintf("Binned~~Pearson~~italic(r)==%.2f", pearson_bin)),
  data.table(line = 4, x = x_pos, y = y_pos - 3*dy,
             label = sprintf("Binned~~Spearman~~italic(r)[s]==%.2f", spearman_bin))
), use.names = TRUE)


#----------plot the percent change of rho and cum confirmed cases, then add binned-mean points----------#
rho.case.msa$event <- "COVID-19"
rho.case.bin_sum$event <- "COVID-19"
fig_rho_case <- ggplot(rho.case.msa, aes(x = percent_case, y = percent_det_rho)) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size= 3, stroke = 0.5) +
  geom_point(data = rho.case.bin_sum, aes(x = x_mean, y = y_mean, fill = "Bin-averaged", color = "Bin-averaged", shape = "Bin-averaged"), inherit.aes = FALSE, size = 5, stroke = 1, alpha = 1) +
  geom_text(data = ann_cor, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1, 
            size = 4, color = "#009E73", parse = TRUE) + 
  facet_wrap(~ event, ncol = 1, scales = "fixed") + 
  scale_fill_manual(
    name = "MSA",
    breaks = c(msa.info$dname, "Bin-averaged"),
    labels = c(msa.info$dname, "Bin-averaged"),
    values = c(setNames(scales::alpha(msa.info$col, alpha = 0.5), msa.info$dname), "Bin-averaged" = "#69b3a2")
  ) +
  scale_color_manual(
    name = "MSA",
    breaks = c(msa.info$dname, "Bin-averaged"),
    labels = c(msa.info$dname, "Bin-averaged"),
    values = c(setNames(scales::alpha(msa.info$col, alpha = 0.5), msa.info$dname), "Bin-averaged" = "#69b3a2")
  ) +
  scale_shape_manual(
    name = "MSA",
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
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")
ggsave(fig_rho_case, filename = paste(figpath,"/msa/SI_rho_case_",Yname[yindex],".pdf",sep=""), width = 7.2, height = 5)


#-------------------------
# Part 2: SIR results
#-------------------------
#----------data from the spreading_model.R----------#
select_msa <- c("Atlanta","San Francisco")
sir.model <- read.csv(paste(datapath.msa,"/msa/SIR_model_as_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.model$DM <- "Model"
sir.model[which(sir.model$det.rho==0),]$DM<-"Empirical"
sir.model.main <- subset(sir.model,beta==0.3 & mu==0.1 & msa%in%select_msa)
#----------plot infected pop and time----------#
sir.model.time <- subset(sir.model.main, det.rho%in%c(-0.03,0,0.03))
det.rho <- sort(unique(sir.model.time$det.rho),decreasing=FALSE)
rho.breaks <- rho.labels <- paste0(round(det.rho,4)*100,"%")
rho.labels[which(rho.breaks=="0%")] <- "empirical (0%)"
rho.info <- data.frame(det.rho = det.rho, 
                       breaks = rho.breaks,
                       labels = rho.labels,
                       color = col.fit[c(5,3,1)],
                       line = c("dotted","solid","dotdash"))
# color = col.all[1:length(det.rho)])
sir.model.time <- left_join(sir.model.time, rho.info, by="det.rho")
emp.rho <- as.data.table(sir.model.time)[DM == "Empirical", .SD[which.max(infected)], by = msa]
#----------plot infected rate for two msas with (beta=0.3,mu=0.1)----------#
strip_cols <- setNames(msa.info$col, msa.info$dname)[select_msa]
fig_sir_infected <- ggplot(data = sir.model.time, aes(x = time, y = infected/pop)) +
  geom_line(aes(color = labels, linetype = labels, group = rho), linewidth = 0.8) +
  geom_label_repel(data = emp.rho,
                   aes(label = paste("empirical",round(rho, 4),sep=" ")),
                   max.overlaps = Inf,  nudge_x = 40,
                   box.padding = 1, size = 3, fill = 'white',color = "#b56576") +
  # facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  # facet_wrap2(~ msa, ncol = 2, scales = "fixed", 
  #             strip = strip_themed(background_x = elem_list_rect(fill = strip_cols, color = "gray",linewidth = 0.5))) +
  facet_wrap2(~ msa, ncol = 2, strip = strip_themed(text_x = elem_list_text( color = strip_cols, size = 15, face = "plain"))) +
  scale_color_manual(name = TeX('Percentage change in $\\rho$'),
                     breaks = rho.info$labels,
                     labels = rho.info$labels,
                     values = rho.info$col) +
  scale_linetype_manual(name = TeX('Percentage change in $\\rho$'), 
                        breaks = rho.info$labels,
                        labels = rho.info$labels,
                        values = rho.info$line) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  scale_x_continuous(limits = c(0, 180), breaks = seq(0, 180, by = 30)) +
  guides(colour = guide_legend(keywidth = 2), linetype = guide_legend(keywidth = 2)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 15),
        # legend.justification = c(0.5,0.5),
        legend.position = "top")


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
sir.peak.change[metric == "time.change", metric := "Peak time change"][metric == "infected.change",  metric := "Change in peak number of infections"]
sir.peak.change <- sir.peak.change[det.rho%in%c(-0.09,-0.06,-0.03,0.03,0.06,0.09)]
fig_sir_pop <- ggplot(sir.peak.change[metric=="Change in peak number of infections"], aes(x = factor(det.rho), y = change, fill = msa)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "white") +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "gray") +
  geom_text(aes(# label = scales::comma(round(change)), 
                label = scales::comma(round(change, digits = -3)),
                color = msa, 
                hjust = ifelse(msa == "Atlanta", 0.5, 0.5),
                vjust = ifelse(change < 0,  1.5 * ifelse(msa == "Atlanta", 1.5, 1), -0.3 * ifelse(msa == "Atlanta", 1.5, 1))),  
            size = 4, position = position_dodge(width = 0.6), alpha = 1, show.legend = FALSE) +
  scale_x_discrete(labels = function(x) paste0(round(as.numeric(as.character(x)) * 100), "%")) + 
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08)), breaks = scales::pretty_breaks(n = 10), labels = scales::comma) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  scale_fill_manual(name = 'MSA',
                    breaks = msa.info$dname,
                    labels = msa.info$dname,
                    values = scales::alpha(msa.info$col, alpha = 0.7))+
  scale_colour_manual(name = 'MSA',
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
  labs(x = TeX('Percentage change in $\\rho$'), y='Change in the number of infections', fill = "MSA") +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = c(0.162, 0.91),  
        legend.justification = c(1, 1))

fig4 <- (((fig_rho_case|fig_sir_infected) + plot_layout(widths = c(1, 1.35)))/fig_sir_pop) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig4,filename = paste(figpath,"/msa/fig4.pdf",sep=""), width = 14, height = 5*2)






