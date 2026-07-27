# rm(list = ls())

# depend on data from disturb_experiment.R


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
  # labs(x = NULL, y = expression(frac("|" * Delta * rho[disaster] * "|", SD(rho[inter-city])))) +
  labs(x = NULL, y = expression(italic(E)[inter]~"(cross-city variation)")) + 
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
  # labs(x = NULL, y = expression(frac("|" * Delta * rho[disaster] * "|", SD(rho[normal])))) +
  labs(x = NULL, y = expression(italic(E)[normal]~"(within-city variation)")) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        # panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        # panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")
# fig_comp <- (fig_city|fig_normal) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
# ggsave(fig_comp, filename = paste(figpath,"/msa/SI_rho_comp_",Yname[yindex],".pdf",sep=""), width = 7*2, height = 7)




#-------------------------
# Part 2: 螖rho vs disturb (by scenario) with 95% CI
#-------------------------
dset <- which(Dname.set%in%c("Atlanta","San Francisco"))
#----------Run scenarios with multiple perturbation magnitudes----------#
disturb_set <- c(0.01, seq(0.10,0.90,0.10))
disturb_set <- c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
scenario_set <- c(
  "I_scale",
  "II_rewire_strat",
  "III_perm_null"
)
n_rep <- 100   # experiment repeat times

#-----------------------------
# Helpers: read baseline and disturb files
#-----------------------------
read_baseline_rho <- function(city_name, region_name, y_tag, datapath) {
  f0 <- file.path(datapath, region_name, "params",
                  paste0(city_name, "_SLDR_params_", y_tag, ".csv"))
  if (!file.exists(f0)) stop("Baseline file not found: ", f0)
  
  b0 <- fread(f0)
  # keep exponential-kernel rho as baseline
  b0 <- b0[index == "exp", .(day, rho_base = rho)]
  b0
}

read_disturb_rho <- function(city_name, region_name, y_tag, datapath) {
  f1 <- file.path(datapath, region_name, "params",
                  paste0(city_name, "_disturb_", y_tag, ".csv"))
  if (!file.exists(f1)) stop("Disturb file not found: ", f1)
  
  d1 <- fread(f1)
  d1[, msa := city_name]
  d1
}

#-----------------------------
# 1) Load disturb outputs and baseline (empirical) rho
#-----------------------------
region_name <- region[2]     # adjust if you wrote to another region folder
y_tag <- Yname[yindex]

dist_list <- list()
base_list <- list()

for (d in dset) {
  city_name <- Dname.set[d]
  
  dist_list[[city_name]] <- read_disturb_rho(city_name, region_name, y_tag, datapath.msa)
  
  base_list[[city_name]] <- read_baseline_rho(city_name, region_name, y_tag, datapath.msa) %>%
    mutate(msa = city_name)
}
dt_base <- rbindlist(base_list, fill = TRUE)
dt_dist <- rbindlist(dist_list, fill = TRUE)
dt_dist <- dt_dist[scenario %in% scenario_set]
dt_dist <- dt_dist[disturb %in% disturb_set]

#-----------------------------
# 2) Join baseline by (msa, base_day) and compute 螖rho
#-----------------------------
dt2 <- dt_dist %>%
  left_join(dt_base, by = c("msa" = "msa", "base_day" = "day")) %>%
  mutate(
    delta_rho     = rho - rho_base,
    pct_delta_rho = (rho - rho_base) / rho_base
  )
if (any(is.na(dt2$rho_base))) {
  bad <- dt2 %>% filter(is.na(rho_base)) %>% distinct(msa, base_day)
  print(bad)
  stop("Some (city, base_day) did not find baseline rho_base. Check date formats or baseline file.")
}

#-----------------------------
# 3) Summarize mean 卤 95% CI across replicates
#-----------------------------
sum_dt <- dt2 %>%
  group_by(msa, scenario, disturb) %>%
  summarise(
    n          = sum(!is.na(pct_delta_rho)),
    mean_delta = mean(pct_delta_rho, na.rm = TRUE),
    sd_delta   = sd(pct_delta_rho, na.rm = TRUE),
    se_delta   = sd_delta / sqrt(n),
    ci_low     = mean_delta - 1.96 * se_delta,
    ci_high    = mean_delta + 1.96 * se_delta,
    .groups = "drop"
  )

sum_dt$scenario <- factor(sum_dt$scenario, levels = scenario_set)

#-----------------------------
# 4) Plot Fig.1: 螖rho vs disturb with 95% CI ribbon
#-----------------------------
scenario_info <- data.frame(
  breaks = scenario_set,
  labels = c("Uniform scaling",
             "High-Low permutation",
             "Global permutation"),
  col = col.all[c(15,14,13)],
  shape = c( 21,22,24),
  stringsAsFactors = FALSE
)

fig_perturb <- ggplot(sum_dt, aes(x = disturb, y = mean_delta)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high, fill = scenario, group = scenario), alpha = 0.5, color = NA) +
  # scenario lines
  geom_line(aes(color = scenario, group = scenario),linewidth = 1) +
  # scenario points
  geom_point(aes(fill = scenario, shape = scenario),color = "black", size = 4, stroke = 0.5) +
  facet_wrap(~ msa, scales = "free_x", nrow = 2) +
  labs(x = "Perturbation level", y = TeX('Percentage change in $\\rho$ under perturbation')) +
  scale_x_continuous(breaks = disturb_set, labels = as.character(disturb_set)) +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  scale_color_manual(
    name   = "Scenario",
    breaks = scenario_info$breaks,
    labels = scenario_info$labels,
    values = scenario_info$col
  ) +
  scale_fill_manual(
    name   = "Scenario",
    breaks = scenario_info$breaks,
    labels = scenario_info$labels,
    values = scales::alpha(scenario_info$col, alpha = 1)
  ) +
  scale_shape_manual(
    name   = "Scenario",
    breaks = scenario_info$breaks,
    labels = scenario_info$labels,
    values = scenario_info$shape
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE), fill  = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain",size = 15),
        legend.justification = c(1, 1),
        legend.position = c(0.58, 0.75))
# ggsave(fig_perturb, filename = paste(figpath,"/msa/SI_rho_perturb_",Yname[yindex],".pdf",sep=""), width = 4*2, height = 4.5)





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
rho.rank[event_name == "Fire Kincade", event_name := "Kincade Fire"]
event.info$labels <- c(rep("COVID-19",Nmsa),"Hurricane Harvey","Hurricane Dorian","Winter Storm Uri","Kincade Fire")
fig_rho_change <- ggplot(rho.rank, aes(x = msa_rank, y = median_rho, color = event)) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  # geom_errorbar(aes(ymin = min_rho, ymax = max_rho), width = 0.2, linewidth = 0.5) +  
  geom_errorbar(aes(ymin = q25_rho, ymax = q75_rho), width = 0.2, linewidth = 0.5) +  
  geom_point(aes(fill = event, color = event), shape = 21, size = 4, stroke = 1) +
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
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 5, stroke = 0, colour = NA), nrow = 1), color = "none") +
  theme_wy() +
  theme(panel.border = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line.x = axis.arrow,
        axis.line.y = axis.arrow)

fig3 <- ((fig_rho_change/(fig_city|fig_normal))|fig_perturb) + plot_layout(widths = c(2.2, 1)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig3, filename = paste(figpath,"/msa/fig3.pdf",sep=""), width = 16, height = 5*2)







