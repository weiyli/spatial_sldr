# rm(list = ls())


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
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



#-------------------------
# Part 1: |螖rho_COVID|/inter-city variation in rho, |螖rho_COVID|/normal daily variation in rho
#-------------------------
#-------------------------
# 1) City-level mean rho by period and one-week COVID-induced changes (delta rho)
#-------------------------
yindex <- 2   
date.gif <- date.sep(Datevalue,durdate)
period_norm  <- pervalue[1]
period_covid <- pervalue[2]
rho_all <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], ".csv"))
rho_all <- as.data.table(rho_all)[, week := NULL]
rho_all <- rho_all[index == "exp"]
rho_all[, day := as.Date(day)]
Sys.setlocale("LC_TIME", "C") 
rho_all[, wday := weekdays(day)]

#-------------------------
# Analysis 1: |螖rho_COVID| vs inter-city variation in rho
#-------------------------
# |螖rho| per city |COVID - normal|
rho_wide <- dcast(rho_all, msa + wday ~ period, value.var = "rho")
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
# median(|螖rho|)
rho_delta_city <- rho_wide[,  .(abs_delta_rho_covid = median(abs_delta_rho)), by = msa]
# Inter-city SD of median rho in the normal period 
rho_city_norm <- rho_wide[, .(rho_city = median(before, na.rm = TRUE)), by = msa]
sd_across_cities <- sd(rho_city_norm$rho_city, na.rm = TRUE)
iqr_across_cities <- IQR(rho_city_norm$rho_city, na.rm = TRUE)
# Ratio per city (using SD across cities in normal period, as in your note)
rho_delta_city[, ratio_sd := abs_delta_rho_covid / sd_across_cities]
rho_delta_city[, ratio_iqr := abs_delta_rho_covid / iqr_across_cities]
# Plot: a compact bar/point plot for ratio (one number per city)
fig_city <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 4, stroke = 1) +
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
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[inter-city])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

mean(rho_delta_city$ratio_sd)


#-------------------------
# Analysis 2: |螖rho_COVID| vs null distribution of |螖rho| in normal period
#   螖rho_null(t)=rho_{t+1}-rho_t  (within each city, normal period)
#   螖rho_COVID(t)=rho_{t+1}-rho_t (within each city, COVID period)
#-------------------------
# |螖rho| per city |COVID - normal|
rho_wide <- dcast(rho_all, msa + wday ~ period, value.var = "rho")
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
# median(|螖rho|)
rho_delta_city <- rho_wide[,  .(abs_delta_rho_covid = median(abs_delta_rho)), by = msa]
# Within-city SD of rho during normal period 
rho_city_norm_sd <- rho_all[period == period_norm,
                            .(sd_rho_normal = sd(rho, na.rm = TRUE)),
                            by = msa]
# Merge and compute E = |螖rho_COVID| / SD(rho_normal)
rho_delta_city <- merge(rho_delta_city, rho_city_norm_sd, by = "msa", all.x = TRUE)
rho_delta_city[, ratio_sd := abs_delta_rho_covid / sd_rho_normal]
# plot: a compact bar/point plot for ratio (one number per city)
fig_normal <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 4, stroke = 1) +
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
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[normal])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

fig_comp <- (fig_city|fig_normal) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig_comp, filename = paste(figpath,"/msa/SI_rho_comp_",Yname[yindex],".pdf",sep=""), width = 7*2, height = 7)






#-------------------------
# Part 2: uncertainty of delta rho
#-------------------------
#-----------------------------
# MBB functions 
#-----------------------------
mbb_resample_index <- function(n, L) {
  stopifnot(n >= L)
  starts <- 1:(n - L + 1)
  n_blocks <- ceiling(n / L)
  s <- sample(starts, n_blocks, replace = TRUE)
  idx <- unlist(lapply(s, function(a) a:(a + L - 1)))
  idx[1:n]
}

mbb_ci <- function(x, stat_fun, L = 2, B = 2000, seed = 1) {
  # x must be time-ordered
  set.seed(seed)
  n <- length(x)
  if (n < L) stop("Length(x) must be >= L.")
  boots <- replicate(B, {
    idx <- mbb_resample_index(n, L)
    stat_fun(x[idx])
  })
  ci <- stats::quantile(boots, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
  list(ci = ci,
       boot_mean = mean(boots, na.rm = TRUE),
       boot_sd   = stats::sd(boots, na.rm = TRUE))
}


#-------------------------
# 1) City-level mean rho by period and five-weeks COVID-induced changes (delta rho)
#-------------------------
yindex <- 2   
duration <- 7*4 + 6
startdate <- '2020-01-13'
durdate <- '2020-03-30'
enddate <- '2020-05-03'
Datevalue <- c(seq(as.Date(startdate),as.Date(startdate) + duration, by="day"), seq(as.Date(durdate),as.Date(durdate) + duration, by="day"))
Day <- length(Datevalue) 
date.gif <- date.sep(Datevalue,durdate)
date.gif <- date.sep(Datevalue,durdate)
period_norm  <- pervalue[1]
period_covid <- pervalue[2]
rho_all <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], "_half_year.csv"))
rho_all <- as.data.table(rho_all)[, week := NULL]
rho_all <- rho_all[index == "exp"]
rho_all[, day := as.Date(day)]
rho_all <- rho_all[(period == "before" & day >= as.Date(startdate) & day <= (as.Date(startdate) + duration)) | 
                     (period == "during" & day >= as.Date(durdate) & day <= (as.Date(durdate) + duration))]
rho_all <- rho_all[, week_name := format(day, "%a")]

#-------------------------
# build week_id for 5-week matching: 1..5 for both before and during
#-------------------------
# before: week_id 1..5  
# during: week_id 1..5
rho_all[, week_id := floor(as.numeric(day - min(day)) / 7) + 1L, by = .(msa, period)]
rho_b <- rho_all[period == "before", .(msa, week_name, week_b = week_id, rho_b = rho)]
rho_d <- rho_all[period == "during", .(msa, week_name, week_d = week_id, rho_d = rho)]
rho_pair <- rho_b[rho_d, on = .(msa, week_name), allow.cartesian = TRUE]
rho_pair[, delta_rho := (rho_d - rho_b) / rho_b]
setorder(rho_pair, msa, week_name, week_b, week_d)

#-------------------------
# 2) Uncertainty of delta rho: MBB CI over the 5-point delta_rho series
#-------------------------
B_boot  <- 2000
L_block <- 4
seed0   <- 1
delta_ci_city <- rho_pair[, {
  d <- delta_rho
  d <- d[is.finite(d)]
  Tn <- length(d)     # typically 5
  if (Tn < L_block) {
    list(T_weeks = Tn, 
         delta_point = median(d, na.rm = TRUE), 
         ci_low = NA_real_, 
         ci_med = NA_real_, 
         ci_high = NA_real_)
  } else {
    mbb <- mbb_ci(d, stat_fun = function(z) median(z, na.rm = TRUE), L = L_block, B = B_boot, seed = seed0)
    list(T_weeks = Tn, 
         delta_point  = median(d, na.rm = TRUE), 
         ci_low = mbb$ci[1], 
         ci_med = mbb$ci[2],
         ci_high = mbb$ci[3])
  }
}, by = msa]
fwrite(delta_ci_city, paste0(datapath, "/msa/delta_rho_CI_", Yname[yindex], "_half_year.csv"))

#-------------------------
# 3) Figure: city-level 螖rho with 95% CI
#-------------------------
fig_delta_ci <- ggplot(delta_ci_city, aes(x = reorder(msa, ci_med), y = ci_med)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = msa), width = 0.2, linewidth = 0.6) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 4, stroke = 1) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
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
  scale_y_continuous(labels = scales::percent_format(scale = 100)) + 
  labs(x="MSA", y = TeX('Percentage change in $\\rho$ (95\\% CI)')) +
  theme_wy() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.position = "none",
        legend.justification = c(0, 0.5),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line.x = axis.arrow,
        axis.line.y = axis.arrow)
ggsave(fig_delta_ci,filename = paste(figpath,"/msa/SI_delta_rho_CI_",Yname[yindex],".pdf",sep=""), width = 10, height = 5)
# fig_delta_ci <- ggplot(delta_ci_city, aes(x = reorder(msa, ci_med), y = ci_med)) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "gray60") +
#   geom_errorbar(aes(ymin = ci_low, ymax = ci_high, color = msa), width = 0.2, linewidth = 0.8) +
#   geom_point(aes(fill = msa, color = msa, shape = msa), size = 4, stroke = 1) +
#   coord_flip() +
#   scale_fill_manual(name = "MSA",
#                     breaks = msa.info$dname,
#                     labels = msa.info$dname,
#                     values = scales::alpha(msa.info$col, alpha = 0.8)) +
#   scale_colour_manual(name = "MSA",
#                       breaks = msa.info$dname,
#                       labels = msa.info$dname,
#                       values = msa.info$col) +
#   scale_shape_manual(name = "MSA",
#                      breaks = msa.info$dname,
#                      labels = msa.info$dname,
#                      values = msa.info$shape) +
#   labs(x = NULL, y = expression(Delta * rho["COVID-19"]~"(median across 5 matched weeks, 95% CI)")) +
#   theme_wy() +
#   theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
#         panel.background = element_blank(),
#         panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
#         panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
#         legend.position = "none")









#-------------------------
# 1) City-level mean rho by period and five-weeks COVID-induced changes (delta rho)
#-------------------------
yindex <- 2   
duration <- 7*4 + 6
startdate <- '2020-01-13'
durdate <- '2020-03-30'
enddate <- '2020-05-03'
Datevalue <- c(seq(as.Date(startdate),as.Date(startdate) + duration, by="day"), seq(as.Date(durdate),as.Date(durdate) + duration, by="day"))
Day <- length(Datevalue) 
date.gif <- date.sep(Datevalue,durdate)
date.gif <- date.sep(Datevalue,durdate)
period_norm  <- pervalue[1]
period_covid <- pervalue[2]
rho_all <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], "_half_year.csv"))
rho_all <- as.data.table(rho_all)[, week := NULL]
rho_all <- rho_all[index == "exp"]
rho_all[, day := as.Date(day)]
rho_all <- rho_all[(period == "before" & day >= as.Date(startdate) & day <= (as.Date(startdate) + duration)) | 
                   (period == "during" & day >= as.Date(durdate) & day <= (as.Date(durdate) + duration))]
#-------------------------
# (prep) build week_pair for 5-week matching: 1..5 for both before and during
#-------------------------
# before: week_id 1..5  -> week_pair = 1..5
# during: week_id 6..10 -> week_pair = 1..5
rho_all[period == "before", week_id := floor(as.numeric(day - min(day)) / 7) + 1L]
rho_all[period == "during", week_id := floor(as.numeric(day - min(day)) / 7) + 6L]
rho_all[period == "before", week_pair := week_id]
rho_all[period == "during", week_pair := week_id - 5L]

#-------------------------
# Analysis 1: |螖rho_COVID| vs inter-city variation in rho  (5-week summary)
#-------------------------
# |螖rho| per city-weekpair |during - before|
rho_wide <- dcast(rho_all, msa + week_pair ~ period, value.var = "rho", fun.aggregate = median)
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
# median(|螖rho|) across the 5 matched weeks -> one number per city
rho_delta_city <- rho_wide[, .(abs_delta_rho_covid = median(abs_delta_rho, na.rm = TRUE)), by = msa]
# Inter-city SD/IQR of median rho in the normal (before) period
# (here: city-level median of the 5 pre weeks)
rho_city_norm <- rho_wide[, .(rho_city = median(before, na.rm = TRUE)), by = msa]
sd_across_cities  <- sd(rho_city_norm$rho_city, na.rm = TRUE)
iqr_across_cities <- IQR(rho_city_norm$rho_city, na.rm = TRUE)
# Ratio per city (using SD across cities in normal period, as in your note)
rho_delta_city[, ratio_sd  := abs_delta_rho_covid / sd_across_cities]
rho_delta_city[, ratio_iqr := abs_delta_rho_covid / iqr_across_cities]
# Plot: a compact bar/point plot for ratio (one number per city)
fig_city <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 4, stroke = 1) +
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
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[inter-city])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

#-------------------------
# Analysis 2: |螖rho_COVID| vs within-city variation in normal period (5-week summary)
#-------------------------
# 5-week paired |螖rho| (city-level median)
rho_wide <- dcast(rho_all, msa + week_pair ~ period, value.var = "rho", fun.aggregate = median)
rho_wide[, delta_rho := during - before]
rho_wide[, abs_delta_rho := abs(delta_rho)]
rho_delta_city <- rho_wide[, .(abs_delta_rho_covid = median(abs_delta_rho, na.rm = TRUE)), by = msa]
# Within-city SD of rho during normal period (day-level, pre window)
rho_city_norm_sd <- rho_all[period == period_norm,
                            .(sd_rho_normal = sd(rho, na.rm = TRUE)),
                            by = msa]
# Merge and compute E = |螖rho_COVID| / SD(rho_normal)
rho_delta_city <- merge(rho_delta_city, rho_city_norm_sd, by = "msa", all.x = TRUE)
rho_delta_city[, ratio_sd := abs_delta_rho_covid / sd_rho_normal]

fig_normal <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 4, stroke = 1) +
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
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[normal])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

fig_comp <- (fig_city | fig_normal) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig_comp, filename = paste(figpath, "/msa/SI_rho_comp_", Yname[yindex], "_half_year.pdf", sep = ""), width = 7*2, height = 7)






#-------------------------
# Part 3: Effect-size benchmarks using pairwise weekday matching for 5 weeks
#  - E_inter-city = median(|螖rho|) / SD( rho_inter-city baseline )
#  - E_normal     = median(|螖rho|) / SD( rho_normal temporal )
#-------------------------
#-------------------------
# 1) Read and subset rho (5-week before & during windows)
#-------------------------
yindex <- 2
duration  <- 7*4 + 6  # 5 weeks (35 days): 0..34
startdate <- "2020-01-13"
durdate   <- "2020-03-30"
rho_all <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], "_half_year.csv"))
rho_all <- as.data.table(rho_all)[, week := NULL]
rho_all <- rho_all[index == "exp"]
rho_all[, day := as.Date(day)]
rho_all <- rho_all[ (period == "before" & day >= as.Date(startdate) & day <= (as.Date(startdate) + duration)) |
                      (period == "during" & day >= as.Date(durdate)   & day <= (as.Date(durdate)   + duration))]
# weekday label (Mon/Tue/...)
Sys.setlocale("LC_TIME", "C")
rho_all[, week_name := format(day, "%a")]
# week_id: 1..5 within each (msa, period)
rho_all[, week_id := floor(as.numeric(day - min(day)) / 7) + 1L, by = .(msa, period)]

#-------------------------
# 2) Pairwise (cartesian) matching within weekday: before week_b vs during week_d
#-------------------------
rho_b <- rho_all[period == "before", .(msa, week_name, week_b = week_id, rho_b = rho)]
rho_d <- rho_all[period == "during", .(msa, week_name, week_d = week_id, rho_d = rho)]

# cartesian join: for each (msa, weekday), all 5x5 week combinations
rho_pair <- rho_b[rho_d, on = .(msa, week_name), allow.cartesian = TRUE]

# delta rho (relative change)
rho_pair[, delta_rho := (rho_d - rho_b) / rho_b]
rho_pair <- rho_pair[is.finite(delta_rho)]

# city-level shock strength (one number per city): median(|螖rho|)
rho_delta_city <- rho_pair[, .(
  abs_delta_rho_covid = median(abs(delta_rho), na.rm = TRUE)
), by = msa]

#-------------------------
# 3) Denominator 1: inter-city variation of baseline rho (before window)
#-------------------------
# baseline per city: median rho in BEFORE window (day-level)
rho_city_baseline <- rho_all[period == "before", .(rho_city = median(rho, na.rm = TRUE)), by = msa]
sd_across_cities  <- sd(rho_city_baseline$rho_city, na.rm = TRUE)
iqr_across_cities <- IQR(rho_city_baseline$rho_city, na.rm = TRUE)

rho_delta_city[, ratio_sd  := abs_delta_rho_covid / sd_across_cities]
rho_delta_city[, ratio_iqr := abs_delta_rho_covid / iqr_across_cities]

#-------------------------
# 4) Denominator 2: within-city temporal variability in normal period (before window)
#-------------------------
rho_city_norm_sd <- rho_all[period == "before",
                            .(sd_rho_normal = sd(rho, na.rm = TRUE)),
                            by = msa]

rho_delta_city <- merge(rho_delta_city, rho_city_norm_sd, by = "msa", all.x = TRUE)
rho_delta_city[, ratio_sd_normal := abs_delta_rho_covid / sd_rho_normal]

#-------------------------
# 5) Figures (same style as your original)
#-------------------------
fig_city <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd), y = ratio_sd)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 4, stroke = 1) +
  scale_fill_manual(name = "MSA",
                    breaks = msa.info$dname, labels = msa.info$dname,
                    values = scales::alpha(msa.info$col, alpha = 0.8)) +
  scale_colour_manual(name = "MSA",
                      breaks = msa.info$dname, labels = msa.info$dname,
                      values = msa.info$col) +
  scale_shape_manual(name = "MSA",
                     breaks = msa.info$dname, labels = msa.info$dname,
                     values = msa.info$shape) +
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[inter-city])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")

fig_normal <- ggplot(rho_delta_city, aes(x = reorder(msa, ratio_sd_normal), y = ratio_sd_normal)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.8, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 4, stroke = 1) +
  scale_fill_manual(name = "MSA",
                    breaks = msa.info$dname, labels = msa.info$dname,
                    values = scales::alpha(msa.info$col, alpha = 0.8)) +
  scale_colour_manual(name = "MSA",
                      breaks = msa.info$dname, labels = msa.info$dname,
                      values = msa.info$col) +
  scale_shape_manual(name = "MSA",
                     breaks = msa.info$dname, labels = msa.info$dname,
                     values = msa.info$shape) +
  coord_flip() +
  labs(x = NULL, y = expression(frac("|" * Delta * rho[COVID-19] * "|", SD(rho[normal])))) +
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
        legend.position = "none")
fig_comp <- (fig_city | fig_normal) + patchwork::plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
ggsave(fig_comp, filename = paste0(figpath, "/msa/SI_rho_comp_", Yname[yindex], "_half_year.pdf"), width = 7*2, height = 7)

