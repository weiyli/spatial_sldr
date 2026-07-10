# rm(list = ls())

# Dataset description & bias


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
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
date.gif <- date.sep(Datevalue, durdate)%>%as.data.table
date.gif[, day_id := seq_len(.N)]
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day)],"%m-%d"),format(Datevalue[1:(Day)], "%a"),sep=" ")
# x_breaks  <- c(1,3,5,7,8,10,12,14)
x_breaks <- c(1:14)
x_labels  <- id.week[x_breaks]


#=============================#
# 1) Build CBG -> MSA mapping #
#=============================#
cbg2msa_list <- lapply(Dname.msa, function(msa_name){
  f <- file.path(geopath, "msa", paste0(msa_name, ".geojson"))
  block <- sf::read_sf(f)
  dt <- as.data.table(sf::st_drop_geometry(block))
  dt <- unique(
    dt[, .(
      origin_census_block_group = as.numeric(CensusBlockGroup),
      msa = msa_name,
      county_fips = as.character(county_fips)
    )],
    by = "origin_census_block_group"
  )
  dt
})
cbg2msa <- rbindlist(cbg2msa_list)
# sanity
stopifnot(!anyDuplicated(cbg2msa$origin_census_block_group))


#===============================#
# 2) Load ACS (CBG-level) demo  #
#===============================#
demo <- fread(paste0(geopath, "/census/acs_2019_5years_cbg.csv"))
acs_cbg <- as.data.table(demo)[, .(
  origin_census_block_group = as.numeric(cbg_fips),
  total_population = as.numeric(total_population),
  
  median_household_income = as.numeric(median_household_income),
  per_capita_income = as.numeric(per_capita_income),
  
  white_population = as.numeric(white_population),
  black_population = as.numeric(black_population),
  asian_population = as.numeric(asian_population),
  hispanic_population = as.numeric(hispanic_population),
  
  age_under_18 = as.numeric(age_under_18),
  age_65_and_over = as.numeric(age_65_and_over)
)]

acs_cbg[, share_white    := ifelse(total_population > 0, white_population / total_population, NA_real_)]
acs_cbg[, share_black    := ifelse(total_population > 0, black_population / total_population, NA_real_)]
acs_cbg[, share_asian    := ifelse(total_population > 0, asian_population / total_population, NA_real_)]
acs_cbg[, share_hispanic := ifelse(total_population > 0, hispanic_population / total_population, NA_real_)]
acs_cbg[, share_under18  := ifelse(total_population > 0, age_under_18 / total_population, NA_real_)]
acs_cbg[, share_over65   := ifelse(total_population > 0, age_65_and_over / total_population, NA_real_)]


#=====================================================#
# 3) Loop over SafeGraph social-distancing (by Day)    #
#    Outputs: A1, A2, A5 + CBG daily devices           #
#=====================================================#
A1_list  <- vector("list", Day)
A2_list  <- vector("list", Day)
A5_list  <- vector("list", Day)
CBG_list <- vector("list", Day)

for(i in 1:Day){
  
  df <- read.csv(paste(flowpath,"/", year(Datevalue[i]),"/", PDatevalue[i],"/",Datevalue[i], "-social-distancing.csv.gz",sep=""))
  df <- as.data.table(df)[, .(
    origin_census_block_group = as.numeric(origin_census_block_group),
    date_range_start = as.IDate(date_range_start),
    device_count = as.numeric(device_count),
    candidate_device_count = as.numeric(candidate_device_count),
    median_non_home_dwell_time = as.numeric(median_non_home_dwell_time)
  )]
  
  df_msa <- merge(df, cbg2msa, by = "origin_census_block_group", all = FALSE)
  
  #----------A1: Active devices----------#
  A1_list[[i]] <- df_msa[, .(
    active_devices = sum(candidate_device_count, na.rm = TRUE),
    devices = sum(device_count, na.rm = TRUE)
  ), by = .(date_range_start, msa)]
  
  #----------A2: Activity proxy----------#
  # A2_list[[i]] <- df_msa[, .(
  #   total_devices = sum(device_count, na.rm = TRUE),
  #   total_candidates = sum(candidate_device_count, na.rm = TRUE),
  #   dwell_proxy = sum(device_count * median_non_home_dwell_time, na.rm = TRUE)
  # ), by = date_range_start]
  A2_list[[i]] <- df_msa[, .(
    total_candidates = sum(candidate_device_count, na.rm = TRUE),
    total_devices    = sum(device_count, na.rm = TRUE),
    active_ratio     = sum(device_count, na.rm = TRUE) /
      sum(candidate_device_count, na.rm = TRUE),
    total_nonhome_dwell = sum(device_count * median_non_home_dwell_time, na.rm = TRUE)
  ), by = .(date_range_start, msa)]
  
  #----------A5: Presence----------#
  A5_list[[i]] <- unique(df_msa[, .(origin_census_block_group, msa, date_range_start)])
  
  #----------CBG daily (for A3/A4)----------#
  CBG_list[[i]] <- df_msa[, .(
    origin_census_block_group,
    msa,
    date_range_start,
    device_count,
    candidate_device_count
  )]
}

dev_time_msa  <- rbindlist(A1_list)
activity_time <- rbindlist(A2_list)
cbg_presence  <- rbindlist(A5_list)
cbg_daily     <- rbindlist(CBG_list)


#==============================================#
# 4) A1 Plot: active devices time series (MSA) #
#==============================================#
# dev_time_msa[, active_devices_index := active_devices / median(active_devices, na.rm = TRUE), by = msa]
dev_time_msa[, day := as.Date(date_range_start)]
dev_time_msa <- merge(
  dev_time_msa,
  date.gif[, .(day , period, day_id)],
  by = "day",
  all.x = TRUE
)
fig_A1_dev <- ggplot(dev_time_msa, aes(x = day_id, y = active_devices, group = msa)) +
  geom_line(aes(color = msa), linewidth = 0.8) +
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
  labs(x = NULL, y = "The number of active devices") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01,0.01)) +
  scale_y_continuous(labels = scales::comma) + 
  guides(color = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2),
         shape = guide_legend(nrow = 2)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
 
#==============================================#
# 5) A2 Plot: activity proxy time series       #
#==============================================#
activity_time <- as.data.table(activity_time)
activity_time[, day := as.Date(date_range_start)]
activity_time <- merge(
  activity_time,
  date.gif[, .(day, period, day_id)],
  by = "day",
  all.x = TRUE
)
activity_time[, period := str_to_title(period)]
#------------------------------#
# Fig A2-1: total candidates   #
#------------------------------#
fig_A2_candidates <- ggplot(activity_time, aes(x = day_id, y = total_candidates, group = msa)) +
  geom_line(aes(color = msa), linewidth = 0.8) +
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
  labs(x = NULL, y = "Total candidate devices") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01,0.01)) +
  scale_y_continuous(labels = scales::comma) + 
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")

#------------------------------#
# Fig A2-2: active ratio (best proxy for intensity)   #
#------------------------------#
fig_A2_ratio <- ggplot(activity_time, aes(x = day_id, y = active_ratio, group = msa)) +
  geom_line(aes(color = msa), linewidth = 0.8) +
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
  labs(x = NULL, y = "Active-to-candidate ratio (%)") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0.01,0.01)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.3, 0.6)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
   
  # theme(panel.border = element_blank(),
  #       panel.background = element_blank(),
  #       legend.position = "right",
  #       axis.text.x = element_text(angle = 30, hjust = 1),
  #       axis.line.x = axis.arrow,
  #       axis.line.y = axis.arrow)

#====================================================#
# 6) A3/A4 + B1: penetration + representativeness     #
#====================================================#
# CBG mean devices over the whole window (you can restrict to baseline later)
cbg_daily[, day := as.Date(date_range_start)]
cbg_daily <- merge(
  cbg_daily,
  date.gif[, .(day, period, day_id)],
  by = "day",
  all.x = TRUE
)
cbg_daily[, period := str_to_title(period)]
cbg_mean <- cbg_daily[, .(
  mean_devices = mean(device_count, na.rm = TRUE),
  mean_candidates = mean(candidate_device_count, na.rm = TRUE)
), by = .(msa, origin_census_block_group, period)]
df_pen <- merge(cbg_mean, acs_cbg, by = "origin_census_block_group", all.x = TRUE)
df_pen <- df_pen[!is.na(total_population) & total_population > 0]
df_pen[, penetration_per_1000 := mean_devices / total_population * 1000]
# A3(i) penetration distribution by MSA (barplot)
cov_ci_msa <- df_pen[, .(
  pen_med  = median(penetration_per_1000, na.rm = TRUE),
  pen_low  = quantile(penetration_per_1000, 0.25, na.rm = TRUE),  # Q1
  pen_high = quantile(penetration_per_1000, 0.75, na.rm = TRUE)   # Q3
), by = .(msa, period)]
cov_ci_msa[, msa_ord := reorder(msa, pen_med)]
fig_A3_across <- ggplot(cov_ci_msa, aes(x = msa_ord, y = pen_med)) +
  geom_errorbar(aes(ymin = pen_low, ymax = pen_high, color = msa),
                width = 0.2, linewidth = 0.6) +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 4, stroke = 1) +
  facet_wrap(~ period, scales = "fixed", nrow = 3) +
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
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Devices per 1,000 residents (CBG level)") +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )
# A3(ii) within-city heterogeneity
pen_heter_msa <- df_pen[, .(
  pen_mean = mean(penetration_per_1000, na.rm = TRUE),
  pen_sd   = sd(penetration_per_1000, na.rm = TRUE),
  pen_cv   = sd(penetration_per_1000, na.rm = TRUE) /
    mean(penetration_per_1000, na.rm = TRUE),
  pen_iqr  = IQR(penetration_per_1000, na.rm = TRUE),
  pen_med  = median(penetration_per_1000, na.rm = TRUE)
), by = .(msa, period)]

pen_heter_msa <- df_pen[, {
  q1 <- quantile(penetration_per_1000, 0.25, na.rm = TRUE)
  q3 <- quantile(penetration_per_1000, 0.75, na.rm = TRUE)
  x_trim <- penetration_per_1000[
    penetration_per_1000 >= q1 &
      penetration_per_1000 <= q3
  ]
  
  list(
    pen_mean = mean(x_trim, na.rm = TRUE),
    pen_sd   = sd(x_trim, na.rm = TRUE),
    pen_cv  = sd(x_trim, na.rm = TRUE) /
      mean(x_trim, na.rm = TRUE),
    pen_iqr       = q3 - q1,
    pen_med       = median(penetration_per_1000, na.rm = TRUE),
    n       = length(x_trim)
  )
}, by = .(msa, period)]
fig_A3_within <- ggplot(pen_heter_msa, aes(x = reorder(msa, pen_cv), y = pen_cv)) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 4, stroke = 1) +
  facet_wrap(~ period, scales = "fixed", nrow = 3) +
  scale_fill_manual(values = scales::alpha(msa.info$col, 0.8),
                    breaks = msa.info$dname, labels = msa.info$dname) +
  scale_colour_manual(values = msa.info$col,
                      breaks = msa.info$dname, labels = msa.info$dname) +
  scale_shape_manual(values = msa.info$shape,
                     breaks = msa.info$dname, labels = msa.info$dname) +
  labs(x = NULL, y = "Within-city coverage heterogeneity (CV)") +
  scale_y_continuous(expand = c(0.01,0.01)) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

# A4: representativeness correlations
rep_test <- df_pen[, .(
  cor_income = cor(penetration_per_1000, median_household_income, use = "complete.obs"),
  cor_white  = cor(penetration_per_1000, share_white, use = "complete.obs"),
  cor_black  = cor(penetration_per_1000, share_black, use = "complete.obs"),
  cor_hisp   = cor(penetration_per_1000, share_hispanic, use = "complete.obs"),
  cor_under18 = cor(penetration_per_1000, share_under18, use = "complete.obs"),
  cor_over65  = cor(penetration_per_1000, share_over65, use = "complete.obs")
),  by = .(msa, period)]
rep_long <- melt(
  rep_test,
  id.vars = c("msa", "period"),
  variable.name = "metric",
  value.name = "r"
)
metric_map <- data.frame(
  metric = c("cor_income", "cor_white", "cor_black", "cor_hisp", "cor_under18", "cor_over65"),
  metric_label = c("Income", "White", "Black", "Hispanic", "Under 18", "Age 65+")
)
rep_long <- merge(rep_long, metric_map, by = "metric", all.x = TRUE)
rep_long[, metric_label := factor(metric_label, levels = metric_map$metric_label)]
fig_A4_cor <- ggplot(rep_long, aes(x = msa, y = metric_label, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  facet_wrap(~ period, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-0.5, 0.5),
    oob = squish
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "Demographic and socio-economic attributes", fill = "Correlation") +
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
pct_abs_multi <- rep_long[!is.na(r), .(
  pct_lt01 = 100 * mean(abs(r) < 0.1),
  pct_lt02 = 100 * mean(abs(r) < 0.2),
  pct_lt03 = 100 * mean(abs(r) < 0.3),
  n_total  = .N
), by = period]


fig_A1_dev
fig_A2_ratio
fig_A3_across
fig_A3_within
fig_A4_cor

fig.data.bias <- ((fig_A1_dev|fig_A2_ratio)/(fig_A3_across|fig_A3_within)/fig_A4_cor) + plot_layout(heights = c(1.5, 2, 1.5)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
ggsave(fig.data.bias, filename = paste(figpath,"/msa/SI_data_bias_",Yname[2],".pdf",sep=""), width =4*4, height = 4*4)



# # A4: quintile diagnostic
# df_pen[, pen_rank := frank(penetration_per_1000, ties.method = "average"), by = .(msa, period)]
# df_pen[, pen_q := ceiling(5 * pen_rank / max(pen_rank, na.rm = TRUE)), by = .(msa, period)]
# rep_bins <- df_pen[, .(
#   mean_income   = mean(median_household_income, na.rm = TRUE),
#   mean_white    = mean(share_white, na.rm = TRUE),
#   mean_black    = mean(share_black, na.rm = TRUE),
#   mean_hisp     = mean(share_hispanic, na.rm = TRUE),
#   mean_under18  = mean(share_under18, na.rm = TRUE),
#   mean_over65   = mean(share_over65, na.rm = TRUE),
#   n_cbg         = .N
# ), by = .(msa, period, pen_q)]
# fig_rep_income_q <- ggplot(rep_bins, aes(x = pen_q, y = mean_income, group = msa)) +
#   geom_line(aes(color = msa), linewidth = 0.6, alpha = 0.9) +
#   geom_point(aes(fill = msa, color = msa, shape = msa), size = 2.2, stroke = 0.5) +
#   facet_wrap(~ period, nrow = 1) +
#   scale_colour_manual(values = msa.info$col, breaks = msa.info$dname, labels = msa.info$dname) +
#   scale_fill_manual(values = scales::alpha(msa.info$col, 0.8), breaks = msa.info$dname, labels = msa.info$dname) +
#   scale_shape_manual(values = msa.info$shape, breaks = msa.info$dname, labels = msa.info$dname) +
#   scale_x_continuous(breaks = 1:5, labels = paste0("Q", 1:5)) +
#   scale_y_continuous(labels = comma) +
#   labs(x = "Penetration quintile (CBGs)", y = "Mean ACS median household income") +
#   theme_wy() +
#   theme(
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     legend.position = "none"
#   )


#====================================================#
# 7) A5: coverage stability / device churn (CBG presence)    #
#====================================================#
#------------------------------#
# CBG-level coverage stability (churn proxy)
#------------------------------#
cbg_presence[, day := as.Date(date_range_start)]
cbg_presence <- merge(
  cbg_presence,
  date.gif[, .(day, period, day_id)],
  by = "day",
  all.x = TRUE
)
cbg_presence[, period := str_to_title(period)]
# Number of unique observed days per msa×period
period_days <- cbg_presence[, .(
  n_period_days = uniqueN(day)
), by = .(msa, period)]

# CBG presence: how many days each CBG appears within each msa×period
cbg_stability <- cbg_presence[, .(
  n_days = uniqueN(day)
), by = .(msa, period, origin_census_block_group)]

cbg_stability <- merge(
  cbg_stability,
  period_days,
  by = c("msa", "period"),
  all.x = TRUE
)

cbg_stability[, coverage_rate := n_days / n_period_days]

#------------------------------#
# MSA-level summaries for reporting
#------------------------------#
cbg_stability_msa <- cbg_stability[, .(
  mean_days = mean(n_days),
  p10_days  = quantile(n_days, 0.10),
  p50_days  = quantile(n_days, 0.50),
  p90_days  = quantile(n_days, 0.90),
  
  mean_cov  = mean(coverage_rate),
  p10_cov   = quantile(coverage_rate, 0.10),
  p50_cov   = quantile(coverage_rate, 0.50),
  p90_cov   = quantile(coverage_rate, 0.90),
  
  frac_cov_lt80 = mean(coverage_rate < 0.80),
  frac_cov_lt50 = mean(coverage_rate < 0.50)
), by = .(msa, period)]

#------------------------------#
# Candidate device pool stability (optional but recommended)
#------------------------------#
pool_msa <- NULL
pool_msa_sum <- NULL

has_candidate <- "candidate_device_count" %in% names(cbg_presence)
has_device    <- "device_count" %in% names(cbg_presence)

if (has_candidate || has_device) {
  
  pool_msa <- cbg_presence[, .(
    total_candidates = if (has_candidate) sum(candidate_device_count, na.rm = TRUE) else NA_real_,
    total_devices    = if (has_device)    sum(device_count, na.rm = TRUE)          else NA_real_
  ), by = .(msa, period, day)]
  
  pool_msa_sum <- pool_msa[, .(
    mean_candidates = if (has_candidate) mean(total_candidates, na.rm = TRUE) else NA_real_,
    cv_candidates   = if (has_candidate) sd(total_candidates, na.rm = TRUE) / mean(total_candidates, na.rm = TRUE) else NA_real_,
    p10_candidates  = if (has_candidate) quantile(total_candidates, 0.10, na.rm = TRUE) else NA_real_,
    p90_candidates  = if (has_candidate) quantile(total_candidates, 0.90, na.rm = TRUE) else NA_real_,
    
    mean_devices    = if (has_device) mean(total_devices, na.rm = TRUE) else NA_real_,
    cv_devices      = if (has_device) sd(total_devices, na.rm = TRUE) / mean(total_devices, na.rm = TRUE) else NA_real_
  ), by = .(msa, period)]
}

stable_cbg <- cbg_stability[coverage_rate >= 0.80, .(msa, period, origin_census_block_group)]
# You can use stable_cbg to subset df_pen / mobility analysis later.

# (Fig A5-1) CBG presence distribution (n_days), faceted by period
fig_A5_presence <- ggplot(cbg_stability, aes(x = n_days)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ period, nrow = 1, scales = "free_y") +
  labs(x = "Number of observed days (CBG presence)", y = "Number of CBGs") +
  theme_bw()

# (Fig A5-2) Coverage rate distribution
fig_A5_covrate <- ggplot(cbg_stability, aes(x = coverage_rate)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ period, nrow = 1, scales = "free_y") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Coverage rate within period (n_days / n_period_days)", y = "Number of CBGs") +
  theme_bw()

# (Fig A5-3) MSA-level churn summary: median coverage with 10–90% band
fig_A5_msa_cov <- ggplot(cbg_stability_msa, aes(x = reorder(msa, p50_cov), y = p50_cov)) +
  geom_errorbar(aes(ymin = p10_cov, ymax = p90_cov), width = 0.2, linewidth = 0.6) +
  geom_point(size = 2.5) +
  coord_flip() +
  facet_wrap(~ period, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "MSA", y = "CBG coverage rate (median; 10th–90th percentile)") +
  theme_bw()

# (Fig A5-4) Candidate device pool stability (if available)
fig_A5_pool <- NULL
if (!is.null(pool_msa) && has_candidate) {
  fig_A5_pool <- ggplot(pool_msa, aes(x = day, y = total_candidates, group = msa)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ period, nrow = 1, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    labs(x = "Date", y = "Total candidate devices (MSA-level sampling pool)") +
    theme_bw()
}
