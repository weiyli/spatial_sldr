# rm(list = ls())

# Relative uncertainty of rho
# profile loss L(蟻) curve


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(sf)          # read_sf() 
library(spdep)       # poly2nb() 
library(igraph)      # diameter()


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


#----------set the global function pointer used by homo/power/exp.param.est----------#
error_test <<- error_specs[["rmse_rmspe"]]$fun

#----------Load functions----------#
# Build day-level CI from rep-level uncertainty file
make_ci_from_uncert <- function(dt_rep, alpha = 0.05, keep_rep_dist = FALSE) {

  q_lo <- alpha/2
  q_hi <- 1 - alpha/2
  
  # dt_rep must contain: msa, day, rep, eta, rho, error, total_flow, lag
  dt_rep <- as.data.table(dt_rep)
  
  ci <- dt_rep[, .(
    n_rep        = .N,
    rho_median   = median(rho, na.rm=TRUE),
    rho_mean     = mean(rho, na.rm=TRUE),
    rho_sd       = sd(rho, na.rm=TRUE),
    rho_mad      = mad(rho, constant = 1, na.rm=TRUE),
    rho_lo       = as.numeric(quantile(rho, q_lo, na.rm=TRUE, type=7)),
    rho_hi       = as.numeric(quantile(rho, q_hi, na.rm=TRUE, type=7)),
    rho_ci_width = as.numeric(quantile(rho, q_hi, na.rm=TRUE, type=7) - quantile(rho, q_lo, na.rm=TRUE, type=7)),
    rho_iqr      = IQR(rho, na.rm=TRUE, type=7),
    err          = median(error, na.rm=TRUE),
    total_flow   = unique(total_flow)[1],
    lag          = unique(lag)[1]
  ), by = .(msa, day, eta)]
  
  ci[, rho_rel_ci := fifelse(is.finite(rho_median) & rho_median != 0, rho_ci_width / abs(rho_median), NA_real_)]
  setorder(ci, msa, eta, day)
  
  if (keep_rep_dist) return(list(ci = ci, rep = dt_rep))
  ci
}

#----------summarize to day-level CI----------#
region_use <- region[2] 
yname_use <- Yname[2]
msa_pick <- Dname.msa
rep_all <- rbindlist(lapply(msa_pick, function(msa){
  f <- file.path(datapath, region_use, "param", paste0(msa, "_uncertainty_", yname_use, ".csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  dt <- fread(f)
  dt[, day := as.Date(day)]
  dt
}), fill = TRUE)
rep_all <- rep_all[eta == 0.05]
ci_exp <- make_ci_from_uncert(dt_rep = rep_all, alpha = 0.05)
date.gif[, day_id := seq_len(.N)]
ci_exp <- merge(ci_exp, date.gif[, .(day, period, day_id)], by = "day", all.x = TRUE)
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day)],"%m-%d"),format(Datevalue[1:(Day)], "%a"),sep=" ")

#----------Relative uncertainty of rho----------#
fig_ci <- ggplot(ci_exp, aes(x = day_id, y = rho_median, group = msa)) +
  geom_ribbon(aes(ymin = rho_lo, ymax = rho_hi, fill = msa), alpha = 0.20, color = NA) +
  geom_line(aes(color = msa), linewidth = 0.8) +
  geom_point(aes(color = msa), size = 1.2) +
  facet_wrap(~ msa, scales = "free_y", nrow = 3) +
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
  labs(x = "Date", y = TeX("Uncertainty of $\\rho$ (95\\% CI)")) +
  scale_x_continuous(breaks = c(1,3,5,7,8,10,12,14), labels = id.week[c(1,3,5,7,8,10,12,14)]) +
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
# ggsave(fig_ci,filename = paste(figpath,"/msa/SI_rho_CI_",yname_use,".pdf",sep=""), width =4*4, height = 4*2.6)

#----------Relative CI width of rho----------#
fig_ci_rel <- ggplot(ci_exp, aes(x = msa, y = rho_rel_ci, fill = msa)) +
  geom_boxplot(width = 0.2, linewidth = 0.2, outlier.shape = NA, alpha = 0.6) +
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
  labs(x = "MSA", y = TeX("Relative 95\\% CI width of $\\rho$ (\\%)")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.01)) +
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
# ggsave(fig_ci_rel,filename = paste(figpath,"/msa/SI_rho_CI_rel_",yname_use,".pdf",sep=""), width =4*2, height = 4*1)
# fig_ci_all <- (fig_ci/fig_ci_rel) + plot_layout(heights = c(5, 2)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
# ggsave(fig_ci_all, filename = paste(figpath,"/msa/SI_rho_CI_",yname_use,".pdf",sep=""), width = 4*4, height = 4*4.2)
ggsave(fig_ci,filename = paste(figpath,"/msa/SI_rho_CI_",yname_use,".pdf",sep=""), width =4*4, height = 3.6*3)



#----------Min/Mean/Max 95% CI of rho----------#
#========================================================
# Build city-level CI-interval summary using CI-width
# (min / median / max CI-width -> corresponding [rho_lo, rho_hi])
#========================================================
# Safety: ensure needed columns exist
stopifnot(all(c("msa","day","rho_lo","rho_hi","rho_ci_width") %in% names(ci_exp)))

# Keep finite CI widths only
ci_use <- as.data.table(ci_exp)[
  is.finite(rho_ci_width) & !is.na(rho_ci_width) &
    is.finite(rho_lo) & is.finite(rho_hi) & !is.na(rho_lo) & !is.na(rho_hi)
]

# Helper: pick the row in one MSA that matches a target width
# rule: choose the day whose width is closest to target (ties -> earliest day)
pick_by_target_width <- function(dt_msa, target_width) {
  dt_msa[
    order(abs(rho_ci_width - target_width), day),
    .SD[1]
  ]
}

# For each MSA, select narrowest / median-width / widest
ci_interval_summary <- ci_use[, {
  # narrowest & widest are direct
  row_min <- .SD[which.min(rho_ci_width)]
  row_max <- .SD[which.max(rho_ci_width)]
  
  # "median case": choose the day whose CI width is closest to the median width
  w_med <- median(rho_ci_width, na.rm = TRUE)
  row_med <- pick_by_target_width(.SD, w_med)
  
  rbind(
    data.table(case = "Min CI", day = row_min$day, ci_width = row_min$rho_ci_width,
               rho_lo = row_min$rho_lo, rho_hi = row_min$rho_hi),
    data.table(case = "Median CI", day = row_med$day, ci_width = row_med$rho_ci_width,
               rho_lo = row_med$rho_lo, rho_hi = row_med$rho_hi),
    data.table(case = "Max CI", day = row_max$day, ci_width = row_max$rho_ci_width,
               rho_lo = row_max$rho_lo, rho_hi = row_max$rho_hi)
  )
}, by = msa]

# Optional: pretty CI interval string for LaTeX or tables
ci_interval_summary[, ci_95 := sprintf("$[%.4f,\\,%.4f]$", rho_lo, rho_hi)]

#========================================================
# Wide-format (1 row per city, 3 CI-interval columns)
#========================================================
ci_interval_wide <- dcast(
  ci_interval_summary,
  msa ~ case,
  value.var = "ci_95"
)

# If you also want to keep the days used:
ci_day_wide <- dcast(
  ci_interval_summary,
  msa ~ case,
  value.var = "day"
)

#========================================================
# (Optional) check output quickly
#========================================================
print(ci_interval_summary[order(msa, case)])
print(ci_interval_wide[order(msa)])


