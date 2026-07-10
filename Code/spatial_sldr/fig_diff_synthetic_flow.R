# rm(list = ls())

# Interpreting the gap between null and empirical rho
# [ rho vs travel distance ]   [ rho vs travel mode ]   [ rho vs POI entropy ]
# Null and empirical rho from synthetic_flow.R; msa level data from sldr_fit.R; variable from rho_geo_var.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(scales)

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


#----------rho and msa or county----------#
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")


#============================================================
# Part 1. Variables
#============================================================
#----------rho and var----------#
geo_var <- fread(paste0(datapath,"/msa/geo_var.csv",sep=""))

#----------pop and area for county and msa----------#
queen.lag.msa <- fread(file.path(datapath, "msa", "queen_lag.csv"))
queen.lag.msa[, pop  := pop  / 1e6]  # million
queen.lag.msa[, area := area / 1e6]  # km^2
queen.max.msa <- unique(queen.lag.msa[, .(msa, nblock, lag.max, pop, area)])

#----------power-law scaling model: log(rho)=rho0+alpha*log(pop)----------#
yindex <- 2 
#----------msa level: data from sldr_fit.R and od_dist.R----------#
rho.msa <- fread(paste0(datapath,"/msa/SLDR_params_",Yname[yindex],".csv",sep=""))
rho.msa <- rho.msa[index == "exp"&period==pervalue[1]]
rho.msa <- merge(rho.msa, queen.max.msa, by=c('msa'))
var.msa <- rho.msa[, .(msa, day, rho, pop, area)]
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
# aggregate and join by (msa, day) 
dist.msa.all <- as.data.table(dist.msa.all)
dist.msa.all[, weighted_haversine := weighted_haversine / 1000]  # m -> km
dist.msa <- dist.msa.all[, .(
  mean_haversine    = mean(weighted_haversine, na.rm = TRUE),
  mean_destinations = mean(n_destinations, na.rm = TRUE)
), by = .(msa, day)]
rho.dist <- left_join(var.msa, dist.msa, by = c("msa", "day"))
#----------Connecting rho with travel mode----------#
rho_var <- left_join(rho.dist, geo_var, by="msa")
travel.mode <- read.csv(paste(geopath, "/census/acs_2019_5years_travel_mode_msa.csv", sep=""), header=TRUE)
travel.mode <- left_join(travel.mode, data.frame(msa_fips = msa.info$did, msa = msa.info$dname), by = "msa_fips")
rho_var <- left_join(rho_var, travel.mode, by = c("msa"))
vars_norm <- c(
  "total_commute",
  "drive_alone",
  "carpool",
  "public_transit",
  "bicycle",
  "walk",
  "other_means",
  "work_at_home"
)
vars_norm <- vars_norm[vars_norm %in% names(rho_var)]
rho_var[, (vars_norm) := lapply(.SD, function(x) x / pop_total), .SDcols = vars_norm]

rho_var[, c("pop_total", "area_km2") := NULL]
cols_keep <- c(
  "msa",
  "day",
  "county",
  "rho",
  "pop",
  "area",
  "avg_density",
  "mean_haversine",
  # "mean_destinations",
  "radial_slope",
  # "radial_r2",
  "compactness",
  "poly_gini",
  "rank_size_slope",
  
  "carpool",
  "drive_alone",
  "work_at_home",
  "walk",
  "public_transit",
  "bicycle"
)
rho_sub <- rho_var[, ..cols_keep]
rho_long <- melt(
  rho_sub,
  id.vars = c("msa","county", "day", "rho"),
  variable.name = "indicator",
  value.name = "x"
)
rho_dt <- as.data.table(rho_long)
rho_wide <- dcast(
  rho_dt,
  msa + county + day + rho ~ indicator,
  value.var = "x"
)


#----------Analysis of variable contributions----------#
# (1) Structural / geometric controls
vars_A <- c("pop", "area", "avg_density")
# (2) Mobility / interaction variables
vars_B <- c("mean_haversine","drive_alone", "public_transit")
# (3) Metropolitan-scale morphology (internal organization; MSA-only)
vars_C  <- c("radial_slope", "poly_gini", "rank_size_slope") 
# all
vars <- c(vars_A, vars_B, vars_C)
# scale var
safe_scale <- function(x) {
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
to_scale <- vars[vars %in% names(rho_wide)]
rho_wide[, paste0(to_scale, "_z") := lapply(.SD, safe_scale), .SDcols = to_scale]





#============================================================
# Part 2. Null rho and Empirical rho
#============================================================
region_name <- region[2]    
y_tag <- Yname[2]
dset <- 1:Nmsa
# empirical
rho_emp <- fread(paste0(datapath, "/",region_name,"/SLDR_params_", y_tag, ".csv"))
rho_emp <- rho_emp[index == "exp" & msa%in%Dname.set[dset]]
# null
read_null_rho <- function(city_name, region_name, y_tag, datapath,
                          folder = c("param", "params"),
                          suffix = "_SLDR_params_synthetic_",
                          required_cols = c("day", "rho", "model")) {
  folder <- match.arg(folder)
  f1 <- file.path(datapath, region_name, folder, paste0(city_name, suffix, y_tag, ".csv"))
  if (!file.exists(f1)) stop("Null rho file not found: ", f1)
  d1 <- fread(f1)
  d1[, msa := city_name]
  d1[, .(msa, day, model, rho)]   
}
null_list <- list()
for (d in dset) {
  city_name <- Dname.set[d]
  
  null_list[[city_name]] <- read_null_rho(
    city_name = city_name,
    region_name = region_name,
    y_tag = y_tag,
    datapath = datapath,
    folder = "param" 
  )
}
rho_null <- rbindlist(null_list, fill = TRUE)
# merge day-level empirical & null
rho_cmp <- merge(
  rho_emp[, .(msa, day, rho_emp = rho)],
  rho_null[, .(msa, day, model, rho_null = rho)],
  by = c("msa", "day")
)
rho_cmp[, model := stringr::str_to_title(model)]




#============================================================
# Part 3. Interpreting the gap between null and empirical rho
#============================================================
rho_all <- merge(
  rho_cmp,
  rho_wide,
  by = c("msa", "day")
)
rho_all[, delta_rho := rho_emp - rho_null]

var_labels <- c(
  pop = "Population",
  area = "Urban area",
  avg_density = "Population density",
  mean_haversine = "Average travel distance",
  drive_alone = "Drive-alone share",
  public_transit = "Public transit share",
  radial_slope = "Radial density gradient",
  poly_gini = "Polycentricity",
  rank_size_slope = "Rank-size slope"
)

cor_dt <- rbindlist(lapply(vars, function(v) {
  rho_all[, .(
    r = cor(delta_rho, get(v),
            use = "pairwise.complete.obs",
            method = "pearson")
  ), by = .(model)][, variable := v]
}))

cor_dt[, metric_label := var_labels[variable]]

cor_dt[, metric_label := factor(
  metric_label,
  levels = rev(unname(var_labels))
)]


fig_delta_cor <- ggplot(cor_dt, aes(x = model, y = metric_label, fill = r)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = c(-0.5, 0.5),
    oob = squish
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Variables",
    fill = "Correlation"
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





#----------rho diff with each indicator for msa----------#
#============================================================
# 1. Prepare msa-level data for plotting
#============================================================
vars_plot <- c("pop", "public_transit", "poly_gini")
indicator_labels <- c(
  pop            = "Population~(million)",
  public_transit = "Public~transit~share",
  poly_gini = "Polycentricity"
)

rho_plot_msa <- rho_all[, .(
  delta_rho       = median(delta_rho, na.rm = TRUE),
  pop             = median(pop, na.rm = TRUE),
  public_transit  = median(public_transit, na.rm = TRUE),
  poly_gini    = median(poly_gini, na.rm = TRUE)
), by = .(msa, model)]

plot_dt <- melt(
  rho_plot_msa,
  id.vars = c("msa", "model", "delta_rho"),
  measure.vars = vars_plot,
  variable.name = "indicator",
  value.name = "x"
)

plot_dt[, indicator := factor(indicator, levels = vars_plot)]

#============================================================
# 2. Compute R2 / P / slope for each panel
#============================================================

fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "< 0.001",
                sprintf("= %.3f", p)))
}

stats_dt <- plot_dt[
  is.finite(x) & is.finite(delta_rho),
  {
    fit <- lm(delta_rho ~ x)
    s   <- summary(fit)
    
    xr <- range(x, na.rm = TRUE)
    yr <- range(delta_rho, na.rm = TRUE)
    
    list(
      R2    = s$r.squared,
      pval  = coef(s)[2, "Pr(>|t|)"],
      slope = coef(fit)[2],
      x_pos = xr[1] + 0.05 * diff(xr),
      y_pos = yr[1] + 0.10 * diff(yr)
    )
  },
  by = .(model, indicator)
]

stats_dt[, label := sprintf(
  "atop(italic(R)^2 == %.2f*','~~italic(P)~'%s', italic(slope) == %.2f)",
  R2, fmt_p(pval), slope
)]

stats_dt[, indicator := factor(indicator, levels = vars_plot)]

#============================================================
# 3. Plot
#============================================================
fig_delta_var <- ggplot(plot_dt, aes(x = x, y = delta_rho)) +
  geom_abline(slope = 0, intercept = 0,linetype = "dashed", linewidth = 0.6, color = "grey65") +
  geom_point(aes(fill = msa, color = msa, shape = msa), size= 4, stroke = 1) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 1.15,
    color = "#69b3a2",
    fill = "grey80",
    alpha = 0.35
  ) +
  geom_text(
    data = stats_dt,
    aes(label = label),
    x = -Inf,   
    y = Inf,   
    hjust = -0.1,
    vjust = 1.2,
    size = 5,
    parse = TRUE,
    color = "#009E73",
    inherit.aes = FALSE
  ) + 
  # geom_text(
  #   data = stats_dt,
  #   aes(x = x_pos, y = y_pos, label = label),
  #   inherit.aes = FALSE,
  #   hjust = 0, vjust = 0,
  #   size = 5,
  #   parse = TRUE,
  #   color = "#009E73"
  # ) +
  facet_grid(
    model ~ indicator,
    scales = "free_x",
    labeller = labeller(
      indicator = as_labeller(indicator_labels, label_parsed)
    )
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
  labs(x = "Urban metrics", y = TeX("Difference between SLDR and gravity-/radial-based $\\rho$")) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )
#----------rho and var----------#
fig_delta_var <- fig_delta_var + plot_annotation(tag_levels = list(letters[3])) & theme(plot.tag = element_text(size = 20))
# ggsave(fig_delta_var,filename = paste(figpath,"/msa/rho_diff_var_",Yname[yindex],".pdf",sep=""), width = 4.6*3, height = 4*2)
ggsave(fig_delta_var,filename = paste(figpath,"/msa/rho_diff_var_",Yname[yindex],".pdf",sep=""), width = 5*3, height = 4*2)




#-----------------------------
# Figure: rho_emp and rho_null
#-----------------------------
xy.range <- range(c(rho_cmp$rho_emp,rho_cmp$rho_null))
fig_null <- ggplot(rho_cmp, aes(x = rho_emp, y = rho_null)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 3, stroke = 0.5) +
  facet_wrap(~ model, scales = "fixed", nrow = 1) +
  scale_fill_manual(name = "MSA",
                    breaks = msa.info$dname,
                    labels = msa.info$dname,
                    values = scales::alpha(msa.info$col, alpha = 0.5)) +
  scale_colour_manual(name = "MSA",
                      breaks = msa.info$dname,
                      labels = msa.info$dname,
                      values = msa.info$col) +
  scale_shape_manual(name = "MSA",
                     breaks = msa.info$dname,
                     labels = msa.info$dname,
                     values = msa.info$shape) +
  scale_x_continuous(limits = xy.range) +
  scale_y_continuous(limits = xy.range) +
  labs(x = TeX("Daily $\\rho_{Empirical}$"), y = TeX("Daily $\\rho_{Model}$")) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain",size=15),
        legend.position = "right")
fig_null <- fig_null + plot_annotation(tag_levels = list(letters[2])) & theme(plot.tag = element_text(size = 20))
ggsave(fig_null, filename = paste(figpath,"/msa/rho_null_",y_tag,".pdf",sep=""), width = 4.6*2, height = 4)






















#----------rho with each indicator for msa----------#
indicator_labels <- c(
  pop             = "Population~(million)",

  mean_haversine  = "Average~travel~distance~(km)",
  drive_alone = "Drive~alone",
  public_transit = "Public~transit",
  
  radial_slope    = "Radial~density~gradient",
  poly_gini       = "Polycentricity~(Gini)",
  rank_size_slope = "Rank-size~slope"
)
indicator_order <- names(indicator_labels)

rho_dt[, indicator := factor(indicator, levels = indicator_order)]
rho_dt_msa <- rho_dt[county=="msa"]

stats_r2 <- rho_dt_msa[is.finite(x) & is.finite(median_rho),
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
stats_r2[, x_pos := x_pos + 0.03 * (max(rho_dt_msa[indicator==.BY$indicator]$x, na.rm=TRUE) - 
                                      min(rho_dt_msa[indicator==.BY$indicator]$x, na.rm=TRUE)),
         by = indicator]
stats_r2[, y_pos := y_pos - 0.05 * (max(rho_dt_msa$median_rho, na.rm=TRUE) - min(rho_dt_msa$median_rho, na.rm=TRUE))]

#---- 2) Plot with R^2 annotation ----#
if ("indicator" %in% names(stats_r2)) {
  stats_r2[, indicator := factor(indicator, levels = indicator_order)]
}
fig.rho.var <- ggplot(rho_dt_msa, aes(x = x, y = median_rho)) +
  geom_point(aes(fill = msa, color = msa, shape = msa), size = 3, stroke = 0.5) +
  geom_smooth(method = "lm", color = "#69b3a2", se = FALSE) +
  facet_wrap(~ indicator, scales = "free_x", nrow = 3,
             labeller = labeller(indicator = as_labeller(indicator_labels, label_parsed))) +
  geom_text(
    data = stats_r2,
    aes(x = x_pos, y = 0.46, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
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









library(lme4)

model_delta <- lmer(
  delta_rho ~ pop_z + area_z + avg_density_z +
    mean_haversine_z + drive_alone_z + public_transit_z +
    radial_slope_z + poly_gini_z + rank_size_slope_z +
    (1 | msa),
  data = rho_all
)

model_full <- lmer(
  rho_emp ~ rho_null +
    pop_z + area_z + avg_density_z +
    mean_haversine_z + drive_alone_z + public_transit_z +
    radial_slope_z + poly_gini_z + rank_size_slope_z +
    (1 | msa),
  data = rho_all
)


dat_msa <- rho_wide[county == "msa"]
dat_msa[, day_f := factor(day)]  # day fixed effect (stable with only 7 levels)
fit_A_msa <- lmer(rho ~ pop_z + area_z + avg_density_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)
fit_B_msa <- lmer(rho ~ mean_haversine_z + drive_alone_z + public_transit_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)
fit_C_msa <- lmer(rho ~ radial_slope_z + poly_gini_z + rank_size_slope_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)
fit_best_msa <- lmer(rho ~ area_z + drive_alone_z + mean_haversine_z + radial_slope_z + factor(day) + (1 | msa),
                     data = dat_msa, REML = TRUE)
fit_best_AC_msa <- lmer(rho ~ area_z + radial_slope_z + factor(day) + (1 | msa),
                        data = dat_msa, REML = TRUE)
# R2 summary
r2A_msa      <- safe_r2(fit_A_msa)
r2B_msa      <- safe_r2(fit_B_msa)
r2C_msa      <- safe_r2(fit_C_msa)
r2best_msa   <- safe_r2(fit_best_msa)
r2bestAC_msa <- safe_r2(fit_best_AC_msa)
R2_msa <- data.table(
  scale  = "MSA",
  model  = c("A", "B", "C", "Core metrics (ABC)", "Core metrics (AC)"),
  labels = c("Urban size & density", "Mobility structure", "Urban morphology",
             "Core metrics (ABC)", "Core metrics (AC)"),
  R2m    = c(r2A_msa$R2m, r2B_msa$R2m, r2C_msa$R2m, r2best_msa$R2m, r2bestAC_msa$R2m),
  R2c    = c(r2A_msa$R2c, r2B_msa$R2c, r2C_msa$R2c, r2best_msa$R2c, r2bestAC_msa$R2c)
)


# ============================================================
# Part 3A. Final parsimonious model selection: COUNTY level
#   - All-subsets search over 6 covariates (A + B) by default
#   - If you also want C at county level AND those columns exist, set include_C_county=TRUE
#   - Constraint: not(pop & area & density) all together
# ============================================================
options(na.action = "na.fail")  # required for dredge

include_C_county <- FALSE  # set TRUE only if county dataset contains C variables

vars_A <- c("pop_z","area_z","avg_density_z")
vars_B <- c("mean_haversine_z","drive_alone_z","public_transit_z")
vars_C <- c("radial_slope_z","poly_gini_z","rank_size_slope_z")

vars_county <- if (include_C_county) c(vars_A, vars_B, vars_C) else c(vars_A, vars_B)

form_full_cty <- as.formula(
  paste0("rho ~ ", paste(vars_county, collapse = " + "), " + (1|msa)")
)

fit_full_cty <- lmer(form_full_cty, data = dat_county, REML = FALSE)

dd_cty <- dredge(
  fit_full_cty,
  rank   = "AICc",
  trace  = FALSE,
  subset = !(pop_z & area_z & avg_density_z)   # mild collinearity constraint
)

fit_best_cty <- get.models(dd_cty, subset = 1)[[1]]
cand_cty     <- get.models(dd_cty, subset = delta <= 2)

# Print key outputs
summarize_selection(dd_cty, fit_best_cty, cand_cty, tag = "County: parsimonious selection")

# ============================================================
# Part 3B. Final parsimonious model selection: MSA level (daily)
#   - All-subsets search over 9 covariates (A + B + C)
#   - Always include day_f + (1|msa)
#   - Constraint: not(pop & area & density) all together
# ============================================================
vars_msa <- c(vars_A, vars_B, vars_C)

form_full_msa <- as.formula(
  paste0("rho ~ ", paste(vars_msa, collapse = " + "), " + factor(day) + (1|msa)")
)

fit_full_msa <- lmer(form_full_msa, data = dat_msa, REML = FALSE)

dd_msa <- dredge(
  fit_full_msa,
  rank   = "AICc",
  trace  = FALSE,
  subset = !(pop_z & area_z & avg_density_z)
)

fit_best_msa <- get.models(dd_msa, subset = 1)[[1]]
cand_msa     <- get.models(dd_msa, subset = delta <= 2)

summarize_selection(dd_msa, fit_best_msa, cand_msa, tag = "MSA: parsimonious selection")

# Restore NA action if desired
options(na.action = "na.omit")

# ============================================================
# Optional: side-by-side comparison tables (best county vs best msa)
# ============================================================
best_compare <- rbind(
  data.table(
    scale = "County",
    AIC   = AIC(fit_best_cty),
    BIC   = BIC(fit_best_cty),
    R2m   = safe_r2(update(fit_best_cty, REML = TRUE))$R2m,
    R2c   = safe_r2(update(fit_best_cty, REML = TRUE))$R2c,
    formula = deparse(formula(fit_best_cty))
  ),
  data.table(
    scale = "MSA",
    AIC   = AIC(fit_best_msa),
    BIC   = BIC(fit_best_msa),
    R2m   = safe_r2(update(fit_best_msa, REML = TRUE))$R2m,
    R2c   = safe_r2(update(fit_best_msa, REML = TRUE))$R2c,
    formula = deparse(formula(fit_best_msa))
  )
)
best_compare




#===========================
# Combine + plot 
#===========================
R2_all <- rbind(R2_county, R2_msa, fill = TRUE)
R2_all[, labels := factor(labels, levels = c("Core metrics (AC)", "Core metrics (ABC)", 
                                             "A+B", "Urban morphology", "Mobility structure", "Urban size & density"))]
geovar.info <- data.frame(breaks = c("Urban size & density", "Mobility structure", "Urban morphology","A+B",
                                     "Core metrics (ABC)", "Core metrics (AC)"),
                          labels = c("A: Urban size & density", "B: Mobility structure", "C: Urban morphology","A+B",
                                     "Core metrics (ABC)", "Core metrics (AC)"),
                          group = c("A", "B", "C", "A+B","core", "coreAC"),
                          color = col.all[c(15,16,18,17,14,17)])
fig.r2.decomp <- ggplot(R2_all[scale=="MSA" & model %in% c("Core metrics (ABC)", "Core metrics (AC)")], aes(x = R2m, y = labels, fill = labels)) +
  geom_col(width = 0.7) +
  facet_wrap(~ scale, ncol = 1, scales = "free_y") +
  geom_text(aes(x = 0.02, label = sprintf("R^2 == %.3f", R2m)), parse = TRUE,
            color = "black", size = 5, hjust = 0) +
  scale_fill_manual(name="Urban metrics",
                    breaks = geovar.info$breaks,
                    labels = geovar.info$labels,
                    values = scales::alpha(geovar.info$col, alpha = 0.8))+
  scale_colour_manual(name = "Urban metrics",
                      breaks = geovar.info$breaks,
                      labels = geovar.info$labels,
                      values = geovar.info$col) +
  scale_x_continuous(limits = c(0,0.70), expand = 0) +
  labs(x = TeX('$\\R^2$'), y = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")

# plot variables
vars_A <- c("Population", "Urban area (+)", "Population density")
vars_B <- c("Average travel distance (-)", "Drive alone (+)", "Public transit")
vars_C <- c("Radial density gradient (-)", "Polycentricity (Gini)", "Rank-size slope")  
df_tbl <- rbindlist(list(
  data.table(group = "A", var = vars_A),
  data.table(group = "B", var = vars_B),
  data.table(group = "C", var = vars_C)
))
df_tbl[, y := rev(seq_len(.N)), by = group]

core_vars_county <- c("Population", "Average travel distance (-)", "Drive alone (+)", "Public transit")
core_vars_msa <- c("Urban area (+)", "Average travel distance (-)", "Drive alone (+)", "Radial density gradient (-)")
core_county <- copy(df_tbl)[, `:=`(
  scale = "County",
  selected = var %in% core_vars_county
)]
core_msa <- copy(df_tbl)[, `:=`(
  scale = "MSA",
  selected = var %in% core_vars_msa
)]
core_all <- rbind(core_county, core_msa)[scale=="MSA"]
core_all_bg <- merge(core_all, geovar.info, by = "group")
header_info <- copy(geovar.info[-(4:6), ]) %>%as.data.table
header_info[, y := max(core_all$y) + 1]
fig.r2.core <- ggplot(core_all, aes(x = group, y = y)) +
  geom_tile(data = core_all_bg, aes(fill = color),
            alpha = 0.6, width = 0.9, height = 0.9, color = "grey40") +
  geom_text(data = core_all[selected != TRUE], aes(label = var), size = 5) +
  geom_text(data = core_all[selected == TRUE], aes(label = var), size = 5, color = "#ac6894") +
  geom_text(data = header_info[group%in%c("A","C")], color = "#83ab8e",
            aes(x = group, y = y, label = paste0(group, "\n", breaks)),
            inherit.aes = FALSE, fontface = "bold", size = 5) +
  geom_text(data = header_info[group=="B"], aes(x = group, y = y, label = paste0(group, "\n", breaks)),
            inherit.aes = FALSE, fontface = "bold", size = 5) +
  # geom_text(data = core_all[selected == TRUE],
  #           aes(label = "\u2713"), color = "#ac6894",
  #           fontface = "bold", size = 7, vjust = -1.2) +
  # geom_text(data = core_all[selected == TRUE], aes(label = "\u2713"),
  #           color = "#ac6894", fontface = "bold", size = 7, vjust = 0, 
  #           position = position_nudge(x = 0.35)) +
  facet_wrap(~scale, nrow = 2) +
  labs(x = NULL, y = NULL) +
  scale_fill_identity() +
  scale_x_discrete(breaks = NULL, labels = NULL) +
  scale_y_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0.1, 0.2))) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")

#----------Model performance for Atlanta and Boston----------#
fig.rho.var.r2 <- (fig.r2.decomp | fig.r2.core) + plot_layout(widths = c(1, 2)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
fig.rho.var.r2 <- (fig.r2.decomp | fig.r2.core) + plot_layout(widths = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20), legend.position = "bottom", legend.justification = c(0.5, 0))
ggsave(fig.rho.var.r2, filename = paste(figpath,"/msa/SI_rho_var_r2_",Yname[yindex],".pdf",sep=""), width =4*4, height = 5.5*1)


