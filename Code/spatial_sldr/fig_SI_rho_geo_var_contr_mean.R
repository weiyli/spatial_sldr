# rm(list = ls())

# Mean-CBG-area version
# [ rho vs travel distance ]   [ rho vs travel mode ]   [ rho vs area ]
# county level data from urban_sldr_fit.R; msa level data from sldr_fit.R; variable from rho_geo_var.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(sf)          # read_sf() 
library(spdep)       # poly2nb() https://blog.csdn.net/weixin_54000907/article/details/116247097
library(gridExtra)
library(ggtext)
# library(relaimpo)    # 
library(lme4)    # lmer()
library(lmerTest) 
library(MuMIn)   # r.squaredGLMM()



#----------Error-measure: R2----------#
R2_test <- function(y_test,y_predict){
  SS.tot = sum((y_test-mean(y_test))^2)
  SS.err = sum((y_predict-y_test)^2)
  result <- 1-SS.err/SS.tot
  return(result)
}

get.significance <- function(p) {
  if (p < 0.001) {
    return("p<0.001")
  } else if (p < 0.01) {
    return("p<0.01")
  } else if (p < 0.05) {
    return("p<0.05")
  } else {
    return("")
  }
}


# ---------- Helpers ----------
get_R2m <- function(r2_mat) as.numeric(r2_mat[1, "R2m"])
get_R2c <- function(r2_mat) as.numeric(r2_mat[1, "R2c"])

safe_r2 <- function(fit) {
  r2 <- r.squaredGLMM(fit)
  list(R2m = get_R2m(r2), R2c = get_R2c(r2), raw = r2)
}

# A small utility to print core outputs in a consistent way
summarize_selection <- function(dd, best_fit, cand_fits, tag = "MODEL") {
  cat("\n====================", tag, "====================\n")
  cat("Best model (AICc):\n")
  print(formula(best_fit))
  cat("\nSummary(best):\n")
  print(summary(best_fit))
  cat("\nR2 (Nakagawa):\n")
  print(safe_r2(update(best_fit, REML = TRUE)))  # R2 should be reported on REML fit
  cat("\nVariable importance (sum of Akaike weights) within 螖AICc<=2:\n")
  print(sw(cand_fits))
  cat("\nTop 10 models by AICc:\n")
  print(as.data.table(dd)[1:min(10, nrow(dd))])
  invisible(NULL)
}

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
id.week<-paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")

#----------Load msa-county / msa-tract relationship tables----------#
county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

#----------rho and var----------#
geo_var <- fread(paste0(datapath,"/msa/geo_var.csv",sep=""))
geo_var_county <- fread(paste0(datapath,"/county/geo_var.csv",sep=""))

#----------pop and area for county and msa----------#
queen.lag.msa <- fread(file.path(datapath, "msa", "queen_lag.csv"))
queen.lag.msa[, pop  := pop  / 1e6]  # million
queen.lag.msa[, area := area / 1e6]  # km^2
queen.max.msa <- unique(queen.lag.msa[, .(msa, nblock, lag.max, pop, area)])

#----------pop and area for county and msa----------#
queen.lag.county <- fread(file.path(datapath, "county", "queen_lag.csv"))
queen.lag.county[, pop  := pop  / 1e6]  # million
queen.lag.county[, area := area / 1e6]  # km^2
queen.max.county <- unique(queen.lag.county[, .(msa, nblock, lag.max, pop, area,county)])

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
rho_var[, mean_cbg_area := mean_cbg_area_km2]
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
rho_var[, (vars_norm) := lapply(.SD, function(x) x/pop_total), .SDcols = vars_norm]

rho_var[, c("pop_total", "area_km2") := NULL]
cols_keep <- c(
  "msa",
  "day",
  "county",
  "rho",
  "pop",
  "area",
  "mean_cbg_area",
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
rho_dt_msa <- as.data.table(rho_long)


#----------county level: data from sldr_fit_county.R and od_dist.R----------#
rho.county <- fread(paste0(datapath, "/county/SLDR_params_", Yname[yindex], ".csv", sep=""))
rho.county <- rho.county[index == "exp" & period == pervalue[1]]
rho.county <- merge(rho.county, queen.max.county, by = c("msa", "county"))
# daily panel for county: (msa, county, day)
var.county <- rho.county[, .(msa, day, rho, pop, area, county)]
#----------distance and destination (daily, county-level)----------#
dist.county.all <- NULL
for (d in 1:Nmsa) {
  Dname <- Dname.set[d]
  
  county.msa <- subset(county.info, msa_id == d)
  if (nrow(county.msa) == 0) next
  
  for (f in 1:nrow(county.msa)) {
    
    # keep county id as 5-digit string to match your tables
    county.id <- as.numeric(county.msa$county_fips[f])
    # days for this (msa, county)
    days_cty <- unique(var.county[msa == Dname & county == county.id]$day)
    if (length(days_cty) == 0) next
    
    for (dd in 1:length(days_cty)) {
      
      fpath <- paste0(flowpath, "/county/", Dname, "/dist_", days_cty[dd], ".csv")
      if (!file.exists(fpath)) next
      dist <- fread(fpath)
      # ensure county column name is consistent
      if ("county_fips" %in% names(dist)) setnames(dist, "county_fips", "county")
      # enforce same county id format for matching
      dist[, county := as.numeric(county)]
      # filter to the target county only (important: file contains all counties for that MSA/day)
      dist <- dist[county == county.id]
      if (nrow(dist) == 0) next
      dist[, msa := Dname]
      dist[, day := days_cty[dd]]
      dist.county.all <- plyr::rbind.fill(dist.county.all, dist)
    }
  }
}

dist.county.all <- as.data.table(dist.county.all)
dist.county.all[, weighted_haversine := weighted_haversine / 1000]  # m -> km
# aggregate and join by (msa, county, day) 
dist.county <- dist.county.all[, .(
  mean_haversine    = mean(weighted_haversine, na.rm = TRUE),
  mean_destinations = mean(n_destinations, na.rm = TRUE)
), by = .(msa, county, day)]

rho.dist.county <- left_join(var.county, dist.county, by = c("msa", "county", "day"))


#----------Connecting rho with urban indicators----------#
rho_var_county <- left_join(rho.dist.county, geo_var_county, by=c("msa","county"))
rho_var_county[, mean_cbg_area := mean_cbg_area_km2]
# travel mode
travel.mode.county <- read.csv(paste(geopath,"/census/acs_2019_5years_travel_mode_county.csv",sep=""), header=TRUE)
setnames(travel.mode.county, "county_fips", "county")
rho_var_county <- left_join(rho_var_county, travel.mode.county, by=c("county"))
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
rho_var_county[, (vars_norm) := lapply(.SD, function(x) x/pop_total), .SDcols = vars_norm]
rho_var_county[, c("pop_total", "area_km2") := NULL]
cols_keep <- c(
  "msa",
  "day",
  "county",
  "rho",
  "pop",
  "area",
  "mean_cbg_area",
  "mean_haversine",
  "compactness",
  
  "carpool",
  "drive_alone",
  "work_at_home",
  "walk",
  "public_transit",
  "bicycle"
)
rho_sub <- rho_var_county[, ..cols_keep]
rho_long <- melt(
  rho_sub,
  id.vars = c("msa", "county", "day", "rho"),
  variable.name = "indicator",
  value.name = "x"
)
rho_dt_county <- as.data.table(rho_long)
rho_dt_county[, county := as.character(county)]
rho_dt_county[, county := "county"]

#----------row bind msa and county----------#
rho_dt <- rbind(rho_dt_msa, rho_dt_county, fill = TRUE)
rho_wide <- dcast(
  rho_dt,
  msa + county + day + rho ~ indicator,
  value.var = "x"
)



#----------Analysis of variable contributions----------#
# (1) Structural / geometric controls
vars_A <- c("pop", "area", "mean_cbg_area")
# (2) Mobility / interaction variables
vars_B <- c("mean_haversine","drive_alone", "public_transit")
# (3) Metropolitan-scale morphology (internal organization; MSA-only)
vars_C  <- c("radial_slope", "poly_gini", "rank_size_slope") 
# all
vars <- unique(c(vars_A, vars_B, vars_C))
# scale var
safe_scale <- function(x) {
  if (all(!is.finite(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}
to_scale <- vars[vars %in% names(rho_wide)]
rho_wide[, paste0(to_scale, "_z") := lapply(.SD, safe_scale), .SDcols = to_scale]




# ============================================================
# Part 1. Block-wise contributions and AICc core-model selection
# ============================================================
dat_county <- rho_wide[county == "county"]
dat_msa <- rho_wide[county == "msa"]

vars_A_z <- c("pop_z", "area_z", "mean_cbg_area_z")
vars_B_z <- c("mean_haversine_z", "drive_alone_z", "public_transit_z")
vars_C_z <- c("radial_slope_z", "poly_gini_z", "rank_size_slope_z")

sig_p_cut <- 0.01
refit_significant_terms <- function(fit, data, p_cut = 0.01) {
  coef_dt <- as.data.table(coef(summary(fit)), keep.rownames = "term")
  p_col <- "Pr(>|t|)"
  selected_terms <- attr(terms(fit), "term.labels")
  candidate_terms <- setdiff(selected_terms, "factor(day)")
  sig_terms <- coef_dt[
    term %in% candidate_terms & get(p_col) < p_cut,
    term
  ]
  sig_terms <- candidate_terms[candidate_terms %in% sig_terms]
  
  rhs_terms <- c(sig_terms, "factor(day)", "(1 | msa)")
  form_sig <- as.formula(paste("rho ~", paste(rhs_terms, collapse = " + ")))
  fit_sig <- lmer(form_sig, data = data, REML = FALSE)
  attr(fit_sig, "sig_table") <- coef_dt[
    term %in% candidate_terms,
    .(term, p_value = get(p_col), retained = term %in% sig_terms)
  ]
  fit_sig
}

#----------County: block-wise contribution----------#
fit_A_county <- lmer(rho ~ pop_z + area_z + mean_cbg_area_z + factor(day) + (1 | msa),
                     data = dat_county, REML = TRUE)
fit_B_county <- lmer(rho ~ mean_haversine_z + drive_alone_z + public_transit_z + factor(day) + (1 | msa),
                     data = dat_county, REML = TRUE)
fit_AB_county <- lmer(rho ~ pop_z + area_z + mean_cbg_area_z +
                        mean_haversine_z + drive_alone_z + public_transit_z +
                        factor(day) + (1 | msa),
                      data = dat_county, REML = TRUE)

#----------MSA: block-wise contribution----------#
fit_A_msa <- lmer(rho ~ pop_z + area_z + mean_cbg_area_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)
fit_B_msa <- lmer(rho ~ mean_haversine_z + drive_alone_z + public_transit_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)
fit_C_msa <- lmer(rho ~ radial_slope_z + poly_gini_z + rank_size_slope_z + factor(day) + (1 | msa),
                  data = dat_msa, REML = TRUE)

#----------AICc selection: county core model from A+B candidates----------#
options(na.action = "na.fail")

include_C_county <- FALSE
vars_county <- if (include_C_county) c(vars_A_z, vars_B_z, vars_C_z) else c(vars_A_z, vars_B_z)

form_full_cty <- as.formula(
  paste0("rho ~ ", paste(vars_county, collapse = " + "), " + factor(day) + (1 | msa)")
)
fit_full_cty <- lmer(form_full_cty, data = dat_county, REML = FALSE)
dd_cty <- dredge(
  fit_full_cty,
  rank = "AICc",
  trace = FALSE,
  subset = !(pop_z & area_z & mean_cbg_area_z)
)
fit_best_cty <- get.models(dd_cty, subset = 1)[[1]]
cand_cty <- get.models(dd_cty, subset = delta <= 2)
summarize_selection(dd_cty, fit_best_cty, cand_cty, tag = "County: parsimonious selection")
fit_sig_cty <- refit_significant_terms(fit_best_cty, dat_county, sig_p_cut)

#----------AICc selection: MSA Core metrics (AC)----------#
vars_AC_msa <- c(vars_A_z, vars_C_z)
form_full_AC_msa <- as.formula(
  paste0("rho ~ ", paste(vars_AC_msa, collapse = " + "), " + factor(day) + (1 | msa)")
)
fit_full_AC_msa <- lmer(form_full_AC_msa, data = dat_msa, REML = FALSE)
dd_AC_msa <- dredge(
  fit_full_AC_msa,
  rank = "AICc",
  trace = FALSE,
  subset = !(pop_z & area_z & mean_cbg_area_z)
)
fit_best_AC_msa <- get.models(dd_AC_msa, subset = 1)[[1]]
cand_AC_msa <- get.models(dd_AC_msa, subset = delta <= 2)
summarize_selection(dd_AC_msa, fit_best_AC_msa, cand_AC_msa, tag = "MSA: Core metrics (AC)")
fit_sig_AC_msa <- refit_significant_terms(fit_best_AC_msa, dat_msa, sig_p_cut)

#----------AICc selection: MSA Core metrics (ABC)----------#
vars_ABC_msa <- c(vars_A_z, vars_B_z, vars_C_z)
form_full_ABC_msa <- as.formula(
  paste0("rho ~ ", paste(vars_ABC_msa, collapse = " + "), " + factor(day) + (1 | msa)")
)
fit_full_ABC_msa <- lmer(form_full_ABC_msa, data = dat_msa, REML = FALSE)
dd_ABC_msa <- dredge(
  fit_full_ABC_msa,
  rank = "AICc",
  trace = FALSE,
  subset = !(pop_z & area_z & mean_cbg_area_z)
)
fit_best_ABC_msa <- get.models(dd_ABC_msa, subset = 1)[[1]]
cand_ABC_msa <- get.models(dd_ABC_msa, subset = delta <= 2)
summarize_selection(dd_ABC_msa, fit_best_ABC_msa, cand_ABC_msa, tag = "MSA: Core metrics (ABC)")
fit_sig_ABC_msa <- refit_significant_terms(fit_best_ABC_msa, dat_msa, sig_p_cut)

sig_selection <- rbindlist(list(
  cbind(model = "County: core indicators", attr(fit_sig_cty, "sig_table")),
  cbind(model = "MSA: Core metrics (AC)", attr(fit_sig_AC_msa, "sig_table")),
  cbind(model = "MSA: Core metrics (ABC)", attr(fit_sig_ABC_msa, "sig_table"))
), fill = TRUE)
sig_selection

options(na.action = "na.omit")

# Refit AICc + significance-retained core models with REML for R2 summaries.
fit_best_county <- update(fit_sig_cty, REML = TRUE)
fit_best_msa <- update(fit_sig_ABC_msa, REML = TRUE)
fit_best_AC_msa_reml <- update(fit_sig_AC_msa, REML = TRUE)

#----------R2 summary----------#
r2A_cty <- safe_r2(fit_A_county)
r2B_cty <- safe_r2(fit_B_county)
r2AB_cty <- safe_r2(fit_AB_county)
r2best_cty <- safe_r2(fit_best_county)

R2_county <- data.table(
  scale = "County",
  model = c("A", "B", "A+B", "Core indicators"),
  labels = c("Urban size & CBG area", "Mobility structure", "A+B", "Core indicators"),
  R2m = c(r2A_cty$R2m, r2B_cty$R2m, r2AB_cty$R2m, r2best_cty$R2m),
  R2c = c(r2A_cty$R2c, r2B_cty$R2c, r2AB_cty$R2c, r2best_cty$R2c)
)
R2_county

r2A_msa <- safe_r2(fit_A_msa)
r2B_msa <- safe_r2(fit_B_msa)
r2C_msa <- safe_r2(fit_C_msa)
r2best_msa <- safe_r2(fit_best_msa)
r2bestAC_msa <- safe_r2(fit_best_AC_msa_reml)

R2_msa <- data.table(
  scale = "MSA",
  model = c("A", "B", "C", "Core metrics (ABC)", "Core metrics (AC)"),
  labels = c("Urban size & CBG area", "Mobility structure", "Urban morphology",
             "Core metrics (ABC)", "Core metrics (AC)"),
  R2m = c(r2A_msa$R2m, r2B_msa$R2m, r2C_msa$R2m, r2best_msa$R2m, r2bestAC_msa$R2m),
  R2c = c(r2A_msa$R2c, r2B_msa$R2c, r2C_msa$R2c, r2best_msa$R2c, r2bestAC_msa$R2c)
)
R2_msa

best_compare <- rbind(
  data.table(
    scale = "County",
    model = "Core indicators",
    AIC = AIC(fit_sig_cty),
    BIC = BIC(fit_sig_cty),
    R2m = r2best_cty$R2m,
    R2c = r2best_cty$R2c,
    formula = paste(deparse(formula(fit_sig_cty)), collapse = " ")
  ),
  data.table(
    scale = "MSA",
    model = "Core metrics (AC)",
    AIC = AIC(fit_sig_AC_msa),
    BIC = BIC(fit_sig_AC_msa),
    R2m = r2bestAC_msa$R2m,
    R2c = r2bestAC_msa$R2c,
    formula = paste(deparse(formula(fit_sig_AC_msa)), collapse = " ")
  ),
  data.table(
    scale = "MSA",
    model = "Core metrics (ABC)",
    AIC = AIC(fit_sig_ABC_msa),
    BIC = BIC(fit_sig_ABC_msa),
    R2m = r2best_msa$R2m,
    R2c = r2best_msa$R2c,
    formula = paste(deparse(formula(fit_sig_ABC_msa)), collapse = " ")
  )
)
best_compare




#===========================
# Combine + plot 
#===========================
R2_all <- rbind(R2_county, R2_msa, fill = TRUE)
R2_all[, labels := factor(labels, levels = c("Core metrics (AC)", "Core metrics (ABC)", 
                                             "A+B", "Urban morphology", "Mobility structure", "Urban size & CBG area"))]
geovar.info <- data.frame(breaks = c("Urban size & CBG area", "Mobility structure", "Urban morphology","A+B",
                                     "Core metrics (ABC)", "Core metrics (AC)"),
                          labels = c("A: Urban size & CBG area", "B: Mobility structure", "C: Urban morphology","A+B",
                                     "Core metrics (ABC)", "Core metrics (AC)"),
                          group = c("A", "B", "C", "A+B","core", "coreAC"),
                          color = col.all[c(15,16,18,17,14,17)])
fig.r2.decomp <- ggplot(R2_all[scale=="MSA" & model %in% c("Core metrics (ABC)", "Core metrics (AC)")], aes(x = R2m, y = labels, fill = labels)) +
  geom_col(width = 0.7) +
  facet_wrap(~ scale, ncol = 1, scales = "free_y") +
  geom_text(aes(x = 0.1, label = sprintf("R^2 == %.3f", R2m)), parse = TRUE,
            color = "black", size = 5, hjust = 0) +
  scale_fill_manual(name="Urban metrics",
                    breaks = geovar.info$breaks,
                    labels = geovar.info$labels,
                    values = scales::alpha(geovar.info$col, alpha = 0.8))+
  scale_colour_manual(name = "Urban metrics",
                      breaks = geovar.info$breaks,
                      labels = geovar.info$labels,
                      values = geovar.info$col) +
  scale_x_continuous(limits = c(0,0.80), expand = 0) +
  labs(x = TeX('$\\R^2$'), y = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")

# plot variables
vars_A <- c("Population (+)", "Urban area", "CBG area (+)")
vars_B <- c("Average travel distance", "Drive alone (+)", "Public transit")
vars_C <- c("Radial density gradient (-)", "Polycentricity (Gini)", "Rank-size slope")
df_tbl <- rbindlist(list(
  data.table(group = "A", var = vars_A),
  data.table(group = "B", var = vars_B),
  data.table(group = "C", var = vars_C)
))
df_tbl[, y := rev(seq_len(.N)), by = group]

term_labels <- c(
  pop_z = "Population (+)",
  area_z = "Urban area",
  mean_cbg_area_z = "Mean CBG area (+)",
  mean_haversine_z = "Average travel distance",
  drive_alone_z = "Drive alone (+)",
  public_transit_z = "Public transit",
  radial_slope_z = "Radial density gradient (-)",
  poly_gini_z = "Polycentricity (Gini)",
  rank_size_slope_z = "Rank-size slope"
)
selected_term_labels <- function(fit) {
  selected_terms <- attr(terms(fit), "term.labels")
  selected_terms <- setdiff(selected_terms, "factor(day)")
  selected_labels <- unname(term_labels[selected_terms])
  selected_labels[!is.na(selected_labels)]
}
core_vars_county <- selected_term_labels(fit_best_county)
core_vars_msa <- selected_term_labels(fit_best_msa)
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
# fig.rho.var.r2 <- (fig.r2.decomp | fig.r2.core) + plot_layout(widths = c(1, 2)) +
#   plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
# fig.rho.var.r2 <- (fig.r2.decomp | fig.r2.core) + plot_layout(widths = c(1, 2), guides = "collect") +
#   plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20), legend.position = "bottom", legend.justification = c(0.5, 0))
fig.rho.var.r2 <- (fig.r2.decomp | fig.r2.core) + plot_layout(widths = c(0.7, 2), guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20), legend.position = "right", legend.justification = c(0, 0.5))
ggsave(fig.rho.var.r2, filename = paste(figpath,"/msa/SI_rho_var_r2_",Yname[yindex],"_mean.pdf",sep=""), width =4*4, height = 4.5*1)












