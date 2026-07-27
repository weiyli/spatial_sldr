# rm(list = ls())

#-----------------------------#
# Build Table Data for Kernel Uncertainty (MBB over time)
# Outputs: City, Period, Kernel, T(days), N(CBGs), Mean rho, 95% CI(rho), CI by Moving Block Bootstrap (MBB)
#-----------------------------#

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath  <- 'D:/ood/Data/Geo'
datapath <- 'D:/ood/Data/spatial_sldr'

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
Dname <<- Dname.set[d]
date.gif <- date.sep(Datevalue,durdate)

#----------MBB functions----------#
mbb_resample_index <- function(n, L) {
  stopifnot(n >= L)
  starts <- 1:(n - L + 1)
  n_blocks <- ceiling(n / L)
  s <- sample(starts, n_blocks, replace = TRUE)
  idx <- unlist(lapply(s, function(a) a:(a + L - 1)))
  idx[1:n]
}

mbb_ci <- function(x, stat_fun, L = 7, B = 2000, seed = 1) {
  # x must be time-ordered
  set.seed(seed)
  n <- length(x)
  if (n < L) stop("Length(x) must be >= L.")
  boots <- replicate(B, {
    idx <- mbb_resample_index(n, L)
    stat_fun(x[idx])
  })
  ci <- stats::quantile(boots, probs = c(0.025, 0.5, 0.975), na.rm = TRUE, names = FALSE)
  list(mean = mean(boots, na.rm = TRUE),
       sd   = stats::sd(boots, na.rm = TRUE),
       ci   = ci)
}

#----------User settings----------#
yindex <- 2
B_boot <- 2000
L_block <- 2
seed0 <- 1
kernel_keep <- c("homo","power","exp" )

#----------Read daily rho----------#
rho.daily <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], "_half_year.csv"))
rho.daily <- rho.daily[index %in% kernel_keep]
rho.daily <- rho.daily[period %in% c("before", "during")]
# Ensure day is Date and ordered within groups
rho.daily[, day := as.Date(day)]
setorder(rho.daily, msa, period, index, day)

#----------Kernel naming for table----------#
kernel_map <- c(
  homo = "No-decay",
  power = "Power-law",
  exp  = "Exponential"
)
rho.daily[, kernel := kernel_map[index]]

#----------Compute N(CBGs) per MSA----------#
get_n_cbg <- function(msa_name, region_name) {
  geo_file <- file.path(geopath, region_name, paste0(msa_name, ".geojson"))
  if (!file.exists(geo_file)) {
    return(NA_integer_)
  }
  shp <- suppressWarnings(sf::read_sf(geo_file))
  if (!("CensusBlockGroup" %in% names(shp))) return(NA_integer_)
  length(unique(shp$CensusBlockGroup))
}

ncbg_dt <- rbindlist(lapply(Dname.msa, function(m) {
  data.table(msa = m, n_cbg = get_n_cbg(m, "msa"))
}), fill = TRUE)

#----------Compute MBB CI for each (msa, period, kernel)----------#
stat_fun <- function(z) mean(z, na.rm = TRUE)
tab_dt <- rho.daily[, {
  x <- rho
  # guard: if too few days, fallback to NA
  T_days <- length(x)
  if (T_days < L_block) {
    list(
      T_days   = T_days,
      mean_rho = mean(x, na.rm = TRUE),
      ci_low   = NA_real_,
      ci_med   = NA_real_,
      ci_high  = NA_real_
    )
  } else {
    mbb <- mbb_ci(x, stat_fun = stat_fun, L = L_block, B = B_boot, seed = seed0)
    list(
      T_days   = T_days,
      mean_rho = mbb$mean,
      ci_low   = mbb$ci[1],
      ci_med   = mbb$ci[2],
      ci_high  = mbb$ci[3]
    )
  }
}, by = .(msa, period, index, kernel)]

#----------Merge N(CBGs) and format columns----------#
tab_dt <- merge(tab_dt, ncbg_dt[, .(msa, n_cbg)], by = "msa", all.x = TRUE)
tab_dt[, period_label := fifelse(period == "before", "Before", "During")]
kernel_order <- c("No-decay", "Power-law","Exponential")
tab_dt[, kernel := factor(kernel, levels = kernel_order)]
setorder(tab_dt, msa, period_label, kernel)
tab_out <- tab_dt[, .(
  City      = msa,
  Period    = period_label,
  Kernel    = as.character(kernel),
  `T_days`  = T_days,
  `N_CBGs`  = n_cbg,
  `Mean_rho` = round(mean_rho, 4),
  `CI_rho`   = sprintf("[%.4f, %.4f]", ci_low, ci_high)
)]
out_csv <- file.path(datapath, "msa", paste0("rho_CI_table_", Yname[yindex], ".csv"))
fwrite(tab_out, out_csv)



