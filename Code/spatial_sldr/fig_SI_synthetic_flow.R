# rm(list = ls())

# Data from synthetic_flow.R
# Compare null rho (gravity and radiation models) with empirical rho


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

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
date.gif<-date.sep(Datevalue,durdate)


#--------------------------------------#
#  Null rho and Empirical rho  #
#--------------------------------------#
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
  labs(x = TeX('Daily $\\rho_{Empirical}$'), y = TeX('Daily $\\rho_{Null}$')) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain",size=15),
        legend.position = "right")
ggsave(fig_null, filename = paste(figpath,"/msa/SI_rho_null_",y_tag,".pdf",sep=""), width = 4.6*2, height = 4)


#-----------------------------
# Table
#-----------------------------
# paired test pooled across all city-days
rho_cmp[model == "Gravity", wilcox.test(rho_emp, rho_null, paired = TRUE)]
rho_cmp[model == "Radial", wilcox.test(rho_emp, rho_null, paired = TRUE)]
# 
city_rank <- rho_cmp[, .( rho_emp_med  = median(rho_emp, na.rm=TRUE), rho_null_med = median(rho_null, na.rm=TRUE)),by = .(model, msa)]
city_rank[, cor_spear := cor(rho_emp_med, rho_null_med, method="spearman"), by = model]
unique(city_rank[, .(model, cor_spear)])


#-----------------------------
# Helpers
#-----------------------------
fmt_p <- function(p) {
  out <- character(length(p))
  out[is.na(p)] <- NA_character_
  out[!is.na(p) & p < 2.2e-16] <- "<2.2e-16"
  out[!is.na(p) & p >= 2.2e-16] <- format.pval(p[!is.na(p) & p >= 2.2e-16],
                                               digits = 3, eps = 1e-16)
  out
}

# Paired sign-flip permutation test (distribution-free; exchangeability within pairs)
perm_signflip <- function(d, B = 5000, seed = 1, stat = c("median", "mean")) {
  stat <- match.arg(stat)
  set.seed(seed)
  d <- d[is.finite(d)]
  if (length(d) == 0) return(list(obs = NA_real_, p = NA_real_))
  obs <- if (stat == "median") median(d) else mean(d)
  
  flips <- replicate(B, {
    s <- sample(c(-1, 1), length(d), replace = TRUE)
    dp <- d * s
    if (stat == "median") median(dp) else mean(dp)
  })
  p <- mean(abs(flips) >= abs(obs))
  list(obs = obs, p = p)
}

# Two-sided paired sign test (very low assumptions; uses only direction)
sign_test <- function(d) {
  d <- d[is.finite(d)]
  d <- d[d != 0]
  n <- length(d)
  if (n == 0) return(list(n = 0, p = NA_real_))
  k <- sum(d > 0)
  # exact binomial test under H0: P(d>0)=0.5
  p <- binom.test(k, n, p = 0.5, alternative = "two.sided")$p.value
  list(n = n, p = p)
}

#-----------------------------
# Main: build test table per model
#-----------------------------
build_test_table <- function(rho_cmp, B_perm = 5000, seed_perm = 1,
                             cov_lo = 0.05, cov_hi = 0.95) {
  
  dt <- copy(as.data.table(rho_cmp))
  dt <- dt[is.finite(rho_emp) & is.finite(rho_null)]
  
  # city-level medians for Spearman rank (emp vs null)
  city_rank <- dt[
    , .(
      rho_emp_city = median(rho_emp, na.rm = TRUE),
      rho_null_city = median(rho_null, na.rm = TRUE)
    ),
    by = .(model, msa)
  ]
  
  out <- dt[
    ,
    {
      d <- rho_emp - rho_null
      
      # Effect summaries (assumption-light)
      n_pairs <- .N
      med_delta <- median(d, na.rm = TRUE)
      mean_delta <- mean(d, na.rm = TRUE)
      iqr_delta <- IQR(d, na.rm = TRUE)
      prop_pos <- mean(d > 0, na.rm = TRUE)
      
      # Wilcoxon signed-rank (paired; nonparametric but uses ranks + symmetry-ish)
      p_wil <- tryCatch(wilcox.test(rho_emp, rho_null, paired = TRUE)$p.value, error = function(e) NA_real_)
      
      # Sign test (paired; only direction)
      st <- sign_test(d)
      p_sign <- st$p
      n_sign <- st$n
      
      # Permutation sign-flip (paired; distribution-free given exchangeability)
      perm <- perm_signflip(d, B = B_perm, seed = seed_perm, stat = "median")
      p_perm <- perm$p
      
      # KS test on pooled distributions (nonparametric; ignores pairing but shows distribution mismatch)
      p_ks <- tryCatch(ks.test(rho_emp, rho_null)$p.value, error = function(e) NA_real_)
      
      # Coverage: how often empirical falls outside central envelope of null (assumption-light)
      qlo <- quantile(rho_null, cov_lo, na.rm = TRUE)
      qhi <- quantile(rho_null, cov_hi, na.rm = TRUE)
      prop_emp_outside <- mean(rho_emp < qlo | rho_emp > qhi, na.rm = TRUE)
      
      list(
        n_pairs = n_pairs,
        med_delta = med_delta,
        mean_delta = mean_delta,
        iqr_delta = iqr_delta,
        prop_pos = prop_pos,
        p_wilcoxon = p_wil,
        n_sign = n_sign,
        p_sign = p_sign,
        p_perm_signflip = p_perm,
        p_ks = p_ks,
        null_q05 = qlo,
        null_q95 = qhi,
        prop_emp_outside_05_95 = prop_emp_outside
      )
    },
    by = model
  ]
  
  # attach Spearman city-rank correlation per model
  spear <- city_rank[
    ,
    .(spearman_city = cor(rho_emp_city, rho_null_city, method = "spearman")),
    by = model
  ]
  
  out <- merge(out, spear, by = "model", all.x = TRUE)
  
  # pretty p-values as strings (optional)
  out[
    ,
    `:=`(
      p_wilcoxon = fmt_p(p_wilcoxon),
      p_sign = fmt_p(p_sign),
      p_perm_signflip = fmt_p(p_perm_signflip),
      p_ks = fmt_p(p_ks)
    )
  ]
  
  # reorder columns for a clean table
  setcolorder(out, c(
    "model",
    "n_pairs",
    "med_delta", "iqr_delta", "mean_delta",
    "prop_pos",
    "spearman_city",
    "p_wilcoxon", "p_sign", "p_perm_signflip", "p_ks",
    "null_q05", "null_q95", "prop_emp_outside_05_95"
  ))
  
  out[]
}

#-----------------------------
# Run
#-----------------------------
test_table <- build_test_table(rho_cmp, B_perm = 5000, seed_perm = 1)
test_table

























