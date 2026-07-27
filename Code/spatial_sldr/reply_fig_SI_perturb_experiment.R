# rm(list = ls())

# Data from disturb_experiment.R
# 螖rho vs disturb (by scenario) with 95% CI


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'



#----------Load packages----------#
library(ggtext)     # element_markdown()


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

dset<-which(Dname.set%in%c("Atlanta","San Francisco"))
yindex<-2


#----------Run scenarios with multiple perturbation magnitudes----------#
disturb_set <- c(0.01, seq(0.10,0.90,0.10))
disturb_set <- c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)


# NOTE (EN): scenario_set includes both parameterized and non-parameterized scenarios.
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
  
  dist_list[[city_name]] <- read_disturb_rho(city_name, region_name, y_tag, datapath)
  
  base_list[[city_name]] <- read_baseline_rho(city_name, region_name, y_tag, datapath) %>%
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
  # labels = c(
  #   "Uniform scaling",
  #   "Multiplicative noise",
  #   "Permutation null",
  #   "Shell-localize",
  #   "Shell-delocalize",
  #   "Core downweight",
  #   "Stratified rewiring"
  # ),
  # col = c("#4E79A7", "#59A14F", "#9C755F", "#F28E2B", "#EDC948", "#B07AA1", "#E15759"),
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
  facet_wrap(~ msa, scales = "fixed", nrow = 1) +
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
  guides(color = guide_legend(nrow = 1, byrow = TRUE), fill  = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain",size=15),
        legend.position = "top")
ggsave(fig_perturb, filename = paste(figpath,"/msa/SI_rho_perturb_",Yname[yindex],".pdf",sep=""), width = 4*2, height = 4.5)



