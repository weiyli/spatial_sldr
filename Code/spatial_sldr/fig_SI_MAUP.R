# rm(list = ls())

# Data from sldr_fit_MAUP.R
# Discussion of spatial scale and administrative units


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'

#----------Load packages----------#
library(car)
library(broom)
library(ggtext)


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
# Scatter plot of rho for different scales  #
#--------------------------------------#
yindex <- 2
# scale_info <- data.frame(
#   breaks = c("msa","tract_msa"),
#   labels = c("CBG","Tract"),
#   col    = c('#d25774','#72567a'),
#   shape  = c(21, 24)
# )
scales <- c("msa","tract_msa")
read_one <- function(scale, city, yname){
  candidates <- c(
    file.path(datapath, scale, "param", paste0(city, "_SLDR_params_", yname, ".csv")),
    file.path(datapath, scale, "params", paste0(city, "_SLDR_params_", yname, ".csv"))
  )
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) return(NULL)
  
  x <- fread(f)
  x[, `:=`(msa = city, scale = scale)]
  x
}

rho_all <- rbindlist(lapply(scales, function(sc){
  rbindlist(lapply(Dname.msa, function(ct) read_one(sc, ct, Yname[yindex])), fill = TRUE)
}), fill = TRUE)
rho_all <- rho_all[index == "exp"]
rho_all[, rho := as.numeric(rho)]
rho_all <- rho_all[!is.na(rho)]
# period
rho_all$day <- as.Date(rho_all$day)
rho_all <- left_join(rho_all,date.gif,by='day')%>%as.data.table
rho_scale <- rho_all[, week := NULL]
rho_scale[scale == "msa", scale := "cbg_msa"]
rho_period <- dcast(
  rho_scale,
  msa + day + period ~ scale,
  value.var = "rho"
)
rho_period[, delta_rho := tract_msa - cbg_msa]
rho_pair_scale <- rho_period[
  !is.na(cbg_msa) & !is.na(tract_msa),
  .(msa, day, period, rho_cbg = cbg_msa, rho_tract = tract_msa)
]
rho_pair_scale[, period := stringr::str_to_title(as.character(period))]
xy.range <- range(c(rho_pair_scale$rho_cbg, rho_pair_scale$rho_tract), na.rm = TRUE)


#--------------------------------------#
# A) text annotation for fig_scale (3 lines, parsed) + Spearman (city-level)
#--------------------------------------#
stats_scale <- rho_pair_scale[
  is.finite(rho_cbg) & is.finite(rho_tract),
  {
    #--- (1) day-level linear calibration (within period, pooled) ---#
    fit <- lm(rho_tract ~ rho_cbg)
    
    # identity test: intercept=0, slope=1
    lh <- car::linearHypothesis(fit, c("(Intercept)=0", "rho_cbg=1"))
    p_id <- lh$`Pr(>F)`[2]
    
    # delta test (paired difference, day-level)
    delta <- rho_tract - rho_cbg
    tt <- t.test(delta, mu = 0)
    md <- median(delta, na.rm = TRUE)
    p_d <- tt$p.value
    
    #--- (2) Spearman rank preservation (city-level medians) ---#
    # within this period subset (.SD), aggregate to city medians then Spearman
    city_dt <- .SD[
      , .(
        rho_cbg_med   = median(rho_cbg,   na.rm = TRUE),
        rho_tract_med = median(rho_tract, na.rm = TRUE)
      ),
      by = msa
    ]
    
    sp <- suppressWarnings(
      cor.test(city_dt$rho_cbg_med, city_dt$rho_tract_med,
               method = "spearman", exact = FALSE)
    )
    sp_rho <- unname(sp$estimate)
    sp_p   <- sp$p.value
    
    # anchor (top-left inside facet)
    x_pos <- min(rho_cbg, na.rm = TRUE)
    y_pos <- max(rho_tract, na.rm = TRUE)
    
    list(
      b0 = coef(fit)[1],
      b1 = coef(fit)[2],
      median_delta = md,
      p_delta = p_d,
      p_id = p_id,
      sp_rho = sp_rho,
      sp_p = sp_p,
      x_pos = x_pos,
      y_pos = y_pos
    )
  },
  by = period
]

stats_scale[, `:=`(
  lab1 = sprintf("rho[Tract] == %.2f~rho[CBG] %+ .2f", round(b1, 2), round(b0, 2)),
  lab2 = sprintf("Median(rho[Tract] - rho[CBG]) == %.2f", round(median_delta, 2)),
  lab3 = sprintf("Spearman~~italic(r)[S] == %.2f", round(sp_rho, 2))
)]

ann_scale <- data.table::rbindlist(list(
  stats_scale[, .(period, x = x_pos, y = y_pos, label = lab1)],
  stats_scale[, .(period, x = x_pos, y = y_pos, label = lab2)],
  stats_scale[, .(period, x = x_pos, y = y_pos, label = lab3)]
), use.names = TRUE)

dy_scale <- 0.1 * diff(xy.range)
ann_scale[, line := rep(1:3, each = nrow(stats_scale))]
ann_scale[, y := y - (line - 1) * dy_scale]


#--------------------------------------#
# Figure: rho_CBG vs rho_Tract
#--------------------------------------#
fig_scale <- ggplot(rho_pair_scale, aes(x = rho_cbg, y = rho_tract)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 3, stroke = 0.5) +
  geom_text(
    data = ann_scale,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3,
    color = "#69b3a2",
    parse = TRUE 
  ) +
  facet_wrap(~ period, scales = "fixed", nrow = 1) +
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
  labs(x = TeX('Daily $\\rho_{CBG}$'),
       y = TeX('Daily $\\rho_{Tract}$')) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")

#--------------------------------------#
# Scatter plot of rho for Rook vs Queen  #
#--------------------------------------#
# queen
rho_queen <- fread(paste0(datapath, "/msa/SLDR_params_", Yname[yindex], ".csv"))
rho_queen <- rho_queen[index == "exp", .(day, rho, msa)]
rho_queen[, adj := "Queen"]
# rook
read_one <- function(city, y_tag){
  x <- fread(file.path(datapath, "msa", "params", paste0(city, "_SLDR_params_rook_", y_tag, ".csv")))
  x <- x[index == "exp"]
  x[, `:=`(msa = city)]
  x[, .(day, rho, msa)] 
}
rho_rook <- rbindlist(lapply(Dname.msa, function(ct) read_one(ct, Yname[yindex])), fill = TRUE)
rho_rook <- rho_rook[, adj := "Rook"]
rho_all <- rbindlist(list(rho_queen, rho_rook), fill = TRUE, use.names = TRUE)
rho_all[, rho := as.numeric(rho)]
rho_all <- rho_all[!is.na(rho)]
# period
rho_all$day <- as.Date(rho_all$day)
rho_all <- left_join(rho_all,date.gif,by='day')%>%as.data.table
rho_adj <- as.data.table(rho_all)[, week := NULL]
# rho_all[scale == "msa", scale := "cbg_msa"]
rho_period <- dcast(
  rho_all,
  msa + day + period ~ adj,
  value.var = "rho"
)
rho_period[, delta_rho := Rook - Queen]
rho_pair_adj <- rho_period[
  !is.na(Queen) & !is.na(Rook),
  .(msa, day, period, rho_queen = Queen, rho_rook = Rook)
]
rho_pair_adj[, period := stringr::str_to_title(as.character(period))]
xy.range <- range(c(rho_pair_adj$rho_queen, rho_pair_adj$rho_rook), na.rm = TRUE)

#--------------------------------------#
# Adjacency robustness: Queen vs Rook
#--------------------------------------#
#--------------------------------------#
# B) text annotation for fig_adj (3 lines) + Spearman (city-level)
#--------------------------------------#
stats_adj <- rho_pair_adj[
  is.finite(rho_queen) & is.finite(rho_rook),
  {
    # (1) day-level linear calibration (within period, pooled)
    fit <- lm(rho_rook ~ rho_queen)
    
    # paired median shift (day-level)
    delta <- rho_rook - rho_queen
    md <- median(delta, na.rm = TRUE)
    
    # (2) Spearman rank preservation (city-level medians)
    city_dt <- .SD[
      , .(
        rho_queen_med = median(rho_queen, na.rm = TRUE),
        rho_rook_med  = median(rho_rook,  na.rm = TRUE)
      ),
      by = msa
    ]
    sp <- suppressWarnings(
      cor.test(city_dt$rho_queen_med, city_dt$rho_rook_med,
               method = "spearman", exact = FALSE)
    )
    sp_rho <- unname(sp$estimate)
    
    # anchor (top-left inside facet)
    x_pos <- min(rho_queen, na.rm = TRUE)
    y_pos <- max(rho_rook,  na.rm = TRUE)
    
    list(
      b0 = coef(fit)[1],
      b1 = coef(fit)[2],
      median_delta = md,
      sp_rho = sp_rho,
      x_pos = x_pos,
      y_pos = y_pos
    )
  },
  by = period
]

stats_adj[, `:=`(
  lab1 = sprintf("rho[Rook] == %.2f~rho[Queen] %+ .2f", round(b1,2), round(b0,2)),
  lab2 = sprintf("Median(rho[Rook] - rho[Queen]) == %.2f", round(median_delta,2)),
  lab3 = sprintf("Spearman~~italic(r)[S] == %.2f", round(sp_rho,2))
)]

ann_adj <- data.table::rbindlist(list(
  stats_adj[, .(period, x = x_pos, y = y_pos,      label = lab1)],
  stats_adj[, .(period, x = x_pos, y = y_pos,      label = lab2)],
  stats_adj[, .(period, x = x_pos, y = y_pos,      label = lab3)]
), use.names = TRUE)

dy <- 0.1 * diff(xy.range)   
ann_adj[, line := rep(1:3, each = nrow(stats_adj))]
ann_adj[, y := y - (line - 1) * dy]


#--------------------------------------#
# Figure: rho_Queen vs rho_Rook
#--------------------------------------#
fig_adj <- ggplot(rho_pair_adj, aes(x = rho_queen, y = rho_rook)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              linewidth = 1, color = "#69b3a2") +
  geom_point(aes(fill = msa, color = msa, shape = msa),
             size = 3, stroke = 0.5) +
  geom_text(
    data = ann_adj,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3,
    color = "#69b3a2",
    parse = TRUE 
  ) +
  facet_wrap(~ period, scales = "fixed", nrow = 1) +
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
  labs(x = TeX('Daily $\\rho_{Queen}$'),
       y = TeX('Daily $\\rho_{Rook}$')) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size = 15),
        legend.position = "right")

#--------------------------------------#
# Together  #
#--------------------------------------#
fig_rho <- (fig_scale/fig_adj) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig_rho, filename = paste(figpath,"/msa/SI_rho_MAUP_",Yname[yindex],".pdf",sep=""), width = 4.6*2, height = 4*2)






