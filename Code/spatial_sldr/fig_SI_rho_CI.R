# rm(list = ls())

# rho + 95% bootstrap CI（exp kernel）
# profile loss L(ρ) curve

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
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

#----------set the global function pointer used by homo/power/exp.param.est----------#
error_test <<- error_specs[["rmse_rmspe"]]$fun




#---- choose one region / yindex ----
# 你在 bootstrap 代码里 for(s in 2:Nregion)；这里按你实际的 region 名字填
region_use <- "msa"        # <-- 改成 region[s] 对应的文件夹名，比如 "msa" / "county" / "urban" 等
yname_use  <- Yname[2]     # 你循环 yindex=2:2

#---- choose MSAs to plot (建议选 3~4 个代表) ----
msa_pick <- Dname.msa[-8]  # <-- 你想展示哪些

#---- read and bind CI files ----
ci_all <- rbindlist(lapply(msa_pick, function(msa){
  f <- file.path(datapath, region_use, "param",
                 paste0(msa, "_SLDR_params_CI_", yname_use, "_1.csv"))
  if (!file.exists(f)) stop("Missing file: ", f)
  dt <- fread(f)
  dt[, msa := msa]
  dt
}), fill = TRUE)

#---- keep exp kernel ----
ci_exp <- ci_all[index == "exp"]
ci_exp[, day := as.Date(day)]

#---- plot ----
fig.ci <- ggplot(ci_exp, aes(x = day, y = rho, group = msa)) +
  geom_ribbon(aes(ymin = rho_ci_low, ymax = rho_ci_high, fill = msa), alpha = 0.20, color = NA) +
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
  scale_shape_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$shape
  ) +
  labs(x = "Date", y = TeX("Spatial range exponent $\\rho$ with CI")) +
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
  
#----------Model performance for Atlanta and Boston----------#
yindex<-2
fig.rho.ci <- fig.ci + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 20))
ggsave(fig.ci,filename = paste(figpath,"/msa/SI_rho_CI_",Yname[yindex],".pdf",sep=""), width =4*4, height = 4.5*3)


