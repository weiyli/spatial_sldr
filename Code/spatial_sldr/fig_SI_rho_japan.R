# rm(list = ls())

# Transferability of SLDR to the Japan mobility datasets


#----------Workpath----------#
setwd("D:/ood/")
codepath <- "D:/ood/Code/spatial_sldr/spatial_sldr"
datapath <- "D:/ood/Data/spatial_sldr/japan"
figpath <- "D:/ood/Figure/spatial_sldr/msa"
dir.create(figpath, showWarnings = FALSE, recursive = TRUE)


#----------Load global variables and functions----------#
d <<- 1
source(file.path(codepath, "sldr_global_vars_funs.R"))


#----------Settings----------#
grid.info <- data.frame(
  breaks = c(0.5, 1, 2),
  labels = c("0.5 km", "1 km", "2 km")
)

# Days 0--59 are the 60-day business-as-usual (no-disaster) period;
# days 60--74 are the 15-day Emergency period.
no.disaster.days <- 0:59

#----------Read and prepare results from the no-disaster period----------#
japan.all <- fread(file.path(datapath, "SLDR_params_TO_in_OO.csv"))
japan.all <- japan.all[day %in% no.disaster.days]
japan.all[, grid_size_label := factor(
  paste0(grid_size_km, " km"),
  levels = grid.info$labels
)]
japan.all[,model := index]

# Panel b uses the best-performing exponential model. For each grid
# resolution there are four points (one per city), summarized over the
# 60 business-as-usual days; error bars show the interquartile range.
rho.japan <- japan.all[index == "exp"]
rho.japan.city <- rho.japan[, .(
  median_rho = median(rho, na.rm = TRUE),
  q25_rho = quantile(rho, 0.25, na.rm = TRUE),
  q75_rho = quantile(rho, 0.75, na.rm = TRUE)
), by = .(city, grid_size_km, grid_size_label)]
rho.japan.reference <- rho.japan.city[, .(
  reference_rho = mean(median_rho, na.rm = TRUE)
), by = .(grid_size_km, grid_size_label)]


# Use the baseline 1-km data for the daily model-error comparison.
error.japan.1km <- japan.all[grid_size_km == 1]

#----------a. Daily model error, faceted by city----------#
fig.error.japan <- ggplot(error.japan.1km, aes(x = day, y = error, color = model,shape = model,group = model)) +
  geom_line(aes(color = model), linewidth = 1) +
  geom_point(aes(color = model, fill = model, shape=model), size=3, stroke = 1) +
  facet_wrap(~city, ncol = 4, scales = "fixed") +
  scale_fill_manual(name = "Model",
                    breaks = model.info$breaks,
                    labels = model.info$labels,
                    values = scales::alpha(model.info$col,alpha=1))+ 
  scale_color_manual(name = "Model",
                     breaks = model.info$breaks,
                     labels = model.info$labels,
                     values = model.info$col) + 
  scale_shape_manual(name = "Model",
                     breaks = model.info$breaks,
                     labels = model.info$labels,
                     values = model.info$shape) + 
  scale_x_continuous(breaks = c(0, 15, 30, 45, 59)) +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))+
  labs(x = "Date", y = "RMSE-RMSPE", title = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        panel.spacing = unit(0.1, "lines"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
        strip.text = element_text(face = "plain",size=15),
        axis.text.y = element_text(hjust = 0),
        # axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top",
        legend.key.size = unit(1, "cm")) 


#----------b. City-level rho over the 60-day no-disaster period----------#
city.info <- data.frame(breaks = c("cityA", "cityB", "cityC", "cityD"),
                        labels = c("cityA", "cityB", "cityC", "cityD"),
                        col = col.all[13:16])

fig.rho.japan <- ggplot(rho.japan.city, aes(x = city, y = median_rho, color = city, fill = city)) +
  geom_hline(
    data = rho.japan.reference,
    aes(yintercept = reference_rho),
    color = "#69b3a2",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_errorbar(aes(ymin = q25_rho, ymax = q75_rho), width = 0.2, linewidth = 0.5) +  
  geom_point(aes(fill = city, color = city), shape = 21, size = 4, stroke = 1) +
  facet_wrap(~grid_size_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(name =  'City',
                    breaks = city.info$breaks,
                    labels = city.info$labels,
                    values = scales::alpha(city.info$col,alpha=1))+
  scale_color_manual(name =  'City',
                     breaks = city.info$breaks,
                     labels = city.info$labels,
                     values = city.info$col) +
  labs(x = "City", y = TeX("Spatial range exponent $\\rho$")) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray"),
    strip.text = element_text(face = "plain"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )


#----------Together----------#
fig.SI.rho.japan <- (fig.error.japan / fig.rho.japan) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20))
ggsave(fig.SI.rho.japan, filename = file.path(figpath, "SI_rho_japan.pdf"), width = 16, height = 5.5*2)
