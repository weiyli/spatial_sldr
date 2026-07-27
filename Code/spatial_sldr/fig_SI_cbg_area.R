# rm(list = ls())

# [ CBG area distribution across MSAs ] [ CBG area and rho]
# Save CBG-level area data and draw an SI_cbg_area distribution figure.
# Use CBG area (mean or median) data from rho_geo_var.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(sf)
library(data.table)
library(dplyr)
library(ggplot2)

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


#----------CBG area data for each MSA----------#
cbg.area.msa <- NULL
for(d in 1:Nmsa){
  
  Dname <<- Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  
  block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.dt <- block.msa %>% sf::st_drop_geometry() %>% as.data.table()
  
  if(!("area_m2" %in% names(block.dt))){
    stop(paste("Missing area_m2 in", Dname))
  }
  
  cbg.id.col <- intersect(c("CensusBlockGroup", "cbg_fips", "GEOID", "geoid"), names(block.dt))[1]
  if(is.na(cbg.id.col)){
    block.dt[, cbg_id := as.character(.I)]
  }else{
    block.dt[, cbg_id := as.character(get(cbg.id.col))]
  }
  
  dis.area <- data.table(
    msa = Dname,
    cbg_id = block.dt$cbg_id,
    cbg_area_m2 = as.numeric(block.dt$area_m2),
    cbg_area_km2 = as.numeric(block.dt$area_m2)/10^6
  )
  
  if("pop" %in% names(block.dt)){
    dis.area[, pop := block.dt$pop]
  }
  if("county_fips" %in% names(block.dt)){
    dis.area[, county_fips := as.character(block.dt$county_fips)]
  }
  
  cbg.area.msa <- plyr::rbind.fill(cbg.area.msa, dis.area)
}

cbg.area.msa <- as.data.table(cbg.area.msa)
cbg.area.msa <- cbg.area.msa[is.finite(cbg_area_km2)&cbg_area_km2>0]

#----------Save CBG-level area data and summary----------#
write.csv(cbg.area.msa, file = paste(datapath,"/msa/cbg_area.csv",sep=""), row.names = FALSE)

cbg.area.summary <- cbg.area.msa[, .(
  n_cbg = .N,
  mean_cbg_area_km2 = mean(cbg_area_km2, na.rm = TRUE),
  median_cbg_area_km2 = median(cbg_area_km2, na.rm = TRUE),
  sd_cbg_area_km2 = sd(cbg_area_km2, na.rm = TRUE),
  q25_cbg_area_km2 = quantile(cbg_area_km2, 0.25, na.rm = TRUE),
  q75_cbg_area_km2 = quantile(cbg_area_km2, 0.75, na.rm = TRUE),
  q95_cbg_area_km2 = quantile(cbg_area_km2, 0.95, na.rm = TRUE),
  max_cbg_area_km2 = max(cbg_area_km2, na.rm = TRUE)
), by = msa]

write.csv(cbg.area.summary, file = paste(datapath,"/msa/cbg_area_summary.csv",sep=""), row.names = FALSE)


#----------Read saved CBG area data for plotting----------#
cbg.area.msa <- fread(
  file = paste(datapath,"/msa/cbg_area.csv",sep="")
)
cbg.area.msa <- cbg.area.msa[
  is.finite(cbg_area_km2) & cbg_area_km2 > 0
]
cbg.area.msa[, log_cbg_area := log10(cbg_area_km2)]

cbg.area.mean <- cbg.area.msa[, .(
  mean_cbg_area_km2 = mean(cbg_area_km2, na.rm = TRUE),
  median_cbg_area_km2 = median(cbg_area_km2, na.rm = TRUE)
), by = msa]

cbg.area.mean[, label_area := sprintf(
  "atop(plain(median) == %.2f~plain(km)^2, plain(mean) == %.2f~plain(km)^2)",
  median_cbg_area_km2,
  mean_cbg_area_km2
)]

#----------SI msa level: density plot of CBG polygon area----------#
fig.cbg.area.si <- ggplot(cbg.area.msa, aes(log_cbg_area)) +
  geom_density(aes(fill = msa), alpha = 0.6, color = NA) +
  geom_text(
    data = cbg.area.mean,
    aes(
      x = log10(median_cbg_area_km2),
      y = Inf,
      label = label_area,
      color = msa
    ),
    parse = TRUE,
    vjust = 1.5, hjust = -0.15, size = 4
  ) +
  facet_wrap(~msa, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    breaks = log10(c(0.1, 1, 10, 100)),
    labels = c("0.1", "1", "10", "100")
  ) +
  scale_fill_manual(values = msa.info$col) +
  scale_color_manual(values = msa.info$col) +
  labs(x = TeX('CBG area (km$^2$, log scale)'),
       y = "Density") +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed", linewidth = 0.3),
        panel.grid.major.y = element_line(color = "gray90", linetype = "dashed", linewidth = 0.3),
        panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(face = "plain"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
        legend.position = "none")
ggsave(fig.cbg.area.si, filename = paste(figpath,"/msa/SI_cbg_area.pdf",sep=""),width = 6.8*2, height = 5*2)



#----------Correlation between mean/median CBG area and rho----------#
geo_var <- fread(file = paste(datapath,"/msa/geo_var.csv",sep=""))
rho.msa <- fread(file = paste(datapath,"/msa/SLDR_params_",Yname[2],".csv",sep=""))
rho.msa <- rho.msa[index=="exp"&period==pervalue[1]]
rho.msa <- rho.msa[, .(
  median_rho = median(rho, na.rm = TRUE)
), by = msa]

rho.cbg.area <- merge(
  rho.msa,
  geo_var[, .(
    msa,
    mean_cbg_area_km2,
    median_cbg_area_km2
  )],
  by = "msa"
)

rho.cbg.area.long <- melt(
  rho.cbg.area,
  id.vars = c("msa","median_rho"),
  measure.vars = c("mean_cbg_area_km2","median_cbg_area_km2"),
  variable.name = "cbg_area_stat",
  value.name = "cbg_area_km2"
)

cbg.area.corr <- rho.cbg.area.long[, {
  fit <- lm(median_rho ~ cbg_area_km2)
  pearson_test <- cor.test(cbg_area_km2, median_rho, method = "pearson")
  spearman_test <- cor.test(cbg_area_km2, median_rho, method = "spearman", exact = FALSE)
  .(
    R2 = summary(fit)$r.squared,
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    spearman_r = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value
  )
}, by = cbg_area_stat]

cbg.area.corr[, cbg_area_label := fifelse(
  cbg_area_stat=="mean_cbg_area_km2",
  "Mean CBG area",
  "Median CBG area"
)]
cbg.area.corr[, label := sprintf(
  "R^2 == %.2f*','~~Pearson~~italic(r) == %.2f",
  R2,
  pearson_r
)]

selected.cbg.area.var <- cbg.area.corr[
  which.max(abs(pearson_r)),
  cbg_area_stat
]
message("CBG area variable selected by absolute Pearson correlation: ",
        selected.cbg.area.var)

write.csv(
  cbg.area.corr,
  file = paste(datapath,"/msa/cbg_area_rho_correlation.csv",sep=""),
  row.names = FALSE
)

rho.cbg.area.long[, cbg_area_label := fifelse(
  cbg_area_stat=="mean_cbg_area_km2",
  "Mean CBG area",
  "Median CBG area"
)]
rho.cbg.area.long[, cbg_area_label := factor(
  cbg_area_label,
  levels = c("Mean CBG area","Median CBG area")
)]
cbg.area.corr[, cbg_area_label := factor(
  cbg_area_label,
  levels = c("Mean CBG area","Median CBG area")
)]

fig.cbg.area.rho.si <- ggplot(
  rho.cbg.area.long,
  aes(x = cbg_area_km2, y = median_rho)
) +
  geom_smooth(method = "lm", color = "#69b3a2", se = TRUE) +
  geom_point(
    aes(fill = msa, color = msa, shape = msa),
    size = 3, stroke = 0.5
  ) +
  geom_text(
    data = cbg.area.corr,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.1, vjust = 1.3, size = 4,
    color = "#69b3a2"
  ) +
  facet_wrap(~cbg_area_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col,alpha=0.8)
  ) +
  scale_color_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$col
  ) +
  scale_shape_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$shape
  ) +
  labs(
    x = TeX('CBG area (km$^2$)'),
    y = TeX('Spatial range exponent $\\rho$')
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    legend.position = "right"
  )

ggsave(fig.cbg.area.rho.si, filename = paste(figpath,"/msa/SI_cbg_area_rho.pdf",sep=""), width = 4.6*2, height = 4*1)



#----------Tract area data for each MSA----------#
# Reviewer noted that changing aggregation from CBG to tract reshuffles city ordering.
# We therefore summarize tract-scale polygon areas and relate them to tract-level rho.
tract.area.msa <- NULL
for(d in 1:Nmsa){
  
  Dname <<- Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  
  block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.dt <- block.msa %>% sf::st_drop_geometry() %>% as.data.table()
  
  if(!("area_m2" %in% names(block.dt))){
    stop(paste("Missing area_m2 in", Dname))
  }
  
  cbg.id.col <- intersect(c("CensusBlockGroup", "cbg_fips", "GEOID", "geoid"), names(block.dt))[1]
  if(is.na(cbg.id.col)){
    stop(paste("Missing CBG identifier for tract aggregation in", Dname))
  }
  
  block.dt[, tract_fips := substr(sprintf("%012s", as.character(get(cbg.id.col))), 1, 11)]
  
  dis.area <- block.dt[, .(
    tract_area_m2 = sum(as.numeric(area_m2), na.rm = TRUE),
    n_cbg = .N
  ), by = tract_fips]
  dis.area[, `:=`(
    msa = Dname,
    tract_area_km2 = tract_area_m2/10^6
  )]
  
  tract.area.msa <- plyr::rbind.fill(tract.area.msa, dis.area)
}

tract.area.msa <- as.data.table(tract.area.msa)
tract.area.msa <- tract.area.msa[is.finite(tract_area_km2)&tract_area_km2>0]

#----------Save tract-level area data and summary----------#
write.csv(tract.area.msa, file = paste(datapath,"/tract_msa/tract_area.csv",sep=""), row.names = FALSE)

tract.area.summary <- tract.area.msa[, .(
  n_tract = .N,
  mean_tract_area_km2 = mean(tract_area_km2, na.rm = TRUE),
  median_tract_area_km2 = median(tract_area_km2, na.rm = TRUE),
  sd_tract_area_km2 = sd(tract_area_km2, na.rm = TRUE),
  q25_tract_area_km2 = quantile(tract_area_km2, 0.25, na.rm = TRUE),
  q75_tract_area_km2 = quantile(tract_area_km2, 0.75, na.rm = TRUE),
  q95_tract_area_km2 = quantile(tract_area_km2, 0.95, na.rm = TRUE),
  max_tract_area_km2 = max(tract_area_km2, na.rm = TRUE)
), by = msa]

write.csv(tract.area.summary, file = paste(datapath,"/tract_msa/tract_area_summary.csv",sep=""), row.names = FALSE)

#----------Read saved tract area data for plotting----------#
tract.area.msa <- fread(
  file = paste(datapath,"/tract_msa/tract_area.csv",sep="")
)
tract.area.msa <- tract.area.msa[
  is.finite(tract_area_km2) & tract_area_km2 > 0
]
tract.area.msa[, log_tract_area := log10(tract_area_km2)]

tract.area.mean <- tract.area.msa[, .(
  mean_tract_area_km2 = mean(tract_area_km2, na.rm = TRUE),
  median_tract_area_km2 = median(tract_area_km2, na.rm = TRUE)
), by = msa]

tract.area.mean[, label_area := sprintf(
  "atop(plain(median) == %.2f~plain(km)^2, plain(mean) == %.2f~plain(km)^2)",
  median_tract_area_km2,
  mean_tract_area_km2
)]

#----------SI msa level: density plot of tract polygon area----------#
fig.tract.area.si <- ggplot(tract.area.msa, aes(x = log_tract_area)) +
  geom_density(aes(fill = msa), alpha = 0.6, color = NA) +
  geom_text(
    data = tract.area.mean,
    aes(
      x = log10(median_tract_area_km2) + 0.08,
      y = Inf,
      label = label_area,
      color = msa
    ),
    parse = TRUE,
    vjust = 1.2, hjust = 0, size = 4
  ) +
  facet_wrap(~msa, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    breaks = log10(c(0.1, 1, 10, 100, 1000)),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_fill_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col, alpha = 0.8)
  ) +
  scale_color_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$col
  ) +
  labs(
    x = TeX('Tract area (km$^2$, log scale)'),
    y = "Density"
  ) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed", linewidth = 0.3),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed", linewidth = 0.3),
        panel.spacing = unit(0.1, "lines"),
        strip.text = element_text(face = "plain"),
        axis.line.x.bottom = element_line(color = "black", linewidth = 0.5),
        legend.position = "none")

ggsave(fig.tract.area.si, filename = paste(figpath,"/msa/SI_tract_area.pdf",sep=""), width = 6.8*2, height = 5*2)


#----------Tract area and tract-level rho across MSAs----------#
date.gif <- date.sep(Datevalue, durdate)
date.gif$day <- as.Date(date.gif$day)

rho.tract <- NULL
for(d in 1:Nmsa){
  Dname <- Dname.set[d]
  fpath <- paste(datapath, "/tract_msa/param/", Dname, "_SLDR_params_TO_in_OO.csv", sep = "")
  if(!file.exists(fpath)) next
  
  sldr.params <- fread(fpath)
  sldr.params[, day := as.Date(day)]
  sldr.params <- left_join(
    sldr.params,
    date.gif,
    by = "day"
  ) %>% as.data.table()
  sldr.params <- sldr.params[index=="exp" & period==pervalue[1]]
  sldr.params[, msa := Dname]
  rho.tract <- plyr::rbind.fill(rho.tract, sldr.params)
}
rho.tract <- as.data.table(rho.tract)
rho.tract <- rho.tract[, .(
  median_rho = median(rho, na.rm = TRUE)
), by = msa]

rho.tract.area <- merge(
  rho.tract,
  tract.area.summary[, .(
    msa,
    mean_tract_area_km2,
    median_tract_area_km2
  )],
  by = "msa"
)

rho.tract.area.long <- melt(
  rho.tract.area,
  id.vars = c("msa","median_rho"),
  measure.vars = c("mean_tract_area_km2","median_tract_area_km2"),
  variable.name = "tract_area_stat",
  value.name = "tract_area_km2"
)

tract.area.corr <- rho.tract.area.long[, {
  fit <- lm(median_rho ~ tract_area_km2)
  pearson_test <- cor.test(tract_area_km2, median_rho, method = "pearson")
  spearman_test <- cor.test(tract_area_km2, median_rho, method = "spearman", exact = FALSE)
  .(
    R2 = summary(fit)$r.squared,
    pearson_r = unname(pearson_test$estimate),
    pearson_p = pearson_test$p.value,
    spearman_r = unname(spearman_test$estimate),
    spearman_p = spearman_test$p.value
  )
}, by = tract_area_stat]

tract.area.corr[, tract_area_label := fifelse(
  tract_area_stat=="mean_tract_area_km2",
  "Mean Tract area",
  "Median Tract area"
)]
tract.area.corr[, label := sprintf(
  "R^2 == %.2f*','~~Pearson~~italic(r) == %.2f",
  R2,
  pearson_r
)]

write.csv(
  tract.area.corr,
  file = paste(datapath,"/tract_msa/tract_area_rho_correlation.csv",sep=""),
  row.names = FALSE
)

rho.tract.area.long[, tract_area_label := fifelse(
  tract_area_stat=="mean_tract_area_km2",
  "Mean Tract area",
  "Median Tract area"
)]
rho.tract.area.long[, tract_area_label := factor(
  tract_area_label,
  levels = c("Mean Tract area","Median Tract area")
)]
tract.area.corr[, tract_area_label := factor(
  tract_area_label,
  levels = c("Mean Tract area","Median Tract area")
)]

fig.tract.area.rho.si <- ggplot(
  rho.tract.area.long,
  aes(x = tract_area_km2, y = median_rho)
) +
  geom_smooth(method = "lm", color = "#69b3a2", se = TRUE) +
  geom_point(
    aes(fill = msa, color = msa, shape = msa),
    size = 3, stroke = 0.5
  ) +
  geom_text(
    data = tract.area.corr,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.1, vjust = 1.3, size = 4,
    color = "#69b3a2"
  ) +
  facet_wrap(~tract_area_label, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col,alpha=0.8)
  ) +
  scale_color_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$col
  ) +
  scale_shape_manual(
    name = 'MSA',
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = msa.info$shape
  ) +
  labs(
    x = TeX('Tract area (km$^2$)'),
    y = TeX('Spatial range exponent $\\rho$')
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(color = "gray", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
    strip.text = element_text(face = "plain", size = 15),
    legend.position = "right"
  )

ggsave(fig.tract.area.rho.si, filename = paste(figpath,"/msa/SI_tract_area_rho.pdf",sep=""), width = 4.6*2, height = 4*1)



#----------CBG rho versus tract rho: Before and During----------#
# Follow the MAUP figure exactly: daily paired observations are shown, while
# Spearman rank correlations are calculated from city medians within each period.
date.gif <- date.sep(Datevalue, durdate)
yindex <- 2
aggregation.scales <- c("msa","tract_msa")

read.scale.rho <- function(scale, city, yname){
  candidates <- c(
    file.path(
      datapath, scale, "param",
      paste0(city, "_SLDR_params_", yname, ".csv")
    ),
    file.path(
      datapath, scale, "params",
      paste0(city, "_SLDR_params_", yname, ".csv")
    )
  )
  f <- candidates[file.exists(candidates)][1]
  if(is.na(f)) return(NULL)
  
  x <- fread(f)
  x[, `:=`(msa = city, scale = scale)]
  x
}

rho.all.period <- rbindlist(lapply(aggregation.scales, function(sc){
  rbindlist(
    lapply(
      Dname.msa,
      function(ct) read.scale.rho(sc, ct, Yname[yindex])
    ),
    fill = TRUE
  )
}), fill = TRUE)

rho.all.period <- rho.all.period[index=="exp"]
rho.all.period[, rho := as.numeric(rho)]
rho.all.period <- rho.all.period[!is.na(rho)]
rho.all.period[, day := as.Date(day)]
rho.all.period <- left_join(
  rho.all.period,
  date.gif,
  by = "day"
) %>% as.data.table()

rho.scale.period <- copy(rho.all.period)
if("week" %in% names(rho.scale.period)){
  rho.scale.period[, week := NULL]
}
rho.scale.period[scale=="msa", scale := "cbg_msa"]

rho.period.wide <- dcast(
  rho.scale.period,
  msa+day+period~scale,
  value.var = "rho"
)

rho.pair.scale <- rho.period.wide[
  !is.na(cbg_msa)&!is.na(tract_msa),
  .(
    msa,
    day,
    period,
    rho_cbg = cbg_msa,
    rho_tract = tract_msa
  )
]
rho.pair.scale[, period := stringr::str_to_title(as.character(period))]
rho.pair.scale[, period := factor(
  period,
  levels = c("Before","During")
)]

rho.pair.scale.median <- rho.pair.scale[, .(
  rho_cbg = median(rho_cbg, na.rm = TRUE),
  rho_tract = median(rho_tract, na.rm = TRUE)
), by = .(msa, period)]

xy.range <- range(
  c(rho.pair.scale.median$rho_cbg, rho.pair.scale.median$rho_tract),
  na.rm = TRUE
)

stats.scale <- rho.pair.scale.median[
  is.finite(rho_cbg)&is.finite(rho_tract),
  {
    fit <- lm(rho_tract~rho_cbg)
    delta <- rho_tract-rho_cbg

    sp <- suppressWarnings(cor.test(
      rho_cbg,
      rho_tract,
      method = "spearman",
      exact = FALSE
    ))
    
    list(
      b0 = coef(fit)[1],
      b1 = coef(fit)[2],
      median_delta = median(delta, na.rm = TRUE),
      sp_rho = unname(sp$estimate),
      sp_p = sp$p.value,
      x_pos = min(rho_cbg, na.rm = TRUE),
      y_pos = max(rho_tract, na.rm = TRUE)
    )
  },
  by = period
]

stats.scale[, `:=`(
  lab1 = sprintf(
    "rho[Tract] == %.2f~rho[CBG] %+ .2f",
    round(b1,2),
    round(b0,2)
  ),
  lab2 = sprintf(
    "Median(rho[Tract]-rho[CBG]) == %.2f",
    round(median_delta,2)
  ),
  lab3 = sprintf(
    "Spearman~~italic(r)[S] == %.2f",
    round(sp_rho,2)
  )
)]

ann.scale <- data.table::rbindlist(list(
  stats.scale[, .(period, x = x_pos, y = y_pos, label = lab1)],
  stats.scale[, .(period, x = x_pos, y = y_pos, label = lab2)],
  stats.scale[, .(period, x = x_pos, y = y_pos, label = lab3)]
), use.names = TRUE)

dy.scale <- 0.1*diff(xy.range)
ann.scale[, line := rep(1:3, each = nrow(stats.scale))]
ann.scale[, y := y-(line-1)*dy.scale]

write.csv(
  stats.scale,
  file = paste(
    datapath,
    "/tract_msa/cbg_tract_rho_period.csv",
    sep = ""
  ),
  row.names = FALSE
)

fig.cbg.tract.aggregation.si <- ggplot(
  rho.pair.scale.median,
  aes(x = rho_cbg, y = rho_tract)
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 1,
    color = "#69b3a2"
  ) +
  geom_point(
    aes(fill = msa, color = msa, shape = msa),
    size = 3,
    stroke = 0.5
  ) +
  geom_text(
    data = ann.scale,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3,
    color = "#69b3a2",
    parse = TRUE
  ) +
  facet_wrap(~period, scales = "fixed", nrow = 1) +
  scale_fill_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col, alpha = 0.5)
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
  scale_x_continuous(limits = xy.range) +
  scale_y_continuous(limits = xy.range) +
  labs(
    x = TeX('Median $\\rho_{CBG}$'),
    y = TeX('Median $\\rho_{Tract}$')
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(
      color = "gray",
      fill = NA,
      linewidth = 1
    ),
    strip.background = element_rect(
      fill = "#f0f0f0",
      color = "gray",
      linewidth = 1
    ),
    strip.text = element_text(face = "plain", size = 15),
    legend.position = "right"
  )

ggsave(
  fig.cbg.tract.aggregation.si,
  filename = paste(
    figpath,
    "/msa/SI_cbg_tract_median_rho.pdf",
    sep = ""
  ),
  width = 4.6*2,
  height = 4
)


#----------Aggregation intensity and period-specific changes in rho----------#
aggregation.period <- merge(
  rho.pair.scale.median,
  cbg.area.summary[, .(
    msa,
    mean_cbg_area_km2,
    median_cbg_area_km2
  )],
  by = "msa"
)
aggregation.period <- merge(
  aggregation.period,
  tract.area.summary[, .(
    msa,
    mean_tract_area_km2,
    median_tract_area_km2
  )],
  by = "msa"
)

aggregation.period[, `:=`(
  rank_rho_cbg = frank(
    rho_cbg,
    ties.method = "average"
  ),
  rank_rho_tract = frank(
    rho_tract,
    ties.method = "average"
  )
), by = period]

aggregation.rank.long <- rbindlist(list(
  aggregation.period[, .(
    msa, period, rho_cbg, rho_tract,
    rank_rho_cbg, rank_rho_tract,
    area_stat = "Mean area",
    cbg_area_km2 = mean_cbg_area_km2,
    tract_area_km2 = mean_tract_area_km2
  )],
  aggregation.period[, .(
    msa, period, rho_cbg, rho_tract,
    rank_rho_cbg, rank_rho_tract,
    area_stat = "Median area",
    cbg_area_km2 = median_cbg_area_km2,
    tract_area_km2 = median_tract_area_km2
  )]
))
aggregation.rank.long[, `:=`(
  rank_cbg_area = frank(cbg_area_km2, ties.method = "average"),
  rank_tract_area = frank(tract_area_km2, ties.method = "average")
), by = .(period, area_stat)]
aggregation.rank.long[, `:=`(
  area_rank_shift = rank_tract_area-rank_cbg_area,
  rho_rank_shift = rank_rho_tract-rank_rho_cbg,
  area_stat = factor(area_stat, levels = c("Mean area", "Median area"))
)]

write.csv(
  aggregation.rank.long,
  file = paste(
    datapath,
    "/tract_msa/cbg_tract_period.csv",
    sep = ""
  ),
  row.names = FALSE
)

# Area-rank shift versus rho-rank shift.
stats.rank.shift <- aggregation.rank.long[
  is.finite(area_rank_shift)&is.finite(rho_rank_shift),
  {
    fit <- lm(rho_rank_shift~area_rank_shift)
    spearman.test <- suppressWarnings(cor.test(
      area_rank_shift,
      rho_rank_shift,
      method = "spearman",
      exact = FALSE
    ))
    list(
      R2 = summary(fit)$r.squared,
      spearman_r = unname(spearman.test$estimate),
      spearman_p = spearman.test$p.value
    )
  },
  by = .(area_stat, period)
]
stats.rank.shift[, `:=`(
  label = sprintf(
    "R^2 == %.3f*','~~Spearman~~italic(r)[S] == %.2f",
    R2,
    spearman_r
  )
)]

rank.shift.y.range <- range(
  aggregation.rank.long$rho_rank_shift,
  na.rm = TRUE
)
rank.shift.x.pos <- min(
  aggregation.period$area_rank_shift,
  na.rm = TRUE
)
rank.shift.y.pos <- max(rank.shift.y.range) +
  0.18*diff(rank.shift.y.range)

ann.rank.shift <- stats.rank.shift[, .(
  area_stat,
  period,
  x = -Inf,
  y = Inf,
  label
)]

write.csv(
  stats.rank.shift,
  file = paste(
    datapath,
    "/tract_msa/area_rho_rank_shift.csv",
    sep = ""
  ),
  row.names = FALSE
)

fig.rank.shift <- ggplot(
  aggregation.rank.long,
  aes(x = area_rank_shift, y = rho_rank_shift)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "gray60"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "gray60"
  ) +
  geom_smooth(
    method = "lm",
    color = "#69b3a2",
    se = TRUE
  ) +
  geom_point(
    aes(fill = msa, color = msa, shape = msa),
    size = 3,
    stroke = 0.5
  ) +
  geom_text(
    data = ann.rank.shift,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.25,
    size = 4,
    color = "#69b3a2",
    parse = TRUE
  ) +
  facet_grid(area_stat~period, scales = "fixed") +
  scale_fill_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col, alpha = 0.5)
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
  labs(
    x = "Change in area rank (Tract - CBG)",
    y = TeX('Change in $\\rho$ rank (Tract - CBG)')
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(
      color = "gray",
      fill = NA,
      linewidth = 1
    ),
    strip.background = element_rect(
      fill = "#f0f0f0",
      color = "gray",
      linewidth = 1
    ),
    strip.text = element_text(face = "plain", size = 15),
    legend.position = "right"
  )

ggsave(
  fig.rank.shift,
  filename = paste(
    figpath,
    "/msa/SI_cbg_tract_rank_shift.pdf",
    sep = ""
  ),
  width = 4.6*2,
  height = 4
)

#----------Change in rho ranking from CBG to Tract----------#
rho.rank.transition <- rbindlist(list(
  aggregation.period[, .(
    msa, period, aggregation = "CBG", rho_rank = rank_rho_cbg
  )],
  aggregation.period[, .(
    msa, period, aggregation = "Tract", rho_rank = rank_rho_tract
  )]
))
rho.rank.transition[, aggregation := factor(
  aggregation,
  levels = c("CBG", "Tract")
)]

write.csv(
  rho.rank.transition,
  file = paste(
    datapath,
    "/tract_msa/cbg_tract_rho_rank_transition.csv",
    sep = ""
  ),
  row.names = FALSE
)

fig.rho.rank.transition <- ggplot(
  rho.rank.transition,
  aes(
    x = aggregation,
    y = rho_rank,
    group = msa,
    color = msa
  )
) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(
    aes(fill = msa, shape = msa),
    size = 3,
    stroke = 0.5
  ) +
  facet_wrap(~period, nrow = 1) +
  scale_y_reverse(
    breaks = seq_len(Nmsa),
    minor_breaks = NULL
  ) +
  scale_fill_manual(
    name = "MSA",
    breaks = msa.info$dname,
    labels = msa.info$dname,
    values = scales::alpha(msa.info$col, alpha = 0.7)
  ) +
  scale_color_manual(
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
  labs(
    x = "Spatial aggregation",
    y = TeX('$\\rho$ rank')
  ) +
  theme_wy() +
  theme(
    panel.background = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray90",
      linetype = "dashed",
      linewidth = 0.3
    ),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(
      color = "gray",
      fill = NA,
      linewidth = 1
    ),
    strip.background = element_rect(
      fill = "#f0f0f0",
      color = "gray",
      linewidth = 1
    ),
    strip.text = element_text(face = "plain", size = 15),
    legend.position = "right"
  )

ggsave(
  fig.rho.rank.transition,
  filename = paste(
    figpath,
    "/msa/SI_cbg_tract_rho_rank_transition.pdf",
    sep = ""
  ),
  width = 4.6*2,
  height = 4
)


#--------------------------------------#
# Together  #
#--------------------------------------#
fig.cbg.tract.area <- (
  fig.tract.area.rho.si/
    fig.rho.rank.transition/
    fig.rank.shift
) +
  plot_layout(heights = c(1, 1, 1.8)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 20))
ggsave(
  fig.cbg.tract.area,
  filename = paste(figpath,"/msa/SI_cbg_tract_area.pdf",sep=""),
  width = 4.6*2,
  height = 4*3
)


