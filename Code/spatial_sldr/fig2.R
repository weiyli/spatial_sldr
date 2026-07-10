# rm(list = ls())

# Interpreting the gap between null and empirical rho
# Null and empirical rho, correlation between Tij and dij from synthetic_flow.R; msa level data from sldr_fit.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(scales)
library(ggtext)

library(sf)          # read_sf() 
library(spdep)       # poly2nb() 

library(igraph)
library(ggraph)
library(tidygraph)
library(grid)


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
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")


#============================================================
# Part 1. Null rho and Empirical rho
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
#-----------------------------
# Define Low / High for legend text
#-----------------------------
cate.info <- data.frame(breaks = c("Low", "High"), 
                        labels = c("Low", "High"), 
                        col = c('#d25774', '#72567a'))
#-----------------------------
# Define Low / High based on empirical-null rho gap
#-----------------------------
# # threshold can be modified if needed
# msa_lower <- rho_cmp[rho_null < 0.35, unique(msa)]
rho_cmp[, rho_gap := rho_emp - rho_null]
# MSA-level average gap across days and null models
msa_gap <- rho_cmp[
  ,
  .(
    mean_gap = mean(rho_gap, na.rm = TRUE),
    sd_gap   = sd(rho_gap, na.rm = TRUE),
    n        = .N
  ),
  by = msa
]
# data-driven two-group classification
set.seed(123)
km <- kmeans(msa_gap[, .(mean_gap)], centers = 2, nstart = 50)
# Extract the two one-dimensional cluster centers
centers <- sort(as.numeric(km$centers[, "mean_gap"]))
# Threshold is the midpoint between the two centers
gap_threshold <- mean(centers)
msa_gap[, gap_group_raw := km$cluster]
# label groups by gap size
cluster_order <- msa_gap[,.(cluster_mean = mean(mean_gap)), by = gap_group_raw][order(cluster_mean)]
small_gap_cluster <- cluster_order$gap_group_raw[1]
large_gap_cluster <- cluster_order$gap_group_raw[2]
msa_gap[,PModel := ifelse(gap_group_raw == small_gap_cluster, "High", "Low")] # Same as: msa_gap[,PModel := ifelse(mean_gap<gap_threshold, "High","Low")]
# merge back to daily comparison data
rho_cmp <- merge(
  rho_cmp,
  msa_gap[, .(msa, mean_gap, PModel)],
  by = "msa",
  all.x = TRUE
)
# check classification
msa_gap[order(mean_gap)]
msa_lower <- rho_cmp[PModel == "Low", unique(msa)]

cate_color_map <- setNames(cate.info$col, cate.info$breaks)
msa_cate <- data.table(msa = msa.info$dname)
msa_cate[, cate := ifelse(msa %in% msa_lower, "Low", "High")]
# color only the legend text
msa_cate[, label_col := paste0(
  "<span style='color:",
  cate_color_map[cate],
  "'>", msa, "</span>"
)]
label_vec <- setNames(msa_cate$label_col, msa_cate$msa)
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
                    labels = label_vec,
                    values = scales::alpha(msa.info$col, alpha = 0.5)) +
  scale_color_manual(name = "MSA",
                     breaks = msa.info$dname,
                     labels = label_vec,
                     values = msa.info$col) +
  scale_shape_manual(name = "MSA",
                     breaks = msa.info$dname,
                     labels = label_vec,
                     values = msa.info$shape) +
  scale_x_continuous(limits = xy.range) +
  scale_y_continuous(limits = xy.range) +
  labs(x = TeX("Daily $\\rho_{Empirical}$"), y = TeX("Daily $\\rho_{Model}$")) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain",size = 15),
        legend.text = ggtext::element_markdown(),
        legend.position = "right")


#============================================================
# Part 2. Interpreting the gap between null and empirical rho
#============================================================
corr <- fread(paste0(datapath, "/",region_name,"/flow_dist_corr.csv"))
# rho_null < 0.35
corr[, cate := ifelse(msa %in% msa_lower, "Low", "High")]
#-----------------------------
# Figure: tij and dij for all msas
#-----------------------------
corr_summary <- corr[, .(
  med = median(corr_value),
  lower_ci = quantile(corr_value, 0.25),
  upper_ci = quantile(corr_value, 0.75)
), by = .(msa,corr_name)]
corr_summary[, cate := ifelse(msa %in% msa_lower, "Low", "High")][, cate := factor(cate, levels = c("Low", "High"))]
cate_color_map <- setNames(cate.info$col, cate.info$breaks)
corr_summary[, msa_col := paste0(
  "<span style='color:",
  cate_color_map[cate],
  "'>", msa, "</span>"
)]
# # order by cate
# ordered_labels <- corr_summary[corr_name == "spearman_adjusted"][
#   order(cate, msa),
#   unique(msa_col)
# ]
# corr_summary[, msa_col := factor(msa_col, levels = ordered_labels)]

# order by msa
ordered_labels <- corr_summary[corr_name == "spearman_adjusted"][
  order(msa),
  unique(msa_col)
]
corr_summary[, msa_col := factor(msa_col, levels = ordered_labels)]

corr_summary[corr_name=="spearman_adjusted"]

fig_corr <- ggplot(corr_summary[corr_name=="spearman_adjusted"], aes(x = msa_col, y = med, fill = cate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "gray",
                position = position_dodge(width = 0.15), width = 0.25) +
  scale_fill_manual(name = TeX("$\\rho_{Model}$"),
                    breaks = cate.info$breaks,
                    labels = cate.info$labels,
                    values = scales::alpha(cate.info$col, alpha = 0.8)) +
  scale_color_manual(name = TeX("$\\rho_{Model}$"),
                     breaks = cate.info$breaks,
                     labels = cate.info$labels,
                     values = cate.info$col) +
  # scale_fill_manual(name = TeX("$\\rho_{Model}$ category"),
  #                   breaks = c("lower", "higher"),
  #                   labels = c(TeX("Low $\\rho_{Model}$"), TeX("High $\\rho_{Model}$")),
  #                   values = scales::alpha(cate.info$col, alpha = 0.8)) +
  # scale_color_manual(name = TeX("$\\rho_{Model}$ category"),
  #                    breaks = c("lower", "higher"),
  #                    labels = c(TeX("Low $\\rho_{Model}$"), TeX("High $\\rho_{Model}$")),
  #                    values = cate.info$col) +
  labs(x = "MSA", y = TeX('Spearman $\\textit{r_s}$ (flow vs. distance)')) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        # axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.x = ggtext::element_markdown(angle = 30, hjust = 1, size = rel(1)),
        strip.text = element_text(face = "plain",size = 15),
        legend.position = "right")
#-----------------------------
# Figure: tij and dij for two categories
#-----------------------------
corr_summary <- corr[, .(
  mea = mean(corr_value),
  med = median(corr_value),
  lower_ci = quantile(corr_value, 0.25),
  upper_ci = quantile(corr_value, 0.75)
), by = .(cate,corr_name)]
fig_corr_cate <- ggplot(corr_summary[corr_name=="spearman_adjusted"], aes(x = cate, y = med, fill = cate)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), color = "gray",
                position = position_dodge(width = 0.25), width = 0.25) +
  scale_fill_manual(name = TeX("$\\rho_{Model}$"),
                    breaks = cate.info$breaks,
                    labels = cate.info$labels,
                    values = scales::alpha(cate.info$col, alpha = 0.8)) +
  scale_color_manual(name = TeX("$\\rho_{Model}$"),
                     breaks = cate.info$breaks,
                     labels = cate.info$labels,
                     values = cate.info$col) +
  labs(x = TeX("$\\rho_{Model}$"), y = TeX('Spearman $\\textit{r_s} (flow vs. distance)$')) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain",size=15),
        legend.position = "none")

fig_null_corr <- (fig_null/fig_corr) + plot_annotation(tag_levels = list(letters[3:4])) & theme(plot.tag = element_text(size = 20))
ggsave(fig_null_corr, filename = paste(figpath,"/msa/fig2.pdf",sep=""), width = 5*2, height = 4.5*2)



#============================================================
# Part 3. Edge and node coherence
#============================================================
#----------Plot the fitting results: spearman rank correlation and map----------#
Dname<<-Dname.set[1]
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif <- date.sep(Datevalue,durdate)
yindex <- 2
#----------data for plotting Pearson's correlation: data form sldr_fit.R----------#
dis.cor <- read.csv(paste(datapath,"/msa/SLDR_rank_cor_",Yname[yindex],".csv",sep=""), header=TRUE)
dis.cor$day <- as.Date(dis.cor$day)
dis.cor <- left_join(dis.cor,date.gif,by='day')
dis.cor.before <- subset(dis.cor,period==pervalue[1])
dis.cor.before <- left_join(dis.cor.before,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')

#----------data for plotting the map with max freq R2----------#
dis.r2 <- read.csv(paste(datapath,"/msa/SLDR_r2_",Yname[yindex],".csv",sep=""), header=TRUE)
dis.r2$day <- as.Date(dis.r2$day)
dis.r2 <- left_join(dis.r2,date.gif,by='day')
max.r2.day <- setDT(dis.r2)[, .SD[which.max(r2_exp)], by = .(msa, period), .SDcols = c("day", "r2_exp")]
max.r2 <- max.r2.day[, .N, by = .(period, day)][order(-N)][, .SD[1], by = period]
max.date <- c(as.Date(max.r2[period==pervalue[1],day]),as.Date(max.r2[period==pervalue[2],day]))
print(max.date[1])
max.day <- length(max.date)

#============================================================
# Part 3. Edge and node coherence
#============================================================
d <- 1
Dname <<- Dname.set[d]
source(paste(codepath, "/sldr_global_vars_funs.R", sep = ""))
date.gif <- date.sep(Datevalue, durdate)
#----------MSA BlockID, NBlock----------#
block.msa <- sf::read_sf(paste(geopath, "/msa/", Dname, ".geojson", sep = ""))
BlockID <- unique(block.msa$CensusBlockGroup)
NBlock <- length(BlockID)
#----------Center point of each block----------#
block.msa$centroid <- st_centroid(block.msa$geometry, of_largest_polygon = TRUE)
s1 <- st_coordinates(block.msa$centroid)
block.msa$lon <- s1[, 1]
block.msa$lat <- s1[, 2]
#----------Generate edges of network and W: 1 lag queen adjacency matrices----------#
W <- spdep::poly2nb(block.msa, queen = TRUE)
edges <- NULL
for (i in 1:length(W)) {
  if (length(W[[i]]) > 0) {
    edges <- plyr::rbind.fill(edges, data.frame(from = i, to = W[[i]]))
  }
}
dis.edges <- distinct(edges)
colnames(dis.edges) <- c("from", "to")
s1_from <- data.frame(from = 1:NBlock, Olon = block.msa$lon, Olat = block.msa$lat)
s1_to   <- data.frame(to   = 1:NBlock, Dlon = block.msa$lon, Dlat = block.msa$lat)
dis.edges <- left_join(dis.edges, s1_from, by = "from")
dis.edges <- left_join(dis.edges, s1_to,   by = "to")
#----------Data from sldr_fit.R: SLDR_fit----------#
sldr.fit <- read.csv(paste(datapath, "/msa/fit/", Dname, "_SLDR_fit_", Yname[yindex], ".csv", sep = ""),header = TRUE) %>% setDT()
sldr.fit$day <- as.Date(sldr.fit$day)
sldr.fit <- left_join(sldr.fit, date.gif, by = "day")
sldr.fit$period <- left_join(sldr.fit, data.frame(period = pervalue, ID = ID),by = "period")$ID
sldr.fit$res <- sldr.fit$empirical - sldr.fit$exp
sldr.fit <- left_join(sldr.fit, data.frame(day = Datevalue, day.id = factor(c(1:Day))),by = "day")
#----------Select one day----------#
land.fit <- subset(sldr.fit, day == max.date[1])
land.fit$lon <- block.msa$lon
land.fit$lat <- block.msa$lat
#----------Empirical flows: nodes----------#
df.nodes <- data.frame(
  lon = land.fit$lon,
  lat = land.fit$lat,
  mob = land.fit$empirical
)
df.nodes$mob <- df.nodes$mob / max(df.nodes$mob, na.rm = TRUE)

#========================================================#
#              Set local plotting extent here            #
#========================================================#
center_lon <- mean(block.msa$lon, na.rm = TRUE)
center_lat <- mean(block.msa$lat, na.rm = TRUE)
dx <- 0.022
dy <- 0.016
xmin <- center_lon - dx - 0.0012
xmax <- center_lon + dx - 0.002
ymin <- center_lat - dy + 0.0002
ymax <- center_lat + dy - 0.002
#----------Crop local polygon basemap----------#
bbox_local <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = st_crs(block.msa))
block.local <- st_crop(block.msa, bbox_local)

#----------Keep local nodes only----------#
df.nodes.local <- df.nodes %>% filter(lon >= xmin, lon <= xmax, lat >= ymin, lat <= ymax)
df.nodes.local[which.max(df.nodes.local$lat),]$mob<-0.015
df.nodes.local[which.min(df.nodes.local$lat),]$mob<-0.01
df.nodes.local[which.max(df.nodes.local$mob),]$mob <- 0.06
df.nodes.local[which.min(abs(df.nodes.local$lon+84.35745)),]$mob <- 0.1

#----------Keep local edges only----------#
dis.edges.local <- dis.edges %>%
  filter(
    Olon >= xmin, Olon <= xmax,
    Olat >= ymin, Olat <= ymax,
    Dlon >= xmin, Dlon <= xmax,
    Dlat >= ymin, Dlat <= ymax
  )

#-----------------------------
# plot
#-----------------------------
col_edge_base <- "#d25774"   # base colour for edge/node gradient
col_blue <- "#4b6aa8"        # nodes in panel a; edges in panel b
grad_cols <- c("#f7e8ee", "#d25774")
#----------Keep local edges only----------#
fig_map <- ggplot() +
  #----------Basemap polygons----------#
  geom_sf(data = block.local, fill = "white", color = "gray70",linewidth = 0.15) +
  #----------Network edges----------#
  # geom_segment(data = na.omit(dis.edges.local),
  #   aes(x = Olon, y = Olat, xend = Dlon, yend = Dlat),
  #   arrow = arrow(length = unit(0, "inches")),
  #   color = col_blue,
  #   linewidth = 1.2,
  #   linetype = "dashed",
  #   alpha = 1,
  #   lineend = "round"
  # ) +
  #----------Nodes----------#
  geom_point(data = df.nodes.local,
    aes(x = lon, y = lat, color = mob, size = mob),
    alpha = 1
  ) +
  scale_colour_gradientn(
    name = Dname,
    colors = grad_cols
    # colors = col.fit
  ) +
  scale_size_continuous(
    name = "Flows",
    range = c(2, 15),
    breaks = c(0.1, 0.5, 1.0),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  coord_sf(
    xlim = c(xmin, xmax),
    ylim = c(ymin, ymax),
    expand = FALSE
  ) +
  labs(title = "Node-level coherence") +
  theme_wy() +
  theme(
    panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5, linetype = "solid"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2, override.aes = list(size = 3))
  )

#-----------------------------
# 4. Combine two panels
#-----------------------------
# fig_edge_node <- (p1 / p2) + plot_annotation(tag_levels = list(letters[1:2])) & theme(plot.tag = element_text(size = 20))
ggsave(fig_map, filename = paste(figpath, "/msa/fig2b.pdf", sep = ""), width = 4, height = 3.5)








#============================================================
# Part 3. Edge and node coherence
#============================================================
#-----------------------------
# 0. Main colours
#-----------------------------
col_edge_base <- "#d25774"   # base colour for edge/node gradient
col_blue <- "#4b6aa8"        # nodes in panel a; edges in panel b

# low -> high gradient based on #d25774
grad_cols <- c("#f7e8ee", "#e9a9bc", "#d25774")

#-----------------------------
# 1. Build a star network
#-----------------------------
n_leaf <- 6
g <- make_star(n = n_leaf + 1, mode = "undirected", center = 1)
V(g)$name <- as.character(1:vcount(g))

angles <- seq(0, 2*pi, length.out = n_leaf + 1)[-1]
node_pos <- data.frame(
  name = as.character(1:(n_leaf + 1)),
  x = c(0, cos(angles) * 1.2),
  y = c(0, sin(angles) * 1.2)
)

#-----------------------------
# 2. Panel A: Edge-level decay
#-----------------------------
graph_edge <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  left_join(node_pos, by = "name") %>%
  activate(edges) %>%
  mutate(
    edge_r = c(0.85, 0.55, 0.20, 0.20, 0.20, 0.20),
    edge_width = c(5.5, 4.2, 1.6, 3.6, 2.0, 1.0)
  )

p1 <- ggraph(graph_edge, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(edge_colour = edge_r, edge_width = edge_width),
    lineend = "round",
    alpha = 0.95,
    show.legend = TRUE
  ) +
  scale_edge_width_identity() +
  scale_edge_colour_gradientn(
    colours = grad_cols,
    limits = c(0, 1),
    breaks = c(0.2, 0.4, 0.6, 0.8),
    labels = c("0.2", "0.4", "0.6", "0.8"),
    name = "Edge correlation"
  ) +
  geom_node_point(
    data = node_pos %>% filter(name != "1"),
    aes(x = x, y = y),
    size = 14,                  # slightly smaller than before
    shape = 21,
    alpha = 1,
    fill = col_blue,
    colour = "black",
    stroke = 1
  ) +
  geom_node_point(
    data = node_pos %>% filter(name == "1"),
    aes(x = x, y = y),
    size = 16,                  # slightly smaller than before
    shape = 21,
    alpha = 1,
    fill = col_blue,
    colour = "black",
    stroke = 1
  ) +
  labs(title = "Edge-level coherence") +
  coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

#-----------------------------
# 3. Panel B: Node-level coherence
#-----------------------------
node_df2 <- node_pos %>%
  mutate(
    node_group = c("Center", "Low", "High", "Very low", "Low", "Medium", "Low"),
    node_r = c(0.50, 0.30, 0.80, 0.10, 0.30, 0.60, 0.30),  # center also gets a value
    size_val = c(20, 12, 24, 13, 12, 21, 12)
  )

graph_node <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  left_join(node_df2, by = "name")

p2 <- ggraph(graph_node, layout = "manual", x = x, y = y) +
  geom_edge_link(
    colour = col_blue,
    width = 1.2,
    linetype = "dashed",
    alpha = 0.8,
    lineend = "round"
  ) +
  geom_node_point(
    aes(x = x, y = y, size = size_val, fill = node_r),
    shape = 21,
    colour = "black",
    stroke = 1
  ) +
  scale_size_identity() +
  scale_fill_gradientn(
    colours = grad_cols,
    limits = c(0, 1),
    breaks = c(0.1, 0.3, 0.5, 0.6, 0.8),
    labels = c("0.1", "0.3", "0.5", "0.6", "0.8"),
    name = "Node correlation"
  ) +
  labs(title = "Node-level coherence") +
  coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
  theme_void() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

#-----------------------------
# 4. Combine two panels
#-----------------------------
fig_edge_node <- (p1 / p2) + plot_annotation(tag_levels = list(letters[1:2])) & theme(plot.tag = element_text(size = 20))
fig_edge_node
# save
ggsave(plot = fig_edge_node, filename = paste(figpath, "/msa/fig2_frame.pdf", sep = ""), width = 4, height = 8)





