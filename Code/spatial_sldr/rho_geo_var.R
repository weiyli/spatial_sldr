# rm(list = ls())

# Urban spatial structure

# polycentricity and compactness

# We focus on two key dimensions of urban spatial structure: the monocentricity-polycentricity dimension and the concentration-dispersion dimension (Anas et al., 1998).

# The rank-size distribution as an indicator for mono/polycentricity.


#---------- Work paths ----------#
setwd("D:/ood/")
codepath <- "D:/ood/Code/spatial_sldr/spatial_sldr"
geopath  <- "D:/ood/Data/Geo"
flowpath <- "D:/ood/Data/Flow"
datapath <- "D:/ood/Data/spatial_sldr"


#---------- Load packages ----------#
library(sf)        # st_read(), st_transform(), geometry operations
library(spdep)     # poly2nb() etc. (used in SLDR pipeline)
library(igraph)    # network metrics (e.g., diameter)
library(geosphere) # distm(), distHaversine


#---------- Load MSAs ----------#
Dname.msa <- c("Atlanta",
               "Boston",
               "Chicago",
               "Dallas",
               "Houston",
               "Los Angeles",
               "Miami",
               "New York",
               "Philadelphia",
               "San Francisco",
               "Seattle",
               "Washington, D.C.")
Dname.dis <- c("Houston",
               "Jacksonville",
               "Houston",
               "Santa Rosa-Petaluma")
Dname.set <- c(Dname.msa, Dname.dis)

#---------- Load global SLDR vars / functions ----------#
d <<- 1
source(paste0(codepath, "/sldr_global_vars_funs.R"))

Flow.name <- "Intra_Flow"

#---------- Set the global function pointer used by homo/power/exp.param.est ----------#
# error_specs is assumed to be loaded from sldr_global_vars_funs.R
error_test <<- error_specs[["rmse_rmspe"]]$fun


#======================================================================
# Helper functions for geography/density controls
# Key design:
# - sf_msa is EPSG:4326 (lon/lat) and already contains BOTH pop and area_m2
# - Radial distance uses Haversine (meters), so no projection is needed
# - Density uses area_m2 attribute (not st_area on lon/lat)
# - BUT: polycentricity clustering via st_buffer(1500) requires a metric CRS,
#        so we project ONLY for that step (recommended: EPSG:3857 or 5070)
#======================================================================

#----1. Population-weighted center (lon/lat; returns c(lon,lat))----#
pop_weighted_center <- function(sf_obj, pop_col = "pop") {
  cen <- sf::st_centroid(sf::st_geometry(sf_obj))
  coords <- sf::st_coordinates(cen)  # X=lon, Y=lat
  
  w <- sf_obj[[pop_col]]
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) == 0) w <- rep(1, length(w))
  
  c(sum(coords[, 1] * w) / sum(w),
    sum(coords[, 2] * w) / sum(w))
}

#----2. Gini coefficient (lower = more even)----#
gini <- function(x) {
  x <- x[is.finite(x) & x >= 0]
  if (length(x) <= 1) return(NA_real_)
  if (sum(x) == 0) return(0)
  x <- sort(x)
  n <- length(x)
  (2 * sum(x * seq_len(n)) / (n * sum(x))) - (n + 1) / n
}

#----3. Radial gradient: log(density) ~ distance_to_center (meters)----#
# Uses:
# - density = pop / area_m2 (attribute)
# - distance = Haversine on lon/lat centroids (meters)
radial_gradient <- function(sf_obj, pop_col = "pop", area_col = "area_m2") {
  if (!(pop_col %in% names(sf_obj))) stop("Missing pop column: ", pop_col)
  if (!(area_col %in% names(sf_obj))) stop("Missing area column: ", area_col)
  
  # density (persons/m^2)
  P <- sf_obj[[pop_col]]
  A_m2 <- sf_obj[[area_col]]
  A_m2[!is.finite(A_m2) | A_m2 <= 0] <- NA_real_
  dens <- P / A_m2
  dens[!is.finite(dens) | dens <= 0] <- NA_real_
  
  # centroid coords (lon/lat)
  cen <- sf::st_centroid(sf::st_geometry(sf_obj))
  coords <- sf::st_coordinates(cen)  # X=lon, Y=lat
  
  # population-weighted center (lon/lat)
  ctr <- pop_weighted_center(sf_obj, pop_col = pop_col)
  
  # haversine distance to center (meters)
  rr <- geosphere::distHaversine(coords, matrix(ctr, nrow = nrow(coords), ncol = 2, byrow = TRUE))
  rr[!is.finite(rr)] <- NA_real_
  
  # fit log(density) ~ distance
  dt <- data.table(d = rr, logdens = log(dens))
  dt <- dt[is.finite(d) & is.finite(logdens)]
  if (nrow(dt) < 30) return(list(slope = NA_real_, r2 = NA_real_))
  
  fit <- lm(logdens ~ d, data = dt)
  list(slope = coef(fit)[2], r2 = summary(fit)$r.squared)
}

#----4. Polsby-Popper compactness: 4*pi*A / P^2----#
# IMPORTANT: compactness needs planar metric CRS (not lon/lat)
compactness_pp <- function(poly_sf_proj) {
  poly_u <- sf::st_union(sf::st_make_valid(poly_sf_proj))
  A <- as.numeric(sf::st_area(poly_u))
  P <- as.numeric(sf::st_length(sf::st_cast(poly_u, "MULTILINESTRING")))
  if (!is.finite(A) || !is.finite(P) || P == 0) return(NA_real_)
  (4 * pi * A) / (P^2)
}



#----------Part 1: msa----------#
#======================================================================
# Extract city-level geography/morphology metrics for a single MSA
# Input sf_msa:
# - EPSG:4326 lon/lat
# - has columns: pop, area_m2
#======================================================================
extract_msa_metrics <- function(sf_msa, msa_name,
                                pop_col = "pop",
                                area_col = "area_m2",
                                proj_crs = 3857) {
  
  #-----------------------------
  # Basic totals (use area_m2 attribute)
  #-----------------------------
  if (!(pop_col %in% names(sf_msa))) stop("Missing pop column: ", pop_col)
  if (!(area_col %in% names(sf_msa))) stop("Missing area column: ", area_col)
  
  pop_total  <- sum(sf_msa[[pop_col]], na.rm = TRUE)
  area_total <- sum(sf_msa[[area_col]], na.rm = TRUE)  # m^2 (from attribute)
  avg_density <- pop_total / area_total                # persons / m^2
  
  #-----------------------------
  # CBG area statistics
  #-----------------------------
  cbg_area_km2 <- sf_msa[[area_col]] / 1e6
  mean_cbg_area_km2   <- mean(cbg_area_km2, na.rm = TRUE)
  median_cbg_area_km2 <- median(cbg_area_km2, na.rm = TRUE)
  n_cbg <- nrow(sf_msa)
  
  #-----------------------------
  # Radial density gradient (no projection; Haversine distance in meters)
  #-----------------------------
  rg <- radial_gradient(sf_msa, pop_col = pop_col, area_col = area_col)
  
  #-----------------------------
  # For compactness + polycentric clustering, we need a metric CRS
  #-----------------------------
  sf_msa_proj <- sf_msa
  if (sf::st_is_longlat(sf_msa_proj)) {
    sf_msa_proj <- sf::st_transform(sf_msa_proj, proj_crs)
  }
  
  # metro boundary polygon (projected)
  poly_proj <- sf::st_as_sf(sf::st_union(sf::st_make_valid(sf_msa_proj)))
  
  # compactness
  comp_pp <- compactness_pp(poly_proj)
  
  #-----------------------------
  # Morphological polycentricity (simple proxy)
  # - density computed from pop/area_m2 attribute (not st_area)
  # - cluster high-density units using buffer in projected CRS
  #-----------------------------
  dens_attr <- sf_msa[[pop_col]] / sf_msa[[area_col]]
  thr <- quantile(dens_attr, 0.95, na.rm = TRUE)  # top 5% density as "centers"
  centers_idx <- which(is.finite(dens_attr) & dens_attr >= thr)
  centers_proj <- sf_msa_proj[centers_idx, ]
  
  if (nrow(centers_proj) >= 5) {
    # buffer 1.5 km in projected CRS
    buf <- sf::st_buffer(sf::st_make_valid(centers_proj), dist = 1000)
    clus_geom <- sf::st_union(buf)
    clus_poly <- sf::st_as_sf(sf::st_cast(clus_geom, "POLYGON"))
    clus_poly$id <- seq_len(nrow(clus_poly))
    
    # intersect projected tessellation with clusters
    inter <- sf::st_intersection(sf_msa_proj, clus_poly["id"])
    
    center_sizes <- inter %>%
      sf::st_drop_geometry() %>%
      group_by(id) %>%
      summarise(size = sum(.data[[pop_col]], na.rm = TRUE), .groups = "drop") %>%
      pull(size)
    
    poly_gini <- gini(center_sizes)
    
    rank_slope_k <- function(x, k) {
      x <- sort(x, decreasing = TRUE)
      if (length(x) < k) return(NA_real_)
      x <- x[1:k]
      r <- 1:k
      abs(coef(lm(log(x) ~ log(r)))[2])
    }
    
    rs <- mean(c(rank_slope_k(center_sizes, 2),
                 rank_slope_k(center_sizes, 3),
                 rank_slope_k(center_sizes, 4)), na.rm = TRUE)
  } else {
    poly_gini <- NA_real_
    rs <- NA_real_
  }
  
  #-----------------------------
  # Output
  #-----------------------------
  data.table(
    msa = msa_name,
    pop_total = pop_total,
    area_km2 = area_total / 1e6,
    avg_density = avg_density * 1e6,  # persons per km^2
    mean_cbg_area_km2 = mean_cbg_area_km2,
    median_cbg_area_km2 = median_cbg_area_km2,
    n_cbg = n_cbg,
    
    radial_slope = rg$slope,
    radial_r2 = rg$r2,
    compactness = comp_pp,
    poly_gini = poly_gini,
    rank_size_slope = rs
  )
}



#======================================================================
# Batch compute geography/morphology metrics via a for-loop
#======================================================================
metrics_list <- vector("list", Nmsa)
for (d in 1:Nmsa) {
  
  #---------- Select MSA ----------#
  Dname <<- Dname.set[d]
  source(paste(codepath, "/sldr_global_vars_funs.R", sep=""))
  
  #---------- Read MSA geojson ----------#
  block.msa <- sf::read_sf(paste0(geopath, "/msa/", Dname, ".geojson"))
  
  pop_col  <- "pop"
  area_col <- "area_m2"
  
  if (!(pop_col %in% names(block.msa))) {
    message("[Skip] Missing pop column '", pop_col, "' in: ", Dname)
    next
  }
  if (!(area_col %in% names(block.msa))) {
    message("[Skip] Missing area column '", area_col, "' in: ", Dname)
    next
  }
  
  #---------- Extract metrics ----------#
  metrics_list[[d]] <- extract_msa_metrics(
    sf_msa   = block.msa,
    msa_name = Dname,
    pop_col  = pop_col,
    area_col = area_col,
    proj_crs = 3857   # you can change to 5070 for CONUS equal-area
  )
  
  message("[Done] ", Dname)
}

msa_metrics <- data.table::rbindlist(metrics_list, fill = TRUE)
write.csv(msa_metrics, file = paste(datapath,"/msa/geo_var.csv",sep=""), row.names = FALSE)







#----------Part 2: county_msa----------#
extract_county_metrics <- function(sf_county, msa_name, county_id,
                                   pop_col = "pop",
                                   area_col = "area_m2",
                                   proj_crs = 3857) {
  
  if (!(pop_col %in% names(sf_county))) stop("Missing pop column: ", pop_col)
  if (!(area_col %in% names(sf_county))) stop("Missing area column: ", area_col)
  
  pop_total  <- sum(sf_county[[pop_col]], na.rm = TRUE)
  area_total <- sum(sf_county[[area_col]], na.rm = TRUE)  # m^2
  avg_density <- pop_total / area_total                   # persons/m^2
  
  #-----------------------------
  # CBG area statistics within county
  #-----------------------------
  cbg_area_km2 <- sf_county[[area_col]] / 1e6
  mean_cbg_area_km2   <- mean(cbg_area_km2, na.rm = TRUE)
  median_cbg_area_km2 <- median(cbg_area_km2, na.rm = TRUE)
  n_cbg               <- nrow(sf_county)
  
  
  rg <- radial_gradient(sf_county, pop_col = pop_col, area_col = area_col)
  
  sf_proj <- sf_county
  if (sf::st_is_longlat(sf_proj)) sf_proj <- sf::st_transform(sf_proj, proj_crs)
  
  poly_proj <- sf::st_as_sf(sf::st_union(sf::st_make_valid(sf_proj)))
  comp_pp <- compactness_pp(poly_proj)
  
  dens_attr <- sf_county[[pop_col]] / sf_county[[area_col]]
  thr <- quantile(dens_attr, 0.95, na.rm = TRUE)
  centers_idx <- which(is.finite(dens_attr) & dens_attr >= thr)
  centers_proj <- sf_proj[centers_idx, ]
  
  if (nrow(centers_proj) >= 5) {
    buf <- sf::st_buffer(sf::st_make_valid(centers_proj), dist = 1000)  # 1 km
    clus_geom <- sf::st_union(buf)
    clus_poly <- sf::st_as_sf(sf::st_cast(clus_geom, "POLYGON"))
    clus_poly$id <- seq_len(nrow(clus_poly))
    
    inter <- sf::st_intersection(sf_proj, clus_poly["id"])
    
    center_sizes <- inter %>%
      sf::st_drop_geometry() %>%
      group_by(id) %>%
      summarise(size = sum(.data[[pop_col]], na.rm = TRUE), .groups = "drop") %>%
      pull(size)
    
    poly_gini <- gini(center_sizes)
    
    rank_slope_k <- function(x, k) {
      x <- sort(x, decreasing = TRUE)
      if (length(x) < k) return(NA_real_)
      x <- x[1:k]
      r <- 1:k
      abs(coef(lm(log(x) ~ log(r)))[2])
    }
    
    rs <- mean(c(rank_slope_k(center_sizes, 2),
                 rank_slope_k(center_sizes, 3),
                 rank_slope_k(center_sizes, 4)), na.rm = TRUE)
  } else {
    poly_gini <- NA_real_
    rs <- NA_real_
  }
  
  data.table(
    msa = msa_name,
    pop_total = pop_total,
    area_km2 = area_total / 1e6,
    avg_density = avg_density * 1e6,  # persons per km^2
    
    mean_cbg_area_km2 = mean_cbg_area_km2,
    median_cbg_area_km2 = median_cbg_area_km2,
    n_cbg = n_cbg,
    
    radial_slope = rg$slope,
    radial_r2 = rg$r2,
    compactness = comp_pp,
    poly_gini = poly_gini,
    rank_size_slope = rs,
    county = county_id
  )
}

#----------Load msa-county / msa-tract relationship tables----------#
county <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county$county_fips <- sprintf("%05d", as.numeric(county$county_fips))

#======================================================================
# Batch compute COUNTY-level metrics for all MSA-counties
#======================================================================
metrics_list <- list()
k <- 0
for (d in 1:Nmsa) {
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  county.msa <- subset(county,msa_id==d)
  Ncounty <- nrow(county.msa)
  
  block.msa.all <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.msa.all <- block.msa.all %>% mutate(pop = ifelse(is.na(pop), 0, pop))
  
  for(f in 1:Ncounty){
    
    #----------BlockID, NBlock----------#
    county.id <- county.msa$county_fips[f]
    block.county <- block.msa.all[which(block.msa.all$county_fips==county.msa$county_fips[f]),]
    
    # align area column name with your files:
    pop_col <- "pop"
    if ("area_m2" %in% names(block.county)) {
      area_col <- "area_m2"
    } else if ("area" %in% names(block.county)) {
      # your rho-fitting script uses area <- block.msa$area
      area_col <- "area"
      # if area is km^2 or m^2 is ambiguous, you should standardize before using!
      # Here we assume area is m^2; change if needed.
    } else {
      message("[Skip] Missing area column (area_m2/area) in: ", Dname, " ", county.id)
      next
    }
    
    if (!(pop_col %in% names(block.county))) {
      message("[Skip] Missing pop column in: ", Dname, " ", county.id)
      next
    }
    
    k <- k + 1
    metrics_list[[k]] <- extract_county_metrics(
      sf_county = block.county,
      msa_name  = Dname,
      county_id = county.id,
      pop_col   = pop_col,
      area_col  = area_col,
      proj_crs  = 3857
    )
    
    message("[Done] ", Dname, " / county ", county.id)
  }
}

county_metrics <- data.table::rbindlist(metrics_list, fill = TRUE)
write.csv(county_metrics, file = paste(datapath,"/county/geo_var.csv",sep=""), row.names = FALSE)




