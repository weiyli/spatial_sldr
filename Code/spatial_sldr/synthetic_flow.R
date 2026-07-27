# rm(list = ls())

# Null flow models (gravity model and radiation model) -> synthetic rho
# Correlation between flow (Tij) and distance (dij)

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'

#----------Load packages----------#
library(sf)          # read_sf() 
library(spdep)       # poly2nb() 
library(igraph)      # diameter()
library(geosphere)

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

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'

#----------Load functions----------#
#----1. Distance matrix D (meters)----#
compute_D <- function(tess) {
  coords <- sf::st_coordinates(sf::st_centroid(tess$geometry))
  D <- geosphere::distm(coords, fun = geosphere::distHaversine)
  diag(D) <- Inf
  D
}

#----2. Population-weighted center (lon, lat)----#
compute_center <- function(tess, P) {
  coords <- sf::st_coordinates(sf::st_centroid(tess$geometry))
  w <- P
  w[is.na(w) | w < 0] <- 0
  if (sum(w) == 0) w <- rep(1, length(w))
  c(sum(coords[,1]*w)/sum(w), sum(coords[,2]*w)/sum(w))
}

#----3. Safe numeric cleanup: replace non-finite with 0, diag=0----#
sanitize_kernel <- function(K) {
  K[!is.finite(K)] <- 0
  diag(K) <- 0
  K
}

#----Global scaling: sum(K_scale)=T_total_observed----#
scale_global <- function(K, T_total_observed) {
  denom <- sum(K)
  if (!is.finite(denom) || denom <= 0) stop("Global scaling failed: sum(K)<=0.")
  S <- T_total_observed / denom
  K_scale <- S * K
  diag(K_scale) <- 0
  K_scale 
}

#----Optional: Origin constrained scaling: rowSums(K_scale)=TO----#
scale_origin <- function(K, TO) {
  row_den <- rowSums(K)
  row_den[row_den <= 0] <- 1
  K_scale <- K / row_den
  K_scale <- K_scale * TO   # recycle along columns
  diag(K_scale) <- 0
  K_scale 
}

#----Gravity kernel: K_ij = (P_i^alpha * P_j^gamma) / D_ij^beta----#
gravity_kernel <- function(P, D, alpha=1, beta=1, gamma=1, eps=1e-6) {
  N <- length(P)
  Pi <- matrix(P^alpha, nrow=N, ncol=N, byrow=FALSE)
  Pj <- matrix(P^gamma, nrow=N, ncol=N, byrow=TRUE)
  K <- (Pi * Pj) / ((D + eps)^beta)
  sanitize_kernel(K)
}

#----Radiation kernel (parameter-free)----#
radiation_kernel <- function(P, D) {
  N <- length(P)
  K <- matrix(0, nrow=N, ncol=N)
  
  diag(D) <- Inf
  
  for (i in 1:N) {
    ord <- order(D[i, ], na.last = TRUE)
    ord <- ord[ord != i]
    
    mi <- P[i]
    mj <- P[ord]
    
    cumP <- cumsum(P[ord])
    s_vec <- c(0, cumP[-length(cumP)])  # strictly closer than j
    
    denom1 <- (mi + s_vec)
    denom2 <- (mi + mj + s_vec)
    
    kij <- (mi * mj) / (denom1 * denom2)
    kij[!is.finite(kij)] <- 0
    K[i, ord] <- kij
  }
  
  sanitize_kernel(K)
}

#----Radial kernel (extreme geography null)----#
radial_kernel <- function(P, tess, radial_type=c("power","exp"), eps=1) {
  radial_type <- match.arg(radial_type)
  N <- length(P)
  
  coords <- sf::st_coordinates(sf::st_centroid(tess$geometry))
  center <- compute_center(tess, P)
  
  rr <- geosphere::distHaversine(coords, matrix(center, nrow=N, ncol=2, byrow=TRUE))
  rr[!is.finite(rr)] <- 0
  
  if (radial_type == "power") {
    A <- P / (rr + eps)
  } else {
    eta <- median(rr[rr > 0], na.rm = TRUE)
    if (!is.finite(eta) || eta <= 0) eta <- 1
    A <- P * exp(-rr / eta)
  }
  A[!is.finite(A)] <- 0
  
  K <- tcrossprod(A, A)   # K_ij = A_i * A_j
  sanitize_kernel(K)
}


#----Main: generate null mobility vector from specified model----#
generate_null_mob <- function(model,
                              P,
                              tess,
                              T_total_observed,
                              radial_type="exp",
                              gravity_par=list(alpha=1,beta=1,gamma=1),
                              drop_na_pop=TRUE) {
  
  # handle NA pop
  idx_keep <- rep(TRUE, length(P))
  if (drop_na_pop) {
    idx_keep <- is.finite(P) & !is.na(P) & (P > 0)
    # keep only valid pop units
    tess <- tess[idx_keep, ]
    P <- P[idx_keep]
  } else {
    P[is.na(P) | !is.finite(P)] <- 0
  }
  
  D <- compute_D(tess)
  
  # choose kernel
  if (model == "gravity") {
    K <- gravity_kernel(P, D,
                        alpha=gravity_par$alpha,
                        beta=gravity_par$beta,
                        gamma=gravity_par$gamma)
  } else if (model == "radiation") {
    K <- radiation_kernel(P, D)
  } else if (model == "radial") {
    K <- radial_kernel(P, tess, radial_type=radial_type)
  } else {
    stop("model must be one of: 'gravity','radiation','radial'")
  }
  
  # global scaling to match total observed flow (your current choice)
  K_scale <- scale_global(K, T_total_observed)
  
  mob <- rowSums(K_scale)
  
  # if dropped, return with original length + mapping (optional)
  if (drop_na_pop) {
    out <- rep(0, length(idx_keep))
    out[idx_keep] <- mob
    return(out)
  } else {
    return(mob)
  }
}




#----------SLDR fitting for synthetic flow----------#
# models <- c("gravity","radiation","radial")
models <- c("gravity","radial")
for(s in 2:Nregion){
  for(yindex in 2:2){
    # for(d in 8:8){
    for(d in 1:Nmsa){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # date.gif<-date.sep(Datevalue,durdate)
      
      #----------1.pop and total flow----------#
      #----------BlockID, NBlock----------#
      block.msa <- sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID <- unique(block.msa$CensusBlockGroup)
      NBlock <- length(BlockID)
      #----------W: 1 to lag.max lag queen adjacency matrices----------#
      Wd <- spdep::poly2nb(block.msa, queen = T) 
      W1 <- nb.mat(nb=Wd)   # W1 <- spdep::nb2mat(Wd, style = "B") 
      lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed")%>%igraph::diameter()
      W <<- queen.adj(W1, lag.max)
      slag <<- seq(1,lag.max,1)
      
      Mob_list <- lapply(models, function(x) matrix(0, nrow=NBlock, ncol=Day))
      names(Mob_list) <- models
      pop <- block.msa$pop
      
      for (i in 1:Day) {
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        T_total_observed <- sum(s2$mob, na.rm=TRUE)
        
        Mob_list[["gravity"]][, i] <- generate_null_mob(
          model="gravity",
          P=pop,
          tess=block.msa,
          T_total_observed=T_total_observed,
          gravity_par=list(alpha=1,beta=1,gamma=1),
          drop_na_pop=TRUE
        )
        
        # Mob_list[["radiation"]][, i] <- generate_null_mob(
        #   model="radiation",
        #   P=pop,
        #   tess=block.msa,
        #   T_total_observed=T_total_observed,
        #   drop_na_pop=TRUE
        # )
        
        Mob_list[["radial"]][, i] <- generate_null_mob(
          model="radial",
          P=pop,
          tess=block.msa,
          T_total_observed=T_total_observed,
          radial_type="exp",   # or "power"
          drop_na_pop=TRUE
        )
      }
      
      
      fit_one_model <- function(Mob_mat, model_name) {
        Mob <<- Mob_mat
        exp_list <- lapply(1:Day, function(g) exp.day.est(g))
        out <- do.call(rbind, exp_list) %>% as.data.frame()
        out$index <- "exp"
        out$model <- model_name
        out
      }
      
      lower_bound<<-0
      upper_bound<<-5
      es_name <- "rmse_rmspe"
      error_test <<- error_specs[[es_name]]$fun
      
      res_all <- do.call(rbind, lapply(names(Mob_list), function(m) {
        fit_one_model(Mob_list[[m]], m)
      }))
      
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=res_all[,1], error=res_all[,2], index=res_all[,3], model=res_all[,4])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/param/",Dname,"_SLDR_params_synthetic_",Yname[yindex],".csv",sep=""), row.names = FALSE)
      
      print(d)
 
    } # dis
  }# yindex
}# region





#----------Part2: correlation between Tij and dij----------#
flowpath <- 'D:/ood/Data/Flow/daily_flow'
datapath <- 'D:/ood/Data/spatial_sldr'

compute_OD <- function(tess) {
  coords <- sf::st_coordinates(sf::st_centroid(tess$geometry))
  D <- geosphere::distm(coords, fun = geosphere::distHaversine)
  return(D)
}

flow_dist_corr <- NULL
for(d in 1:Nmsa){
  
  # Read the MSA geographical data
  Dname <<- Dname.set[d]
  block.msa <- sf::read_sf(paste(geopath, "/msa/", Dname, ".geojson", sep=""))
  BlockID <- unique(block.msa$CensusBlockGroup)
  NBlock <- length(BlockID)
  pairs <- choose(NBlock, 2)
  pop <- block.msa$pop/1e6
  # Compute the distance matrix dij using Haversine formula: m->km
  D <- compute_OD(block.msa)/1000  

  for(i in 1:Day){
    
    # Filter the OD flow data to only include blocks in the current MSA
    df <- fread(paste0(flowpath, "/", Dname, "/Daily_Flow_", Datevalue[i], ".csv"))
    OD.FTH <- df[(orig_bg != dest_bg) & (no_visits>2)]
    # OD.FTH <- df[no_visits>2]
    # OD.FTH <- df[orig_bg != dest_bg]
    
    OD_pairs <- unique(OD.FTH[, c("orig_bg", "dest_bg")])
    pairs_ratio <- length(OD_pairs)/pairs

    # Extract Tij (flow) and dij (distance) for each pair of blocks
    Tij <- OD.FTH$no_visits 
    orig_bg <- match(OD.FTH$orig_bg, BlockID)
    dest_bg <- match(OD.FTH$dest_bg, BlockID)
    dij <- D[cbind(orig_bg, dest_bg)]  
    
    # Remove the population effect (mi, mj) from Tij     
    valid_rows <- !(pop[orig_bg] == 0 | pop[dest_bg] == 0 | is.na(pop[orig_bg])| is.na(pop[dest_bg]))
    Tij_valid <- Tij[valid_rows]
    orig_bg_valid <- orig_bg[valid_rows]
    dest_bg_valid <- dest_bg[valid_rows]
    Tij_adjusted <- Tij_valid / (pop[orig_bg_valid] * pop[dest_bg_valid])  # Remove population effect
    dij_adjusted <- D[cbind(orig_bg_valid, dest_bg_valid)]  
    
    # Calculate the correlation between Tij and dij 
    pearson_cor <- cor(Tij, dij, method = "pearson")  
    spearman_cor <- cor(Tij, dij, method = "spearman")
    # Calculate the correlation between adjusted Tij and dij
    pearson_cor_adjusted <- cor(Tij_adjusted, dij_adjusted, method = "pearson")
    spearman_cor_adjusted <- cor(Tij_adjusted, dij_adjusted, method = "spearman")
    
    # Calculate the correlation between log(Tij) and log(dij)
    log_Tij <- log(Tij)
    log_dij <- log(dij)
    log_pearson_cor <- cor(log_Tij, log_dij, method = "pearson")  
    log_spearman_cor <- cor(log_Tij, log_dij, method = "spearman")
    
    log_Tij_adjusted <- log(Tij_adjusted)
    log_dij_adjusted <- log(dij_adjusted)
    log_pearson_cor_adjusted <- cor(log_Tij_adjusted, log_dij_adjusted, method = "pearson")  
    log_spearman_cor_adjusted <- cor(log_Tij_adjusted, log_dij_adjusted , method = "spearman")
    
    # Store the correlation results in a data frame
    corr <- data.frame(
      msa = Dname,
      day = Datevalue[i],
      pairs_ratio = pairs_ratio,
      corr_value = c(pearson_cor, spearman_cor, pearson_cor_adjusted, spearman_cor_adjusted,
                     log_pearson_cor, log_spearman_cor, log_pearson_cor_adjusted, log_spearman_cor_adjusted),
      corr_name = c('pearson', 'spearman', 'pearson_adjusted', 'spearman_adjusted', 
                    'log_pearson', 'log_spearman', 'log_pearson_adjusted', 'log_spearman_adjusted')
    )
    # corr <- data.frame(
    #   msa = Dname,
    #   day = Datevalue[i],
    #   pairs_ratio = pairs_ratio,
    #   corr_value = c(pearson_cor,spearman_cor,log_pearson_cor,log_spearman_cor),
    #   corr_name = c('pearson', 'spearman','log_pearson','log_spearman')
    # )

    flow_dist_corr <- plyr::rbind.fill(flow_dist_corr, corr)
    
  } 
  print(d)  
}
write.csv(flow_dist_corr, file = paste(datapath,"/msa/flow_dist_corr.csv",sep=""), row.names = FALSE)

  
  








