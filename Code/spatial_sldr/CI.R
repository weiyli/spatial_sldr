# rm(list = ls())

# SLDR fitting and higer-order queen weight matrix at the msa level for rho estimation


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'


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

#----------county in each msa----------#
county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'



#----------Bootstrap settings (robustness / CI)----------#
B <- 100                 # bootstrap 
keep_ratio <- 0.90       # keep 90% CBGs
set.seed(1)

#----------Helper: fit one day for three kernels under current globals----------#
fit_one_day_allkernels <- function(g){
  out_h <- homo.day.est(g)   # returns c(rho, error) or similar
  out_p <- power.day.est(g)
  out_e <- exp.day.est(g)
  rbind(
    data.frame(day = Datevalue[g], index = "homo",  rho = out_h[1], error = out_h[2]),
    data.frame(day = Datevalue[g], index = "power", rho = out_p[1], error = out_p[2]),
    data.frame(day = Datevalue[g], index = "exp",   rho = out_e[1], error = out_e[2])
  )
}


#----------SLDR fitting + CI----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in 8:8){
      
      Dname <<- Dname.set[d]
      source(paste0(codepath, "/sldr_global_vars_funs.R"))
      
      #------------------------------------------------------------
      # Read spatial tessellation (CBG polygons) for a given MSA
      #------------------------------------------------------------
      block_full <- sf::read_sf(paste0(geopath, "/", region[s], "/", Dname, ".geojson"))
      
      # IMPORTANT: Fix the ordering of spatial units to guarantee that
      # the adjacency matrix and mobility vectors are consistently aligned
      # across point estimation and bootstrap replicates.
      block_full <- block_full[order(block_full$CensusBlockGroup), ]
      BlockID_full <- unique(block_full$CensusBlockGroup)
      BlockID_full <- as.character(BlockID_full)
      
      #------------------------------------------------------------
      # Pre-load daily mobility data (avoid repeated disk I/O)
      #------------------------------------------------------------
      flow_list <- vector("list", Day)
      for(i in 1:Day){
        tmp <- read.csv(
          paste0(flowpath, "/", region[s], "/", Dname, "/Intra_Flow_", Datevalue[i], ".csv"),
          header = TRUE
        )
        # Ensure consistent ID type for fast matching
        tmp$CensusBlockGroup <- as.character(tmp$CensusBlockGroup)
        flow_list[[i]] <- tmp
      }
      
      #------------------------------------------------------------
      # Compute queen contiguity adjacency ONCE on the full sample
      #------------------------------------------------------------
      # Rationale: geometric contiguity (poly2nb) is computationally expensive.
      # For delete-p bootstrap over spatial units, we can form the induced subgraph
      # by subsetting the full adjacency matrix, which is equivalent to recomputing
      # contiguity after node deletion (node removal cannot create new edges).
      Wd_full <- spdep::poly2nb(block_full, queen = TRUE)
      W1_full <- nb.mat(nb = Wd_full)
      rm(Wd_full)
      
      #------------------------------------------------------------
      # Externalize higher-order queen adjacency construction
      #------------------------------------------------------------
      # Build the full higher-order adjacency list only once:
      # W_full_list[[l]] encodes the l-th order neighborhood structure
      # derived from the full contiguity graph.
      g_full <- igraph::graph_from_adjacency_matrix(W1_full, mode = "undirected", diag = FALSE)
      lag.max_full <- igraph::diameter(g_full)
      
      W_full_list <- queen.adj(W1_full, lag.max_full)  # computed once per MSA
      slag_full <- seq_len(lag.max_full)
      
      #------------------------------------------------------------
      # Helper: build global variables given a retained set of CBG IDs
      # (delete-p bootstrap via induced subgraph of precomputed W_full_list)
      #------------------------------------------------------------
      build_globals_from_keep <- function(keep_ids){
        
        keep_ids <- as.character(keep_ids)
        
        # Map retained IDs to indices of the full ordering and keep them sorted.
        keep_idx <- match(keep_ids, BlockID_full)
        keep_idx <- keep_idx[!is.na(keep_idx)]
        keep_idx <- sort(keep_idx)
        
        # Corresponding spatial unit IDs (aligned with W / Mob_mat)
        BlockID <- BlockID_full[keep_idx]
        NBlock_local <- length(BlockID)
        
        # Induced subgraph for the 1st-order adjacency (useful for diagnostics)
        W1_sub <- W1_full[keep_idx, keep_idx, drop = FALSE]
        
        # (Optional) diameter on the induced subgraph (for reporting / diagnostics)
        # NOTE: we keep the kernel definition anchored to the full lag structure,
        # and represent node deletion by subsetting the precomputed higher-order matrices.
        g_sub <- igraph::graph_from_adjacency_matrix(W1_sub, mode = "undirected", diag = FALSE)
        lag.max_local <- igraph::diameter(g_sub)
        
        #------------------------------------------------------------
        # Induced subgraph of the FULL higher-order adjacency list
        #------------------------------------------------------------
        # Subset each lag matrix to the retained nodes. This implements structural
        # perturbation by node deletion while keeping the lag definition fixed.
        W <<- lapply(W_full_list, function(Wl) Wl[keep_idx, keep_idx, drop = FALSE])
        slag <<- slag_full
        
        #------------------------------------------------------------
        # Build the mobility matrix (rows aligned with the adjacency matrix)
        #------------------------------------------------------------
        Mob_mat <- matrix(0, nrow = NBlock_local, ncol = Day)
        for(i in 1:Day){
          s1 <- flow_list[[i]]
          idx <- match(BlockID, s1$CensusBlockGroup)
          mob_vec <- s1[idx, yindex + 1]
          mob_vec[is.na(mob_vec)] <- 0
          Mob_mat[, i] <- mob_vec
        }
        
        # Assign globals expected by the estimator functions (homo/power/exp)
        Mob   <<- Mob_mat
        NBlock <<- NBlock_local
        
        # For logging / output consistency:
        # - lag.max: keep induced-subgraph diameter (reflects structural disruption)
        # - slag: fixed to full lag definition used by W_full_list
        lag.max <<- lag.max_local
        
        invisible(TRUE)
      }
      
      #------------------------------------------------------------
      # (1) Point estimate on the full sample
      #------------------------------------------------------------
      build_globals_from_keep(BlockID_full)
      
      # Parameter bounds for different kernels
      lower_bound <<- 0
      upper_bound <<- 5
      lower_power <<- 40
      upper_power <<- 80
      
      # Fit all kernels day-by-day
      point_list <- lapply(1:Day, fit_one_day_allkernels)
      point_df <- do.call(rbind, point_list)
      point_dt <- as.data.table(point_df)
      point_dt[, `:=`(msa = Dname, region = region[s], yindex = yindex)]
      
      #------------------------------------------------------------
      # (2) Delete-p bootstrap for uncertainty quantification (structural robustness CI)
      #------------------------------------------------------------
      boot_res <- vector("list", B)
      n_fail <- 0L
      
      for(b in 1:B){
        keep_ids <- sample(
          BlockID_full,
          size = floor(length(BlockID_full) * keep_ratio),
          replace = FALSE
        )
        
        tryCatch({
          build_globals_from_keep(keep_ids)
          
          tmp_list <- lapply(1:Day, fit_one_day_allkernels)
          tmp <- as.data.table(do.call(rbind, tmp_list))
          
          tmp[, `:=`(
            boot = b,
            nblock = NBlock,
            lag = lag.max,         # induced-subgraph diameter (diagnostic)
            msa = Dname,
            region = region[s],
            yindex = yindex,
            keep_ratio = keep_ratio
          )]
          
          boot_res[[b]] <- tmp
        }, error = function(e){
          n_fail <<- n_fail + 1L
          boot_res[[b]] <- NULL
        })
      }
      
      boot_dt <- rbindlist(boot_res, fill = TRUE)
      
      #------------------------------------------------------------
      # Save full bootstrap trajectories (rho estimates per day and kernel)
      #------------------------------------------------------------
      boot_file <- paste0(datapath, "/", region[s], "/param/", Dname, "_SLDR_bootstrap_rho_", Yname[yindex], "_1.csv")
      boot_dt[, `:=`(B = B, n_fail = n_fail)]
      fwrite(boot_dt, boot_file)
      
      #------------------------------------------------------------
      # Summarize bootstrap confidence intervals by (day, kernel)
      #------------------------------------------------------------
      ci_dt <- boot_dt[, .(
        rho_ci_low    = quantile(rho, 0.025, na.rm = TRUE),
        rho_ci_high   = quantile(rho, 0.975, na.rm = TRUE),
        err_ci_low    = quantile(error, 0.025, na.rm = TRUE),
        err_ci_high   = quantile(error, 0.975, na.rm = TRUE),
        nblock_median = median(nblock, na.rm = TRUE),
        lag_median    = median(lag, na.rm = TRUE),
        n_boot_used   = sum(!is.na(rho))
      ), by = .(day, index)]
      
      #------------------------------------------------------------
      # Merge point estimates with bootstrap CIs and write output
      #------------------------------------------------------------
      out_dt <- merge(point_dt, ci_dt, by = c("day", "index"), all.x = TRUE)
      out_dt[, `:=`(
        msa = Dname,
        region = region[s],
        yindex = yindex,
        keep_ratio = keep_ratio,
        B = B,
        n_fail = n_fail
      )]
      
      out_file <- paste0(datapath, "/", region[s], "/param/", Dname, "_SLDR_params_CI_", Yname[yindex], "_1.csv")
      fwrite(out_dt, out_file)
      
      message("Done: ", Dname, " (bootstrap saved: ", boot_file, ")")
      print(d)
    }
  }
}





 




#----------SLDR fitting + CI----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in 1:7){
      
      Dname <<- Dname.set[d]
      source(paste0(codepath, "/sldr_global_vars_funs.R"))
      
      #------------------------------------------------------------
      # Read spatial tessellation (CBG polygons) for a given MSA
      #------------------------------------------------------------
      block_full <- sf::read_sf(paste0(geopath, "/", region[s], "/", Dname, ".geojson"))
      
      # IMPORTANT: Fix the ordering of spatial units to guarantee that
      # the adjacency matrix and mobility vectors are consistently aligned
      # across point estimation and bootstrap replicates.
      block_full <- block_full[order(block_full$CensusBlockGroup), ]
      BlockID_full <- unique(block_full$CensusBlockGroup)
      BlockID_full <- as.character(BlockID_full)
      
      #------------------------------------------------------------
      # Pre-load daily mobility data (avoid repeated disk I/O)
      #------------------------------------------------------------
      flow_list <- vector("list", Day)
      for(i in 1:Day){
        tmp <- read.csv(
          paste0(flowpath, "/", region[s], "/", Dname, "/Intra_Flow_", Datevalue[i], ".csv"),
          header = TRUE
        )
        # Ensure consistent ID type for fast matching
        tmp$CensusBlockGroup <- as.character(tmp$CensusBlockGroup)
        flow_list[[i]] <- tmp
      }
      
      #------------------------------------------------------------
      # Compute queen contiguity adjacency ONCE on the full sample
      #------------------------------------------------------------
      # Rationale: geometric contiguity (poly2nb) is computationally expensive.
      # For delete-p bootstrap over spatial units, we can form the induced subgraph
      # of the full contiguity graph by subsetting the full adjacency matrix,
      # which is equivalent to recomputing contiguity after node deletion because
      # removing units can only remove edges, not create new contiguity relations.
      Wd_full <- spdep::poly2nb(block_full, queen = TRUE)
      W1_full <- nb.mat(nb = Wd_full)
      rm(Wd_full)
      
      #------------------------------------------------------------
      # Helper: build global variables given a retained set of CBG IDs
      # (delete-p bootstrap via induced subgraph)
      #------------------------------------------------------------
      build_globals_from_keep <- function(keep_ids){
        
        keep_ids <- as.character(keep_ids)
        
        # Map retained IDs to indices of the full ordering and keep them sorted.
        # Sorting ensures that W1 and Mob_mat share the same row/column ordering.
        keep_idx <- match(keep_ids, BlockID_full)
        keep_idx <- keep_idx[!is.na(keep_idx)]
        keep_idx <- sort(keep_idx)
        
        # Induced subgraph adjacency matrix after deleting spatial units
        W1 <- W1_full[keep_idx, keep_idx, drop = FALSE]
        
        # Corresponding spatial unit IDs (aligned with W1)
        BlockID <- BlockID_full[keep_idx]
        NBlock_local <- length(BlockID)
        
        # Determine maximum spatial lag based on the induced subgraph diameter
        # (used to set the highest order of spatial dependence)
        g_sub <- igraph::graph_from_adjacency_matrix(W1, mode = "undirected", diag = FALSE)
        lag.max_local <- igraph::diameter(g_sub)
        
        # Generate higher-order queen adjacency matrices up to lag.max_local
        W <<- queen.adj(W1, lag.max_local)
        slag <<- seq_len(lag.max_local)
        
        #------------------------------------------------------------
        # Build the mobility matrix (rows aligned with the adjacency matrix)
        #------------------------------------------------------------
        # Using match() instead of join operations is faster and preserves ordering.
        Mob_mat <- matrix(0, nrow = NBlock_local, ncol = Day)
        for(i in 1:Day){
          s1 <- flow_list[[i]]
          idx <- match(BlockID, s1$CensusBlockGroup)
          mob_vec <- s1[idx, yindex + 1]
          mob_vec[is.na(mob_vec)] <- 0
          Mob_mat[, i] <- mob_vec
        }
        
        # Assign globals expected by the estimator functions (homo/power/exp)
        Mob   <<- Mob_mat
        NBlock <<- NBlock_local
        lag.max <<- lag.max_local
        
        invisible(TRUE)
      }
      
      #------------------------------------------------------------
      # (1) Point estimate on the full sample
      #------------------------------------------------------------
      build_globals_from_keep(BlockID_full)
      
      # Parameter bounds for different kernels
      lower_bound <<- 0
      upper_bound <<- 5
      lower_power <<- 10
      upper_power <<- 100
      
      # Fit all kernels day-by-day
      point_list <- lapply(1:Day, fit_one_day_allkernels)
      point_df <- do.call(rbind, point_list)
      point_dt <- as.data.table(point_df)
      point_dt[, `:=`(msa = Dname, region = region[s], yindex = yindex)]
      
      #------------------------------------------------------------
      # (2) Delete-p bootstrap for uncertainty quantification (structural robustness CI)
      #------------------------------------------------------------
      # We randomly retain keep_ratio of spatial units (CBGs), rebuild the induced
      # contiguity subgraph, and re-estimate rho for each day and each kernel.
      boot_res <- vector("list", B)
      n_fail <- 0L
      
      for(b in 1:B){
        keep_ids <- sample(
          BlockID_full,
          size = floor(length(BlockID_full) * keep_ratio),
          replace = FALSE
        )
        
        tryCatch({
          build_globals_from_keep(keep_ids)
          
          tmp_list <- lapply(1:Day, fit_one_day_allkernels)
          tmp <- as.data.table(do.call(rbind, tmp_list))
          
          tmp[, `:=`(
            boot = b,
            nblock = NBlock,
            lag = lag.max,
            msa = Dname,
            region = region[s],
            yindex = yindex,
            keep_ratio = keep_ratio
          )]
          
          boot_res[[b]] <- tmp
        }, error = function(e){
          n_fail <<- n_fail + 1L
          boot_res[[b]] <- NULL
        })
      }
      
      boot_dt <- rbindlist(boot_res, fill = TRUE)
      
      #------------------------------------------------------------
      # Save full bootstrap trajectories (rho estimates per day and kernel)
      #------------------------------------------------------------
      boot_file <- paste0(datapath, "/", region[s], "/param/",
                          Dname, "_SLDR_bootstrap_rho_", Yname[yindex], ".csv")
      boot_dt[, `:=`(B = B, n_fail = n_fail)]
      fwrite(boot_dt, boot_file)
      
      #------------------------------------------------------------
      # Summarize bootstrap confidence intervals by (day, kernel)
      #------------------------------------------------------------
      ci_dt <- boot_dt[, .(
        rho_ci_low    = quantile(rho, 0.025, na.rm = TRUE),
        rho_ci_high   = quantile(rho, 0.975, na.rm = TRUE),
        err_ci_low    = quantile(error, 0.025, na.rm = TRUE),
        err_ci_high   = quantile(error, 0.975, na.rm = TRUE),
        nblock_median = median(nblock, na.rm = TRUE),
        lag_median    = median(lag, na.rm = TRUE),
        n_boot_used   = sum(!is.na(rho))
      ), by = .(day, index)]
      
      #------------------------------------------------------------
      # Merge point estimates with bootstrap CIs and write output
      #------------------------------------------------------------
      out_dt <- merge(point_dt, ci_dt, by = c("day", "index"), all.x = TRUE)
      out_dt[, `:=`(
        msa = Dname,
        region = region[s],
        yindex = yindex,
        keep_ratio = keep_ratio,
        B = B,
        n_fail = n_fail
      )]
      
      out_file <- paste0(datapath, "/", region[s], "/param/",
                         Dname, "_SLDR_params_CI_", Yname[yindex], ".csv")
      fwrite(out_dt, out_file)
      
      message("Done: ", Dname, " (bootstrap saved: ", boot_file, ")")
      print(d)
    }
  }
}










#----------SLDR fitting + CI----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in 1:7){
      
      Dname <<- Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      
      #----------Read geometry once----------#
      block_full <- sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID_full <- unique(block_full$CensusBlockGroup)
      
      #----------Pre-read daily flows once (avoid repeated IO)----------#
      flow_list <- vector("list", Day)
      for(i in 1:Day){
        flow_list[[i]] <- read.csv(
          paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""),
          header = TRUE
        )
      }
      
      #----------Function: build globals given a kept ID set----------#
      build_globals_from_keep <- function(keep_ids){
        
        block.msa <- block_full[block_full$CensusBlockGroup %in% keep_ids, ]
        BlockID <- unique(block.msa$CensusBlockGroup)
        NBlock_local <- length(BlockID)  # NEW: avoid name clash
        
        # W: higher-order queen adjacency
        Wd <- spdep::poly2nb(block.msa, queen = TRUE)
        W1 <- nb.mat(nb = Wd)
        lag.max_local <- igraph::graph_from_adjacency_matrix(W1, mode = "directed") %>% igraph::diameter()
        
        W <<- queen.adj(W1, lag.max_local)
        slag <<- seq(1, lag.max_local, 1)
        
        # Mob: NBlock x Day
        Mob_mat <- matrix(0, nrow = NBlock_local, ncol = Day)
        for(i in 1:Day){
          s1 <- flow_list[[i]]
          s2 <- left_join(
            data.frame(CensusBlockGroup = BlockID),
            data.frame(CensusBlockGroup = s1$CensusBlockGroup, mob = s1[, yindex + 1]),
            by = "CensusBlockGroup"
          )
          s2$mob[is.na(s2$mob)] <- 0
          Mob_mat[, i] <- s2$mob
        }
        
        # assign globals expected by estimator functions
        Mob   <<- Mob_mat
        NBlock <<- NBlock_local
        lag.max <<- lag.max_local
        invisible(TRUE)
      }
      
      #----------(1) Point estimate on full sample----------#
      build_globals_from_keep(BlockID_full)
      
      lower_bound <<- 0
      upper_bound <<- 5
      lower_power <<- 10
      upper_power <<- 100
      
      point_list <- lapply(1:Day, fit_one_day_allkernels)
      point_df <- do.call(rbind, point_list)
      point_dt <- as.data.table(point_df)
      point_dt[, `:=`(msa = Dname, region = region[s], yindex = yindex)]  
      
      #----------(2) Bootstrap / delete-p% for CI----------#
      boot_res <- vector("list", B)
      n_fail <- 0L
      
      for(b in 1:B){
        keep_ids <- sample(BlockID_full, size = floor(length(BlockID_full) * keep_ratio), replace = FALSE)
        
        tryCatch({
          build_globals_from_keep(keep_ids)
          
          tmp_list <- lapply(1:Day, fit_one_day_allkernels)
          tmp <- as.data.table(do.call(rbind, tmp_list))
          
          tmp[, `:=`(
            boot = b,
            nblock = NBlock,
            lag = lag.max,
            msa = Dname,              
            region = region[s],       
            yindex = yindex,         
            keep_ratio = keep_ratio 
          )]
          
          boot_res[[b]] <- tmp
        }, error = function(e){
          n_fail <<- n_fail + 1L
          boot_res[[b]] <- NULL
        })
      }
      
      boot_dt <- rbindlist(boot_res, fill = TRUE)
      
      #----------Save full bootstrap trajectories (each rho)----------#
      boot_file <- paste(datapath, "/", region[s], "/param/", Dname, "_SLDR_bootstrap_rho_", Yname[yindex], ".csv", sep = "")
      boot_dt[, `:=`(B = B, n_fail = n_fail)] 
      fwrite(boot_dt, boot_file)
      
      #----------Summarize CI by (day, index)----------#
      ci_dt <- boot_dt[,.(
          rho_ci_low    = quantile(rho, 0.025, na.rm = TRUE),
          rho_ci_high   = quantile(rho, 0.975, na.rm = TRUE),
          err_ci_low    = quantile(error, 0.025, na.rm = TRUE),
          err_ci_high   = quantile(error, 0.975, na.rm = TRUE),
          nblock_median = median(nblock, na.rm = TRUE),
          lag_median    = median(lag, na.rm = TRUE),
          n_boot_used   = sum(!is.na(rho))     
        ), by = .(day, index)]
      
      #----------Merge point + CI----------#
      out_dt <- merge(point_dt, ci_dt, by = c("day","index"), all.x = TRUE)
      out_dt[, `:=`(
        msa = Dname,
        region = region[s],   # NEW
        yindex = yindex,      # NEW
        keep_ratio = keep_ratio,
        B = B,
        n_fail = n_fail
      )]
      
      #----------Write summary----------#
      out_file <- paste(
        datapath, "/", region[s], "/param/",
        Dname, "_SLDR_params_CI_", Yname[yindex], ".csv",
        sep = ""
      )
      fwrite(out_dt, out_file)
      
      message("Done: ", Dname, " (bootstrap saved: ", boot_file, ")")
      print(d)
    }
  }
}





#----------SLDR fitting----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in 1:Nmsa){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # date.gif<-date.sep(Datevalue,durdate)
      
      #----------BlockID, NBlock----------#
      block.msa<-sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID<-unique(block.msa$CensusBlockGroup)
      NBlock<-length(BlockID)
      #----------W: 1 to lag.max lag queen adjacency matrices----------#
      Wd <- spdep::poly2nb(block.msa, queen = T) 
      W1 <- nb.mat(nb=Wd)   # W1 <- spdep::nb2mat(Wd, style = "B") 
      lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed")%>%igraph::diameter()
      W <<- queen.adj(W1, lag.max)
      slag <<- seq(1,lag.max,1)
      
      #----------1.parameters----------#
      #----------Mobility: Mob----------#
      Mob<<-matrix(0, nrow=NBlock, ncol=Day)
      for(i in 1:Day){
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        Mob[,i] <- s2$mob
      } # day
      #----------Optimal parameters: rho----------#
      lower_bound<<-0
      upper_bound<<-5
      lower_power<<-10
      upper_power<<-100
      homo_list <- lapply(1:Day, function(g) {homo.day.est(g)})
      homo_result <- do.call(rbind, homo_list)%>%as.data.frame
      homo_result$index <- 'homo'
      power_list <- lapply(1:Day, function(g) {power.day.est(g)})
      power_result <- do.call(rbind, power_list)%>%as.data.frame
      power_result$index <- 'power'
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list)%>%as.data.frame
      exp_result$index <- 'exp'
      result <- rbind(homo_result, power_result, exp_result)
      #----------rho, error, convergence----------#
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/param/",Dname,"_SLDR_params_CI_",Yname[yindex],".csv",sep=""), row.names = FALSE)
      
      
      print(d)
    } # dis
  }# yindex
}# region



#----------SLDR fitting results for all msa_county----------#
for(s in 2:Nregion){
  for(yindex in 1:Ynum){
    dis.rho <- dis.r2 <- dis.cor <- NULL
    for(d in 1:Nmsa){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      date.gif<-date.sep(Datevalue,durdate)
      
      #----------params----------#
      sldr.params <- read.csv(paste(datapath,"/",region[s],"/params/",Dname,"_SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)
      sldr.params$event <- event
      sldr.params$msa <- Dname
      sldr.params$day <- as.Date(sldr.params$day)
      sldr.params <- left_join(sldr.params, date.gif, by='day')
      dis.rho <- plyr::rbind.fill(dis.rho, sldr.params)
      
      #----------r2----------# 
      sldr.fit <- read.csv(paste(datapath,"/",region[s],"/fit/",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
      r2 <- sldr.fit[, .(r2_homo = R2(empirical, homo),
                         r2_power = R2(empirical, power),
                         r2_exp = R2(empirical, exp)), by = day]
      r2$event <- event
      r2$msa <- Dname
      dis.r2 <- plyr::rbind.fill(dis.r2, r2)
      
      #----------rank cor: apply the rank correlations for each day, sorting by 'empirical'----------# 
      # sldr.fit <- sldr.fit[order(day, empirical)]
      rcor <- sldr.fit[, {
        # Calculate Pearson correlation
        pearson_homo <- cor(empirical, homo, method = "pearson", use = "complete.obs")
        pearson_power <- cor(empirical, power, method = "pearson", use = "complete.obs")
        pearson_exp <- cor(empirical, exp, method = "pearson", use = "complete.obs")
        # Calculate spearman correlation on the ranks
        spearman_homo <- cor(empirical, homo, method = "spearman", use = "complete.obs")
        spearman_power <- cor(empirical, power, method = "spearman", use = "complete.obs")
        spearman_exp <- cor(empirical, exp, method = "spearman", use = "complete.obs")
        # Calculate Kendall correlation on the ranks
        kendall_homo <- cor(empirical, homo, method = "kendall", use = "complete.obs")
        kendall_power <- cor(empirical, power, method = "kendall", use = "complete.obs")
        kendall_exp <- cor(empirical, exp, method = "kendall", use = "complete.obs")
        # Return the correlations for this day
        .(pearson_homo = pearson_homo, pearson_power = pearson_power, pearson_exp = pearson_exp,
          spearman_homo = spearman_homo, spearman_power = spearman_power, spearman_exp = spearman_exp,
          kendall_homo = kendall_homo, kendall_power = kendall_power, kendall_exp = kendall_exp)
      }, by = day]
      rcor$event <- event
      rcor$msa <- Dname
      dis.cor <- plyr::rbind.fill(dis.cor, rcor)
      print(d)
    } # dis
    write.csv(dis.rho, file=paste(datapath,"/",region[s],"/SLDR_params_",Yname[yindex],".csv",sep=""), row.names = FALSE)
    write.csv(dis.r2, file=paste(datapath,"/",region[s],"/SLDR_r2_",Yname[yindex],".csv",sep=""), row.names = FALSE)
    write.csv(dis.cor, file=paste(datapath,"/",region[s],"/SLDR_rank_cor_",Yname[yindex],".csv",sep=""), row.names = FALSE)
    
  } # yindex
} # region



