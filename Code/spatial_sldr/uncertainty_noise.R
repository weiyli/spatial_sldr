# rm(list = ls())

# Scale-invariance and structure-sensitivity tests for rho


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
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

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
dset <- c(8)


#----------Experiment II(a): multiplicative noise (total-preserving)----------#
rescale_to_total <- function(x, target_total) {
  s_new <- sum(x)
  if (is.finite(s_new) && s_new > 0) return(x * (target_total / s_new))
  x
}

perturb_mult_noise <- function(x_base, eta = 0.03, target_total) {
  eps <- runif(length(x_base), -eta, eta)
  x_new <- x_base * (1 + eps)
  x_new[x_new < 0] <- 0
  rescale_to_total(x_new, target_total)
}

#----------Perturbation magnitudes----------#
# eta_set <- c(0.05, 0.10)
eta_set <- c(0.05)

#----------MC settings----------#
n_rep     <- 100
base_seed <- 12345

#----------Main loop----------#
for (s in 2:Nregion) {
  for (yindex in 2:2) {
    for (d in dset) {
      
      Dname <<- Dname.set[d]
      source(paste0(codepath, "/sldr_global_vars_funs.R"))
      
      #----------BlockID, NBlock----------#
      block.msa <- sf::read_sf(paste0(geopath, "/", region[s], "/", Dname, ".geojson"))
      BlockID   <- unique(block.msa$CensusBlockGroup)
      NBlock    <- length(BlockID)
      
      #----------W (queen adjacency up to lag.max)----------#
      Wd <- spdep::poly2nb(block.msa, queen = TRUE)
      W1 <- nb.mat(nb = Wd)
      lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed") %>% igraph::diameter()
      W    <<- queen.adj(W1, lag.max)
      slag <<- seq(1, lag.max, 1)

      #----------Mobility: Mob----------#
      Mob_day <- matrix(0, nrow=NBlock, ncol=Day)
      for(i in 1:Day){
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        Mob_day[,i] <- s2$mob
      } # day
      
      # day-specific totals (total-preserving)
      target_total_vec <- colSums(Mob_day)
      
      #----------Rho bounds----------#
      lower_bound <<- 0
      upper_bound <<- 5
      
      #----------Build task grid = Day * n_rep----------#
      n_task <- Day * n_rep
      task_id <- seq_len(n_task)
      
      # labels for each task (day varies slow; rep varies fast)
      day_idx_vec <- rep(1:Day, each = n_rep)       
      rep_idx_vec <- rep(1:n_rep, times = Day)       
      
      #----------Run eta grid----------#
      out_all <- vector("list", length(eta_set))
      for (j in seq_along(eta_set)) {
        
        eta <- eta_set[j]
        
        # 1) generate perturbed x vectors for ALL (day,rep) tasks
        x_list <- lapply(task_id, function(k) {
          
          g   <- day_idx_vec[k]     # day index
          rep <- rep_idx_vec[k]     # replicate id
          
          x_base <- Mob_day[, g]
          target_total <- target_total_vec[g]
          
          # reproducible seed per (s,d,yindex,eta,day,rep)
          set.seed(
            base_seed +
              1e5 * s +
              1e3 * d +
              100 * yindex +
              round(eta * 1000) +
              10 * g +
              rep
          )
          
          perturb_mult_noise(x_base, eta = eta, target_total = target_total)
        })

        Mob <<- do.call(cbind, x_list)  # NBlock x (Day*n_rep)
        
        # 2) estimate rho for each task column
        exp_list   <- lapply(task_id, exp.day.est)
        exp_result <- do.call(rbind, exp_list) %>% as.data.frame()
        
        rho_vec <- exp_result[, 1]
        err_vec <- exp_result[, 2]
        
        # 3) build output with your preferred day labeling
        dt_task <- data.table::data.table(
          msa       = Dname,
          day       = Datevalue[day_idx_vec],   
          rep       = rep_idx_vec,
          lag       = lag.max,
          eta       = eta,
          total_flow = target_total_vec[day_idx_vec],
          rho       = rho_vec,
          error     = err_vec
        )
        
        out_all[[j]] <- dt_task
      } # eta
      
      SLDR_all <- do.call(rbind, out_all)
     
      out_file <- paste0(datapath, "/", region[s], "/param/", Dname, "_uncertainty_", Yname[yindex], ".csv")
      write.csv(SLDR_all, file = out_file, row.names = FALSE)
      
      message("Done: ", Dname, " | region=", region[s], " | tasks=", n_task, " | file=", out_file)
      
    }
  }
}




