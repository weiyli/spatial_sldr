# rm(list = ls())


# Scale-invariance and structure-sensitivity tests for rho


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

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
dset<-which(Dname.set%in%c("Atlanta","San Francisco"))


#----------Load function----------#
#----------Preserve total outflow----------#
# Rescale x so that sum(x) matches a target total
rescale_to_total <- function(x, target_total) {
  s_new <- sum(x)
  if (is.finite(s_new) && s_new > 0) return(x * (target_total / s_new))
  x
}

#----------Experiment I: uniform scaling (NO total preservation)----------#
# NOTE (EN): Intentionally changes total volume (for scale invariance test).
perturb_scale <- function(x_base, alpha, target_total) {
  x_new <- x_base * alpha
  rescale_to_total(x_new, target_total)  
}

#----------Experiment II: rewiring (stratified permutation)----------#
# NOTE (EN): Permute within core/periphery groups separately; preserves totals and group marginals.
perturb_rewire_strat <- function(x_base, core_ids, rp = 0.3, target_total) {
  ids_all <- seq_along(x_base)
  peri_ids <- setdiff(ids_all, core_ids)
  
  x_new <- x_base
  
  k1 <- max(1, floor(length(core_ids) * rp))
  k2 <- max(1, floor(length(peri_ids) * rp))
  
  idx1 <- sample(core_ids, k1)
  idx2 <- sample(peri_ids, k2)
  
  x_new[idx1] <- sample(x_base[idx1], k1, replace = FALSE)
  x_new[idx2] <- sample(x_base[idx2], k2, replace = FALSE)
  
  rescale_to_total(x_new, target_total)
}

#----------Experiment III: permutation null (total-preserving)----------#
# NOTE (EN): Shuffles values across CBGs; preserves distribution and total, destroys spatial alignment.
perturb_permute <- function(x_base, eta = 0.3, target_total) {
  n <- length(x_base)
  k <- floor(n * eta)
  idx <- sample(seq_len(n), k)
  
  x_new <- x_base
  x_new[idx] <- sample(x_base[idx], k, replace = FALSE)
  
  rescale_to_total(x_new, target_total)
}


#----------Run scenarios with multiple perturbation magnitudes----------#
disturb_set <- c(0.01, seq(0.10,0.90,0.10))

disturb_set <- c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)

# NOTE (EN): scenario_set includes both parameterized and non-parameterized scenarios.
scenario_set <- c(
  "I_scale",
  "II_rewire_strat",
  "III_perm_null"
)

n_rep <- 100              # experiment repeat times (replicates)
base_seed <- 12345        # global reproducibility seed

#----------SLDR fitting----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in dset){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # date.gif<-date.sep(Datevalue,durdate)
      
      #----------BlockID, NBlock----------#
      block.msa<-sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID<-unique(block.msa$CensusBlockGroup)
      NBlock<-length(BlockID)
      
      #----------W: 1 to lag.max lag queen adjacency matrices----------#
      # NOTE (EN): This block is expensive; compute ONCE per (s,d).
      Wd <- spdep::poly2nb(block.msa, queen = TRUE)
      W1 <- nb.mat(nb=Wd)
      lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed") %>% igraph::diameter()
      W <<- queen.adj(W1, lag.max)
      slag <<- seq(1, lag.max, 1)
      
      #----------Select one base day----------#
      # NOTE (EN): sample one day as baseline (same baseline used across all scenarios and disturb).
      set.seed(base_seed + 1000*s + 10*d + yindex)
      base_idx <- sample(1:(Day/2), 1)
      base_date <- Datevalue[base_idx]
      
      s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_", base_date, ".csv", sep = ""), header=TRUE)
      s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                      data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                      by="CensusBlockGroup")
      s2$mob[is.na(s2$mob)] <- 0
      x_base <- s2$mob
      
      #----------Define core IDs for structure experiments----------#
      core_frac <- 0.1
      k <- max(1, floor(length(x_base) * core_frac))
      core_ids <- order(x_base, decreasing = TRUE)[1:k]
      
      #----------Optimal parameters: rho----------#
      lower_bound<<-0
      upper_bound<<-5
      
      #----------Run all scenarios and perturbation magnitudes----------#
      out_all <- list()
      idx_out <- 1
      
      for(scenario in scenario_set){
        for(disturb in disturb_set){
          
          target_total <- sum(x_base) * (1-disturb)
          
          #----------Mobility: Mob (replicates for this scenario & disturb)----------#
          Mob<<-matrix(0, nrow=NBlock, ncol=n_rep)
          
          # NOTE (EN): Use a scenario- and disturb-specific seed, so results are reproducible
          # but different across scenarios/magnitudes.
          set.seed(base_seed + 1e5*s + 1e3*d + 100*yindex + 10*match(scenario, scenario_set) + round(disturb*1000))
          
          for(i in 1:n_rep){
            #----------experiment----------#
            
            if(scenario == "I_scale"){
              # disturb is the scaling factor alpha (0.01, 0.1, 0.5, 0.9)
              alpha_use <- max(0, 1 - disturb)
              Mob[,i] <- perturb_scale(x_base, alpha = alpha_use, target_total = target_total)
              
            } else if(scenario == "II_rewire_strat"){
              # NOTE (EN): rewiring has no continuous magnitude; repeat with different permutations.
              Mob[,i] <- perturb_rewire_strat(x_base, core_ids, rp = disturb, target_total = target_total)
              
            } else if(scenario == "III_perm_null"){
              # NOTE (EN): permutation has no continuous magnitude; we still repeat under each disturb for consistency.
              Mob[,i] <- perturb_permute(x_base, eta = disturb, target_total = target_total)
              
            } else {
              stop("Unknown scenario: ", scenario)
            }
          } # rep
          
          #----------Estimate rho for each replicate----------#
          exp_list <- lapply(1:n_rep, function(g) {exp.day.est(g)})
          exp_result <- do.call(rbind, exp_list) %>% as.data.frame
          
          #----------Assemble output----------#
          tmp <- data.frame(
            base_day = base_date,
            lag = lag.max,
            rep = 1:n_rep,
            scenario = scenario,
            disturb = disturb,
            rho = exp_result[,1],
            error = exp_result[,2]
          )
          
          
          out_all[[idx_out]] <- tmp
          idx_out <- idx_out + 1
        } # disturb
      } # scenario
      
      SLDR_all <- do.call(rbind, out_all)
      
      #----------Write output----------#
      write.csv(SLDR_all,file=paste(datapath,"/",region[s],"/params/",Dname,"_disturb_",Yname[yindex],".csv",sep=""),
                row.names = FALSE)
      
      print(d)
    } # dis
  } # yindex
} # region

