# rm(list = ls())

# simulate the spatially correlated SIR (Susceptible–Infected–Recovered) model on a graph while varying the spatial range exponent rho to analyze its impact on disease spread dynamics.


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Load packages----------#
library(igraph)
library(doParallel)
library(foreach)


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


#-----------------row_normalize: Row normalization--------------------#
#' @description Normalize each row of a matrix so that the row sum equals 1.
#' @param weight_matrix A square matrix (e.g., OD matrix or spatial weight matrix).
#' @return Row-normalized weight matrix.
row_normalize <- function(weight_matrix) {
  row_sums <- rowSums(weight_matrix)  # Compute row sums
  row_sums[row_sums == 0] <- 1        # Avoid division by zero for empty rows
  normalized_matrix <- weight_matrix / row_sums  # Element-wise division
  return(normalized_matrix)
}


#----------------- update_SIR: Efficiently update SIR states using matrices --------------------#
#' @description Update SIR model states based on adjacency and spatial weights using matrix operations.
#' @param node_pop      Population sizes for each node (block group).
#' @param node_S        Number of susceptible individuals per node.
#' @param node_I        Number of infected individuals per node.
#' @param node_R        Number of recovered individuals per node.
#' @param weight_matrix Spatial decay weight matrix.
#' @param beta          Infection rate (probability of transmission per contact).
#' @param mu            Recovery rate (probability of an infected individual recovering per step).
#' @return Updated node states (S, I, R counts for each node).
update_SIR <- function(node_pop, node_S, node_I, node_R, weight_matrix, beta, mu) {
  
  N <- length(node_pop)
  
  # Compute infected ratios safely, avoiding division by zero
  infected_ratio <- ifelse(node_pop > 0, node_I / node_pop, 0)
  
  # Compute spatial exposure risk (original weight-based transmission)
  exposure_risk <- beta * (weight_matrix %*% infected_ratio) 
  
  # Compute effective infection probability based on weighted exposure
  # infection_prob <- 1 - exp(-exposure_risk)  # Convert to probability of infection
  infection_prob <- exposure_risk
  
  # Simulate new infections (S -> I) using binomial draws
  new_infections <- rbinom(N, node_S, infection_prob)
  
  # Simulate recoveries (I -> R) using binomial draws
  new_recoveries <- rbinom(N, node_I, mu)
  
  # Update S, I, R values
  node_S <- node_S - new_infections
  node_I <- node_I + new_infections - new_recoveries
  node_R <- node_R + new_recoveries
  
  return(list(node_S = node_S, node_I = node_I, node_R = node_R))
}


#-----------------run_SIR_simulation: SIR simulation with spatial correlation--------------------#
#' @description Run an SIR (Susceptible-Infected-Recovered) model simulation considering spatial correlation.
#' @param spatial_pop_data  Spatial polygons (e.g., census blocks, counties, MSAs) and population data.
#' @param beta          Infection rate (probability of transmission per contact).
#' @param mu            Recovery rate (probability of an infected individual recovering per step).
#' @param rho           Spatial lag decay rate (controls the effect of spatial correlation).
#' @param N0            Initial number of infected nodes (ensures it does not exceed total nodes).
#' @param I0            Initial proportion of the population in each initially infected nodes.
#' @param time_steps    Number of simulation steps (iterations).
#' @param simulations   Number of times the SIR model simulation is repeated
#' @return A data.table containing the number of susceptible, infected, and recovered individuals at each time step.
#' @examples result_data <- run_SIR_simulation(spatial_pop_data, rho = 0.2, beta = 0.3, mu = 0.1, time_steps = 100, N0 = 5, I0 = 0.01)
run_SIR_simulation <- function(spatial_pop_data, rho = 0.2, beta = 0.3, mu = 0.1, N0 = 5, I0 = 0.01, 
                               time_steps = 1000, simulations = 100) {
  # Compute adjacency & distance matrices
  Wd <- spdep::poly2nb(spatial_pop_data, queen = TRUE)  
  adj_matrix <- nb.mat(nb = Wd)  
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed")
  dist_matrix <- igraph::distances(g)  
  
  # Compute spatial lag decay weights
  weight_matrix <- exp(-dist_matrix / rho)
  weight_matrix[is.infinite(weight_matrix)] <- 0  
  # weight_matrix <- weight_matrix / rowSums(weight_matrix)  # Row normalization
  weight_matrix <- row_normalize(weight_matrix)
  
  # Function to run a single SIR simulation
  single_run <- function() {
    
    node_pop <- spatial_pop_data$pop  
    N <- length(node_pop)  
    
    # Start with all individuals susceptible
    node_S <- node_pop  
    node_I <- rep(0, N)  
    node_R <- rep(0, N)  
    
    # Ensure initial infections do not exceed total population
    N0 <- min(N0, N)
    initial_infected <- sample(which(node_pop != 0), N0)
    node_I[initial_infected] <- round(node_pop[initial_infected] * I0)  
    node_S[initial_infected] <- node_pop[initial_infected] - node_I[initial_infected]
    
    # Simulation loop
    results <- data.table(time = integer(), susceptible = integer(), infected = integer(), recovered = integer())
    no_infection_count <- 0  
    
    for (t in 1:time_steps) {
      updated_states <- update_SIR(node_pop, node_S, node_I, node_R, weight_matrix, beta, mu)
      
      node_S <- updated_states$node_S
      node_I <- updated_states$node_I
      node_R <- updated_states$node_R
      
      total_infected <- sum(node_I)  
      
      # Store results
      results <- rbind(results, data.table(
        time = t,
        susceptible = sum(node_S),
        infected = sum(node_I),
        recovered = sum(node_R)
      ))
      
      # Early stopping if no new infections for 5 consecutive steps
      if (total_infected == 0) {
        no_infection_count <- no_infection_count + 1
      } else {
        no_infection_count <- 0  
      }
      
      if (no_infection_count >= 5) {
        message("Simulation stopped early at t=", t, " due to no new infections for 5 consecutive time steps.")
        break  
      }
    }  
    return(results)
  }
  
  # Run multiple simulations and average the results
  all_runs <- lapply(1:simulations, function(x) single_run())
  
  # Combine all simulations and compute means
  avg_results <- rbindlist(all_runs)[, .(
    susceptible = mean(susceptible),
    infected = mean(infected),
    recovered = mean(recovered)
  ), by = time]
  
  avg_results[, pop := sum(spatial_pop_data$pop)]
  
  return(avg_results)
}



#----------rho and msa sorted by median_rho: data from sldr_fit.R----------#
dis.rho <- read.csv(paste(datapath,"/msa/SLDR_params_",Yname[2],".csv",sep=""), header=TRUE)%>%setDT
rho.msa <- dis.rho[, .(
  max_rho = max(rho),
  max_day = day[which(rho==max(rho))]
), by = .(msa, index)]
rho.msa <- rho.msa[index == "exp"]
rho.rank <- rho.msa[order(max_rho)]
rho.rank[, msa_rank := .I]
# select.msa <- rho.rank[msa_rank %in% c(1, 12)]
select.msa <- rho.rank
dset <- match(select.msa$msa, Dname.set)
#----------initial parameters----------#
# rho.set <- seq(0.2,0.6,0.05)    # spatial decay rate
# beta.set <- c(0.3,0.5,0.8)      # infection rate
# mu.set <- c(0.1,0.15,0.2)       # recovery rate
beta.set <- c(0.3)      # infection rate
mu.set <- c(0.1)       # recovery rate
time_steps <- 1000     # number of simulation steps
simulations <- 10      # number of repeated simulations
N0 <- 5                # number of initially infected nodes
I0 <- 0.01             # initial infected probability
sir.all <- NULL
for(i in 1:length(dset)){
  
  #----------spatial_pop_data----------#
  d <- dset[i]
  Dname <<- Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  block.msa <- block.msa %>% mutate(pop = ifelse(is.na(pop), 0, pop))
  # block.msa <- block.msa[!is.na(block.msa$pop),]
  
  #----------Generate all parameter combinations----------#
  day.rho <- select.msa[i,]
  max.rho <- day.rho$max_rho
  max.day <- day.rho$max_day
  # det.rho <- seq(0, 0.9, 0.1)
  det.rho <- seq(0,0.5,0.25)
  rho.set <- c(max.rho * (1 - det.rho), max.rho * (1 + det.rho))
  day.set <- rep(max.day,length(rho.set))
  # c(rep(max.day,length(det.rho)),rep(max.day,length(det.rho)))
  det.rho.set <- c(-det.rho,det.rho)
  rho.day.msa <- data.frame(rho=rho.set, day=day.set, det.rho=det.rho.set)%>%unique
  # Generate all parameter combinations
  param.set <- expand.grid(rho = rho.day.msa$rho, beta = beta.set, mu = mu.set, stringsAsFactors = FALSE)
  
  # #----------Data from sldr_fit.R: SLDR_fit----------# 
  # sldr.fit <- read.csv(paste(datapath,"/msa/fit/",Dname,"_SLDR_fit_",Yname[2],".csv",sep=""), header=TRUE)
  # sldr.fit.day <- subset(sldr.fit,day==max.day)
  # sldr.fit$res<-sldr.fit$empirical-sldr.fit$exp
  
  #----------one simulation over all parameter combinations----------#
  sir.msa.list <- lapply(seq_len(nrow(param.set)), function(j) {
    params <- param.set[j, ]  
    sir.msa.rho <- run_SIR_simulation(spatial_pop_data = block.msa, 
                                      rho = params$rho, 
                                      beta = params$beta, 
                                      mu = params$mu, 
                                      N0, I0,
                                      time_steps, 
                                      simulations)
    sir.msa.rho$beta <- params$beta
    sir.msa.rho$mu <- params$mu
    sir.msa.rho$rho <- params$rho
    return(sir.msa.rho)
  })
  sir.msa <- do.call(rbind, sir.msa.list)
  sir.msa <- left_join(sir.msa,rho.day.msa,by=c("rho"))
  sir.msa$msa <- Dname
  sir.all <- plyr::rbind.fill(sir.all, sir.msa)
  print(Dname)
} # msa
# write.csv(sir.all, file=paste(datapath,"/msa/SIR_model_",Yname[2],".csv",sep=""), row.names = FALSE)
# write.csv(sir.all, file=paste(datapath,"/msa/SIR_model_all_",Yname[2],".csv",sep=""), row.names = FALSE) 0-0.9
write.csv(sir.all, file=paste(datapath,"/msa/SIR_model_rho_",Yname[2],".csv",sep=""), row.names = FALSE)





