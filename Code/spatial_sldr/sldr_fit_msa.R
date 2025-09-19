# rm(list = ls())

# SLDR fitting and higer-order queen weight matrix at the msa level for rho estimation


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


#----------county in each msa----------#
county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'

#----------W: Calculate the max lag h with sum(W[[h]])>0 for msa or county----------#
for(s in 2:Nregion){
  dis.W<-NULL
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    
    #----------BlockID, NBlock----------#
    block.msa <- sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
    block.msa <- block.msa %>% mutate(pop = ifelse(is.na(pop), 0, pop))
    BlockID<-unique(block.msa$CensusBlockGroup)
    NBlock<-length(BlockID)
    pop <- block.msa$pop
    area <- block.msa$area_m2
    
    #----------W: 1 to lag.max lag queen adjacency matrices----------#
    Wd <- spdep::poly2nb(block.msa, queen = T) 
    W1 <- nb.mat(nb=Wd)   # W1 <- spdep::nb2mat(Wd, style = "B") 
    lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed")%>%igraph::diameter()
    W <- queen.adj(W1, lag.max)
    
    #----------dis.W: queen matrices for all msas and counties----------#
    region.W <- data.frame(
      msa = Dname,
      nblock = NBlock,
      lag.max = lag.max,
      lag = 1:lag.max,
      nozero = sapply(1:lag.max, function(h) sum(W[[h]])),
      pop.adj = sapply(1:lag.max, function(h) sum(pop[which(W[[h]]==1, arr.ind = TRUE)[,2]])),
      area.adj = sapply(1:lag.max, function(h) sum(area[which(W[[h]]==1, arr.ind = TRUE)[,2]])),
      pop = sum(pop),
      area = sum(area)
    )
    dis.W <- plyr::rbind.fill(dis.W, region.W)
    print(d)
    
  } # dis
  write.csv(dis.W,file=paste(datapath,"/",region[s],"/queen_lag.csv",sep=""),row.names = FALSE)
} # region


#----------SLDR fitting----------#
for(s in 2:Nregion){
  for(yindex in 1:Ynum){
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
      lower_bound<<--5
      upper_bound<<-5
      homo_list <- lapply(1:Day, function(g) {homo.day.est(g)})
      homo_result <- do.call(rbind, homo_list)%>%as.data.frame
      homo_result$index <- 'homo'
      # linear_list <- lapply(1:Day, function(g) {linear.day.est(g)})
      # linear_result <- do.call(rbind, linear_list)%>%as.data.frame
      # linear_result$index <- 'linear'
      power_list <- lapply(1:Day, function(g) {power.day.est(g)})
      power_result <- do.call(rbind, power_list)%>%as.data.frame
      power_result$index <- 'power'
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list)%>%as.data.frame
      exp_result$index <- 'exp'
      result <- rbind(homo_result, power_result, exp_result)
      #----------rho, error, convergence----------#
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/params/",Dname,"_SLDR_params_",Yname[yindex],".csv",sep=""), row.names = FALSE)
      
      
      #----------2.fitting----------#
      SLDR.fit<-NULL
      SLDR.homo <- subset(SLDR,index=='homo')
      SLDR.power <- subset(SLDR,index=='power')
      SLDR.exp <- subset(SLDR,index=='exp')
      for(i in 1:Day){
        # empirical
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        M.day <- s2$mob
        # homo
        rho<-SLDR.homo$rho[i]
        lag.h <- matrix(rho, nrow = NBlock, ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.homo <- rowSums(lag.h*WM.h)
        # power
        rho<-SLDR.power$rho[i]
        lag.h <- matrix(rep(slag^(-rho), NBlock), ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.power <- rowSums(lag.h*WM.h)
        # exp
        rho<-SLDR.exp$rho[i]
        lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.exp <- rowSums(lag.h*WM.h)
        SLDR.fit <- plyr::rbind.fill(SLDR.fit, data.frame(day=Datevalue[i], 
                                                          empirical=M.day, 
                                                          homo=M.day.sim.homo,
                                                          power=M.day.sim.power,
                                                          exp=M.day.sim.exp))
      } # day
      write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fit/",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""),row.names = FALSE)
      
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




#----------Part2: disaster----------#
flowpath <- 'D:/ood/Data/Flow/disaster'
datapath <- 'D:/ood/Data/spatial_sldr/disaster'

#----------W: Calculate the max lag h with sum(W[[h]])>0 for msa or county----------#
for(s in 2:Nregion){
  dis.W<-NULL
  for(d in Nmsa+(1:Ndis)){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    
    #----------BlockID, NBlock----------#
    block.msa<-sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
    block.msa <- block.msa %>% mutate(pop = ifelse(is.na(pop), 0, pop))
    # nrow(block.msa)
    # block.msa <- block.msa[st_geometry_type(block.msa) %in% c("POLYGON", "MULTIPOLYGON"), ]
    # nrow(block.msa)
    BlockID<-unique(block.msa$CensusBlockGroup)
    NBlock<-length(BlockID)
    pop <- block.msa$pop
    area <- block.msa$area_m2
    
    #----------W: 1 to lag.max lag queen adjacency matrices----------#
    Wd <- spdep::poly2nb(block.msa, queen = T) 
    W1 <- nb.mat(nb=Wd)   # W1 <- spdep::nb2mat(Wd, style = "B") 
    lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed")%>%igraph::diameter()
    W <- queen.adj(W1, lag.max)
    
    #----------dis.W: queen matrices for all msas and counties----------#
    region.W <- data.frame(
      msa = Dname,
      nblock = NBlock,
      lag.max = lag.max,
      lag = 1:lag.max,
      nozero = sapply(1:lag.max, function(h) sum(W[[h]])),
      pop.adj = sapply(1:lag.max, function(h) sum(pop[which(W[[h]]==1, arr.ind = TRUE)[,2]])),
      area.adj = sapply(1:lag.max, function(h) sum(area[which(W[[h]]==1, arr.ind = TRUE)[,2]])),
      pop = sum(pop),
      area = sum(area)
    )
    dis.W <- plyr::rbind.fill(dis.W, region.W)
    print(d)
    
  } # dis
  write.csv(dis.W,file=paste(datapath,"/",region[s],"/queen_lag.csv",sep=""),row.names = FALSE)
} # region



#----------SLDR fitting----------#
for(s in 2:Nregion){
  for(yindex in 1:Ynum){
    for(d in Nmsa+(1:Ndis)){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # date.gif<-date.sep(Datevalue,durdate)
      
      #----------BlockID, NBlock----------#
      block.msa<-sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      # block.msa <- block.msa[st_geometry_type(block.msa) %in% c("POLYGON", "MULTIPOLYGON"), ]
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
      lower_bound<<--5
      upper_bound<<-5
      homo_list <- lapply(1:Day, function(g) {homo.day.est(g)})
      homo_result <- do.call(rbind, homo_list)%>%as.data.frame
      homo_result$index <- 'homo'
      # linear_list <- lapply(1:Day, function(g) {linear.day.est(g)})
      # linear_result <- do.call(rbind, linear_list)%>%as.data.frame
      # linear_result$index <- 'linear'
      power_list <- lapply(1:Day, function(g) {power.day.est(g)})
      power_result <- do.call(rbind, power_list)%>%as.data.frame
      power_result$index <- 'power'
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list)%>%as.data.frame
      exp_result$index <- 'exp'
      result <- rbind(homo_result, power_result, exp_result)
      #----------rho, error, convergence----------#
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/params/",event,"_",Dname,"_SLDR_params_",Yname[yindex],".csv",sep=""), row.names = FALSE)
      
      #----------2.fitting----------#
      SLDR.fit<-NULL
      SLDR.homo <- subset(SLDR,index=='homo')
      SLDR.power <- subset(SLDR,index=='power')
      SLDR.exp <- subset(SLDR,index=='exp')
      for(i in 1:Day){
        # empirical
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        M.day <- s2$mob
        # homo
        rho<-SLDR.homo$rho[i]
        lag.h <- matrix(rho, nrow = NBlock, ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.homo <- rowSums(lag.h*WM.h)
        # power
        rho<-SLDR.power$rho[i]
        lag.h <- matrix(rep(slag^(-rho), NBlock), ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.power <- rowSums(lag.h*WM.h)
        # exp
        rho<-SLDR.exp$rho[i]
        lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.exp <- rowSums(lag.h*WM.h)
        SLDR.fit <- plyr::rbind.fill(SLDR.fit, data.frame(day=Datevalue[i], 
                                                          empirical=M.day, 
                                                          homo=M.day.sim.homo,
                                                          power=M.day.sim.power,
                                                          exp=M.day.sim.exp))
      } # day
      write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fit/",event,"_",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""),row.names = FALSE)
      
      print(d)
    } # dis
  }# yindex
}# region


#----------SLDR fitting results for all msa_county----------#
for(s in 2:Nregion){
  for(yindex in 1:Ynum){
    dis.rho <- dis.r2 <- dis.cor <- NULL
    for(d in Nmsa+(1:Ndis)){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      date.gif<-date.sep.dis(Datevalue,startdate,enddate)
      
      #----------params----------#
      sldr.params <- read.csv(paste(datapath,"/",region[s],"/params/",event,"_",Dname,"_SLDR_params_",Yname[yindex],".csv",sep=""), header=TRUE)
      sldr.params$event <- event
      sldr.params$msa <- Dname
      sldr.params$day <- as.Date(sldr.params$day)
      sldr.params <- left_join(sldr.params, date.gif, by='day')
      dis.rho <- plyr::rbind.fill(dis.rho, sldr.params)
      
      #----------r2----------# 
      sldr.fit <- read.csv(paste(datapath,"/",region[s],"/fit/",event,"_", Dname, if (region[s]=='county') paste("_", county.name, sep = ""),"_SLDR_fit_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
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




