# rm(list = ls())

# SLDR fitting and higer-order queen weight matrix for msa or county level for rho


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo/cbg'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


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
county<-read.csv(paste(geopath,"/urban/select_county.csv",sep=""), header=TRUE)
county$state.fips <- sprintf("%02d", as.numeric(county$state.fips))
county$county.fips <- sprintf("%03d", as.numeric(county$county.fips))


#----------W: Calculate the max lag h with sum(W[[h]])>0 for msa or county----------#
dis.W<-NULL
for(d in 1:Nmsa){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  county.msa <- subset(county,msa.id==d)
  Ncounty <- nrow(county.msa)
  
  for(f in 1:Ncounty){
    
    county.id <- paste(county.msa$state.fips[f],county.msa$county.fips[f],sep='')
    block.msa<-sf::read_sf(paste(geopath,"/urban/",Dname,"_",county.id,".geojson",sep=""))
    BlockID<-unique(block.msa$CensusBlockGroup)
    NBlock<-length(BlockID)
    pop <- block.msa$pop
    area <- block.msa$area
    
    #----------W: 1 to lag.max lag queen adjacency matrices----------#
    Wd <- spdep::poly2nb(block.msa, queen = T) 
    W1 <- nb.mat(nb=Wd)   # W1 <- spdep::nb2mat(Wd, style = "B") 
    lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed")%>%igraph::diameter()
    W <- queen.adj(W1, lag.max)
    
    #----------dis.W: queen matrices for all msas and counties----------#
    region.W <- data.frame(
      msa = Dname,
      county = county.id,
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
  } # county
  print(d)
} # dis
write.csv(dis.W,file=paste(datapath,"/urban/queen_lag.csv",sep=""),row.names = FALSE)


#----------SLDR fitting----------#
for(yindex in 1:Ynum){
  
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    county.msa <- subset(county,msa.id==d)
    Ncounty <- nrow(county.msa)
    
    for(f in 1:Ncounty){
      
      #----------BlockID, NBlock----------#
      county.id <- paste(county.msa$state.fips[f],county.msa$county.fips[f],sep='')
      block.msa<-sf::read_sf(paste(geopath,"/urban/",Dname,"_",county.id,".geojson",sep=""))
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
        s1 <- read.csv(paste(flowpath,"/msa/", Dname, "/Intra_Flow_", Datevalue[i], ".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        Mob[,i] <- s2$mob
      } # day
      #----------Optimal parameters: rho----------#
      lower_bound<<-0
      upper_bound<<-1
      homo_list <- lapply(1:Day, function(g) {homo.day.est(g)})
      homo_result <- do.call(rbind, homo_list)%>%as.data.frame
      homo_result$index <- 'homo'
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list)%>%as.data.frame
      exp_result$index <- 'exp'
      lower_bound<<-0
      upper_bound<<-100
      power_list <- lapply(1:Day, function(g) {power.day.est(g)})
      power_result <- do.call(rbind, power_list)%>%as.data.frame
      power_result$index <- 'power'
      result <- rbind(homo_result, power_result, exp_result)
      #----------rho, error, convergence----------#
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      folder_path <- file.path(datapath, "urban", "params", Dname)
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      write.csv(SLDR, file=file.path(folder_path, paste(county.id, "_SLDR_params_", Yname[yindex], ".csv", sep = "")), row.names = FALSE)
      
      #----------2.fitting----------#
      SLDR.fit<-NULL
      SLDR.homo <- subset(SLDR,index=='homo')
      SLDR.power <- subset(SLDR,index=='power')
      SLDR.exp <- subset(SLDR,index=='exp')
      for(i in 1:Day){
        # empirical
        s1 <- read.csv(paste(flowpath,"/msa/", Dname,"/Intra_Flow_", Datevalue[i], ".csv", sep = ""), header=TRUE)
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
      folder_path <- file.path(datapath, "urban", "fit", Dname)
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      write.csv(SLDR.fit, file=file.path(folder_path, paste(county.id, "_SLDR_fit_", Yname[yindex], ".csv", sep = "")), row.names = FALSE)
      
    } # county
    print(d)
  } # dis
}# yindex




#----------SLDR fitting results for all msa_county----------#
for(yindex in 1:Ynum){
  
  dis.rho <- NULL
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    date.gif<-date.sep(Datevalue,durdate)
    county.msa <- subset(county,msa.id==d)
    Ncounty <- nrow(county.msa)
    
    for(f in 1:Ncounty){
      county.id <- paste(county.msa$state.fips[f],county.msa$county.fips[f],sep='')
      SLDR <- read.csv(paste(datapath,"/urban/params/", Dname,"/",county.id, "_SLDR_params_", Yname[yindex], ".csv", sep = ""), header=TRUE)
      SLDR$county <- county.id
      SLDR$msa <- Dname
      SLDR$day <- as.Date(SLDR$day)
      SLDR <- subset(SLDR, day%in%Datevalue)
      SLDR <- left_join(SLDR,date.gif,by='day')
      dis.rho <- plyr::rbind.fill(dis.rho,SLDR)
    } # county
    
  } # msa
  write.csv(dis.rho, file = paste(datapath,"/urban/SLDR_params_",Yname[yindex],".csv",sep=""), row.names = FALSE)
} # yindex



