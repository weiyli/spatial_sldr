# rm(list = ls())

# SLDR fitting at the MSA level for estimating rho over a 10-week period.


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'


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
dset <- 1:Nmsa

#----------SLDR fitting for queen----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in dset){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      Datevalue <- seq(as.Date('2020-01-01'), as.Date('2020-06-30'), by="day")
      Day <- length(Datevalue)
      
      #----------BlockID, NBlock----------#
      block.msa <- sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID <- unique(block.msa$CensusBlockGroup)
      NBlock <-length(BlockID)
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
      lower_power<<-40
      upper_power<<-80
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
      SLDR <- data.frame(day = Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/params/",Dname,"_SLDR_params_",Yname[yindex],"_half_year.csv",sep=""), row.names = FALSE)
      
      
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
      write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fit/",Dname,"_SLDR_fit_",Yname[yindex],"_half_year.csv",sep=""),row.names = FALSE)
      
      print(d)
    } # dis
  }# yindex
}# region





#----------SLDR fitting results for all msa_county----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    dis.rho <- NULL
    for(d in 1:Nmsa){
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # new date: 5 weeks before + 5 weeks during
      duration <- 7*4 + 6
      startdate <- '2020-01-13'
      durdate <- '2020-03-30'
      enddate <- '2020-05-03'
      Datevalue <- c(seq(as.Date(startdate),as.Date(startdate) + duration, by="day"), seq(as.Date(durdate),as.Date(durdate) + duration, by="day"))
      Day <- length(Datevalue) 
      date.gif <- date.sep(Datevalue,durdate)
      #----------params----------#
      sldr.params <- read.csv(paste(datapath,"/",region[s],"/params/",Dname,"_SLDR_params_",Yname[yindex],"_half_year.csv",sep=""), header=TRUE)
      sldr.params$event <- event
      sldr.params$msa <- Dname
      sldr.params$day <- as.Date(sldr.params$day)
      sldr.params <- left_join(sldr.params, date.gif, by='day')
      dis.rho <- plyr::rbind.fill(dis.rho, sldr.params)
    } # dis
    write.csv(dis.rho, file=paste(datapath,"/",region[s],"/SLDR_params_",Yname[yindex],"_half_year.csv",sep=""), row.names = FALSE)
  } # yindex
} # region


