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


# #----------county in each msa----------#
# county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
# county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'

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
      lower_power<<-30
      upper_power<<-80
      
      for (es_name in names(error_specs)) {
        
        err_tag <- error_specs[[es_name]]$tag
        
        # set the global function pointer used by homo/power/exp.param.est
        error_test <<- error_specs[[es_name]]$fun
        
        #----------(A) estimate rho under this loss----------#
        homo_list  <- lapply(1:Day, function(g) homo.day.est(g))
        homo_result <- as.data.frame(do.call(rbind, homo_list))
        homo_result$index <- "homo"
        
        power_list  <- lapply(1:Day, function(g) power.day.est(g))
        power_result <- as.data.frame(do.call(rbind, power_list))
        power_result$index <- "power"
        
        exp_list  <- lapply(1:Day, function(g) exp.day.est(g))
        exp_result <- as.data.frame(do.call(rbind, exp_list))
        exp_result$index <- "exp"
        
        result <- rbind(homo_result, power_result, exp_result)
        
        #----------rho, error, convergence----------#
        SLDR <- data.frame(
          day   = Datevalue,
          lag   = lag.max,
          rho   = result[,1],
          error = result[,2],
          index = result[,3],
          err_tag = err_tag
        )
        write.csv(SLDR, file=paste(datapath,"/",region[s],"/param/",Dname,"_SLDR_params_",err_tag,"_",Yname[yindex],".csv",sep=""), row.names = FALSE)
        
        
        #----------(B) Fitting outputs under this loss (use the rho estimated above)----------#
        SLDR.fit <- NULL
        SLDR.homo  <- subset(SLDR, index == "homo")
        SLDR.power <- subset(SLDR, index == "power")
        SLDR.exp   <- subset(SLDR, index == "exp")
        
        for(i in 1:Day){
          s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
          s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                          data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                          by="CensusBlockGroup")
          s2$mob[is.na(s2$mob)] <- 0
          M.day <- s2$mob
          
          # homo
          rho <- SLDR.homo$rho[i]
          lag.h <- matrix(rho, nrow=NBlock, ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.homo <- rowSums(lag.h * WM.h)
          
          # power
          rho <- SLDR.power$rho[i]
          lag.h <- matrix(rep(slag^(-rho), NBlock), ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.power <- rowSums(lag.h * WM.h)
          
          # exp
          rho <- SLDR.exp$rho[i]
          lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.exp <- rowSums(lag.h * WM.h)
          
          SLDR.fit <- plyr::rbind.fill(
            SLDR.fit,
            data.frame(
              day = Datevalue[i],
              empirical = M.day,
              homo = M.day.sim.homo,
              power = M.day.sim.power,
              exp = M.day.sim.exp,
              err_tag = err_tag
            )
          )
        } # day
       
        write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fitting/",Dname,"_SLDR_fit_",err_tag,"_",Yname[yindex],".csv",sep=""),row.names = FALSE)
        
        message("Finished: ", Dname, " | ", err_tag)
      } # err
      
    } # dis
  }# yindex
}# region




#----------Part2: disaster----------#
flowpath <- 'D:/ood/Data/Flow/disaster'
datapath <- 'D:/ood/Data/spatial_sldr/disaster'

#----------SLDR fitting----------#
for(s in 2:Nregion){
  for(yindex in 2:Ynum){
    for(d in Nmsa+(1:Ndis)){
      
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
      lower_power<<-30
      upper_power<<-80
      
      for (es_name in names(error_specs)) {
        
        err_tag <- error_specs[[es_name]]$tag
        
        # set the global function pointer used by homo/power/exp.param.est
        error_test <<- error_specs[[es_name]]$fun
        
        #----------(A) estimate rho under this loss----------#
        homo_list  <- lapply(1:Day, function(g) homo.day.est(g))
        homo_result <- as.data.frame(do.call(rbind, homo_list))
        homo_result$index <- "homo"
        
        power_list  <- lapply(1:Day, function(g) power.day.est(g))
        power_result <- as.data.frame(do.call(rbind, power_list))
        power_result$index <- "power"
        
        exp_list  <- lapply(1:Day, function(g) exp.day.est(g))
        exp_result <- as.data.frame(do.call(rbind, exp_list))
        exp_result$index <- "exp"
        
        result <- rbind(homo_result, power_result, exp_result)
        
        #----------rho, error, convergence----------#
        SLDR <- data.frame(
          day   = Datevalue,
          lag   = lag.max,
          rho   = result[,1],
          error = result[,2],
          index = result[,3],
          err_tag = err_tag
        )
        write.csv(SLDR, file=paste(datapath,"/",region[s],"/param/",event,"_",Dname,"_SLDR_params_",err_tag,"_",Yname[yindex],".csv",sep=""), row.names = FALSE)
        
        
        #----------(B) Fitting outputs under this loss (use the rho estimated above)----------#
        SLDR.fit <- NULL
        SLDR.homo  <- subset(SLDR, index == "homo")
        SLDR.power <- subset(SLDR, index == "power")
        SLDR.exp   <- subset(SLDR, index == "exp")
        
        for(i in 1:Day){
          s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
          s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                          data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                          by="CensusBlockGroup")
          s2$mob[is.na(s2$mob)] <- 0
          M.day <- s2$mob
          
          # homo
          rho <- SLDR.homo$rho[i]
          lag.h <- matrix(rho, nrow=NBlock, ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.homo <- rowSums(lag.h * WM.h)
          
          # power
          rho <- SLDR.power$rho[i]
          lag.h <- matrix(rep(slag^(-rho), NBlock), ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.power <- rowSums(lag.h * WM.h)
          
          # exp
          rho <- SLDR.exp$rho[i]
          lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol=length(slag), byrow=TRUE)
          WM.h <- sapply(W, "%*%", M.day)
          M.day.sim.exp <- rowSums(lag.h * WM.h)
          
          SLDR.fit <- plyr::rbind.fill(
            SLDR.fit,
            data.frame(
              day = Datevalue[i],
              empirical = M.day,
              homo = M.day.sim.homo,
              power = M.day.sim.power,
              exp = M.day.sim.exp,
              err_tag = err_tag
            )
          )
        } # day
        
        write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fitting/",event,"_",Dname,"_SLDR_fit_",err_tag,"_",Yname[yindex],".csv",sep=""),row.names = FALSE)
        
        message("Finished: ", Dname, " | ", err_tag)
      } # err
      
    } # dis
  }# yindex
}# region



