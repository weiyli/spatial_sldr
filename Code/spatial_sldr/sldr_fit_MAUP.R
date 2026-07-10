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


#----------Load tract/county polygons (for fitting at those scales)----------#
msa.boundary.all   <- sf::read_sf(paste(geopath,"/census/tigris_msa_boundary_2010_2019.geojson",sep=""))
county.boundary.all<- sf::read_sf(paste(geopath,"/census/tigris_county_boundary_2010_2019.geojson",sep=""))
tract.boundary.all <- sf::read_sf(paste(geopath,"/census/tigris_tract_boundary_2010_2019.geojson",sep=""))

#----------Load msa-county / msa-tract relationship tables----------#
county.info <- read.csv(paste(geopath,"/msa/selected_msa_county.csv",sep=""), header=TRUE)
county.info$county_fips <- sprintf("%05d", as.numeric(county.info$county_fips))

tract.info <- read.csv(paste(geopath,"/msa/selected_msa_tract.csv",sep=""), header=TRUE)
tract.info$tract_fips <- sprintf("%011s", as.numeric(tract.info$tract_fips)) 

#----------set the global function pointer used by homo/power/exp.param.est----------#
error_test <<- error_specs[["rmse_rmspe"]]$fun



#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'

#----------Create output dirs----------#
for(rg in region_add){
  dir.create(paste0(datapath,"/",rg,"/param"), recursive=TRUE, showWarnings=FALSE)
  dir.create(paste0(datapath,"/",rg,"/fitting"),    recursive=TRUE, showWarnings=FALSE)
}

#----------SLDR fitting: county_msa or tract_msa----------#
for(rg in region_add){
  for(yindex in 2:2){
    for(d in 1:Nmsa){
            
      Dname <<- Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
      
      #====================================================#
      #  A) Load spatial units (polygons) at the target scale
      #====================================================#
      if(rg == "tract_msa"){
        # tracts within this MSA (polygon units)
        block.msa$tract_fips <- stringr::str_sub(sprintf("%012s", as.numeric(block.msa$CensusBlockGroup)), 1, -2)
        unit.sf <- block.msa %>% group_by(tract_fips) %>% summarise(geometry = sf::st_union(geometry),.groups = "drop")
        unit.sf <- unit.sf %>% sf::st_make_valid() %>% sf::st_cast("MULTIPOLYGON", warn = FALSE)
        NodeID <- unique(unit.sf$tract_fips)
      } else if(rg == "county_msa"){
        # counties within this MSA (polygon units)
        block.msa$county_fips <- sprintf("%05d", as.numeric(block.msa$county_fips))
        unit.sf <- block.msa %>% group_by(county_fips) %>% summarise(geometry = sf::st_union(geometry),.groups = "drop")
        unit.sf <- unit.sf %>% sf::st_make_valid() %>% sf::st_cast("MULTIPOLYGON", warn = FALSE)
        NodeID <- unique(unit.sf$county_fips)
      }
      # if(rg == "county_msa"){
      #   # counties within this MSA (polygon units)
      #   msa.boundary <- msa.boundary.all[msa.boundary.all$msa_fips == Dname.id, ]
      #   unit.sf <- sf::st_intersection(county.boundary.all, msa.boundary)
      #   unit.sf$county_fips <- sprintf("%05d", as.numeric(unit.sf$county_fips))
      #   NodeID <- unique(unit.sf$county_fips)
      #   
      # } else if(rg == "tract_msa"){
      #   # tracts within this MSA (polygon units)
      #   msa.boundary <- msa.boundary.all[msa.boundary.all$msa_fips == Dname.id, ]
      #   unit.sf <- sf::st_intersection(tract.boundary.all, msa.boundary)
      #   unit.sf$tract_fips <- sprintf("%011s", as.character(unit.sf$tract_fips))
      #   NodeID <- unique(unit.sf$tract_fips)
      # }
      BlockID <<- NodeID
      NBlock <<- length(BlockID)
      if (is.na(NBlock) || NBlock < 5) {
        message("Skip: ", rg, " | ", Dname, " (NBlock=", NBlock, ")")
        next  # go to next d
      }
      
      #====================================================#
      #  B) Build W (queen contiguity) and lag shells
      #====================================================#
      Wd <- spdep::poly2nb(unit.sf, queen = TRUE)
      W1 <- nb.mat(nb = Wd)
      lag.max <- igraph::graph_from_adjacency_matrix(W1, mode = "directed") %>% igraph::diameter()
      W <<- queen.adj(W1, lag.max)
      slag <<- seq(1, lag.max, 1)
      
      #====================================================#
      #  C) Load mobility series Mob (node-by-day)
      #====================================================#
      Mob <<- matrix(0, nrow = NBlock, ncol = Day)
      for(i in 1:Day){
        # flow data from od_flow.R
        s1 <- read.csv(paste(flowpath,"/", rg, "/", Dname, "/Intra_Flow_", Datevalue[i], ".csv", sep=""), header=TRUE)
        if (rg == "tract_msa") {
          s1$tract_fips <- sprintf("%011s", as.numeric(s1$tract_fips))
          s2 <- dplyr::left_join(data.frame(ID = BlockID),
                                 data.frame(ID = s1$tract_fips, mob = s1[, yindex+1]),
                                 by = "ID")
          s2$mob[is.na(s2$mob)] <- 0
          Mob[, i] <- s2$mob
        } else if (rg == "county_msa") {
          s1$county_fips <- sprintf("%05s",  as.numeric(s1$county_fips))
          s2 <- dplyr::left_join(data.frame(ID = BlockID),
                                 data.frame(ID = s1$county_fips, mob = s1[, yindex+1]),
                                 by = "ID")
          s2$mob[is.na(s2$mob)] <- 0
          Mob[, i] <- s2$mob
        }
      } # day
      
      #====================================================#
      #  D) Estimate rho day-by-day (your existing functions)
      #====================================================#
      lower_bound <<- 0
      upper_bound <<- 5
      lower_power <<- 40
      upper_power <<- 80
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list) %>% as.data.frame
      exp_result$index <- "exp"
      result <- exp_result
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR,file = paste(datapath,"/", rg, "/param/", Dname, "_SLDR_params_", Yname[yindex], ".csv", sep=""),row.names = FALSE)

      #====================================================#
      #  E) Fitting series (empirical vs simulated)
      #====================================================#
      SLDR.fit <- NULL
      SLDR.exp   <- subset(SLDR, index=="exp")
      for(i in 1:Day){
        # empirical M.day
        s1 <- read.csv(paste(flowpath,"/", rg, "/", Dname, "/Intra_Flow_", Datevalue[i], ".csv", sep=""), header=TRUE)
        if (rg == "tract_msa") {
          s1$tract_fips <- sprintf("%011s", as.numeric(s1$tract_fips))
          s2 <- dplyr::left_join(data.frame(ID = BlockID),
                                 data.frame(ID = s1$tract_fips, mob = s1[, yindex+1]),
                                 by = "ID")
          s2$mob[is.na(s2$mob)] <- 0
          M.day <- s2$mob
        } else if (rg == "county_msa") {
          s1$county_fips <- sprintf("%05s",  as.numeric(s1$county_fips))
          s2 <- dplyr::left_join(data.frame(ID = BlockID),
                                 data.frame(ID = s1$county_fips, mob = s1[, yindex+1]),
                                 by = "ID")
          s2$mob[is.na(s2$mob)] <- 0
          M.day <- s2$mob
        }
        # exp
        rho <- SLDR.exp$rho[i]
        lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol=length(slag), byrow=TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.exp <- rowSums(lag.h * WM.h)
        SLDR.fit <- plyr::rbind.fill(SLDR.fit, data.frame(day=Datevalue[i], 
                                                          empirical=M.day,
                                                          exp=M.day.sim.exp))

      } # day
      
      write.csv(SLDR.fit, file = paste(datapath,"/", rg, "/fitting/", Dname, "_SLDR_fit_", Yname[yindex], ".csv", sep=""), row.names = FALSE)
      message("Finished: ", rg, " | ", Dname)
      
    } # d
  } # yindex
} # region_add





#----------SLDR fitting for rook adjacency matrices----------#
for(s in 2:Nregion){
  for(yindex in 2:2){
    for(d in 1:Nmsa){
    # for(d in 8:8){
      
      Dname<<-Dname.set[d]
      source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
      # date.gif<-date.sep(Datevalue,durdate)
      
      #----------BlockID, NBlock----------#
      block.msa<-sf::read_sf(paste(geopath,"/",region[s],"/",Dname,".geojson",sep=""))
      BlockID<-unique(block.msa$CensusBlockGroup)
      NBlock<-length(BlockID)
      #----------W: 1 to lag.max lag rook adjacency matrices----------#
      Wd <- spdep::poly2nb(block.msa, queen = FALSE) 
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
      exp_list <- lapply(1:Day, function(g) {exp.day.est(g)})
      exp_result <- do.call(rbind, exp_list)%>%as.data.frame
      exp_result$index <- 'exp'
      result <- exp_result
      #----------rho, error, convergence----------#
      SLDR <- data.frame(day=Datevalue, lag=lag.max, rho=result[,1], error=result[,2], index=result[,3])
      write.csv(SLDR, file=paste(datapath,"/",region[s],"/params/",Dname,"_SLDR_params_rook_",Yname[yindex],".csv",sep=""), row.names = FALSE)
      
      
      #----------2.fitting----------#
      SLDR.fit<-NULL
      SLDR.exp <- subset(SLDR,index=='exp')
      for(i in 1:Day){
        # empirical
        s1 <- read.csv(paste(flowpath,"/", region[s], "/", Dname, "/Intra_Flow_",Datevalue[i],".csv", sep = ""), header=TRUE)
        s2 <- left_join(data.frame(CensusBlockGroup=BlockID),
                        data.frame(CensusBlockGroup=s1$CensusBlockGroup, mob=s1[,yindex+1]),
                        by="CensusBlockGroup")
        s2$mob[is.na(s2$mob)] <- 0
        M.day <- s2$mob
        # exp
        rho<-SLDR.exp$rho[i]
        lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol = length(slag), byrow = TRUE)
        WM.h <- sapply(W, "%*%", M.day)
        M.day.sim.exp <- rowSums(lag.h*WM.h)
        SLDR.fit <- plyr::rbind.fill(SLDR.fit, data.frame(day=Datevalue[i], 
                                                          empirical=M.day, 
                                                          exp=M.day.sim.exp))
      } # day
      write.csv(SLDR.fit,file=paste(datapath,"/",region[s],"/fit/",Dname,"_SLDR_fit_rook_",Yname[yindex],".csv",sep=""),row.names = FALSE)
      
      print(d)
    } # dis
  }# yindex
}# region


