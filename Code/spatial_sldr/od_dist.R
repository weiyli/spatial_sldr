# rm(list = ls())

# calculate the distance (queen or other) of OD flows


#----------Workpath----------#
setwd("/home/weiy.li/")
codepath <- '/home/weiy.li/Code/spatial_sldr'
geopath <- '/home/weiy.li/Data/Geo'
#----------Part1: OD flow for COVID-19----------#
flowpath <- '/work/rwuirlab/Xinhua/WeiyuLi/Data/Flow'


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
datapath <- '/home/weiy.li/Data/Flow'
#----------Mobility----------#
for(i in 1:Day){
  
  # home (census) cbg -> cbg
  df <- read.csv(paste(flowpath,"/", year(Datevalue[i]),"/", PDatevalue[i],"/",Datevalue[i], "-social-distancing.csv.gz",sep=""))
  OD<-OD.fun(orig_cbgs=df$origin_census_block_group, dest_cbgs=df$destination_cbgs)
  
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
    BlockID <- unique(block.msa$CensusBlockGroup)
    # OD <- data.frame(from=rep(BlockID,NBlock), to=rep(BlockID,each=NBlock))
    unique(block.msa$State)
    
    #----------compute adjacency & distance matrices----------#
    # haversine distance
    coords <- sf::st_centroid(block.msa$geometry)%>%sf::st_coordinates()
    dist_haversine <- geosphere::distm(coords, fun = geosphere::distHaversine) 
    # queen distance
    Wd <- spdep::poly2nb(block.msa, queen = TRUE)  
    adj_matrix <- nb.mat(nb = Wd)  
    g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "directed")
    dist_queen <- igraph::distances(g)  
    
    for(s in 1:Nregion){
      if(region[s]=="msa"){
        
        #----------OD ----------#
        OD1 <- OD[OD$orig_bg %in% BlockID & OD$dest_bg %in% BlockID, ]
        OD.FTH <- OD1[OD1$no_visits>FTH,]
        orig_indices <- match(OD.FTH$orig_bg, BlockID)
        dest_indices <- match(OD.FTH$dest_bg, BlockID)
        OD.FTH$haversine <- dist_haversine[cbind(orig_indices, dest_indices)]
        OD.FTH$queen <- dist_queen[cbind(orig_indices, dest_indices)]
        
        #----------distance----------#
        result <- as.data.table(OD.FTH)[, .(
          weighted_haversine = sum(haversine * no_visits) / sum(no_visits),
          weighted_queen     = sum(queen * no_visits) / sum(no_visits),
          total_visits       = sum(no_visits),
          n_destinations     = .N
        ), by = orig_bg]
        colnames(result)[1]<-"CensusBlockGroup"
        dist.msa <- left_join(data.frame(CensusBlockGroup=BlockID), result, by="CensusBlockGroup")
        
        output_dir <- paste(datapath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/dist_", Datevalue[i], ".csv", sep="")
        write.csv(dist.msa, file=output_file, row.names=FALSE)

      }else{
        
        county.msa.all <- subset(county.info,msa_id==d)
        dist.county<-NULL
        for(f in 1:nrow(county.msa.all)){
          # keep the od with BlockID in the county and the no_visits>FTH(0)
          county.msa <- subset(block.msa,county_fips ==county.msa.all$county_fips[f])
          CountyID <- unique(county.msa$CensusBlockGroup)
          OD1 <- OD[OD$orig_bg %in% CountyID & OD$dest_bg %in% CountyID, ]
          OD.FTH <- OD1[OD1$no_visits>FTH,]
          orig_indices <- match(OD.FTH$orig_bg, BlockID)
          dest_indices <- match(OD.FTH$dest_bg, BlockID)
          OD.FTH$haversine <- dist_haversine[cbind(orig_indices, dest_indices)]
          OD.FTH$queen <- dist_queen[cbind(orig_indices, dest_indices)]
          
          #----------distance----------#
          result <- as.data.table(OD.FTH)[, .(
            weighted_haversine = sum(haversine * no_visits) / sum(no_visits),
            weighted_queen     = sum(queen * no_visits) / sum(no_visits),
            total_visits       = sum(no_visits),
            n_destinations     = .N
          ), by = orig_bg]
          colnames(result)[1]<-"CensusBlockGroup"
          result <- left_join(data.frame(CensusBlockGroup=CountyID), result, by="CensusBlockGroup")
          result$county_fips<-county.msa.all$county_fips[f]
          result$day <- Datevalue[i]
          dist.county <- plyr::rbind.fill(dist.county, result)
        } # county
        
        # save od dist 
        output_dir <- paste(datapath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/dist_", Datevalue[i], ".csv", sep="")
        write.csv(dist.county, file=output_file, row.names=FALSE)
        
      } 
    } # region: msa or county
  } # dis
  print(i)
} # day











#----------Part2: disaster----------#
datapath <- '/home/weiy.li/Data/Flow/disaster'
#----------Mobility----------#
for(d in c(Nmsa+(1:Ndis))){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
  BlockID<-unique(block.msa$CensusBlockGroup)
  
  for(i in 1:Day){
    
    # home (census) cbg -> cbg
    if(d==Nmsa+1){
      OD <- read.csv(paste(flowpath, "/2017 Harvey flow/",Datevalue[i],".csv",sep=""),header=FALSE)
      names(OD)<-c("orig_bg","dest_bg","no_visits")
    }else{
      df <- read.csv(paste(flowpath, "/", year(Datevalue[i]),"/", PDatevalue[i],"/",Datevalue[i], "-social-distancing.csv.gz",sep=""))
      OD<-OD.fun(orig_cbgs=df$origin_census_block_group, dest_cbgs=df$destination_cbgs)
    }
    
    for(s in 1:Nregion){
      if(region[s]=="msa"){
        # keep the od with BlockID in the msa and the no_visits>FTH(0)
        OD1 <- OD[OD$orig_bg %in% BlockID & OD$dest_bg %in% BlockID, ]
        OD.FTH <- OD1[OD1$no_visits>FTH,]
        
        #----------Part1: Intra_Flow, calculate the inflow and outflow using flow.fun(): TI,TO----------#
        # OD include OO
        result<-flow.fun(ID=BlockID,OD=OD.FTH)
        colnames(result)[1]<-"CensusBlockGroup"
        result[is.na(result)] <- 0
        OD_in_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                             TI_in_OO=result$TI,
                             TO_in_OO=result$TO)
        # OD exclude OO
        result<-flow.fun(ID=BlockID,OD=OD.FTH[OD.FTH$orig_bg != OD.FTH$dest_bg, ])
        colnames(result)[1]<-"CensusBlockGroup"
        result[is.na(result)] <- 0
        OD_ex_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                             TI_ex_OO=result$TI,
                             TO_ex_OO=result$TO)
        result<-left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
        # save od flows in each month
        output_dir <- paste(datapath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(result, file=output_file, row.names=FALSE)
        
        #----------Part2: Intra_Degree, calculate the indegree and outdree using degree.fun(): TI,TO----------#
        # OD include OO
        result<-degree.fun(ID=BlockID,OD=OD.FTH)
        colnames(result)[1]<-"CensusBlockGroup"
        result[is.na(result)] <- 0
        OD_in_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                             TI_in_OO=result$TI,
                             TO_in_OO=result$TO)
        # OD exclude OO
        result<-degree.fun(ID=BlockID,OD=OD.FTH[OD.FTH$orig_bg != OD.FTH$dest_bg, ])
        colnames(result)[1]<-"CensusBlockGroup"
        result[is.na(result)] <- 0
        OD_ex_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                             TI_ex_OO=result$TI,
                             TO_ex_OO=result$TO)
        result<-left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
        # save od degree 
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(result, file=output_file, row.names=FALSE)
        
      }else{
        
        county.msa.all <- subset(county.info,msa_id==d)
        flow.county <- degree.county <-NULL
        for(f in 1:nrow(county.msa.all)){
          # keep the od with BlockID in the county and the no_visits>FTH(0)
          county.msa <- subset(block.msa,county_fips ==county.msa.all$county_fips[f])
          CountyID <- unique(county.msa$CensusBlockGroup)
          OD1 <- OD[OD$orig_bg %in% CountyID & OD$dest_bg %in% CountyID, ]
          OD.FTH <- OD1[OD1$no_visits>FTH,]
          
          #----------Part1: Intra_Flow, calculate the inflow and outflow using flow.fun(): TI,TO----------#
          # OD include OO
          result<-flow.fun(ID=CountyID,OD=OD.FTH)
          colnames(result)[1]<-"CensusBlockGroup"
          result[is.na(result)] <- 0
          OD_in_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                               TI_in_OO=result$TI,
                               TO_in_OO=result$TO)
          # OD exclude OO
          result<-flow.fun(ID=CountyID,OD=OD.FTH[OD.FTH$orig_bg != OD.FTH$dest_bg, ])
          colnames(result)[1]<-"CensusBlockGroup"
          result[is.na(result)] <- 0
          OD_ex_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                               TI_ex_OO=result$TI,
                               TO_ex_OO=result$TO)
          result<-left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
          result$county_fips<-county.msa.all$county_fips[f]
          flow.county <- plyr::rbind.fill(flow.county, result)
          
          #----------Part2: Intra_Degree, calculate the indegree and outdree using degree.fun(): TI,TO----------#
          # OD include OO
          result<-degree.fun(ID=CountyID,OD=OD.FTH)
          colnames(result)[1]<-"CensusBlockGroup"
          result[is.na(result)] <- 0
          OD_in_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                               TI_in_OO=result$TI,
                               TO_in_OO=result$TO)
          # OD exclude OO
          result<-degree.fun(ID=CountyID,OD=OD.FTH[OD.FTH$orig_bg != OD.FTH$dest_bg, ])
          colnames(result)[1]<-"CensusBlockGroup"
          result[is.na(result)] <- 0
          OD_ex_OO<-data.frame(CensusBlockGroup=result$CensusBlockGroup,
                               TI_ex_OO=result$TI,
                               TO_ex_OO=result$TO)
          result<-left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
          result$county_fips<-county.msa.all$county_fips[f]
          degree.county <- plyr::rbind.fill(degree.county, result)
        } # county
        
        # save od flows
        output_dir <- paste(datapath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(flow.county, file=output_file, row.names=FALSE)
        
        # save od degree
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(degree.county, file=output_file, row.names=FALSE)
      } 
    } # region: msa or county
  } # day 
  print(d)
} # dis
