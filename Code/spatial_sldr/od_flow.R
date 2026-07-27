# rm(list = ls())

# calculate the intra-od flow

#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

# #----------unzip the raw mobility data----------#
# year<-c(2019,2020,2021)
# for(i in 2:length(year)){
#   zip_file <- paste("/work/rwuirlab/Safegraph/Safegraph_Daily/", year[i],".zip",sep="")
#   target_file <- paste(datapath,year[i],sep="")
#   if (!dir.exists(target_file)) {
#     dir.create(target_file, recursive = TRUE)
#   }
#   unzip(zip_file, exdir = target_file)
# }

# #----------unzip the raw mobility data----------#
# year<-c(2019,2020,2021)
# for(i in 2:length(year)){
#   zip_file <- paste(flowpath,"/", year[i], ".zip", sep="")
#   target_file <- paste(flowpath, "/", year[i], sep="")
#   if (!dir.exists(target_file)) {
#     dir.create(target_file, recursive = TRUE)
#   }
#   unzip(zip_file, exdir = target_file)
# }


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

#----------tract in each msa----------#
tract.info <- read.csv(paste(geopath,"/msa/selected_msa_tract.csv",sep=""), header=TRUE)
tract.info$tract_fips <- sprintf("%011s", as.numeric(tract.info$tract_fips))

#----------Part1: COVID-19----------#
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'

#----------Mobility: CBG_MSA and CBG_County----------#
for(i in 1:Day){
  
  # home (census) cbg -> cbg
  df <- read.csv(paste(flowpath,"/", year(Datevalue[i]),"/", PDatevalue[i],"/",Datevalue[i], "-social-distancing.csv.gz",sep=""))
  OD <- OD.fun(orig_cbgs=df$origin_census_block_group, dest_cbgs=df$destination_cbgs)

  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    block.msa <- sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
    BlockID <- unique(block.msa$CensusBlockGroup)

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
        flow.msa <- left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
        output_dir <- paste(flowpath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(flow.msa, file=output_file, row.names=FALSE)
        
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
        degree.msa <- left_join(OD_in_OO,OD_ex_OO,by="CensusBlockGroup")
        # save od degree in each month
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(degree.msa, file=output_file, row.names=FALSE)
        
      }else{
        
        county.msa.all <- subset(county.info,msa_id==d)
        flow.county<-degree.county<-NULL
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
          result$day <- Datevalue[i]
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
        output_dir <- paste(flowpath, region[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(flow.county, file=output_file, row.names=FALSE)
        
        # save od degree in each month
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(degree.county, file=output_file, row.names=FALSE)
      } 
    } # region: msa or county
  } # dis
  print(i)
} # day


#----------Mobility data for sldr_fit_MAUP.R: County_MSA and Tract_MSA----------#
for(i in 1:Day){
  
  # home (census) cbg -> cbg
  df <- read.csv(paste(flowpath,"/", year(Datevalue[i]),"/", PDatevalue[i],"/",Datevalue[i], "-social-distancing.csv.gz",sep=""))
  OD <- OD.fun(orig_cbgs=df$origin_census_block_group, dest_cbgs=df$destination_cbgs)
  
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
    BlockID<-unique(block.msa$CensusBlockGroup)
    
    for(s in 1:length(region_add)){
      if(region_add[s]=="tract_msa"){
        
        block.msa$tract_fips <- stringr::str_sub(sprintf("%012s", as.numeric(block.msa$CensusBlockGroup)), 1, -2)
        #----------tract list within this MSA (from selected_msa_tract.csv / selected.tract)----------#
        tract.msa.all <- subset(tract.info, msa_id == d)  
        TractID <- sprintf("%011s", as.numeric(unique(tract.msa.all$tract_fips)))%>%unique()
        tract.msa <- subset(block.msa, tract_fips%in%TractID)
        cbg.tract <- unique(tract.msa$CensusBlockGroup)
        OD1 <- OD[OD$orig_bg %in%  cbg.tract & OD$dest_bg %in% cbg.tract, ]
        OD.FTH <- OD1[OD1$no_visits>FTH,]
        
        cbg2tract <- tract.msa %>% sf::st_drop_geometry() %>% dplyr::select(CensusBlockGroup, tract_fips) %>% dplyr::distinct()
        cbg2tract$tract_fips <- sprintf("%011s", as.numeric(cbg2tract$tract_fips))
        OD.tract <- dplyr::left_join(OD.FTH, cbg2tract, by = c("orig_bg" = "CensusBlockGroup")) %>% dplyr::rename(orig_tract = tract_fips)
        OD.tract <- dplyr::left_join(OD.tract, cbg2tract, by = c("dest_bg" = "CensusBlockGroup")) %>% dplyr::rename(dest_tract = tract_fips)
        OD.tract <- OD.tract %>% dplyr::filter(!is.na(orig_tract) & !is.na(dest_tract))
        OD.tract <- OD.tract %>% dplyr::group_by(orig_tract, dest_tract) %>% dplyr::summarise(no_visits = sum(no_visits, na.rm = TRUE), .groups="drop")
        
        #----------Part1: Intra_Flow (tract)----------#
        result <- flow.fun.tract(ID = TractID, OD = OD.tract)
        colnames(result)[1] <- "tract_fips"
        result[is.na(result)] <- 0
        OD_in_OO <- data.frame(tract_fips = result$tract_fips,
                               TI_in_OO = result$TI,
                               TO_in_OO = result$TO)
        
        result <- flow.fun.tract(ID = TractID, OD = OD.tract[OD.tract$orig_tract != OD.tract$dest_tract, ])
        colnames(result)[1] <- "tract_fips"
        result[is.na(result)] <- 0
        OD_ex_OO <- data.frame(tract_fips = result$tract_fips,
                               TI_ex_OO = result$TI,
                               TO_ex_OO = result$TO)
        
        flow.tract.msa <- left_join(OD_in_OO, OD_ex_OO, by = "tract_fips")
        
        output_dir <- paste(flowpath, region_add[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(flow.tract.msa, file = output_file, row.names = FALSE)
        
        #----------Part2: Intra_Degree (tract)----------#
        result <- degree.fun.tract(ID = TractID, OD = OD.tract)
        colnames(result)[1] <- "tract_fips"
        result[is.na(result)] <- 0
        OD_in_OO <- data.frame(tract_fips = result$tract_fips,
                               TI_in_OO = result$TI,
                               TO_in_OO = result$TO)
        
        result <- degree.fun.tract(ID = TractID, OD = OD.tract[OD.tract$orig_tract != OD.tract$dest_tract, ])
        colnames(result)[1] <- "tract_fips"
        result[is.na(result)] <- 0
        OD_ex_OO <- data.frame(tract_fips = result$tract_fips,
                               TI_ex_OO = result$TI,
                               TO_ex_OO = result$TO)
        
        degree.tract.msa <- left_join(OD_in_OO, OD_ex_OO, by = "tract_fips")
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(degree.tract.msa, file = output_file, row.names = FALSE)
        
      } else if(region_add[s]=="county_msa"){
        
        #----------keep the od with BlockID in the county and the no_visits>FTH(0)----------#
        county.msa.all <- subset(county.info,msa_id==d)
        CountyID <- sprintf("%05d", as.numeric(unique(county.msa.all$county_fips)))%>%unique
        
        county.msa <- subset(block.msa, county_fips%in%CountyID)
        cbg.county <- unique(county.msa$CensusBlockGroup)
        OD1 <- OD[OD$orig_bg %in%  cbg.county & OD$dest_bg %in%  cbg.county, ]
        OD.FTH <- OD1[OD1$no_visits>FTH,]
        cbg2county <- county.msa %>% sf::st_drop_geometry() %>% dplyr::select(CensusBlockGroup, county_fips) %>% dplyr::distinct()
        cbg2county$county_fips <- sprintf("%05d", as.numeric(cbg2county$county_fips))
        OD.county <- dplyr::left_join(OD.FTH, cbg2county, by = c("orig_bg" = "CensusBlockGroup")) %>% dplyr::rename(orig_county = county_fips)
        OD.county <- dplyr::left_join(OD.county, cbg2county, by = c("dest_bg" = "CensusBlockGroup")) %>% dplyr::rename(dest_county = county_fips)
        OD.county <- OD.county %>% dplyr::filter(!is.na(orig_county) & !is.na(dest_county))
        OD.county <- OD.county %>% dplyr::group_by(orig_county, dest_county) %>% dplyr::summarise(no_visits = sum(no_visits, na.rm = TRUE), .groups="drop")
        
        #----------Part1: Intra_Flow, calculate the inflow and outflow using flow.fun.county(): TI,TO----------#
        # OD include OO
        result <- flow.fun.county(ID = CountyID, OD = OD.county)
        colnames(result)[1] <- "county_fips"
        result[is.na(result)] <- 0
        OD_in_OO<-data.frame(county_fips=result$county_fips,
                             TI_in_OO=result$TI,
                             TO_in_OO=result$TO)
        # OD exclude OO
        result <- flow.fun.county(ID=CountyID, OD=OD.county[OD.county$orig_county != OD.county$dest_county, ])
        colnames(result)[1] <- "county_fips"
        result[is.na(result)] <- 0
        OD_ex_OO<-data.frame(county_fips=result$county_fips,
                             TI_ex_OO=result$TI,
                             TO_ex_OO=result$TO)
        flow.county.msa <- left_join(OD_in_OO, OD_ex_OO, by="county_fips")
        output_dir <- paste(flowpath, region_add[s], Dname, sep="/")
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        output_file <- paste(output_dir, "/Intra_Flow_", Datevalue[i], ".csv", sep="")
        write.csv(flow.county.msa, file=output_file, row.names=FALSE)
        
        #----------Part2: Intra_Degree, calculate the indegree and outdree using degree.fun.county(): TI,TO----------#
        # OD include OO
        result <- degree.fun.county(ID = CountyID, OD = OD.county)
        colnames(result)[1] <- "county_fips"
        result[is.na(result)] <- 0
        OD_in_OO <- data.frame(county_fips = result$county_fips,
                               TI_in_OO = result$TI,
                               TO_in_OO = result$TO)
        # OD exclude OO
        result <- degree.fun.county(ID = CountyID,
                                    OD = OD.county[OD.county$orig_county != OD.county$dest_county, ])
        colnames(result)[1] <- "county_fips"
        result[is.na(result)] <- 0
        OD_ex_OO <- data.frame(county_fips = result$county_fips,
                               TI_ex_OO = result$TI,
                               TO_ex_OO = result$TO)
        degree.county.msa <- left_join(OD_in_OO, OD_ex_OO, by = "county_fips")
        output_file <- paste(output_dir, "/Intra_Degree_", Datevalue[i], ".csv", sep="")
        write.csv(degree.county.msa, file = output_file, row.names = FALSE)
      }
    
    } # region: county_msa or tract_msa
  } # dis
  print(i)
} # day



#----------Part2: disaster----------#
flowpath <- 'D:/ood/Data/Flow/disaster'
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



