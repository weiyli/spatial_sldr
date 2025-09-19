# rm(list = ls())

# select the region: msa and county

#----------Workpath----------#
setwd("/home/weiy.li/")
codepath <- '/home/weiy.li/Code/spatial_sldr'
datapath <- '/home/weiy.li/Data/spatial_sldr'
figpath <- '/home/weiy.li/Figure/spatial_sldr'
geopath <- '/home/weiy.li/Data/Geo'


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'

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

#----------Load packages----------#
library(sf)          # poly2nb()

#----------demographic data: race and income----------#
# data ref (acs_7_categories_geodata.R): ACS 2015-2019 5 years data, with 1055 cbgs with 0 pop and NA income, the total NA income is 8065
demo <- read.csv(paste(geopath,"/census/acs_2019_5years_cbg.csv",sep=""), header=TRUE)
dis.demo <- setDT(demo)[, .(
  CensusBlockGroup = cbg_fips,
  pop = total_population,
  p_black = black_population / total_population,
  p_white = white_population / total_population,
  p_asian = asian_population / total_population,
  p_hispanic = hispanic_population / total_population,
  median_household_income = median_household_income
)]
# pop = 0 -> NA
dis.demo[pop == 0, `:=`(
  pop = NA,
  p_black = NA,
  p_white = NA,
  p_asian = NA,
  p_hispanic = NA)]

#----------geo data----------#
# msa, county, and cbg from 2010_2019 census
msa <- sf::read_sf(paste(geopath,"/census/tigris_msa_boundary_2010_2019.geojson",sep=""))
county <- sf::read_sf(paste(geopath,"/census/tigris_county_boundary_2010_2019.geojson",sep=""))

# cbg from safegraph
block <- sf::read_sf(paste(geopath,"/census/cbg.geojson",sep=""))
block$county_fips <- paste0(block$StateFIPS,block$CountyFIPS)
block$CensusBlockGroup<-as.numeric(block$CensusBlockGroup)
BlockID <- unique(block$CensusBlockGroup)
NBlock <- length(BlockID)

selected.county <- county.num <- NULL
for(d in 1:(Nmsa+Ndis)){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  
  #----------counties within msa from 2010_2019 census----------#
  msa.boundary <- msa[which(msa$msa_fips==Dname.id),]
  county.msa <- st_intersection(county, msa.boundary)
  county.info <- data.frame(county_fips=county.msa$county_fips, 
                            county_name=county.msa$county_name, 
                            msa_id=d, 
                            msa_fips=county.msa$msa_fips,
                            msa_name=Dname)%>%unique
  selected.county <- plyr::rbind.fill(selected.county, county.info)
  county.num <- plyr::rbind.fill(county.num,data.frame(county_num=nrow(county.info), msa_name=Dname))
  
  #----------cbgs within counties of msa from 2010_2019 census----------#
  # cbg.msa <- cbg[cbg$county_fips %in% county.info$county_fips,]
  cbg.msa <- block[block$county_fips %in% county.info$county_fips,]
  cbg.msa$area_m2 <- sf::st_area(cbg.msa)%>%as.numeric
  cbg.msa <- left_join(cbg.msa,dis.demo,by="CensusBlockGroup")
  selected.columns <- c("StateFIPS", "CountyFIPS", "county_fips", "CensusBlockGroup", 
                        "State", "County", 
                        "geometry", "pop", "median_household_income", "area_m2")
  cbg.msa <- cbg.msa[, selected.columns]
  
  output.path <- paste(geopath,"/msa/", Dname, ".geojson", sep = "")
  if (!file.exists(output.path)) {
    st_write(cbg.msa, output.path, driver = "GeoJSON")
    message("File saved to: ", output.path)
  } else {
    message("File already exists at: ", output.path)
  }
  
  print(d)
  
}
write.csv(selected.county, file=paste(geopath,"/msa/selected_msa_county.csv",sep=""),row.names = FALSE)





A.panel <- list()
for(d in 1:(Nmsa+Ndis)){
  
  Dname<<-Dname.set[d]
  source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
  
  #----------counties within msa from 2010_2019 census----------#
  msa.boundary <- msa[which(msa$msa_fips==Dname.id),]
  county.msa <- st_intersection(county, msa.boundary)
  A.panel[[d]] <- ggplot() +
    geom_sf(data = county.msa, color = "black") +
    theme_minimal() +
    ggtitle(paste0("County: ",Dname))
  
}
fig.msa.map <- (A.panel[[1]]|A.panel[[2]]|A.panel[[3]]|A.panel[[4]])/
  (A.panel[[5]]|A.panel[[6]]|A.panel[[7]]|A.panel[[8]])/
  (A.panel[[9]]|A.panel[[10]]|A.panel[[11]]|A.panel[[12]])/
  (A.panel[[13]]|A.panel[[14]]|A.panel[[15]]|A.panel[[16]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.msa.map, filename = paste(figpath,"/msa/msa_county_boundary.pdf",sep=""), width =5*4, height = 5*4)




f0 <- ggplot() +
  geom_sf(data = boundary.msa, color = "black") +
  # geom_sf(data = county.msa[which(county.msa$county_fips==51113),], color = "red") +
  theme_minimal() +
  ggtitle(paste0("MSA: ",Dname))

f1 <- ggplot() +
  geom_sf(data = county.msa, color = "black") +
  theme_minimal() +
  ggtitle(paste0("County: ",Dname))

f2 <- ggplot() +
  geom_sf(data = cbg.msa, color = "black") +
  theme_minimal() +
  ggtitle(paste0("CBG: ",Dname))

f0|f1|f2





