# rm(list = ls())

# obtain the data from the census.


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'


#----------Load packages----------#
library(tidycensus)
library(dplyr)
library(tigris)
library(sf)


# Step 0: Census API and states
census_api_key("0fe438ad561905e0e5050f8b2146396c63af06ba",install = TRUE)
state_county<-tigris::fips_codes
names(state_county)<-c("state_usps","state_fips","state_name","county_fips","county_name")
# write.csv(state_county,paste(geopath,"/us_county_fips.csv",sep=''),row.names = FALSE,quote = TRUE)
states <- unique(tigris::fips_codes$state_code)[1:51]   

#----------Tract boundary----------#
options(tigris_use_cache = TRUE, tigris_class = "sf")
sf::sf_use_s2(FALSE) 
tract_2019 <- tracts(year = 2019, cb = TRUE)     
tract_2019_wgs84 <- st_transform(tract_2019, crs = 4326)
tract_2019_boundaries <- tract_2019_wgs84 %>% select(GEOID, NAME, LSAD, geometry)
names(tract_2019_boundaries) <- c("tract_fips","tract_name","tract_namelsad","geometry")
st_write(tract_2019_boundaries,paste(geopath, "/census/tigris_tract_boundary_2010_2019.geojson", sep=""), delete_dsn = TRUE)

#----------MSA boundary----------#
options(tigris_use_cache = TRUE, tigris_class = "sf")
msa_2019 <- core_based_statistical_areas(year = 2019)
msa_2019_wgs84 <- st_transform(msa_2019, crs = 4326)
msa_2019_boundaries <- msa_2019_wgs84  %>% select(GEOID, NAME, NAMELSAD, geometry) 
names(msa_2019_boundaries)<-c("msa_fips","msa_name","msa_namelsad","geometry")
st_write(msa_2019_boundaries, paste(geopath, "/census/tigris_msa_boundary_2010_2019.geojson", sep=""), delete_dsn = TRUE)

#----------County boundary----------#
options(tigris_use_cache = TRUE, tigris_class = "sf")
county_2019 <- counties(year = 2019)
county_2019_wgs84 <- st_transform(county_2019, crs = 4326)
county_2019_boundaries <- county_2019_wgs84  %>% select(GEOID, NAME, NAMELSAD, geometry) 
names(county_2019_boundaries)<-c("county_fips","county_name","county_namelsad","geometry")
st_write(county_2019_boundaries, paste(geopath, "/census/tigris_county_boundary_2010_2019.geojson", sep=""), delete_dsn = TRUE)






