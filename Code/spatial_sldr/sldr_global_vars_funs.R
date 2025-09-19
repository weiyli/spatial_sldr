# rm(list = ls())


#----------Load data packages----------#
library(data.table) # setDT
library(plyr)       # rbind.fill()
library(dplyr)      # left_join()
#----------Load plot packages----------#
library(RColorBrewer)  #set color 
library(scales)     # date_format()  viridis
library(patchwork)  # combine figs
library(grid)
library(ggplot2)
library(latex2exp)  # TeX()
library(ggthemes)


#----------Part1: Data process----------#
#----------Disaster, Area, and Date----------#
# msa name
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
# number of msa and dis
Nmsa<-length(Dname.msa)
Ndis<-length(Dname.dis)
# disaster name
event.msa <- rep(c('COVID-19'),Nmsa)
event.dis <- c("Hurricane_Harvey", "Hurricane_Dorian", "Storm_Texas", "Fire_Kincade")
# event.dis <- c("2017 Hurricane Harvey", "2019 Hurricane Dorian", "2021 Storm Texas", "2019 Fire Kincade")
event.set <- c(event.msa,event.dis)
# msa code
Dname.id.msa<-c(12060,14460,16980,19100,26420,31080,33100,35620,37980,41860,42660,47900)
Dname.id.dis <- c(26420, 27260, 26420, 42220)
Dname.id.set<-c(Dname.id.msa,Dname.id.dis)
# msa state
state.msa <- c('GA',
               'MA', 
               'IL', 
               'TX', 
               'TX', 
               'CA', 
               'FL', 
               'NY', 
               'PA', 
               'CA', 
               'WA', 
               'DC')
state.dis <- c('TX', 
               'FL', 
               'TX', 
               'CA')
state.set <- c(state.msa, state.dis)

# msa county name
county.name.msa <- c('Fulton', 'Middlesex', 'Cook', 'Dallas', 'Harris', 'Los Angeles', 'Miami-Dade', 'Kings', 'Philadelphia', 'Alameda', 'King', 'Montgomery')
county.name.dis <- c('Harris', 'Duval', 'Harris', 'Sonoma')
county.name.set <- c(county.name.msa, county.name.dis)

# msa county fips
county.id.msa <- c('13121', '25017', '17031', '48113', '48201', '06037', '12086', '36047', '42101', '06001', '53033', '24031')
county.id.dis <- c('48201', '12031', '48201', '06097')
county.id.set <- c(county.id.msa, county.id.dis)


last.day <- c("2020-01-31", "2020-02-29", "2020-03-31", "2020-04-30", 
              "2020-05-31", "2020-06-30","2020-07-31", "2020-08-31", 
              "2020-09-30", "2020-10-31", "2020-11-30", "2020-12-31")%>%as.Date

# msa date
duration<-6
startdate.msa<-rep(c('2020-02-10'),Nmsa)
durdate.msa<-rep(c('2020-03-30'),Nmsa)
enddate.msa<-rep(c('2020-04-05'),Nmsa)
