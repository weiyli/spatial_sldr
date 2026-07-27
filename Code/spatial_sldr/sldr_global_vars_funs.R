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
# county.name.msa <- c('Fulton', 'Middlesex', 'Cook', 'Dallas', 'Harris', 'Los Angeles', 'Miami-Dade', 'New York', 'Philadelphia', 'San Francisco', 'King', 'Montgomery')
county.name.msa <- c('Fulton', 'Middlesex', 'Cook', 'Dallas', 'Harris', 'Los Angeles', 'Miami-Dade', 'Kings', 'Philadelphia', 'Alameda', 'King', 'Montgomery')
county.name.dis <- c('Harris', 'Duval', 'Harris', 'Sonoma')
county.name.set <- c(county.name.msa, county.name.dis)

# msa county fips
# county.id.msa <- c('13121', '25017', '17031', '48113', '48201', '06037', '12086', '36061', '42101', '06075', '53033', '24031')
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

# duration<-13
# startdate.msa<-rep(c('2020-02-03'),Nmsa)
# durdate.msa<-rep(c('2020-03-30'),Nmsa)
# enddate.msa<-rep(c('2020-04-12'),Nmsa) 

# 2019 Hurricane Dorian: August 28, 2019–September 2, 2019;
# 
# 2019 Tropical Storm Imelda: September 17, 2019-September 21, 2019;
# 	
# 2019 Saddleridge Wildfire: October 10, 2019-October 31, 2019;
# 
# 2019 Kincade Wildfire: October 23, 2019-November 6, 2019;
# 
# 2021 Texas Winter Freeze: February 10, 2021-February 20, 2021.

startdate.dis <- c("2017-08-25", "2019-09-02", "2021-02-10", "2019-10-23")
durdate.dis <- c("2017-08-28", "2019-09-04", "2021-02-17", "2019-10-28")
enddate.dis <- c("2017-09-01", "2019-09-05", "2021-02-20", "2019-11-06")

startdate.set <- c(startdate.msa, startdate.dis)
durdate.set <- c(durdate.msa, durdate.dis)
enddate.set <- c(enddate.msa, enddate.dis)


#-----------------Other variables--------------------#
region <- c('county','msa')
Nregion <- length(region)
region_add <- c('tract_msa','county_msa')

Yname<-c("TI_in_OO","TO_in_OO","TI_ex_OO","TO_ex_OO")
Ynum<-length(Yname)

pervalue<-c("before", "during")
pernum<-length(pervalue)
ID<-factor(c(1:2))   
period.info <- data.frame(id=ID,
                          breaks=c("before","during"),
                          labels=c("Before","During"),
                          labels.new=c("Non-disaster","Disaster"),
                          col=c("#984EA3","#5F4B54"),  
                          shape=c(21, 22)) 
                     

FTH<-0  # flow intensity threshold
DM<-c('Empirical','Model')
DM.info<-data.frame(breaks=DM,
                    labels=DM,
                    line=c("solid","dotdash"))


#----------Select MSA, Area, and Date----------#
# msa
# dindex<-which(Dname.set==Dname)
dindex<-d
Dname<-Dname.set[dindex]
Dname.id<-Dname.id.set[dindex]
event<-event.set[dindex]
state<-state.set[dindex]
county.name<-county.name.set[dindex]
county.id<-county.id.set[dindex]
# date
startdate<-startdate.set[dindex]
durdate<-durdate.set[dindex]
enddate<-enddate.set[dindex]
if(dindex>Nmsa){
  if(dindex==Nmsa+1){
    duration.dis<-as.Date(enddate)-as.Date(startdate)
  }else{
    duration.dis<-as.Date(enddate)-as.Date(startdate)+1
  }
  duration.dis<-ifelse(duration.dis<7,7,duration.dis)  # less than 7 days are calculated as 7 days
  Datevalue <- seq(as.Date(startdate)-duration.dis,as.Date(enddate)+duration.dis, by="day") 
}else{
  # for model
  Datevalue<-c(seq(as.Date(startdate),as.Date(startdate) + duration, by="day"),
               seq(as.Date(durdate),as.Date(durdate) + duration, by="day"))
  # # for data 1
  # Datevalue<-seq(as.Date(startdate),as.Date(enddate), by="day")
  
  # # for data 2
  # Datevalue <- seq(as.Date('2020-01-01'),as.Date('2020-06-30'), by="day")
  
  # Datevalue <- c(seq(as.Date('2020-02-10'),as.Date('2020-02-10')+duration, by="day"),
  #                seq(as.Date('2020-03-02'),as.Date('2020-03-02')+duration, by="day"),
  #                seq(as.Date('2020-03-30'),as.Date('2020-03-30')+duration, by="day"),
  #                seq(as.Date('2020-04-27'),as.Date('2020-04-27')+duration, by="day"))
}
Day<-length(Datevalue)
PDatevalue <- gsub('-', '/', Datevalue)  



# Datevalue<-seq(as.Date('2020-01-01'),as.Date('2020-06-30'), by="day")
# Day<-length(Datevalue)
# PDatevalue <- gsub('-', '/', Datevalue)

# Datevalue<-seq(as.Date('2020-01-01'),as.Date('2020-02-09'), by="day")
# Day<-length(Datevalue)
# PDatevalue <- gsub('-', '/', Datevalue)


# Datevalue<-seq(as.Date('2020-02-17'),as.Date('2020-03-29'), by="day")
# Day<-length(Datevalue)
# PDatevalue <- gsub('-', '/', Datevalue)



#----------Part2: Load plot variables----------#
col.all <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
             '#53738c','#a5a9b0','#a78982','#696a6c','#92699e',
             '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
             '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
             '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
             '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
             '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
             '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
             '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
             '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
             '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
             '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
shape_set <- c(1,10,13,16,19,21, 0,7,12,14,15,22)
# col.map <- c('#011638', '#2e294e', '#9055a2','#ffddd2', '#e29578', '#e07a5f')
col.lag <- c('#D3D3E7','#475F78','#D4A1C5')

# msa color for 1:Nmsa
col.msa <- col.all[1:Nmsa]
shape.msa <- c(1,10,13,16,19,21,  0,7,12,14,15,22)
msa.info<-data.frame(dname = Dname.msa, 
                     did = Dname.id.msa, 
                     breaks = "COVID-19",
                     labels = "COVID-19",
                     state = state.msa, 
                     col = col.msa,
                     shape = shape.msa)

# dis color for 1:Ndis
col.dis <- col.all[Nmsa+(1:Ndis)]
shape.dis <- c(1,10,13,16,19,21,  0,7,12,14,15,22)[1:Ndis]
dis.info <- data.frame(dname = Dname.dis, 
                       did = Dname.id.dis, 
                       rank = 1:Ndis,
                       breaks = c("Hurricane_Harvey", "Hurricane_Dorian", "Storm_Texas", "Fire_Kincade"),
                       # labels = c("Hurricane Harvey", "Hurricane Dorian", "Storm Texas", "Fire Kincade"),
                       labels = c("Hurricane Harvey", "Hurricane Dorian", "Winter Storm Uri","Fire Kincade"),
                       state = state.dis, 
                       col = col.dis,
                       shape = shape.dis)

event.info <- data.frame(dname = Dname.set, 
                         did = Dname.id.set, 
                         breaks = c(rep("COVID-19",Nmsa),
                                    "Hurricane_Harvey", "Hurricane_Dorian", "Storm_Texas", "Fire_Kincade"),
                         # labels = c(rep("COVID-19",Nmsa),
                         #            "Hurricane Harvey", "Hurricane Dorian", "Storm Texas", "Fire Kincade"),
                         labels = c(rep("COVID-19",Nmsa),
                                    "Hurricane Harvey","Hurricane Dorian","Winter Storm Uri","Fire Kincade"),
                         col = c(rep("#474747",Nmsa),c("#698CC7","#6AC7A4","#C7A469","#C7698C"))
                        )

col.fit <- c('#355070', '#6d597a', '#b56576', '#e56b6f', '#eaac8b')
col.map <- c("#D4A1C5","#475F78","#D9EDF3","#D3D3E7","white")
col.map <- c('#22223b', '#4a4e69', '#9a8c98', '#c9ada7', '#f2e9e4')

model.info <- data.frame(breaks=c("homo","power","exp"),
                         labels=c("No-decay","Power-law","Exponential"),
                         col=col.all[1:3],
                         shape=c(21,22,24))

axis.arrow <- element_line(arrow = arrow(angle = 20,length = unit(0.15, "inches"), ends = "last", type = "closed"))

#----------Part3: Load plot functions----------#
theme_wy <- function(base_size=15, base_family="Helvetica Neue") {
  (theme_foundation(base_size=base_size)
   + theme(plot.margin=unit(c(3,3,3,3),"mm"),  #t,r,b,l
           #plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
           plot.title = element_text(size = rel(1.2), hjust = 0.5),
           text = element_text(),
           # panel.background = element_rect(colour = NA),
           panel.background = element_blank(),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour=NA),
           axis.title = element_text(size = rel(1)),
           axis.title.x = element_text(vjust = 0,hjust = 0.5),
           axis.title.y = element_text(angle = 90,vjust =0.5),
           axis.text = element_text(), 
           # axis.line = element_line(colour="black",linewidth=0.5),
           axis.ticks = element_line(colour="black",linewidth=0.5),
           axis.ticks.length = unit(1.2, "mm"),
           axis.text.x = element_text(margin = unit(c(t = 1.5, r = 0, b = 0, l = 0), "mm")),
           axis.text.y = element_text(margin = unit(c(t = 0, r = 1.5, b = 0, l = 0), "mm")),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.margin = margin(t=0,unit="cm"),
           # legend.title = element_text(face="italic"),
           legend.title = element_text(face = "plain"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}






#----------Part4: Load calculation functions----------#
#----------Part4.1: od_flow.R----------#
#-----------------OD.fun: Calculate the no_visits between OD--------------------#
#' @param orig_cbgs    The original block group
#' @param dest_cbgs    The destination block group
#' @return the no_visits between OD
OD.fun<- function(orig_cbgs,dest_cbgs){
  result<-lapply(1:length(dest_cbgs), 
                 function(i){
                   s1<-strsplit(dest_cbgs[[i]],split = ",")
                   s2<-gsub("[^0-9]", "", s1[[1]])
                   s3<-data.frame(orig_bg=orig_cbgs[i],
                                  dest_bg=substr(s2, 1, 12),
                                  no_visits=substr(s2, 13, nchar(s2)))
                 }) %>% bind_rows()
  result$dest_bg <- as.numeric(result$dest_bg)
  result$no_visits <- as.numeric(result$no_visits)
  return(result)
}
#-----------------flow.fun: Calculate the Total in flow and Total out flow--------------------#
#' @param ID      Block IDs
#' @param OD      The no_visits between OD
#' @return        data.frame(ID, TI, TO)
flow.fun <- function(ID,OD){
  
  OD <- as.data.table(OD)
  inflow <- OD[, .(TI = sum(no_visits)), by = dest_bg]
  outflow <- OD[, .(TO = sum(no_visits)), by = orig_bg]
  
  df1<-data.frame(ID=ID)
  df2<-data.frame(ID=inflow$dest_bg,TI=inflow$TI)
  df3<-left_join(df1,df2,by="ID")
  df4<-data.frame(ID=outflow$orig_bg,TO=outflow$TO)
  df5<-left_join(df3,df4,by="ID")
  
  return(df5)
}

#-----------------degree.fun: Calculate the Total in degree and Total out degree--------------------#
#' @param ID        The block 
#' @param OD        The no_visits between OD
#' @return the ID,TI,TO
degree.fun<- function(ID,OD){
  
  OD <- as.data.table(OD)
  inflow <- OD[, .(TI = .N), by = dest_bg]
  outflow <- OD[, .(TO = .N), by = orig_bg]
  
  df1<-data.frame(ID=ID)
  df2<-data.frame(ID=inflow$dest_bg,TI=inflow$TI)
  df3<-left_join(df1,df2,by="ID")
  df4<-data.frame(ID=outflow$orig_bg,TO=outflow$TO)
  df5<-left_join(df3,df4,by="ID")
  
  return(df5)
}

#-----------------flow.fun.county: County-level TI / TO--------------------#
#' @param ID        County IDs (e.g. county_fips, character)
#' @param OD        data.table or data.frame with: orig_county, dest_county, no_visits
#' @return          data.frame(ID, TI, TO)
flow.fun.county <- function(ID, OD){
  
  OD <- as.data.table(OD)
  inflow <- OD[, .(TI = sum(no_visits)), by = dest_county]
  outflow <- OD[, .(TO = sum(no_visits)), by = orig_county]
  
  df1 <- data.frame(ID = ID)
  df2 <- data.frame(ID = inflow$dest_county,TI = inflow$TI)
  df3 <- left_join(df1, df2, by = "ID")
  df4 <- data.frame(ID = outflow$orig_county, TO = outflow$TO)
  df5 <- left_join(df3, df4, by = "ID")

  return(df5)
}

#-----------------degree.fun.county: County-level in/out degree--------------------#
#' @param ID        County IDs (e.g., county_fips, character)
#' @param OD        data.table or data.frame with: orig_county, dest_county, no_visits
#' @return          data.frame(ID, TI, TO)  where TI/TO are in/out degree counts
degree.fun.county <- function(ID, OD){
  
  OD <- as.data.table(OD)
  inflow <- OD[, .(TI = .N), by = dest_county]
  outflow <- OD[, .(TO = .N), by = orig_county]
  
  df1 <- data.frame(ID = ID)
  df2 <- data.frame(ID = inflow$dest_county, TI = inflow$TI)
  df3 <- dplyr::left_join(df1, df2, by = "ID")
  df4 <- data.frame(ID = outflow$orig_county,TO = outflow$TO)
  df5 <- dplyr::left_join(df3, df4, by = "ID")
  
  return(df5)
}

#-----------------flow.fun.tract: Tract-level TI / TO--------------------#
#' @param ID        County IDs (e.g. tract_fips, character)
#' @param OD        data.table or data.frame with: orig_tract, dest_tract, no_visits
#' @return          data.frame(ID, TI, TO)
flow.fun.tract <- function(ID, OD){
  OD <- data.table::as.data.table(OD)
  inflow <- OD[, .(TI = sum(no_visits)), by = dest_tract]
  outflow <- OD[, .(TO = sum(no_visits)), by = orig_tract]
  
  df1 <- data.frame(ID = as.character(ID))
  df2 <- data.frame(ID = as.character(inflow$dest_tract), TI = inflow$TI)
  df3 <- dplyr::left_join(df1, df2, by = "ID")
  df4 <- data.frame(ID = as.character(outflow$orig_tract), TO = outflow$TO)
  df5 <- dplyr::left_join(df3, df4, by = "ID")
  df5
}
#-----------------degree.fun.tract: County-level in/out degree--------------------#
#' @param ID        County IDs (e.g., tract_fips, character)
#' @param OD        data.table or data.frame with: orig_tract, dest_tract, no_visits
#' @return          data.frame(ID, TI, TO)  where TI/TO are in/out degree counts
degree.fun.tract <- function(ID, OD){
  OD <- data.table::as.data.table(OD)
  inflow <- OD[, .(TI = .N), by = dest_tract]
  outflow <- OD[, .(TO = .N), by = orig_tract]
  
  df1 <- data.frame(ID = as.character(ID))
  df2 <- data.frame(ID = as.character(inflow$dest_tract), TI = inflow$TI)
  df3 <- dplyr::left_join(df1, df2, by = "ID")
  df4 <- data.frame(ID = as.character(outflow$orig_tract), TO = outflow$TO)
  df5 <- dplyr::left_join(df3, df4, by = "ID")
  df5
}




#-----------------date.sep: week and period for COVID-19--------------------#
#' @param Datevalue  Period
#' @param durdate    The start of disaster
#' @return data for week and period
date.sep<- function(Datevalue,durdate){
  s1<-format(Datevalue,'%w')
  s2<-which(s1==0|s1==6)
  s3<-1:Day
  s3[s2]<-"weekend"
  s3[-s2]<-"weekday"
  # before disaster, or during disaster: s6
  s4 <- which(Datevalue==durdate)
  s5 <- 1:Day
  if(length(s4)>0){
    s5[1:(s4-1)] <- "before"
    s5[s4:Day] <- "during"
  }else{
    s5[1:Day] <- "before"
  }
  data<-data.frame(day=Datevalue,
                   week=s3,
                   period=s5)
  return(data)
}


#-----------------data.sep.dis: week and period for disaster--------------------#
#' @param Datevalue  Period
#' @param startdate  The start of disaster
#' @param enddate    The end of disaster
#' return data for week and period
date.sep.dis <- function(Datevalue,startdate,enddate){
  s1<-format(Datevalue,'%w')
  s2<-which(s1==0|s1==6)
  s3<-1:Day
  s3[s2]<-"weekend"
  s3[-s2]<-"weekday"
  #before flood, during flood, or after flood: s6
  s4<-which(Datevalue==startdate)
  s5<-which(Datevalue==enddate)
  s6<-1:Day
  s6[1:(s4-1)]<-"before"
  s6[s4:s5]<-"during"
  s6[(s5+1):Day]<-"after"
  data<-data.frame(day=Datevalue,
                   week=s3,
                   period=s6)
  return(data)
}


#----------Generates a weights matrix for a neighbours list: nb2mat----------#
#' @param nb  1 lag queen neighborhood list
#' return matrix of 1 lag queen-type adjacency matrix--W1
nb.mat <- function(nb) {
  
  n <- length(nb)
  mat <- matrix(0, n, n)
  
  # Get the number of neighbors for each node
  cards <- card(nb)
  
  # Change elements at positions with zero neighbors with 1
  cards0<-which(cards==0)
  if(length(cards0)>0){
    cards[cards0] <- 1
    nb[[cards0]] <- 1
  }
  
  # Create the adjacency matrix using matrix indexing
  rows <- rep(1:n, cards)
  cols <- unlist(nb)
  mat[cbind(rows, cols)] <- 1
  
  # Change elements at positions with zero neighbors with 0
  if(length(cards0)>0){
    mat[cards0,] <- 0
  }
  
  return(mat)
}


#----------Compute higher-order queen adjacency matrices: queen.adj----------#
#' @param W1        1 lag queen matrix
#' @param lag.max   the maximum lag
#' return list of 1-n lag queen-type adjacency matrix--W
queen.adj <- function(W1, lag.max) {
  # Initialize a list to store adjacency matrices of different orders
  W<-list()
  W[[1]] <- W1
  
  # Track already marked adjacency relationships
  already <- (W1 == 1)
  
  for (h in 2:lag.max) {
    # Compute the current order adjacency matrix
    Wh <- ((W[[h-1]] %*% W1) > 0) * 1
    diag(Wh) <- 0
    
    # Avoid duplicate adjacency relationships
    Wh[already & (Wh == 1)] <- 0
    W[[h]] <- Wh
    
    # Update already marked adjacency relationships
    already <- already | (W[[h]] == 1)
  }
  
  return(W)
}

#----------Error-measure: R2----------#
R2 <- function(y_test,y_predict){
  SS.tot = sum((y_test-mean(y_test))^2)
  SS.err = sum((y_predict-y_test)^2)
  result <- 1-SS.err/SS.tot
  return(result)
}


#----------SSI (Sørensen similarity index)----------#
calc_ssi <- function(y_test, y_predict, eps = 1e-12) {
  stopifnot(length(y_test) == length(y_predict))
  den <- sum(y_test, na.rm = TRUE) + sum(y_predict, na.rm = TRUE)
  if (den <= eps) return(NA_real_)
  2 * sum(pmin(y_test, y_predict), na.rm = TRUE) / den
}
#----------CPC (asymmetric, w.r.t observed flows)----------#
calc_cpc <- function(y_test, y_predict, eps = 1e-12) {
  stopifnot(length(y_test) == length(y_predict))
  den <- sum(y_test, na.rm = TRUE)
  if (den <= eps) return(NA_real_)
  sum(pmin(y_test, y_predict), na.rm = TRUE) / den
}
#----------PCC (Pearson correlation coefficient)----------#
calc_pcc <- function(y_test, y_predict) {
  stopifnot(length(y_test) == length(y_predict))
  if (sd(y_test, na.rm = TRUE) == 0 || sd(y_predict, na.rm = TRUE) == 0)
    return(NA_real_)
  cor(y_test, y_predict, use = "complete.obs", method = "pearson")
}
#----------SRC (Spearman rank correlation)----------#
calc_src <- function(y_test, y_predict) {
  stopifnot(length(y_test) == length(y_predict))
  if (sd(y_test, na.rm = TRUE) == 0 || sd(y_predict, na.rm = TRUE) == 0)
    return(NA_real_)
  cor(y_test, y_predict, use = "complete.obs", method = "spearman")
}


#----------Error-measure: error_test----------#
#' #'Root-mean-squared-error (RMSE) between observed dynamics and model predictions
# error_test <- function(y_test,y_predict){
# index0 <- which(y_test==0)
# if(length(index0)==0){
#   rmse_test = (sum((y_predict-y_test)^2)/length(y_test))^0.5
#   rmspe_test = (sum(((y_predict-y_test)/y_test)^2*100)/length(y_test))^0.5
#   result <- (rmse_test+rmspe_test)/2
# }else{
#   #y_test==0
#   rmse_test0 <- (sum((y_predict[index0]-y_test[index0])^2)/length(y_test))^0.5
#   #y_test!=0
#   rmse_test = (sum((y_predict[-index0]-y_test[-index0])^2)/length(y_test))^0.5
#   rmspe_test = (sum(((y_predict[-index0]-y_test[-index0])/y_test[-index0])^2*100)/length(y_test))^0.5
#   result <- (rmse_test+rmspe_test)/2+rmse_test0
# }
# return(result)
# }

error_rmse_rmspe <- function(y_test,y_predict){
  index0 <- which(y_test==0)
  if(length(index0)==0){
    rmse_test = (sum((y_predict-y_test)^2)/length(y_test))^0.5
    rmspe_test = (sum(((y_predict-y_test)/y_test)^2*100)/length(y_test))^0.5
    result <- (rmse_test+rmspe_test)/2
  }else{
    #y_test==0
    rmse_test0 <- (sum((y_predict[index0]-y_test[index0])^2)/length(y_test))^0.5
    #y_test!=0
    rmse_test = (sum((y_predict[-index0]-y_test[-index0])^2)/length(y_test))^0.5
    rmspe_test = (sum(((y_predict[-index0]-y_test[-index0])/y_test[-index0])^2*100)/length(y_test))^0.5
    result <- (rmse_test+rmspe_test)/2+rmse_test0
  }
  return(result)
}
error_rmse <- function(y_test, y_predict){
  y_predict <- pmax(y_predict, 0)  # non-negativity
  sqrt(mean((y_predict - y_test)^2))
}
error_rmspe <- function(y_test, y_predict, c = 1){
  y_predict <- pmax(y_predict, 0)
  sqrt(mean(((y_predict - y_test) / (y_test + c))^2))
}
error_logrmse <- function(y_test, y_predict){
  y_predict <- pmax(y_predict, 0) 
  sqrt(mean((log1p(y_predict) - log1p(y_test))^2))
}

#----------error metric configs (choose what you want to run)----------#
error_specs <- list(
  rmse_rmspe = list(
    tag = "rmse_rmspe",
    fun = error_rmse_rmspe
  ),
  rmse = list(
    tag = "rmse",
    fun = error_rmse
  ),
  rmspe = list(
    tag = "srmspe",
    fun = function(y_test, y_predict) error_rmspe(y_test, y_predict, c = 1)
  ),
  logrmse = list(
    tag = "logrmse",
    fun = error_logrmse
  )
)


#----------constant decay: rho----------#
#----------Local Parameter Estimation: homo.param.est----------#
#' @param rho     spatial correlation length parameters 
#' return optimal local parameters--(rho)
homo.param.est <- function (rho) {
  lag.h <- matrix(rho, nrow = NBlock, ncol = length(slag), byrow = TRUE)
  WM.h <- sapply(W, "%*%", M.day)
  M.day.sim<-rowSums(lag.h*WM.h)
  return(error_test(M.day,M.day.sim))
}
#----------Return the optimal phi, rho, error, fitting results and convergence----------#
homo.day.est <- function (g) {
  M.day<<-Mob[,g]
  fit <- optimize(f = homo.param.est, interval = c(lower_bound, upper_bound))
  return(c(fit$minimum, fit$objective))
}

#----------power law decay: rho----------#
#----------Local Parameter Estimation: power.param.est----------#
#' @param rho     spatial correlation length parameters 
#' return optimal local parameters--(rho)
power.param.est <- function (rho) {
  lag.h <- matrix(rep(slag^(-rho), NBlock), ncol = length(slag), byrow = TRUE)
  WM.h <- sapply(W, "%*%", M.day)
  M.day.sim<-rowSums(lag.h*WM.h)
  return(error_test(M.day,M.day.sim))
}
#----------Return the optimal phi, rho, error, fitting results and convergence----------#
power.day.est <- function (g) {
  M.day<<-Mob[,g]
  fit <- optimize(f = power.param.est, interval = c(lower_power, upper_power))
  return(c(fit$minimum, fit$objective))
}

#----------exp decay: rho----------#
#----------Local Parameter Estimation: exp.param.est----------#
#' @param rho     spatial correlation length parameters 
#' return optimal local parameters--(rho)
exp.param.est <- function (rho) {
  lag.h <- matrix(rep(exp(-slag/rho), NBlock), ncol = length(slag), byrow = TRUE)
  WM.h <- sapply(W, "%*%", M.day)
  M.day.sim<-rowSums(lag.h*WM.h)
  return(error_test(M.day,M.day.sim))
}
#----------Return the optimal phi, rho, error, fitting results and convergence----------#
exp.day.est <- function (g) {
  M.day<<-Mob[,g]
  fit <- optimize(f = exp.param.est, interval = c(lower_bound, upper_bound))
  return(c(fit$minimum, fit$objective))
}



#----------power law decay: con and rho----------#
#----------Local Parameter Estimation: npower.param.est----------#
#' @param params  con is lag-specific normalisation constant; rho is spatial correlation length parameters 
#' @references    Li, Y., Schubert, S., Kropp, J.P. et al. On the influence of density and morphology on the Urban Heat Island intensity. Nat Commun 11, 2647 (2020). https://doi.org/10.1038/s41467-020-16461-9
#' return optimal local parameters--(rho)
npower.param.est <- function (params) {
  con <- params[1]
  rho <- params[2]
  norms <- con*sum(slag^(-1/rho))
  lag.h <- matrix(rep(slag^(-1/rho)/norms, NBlock), ncol = length(slag), byrow = TRUE)
  WM.h <- sapply(W, "%*%", M.day)
  M.day.sim <- rowSums(lag.h*WM.h)
  return(error_test(M.day,M.day.sim))
}
#----------Return the optimal phi, rho, error, fitting results and convergence----------#
npower.day.est <- function (g) {
  M.day <<- Mob[,g]
  fit <- optim(par = init.params, fn = npower.param.est,
               method = Method, control = list(maxit=5000))
  return(c(fit$par, fit$value, fit$convergence))
}

#----------exp decay: con and rho----------#
#----------Local Parameter Estimation: nexp.param.est----------#
#' @param params  con is lag-specific normalisation constant; rho is spatial correlation length parameters 
#' return optimal local parameters--(con,rho)
nexp.param.est <- function (params) {
  con <- params[1]
  rho <- params[2]
  norms <- con*sum(exp(-slag/rho))
  lag.h <- matrix(rep(exp(-slag/rho)/norms, NBlock), ncol = length(slag), byrow = TRUE)
  WM.h <- sapply(W, "%*%", M.day)
  M.day.sim<-rowSums(lag.h*WM.h)
  return(error_test(M.day,M.day.sim))
}
#----------Return the optimal phi, rho, error, fitting results and convergence----------#
nexp.day.est <- function (g) {
  M.day<<-Mob[,g]
  fit <- optim(par = init.params, fn = nexp.param.est,
               method = Method, control = list(maxit=5000))
  return(c(fit$par,fit$value,fit$convergence))
}


