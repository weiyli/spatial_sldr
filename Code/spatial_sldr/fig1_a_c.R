# rm(list = ls())

# Model performance


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Plot the spatial map of TI,TO,TI-TO----------#
library(sf)          # poly2nb()
library(spdep)       # poly2nb() 
library(igraph)
library(patchwork)  


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


#----------rho and msa or county----------#
A.msa <- B.msa <- C.msa <- D.msa <- list()

#----------Map----------#
d=4 
d=5
Dname<<-Dname.set[d]
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif<-date.sep(Datevalue,durdate)
#----------MSA BlockID, NBlock----------#
block.msa<-sf::read_sf(paste(geopath,"/msa/",Dname,".geojson",sep=""))
BlockID<-unique(block.msa$CensusBlockGroup)
NBlock<-length(BlockID)
#----------Center point of each block----------#
block.msa$centroid <- st_centroid(block.msa$geometry,of_largest_polygon = TRUE)
s1 <- st_coordinates(block.msa$centroid)
block.msa$lon<-s1[,1]
block.msa$lat<-s1[,2]

#----------Generate a queen-type spatial weight matrix: W with lag h----------#
Wd<-poly2nb(block.msa,queen = T) 
W<<-list()
W[[1]] <- nb2mat(Wd,style = "B")  # 1 lag queen-type spatial weight matrix: w0
lag.max <- igraph::graph.adjacency(W[[1]], mode="directed", weight="W")%>%diameter()
print(lag.max)
lag.max<-2
slag<<-seq(1,lag.max,1)
for(h in 2:length(slag)){  # h lag queen-type spatial weight matrix
  if(h==2){
    W1 <- nblag(Wd, maxlag=slag[h])%>%nblag_cumul()
    W[[h]] <- nb2mat(W1,style = "B")-W[[1]]
  }else{
    W1 <- nblag(Wd, maxlag=slag[h])%>%nblag_cumul()
    W2 <- nblag(Wd, maxlag=slag[h-1])%>%nblag_cumul()
    W[[h]] <- nb2mat(W1,style = "B")-nb2mat(W2,style = "B")
  }
}
#----------Generate edges of network----------#
edges<-NULL
for(i in c(1,2)){
  edges_info <- graph_from_adjacency_matrix(W[[i]], mode = "undirected")%>%get.edgelist()
  edges<-plyr::rbind.fill(edges,data.frame(from=edges_info[,1],to=edges_info[,2],lag=i))
}
dis.edges <- distinct(edges)
colnames(dis.edges) <- c("from", "to","lag")
s1<-data.frame(from=1:NBlock,to=1:NBlock,lon=block.msa$lon,lat=block.msa$lat)
s2<-left_join(dis.edges,s1,by="from")
s3<-left_join(dis.edges,s1,by="to")
dis.edges$Olon<-s2$lon
dis.edges$Olat<-s2$lat
dis.edges$Dlon<-s3$lon
dis.edges$Dlat<-s3$lat
#----------Generate edges of network----------#
col.lag<-c("#D3D3E7","#475F78","#D4A1C5")
col.DM <- c("#FEC98DFF","#721F81FF")
for(d in 1:2){
  df<-subset(dis.edges,lag==d)
  A.msa[[d]] <- ggplot() +
    geom_segment(data = df, aes(x = Olon, y = Olat, xend = Dlon, yend = Dlat),
                 color=col.lag[d],size=0.6,
                 arrow = arrow(length = unit(0, "inches")))+
    geom_point(data = block.msa, aes(x = lon, y = lat), 
               fill = scales::alpha(col.lag[3],alpha=0.5), 
               color=scales::alpha(col.lag[3],alpha=0.8),
               size= 1.5, stroke = 1) +
    theme_wy()+
    theme(panel.border = element_blank(),
          plot.title = element_text(color = col.DM[1]), 
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "left")
}
#----------Generate queen spatial neighborhood----------#
rows <- 5
cols <- 5
grid_data <- expand.grid(x = 1:cols, y = 1:rows)
grid_data$xmid <- grid_data$x - 0.5
grid_data$ymid <- rows - grid_data$y + 0.5
queen <- matrix(c(c(2,2,2,2,2),
                  c(2,1,1,1,2),
                  c(2,1,0,1,2),
                  c(2,1,1,1,2),
                  c(2,2,2,2,2)),
                nrow = 5, ncol = 5,byrow=TRUE)
A.msa[[3]] <- ggplot(grid_data, aes(x = xmid, y = ymid)) +
  geom_tile(aes(fill = as.factor(queen[cbind(grid_data$y, grid_data$x)])), 
            color = "black", width = 1, height = 1) +
  scale_fill_manual(name ='Lag', values = c("0" = col.lag[3], "1" = col.lag[1],"2"=col.lag[2])) +
  coord_fixed(ratio = 1) +
  theme_wy()+
  theme(panel.border = element_blank(),
        plot.title = element_text(color = col.DM[1]), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right")
ggsave(A.msa[[1]],filename = paste(figpath,"/msa/lag1.pdf",sep=""),width = 4*1, height = 4*1)  
ggsave(A.msa[[2]],filename = paste(figpath,"/msa/lag2.pdf",sep=""),width = 4*1, height = 4*1)
ggsave(A.msa[[3]],filename = paste(figpath,"/msa/queen.pdf",sep=""),width = 4*1, height = 4*1)

g<- (A.msa[[1]]|A.msa[[2]]) +
  # plot_annotation(tag_levels = list('b','')) & 
  theme(plot.tag = element_text(size = 20))
ggsave(g,filename = paste(figpath,"/msa/lag.h.pdf",sep=""),width = 4*2, height = 4*1)   



#----------Spatial correlation length----------#
# exp_values <- seq(1,2, length.out = Nexp)
# Nexp <- length(exp_values)
# col.exp <- colorRampPalette(c("#984EA3","#5F4B54"))(Nexp)
exp_values <- 2
Nexp <- length(exp_values)
col.exp <- "#5F4B54"
x_values <- seq(0, 15, length.out = 100)
data<-NULL
for(i in 1:Nexp){
  y_values <- exp(-x_values / exp_values[i])
  df <- data.frame(x = x_values, y = y_values, index = exp_values[i])
  data<-plyr::rbind.fill(data,df)
}

fig.exp <- ggplot(data, aes(x = x, y = y)) +
  geom_line(aes(color = factor(index)), linetype = "solid", size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = col.exp) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(fig.exp, filename = paste(figpath,"/msa/brief_lag_exp.pdf",sep=""),width = 4*1, height = 4*1)   






















