# rm(list = ls())

# SLDR fitting and higher-order queen weight matrix for Japan 500 m grids


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
datapath <- 'D:/ood/Data/spatial_sldr/japan'


#----------Load packages----------#
d<<-1
source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))


#----------Load cities----------#
Dname.japan <- c('cityA',
                 'cityB',
                 'cityC',
                 'cityD')
Dfile.japan <- c('cityA_groundtruthdata.csv.gz',
                 'cityB_challengedata.csv.gz',
                 'cityC_challengedata.csv.gz',
                 'cityD_challengedata.csv.gz')
Dname.set <- Dname.japan
Ncity <- length(Dname.set)

yindex <- 2
lag.max.user <- 20     # NA: use the full grid extent
save.od <- TRUE
save.fit <- TRUE

#----------Set the global function pointer used by parameter estimation----------#
error_test <<- error_specs[["rmse_rmspe"]]$fun
lower_bound <<- 0
upper_bound <<- 5
lower_power <<- 40
upper_power <<- 80


#----------Output paths----------#
dir.create(paste(datapath,"/processed",sep=""), showWarnings=FALSE, recursive=TRUE)
dir.create(paste(datapath,"/param",sep=""), showWarnings=FALSE, recursive=TRUE)
dir.create(paste(datapath,"/fitting",sep=""), showWarnings=FALSE, recursive=TRUE)


#----------Read csv.gz without R.utils----------#
read.csv.gz <- function(file){
  con <- gzfile(file,open="rt")
  on.exit(close(con))
  return(read.csv(con,header=TRUE)%>%as.data.table)
}


#----------Write csv.gz----------#
write.csv.gz <- function(x,file){
  con <- gzfile(file,open="wt")
  on.exit(close(con))
  write.csv(x,con,row.names=FALSE)
}


#----------Square-window sum for regular grids----------#
window.sum <- function(M,radius){
  if(radius==0) return(M)
  
  nr <- nrow(M)
  nc <- ncol(M)
  M.pad <- matrix(0,nrow=nr+2*radius,ncol=nc+2*radius)
  M.pad[(radius+1):(radius+nr),(radius+1):(radius+nc)] <- M
  
  M.sum <- rbind(0,cbind(0,apply(apply(M.pad,2,cumsum),1,cumsum)))
  row.id <- 1:nr
  col.id <- 1:nc
  
  return(M.sum[row.id+2*radius+1,col.id+2*radius+1,drop=FALSE]-
           M.sum[row.id,col.id+2*radius+1,drop=FALSE]-
           M.sum[row.id+2*radius+1,col.id,drop=FALSE]+
           M.sum[row.id,col.id,drop=FALSE])
}


#----------Queen-ring mobility sums----------#
queen.ring.sum <- function(M,lag.max){
  WM <- matrix(0,nrow=length(M),ncol=lag.max)
  M.previous <- M
  for(h in 1:lag.max){
    M.current <- window.sum(M,h)
    WM[,h] <- as.vector(M.current-M.previous)
    M.previous <- M.current
  }
  return(WM)
}


#----------Trajectory to OD flows----------#
trajectory.od <- function(trajectory){
  setorderv(trajectory,c("uid","d","t"))
  trajectory[,`:=`(
    dest.x=shift(x,type="lead"),
    dest.y=shift(y,type="lead")
  ),by=.(uid,d)]
  
  od <- trajectory[
    !is.na(dest.x)&!is.na(dest.y),
    .(day=d,
      orig.x=x,
      orig.y=y,
      dest.x=as.integer(dest.x),
      dest.y=as.integer(dest.y))
  ]
  od <- od[,.(no.visits=.N),by=.(day,orig.x,orig.y,dest.x,dest.y)]
  return(od)
}


#----------Mobility matrix----------#
mobility.matrix <- function(od,grid,day.value,Yname.value){
  if(Yname.value=='TO_in_OO'){
    mob <- od[day==day.value,
              .(mob=sum(no.visits)),
              by=.(x=orig.x,y=orig.y)]
  }
  if(Yname.value=='TI_in_OO'){
    mob <- od[day==day.value,
              .(mob=sum(no.visits)),
              by=.(x=dest.x,y=dest.y)]
  }
  
  mob <- left_join(grid,mob,by=c("x","y"))%>%as.data.table
  mob[is.na(mob),mob:=0]
  setorder(mob,y,x)
  return(matrix(mob$mob,
                nrow=length(unique(grid$x)),
                ncol=length(unique(grid$y))))
}


#----------SLDR fitting for one day----------#
sldr.day.est <- function(M.day,WM,slag){
  homo.predict <- function(rho){
    return(as.numeric(rho*rowSums(WM)))
  }
  power.predict <- function(rho){
    return(as.numeric(WM%*%(slag^(-rho))))
  }
  exp.predict <- function(rho){
    return(as.numeric(WM%*%exp(-slag/rho)))
  }
  
  homo.fit <- optimize(function(rho){
    error_test(M.day,homo.predict(rho))
  },interval=c(lower_bound,upper_bound))
  power.fit <- optimize(function(rho){
    error_test(M.day,power.predict(rho))
  },interval=c(lower_power,upper_power))
  exp.fit <- optimize(function(rho){
    error_test(M.day,exp.predict(rho))
  },interval=c(max(0.01,lower_bound),upper_bound))
  
  params <- data.frame(
    rho=c(homo.fit$minimum,power.fit$minimum,exp.fit$minimum),
    error=c(homo.fit$objective,power.fit$objective,exp.fit$objective),
    index=c('homo','power','exp')
  )
  fitting <- data.frame(
    empirical=M.day,
    homo=homo.predict(homo.fit$minimum),
    power=power.predict(power.fit$minimum),
    exp=exp.predict(exp.fit$minimum)
  )
  return(list(params=params,fitting=fitting))
}


#----------SLDR fitting for Japan----------#
dis.rho <- dis.r2 <- NULL
for(d in 1:Ncity){
  
  Dname <- Dname.set[d]
  city.file <- paste(datapath,"/",Dfile.japan[d],sep="")
  message("Reading ",Dname,": ",city.file)
  
  #----------Trajectory data----------#
  trajectory <- read.csv.gz(city.file)
  trajectory <- trajectory[
    !is.na(uid)&!is.na(d)&!is.na(t)&!is.na(x)&!is.na(y)&
      x!=999&y!=999
  ]
  trajectory[,`:=`(
    uid=as.character(uid),
    d=as.integer(d),
    t=as.integer(t),
    x=as.integer(x),
    y=as.integer(y)
  )]
  
  #----------Grid cells----------#
  grid <- CJ(
    x=seq(min(trajectory$x),max(trajectory$x),1),
    y=seq(min(trajectory$y),max(trajectory$y),1)
  )
  setorder(grid,y,x)
  grid[,cell.index:=.I]
  grid[,cell.id:=paste(x,y,sep="_")]
  write.csv(grid,
            file=paste(datapath,"/processed/",Dname,"_grid_cells.csv",sep=""),
            row.names=FALSE)
  
  NBlock <- nrow(grid)
  nx <- length(unique(grid$x))
  ny <- length(unique(grid$y))
  lag.max <- ifelse(is.na(lag.max.user),
                    max(nx-1,ny-1),
                    min(lag.max.user,max(nx-1,ny-1)))
  slag <- seq(1,lag.max,1)
  
  #----------OD flows----------#
  message("Constructing OD for ",Dname)
  od <- trajectory.od(trajectory)
  rm(trajectory)
  gc()
  if(save.od){
    write.csv.gz(
      od,
      paste(datapath,"/processed/",Dname,"_OD.csv.gz",sep="")
    )
  }
  Datevalue <- sort(unique(od$day))
  Datevalue <- tail(Datevalue, 30)
  Day <- length(Datevalue)
  
    #----------1.parameters----------#
    SLDR <- SLDR.fit <- NULL
    for(i in 1:Day){
      message("Fitting ",Dname," ",Yname[yindex],", day ",Datevalue[i])
      
      M.matrix <- mobility.matrix(od,grid,Datevalue[i],Yname[yindex])
      M.day <- as.vector(M.matrix)
      WM <- queen.ring.sum(M.matrix,lag.max)
      result <- sldr.day.est(M.day,WM,slag)
      
      params.day <- result$params
      params.day$day <- Datevalue[i]
      params.day$lag <- lag.max
      SLDR <- plyr::rbind.fill(SLDR,params.day)
      
      #----------2.fitting----------#
      fitting.day <- result$fitting
      fitting.day$day <- Datevalue[i]
      fitting.day$cell.index <- grid$cell.index
      fitting.day$x <- grid$x
      fitting.day$y <- grid$y
      SLDR.fit <- plyr::rbind.fill(SLDR.fit,fitting.day)
      
      rm(M.matrix,M.day,WM,result)
      gc()
    } # day
    
    write.csv(
      SLDR,
      file=paste(datapath,"/param/",Dname,
                 "_SLDR_params_",Yname[yindex],".csv",sep=""),
      row.names=FALSE
    )
    if(save.fit){
      write.csv.gz(
        SLDR.fit,
        paste(datapath,"/fitting/",Dname,
              "_SLDR_fit_",Yname[yindex],".csv.gz",sep="")
      )
    }
    
    #----------R2 and rank correlations----------#
    SLDR.fit <- as.data.table(SLDR.fit)
    r2 <- SLDR.fit[,.(r2_homo=R2(empirical,homo),
                      r2_power=R2(empirical,power),
                      r2_exp=R2(empirical,exp)),
                   by=day]
    rcor <- SLDR.fit[,.(pearson_homo=calc_pcc(empirical,homo),
                        pearson_power=calc_pcc(empirical,power),
                        pearson_exp=calc_pcc(empirical,exp),
                        spearman_homo=calc_src(empirical,homo),
                        spearman_power=calc_src(empirical,power),
                        spearman_exp=calc_src(empirical,exp),
                        kendall_homo=cor(empirical,homo,method="kendall"),
                        kendall_power=cor(empirical,power,method="kendall"),
                        kendall_exp=cor(empirical,exp,method="kendall")),
                     by=day]
    r2$city <- rcor$city <- Dname
    r2$yname <- rcor$yname <- Yname[yindex]
    
    SLDR$city <- Dname
    SLDR$yname <- Yname[yindex]
    dis.rho <- plyr::rbind.fill(dis.rho,SLDR)
    dis.r2 <- plyr::rbind.fill(dis.r2,left_join(r2,rcor,by=c("day","city","yname")))
    
  rm(od,grid,SLDR.fit)
  gc()
  print(d)
} # city


#----------SLDR fitting results for all Japan cities----------#
write.csv(
  dis.rho,
  file=paste(datapath,"/SLDR_params_",Yname[yindex],".csv",sep=""),
  row.names=FALSE
)
write.csv(
  dis.r2,
  file=paste(datapath,"/SLDR_r2_",Yname[yindex],".csv",sep=""),
  row.names=FALSE
)
