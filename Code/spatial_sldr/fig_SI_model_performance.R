# rm(list = ls())

# Data from sldr_fit_msa.R and sldr_fit_msa_error.R
# Model selection


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
#----------Part1: COVID-19----------#
datapath <- 'D:/ood/Data/spatial_sldr'
figpath <- 'D:/ood/Figure/spatial_sldr'


#----------Load packages----------#
library(ggsignif)   # gemo_signif()
library(ggtext)     # element_markdown()


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
# Dname<<-Dname.set[1]
# source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
date.gif<-date.sep(Datevalue,durdate)
#----------set the global function pointer used by homo/power/exp.param.est----------#
error_test <<- error_specs[["rmse_rmspe"]]$fun

#----------SMAPE----------#
for(yindex in 1:Ynum){
  
  SMAPE<-NULL
  for(d in 1:Nmsa){
    
    Dname<<-Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    date.gif<-date.sep(Datevalue,durdate)
    
    #----------Data from sldr_fit.R: SLDR_fit----------# 
    sldr.fit <- read.csv(paste(datapath,"/msa/fit/",Dname,"_SLDR_fit_",Yname[yindex],".csv",sep=""), header=TRUE)%>%setDT
    sldr.fit <- melt(sldr.fit, id.vars = c("day", "empirical"), 
                     measure.vars = c("homo", "exp", "power"),
                     variable.name = "model", value.name = "predicted")
    error.msa <- sldr.fit[, .(error = error_test(empirical, predicted)), by = .(model,day)]
    error.msa$msa <- Dname
    SMAPE <- plyr::rbind.fill(SMAPE, error.msa)
  } # msa
  write.csv(SMAPE, file=paste(datapath,"/msa/SMAPE_",Yname[yindex],".csv",sep=""), row.names = FALSE)
} # yindex


#----------SSI, CPC, PCC, SRC----------#
for (yindex in 1:Ynum) {
  
  SSI_out   <- NULL
  CPC_out   <- NULL
  PCC_out   <- NULL
  SRC_out   <- NULL
  SMAPE_out <- NULL
  
  for (d in 1:Nmsa) {
    
    Dname <<- Dname.set[d]
    source(paste(codepath,"/sldr_global_vars_funs.R",sep=""))
    date.gif <- date.sep(Datevalue, durdate)   # should include: day, period
    
    #----------Data from sldr_fit.R: SLDR_fit----------#
    sldr.fit <- fread(paste(datapath, "/msa/fit/", Dname, "_SLDR_fit_", Yname[yindex], ".csv", sep = ""))
    
    sldr.fit <- melt(sldr.fit,
                     id.vars = c("day", "empirical"),
                     measure.vars = c("homo", "exp", "power"),
                     variable.name = "model", value.name = "predicted")
    
    sldr.fit <- merge(sldr.fit, as.data.table(date.gif), by = "day", all.x = TRUE)
    
    #----------SMAPE----------#
    smape.msa <- sldr.fit[, .(value = error_test(empirical, predicted)), by = .(model, period)]
    smape.msa$msa <- Dname
    SMAPE_out <- plyr::rbind.fill(SMAPE_out, smape.msa)
    
    #----------SSI----------#
    ssi.msa <- sldr.fit[, .(value = calc_ssi(empirical, predicted)), by = .(model, period)]
    ssi.msa$msa <- Dname
    SSI_out <- plyr::rbind.fill(SSI_out, ssi.msa)
    
    #----------CPC----------#
    cpc.msa <- sldr.fit[, .(value = calc_cpc(empirical, predicted)), by = .(model, period)]
    cpc.msa$msa <- Dname
    CPC_out <- plyr::rbind.fill(CPC_out, cpc.msa)
    
    #----------PCC----------#
    pcc.msa <- sldr.fit[, .(value = calc_pcc(empirical, predicted)), by = .(model, period)]
    pcc.msa$msa <- Dname
    PCC_out <- plyr::rbind.fill(PCC_out, pcc.msa)
    
    #----------SRC (Spearman rank correlation)----------#
    src.msa <- sldr.fit[, .(value = calc_src(empirical, predicted)), by = .(model, period)]
    src.msa$msa <- Dname
    SRC_out <- plyr::rbind.fill(SRC_out, src.msa)
    
  } # msa
  
  #----------Write----------#
  SMAPE_out$metric <- "SMAPE"
  SSI_out$metric   <- "SSI"
  CPC_out$metric   <- "CPC"
  PCC_out$metric   <- "PCC"
  SRC_out$metric   <- "SRC"
  
  metrics <- plyr::rbind.fill(SMAPE_out, SSI_out, CPC_out, PCC_out, SRC_out)
  
  write.csv(metrics, file = paste(datapath, "/msa/metrics_", Yname[yindex], ".csv", sep = ""), row.names = FALSE)
  
} # yindex


#----------Plot----------# 
A.msa <- B.msa <- C.msa <- D.msa <- list()
Sys.setlocale("LC_TIME", "C")
id.week <- paste(format(Datevalue[1:(Day/2)],"%m-%d"),format(Datevalue[1:(Day/2)], "%a"),sep=" ")
id.week.during <- paste(format(Datevalue[Day/2+(1:(Day/2))],"%m-%d"),format(Datevalue[Day/2+(1:(Day/2))], "%a"),sep=" ")
for(yindex in 2:2){
  
  SMAPE <- read.csv(paste(datapath,"/msa/SMAPE_",Yname[yindex],".csv",sep=""), header=TRUE)
  SMAPE$day <- as.Date(SMAPE$day)
  SMAPE <- left_join(SMAPE,date.gif,by="day")
  SMAPE.before <- subset(SMAPE,period==pervalue[1])
  SMAPE.before <- left_join(SMAPE.before,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')
  
  SMAPE.during <- subset(SMAPE,period==pervalue[2])
  SMAPE.during <- left_join(SMAPE.during,data.frame(day=Datevalue,day.id=c(1:Day)),by='day')
  
  metrics <- fread(paste0(datapath, "/msa/metrics_", Yname[yindex], ".csv"))
  metrics.before <- metrics[period == pervalue[1] & metric!="SMAPE"]
  metrics.during <- metrics[period == pervalue[2] & metric!="SMAPE"]
  
  #----------Before----------# 
  # daily SMAPE
  A.msa[[yindex]] <- ggplot(SMAPE.before, aes(x = day.id, y = error)) +
    geom_line(aes(color = model), linewidth = 1) +
    geom_point(aes(color = model, fill = model, shape=model), size=3, stroke = 1) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    scale_x_continuous(breaks = c(1:(Day/2)),labels = id.week)+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = "Date", y = "RMSE-RMSPE", title = NULL) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain",size=15),
          axis.text.y = element_text(hjust = 0),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top",
          legend.key.size = unit(1, "cm")) 
  # all metrics: CPC, PCC, SSI, SRC
  C.msa[[yindex]] <- ggplot(metrics.before, aes(x = msa, y = value)) +
    geom_line(aes(color = model, group = model), linewidth = 1) +
    geom_point(aes(fill = model, color = model, shape = model),
               size = 3, stroke = 1) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    guides(fill = guide_legend(override.aes = list(size=3)))+
    labs(x = "MSA", y = "Metric value") +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
          strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
          strip.text = element_text(face = "plain",size=15),
          legend.position = "top",
          legend.key.size = unit(1, "cm"),
          axis.title.x = element_text(margin = margin(b = 4)),
          axis.text.x = ggtext::element_markdown(size = rel(1), angle = 40, hjust = 1))
  
  fig_model_before <- (A.msa[[yindex]]/C.msa[[yindex]]) + plot_layout(heights = c(3, 1)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig_model_before, filename = paste(figpath,"/msa/SI_fit_model_before_",Yname[yindex],".pdf",sep=""), width = 4*4, height = 4.3*4)
 
  
  #----------During----------# 
  # daily SMAPE
  B.msa[[yindex]] <- ggplot(SMAPE.during, aes(x = day.id, y = error)) +
    geom_line(aes(color = model), linewidth = 1) +
    geom_point(aes(color = model, fill = model, shape=model), size=3, stroke = 1) +
    facet_wrap(~ msa, ncol = 4, scales = "fixed") + 
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    scale_x_continuous(breaks = c(Day/2+(1:(Day/2))),labels = id.week.during)+
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x)))+
    labs(x = "Date", y = "RMSE-RMSPE", title = NULL) +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
          panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
          panel.spacing = unit(0.1, "lines"),
          axis.line.x.bottom = element_line(color = "black", linewidth = 0.5), 
          strip.text = element_text(face = "plain",size=15),
          axis.text.y = element_text(hjust = 0),
          axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top",
          legend.key.size = unit(1, "cm"))
  # all metrics: CPC, PCC, SSI, SRC
  D.msa[[yindex]] <- ggplot(metrics.during, aes(x = msa, y = value)) +
    geom_line(aes(color = model, group = model), linewidth = 1) +
    geom_point(aes(fill = model, color = model, shape = model),
               size = 3, stroke = 1) +
    facet_wrap(~ metric, nrow = 1, scales = "free_y") +
    scale_fill_manual(name = "Model",
                      breaks = model.info$breaks,
                      labels = model.info$labels,
                      values = scales::alpha(model.info$col,alpha=1))+ 
    scale_color_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$col) + 
    scale_shape_manual(name = "Model",
                       breaks = model.info$breaks,
                       labels = model.info$labels,
                       values = model.info$shape) + 
    guides(fill = guide_legend(override.aes = list(size=3)))+
    labs(x = "MSA", y = "Metric value") +
    theme_wy() +
    theme(panel.background = element_blank(),
          panel.spacing = unit(0.2, "lines"),
          panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
          strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
          strip.text = element_text(face = "plain",size=15),
          legend.position = "top",
          legend.key.size = unit(1, "cm"),
          axis.title.x = element_text(margin = margin(b = 4)),
          axis.text.x = ggtext::element_markdown(size = rel(1), angle = 40, hjust = 1))
  fig_model_during <- (B.msa[[yindex]]/D.msa[[yindex]]) + plot_layout(heights = c(3, 1)) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
  ggsave(fig_model_during, filename = paste(figpath,"/msa/SI_fit_model_during_",Yname[yindex],".pdf",sep=""), width = 4*4, height = 4.3*4)
  
} # Yname










#----------Robustness across loss functions: exponential kernel win rate----------#
param_all <- NULL
for (s in 2:Nregion) {
  for (yindex in 2:2) {    
    for (d in 1:Nmsa) {
      
      Dname <- Dname.msa[d]
      source(paste(codepath, "/sldr_global_vars_funs.R", sep=""))
      date.gif<-date.sep(Datevalue,durdate)
      
      for (es_name in names(error_specs)) {
        
        err_tag <- error_specs[[es_name]]$tag
        tmp <- fread(paste0(datapath, "/", region[s], "/param/",Dname, "_SLDR_params_", err_tag, "_", Yname[yindex], ".csv"))
        tmp$msa     <- Dname
        tmp$err_tag <- err_tag
        tmp$yindex  <- Yname[yindex]
        tmp$region  <- region[s]
        
        # attach period
        tmp$day <- as.Date(tmp$day)
        tmp <- left_join(tmp, date.gif, by = "day")
        
        param_all <- rbind(param_all, tmp)
      }
    }
  }
}
param_all$period <- tools::toTitleCase(as.character(param_all$period))

#----------Plot win-rate of Exponential across losses----------#
# Aggregate daily errors within each (msa, period, kernel, err_tag)
sum_dt <- param_all[, .(err_med = median(error, na.rm = TRUE), err_mean = mean(error, na.rm = TRUE)),
                    by = .(msa, period, err_tag, index)]
# Winner kernel per (msa, period, err_tag)
winner_dt <- sum_dt[ , .SD[which.min(err_med)], by = .(msa, period, err_tag)]
winner_dt[, exp_win := (index == "exp")]
winrate_dt <- winner_dt[, .(n_cases = .N, exp_win_rate = mean(exp_win)), by = .(period, err_tag)]
loss_label_map <- c(
  rmse         = "RMSE",
  logrmse      = "Log-RMSE",
  srmspe       = "sRMSPE",
  rmse_rmspe   = "RMSE-RMSPE"
)
winrate_dt[, err_label := factor(
  err_tag,
  levels = names(loss_label_map),
  labels = loss_label_map
)]

fig_winrate <- ggplot(winrate_dt, aes(x = err_label, y = exp_win_rate)) +
  geom_col(fill = "#ad98c3", width = 0.65) +
  facet_wrap(~ period, scales = "fixed", nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Loss function", y = "Exponential win rate") +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 1), 
        strip.background = element_rect(fill = "#f0f0f0", color = "gray", linewidth = 1),
        strip.text = element_text(face = "plain", size=15),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right")
ggsave(fig_winrate, filename = paste(figpath,"/msa/SI_loss_winrate_",Yname[yindex],".pdf",sep=""), width = 4*2, height = 4*1.2)


# cross-city and intra city 
sum_dt[, kernel_rank := frank(err_med, ties.method = "average"), by = .(msa, period, err_tag)]
rank_wide <- dcast(
  sum_dt,
  msa + period + index ~ err_tag,
  value.var = "kernel_rank"
)%>%as.data.table
loss_tags <- setdiff(colnames(rank_wide), c("msa","period","index"))

u <- unique(rank_wide[, .(msa, period)])

rank_corr_dt <- rbindlist(
  lapply(seq_len(nrow(u)), function(i){
    
    msa_i    <- u[i, msa]
    period_i <- u[i, period]
    
    sub <- rank_wide[msa == msa_i & period == period_i]
    
    combs <- combn(loss_tags, 2, simplify = FALSE)
    
    rbindlist(lapply(combs, function(p){
      data.table(
        msa    = msa_i,
        period = period_i,
        loss_1 = p[1],
        loss_2 = p[2],
        spearman_rho = suppressWarnings(
          cor(sub[[p[1]]], sub[[p[2]]],
              method = "spearman",
              use = "complete.obs")
        )
      )
    }))
  }),
  fill = TRUE
)

rank_corr_dt[, loss_1 := factor(loss_1, levels = names(loss_label_map), labels = loss_label_map)]
rank_corr_dt[, loss_2 := factor(loss_2, levels = names(loss_label_map), labels = loss_label_map)]

rank_corr_summary <- rank_corr_dt[
  , .(
    median_rho = median(spearman_rho, na.rm = TRUE),
    IQR_low    = quantile(spearman_rho, 0.25, na.rm = TRUE),
    IQR_high   = quantile(spearman_rho, 0.75, na.rm = TRUE)
  ),
  by = .(period, loss_1, loss_2)
]




rho_loss_dt <- param_all[index == "exp",
                            .(rho_mean = mean(rho, na.rm = TRUE),
                              rho_med  = median(rho, na.rm = TRUE)),
                            by = .(msa, period, err_tag)
]
p_rho_loss <- ggplot(rho_loss_dt, aes(x = err_tag, y = rho_med, group = err_tag)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "gray30", linewidth = 0.3) +
  facet_wrap(~ period, nrow = 1) +
  labs(x = "Loss function", y = expression("Median of daily " * rho * " (Exponential kernel)")) +
  theme_wy() +
  theme(panel.border = element_rect(color="gray50", fill=NA, linewidth=0.5),
        strip.background = element_rect(fill="#f0f0f0", color="gray50", linewidth=0.5),
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(p_rho_loss, filename = file.path(figpath, paste0("robust_loss_rho_exp_", Yname[yindex], ".pdf")),
       width = 10, height = 4)



rho_delta_dt <- dcast(rho_loss_dt, msa + err_tag ~ period, value.var = "rho_med")
rho_delta_dt[, delta := During - Before]

p_delta <- ggplot(rho_delta_dt, aes(x = err_tag, y = delta)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "gray30", linewidth = 0.3) +
  labs(x = "Loss function", y = expression(Delta * rho~"(During - Before), Exponential kernel")) +
  theme_wy() +
  theme(panel.border = element_rect(color="gray50", fill=NA, linewidth=0.5),
        axis.text.x = element_text(angle = 20, hjust = 1))

