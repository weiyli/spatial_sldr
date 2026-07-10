# rm(list = ls())

# fig2_brief.R + fig4.R


#----------Workpath----------#
setwd("D:/ood/")
codepath <- 'D:/ood/Code/spatial_sldr/spatial_sldr'
geopath <- 'D:/ood/Data/Geo'
flowpath <- 'D:/ood/Data/Flow'
datapath.msa <- 'D:/ood/Data/spatial_sldr'
datapath.dis <- 'D:/ood/Data/spatial_sldr/disaster'
figpath <- 'D:/ood/Figure/spatial_sldr'

#----------Load packages----------#
library(ggrepel)     # geom_text_repel()
library(stringr)

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
Dname<<-Dname.set[d]
date.gif<-date.sep(Datevalue, durdate)


#----------plot SIR results----------#
E.msa <- F.msa <- list()
yindex <- 2
#----------data from the spreading_model.R----------#
sir.model <- read.csv(paste(datapath.msa,"/msa/SIR_model_",Yname[yindex],".csv",sep=""), header=TRUE)
sir.model$DM <- "Model"
sir.model[which(sir.model$det.rho==0),]$DM<-"Empirical"
sir.model.main <- subset(sir.model,beta==0.3&mu==0.1&msa%in%c("Atlanta","San Francisco"))
#----------plot infection peak and rho----------#
sir.model.peak <- as.data.table(sir.model.main)[, .SD[which.max(infected)], by = .(msa, day, rho, det.rho)]
baseline <- sir.model.peak[det.rho == 0, .(msa, day, time.base = time, infected.base = infected)]
sir.model.peak <- merge(sir.model.peak, baseline, by = c("msa", "day"), all.x = TRUE)
sir.model.peak[, ':='(
  time.pct =  (time - time.base) / time.base,
  infected.pct = (infected - infected.base) / infected.base
)]
E.msa[[yindex]] <- ggplot(sir.model.peak, aes(det.rho,round(infected.pct,2))) +  
  geom_line(aes(color = msa, group = msa), linetype = "dotted", na.rm = TRUE) +
  geom_point(aes(fill= msa, color = msa, shape = msa),size=3,stroke = 0.5) +
  geom_abline(intercept = 0, slope = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1, color = "#69b3a2") +
  scale_fill_manual(name = 'MSA',
                    breaks = msa.info$dname,
                    labels = msa.info$dname,
                    values = scales::alpha(msa.info$col,alpha=0.8))+
  scale_color_manual(name = 'MSA',
                     breaks = msa.info$dname,
                     labels = msa.info$dname,
                     values = msa.info$col) +
  scale_shape_manual(name = 'MSA',
                     breaks = msa.info$dname,
                     labels = msa.info$dname,
                     values = msa.info$shape) +
  scale_x_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.2),
                     labels = scales::percent_format(scale = 100)) +
  scale_y_continuous(limits = c(-1, 1),
                     breaks = seq(-1, 1, 0.2),
                     labels = scales::percent_format(scale = 100)) +
  labs(x=TeX('Percentage change in $\\rho$'), y='Percentage change in infection peak')+
  theme_wy() +
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.title = element_blank(), 
        legend.key.size = unit(0.5, "cm"),
        legend.position = "inside",
        legend.position.inside = c(0.05,0.85),
        legend.text = element_text(hjust=0),
        legend.justification = c(0,0.5),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank())
#----------plot infected pop and time----------#
sir.model.time <- subset(sir.model.main,det.rho%in%c(-0.5,0,0.5))
det.rho <- sort(unique(sir.model.time$det.rho),decreasing=TRUE)
rho.breaks <- rho.labels <- paste0(round(det.rho,4)*100,"%")
rho.labels[which(rho.breaks=="0%")] <- "0% empirical"
rho.info <- data.frame(det.rho = det.rho, 
                       breaks = rho.breaks,
                       labels = rho.labels,
                       color = col.fit[c(5,3,1)],
                       line = c("dotted","solid","dotdash"))
# color = col.all[1:length(det.rho)])
sir.model.time <- left_join(sir.model.time, rho.info, by="det.rho")
emp.rho <- as.data.table(sir.model.time)[DM == "Empirical", .SD[which.max(infected)], by = msa]
#----------plot infected rate for two msas with (beta=0.3,mu=0.1)----------#
F.msa[[yindex]] <- ggplot(data = sir.model.time, aes(x = time, y = infected/pop)) +
  geom_line(aes(color = labels, linetype = labels, group = rho), linewidth = 0.8) +
  geom_label_repel(data = emp.rho,
                   aes(label = paste("empirical",round(rho, 4),sep=" ")), 
                   # aes(label = sprintf("rho == %.4f", round(rho, 4))), 
                   max.overlaps = Inf,
                   box.padding = 1, size = 3, fill = 'white',color="black") +
  facet_wrap(~ msa, ncol = 2, scales = "fixed") + 
  scale_colour_manual(name = TeX('% change in $\\rho$'),
                      breaks = rho.info$labels,
                      labels = rho.info$labels,
                      values = rho.info$col) +
  scale_linetype_manual(name = TeX('% change in $\\rho$'), 
                        breaks = rho.info$labels,
                        labels = rho.info$labels,
                        values = rho.info$line) +
  labs(x = "Time", y = "Infected fraction", title = NULL) +
  theme_wy() +
  theme(panel.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        panel.border = element_rect(color = "gray", fill = NA, linewidth = 0.5), 
        panel.grid.major.x = element_line(color = "gray90", linetype = "dashed"),
        panel.grid.major.y =  element_line(color = "gray90", linetype = "dashed"),
        strip.background = element_rect(fill = "gray90", color = "gray", linewidth = 0.5),
        strip.text = element_text(face = "plain", size = 15, color = event.info$col[d]),
        legend.position = "right",
        legend.justification = c(0,0.5))
fig.rho.sir <- (E.msa[[yindex]]|F.msa[[yindex]]) + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 20))
ggsave(fig.rho.sir,filename = paste(figpath,"/msa/SI_sir_rho_pop_",Yname[yindex],".pdf",sep=""), width = 6.8*2, height = 5)






