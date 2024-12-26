library(ggplot2)
library(see)
library(ggpubr)
library(ape)
library(phytools)
library(rr2)
library(nlme)
library(caper)
library(geiger)
library(grid)
library(Cairo)
library(qpdf)
library(grid)
library(dplyr)
library(cowplot)
library(magick)
library(grid)
library(gridExtra)
library(deeptime)
library(devtools)
library(ggeasy)
library("ggsci")
library(vegan)
library(factoextra)
library(FactoMineR)
library(PCAmixdata)
library(png)
library(RPANDA)
library(geoscale)
library(patchwork)
library(pairwiseAdonis)
library(ecole)
library(magrittr) # for piping
library(ggtree)
library(paleotree)
library(ggstatsplot)


# PLEASE CHECK THE FOLLOWING TWO LINES
# if you are using the "Rscript" command in linux/mac terminal
setwd(paste0(getwd(),"/../data/all_trees"))

# if you are using RStudio IDE
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../data/"))

# OR DIRECTLY PUT THE PATH AS A STRING

data <- read.csv("Final_summary_updated2.csv")

head(data)
nrow(data)


data.taxa.modified <- data
data.taxa.modified$Taxonomic.group <- factor(data.taxa.modified$Taxonomic.group,levels=c("Invertebrate","Herpetofauna","Angiosperm"))

geological_events <- c("GS",
                       "DV+SIS",
                       "LC-w-Asia",
                       "HCA+EOC",
                       "Peak H-T-O",
                       "IA+IM+\nC4-Ex")
geological_events_windows <- matrix(c(170,150,
                                      66,59.9,
                                      48,39.9,
                                      39,31.9,
                                      19,15,
                                      11,2.9), nrow=length(geological_events),ncol=2,byrow=TRUE)
geological_events_windows
geological_events_ends <- geological_events_windows[,2]
geological_events_starts <- geological_events_windows[,1]
geological_events_dates <- (geological_events_starts+geological_events_ends)/2


######################################## Representation of data ###########################################################
################ Lineages, Clade ages, Endemicity, Habitat types, Biogeographic orgins, Peninsular India ##################
data.modified <- c()

tax.grp <- c("Invertebrate","Herpetofauna","Angiosperm")
num.tax.grp <- length(tax.grp)
for(i in 1:num.tax.grp){
  #i=1
  data.tax <- data[which(data$Taxonomic.group==tax.grp[i]),]
  data.tax <- data.tax[order(data.tax$Stem.age..Mya.,decreasing = TRUE),]
  data.modified <- rbind(data.modified,data.tax)
}
data.modified$Lineage <- c(1:length(data.modified$Lineage))

data.modified$Lineage <- factor(data.modified$Lineage, levels = data.modified$Lineage)


pl1 <- ggplot(data.modified, aes(x=Lineage,y=Crown.age..Mya.,fill=Taxonomic.group))+
  scale_x_discrete(position = "top",limits = as.factor(c(1:40))#,breaks=as.factor(c(5,8,16,28,29,34))
                   )+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[1],ymax=geological_events_ends[1]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[1],label=geological_events[1]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[2],ymax=geological_events_ends[2]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[2],label=geological_events[2]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[3],ymax=geological_events_ends[3]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[3]),label=geological_events[3],
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[4],ymax=geological_events_ends[4]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[4],label=geological_events[4]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[5],ymax=geological_events_ends[5]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[5],label=geological_events[5]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[6],ymax=geological_events_ends[6]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[6],label=geological_events[6]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  scale_fill_manual(values = c("#77dd77","#898989","#00c7e7"))+
  scale_color_manual(values = c("#77dd77","#898989","#00c7e7"))+
  geom_bar(stat="identity",width=0.8,alpha=0.95,
           #color="black",
           size=0.2)+geom_linerange(aes(ymin=Crown.age..Mya.,ymax=Stem.age..Mya.,colour=Taxonomic.group))+
  theme_classic()+
  geom_hline(yintercept = 0,size=0.25)+
  scale_y_reverse(breaks=c(200,175,150,125,100,75,50,25,0))+theme(axis.text.x=element_text(size=14),
                          legend.position = "none",
                          plot.title = element_blank(),
                        axis.title.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.line.y = element_blank(),
                        axis.line.x = element_blank(),
                        axis.text.y=element_blank(),
                        axis.title.x=element_text(size=16))+xlab("Lineages")+ylab("Age (Mya)")+ggtitle("Clade ages")+
  geom_text(aes(x = as.factor(37),y=-18,label="Endemicity\n(Endemic Species\n/Species Richness)"),#,"Endemic\nSpecies","Biogeographic\nOrigins","Habitat\nType")),
            size = 3.5,fontface="bold",check_overlap = TRUE)+
  geom_text(aes(x = as.factor(37),y=-48,label="Biogeographic\nOrigin"),
            size = 3.5,fontface="bold",check_overlap = TRUE)+
  geom_text(aes(x = as.factor(37),y=-78,label="Habitat\nType"),
            size = 3.5,fontface="bold",check_overlap = TRUE)+
  geom_text(aes(x=Lineage,y=-18,
                  #(Stem.age..Mya.+8),  
                label=paste0(Endemic.species,"/",No..of.species)),size = 3.0,color="#000000")+
  geom_text(aes(x=Lineage,y=-48,
                #(Stem.age..Mya.+8),  
                label=Biogeographic.origin),size = 3.0,color="#000000")+
  geom_text(aes(x=Lineage,y=-78,
                #(Stem.age..Mya.+8),  
                label=Habitat),size = 3.0,color="#000000")+coord_flip()

plg <- gggeo_scale(pl1)
CairoPDF("../graphs/1_Data_biogeo_clade_age_habitat.pdf",11.70,8.30,bg="transparent")
plg
dev.off()

###### Figure 2: PhyloANOVA using clade rates #######

############################################################################################################################
############################################## DIVERSIFICATION PATTERNS #########################################################
data.modified <- c()

tax.grp <- c("Invertebrate","Herpetofauna","Angiosperm")
num.tax.grp <- length(tax.grp)
for(i in 1:num.tax.grp){
  #i=1
  data.tax <- data[which(data$Taxonomic.group==tax.grp[i]),]
  data.tax <- data.tax[order(data.tax$Crown.age..Mya.,decreasing = TRUE),]
  data.modified <- rbind(data.modified,data.tax)
}
data.modified$Lineage <- c(1:length(data.modified$Lineage))
#View(data.modified)
data.modified$Lineage <- factor(data.modified$Lineage, levels = data.modified$Lineage)

data.modified$Driver[which(data.modified$Driver=="")] <- "Other"
data.modified$Final.Results[which(data.modified$Final.Results == "SC1")] <- "SC1\n[Gradual Accumulation]"
data.modified$Final.Results[which(data.modified$Final.Results == "SC2")] <- "SC2\n[Saturated Accumulation]"
data.modified$Final.Results[which(data.modified$Final.Results == "SC2+Shift")] <- "SC2+Shift\n[Saturated Accumulation\n+Episodic Rate Shift(s)]"
data.modified$Final.Results[which(data.modified$Final.Results == "SC3+Shift")] <- "SC3+Shift\n[Waxing and Waning\n+Episodic Rate Shift(s)]"
data.modified$Final.Results[which(data.modified$Final.Results == "Shift")] <- "Shift\n[Gradual Accumulation\n+Episodic Rate Shift(s)]"

data.modified$Final.Results <- factor(data.modified$Final.Results, levels = c("Shift\n[Gradual Accumulation\n+Episodic Rate Shift(s)]",
                                                                              "SC3+Shift\n[Waxing and Waning\n+Episodic Rate Shift(s)]",
                                                                              "SC2+Shift\n[Saturated Accumulation\n+Episodic Rate Shift(s)]",
                                                                              "SC2\n[Saturated Accumulation]",
                                                                              "SC1\n[Gradual Accumulation]"))
pl_scenarios_drivers <- ggplot(data.modified, aes(x=Final.Results,fill=Driver))+
  scale_fill_manual(values = c("darkorange","black","darkgreen","gray"))+
  geom_bar(width=0.5)+theme_classic()+
  ylim(c(0,25))+
  ggtitle("Diversification scenarios across lineages")+xlab("Diversification Scenarios")+ylab("No. of lineages")+
  theme(axis.text = element_text(colour = "black"),
        axis.text.y=element_text(size=8,hjust=0.5),
        axis.text.x=element_text(size=8),
        axis.title=element_text(size=12,face="bold"),legend.text=element_text(size=6),
        legend.position = c(0.8,0.2),legend.key.size = unit(4, 'mm'),legend.title=element_text(size=7),
        legend.box.background = element_rect(colour = "black",size=0.1),
        legend.box.margin = margin(6, 10, 6, 6),plot.title = element_text(hjust = 0.5,face="bold"))+coord_flip()

data.for.episodic.rate.shift <- data[union(union(which(!is.na(data$Rate_Shift_Time.CoMET.)),which(!is.na(data$Peak_or_Dip_Time.SES.))),which(!is.na(data$Rate_Shift_Time_2.CoMET.))),]
data.for.episodic.rate.shift$type <- rep(NA,nrow(data.for.episodic.rate.shift))
data.for.episodic.rate.shift$type[which(!is.na(data.for.episodic.rate.shift$Peak_or_Dip_Time.SES.))] <- "SES\nmodel"
ind.comet <- union(which(!is.na(data.for.episodic.rate.shift$Rate_Shift_Time.CoMET.)),
                   which(!is.na(data.for.episodic.rate.shift$Rate_Shift_Time_2.CoMET.)))
data.for.episodic.rate.shift$type[ind.comet[which(ind.comet %in% which(!is.na(data.for.episodic.rate.shift$Peak_or_Dip_Time.SES.)))]] <- "CoMET+\nSES model"
data.for.episodic.rate.shift$type[ind.comet[which(!(ind.comet %in% which(!is.na(data.for.episodic.rate.shift$Peak_or_Dip_Time.SES.))))]] <- "CoMET\nanalysis"

data.for.episodic.rate.shift <- data.for.episodic.rate.shift[,c("Rate_Shift_Time_2.CoMET.","Rate_Shift_Time.CoMET.","Peak_or_Dip_Time.SES.","type","Taxonomic.group")]

data.for.episodic.rate.shift %>% 
  pivot_longer(c(Rate_Shift_Time.CoMET.,Peak_or_Dip_Time.SES.,Rate_Shift_Time_2.CoMET.),names_to = "rate_shift_type",values_to = "rate_shift_time") -> data.for.episodic.rate.shift
unique(data.for.episodic.rate.shift$rate_shift_type)
data.for.episodic.rate.shift$rate_shift_type[which(data.for.episodic.rate.shift$rate_shift_type == "Rate_Shift_Time.CoMET.")] <- "1st Rate Shift"
data.for.episodic.rate.shift$rate_shift_type[which(data.for.episodic.rate.shift$rate_shift_type == "Rate_Shift_Time_2.CoMET.")] <- "2nd Rate Shift"
data.for.episodic.rate.shift$rate_shift_type[which(data.for.episodic.rate.shift$rate_shift_type == "Peak_or_Dip_Time.SES.")] <- "Peak or Dip"

p.episodic <- ggplot(data.for.episodic.rate.shift,aes(x=type,y=rate_shift_time))+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[6],ymax=geological_events_ends[6]),
            alpha = 0.06,
            fill = "#cce7e8")+
  ggtitle("Estimated times of strong episodic rate shifts")+
  geom_text(aes(x = as.factor(""),y=geological_events_dates[6],label=geological_events[6]),
            size = 2.5,angle=90,check_overlap = TRUE)+
  geom_boxplot(outliers = F,width=0.7,fill=NA,lwd=0.2)+
  xlab("Method")+ylab("Time (Mya)")+
  # geom_hline(yintercept = 0,size=0.35)+
  geom_jitter(aes(fill=Taxonomic.group,shape=rate_shift_type),position=position_jitter(0.4),size=2,stroke=0)+
  scale_shape_manual(values=c(24,25,21))+
  scale_y_reverse(limits=c(40,0),breaks=c(40,30,20,10,0))+
  #scale_x_discrete(position = "top")+
  #geom_vline(xintercept = 0,size=0.25)+
  scale_fill_manual(values = c("#77dd77","#898989","#00c7e7"))+
  guides(fill=guide_legend(title="Taxonomic groups",override.aes = list(color=c("#77dd77","#898989","#00c7e7"))),
         shape=guide_legend(title="Rate-shift type",override.aes = list(fill="black")))+
  coord_flip()+theme_classic()+theme(legend.position = c(0.25,0.18),
                                     legend.direction = "vertical",
                                     legend.background = element_blank(),
                                     legend.box = "horizontal",
                                     legend.key.size = unit(4, 'mm'),legend.title=element_text(size=7),
                                     legend.text = element_text(size=6),
                                     legend.box.background = element_rect(colour = "black",size=0.1),
                                     legend.box.margin = margin(6, 10, 6, 6),
                                     axis.title = element_text(size=12,face="bold"),
                                     # axis.ticks.y = element_blank(),
                                     axis.title.y =  element_text(size=12,margin = margin(r = 15),colour = "black",face="bold"),
                                     axis.text.y = element_text(size=8,hjust = 1,vjust=0.5,colour = "black"),
                                     axis.text.x = element_text(size=8,colour = "black"),
                                     axis.line.x = element_blank(),plot.title = element_text(hjust = 0.5,face="bold"))
p.episodic <- gggeo_scale(p.episodic,abbrv = F,size = 2.5,height = unit(1, "line"))
# pdf("../graphs/3_Rate_shifts2.pdf",11.7,8.3,bg="transparent")
# ggdraw()+draw_plot(p.episodic)+draw_plot(pl_scenarios_drivers,height=1,x=-0.265,y=0.23,scale=0.4)
# dev.off()

rownames(geological_events_windows) <- geological_events
colnames(geological_events_windows) <- c("older bound","newer bound")
geological_events_windows <- as.data.frame(geological_events_windows)

tax.grp.SC1 <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$Final.Results=="SC1")])
tax.grp.SC2 <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$Final.Results=="SC2")])
tax.grp.SC2_Shift <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$Final.Results=="SC2+Shift")])
tax.grp.SC3_Shift <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$Final.Results=="SC3+Shift")])
tax.grp.Shift <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$Final.Results=="Shift")])
tax.grp.TDD <- as.vector(data.taxa.modified$Taxonomic.group[which(data.taxa.modified$TDD=="TDD")])


make_beautiful_RTT <- function(filenames,tax.grp,text.title,ylb=T){
  # filenames <- c(filenames.SC2,filenames.SC2_Shift,filenames.SC3_Shift,filenames.Shift)
  
  
  # tax.grp <- c(tax.grp.SC2,tax.grp.SC2_Shift,tax.grp.SC3_Shift,tax.grp.Shift)
  # text.title <- "Time-Varying Rates"
  # ylb = F
  # filenames = filenames.SC1;tax.grp = tax.grp.SC1;text.title="Constant Rates"
  cat("Started\n\nGetting geological information for plotting.....\n\n")
  
  cols <- c("#00ff00","#000000","#0000ff")
  names(cols) <- c("Angiosperm","Herpetofauna","Invertebrate")
  
  cols <- cols[which(names(cols) %in% unique(tax.grp))]
  
  long.format.data <- matrix(ncol=4)
  colnames(long.format.data) <- c("Lineage ID","Net-Diversification Rate (events/Myr/lineage)","Time (Mya)","Taxonomic group")
  for(i in 1:length(filenames)){
    #i <- 6
    rates <- as.vector(unlist(read.csv(filenames[i])[1,-1]))
    mean.rate <- mean(rates)
    rates <- rates/mean.rate
    long.format.data <- rbind(long.format.data,
                              cbind(i,
                                    c(rates,rates[length(rates)]),
                                    c(as.numeric(substr(colnames(read.csv(filenames[i])),2,nchar(colnames(read.csv(filenames[i])))))[-1],0), 
                                    tax.grp[i])
                              )
  }
  
  long.format.data <- as_tibble(long.format.data[-1,])
  
  
  if(text.title == "Extinction rates"){
    plot_title <- text.title
  }else{
    plot_title <- paste0("RTT for Lineages with ",text.title)
  }
  
  timebins <- seq(0,150,5)
  long.format.data$time.bin <- rep(NA,nrow(long.format.data))
  for(i in 2:(length(timebins))){
    #i <- 2
    ind.less.than.max.of.this.time.bin <- which(as.numeric(long.format.data$`Time (Mya)`)<timebins[i])
    ind.great.eq.than.min.of.this.time.bin <- which(as.numeric(long.format.data$`Time (Mya)`)>=timebins[i-1])
    ind.for.this.time.bin <- intersect(ind.less.than.max.of.this.time.bin,ind.great.eq.than.min.of.this.time.bin)
    long.format.data$time.bin[ind.for.this.time.bin] <- mean(c(timebins[i],timebins[i-1]))
  }
  
  #usinf ggplot
  #y.lab <- ifelse(ylb,"Net-Diversification Rate\n(events/Myr/lineage)","")
  y.lab <- ifelse(ylb,"Net-Diversification Rate\nnormalised by lineage mean","")
  
  pl1 <- ggplot(long.format.data, aes(x = time.bin, y = as.numeric(`Net-Diversification Rate (events/Myr/lineage)`), group = time.bin)) + 
    ggtitle(plot_title)+xlab('Time (Mya)')+ylab(y.lab)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[1],xmax=geological_events_ends[1]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[1],label=geological_events[1]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[2],xmax=geological_events_ends[2]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[2],label=geological_events[2]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[3],xmax=geological_events_ends[3]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[3]),label=geological_events[3],
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[4],xmax=geological_events_ends[4]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[4],label=geological_events[4]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[5],xmax=geological_events_ends[5]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[5],label=geological_events[5]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[6],xmax=geological_events_ends[6]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 4.2,x=geological_events_dates[6],label=geological_events[6]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_text(aes(y = 0.3,x=-5,label="PRESENT DAY"),
              size = 2.5,angle=90,check_overlap = TRUE,color="black")+
    geom_boxplot(width=4.2,outliers = T,outlier.size = 0.1, lwd=0.1,alpha=0.5)+
    #geom_jitter(shape=16,size=0.05,aes(colour = `Taxonomic group`),alpha=0.5)+
    theme_classic()+
    #scale_fill_manual("Taxonomic group",values = cols)+
    scale_color_manual("Taxonomic group",values = cols)+
    #geom_vline(xintercept = 0,size=0.25)+
    # guides(color = guide_legend( 
    #   override.aes=list(shape = 20,size=3,stroke=0)))+
    # geom_segment(aes(x = 0, xend = Inf, y = Inf), 
    #              linetype = "solid", color = "black", size = 0.1)+
    #ylim(c(min(as.numeric(long.format.data$`Net-Diversification Rate (events/Myr/lineage)`)),0.43))+
    ylim(c(-0.02),5)+
    scale_x_reverse(breaks=c(150,125,100,75,50,25,0),limits=c(150,0))+theme(
      # legend.position = c(0.2,0.55),legend.key.size = unit(0.1, 'cm'),
      #                                                                       legend.text = element_text(size = 7,colour = "black"),
      #                                                                       legend.title = element_text(size = 8,colour = "black"),
                                                                            axis.line.x = element_blank(),
                                                                            axis.title = element_text(size=12, face="bold"),
                                                                            axis.text = element_text(color="black",size=8),
                                                                            plot.title = element_text(hjust = 0.5,face="bold"))
  
  gggeo_scale(pl1,abbrv = F,size = 2.5,height = unit(1, "line"))
  
  
}

# filenames for RTT Plots 
filenames.SC1 <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC1")],"_nd_rates.csv")
filenames.SC2 <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2")],"_nd_rates.csv")
filenames.SC2_Shift <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2+Shift")],"_nd_rates.csv")
filenames.SC3_Shift <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC3+Shift")],"_nd_rates.csv")
filenames.Shift <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="Shift")],"_nd_rates.csv")
filenames.TDD <- paste0("./estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$TDD=="TDD")],"_nd_rates.csv")



p1 <- make_beautiful_RTT(c(filenames.SC2,filenames.SC2_Shift,filenames.SC3_Shift,filenames.Shift),
                   c(tax.grp.SC2,tax.grp.SC2_Shift,tax.grp.SC3_Shift,tax.grp.Shift),
                   "Time-Varying Rates")

p2 <- make_beautiful_RTT(filenames.SC1,tax.grp.SC1,"Constant Rates")

CairoPDF("../graphs/3_Diversification_patterns.pdf",11.7,8.3,bg="transparent")
grid.arrange(pl_scenarios_drivers,p.episodic,p2,p1,nrow=2,ncol=2)
dev.off()
#########################################################################################################################################



##########################################################################################################################
########################################## Multivariate-analysis using FAMD ##################################################
############################################SCENARIOS-DRIVERS, DISTRIBUTION ##########################################
# input a matrix with - Lineage names, taxonomic affiliations, clade ages, species richness, net diversification rate,
# biogeographic origin, habitat, if clades are habitat-structured, diversification patterns, 
# if it shows Diversity-dependent diversification, if it shows Temperature-dependent diversification
# but make sure the lineage names column comes first in the matrix

mk_FAMD_plot <- function(df,list.variables,palette){
  #df=dt;list.variables= c("Habitat Types");palette=c("grey10","grey40","grey70")
  lin.names <- df$Lineage
  df <- df[,-which(colnames(df)=="Lineage")]
  rownames(df) <- lin.names
  for(i in c("Himalayan-Orogeny Dependence","Temperature Dependence","Diversity Dependence","Scenarios","Habitat Types","Biogeographic Origins","Taxonomic groups")){
    df[,i] <- as.factor(df[,i])
  }
  
  diss.matrix <- kmed::distmix(df,method = "gower",idnum = c(2,3,4),idbin = NULL,idcat = c(1,5:10))
  
  res.famd <- FAMD(df,graph = FALSE)
  famd_coordinates <- res.famd$ind$coord
  
  # Calculate dissimilarity matrix. we are using mixed data and hence we are using Gower distance
  #diss.matrix <- vegdist(famd_coordinates, method = "gower")
  mat <- as.matrix(diss.matrix)
  #View(mat)
  tree <- read.tree("SUPER_TREE.nwk")
  tree <- drop.tip(tree,"Montecincla")
  phylogenetic.diss.mat <-cophenetic(tree)
  phylogenetic.diss.mat <- phylogenetic.diss.mat/max(phylogenetic.diss.mat)
  mantel.test <- vegan::mantel(xdis = mat,
                ydis = phylogenetic.diss.mat,
                method = "pearson",
                permutations = 999)
  cat("\nSignificance of correlation between phylogenetic distance and trait distance: ",mantel.test$signif,"\n")
  weighted.diss.mat <- mat*phylogenetic.diss.mat
  #View(weighted.diss.mat)
  weighted.diss.mat <- as.dist(weighted.diss.mat)
  
  p_scree <- fviz_screeplot(res.famd)
  
  # contribution to the 1st Dimension
  pc1 <- fviz_contrib(res.famd,"var",axes=1)
  # contribution to the 2nd Dimension
  pc2 <- fviz_contrib(res.famd,"var",axes=2)
  # visualize quantitative variables
  plot.quanti <- fviz_famd_var(res.famd,"quanti.var",repel=TRUE,col.var = "contrib",gradient.cols=c("red","white","blue"))
  
  groups.variables <- c(rep("A.Taxonomic group",3),
                        rep("B. Biogeographic origin",2),
                        rep("C. Habitat Type",3),
                        #rep("D. Habitat structuring",2),
                        rep("E. Scenarios",5),
                        rep("F. Drivers",6))
  
  # visualize qualitative variables
  plot.quali <- fviz_famd_var(res.famd,geom=c("arrow","text"),"quali.var",col.var = factor(groups.variables),
                              palette="jco",repel = TRUE)
  
  CairoPDF("../graphs/supplementary_FAMD.pdf",8.3,11.7,bg="transparent")
  grid.arrange(p_scree,pc1,pc2,nrow=3)
  grid.arrange(plot.quanti,plot.quali,nrow=2)
  dev.off()
  
  if(is.null(list.variables)){
    # plot of variables
    fviz_famd_var(res.famd,geom=c("arrow","text"),repel=TRUE)+labs(title = NULL)+
      theme(axis.title = element_text(size=12,face="bold"),axis.text = element_text(color="black",size=8))
  }else{
    perm.var <- permanova_pairwise(diss.matrix,df[,list.variables[1]],padj = "bonferroni")
    perm.var.phylo <- permanova_pairwise(weighted.diss.mat,df[,list.variables[1]],padj = "bonferroni")
    
    name <- unlist(strsplit(list.variables[1],"\n"))
    
    if(!("PermANOVA_results" %in% list.files("../graphs/"))){
      dir.create("../graphs/PermANOVA_results")
    }
    write.csv(perm.var,paste0("../graphs/PermANOVA_results/Permanova_",name,".csv"))
    write.csv(perm.var.phylo,paste0("../graphs/PermANOVA_results/Phylogeny_weighted_Permanova_",name,".csv"))
    # graph of lineages/individual points
    return(fviz_ellipses(res.famd,list.variables,
                  ellipse.level = 0.95, palette = palette,ellipse.type = c("convex"),repel=TRUE,geom="point",labelsize = 1,mean.point.size = 0)+
      #annotation_custom(permanova.table, xmin=-3.8, xmax=-1, ymin=3, ymax=7.5)+
      theme(axis.title = element_text(size=12,face="bold"),axis.text = element_text(color="black",size=8)))
    
    
    
    
    
  }
  
}

dt <- data.taxa.modified[,c("Lineage","Taxonomic.group","Crown.age..Mya.",
                            "No..of.species","Net.diversification.rate..ClaDS.",
                            "Biogeographic.origin","Habitat",#"Habitat.differentiation.among.clades.",
                            "Final.Results",
                            "DDD",
                            "TDD",
                            "Him.oro"
                            # ,"C4.Exp"
)]


colnames(dt) <- c("Lineage","Taxonomic groups","Clade Age",
                  "Species Richness","Net-diversification Rate",
                  "Biogeographic Origins","Habitat Types",#"Habitat Structuring",
                  "Scenarios","Diversity Dependence","Temperature Dependence","Himalayan-Orogeny Dependence"
                  # ,"C4-Plants-Expansion\nDependence"
)


# Primary plot defining the diversification space in terms of all variables and how variables are correlated
pl_all <- mk_FAMD_plot(dt,list.variables = NULL,NULL)
pl_FAMD_SC <- mk_FAMD_plot(dt,c("Scenarios"),"jco")+guides(color=guide_legend(ncol=3,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
pl_FAMD_Tax <- mk_FAMD_plot(dt,c("Taxonomic groups"),
                            c("#77dd77","#898989","#00c7e7"))+
  guides(color=guide_legend(ncol=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
pl_FAMD_Hab <- mk_FAMD_plot(dt,c("Habitat Types"),c("grey10","grey40","grey70"))+
  guides(color=guide_legend(ncol=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
pl_FAMD_Biog <- mk_FAMD_plot(dt,c("Biogeographic Origins"),c("grey10","grey40"))+
  guides(color=guide_legend(nrow=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
# pl_FAMD_HS <- mk_FAMD_plot(dt,c("Habitat Structuring"),c("grey10","grey40"))+
#   guides(color=guide_legend(ncol=2,title.position = "top"))+
#   theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
#         plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
#                            "inches"))+labs(title = "")
pl_FAMD_DDD <- mk_FAMD_plot(dt,c("Diversity Dependence"),"jco")+
  guides(color=guide_legend(nrow=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
pl_FAMD_TDD <- mk_FAMD_plot(dt,c("Temperature Dependence"),"jco")+
  guides(color=guide_legend(nrow=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
# pl_FAMD_HO_DD <- mk_FAMD_plot(dt,c("Himalayan-Orogeny Dependence"),"jco")+
#   guides(color=guide_legend(nrow=2,title.position = "top"))+
#   theme(legend.position="top",legend.text=element_text(size=8),legend.title=element_text(size=9),legend.title.align=0.5)+labs(title = "")
# pl_FAMD_C4_DD <- mk_FAMD_plot(dt,c("C4-Plants-Expansion\nDependence"),"jco")+
#   guides(color=guide_legend(nrow=2))+
#   theme(legend.position="top",legend.text=element_text(size=8),legend.title=element_text(size=9))+labs(title = "")
# Plot to show how patterns and other variables are distributed and structured across the diversification space
pdf("../graphs/4_scenarios_drivers.pdf",8.3,11.7,bg="transparent")
grid.arrange(pl_all,arrangeGrob(pl_FAMD_SC,pl_FAMD_Tax,pl_FAMD_Hab,
                                pl_FAMD_Biog,pl_FAMD_DDD,pl_FAMD_TDD,nrow=2,ncol=3),nrow = 2)
dev.off()
write.csv(dt,"../graphs/FAMD_input_matrix.csv",row.names = F)
##############################################################

mk_pgls_plot <- function(formula,phy,df.pgls,xlab,ylab,xlim.min,ylim.min){
  # formula=Speciation.rate..RPANDA.~Crown.age..Mya.;phy = tree;df.pgls = data.taxa.modified;xlab = "Clade age (Mya)";ylab="";xlim.min = 0;ylim.min = 0
  # formula=Extinction.rate..ClaDS.~Crown.age..Mya.;phy = tree;df.pgls = data.taxa.modified;xlab = "Clade age (Mya)";ylab="";xlim.min = 0;ylim.min = 0
  # formula=Pulled.Diversification.Rate~Crown.age..Mya.;phy = tree;df.pgls = data.taxa.modified3;xlab = "";ylab="";xlim.min = 0;ylim.min = -0.03
  Response <- as.character(formula)[2]
  Predictor <- as.character(formula)[3]
  
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  
  df.pgls <- df.pgls[match(phy$tip.label,df.pgls$Lineage),] # order data in the same format as tip labels of tree
  
  comp.data <- comparative.data(phy=phy,data = df.pgls,names.col = Lineage,vcv = TRUE, 
                                na.omit = FALSE, warn.dropped = TRUE)
  pgls.model <- pgls(Response~Predictor,data = comp.data,lambda = 'ML')
  
  summ.pgls <- summary(pgls.model)
  lambda.r <- summ.pgls$param["lambda"]
  coeff.r <- as.data.frame(summ.pgls$coefficients)
  p.val.r <- coeff.r$`Pr(>|t|)`[-1]
  names(p.val.r) <- rownames(coeff.r)[-1]
  intercept.r <- coeff.r[1,"Estimate"]
  slope.r <- coeff.r[2,"Estimate"]
  
  r2.r <- summ.pgls$r.squared
  stats <- as.character(as.expression(substitute(italic("R")^2~"="~r2*","~italic("p")~"="~p,
                                                 list(r2 = round(r2.r,4),
                                                      p = round(unname(p.val.r[1]),5)
                                                 ))))
  if(is.na(lambda.r)){
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = "NA"))))
  }else{
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = round(unname(lambda.r)[1],3)))))
  }
  
  if(p.val.r > 0.05){
    lty<-3
  }else{lty<-1}
  colnames(df.pgls)[which(colnames(df.pgls) == Response)] <- "Response"
  colnames(df.pgls)[which(colnames(df.pgls) == Predictor)] <- "Predictor"
  
  get_labels <- function(x,lab) {
    if(!(-1 %in% unlist(gregexpr("log",lab)))){
      exp(x)
    }else{
      x
    }
    
  }
  get_breaks_custom <- function(vals,lab){
    # vals <- df.pgls$Response;lab=""
    if(!(-1 %in% unlist(gregexpr("log",lab)))){
      log(c(1,2,5,10,20,50,100,200))
    }else{
      max.val <- max(vals)
      min.val <- min(vals)
      interval <- max.val/5
      if(min.val>1){
        max.val <- round(max.val/5,0)*5
        interval <- round(interval,0)
      }else{
        if(interval < 0.05){
          max.val <- round(max.val/5,2)*5
          interval <- round(interval,2)
        }else{
          max.val <- round(max.val/5,1)*5
          interval <- round(interval,1)
        }
      }
      seq(0,max.val,interval)
    }
  }
  
  breaks.x <- get_breaks_custom(df.pgls$Predictor,xlab)
  breaks.y <- get_breaks_custom(df.pgls$Response,ylab)
  plot <- ggplot(df.pgls,aes(x=Predictor,y=Response)) + 
    geom_point(
      # aes(color=Biogeographic.origin,
      #              fill=Taxonomic.group,
      #              shape=Habitat),
               alpha=0.6,size=1,stroke=0.7) +
    # scale_fill_manual(values = c(alpha("#00ff00",0.8),alpha("#000000",0.8),alpha("#0000ff",0.8)),
    #                   guide=guide_legend(override.aes = list(shape = 21,stroke=0,
    #                                                          color=c(alpha("#77dd77",0.8),alpha("#898989",0.8),alpha("#00c7e7",0.8)))))+
    # scale_shape_manual(values = c(21,22,24))+
    # scale_color_manual(
    #   values = c("#aa4231", "#222233"),
    #   guide = guide_legend(override.aes = list(shape = 21))
    # )+
    scale_x_continuous(labels = get_labels(breaks.x,lab = xlab), breaks = breaks.x,limits = c(xlim.min,(max(df.pgls$Predictor)*1.1)))+
    scale_y_continuous(labels = get_labels(breaks.y,lab = ylab), breaks = breaks.y,limits = c(ylim.min,(max(df.pgls$Response)*1.3)))+
    xlab(ifelse(!(-1 %in% unlist(gregexpr("log",xlab))),substr(xlab,(unlist(gregexpr("log",xlab))+4),(nchar(xlab)-1)),xlab))+
    ylab(ifelse(!(-1 %in% unlist(gregexpr("log",ylab))),substr(ylab,(unlist(gregexpr("log",ylab))+4),(nchar(ylab)-1)),ylab))+
    geom_abline(intercept = intercept.r, slope = slope.r,linetype=lty,size=0.25) + 
    # geom_text(data = as.data.frame(eq), aes(label = eq,x = max(df.pgls$Predictor)/1.75,
    #                                         y = (max(df.pgls$Response))*1.45,family = 'serif'), parse = TRUE,size=3.5)+
    geom_text(data = as.data.frame(stats), aes(label = stats,x = max(df.pgls$Predictor)/1.5,
                                               y = (max(df.pgls$Response)),family = 'serif'), inherit.aes = FALSE, parse = TRUE,size=2.5)+
    geom_text(data = as.data.frame(lambda), aes(label = lambda,x = max(df.pgls$Predictor)/1.5,
                                                y = (max(df.pgls$Response))/1.25,family = 'serif'), inherit.aes = FALSE, parse = TRUE,size=2.5)+
    theme_classic()+theme(legend.title = element_blank(),axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"))
  return(plot)

  
}



##################### SUPPLEMENTARY FIGURES ################
# supplementary data representation
CairoPDF("../graphs/supplementary/Species_richness_taxonomic_group.pdf",6,6,bg="transparent")
endemicity <- ggplot(data.taxa.modified,aes(x=(Endemic.species/No..of.species),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(name = "Taxonomic\ngroups", values=c("#00c7e7","#898989","#77dd77"))+
  scale_color_manual(name = "Taxonomic\ngroups",values=c("#00c7e7","#898989","#77dd77"))+
  xlab("Endemicity")+
  geom_density(alpha=0.25)+theme_classic()+theme(legend.position = "inside",
                                                 legend.position.inside = c(0.3,0.75),
                                                 legend.key.size = unit(0.2, 'cm'),
                                                 legend.text = element_text(size = 5,colour = "black"),
                                                 legend.title = element_text(size = 5.5,colour = "black"))
species_richness <- ggplot(data.taxa.modified,aes(x=(No..of.species),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(name = "Taxonomic\ngroups",values=c("#00c7e7","#898989","#77dd77"))+
  scale_color_manual(name = "Taxonomic\ngroups",values=c("#00c7e7","#898989","#77dd77"))+
  xlab("Species Richness")+
  geom_density(alpha=0.25)+theme_classic()+theme(legend.position = "inside",
                                                 legend.position.inside = c(0.5,0.75),
                                                 legend.key.size = unit(0.2, 'cm'),
                                                 legend.text = element_text(size = 5,colour = "black"),
                                                 legend.title = element_text(size = 5.5,colour = "black"),
                                                 axis.title.y=element_blank())
clade_age <- ggplot(data.taxa.modified,aes(x=(Crown.age..Mya.),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(name = "Taxonomic\ngroups",values=c("#00c7e7","#898989","#77dd77"))+
  scale_color_manual(name = "Taxonomic\ngroups",values=c("#00c7e7","#898989","#77dd77"))+
  xlab("Clade age (Mya)")+
  geom_density(alpha=0.25)+theme_classic()+theme(legend.position = "inside",
                                                 legend.position.inside = c(0.75,0.75),
                                                 legend.key.size = unit(0.2, 'cm'),
                                                 legend.text = element_text(size = 5,colour = "black"),
                                                 legend.title = element_text(size = 5.5,colour = "black"),
                                                 axis.title.y=element_blank())
grid.arrange(endemicity,species_richness,clade_age,ncol=3,nrow=2)
dev.off()


# Speciation and Extinction rates for three methods - RPANDA, CoMET, ClaDS
# function to add sample size with categorical predictors
give.n <- function(x){
  data.taxa.modified %>% group_by(!!x) %>% summarise(n=n())
}


table_taxonomic_group <- give.n(data.taxa.modified$Taxonomic.group)


data.taxa.modified2 <-  as_tibble(data.taxa.modified) %>% 
  mutate(Taxonomic.group2 = paste0(Taxonomic.group,"(",
                                   table_taxonomic_group$n[match(data.taxa.modified$Taxonomic.group,table_taxonomic_group$`<fct>`)]
                                   ,")"))
data.taxa.modified2 <- as.data.frame(data.taxa.modified2)

data.taxa.modified2$Taxonomic.group2 <- factor(data.taxa.modified2$Taxonomic.group2,levels = paste0(as.character(table_taxonomic_group$`<fct>`),"(",table_taxonomic_group$n,")"))

n.Angiosperm <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Angiosperm")]
n.Herpetofauna <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Herpetofauna")]
n.Invertebrate <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Invertebrate")]



# habitat
table_habitat_type <- give.n(data.taxa.modified$Habitat)

data.habi.modified <-  as_tibble(data.taxa.modified) %>% 
  mutate(habi = paste0(Habitat,"(",
                       table_habitat_type$n[match(data.taxa.modified$Habitat,table_habitat_type$`<chr>`)]
                       ,")"))
data.habi.modified <- as.data.frame(data.habi.modified)

data.habi.modified$habi <- factor(data.habi.modified$habi,levels = paste0(table_habitat_type$`<chr>`,"(",table_habitat_type$n,")"))


n.FW <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "FW")]
n.Wet <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "Wet")]
n.WetDry <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "Wet + Dry")]



# Biogeography
table_biogeographic_origin <- give.n(data.taxa.modified$Biogeographic.origin)

data.biogeo.modified <-  as_tibble(data.taxa.modified) %>% 
  mutate(biogeo = paste0(Biogeographic.origin,"(",
                         table_biogeographic_origin$n[match(data.taxa.modified$Biogeographic.origin,table_biogeographic_origin$`<chr>`)]
                         ,")"))
data.biogeo.modified <- as.data.frame(data.biogeo.modified)

n.Asian <- table_biogeographic_origin$n[which(table_biogeographic_origin$`<chr>` == "Asian")]
n.Gondwanan <- table_biogeographic_origin$n[which(table_biogeographic_origin$`<chr>` == "Gondwanan")]



# Speciation rates
mean.rate.rpanda <- mean(data.taxa.modified$Speciation.rate..RPANDA.)
mean.rate.CoMET <- mean(data.taxa.modified$Speciation.rate..CoMET.)
mean.rate.ClaDS <- mean(data.taxa.modified$Speciation.rate..ClaDS.)
tree <- read.tree("./SUPER_TREE.nwk")
tree <- drop.tip(tree,"Montecincla")

# comparing across taxonomic groups
pd1<- ggbetweenstats(
  data = data.taxa.modified2,
  x = Taxonomic.group2,
  y = Speciation.rate..RPANDA.,
  color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),
  point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate (RPANDA)")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())



# comparing across habitat types
pd2 <- ggbetweenstats(
  data = data.habi.modified,
  x = habi,y=Speciation.rate..RPANDA.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Speciation rate RPANDA.")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..RPANDA.,color=biogeo,type = "nonparametric")
pd3.data <- extract_stats(pd3)$subtitle_data
pval3 <- round(pd3.data$p.value,5)

pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..RPANDA.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval3<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval3))))),"")
                  #'paste(italic("p"),"=",pval)'
                  ),
            size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd4 <- mk_pgls_plot(Speciation.rate..RPANDA.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd5<- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Speciation.rate..CoMET.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate (CoMET)")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd6 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Speciation.rate..CoMET.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..CoMET.,color=biogeo,type = "nonparametric")
pd7.data <- extract_stats(pd7)$subtitle_data
pval7 <- round(pd7.data$p.value,5)

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..CoMET.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval6<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval7))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd8 <- mk_pgls_plot(Speciation.rate..CoMET.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd9 <- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Speciation.rate..ClaDS.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate (ClaDS)")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd10 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Speciation.rate..ClaDS.,color=habi, 
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..ClaDS.,color=biogeo,type = "nonparametric")
pd11.data <- extract_stats(pd11)$subtitle_data
pval11 <- round(pd11.data$p.value,5)

pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Speciation.rate..ClaDS.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval9<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval11))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd12 <- mk_pgls_plot(Speciation.rate..ClaDS.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "Clade age (Mya)",ylab="",xlim.min = 0,ylim.min = 0)

CairoPDF("../graphs/supplementary/Supplementary_Speciation_Rates.pdf",8.3,11.7,bg="transparent")
# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,nrow=3,ncol=4)
dev.off()

# Extinction rates
mean.rate.rpanda <- mean(data.taxa.modified$Extinction.rate..RPANDA.)
mean.rate.CoMET <- mean(data.taxa.modified$Extinction.rate..CoMET.)
mean.rate.ClaDS <- mean(data.taxa.modified$Extinction.rate..ClaDS.)

# comparing across taxonomic groups
pd1<- ggbetweenstats(
  data = data.taxa.modified2,
  x = Taxonomic.group2,
  y = Extinction.rate..RPANDA.,
  color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),
  point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate (RPANDA)")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())



# comparing across habitat types
pd2 <- ggbetweenstats(
  data = data.habi.modified,
  x = habi,y=Extinction.rate..RPANDA.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Extinction rate RPANDA.")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..RPANDA.,color=biogeo,type = "nonparametric")
pd3.data <- extract_stats(pd3)$subtitle_data
pval3 <- round(pd3.data$p.value,5)

pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..RPANDA.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.25))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval3<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval3))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd4 <- mk_pgls_plot(Extinction.rate..RPANDA.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd5<- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Extinction.rate..CoMET.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate (CoMET)")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd6 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Extinction.rate..CoMET.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..CoMET.,color=biogeo,type = "nonparametric")
pd7.data <- extract_stats(pd7)$subtitle_data
pval7 <- round(pd7.data$p.value,5)

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..CoMET.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.25))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval6<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval7))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd8 <- mk_pgls_plot(Extinction.rate..CoMET.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd9 <- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Extinction.rate..ClaDS.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate (ClaDS)")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd10 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Extinction.rate..ClaDS.,color=habi, 
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(0,0.25))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..ClaDS.,color=biogeo,type = "nonparametric")
pd11.data <- extract_stats(pd11)$subtitle_data
pval11 <- round(pd11.data$p.value,5)

pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Extinction.rate..ClaDS.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(0,0.25))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval9<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval11))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd12 <- mk_pgls_plot(Extinction.rate..ClaDS.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "Clade age (Mya)",ylab="",xlim.min = 0,ylim.min = 0)

CairoPDF("../graphs/supplementary/Supplementary_Extinction_Rates.pdf",8.3,11.7,bg="transparent")
# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,nrow=3,ncol=4)
dev.off()



# Net-diversifcation rates
mean.rate.rpanda <- mean(data.taxa.modified$Net.diversification.rate..RPANDA.)
mean.rate.CoMET <- mean(data.taxa.modified$Net.diversification.rate..CoMET.)
mean.rate.ClaDS <- mean(data.taxa.modified$Net.diversification.rate..ClaDS.)

# comparing across taxonomic groups
pd1<- ggbetweenstats(
  data = data.taxa.modified2,
  x = Taxonomic.group2,
  y = Net.diversification.rate..RPANDA.,
  color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),
  point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-diversification rate (RPANDA)")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())



# comparing across habitat types
pd2 <- ggbetweenstats(
  data = data.habi.modified,
  x = habi,y=Net.diversification.rate..RPANDA.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-diversification rate RPANDA.")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..RPANDA.,color=biogeo,type = "nonparametric")
pd3.data <- extract_stats(pd3)$subtitle_data
pval3 <- round(pd3.data$p.value,5)

pd3 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..RPANDA.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(-0.03,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval3<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval3))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd4 <- mk_pgls_plot(Net.diversification.rate..RPANDA.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd5<- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Net.diversification.rate..CoMET.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-diversification rate (CoMET)")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd6 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Net.diversification.rate..CoMET.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..CoMET.,color=biogeo,type = "nonparametric")
pd7.data <- extract_stats(pd7)$subtitle_data
pval7 <- round(pd7.data$p.value,5)

pd7 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..CoMET.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(-0.03,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval6<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval7))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd8 <- mk_pgls_plot(Net.diversification.rate..CoMET.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "",ylab="",xlim.min = 0,ylim.min = 0)

pd9 <- ggbetweenstats(
  data = data.taxa.modified2, x=Taxonomic.group2,y=Net.diversification.rate..ClaDS.,color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified2$Taxonomic.group2),lwd=0.3,fill=NA)+scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-diversification rate (ClaDS)")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd10 <- ggbetweenstats(
  data = data.habi.modified, x=habi,y=Net.diversification.rate..ClaDS.,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified$habi),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(-0.03,0.5))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..ClaDS.,color=biogeo,type = "nonparametric")
pd11.data <- extract_stats(pd11)$subtitle_data
pval11 <- round(pd11.data$p.value,5)

pd11 <- ggbetweenstats(
  data = data.biogeo.modified,x=biogeo,y=Net.diversification.rate..ClaDS.,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(-0.03,0.5))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval6<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval11))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd12 <- mk_pgls_plot(Net.diversification.rate..ClaDS.~Crown.age..Mya.,phy = tree,df.pgls = data.taxa.modified,xlab = "Clade age (Mya)",ylab="",xlim.min = 0,ylim.min = 0)


CairoPDF("../graphs/supplementary/Supplementary_Net_diversification_Rates.pdf",8.3,11.7,bg="transparent")
# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,pd10,pd11,pd12,nrow=3,ncol=4)
dev.off()


######## Pulled Diversification Rates #################
data.taxa.modified3 <- data.taxa.modified[-which(is.na(data.taxa.modified$Pulled.Diversification.Rate)),]
data.biogeo.modified2 <- data.taxa.modified[-which(is.na(data.taxa.modified$Pulled.Diversification.Rate)),]
data.habi.modified2 <- data.taxa.modified[-which(is.na(data.taxa.modified$Pulled.Diversification.Rate)),]

give.n <- function(x){
  data.taxa.modified3 %>% group_by(!!x) %>% summarise(n=n())
}


table_taxonomic_group <- give.n(data.taxa.modified3$Taxonomic.group)


data.taxa.modified3 <-  as_tibble(data.taxa.modified3) %>% 
  mutate(Taxonomic.group2 = paste0(Taxonomic.group,"(",
                                   table_taxonomic_group$n[match(data.taxa.modified3$Taxonomic.group,table_taxonomic_group$`<fct>`)]
                                   ,")"))
data.taxa.modified3 <- as.data.frame(data.taxa.modified3)

data.taxa.modified3$Taxonomic.group2 <- factor(data.taxa.modified3$Taxonomic.group2,levels = paste0(as.character(table_taxonomic_group$`<fct>`),"(",table_taxonomic_group$n,")"))

n.Angiosperm <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Angiosperm")]
n.Herpetofauna <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Herpetofauna")]
n.Invertebrate <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Invertebrate")]



# habitat
table_habitat_type <- give.n(data.habi.modified2$Habitat)

data.habi.modified2 <-  as_tibble(data.habi.modified2) %>% 
  mutate(habi = paste0(Habitat,"(",
                       table_habitat_type$n[match(data.habi.modified2$Habitat,table_habitat_type$`<chr>`)]
                       ,")"))
data.habi.modified2 <- as.data.frame(data.habi.modified2)

#data.habi.modified2$habi <- factor(data.habi.modified2$habi,levels = paste0(table_habitat_type$`<chr>`,"(",table_habitat_type$n,")"))


n.FW <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "FW")]
n.Wet <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "Wet")]
n.WetDry <- table_habitat_type$n[which(table_habitat_type$`<chr>` == "Wet + Dry")]



# Biogeography
table_biogeographic_origin <- give.n(data.biogeo.modified2$Biogeographic.origin)

data.biogeo.modified2 <-  as_tibble(data.biogeo.modified2) %>% 
  mutate(biogeo = paste0(Biogeographic.origin,"(",
                         table_biogeographic_origin$n[match(data.biogeo.modified2$Biogeographic.origin,table_biogeographic_origin$`<chr>`)]
                         ,")"))
data.biogeo.modified <- as.data.frame(data.biogeo.modified2)

n.Asian <- table_biogeographic_origin$n[which(table_biogeographic_origin$`<chr>` == "Asian")]
n.Gondwanan <- table_biogeographic_origin$n[which(table_biogeographic_origin$`<chr>` == "Gondwanan")]


# pulled rates
mean.rate.rpanda <- mean(data.taxa.modified3$Pulled.Diversification.Rate)
mean.rate.CoMET <- mean(data.taxa.modified3$Pulled.Diversification.Rate)
mean.rate.ClaDS <- mean(data.taxa.modified3$Pulled.Diversification.Rate)

# change tree after dropping the lineages for which PDR could not be estimated
tree.pdr <- drop.tip(tree,data$Lineage[which(!(data$Lineage %in% data.taxa.modified3$Lineage))])

# comparing across taxonomic groups
pd1<- ggbetweenstats(
  data = data.taxa.modified3,
  x = Taxonomic.group2,
  y = Pulled.Diversification.Rate,
  color=Taxonomic.group2,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),
  point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified3$Taxonomic.group2),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("#00c7e7","#898989","#77dd77"))+
  ylim(c(-0.3,0.7))+
  theme_classic()+xlab("Taxonomic group")+ylab("Pulled diversification rate (RPANDA)")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())



# comparing across habitat types
pd2 <- ggbetweenstats(
  data = data.habi.modified2,
  x = habi,y=Pulled.Diversification.Rate,color=habi,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.habi.modified2$habi),lwd=0.3,fill=NA)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  ylim(c(-0.3,0.7))+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Extinction rate RPANDA.")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

#comparing across lineages with varying biogeographic origins
pd3 <- ggbetweenstats(
  data = data.biogeo.modified2,x=biogeo,y=Pulled.Diversification.Rate,color=biogeo,type = "nonparametric")
pd3.data <- extract_stats(pd3)$subtitle_data
pval3 <- round(pd3.data$p.value,5)

pd3 <- ggbetweenstats(
  data = data.biogeo.modified2,x=biogeo,y=Pulled.Diversification.Rate,color=biogeo,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.biogeo.modified2$biogeo),lwd=0.3,fill=NA)+scale_color_manual(values=c("grey10","grey40"))+
  ylim(c(-0.3,0.7))+
  geom_text(aes(y=0.45,x=2,label=ifelse(pval3<0.05,as.character(as.expression((substitute(italic("p")~"="~p,list(p=pval3))))),"")
                #'paste(italic("p"),"=",pval)'
  ),
  size = 3,check_overlap = TRUE,color="black",parse=T)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pd4 <- mk_pgls_plot(Pulled.Diversification.Rate~Crown.age..Mya.,phy = tree.pdr,df.pgls = data.taxa.modified3,xlab = "",ylab="",xlim.min = 0,ylim.min = -0.03)
CairoPDF("../graphs/supplementary/Supplementary_Pulled_diversification_Rates.pdf",8.3,11.7,bg="transparent")
# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,nrow=3,ncol=2)
dev.off()





######## SUPPLEMENTARY scenarios ############

pb0 <- ggplot(data.taxa.modified, aes(x=Final.Results))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+
  theme(
        axis.text.y=element_text(size=12,angle=0,hjust=1,vjust=0.9),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=14))#+coord_flip()

pb1 <-ggplot(data.taxa.modified, aes(x=Final.Results,fill=Taxonomic.group))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Biogeographic origins")+
  scale_fill_manual(name = "Taxonomic\ngroups",values = c("#00c7e7","#898989","#77dd77"))+
  ggtitle(label = "All lineages")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=14),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()

pb2 <- ggplot(data.taxa.modified, aes(x=Final.Results,fill=Habitat))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Habitat types",values = c("grey10","grey40","grey70"))+
  ggtitle(label = "All lineages")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()

pb3 <- ggplot(data.taxa.modified, aes(x=Final.Results,fill=Biogeographic.origin))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Biogeographic\norigins",values = c("grey10","grey40"))+
  ggtitle(label = "All lineages")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()

############# summary with different lots of lineages with different criteria for clade sizes
data.richness.greater.than.10 <- data.taxa.modified[which(data.taxa.modified$No..of.species >= 10),]
data.richness.greater.than.20 <- data.taxa.modified[which(data.taxa.modified$No..of.species >= 20),]
data.richness.greater.than.30 <- data.taxa.modified[which(data.taxa.modified$No..of.species >= 30),]

pb1.geq10 <-ggplot(data.richness.greater.than.10, aes(x=Final.Results,fill=Taxonomic.group))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Biogeographic origins")+
  scale_fill_manual(name = "Taxonomic\ngroups",values = c("#00c7e7","#898989","#77dd77"))+
  ggtitle(label = "\u2265 10 species")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=14),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb2.geq10 <- ggplot(data.richness.greater.than.10, aes(x=Final.Results,fill=Habitat))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Habitat types",values = c("grey10","grey40","grey70"))+
  ggtitle(label = "\u2265 10 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb3.geq10 <- ggplot(data.richness.greater.than.10, aes(x=Final.Results,fill=Biogeographic.origin))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Biogeographic\norigins",values = c("grey10","grey40"))+
  ggtitle(label = "\u2265 10 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()


pb1.geq20 <-ggplot(data.richness.greater.than.20, aes(x=Final.Results,fill=Taxonomic.group))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Biogeographic origins")+
  scale_fill_manual(name = "Taxonomic\ngroups",values = c("#00c7e7","#898989","#77dd77"))+
  ggtitle(label = "\u2265 20 species")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=14),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb2.geq20 <- ggplot(data.richness.greater.than.20, aes(x=Final.Results,fill=Habitat))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Habitat types",values = c("grey40","grey70"))+
  ggtitle(label = "\u2265 20 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb3.geq20 <- ggplot(data.richness.greater.than.20, aes(x=Final.Results,fill=Biogeographic.origin))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Biogeographic\norigins",values = c("grey10","grey40"))+
  ggtitle(label = "\u2265 20 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()



pb1.geq30 <-ggplot(data.richness.greater.than.30, aes(x=Final.Results,fill=Taxonomic.group))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Biogeographic origins")+
  scale_fill_manual(name = "Taxonomic\ngroups",values = c("#898989","#77dd77"))+
  ggtitle(label = "\u2265 30 species")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=14),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb2.geq30 <- ggplot(data.richness.greater.than.30, aes(x=Final.Results,fill=Habitat))+geom_bar(width=0.5)+theme_classic()+xlab("Diversification scenarios")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Habitat types",values = c("grey40","grey70"))+
  ggtitle(label = "\u2265 30 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()

pb3.geq30 <- ggplot(data.richness.greater.than.30, aes(x=Final.Results,fill=Biogeographic.origin))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Biogeographic\norigins",values = c("grey10","grey40"))+
  ggtitle(label = "\u2265 30 species")+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=6,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=6,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_blank(),
        legend.position = "none",
        legend.position.inside = c(0.3,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 2.5,colour = "black"),
        legend.title = element_text(size = 3,colour = "black"))#+coord_flip()


# SUPPLEMENTARY #
CairoPDF("../graphs/supplementary/supplementary_scenarios.pdf",8.3,11.7,bg="transparent")
grid.arrange(pb1,pb2,pb3,
             pb1.geq10,pb2.geq10,pb3.geq10,
             pb1.geq20,pb2.geq20,pb3.geq20,
             pb1.geq30,pb2.geq30,pb3.geq30,
             nrow=4,ncol=3)
dev.off()



pTDD_SR <- ggbetweenstats(
  data = data.taxa.modified,
  x = TDD,y=No..of.species,color=TDD,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 6, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "all",
  centrality.plotting = FALSE,
  results.subtitle = T,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$TDD),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,100))+
  theme_classic()+xlab("Drivers")+ylab("Species Richness")+
  geom_hline(yintercept = mean(data.taxa.modified$No..of.species), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 7),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())


pDDD_SR <- ggbetweenstats(
  data = data.taxa.modified,
  x = DDD,y=No..of.species,color=DDD,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = T,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$DDD),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,100))+
  theme_classic()+xlab("Drivers")+ylab("Species Richness")+
  geom_hline(yintercept = mean(data.taxa.modified$No..of.species), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 7),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pTDD_DR <- ggbetweenstats(
  data = data.taxa.modified,
  x = TDD,y=Net.diversification.rate..ClaDS.,color=TDD,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = T,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$TDD),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Drivers")+ylab("Net-diversification rate")+
  geom_hline(yintercept = mean(data.taxa.modified$Net.diversification.rate..ClaDS.), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 7),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pDDD_DR <- ggbetweenstats(
  data = data.taxa.modified,
  x = DDD,y=Net.diversification.rate..ClaDS.,color=DDD,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = T,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$DDD),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Drivers")+ylab("Net-diversification rate")+
  geom_hline(yintercept = mean(data.taxa.modified$Net.diversification.rate..ClaDS.), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 7),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())


pDS_SR <- ggbetweenstats(
  data = data.taxa.modified,
  x = Final.Results,y=No..of.species,color=Final.Results,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$Final.Results),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,100))+
  theme_classic()+xlab("Scenarios")+ylab("Species Richness")+
  geom_hline(yintercept = mean(data.taxa.modified$No..of.species), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

pDS_DR <-ggbetweenstats(
  data = data.taxa.modified,
  x = Final.Results,y=Net.diversification.rate..ClaDS.,color=Final.Results,
  boxplot.args = list(width = 0, linewidth = 0,col=NA),
  violin.args = list(width = 0, linewidth = 0,col=NA),
  ggsignif.args = list(textsize = 2.5, tip_length = 0.01),   point.args = list(alpha = 0.6, size = 1.5),
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  p.adjust.method = "bonferroni",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  results.subtitle = F,
  bf.message = FALSE
)+geom_boxplot(aes(color=data.taxa.modified$Final.Results),lwd=0.3,fill=NA)+
  scale_color_jco()+
  ylim(c(0,0.5))+
  theme_classic()+xlab("Scenarios")+ylab("Net-diversification rate")+
  geom_hline(yintercept = mean(data.taxa.modified$Net.diversification.rate..ClaDS.), linetype = 2,size=0.5,color="#000000")+
  theme(legend.position="none",
        text = element_text(size = 2),
        axis.title = element_text(size=7.5,color="black"), axis.text = element_text(size=6,color="black"),
        axis.text.x=element_text(angle=20,hjust=1,vjust=0.9),
        axis.title.x=element_blank(),axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank())

CairoPDF("../../graphs/supplementary/supplementary_scenarios_drivers.pdf",8.3,11.7,bg="transparent")
grid.arrange(#arrangeGrob(
  pTDD_SR,
  #pDDD_SR,ncol=2),
            # arrangeGrob(
               pTDD_DR,
              # pDDD_DR,ncol=2),
            # arrangeGrob(
               pDS_SR,
               pDS_DR,
               #ncol=2),
               ncol=2)
dev.off()
pdf_combine(input=c("../graphs/supplementary/supplementary_FAMD.pdf","../graphs/supplementary/supplementary_scenarios_drivers.pdf"),output = "../../graphs/supplementary/supplementary_FAMD2.pdf")

# species richness and div rate box plots comparing TD, DD and Scenarios #
# supplementary tables #
# LTT and RTT supplementary #

############################################################################################################################
############################################################################################################################


t <- read.tree("./SUPER_TREE.nwk")
t <- drop.tip(t,"Montecincla")
t2 <- read.tree("./SUPER_TREE_with_all_tips_ultrametric.nwk")
t2 <- drop.tip(t2,t2$tip.label[grep("Montecincla",t2$tip.label)])

node.numbers <- vector("numeric",Ntip(t))
indices <- c()
names(node.numbers) <- t$tip.label
for(i in 1:Ntip(t)){
  ind <- grep(t$tip.label[i],t2$tip.label)
  node.numbers[i] <- getMRCA(t2,t2$tip.label[ind])
  indices <- c(indices,ind)
}


p <- ggtree(t,ladderize = F) +
  coord_geo(xlim = c(-1650, 0), ylim = c(-2, Ntip(t)), neg = TRUE, abbrv = F,dat = "era",size = 1.5,height = unit(1, "line")) +
  geom_tiplab(as_ylab = TRUE,aes(face="italic",size=0.5))+
  scale_x_continuous(breaks = seq(-1600, 0, 250), labels = abs(seq(-1600, 0, 250))) +
  theme_classic()+theme(axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank())
p2 <- ggtree(t2,size=0.1,ladderize = F) + theme_tree2()+
  coord_geo(xlim = c(-1650, 480), ylim = c(-2, Ntip(t2)), neg = TRUE, abbrv = F,dat = "era",size = 1.5,height = unit(1, "line")) +
  geom_tiplab(geom = "text", as_ylab = TRUE,face="italic",size=0)+
  scale_x_continuous(breaks = seq(-1600, 0, 250), labels = abs(seq(-1600, 0, 250))) +
  theme_classic()+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.line.y = element_blank(),axis.line.x=element_blank())
p2 <- revts(p2)
for(i in 1:Ntip(t)){
  p2 <- p2 + geom_cladelab(node=node.numbers[i], label=paste0(names(node.numbers)[i],"(",Ntip(extract.clade(t2,node.numbers[i])),")"), barsize=0.3, barcolour = c(rep(c("blue","red"),16),"blue")[i],fontsize=1.5,align=T)
} 
CairoPDF("../graphs/supplementary/supplementary_pgls_Tree_clades.pdf",
         8.30,11.70,bg="transparent")
grid.arrange(revts(p),p2,nrow=2)
dev.off()


############################################################################################################################
# CRABS files processing
dir.create("../graphs/supplementary/CRABS")
pdfnames <- list.files("../data/CRABS_results/")
system(paste0("pdftk CRABS_results/",pdfnames[2]," cat 1 output ../graphs/supplementary/CRABS/1_functions_used.pdf"))
for(i in 2:length(pdfnames)){
  #i <- 2
  lineage.name <- substr(pdfnames[i],1,nchar(pdfnames[i])-25)
  
  pdf(paste0("../graphs/supplementary/CRABS/",lineage.name,"_stamp.pdf"),8.3,11.7,bg="transparent")
  print(ggplot(data=as.data.frame(4))+theme_classic()+scale_y_continuous(position = "right")+ylab(lineage.name)+
  theme(axis.title.y = element_text(size=25,hjust=0),panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA)))
  dev.off()
  
  system(paste0("pdftk CRABS_results/",pdfnames[i]," cat 2 output ../graphs/supplementary/CRABS/",lineage.name,"_1.pdf"))
  system(paste0("pdftk ","../graphs/supplementary/CRABS/",lineage.name,"_1.pdf"," multistamp ","../graphs/supplementary/CRABS/",lineage.name,"_stamp.pdf ","output ","../graphs/supplementary/CRABS/",lineage.name,"_stamp1.pdf"))
  system(paste0("rm ../graphs/supplementary/CRABS/",lineage.name,"_1.pdf "))
  
  system(paste0("pdftk CRABS_results/",pdfnames[i]," cat 9 output ../graphs/supplementary/CRABS/",lineage.name,"_2.pdf"))
  system(paste0("pdftk ","../graphs/supplementary/CRABS/",lineage.name,"_2.pdf"," multistamp ","../graphs/supplementary/CRABS/",lineage.name,"_stamp.pdf ","output ","../graphs/supplementary/CRABS/",lineage.name,"_stamp2.pdf"))
  system(paste0("rm ../graphs/supplementary/CRABS/",lineage.name,"_2.pdf "))
  
  system(paste0("pdftk CRABS_results/",pdfnames[i]," cat 16 output ../graphs/supplementary/CRABS/",lineage.name,"_3.pdf"))
  system(paste0("pdftk ","../graphs/supplementary/CRABS/",lineage.name,"_3.pdf"," multistamp ","../graphs/supplementary/CRABS/",lineage.name,"_stamp.pdf ","output ","../graphs/supplementary/CRABS/",lineage.name,"_stamp3.pdf"))
  system(paste0("rm ../graphs/supplementary/CRABS/",lineage.name,"_3.pdf ","../graphs/supplementary/CRABS/",lineage.name,"_stamp.pdf "))
}


pdf_combine(input=paste0("../graphs/supplementary/CRABS/",list.files("../graphs/supplementary/CRABS/")),
            output = "../graphs/supplementary/CRABS_output_combined.pdf")

########### combining PDFs ##################
setwd("../graphs/supplementary/")
pdf_combine(input=c("Caption3.pdf","3. Supplementary_Speciation_Rates.pdf","4. Supplementary_Extinction_Rates.pdf",
                    "5. Supplementary_Net_diversification_Rates.pdf","6. Supplementary_Pulled_diversification_Rates.pdf",
                    "7. supplementary_pgls_Tree_clades.pdf"),
            output = "3_Diversification_rates.pdf")
pdf_combine(input=c("Caption5.pdf","9. CRABS_output_combined_A4.pdf"),
            output = "5_CRABS.pdf")
pdf_combine(input=c("Caption6.pdf","10. supplementary_FAMD2.pdf"),
            output = "6_FAMD.pdf")
pdf_combine(input=c("Caption7.pdf","11. scenario_estimations_RPANDA_models_a5.pdf","12. scenario_estimationsCoMET_a5.pdf",
                    "13. Model_parameter_estimations.pdf","14. times_of_strong_shift_A5.pdf"),
            output = "7_Credibility_through_Simulations.pdf")

pdf_combine(input = c("SI-Fig_Cover_page.pdf","SI Tables.pdf","1_2.pdf","3_Diversification_rates.pdf","4_Scenarios_vs_everything.pdf",
                      "5_CRABS.pdf","6_FAMD.pdf","7_Credibility_through_Simulations.pdf","SI_References.pdf"
                      ),
            output = "SI_Table_and_Figures.pdf")



