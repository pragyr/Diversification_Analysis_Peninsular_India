library(ggplot2)
library(see)
library(ggpubr)
library(ape)
library(phytools)
library(rr2)
library(nlme)
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


# PLEASE CHECK THE FOLLOWING TWO LINES
# if you are using the "Rscript" command in linux/mac terminal
setwd(paste0(getwd(),"/../data/all_trees"))

# if you are using RStudio IDE
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/../data"))

# OR DIRECTLY PUT THE PATH AS A STRING

data <- read.csv("Final_summary_updated.csv")

head(data)
nrow(data)


data.taxa.modified <- data
data.taxa.modified$Taxonomic.group <- factor(data.taxa.modified$Taxonomic.group,levels=c("Arthropod","Mollusc","Amphibian","Reptile","Bird","Plant"))

geological_events <- c("GS",
                       "DV+SIS",
                       "LC-w-Asia",
                       "HCA+EOC",
                       "Peak H-T-O",
                       "M4G+IA+\nIM+C4-Ex")
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

tax.grp <- c("Arthropod","Mollusc","Amphibian","Reptile","Bird","Plant")
num.tax.grp <- length(tax.grp)
for(i in 1:num.tax.grp){
  #i=1
  data.tax <- data[which(data$Taxonomic.group==tax.grp[i]),]
  data.tax <- data.tax[order(data.tax$Stem.age..Mya.,decreasing = TRUE),]
  data.modified <- rbind(data.modified,data.tax)
}
data.modified$Lineage <- c(1:length(data.modified$Lineage))

data.modified$Lineage <- factor(data.modified$Lineage, levels = data.modified$Lineage)

CairoPDF("./Graphs/1_Data_biogeo_clade_age_habitat.pdf",11.70,8.30,bg="transparent")
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
  scale_fill_manual(values = c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093"))+
  scale_color_manual(values = c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093"))+
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
plg
dev.off()
#######################################################################################################################
##################################################### Boxplots, Barplots ###################################################
##################################################### PACE ######################################################
cm <- list( c("Arthropod", "Mollusc"), c("Arthropod", "Amphibian"), c("Arthropod", "Reptile"),c("Arthropod","Plant"),
            c("Mollusc","Amphibian"),c("Mollusc","Reptile"),c("Mollusc","Plant"),
            c("Amphibian","Reptile"),c("Amphibian","Plant"),
            c("Reptile","Plant"))

give.n <- function(x){
  data.taxa.modified %>% group_by(!!x) %>% summarise(n=n())
}


CairoPDF("./Graphs/2_Rates.pdf",8.3,11.7,bg="transparent")

mean.rate.rpanda <- mean(data.taxa.modified$Net.diversification.rate..RPANDA.)
mean.rate.CoMET <- mean(data.taxa.modified$Net.diversification.rate..CoMET.)
mean.rate.ClaDS <- mean(data.taxa.modified$Net.diversification.rate..ClaDS.)


table_taxonomic_group <- give.n(data.taxa.modified$Taxonomic.group)


data.taxa.modified2 <-  as_tibble(data.taxa.modified) %>% 
  mutate(Taxonomic.group2 = paste0(Taxonomic.group,"(",
                                   table_taxonomic_group$n[match(data.taxa.modified$Taxonomic.group,table_taxonomic_group$`<fct>`)]
                                   ,")"))
data.taxa.modified2 <- as.data.frame(data.taxa.modified2)

data.taxa.modified2$Taxonomic.group2 <- factor(data.taxa.modified2$Taxonomic.group2,levels = paste0(as.character(table_taxonomic_group$`<fct>`),"(",table_taxonomic_group$n,")"))

n.Arthropod <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Arthropod")]
n.Amphibian <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Amphibian")]
n.Mollusc <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Mollusc")]
n.Bird <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Bird")]
n.Reptile <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Reptile")]
n.Plant <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Plant")]


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




# comparing across taxonomic groups
cm1 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")), c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")"))
)
pd1<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Net.diversification.rate..RPANDA.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-Diversification rate\n(RPANDA)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm1,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position="none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))



# comparing across habitat types
cm2 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd2 <- ggplot(data.habi.modified, aes(x=habi,y=Net.diversification.rate..RPANDA.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm2,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position="none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins

cm3 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd3 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Net.diversification.rate..RPANDA.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm3,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position = "none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))




cm4 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")), c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")")))
pd4<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Net.diversification.rate..CoMET.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-Diversification rate\n(CoMET)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm4,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position="none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm5 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd5 <- ggplot(data.habi.modified, aes(x=habi,y=Net.diversification.rate..CoMET.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm5,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position="none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm6 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd6 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Net.diversification.rate..CoMET.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm6,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
        legend.position = "none",axis.text=element_text(color="black"),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))






cm7 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")), c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")")))
pd7 <- data.taxa.modified2 %>% 
           ggplot(., aes(x=Taxonomic.group2,y=Net.diversification.rate..ClaDS.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Net-Diversification rate\n(ClaDS)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm7,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))




cm8 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd8 <- ggplot(data.habi.modified, aes(x=habi,y=Net.diversification.rate..ClaDS.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm8,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins
cm9 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd9 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Net.diversification.rate..ClaDS.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Net.diversification.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.45,6),color="red")+
  stat_compare_means(comparisons = cm9,hide.ns = FALSE,label="p.signif",size=7,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))


# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,nrow=3,ncol=3)
dev.off()








############################################################################################################################
############################################## RATE SHIFTS SUMMARY #########################################################
data.modified <- c()

tax.grp <- c("Arthropod","Mollusc","Amphibian","Reptile","Bird","Plant")
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

#data.modified$Driver[which(data.modified$Driver!="")] <- paste0("/",data.modified$Driver[which(data.modified$Driver!="")])
#data.modified$Driver[which(nchar(data.modified$Driver)>4)] <- substr(data.modified$Driver[which(nchar(data.modified$Driver)>4)],1,8)

# CairoPDF("./Graphs/3b_scenarios.pdf",8.3,5.9,bg="transparent")
# ggplot(data.taxa.modified, aes(x=Results.without.shifts))+geom_bar(width=0.5)+xlab("")+ylab("No. of lineages")+
#   theme_classic(base_size = 40)
# dev.off()

CairoPDF("./Graphs/3_Rate_shifts.pdf",11.7,8.3,bg="transparent")
rt_sft_pl <- ggplot(data.modified, aes(x=Lineage,y=Crown.age..Mya.,fill=Taxonomic.group))+
  scale_x_discrete(position = "top",limits = as.factor(c(1:40)))+#,breaks=as.factor(c(5,8,16,28,29,34)))+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[1],ymax=geological_events_ends[1]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[1],label=geological_events[1]),
            size = 3.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[2],ymax=geological_events_ends[2]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[2],label=geological_events[2]),
            size = 3.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[3],ymax=geological_events_ends[3]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[3]),label=geological_events[3],
            size = 3.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[4],ymax=geological_events_ends[4]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[4],label=geological_events[4]),
            size = 3.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[5],ymax=geological_events_ends[5]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[5],label=geological_events[5]),
            size = 3.5,angle=90,check_overlap = TRUE)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[6],ymax=geological_events_ends[6]),
            alpha = 0.06,
            fill = "#cce7e8")+
  geom_text(aes(x = as.factor(37),y=geological_events_dates[6],label=geological_events[6]),
            size = 3.5,angle=90,check_overlap = TRUE)+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[7],ymax=geological_events_ends[7]),
  #           alpha = 0.06,
  #           fill = "#cce7e8",color="gray40",size=0.1)+
  # geom_text(aes(x = as.factor(37),y=geological_events_dates[7],label=geological_events[7]),
  #           size = 2.5,angle=90,check_overlap = TRUE)+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[8],ymax=geological_events_ends[8]),
  #           alpha = 0.06,
  #           fill = "#cce7e8",color="gray40",size=0.1)+
  # geom_text(aes(x = as.factor(37),y=geological_events_dates[8],label=geological_events[8]),
  #           size = 2.5,angle=90,check_overlap = TRUE)+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[9],ymax=geological_events_ends[9]),
  #           alpha = 0.06,
  #           fill = "#cce7e8",color="gray40",size=0.1)+
  # geom_text(aes(x = as.factor(37),y=geological_events_dates[9],label=geological_events[9]),
  #           size = 2.5,angle=90,check_overlap = TRUE)+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[10],ymax=geological_events_ends[10]),
  #           alpha = 0.06,
  #           fill = "#cce7e8",color="gray40",size=0.1)+
  # geom_text(aes(x = as.factor(37),y=geological_events_dates[10],label=geological_events[10]),
  #           size = 2.5,angle=90,check_overlap = TRUE)+
  # geom_rect(aes(xmin = -Inf, xmax = Inf, ymin=geological_events_starts[11],ymax=geological_events_ends[11]),
  #           alpha = 0.06,
  #           fill = "#cce7e8",color="gray40",size=0.1)+
  # geom_text(aes(x = as.factor(37),y=geological_events_dates[11],label=geological_events[11]),
  #           size = 2.5,angle=90,check_overlap = TRUE)+
  scale_fill_manual(values = c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093"))+
  scale_color_manual(values = c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093"))+
  geom_bar(stat="identity",width=0.2,alpha=0.95)+
  #geom_linerange(aes(ymin=Crown.age..Mya.,ymax=Stem.age..Mya.,colour=Taxonomic.group))+
  geom_point(data=data.modified, aes(y=Up_Peak_or_Down_Dip_Time),fill="blue",shape=24,color="white",size=3.5,alpha=0.8)+
  geom_point(data=data.modified, aes(y=Peak_or_Dip_Time),fill="blue",shape=21,color="white",size=3.5,alpha=1)+
  geom_point(data=data.modified, aes(y=Down_Peak_or_Up_Dip_Time),fill="blue",shape=25,color="white",size=3.5,alpha=0.8)+
  theme_classic()+
  geom_hline(yintercept = 0,size=0.35)+
  scale_y_reverse(limits=c(150,-35),breaks=c(150,125,100,75,50,25,0))+theme(axis.text.x=element_text(size=14),
                          legend.position = "none",
                          plot.title = element_blank(),
                          axis.text = element_text(color="black"),
                          #plot.margin = margin(0,0,0,0, "cm"),
                          axis.title.y=element_blank(),
                          axis.ticks.y=element_blank(),
                          axis.line.y = element_blank(),
                          axis.line.x = element_blank(),
                          axis.text.y=element_blank(),
                          axis.title.x=element_text(size=16))+xlab("Lineages")+ylab("Age (Mya)")+ggtitle("Clade ages")+
  geom_text(aes(x = as.factor(37),y=-12.5,label="Scenario"),
            size = 4.5,fontface="bold",check_overlap = TRUE,color="#0000ff")+
  geom_text(aes(x = as.factor(37),y=-32.5,label="Driver"),
            size = 4.5,fontface="bold",check_overlap = TRUE,color="#0000ff")+
  geom_text(aes(x=Lineage,y=-12.5,
                #(Stem.age..Mya.+8),  
                label=Final.Results),size = 3.0,color="#0000ff")+
  geom_text(aes(x=Lineage,y=-32.5,
                #(Stem.age..Mya.+8),  
                label=Driver),size = 3.0,color="#0000ff")+
  coord_flip()

rt_sft_pl <- gggeo_scale(rt_sft_pl,abbrv = F,size = 5)
# rt_sft_pl
# rt_sft_pl + annotation_custom(ggplotGrob(pb1), xmin = 1, xmax = 3, 
#                        ymin = -0.3, ymax = 0.6)
data.modified$Driver[which(data.modified$Driver=="")] <- "Other"

pl_scenarios_drivers <- ggplot(data.modified, aes(x=Final.Results,fill=Driver))+
  scale_fill_manual(values = c("darkorange","black","darkgreen","gray"))+
  geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+
  theme(axis.text = element_text(colour = "black"),axis.text.y=element_text(size=8,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=8,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=9),legend.text=element_text(size=6),
        legend.position = c(0.64,0.75),legend.key.size = unit(4, 'mm'),legend.title=element_text(size=7),
        plot.background = element_rect(colour = "black", fill=NA, size=0.1))
ggdraw()+draw_plot(rt_sft_pl)+draw_plot(pl_scenarios_drivers,height=1,x=-0.3,y=0.2,scale=0.3)
dev.off()
############################################################################################################################


########################################## RTT(TESS) and LTT Plots ################################################
data("InfTemp")

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
  cat("Started\n\nGetting geological information for plotting.....\n\n")
  
  cols <- c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093")
  names(cols) <- c("Amphibian","Arthropod","Bird","Mollusc","Plant","Reptile")
  
  cols <- cols[which(names(cols) %in% unique(tax.grp))]
  
  long.format.data <- matrix(ncol=4)
  colnames(long.format.data) <- c("Lineage ID","Net-Diversification Rate (events/Myr/lineage)","Time (Mya)","Taxonomic group")
  for(i in 1:length(filenames)){
    #i <- 6
    rates <- as.vector(unlist(read.csv(filenames[i])[1,-1]))
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
  
  #usinf ggplot
  y.lab <- ifelse(ylb,"Net-Diversification Rate\n(events/Myr/lineage)","")
  
  pl1 <- ggplot(long.format.data, aes(x = as.numeric(`Time (Mya)`), y = as.numeric(`Net-Diversification Rate (events/Myr/lineage)`), 
                               color = `Taxonomic group`, group = `Lineage ID`)) + 
    ggtitle(plot_title)+xlab('Time (Mya)')+ylab(y.lab)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[1],xmax=geological_events_ends[1]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[1],label=geological_events[1]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[2],xmax=geological_events_ends[2]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[2],label=geological_events[2]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[3],xmax=geological_events_ends[3]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[3]),label=geological_events[3],
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[4],xmax=geological_events_ends[4]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[4],label=geological_events[4]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[5],xmax=geological_events_ends[5]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[5],label=geological_events[5]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[6],xmax=geological_events_ends[6]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = 0.48,x=geological_events_dates[6],label=geological_events[6]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_text(aes(y = 0.3,x=-5,label="PRESENT DAY"),
              size = 2.5,angle=90,check_overlap = TRUE,color="black")+
    scale_fill_manual(values = cols)+
    scale_color_manual("Taxonomic group",values = cols)+
    geom_line(size=0.1) +
    theme_classic()+
    geom_vline(xintercept = 0,size=0.25)+
    guides(color = guide_legend( 
      override.aes=list(shape = 20)))+
    geom_segment(aes(x = 0, xend = Inf, y = Inf), 
                 linetype = "solid", color = "black", size = 0.1)+
    ylim(c(min(as.numeric(long.format.data$`Net-Diversification Rate (events/Myr/lineage)`)),0.6))+
    xlim(c(0,150))+
    scale_x_reverse(breaks=c(150,125,100,75,50,25,0),limits=c(150,0))+theme(legend.position = "none",
                                                                            axis.line.x = element_blank(),
                                                                            axis.text = element_text(color="black",size=6),
                                                                            plot.title = element_text(hjust = 0.5))
  
  gggeo_scale(pl1,abbrv = F,size = 2.5,height = unit(1, "line"))
  
  
}

make_beautiful_LTT <- function(filenames,tax.grp,text.title,ylb=T){
  # filenames <- c(filenames.LTT.SC2,filenames.LTT.SC2_Shift,filenames.LTT.SC3_Shift,filenames.LTT.Shift)
  # tax.grp <- c(tax.grp.SC2,tax.grp.SC2_Shift,tax.grp.SC3_Shift,tax.grp.Shift)
  # text.title <- "Time-Varying Rates"
  # ylb = F
  cat("Started\n\nGetting geological information for plotting.....\n\n")
  
  cols <- c("#243763","#bba40d","#0000ff","#0dc2c7","#3C6255","#db7093")
  names(cols) <- c("Amphibian","Arthropod","Bird","Mollusc","Plant","Reptile")
  
  cols <- cols[which(names(cols) %in% unique(tax.grp))]
  
  long.format.data <- matrix(ncol=4)
  colnames(long.format.data) <- c("Lineage ID","time","N","Taxonomic group")
  
  for(i in 1:length(filenames)){
    #i <- 9
    tree <- read.tree(filenames[i])
    species_and_time <- ltt.plot.coords(tree)
    if(length(which(species_and_time[,2] == 1))>1){
      species_and_time <- species_and_time[-1,]
    }
    species_and_time[,1] <- abs(species_and_time[,1])
    species_and_time[,2] <- species_and_time[,2]+1
    species_and_time <- rbind(c("time"=max(species_and_time[,1])+0.000000001,"N"=1),species_and_time)
    species_and_time[nrow(species_and_time),2] <- Ntip(tree)
    species_and_time[,2] <- log(species_and_time[,2])
    long.format.data <- rbind(long.format.data,
                              cbind("Lineage ID"=i,
                                    species_and_time, 
                                    "Taxonomic group"=tax.grp[i])
    )
  }
  
  long.format.data <- as_tibble(long.format.data[-1,])
  
  
  if(text.title == "Extinction rates"){
    plot_title <- text.title
  }else{
    plot_title <- paste0("LTT for Lineages with ",text.title)
  }
  
  anti_log_labels <- function(x) {
    exp(x)
  }
  
  
  #usinf ggplot
  y.lab <- ifelse(ylb,"log(No. of lineages)\n","")
  
  pl1 <- ggplot(long.format.data, aes(x = as.numeric(time), y = as.numeric(N), 
                                      color = `Taxonomic group`, group = `Lineage ID`)) + 
    ggtitle(plot_title)+xlab('Time (Mya)')+ylab(y.lab)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[1],xmax=geological_events_ends[1]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[1],label=geological_events[1]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[2],xmax=geological_events_ends[2]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[2],label=geological_events[2]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[3],xmax=geological_events_ends[3]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[3]),label=geological_events[3],
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[4],xmax=geological_events_ends[4]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[4],label=geological_events[4]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[5],xmax=geological_events_ends[5]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[5],label=geological_events[5]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[6],xmax=geological_events_ends[6]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = log(120),x=geological_events_dates[6],label=geological_events[6]),
              size = 1.5,angle=90,check_overlap = TRUE,color="black")+
    geom_text(aes(y = log(10),x=-5,label="PRESENT DAY"),
              size = 2.5,angle=90,color="black",check_overlap = TRUE)+
    scale_fill_manual(values = cols)+
    scale_color_manual(values = cols)+
    scale_y_continuous(labels = anti_log_labels,
                       breaks=log(c(1,2,5,10,20,50,100)),
                       limits = log(c(1,200)))+
    geom_step(size=0.1) +
    theme_classic()+
    geom_vline(xintercept = 0,size=0.25)+geom_segment(aes(x = 0, xend = Inf, y = 0), 
                                                      linetype = "solid", color = "black", size = 0.1)+
    geom_segment(aes(x = 0, xend = Inf, y = Inf), 
                 linetype = "solid", color = "black", size = 0.1)+
    scale_x_reverse(breaks=c(150,125,100,75,50,25,0),limits=c(150,0))+theme(legend.position = "none",
                                                                            axis.line.x = element_blank(),
                                                                            axis.text = element_text(color="black",size=6),
                                                                            plot.title = element_text(hjust = 0.5))
  
  gggeo_scale(pl1,abbrv = F, size=2.5,height = unit(1, "line"))
}

# filenames for RTT Plots 
filenames.SC1 <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC1")],"_nd_rates.csv")
filenames.SC2 <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2")],"_nd_rates.csv")
filenames.SC2_Shift <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2+Shift")],"_nd_rates.csv")
filenames.SC3_Shift <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC3+Shift")],"_nd_rates.csv")
filenames.Shift <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="Shift")],"_nd_rates.csv")
filenames.TDD <- paste0("../estimated_rates_through_time/TESS_RTT_data/",data.taxa.modified$Lineage[which(data.taxa.modified$TDD=="TDD")],"_nd_rates.csv")


# filenames for LTT Plots - Scenaro-wise 
filenames.LTT.SC1 <- paste0("../all_trees/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC1")],".tre")
filenames.LTT.SC2 <- paste0("../all_trees/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2")],".tre")
filenames.LTT.SC2_Shift <- paste0("../all_trees/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC2+Shift")],".tre")
filenames.LTT.SC3_Shift <- paste0("../all_trees/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="SC3+Shift")],".tre")
filenames.LTT.Shift <- paste0("../all_trees/",data.taxa.modified$Lineage[which(data.taxa.modified$Final.Results=="Shift")],".tre")


plot_env_var_with_geoscale <- function(env_var_age,ylab,unsmoothened_data,plot_title){
  # env_var_age <- Temp_dat
  # unsmoothened_data <- InfTemp
  # ylab <- "Temperature (\u00B0C)"
  # plot_title <- "Global Paleo-temperature"
  
  colnames(unsmoothened_data) <- colnames(env_var_age)
  
  long.format.data <- as_tibble(rbind(cbind("Line_ID"=rep(1,nrow(unsmoothened_data)),
                                            unsmoothened_data,
                                            "color_var" = rep("A",nrow(unsmoothened_data)),
                                            "width_var" = rep(1,nrow(unsmoothened_data)),
                                            "style_var" = rep("dotted",nrow(unsmoothened_data))
                            ),
                            cbind("Line_ID"=rep(2,nrow(env_var_age)),
                            env_var_age,
                            "color_var" = rep("B",nrow(env_var_age)),
                            "width_var" = rep(2,nrow(env_var_age)),
                            "style_var" = rep("solid",nrow(env_var_age))
                            )
                            
  ))
  
  cat("Started\n\nGetting geological information for plotting.....\n\n")
  #CairoPDF("./Graphs/Temp_trial.pdf")
  pl1 <- ggplot(long.format.data, aes(x = Age, y = Temperature, group = Line_ID)) + 
    ggtitle(paste(plot_title,"through time"))+xlab('Time (Mya)')+ylab(paste0(ylab,"\n"))+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[1],xmax=geological_events_ends[1]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = Inf,hjust=2,x=geological_events_dates[1],label=geological_events[1]),
              size = 2.5,angle=90,check_overlap = TRUE)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[2],xmax=geological_events_ends[2]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = max(long.format.data$Temperature)-2,x=geological_events_dates[2],label=geological_events[2]),
              size = 1.5,angle=90,check_overlap = TRUE)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[3],xmax=geological_events_ends[3]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = max(long.format.data$Temperature)-2,x=geological_events_dates[3]),label=geological_events[3],
              size = 1.5,angle=90,check_overlap = TRUE)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[4],xmax=geological_events_ends[4]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = max(long.format.data$Temperature)-2,x=geological_events_dates[4],label=geological_events[4]),
              size = 1.5,angle=90,check_overlap = TRUE)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[5],xmax=geological_events_ends[5]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = max(long.format.data$Temperature)-2,x=geological_events_dates[5],label=geological_events[5]),
              size = 1.5,angle=90,check_overlap = TRUE)+
    geom_rect(aes(ymin = -Inf, ymax = Inf, xmin=geological_events_starts[6],xmax=geological_events_ends[6]),
              alpha = 0.06,
              fill = "#cce7e8",color=NA)+
    geom_text(aes(y = max(long.format.data$Temperature)-2,x=geological_events_dates[6],label=geological_events[6]),
              size = 1.5,angle=90,check_overlap = TRUE)+
    geom_text(aes(y = max(long.format.data$Temperature)/2,x=-2,label="PRESENT DAY"),
              size = 2.5,angle=90,check_overlap = TRUE)+
    geom_line(aes(color = color_var, size = as.factor(width_var), linetype = style_var)) +
    scale_color_manual(values=c("#A7A6A6","black"))+
    scale_size_manual(values = c(0.1, 0.0005))+
    scale_linetype_manual(values=c("solid","solid"))+
    theme_classic()+
    ylim(c(min(long.format.data$Temperature),25))+
    geom_segment(aes(x = 0, xend = Inf, y = Inf), 
                 linetype = "solid", color = "black", size = 0.1)+
    geom_vline(xintercept = 0,size=0.25)+scale_x_reverse(breaks=c(150,125,100,75,50,25,0),limits=c(150,0))+theme(legend.position = "none",
                                                                                                    axis.line.x = element_blank(),
                                                                                                    axis.text = element_text(color="black",size=6),
                                                                                                    plot.title = element_text(hjust = 0.5))
  
  gggeo_scale(pl1,dat="periods",abbrv=F,size=2.5,height = unit(1, "line"))
  #dev.off()
}

# get paleo-climate data - temperature (here), as we see high prevalence of temperature-dependence
Temp_dat <- read.csv("../paleoclimate_data/Smooth_Inftemp.csv")
Temp_dat <- Temp_dat[,c("Age","Temperature")]


CairoPDF("./Graphs/4_Rates_and_lineages_through_time.pdf",11.7,8.3,bg="transparent")
p2 <- plot_env_var_with_geoscale(env_var_age = Temp_dat,ylab = "Temperature (\u00B0C)",unsmoothened_data = InfTemp,plot_title = "Paleo-temperature")

p1 <- data.taxa.modified %>% 
  mutate(modified_div_results = ifelse(Final.Results=="SC1",
                                       "Constant",
                                       "Time\nVarying"))%>%
  ggplot(aes(x=modified_div_results,y=Crown.age..Mya.))+
  ggtitle("Diversification trends across clades of varying ages")+
  geom_violinhalf(fill="#D2D2D2",color=NA)+#geom_boxplot(fill='#FFFFFF')+
  geom_boxplot(fill="#D2D2D2",varwidth = T,outliers = FALSE,width=0.3,lwd=0.1) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               dotsize = 0.25, aes(fill=Taxonomic.group,color=Taxonomic.group)) +
  ylab("Clade age (Mya)")+xlab("Diversification Trends")+
  scale_y_reverse(breaks=c(150,125,100,75,50,25,0),limits=c(150,0))+
  scale_color_manual(values = c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  scale_fill_manual(values = c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_segment(aes(y = 0, yend = Inf, x = Inf), 
               linetype = "solid", color = "black", size = 0.1)+
  geom_segment(aes(x = 0, xend = Inf, y = 0), 
               linetype = "solid", color = "black", size = 0.1)+
  geom_label(aes(label=paste("n =",table(data.taxa.modified$Final.Results)['SC1']),
                x=1,y=0),size=3,label.size=NA,nudge_x = 0.4,nudge_y = -10,fill="white")+
  geom_label(aes(label=paste("n =",nrow(data.taxa.modified) - table(data.taxa.modified$Final.Results)['SC1']),
                x=2,y=0),size=3,label.size=NA,nudge_x = 0.4,nudge_y = -10,fill="white")+
  coord_flip()+theme_classic()+
  theme(legend.position = "inside",legend.position.inside = c(0.1,0.55),legend.key.size = unit(0.1, 'cm'),
                                     legend.text = element_text(size = 5,colour = "black"),
                                     legend.title = element_text(size = 6,colour = "black"),
                                     axis.line.x = element_blank(),
                                     axis.text = element_text(color="black",size=6),
                                     plot.title = element_text(hjust = 0.5),
                                     plot.margin = unit(c(0.065, 0.08, 0.065, 0.00), #top,right,bottom,left
                                                        "inches")
          )
p4 <- make_beautiful_RTT(c(filenames.SC2,filenames.SC2_Shift,filenames.SC3_Shift,filenames.Shift),
                   c(tax.grp.SC2,tax.grp.SC2_Shift,tax.grp.SC3_Shift,tax.grp.Shift),
                   "Time-Varying Rates",ylb = F)

p3 <- make_beautiful_RTT(filenames.SC1,tax.grp.SC1,"Constant Rates")

p6 <- make_beautiful_LTT(c(filenames.LTT.SC2,filenames.LTT.SC2_Shift,filenames.LTT.SC3_Shift,filenames.LTT.Shift),
                   c(tax.grp.SC2,tax.grp.SC2_Shift,tax.grp.SC3_Shift,tax.grp.Shift),
                   "Time-Varying Rates",ylb=F)

p5 <- make_beautiful_LTT(filenames.LTT.SC1,tax.grp.SC1,"Constant Rates")
grid.arrange(gggeo_scale(p1,abbrv = F,size = 2.5,height = unit(1, "line")),p2,p3,p4,p5,p6,nrow=3,ncol=2)
dev.off()
#########################################################################################################################################



##########################################################################################################################
########################################## Multivariate-analysis using FAMD ##################################################
############################################SCENARIOS-DRIVERS, DISTRIBUTION ##########################################
# input a matrix with - Lineage names, taxonomic affiliations, clade ages, species richness, net diversification rate,
# biogeographic origin, habitat, if clades are habitat-structured, diversification patterns, 
# if it shows Diversity-dependent diversification, if it shows Temperature-dependent diversification
# but make sure the lineage names column comes first in the matrix
# df <- dt
# list.variables <- c("Biogeographic\norigins")
# palette <- c("grey10","grey40","grey70")
mk_FAMD_plot <- function(df,list.variables,palette){
  # df<-dt
  # list.variables <- NULL
  # palette <- c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255")
  lin.names <- df$Lineage
  df <- df[,-1]
  rownames(df) <- lin.names
  
  res.famd <- FAMD(df,graph = FALSE)
  famd_coordinates <- res.famd$ind$coord
  
  # Calculate dissimilarity matrix. we are using mixed data and hence we are using Gower distance
  diss_matrix <- vegdist(famd_coordinates, method = "gower")
  
  CairoPDF("../Summary and Graphs/Graphs/supplementary_FAMD.pdf",8.3,11.7,bg="transparent")
  p_scree <- fviz_screeplot(res.famd)
  
  # contribution to the 1st Dimension
  pc1 <- fviz_contrib(res.famd,"var",axes=1)
  # contribution to the 2nd Dimension
  pc2 <- fviz_contrib(res.famd,"var",axes=2)
  # visualize quantitative variables
  plot.quanti <- fviz_famd_var(res.famd,"quanti.var",repel=TRUE,col.var = "contrib",gradient.cols=c("red","white","blue"))
  
  groups.variables <- c(rep("A.Taxonomic group",6),
                        rep("B. Biogeographic origin",2),
                        rep("C. Habitat Type",3),
                        rep("D. Habitat structuring",2),
                        rep("E. Scenarios",5),
                        rep("F. Drivers",6))
  
  # visualize qualitative variables
  plot.quali <- fviz_famd_var(res.famd,geom=c("arrow","text"),"quali.var",col.var = factor(groups.variables),
                              palette=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"),repel = TRUE)
  
  grid.arrange(p_scree,pc1,pc2,nrow=3)
  grid.arrange(plot.quanti,plot.quali,nrow=2)
  dev.off()
  
  if(is.null(list.variables)){
    # plot of variables
    fviz_famd_var(res.famd,geom=c("arrow","text"),repel=TRUE)+labs(title = NULL)+theme(axis.text = element_text(color="black"))
  }else{
    perm.var <- permanova_pairwise(diss_matrix,df[,list.variables[1]],method = "gower",padj = "holm")
    name <- unlist(strsplit(list.variables[1],"\n"))
    
    write.csv(perm.var,paste0("./Permanova_",paste0(name[1],"_",name[2]),".csv"))
    # graph of lineages/individual points
    fviz_ellipses(res.famd,list.variables,
                  ellipse.level = 0.95, palette = palette,ellipse.type = c("convex"),repel=TRUE,geom="point",labelsize = 1)+
      theme(axis.text = element_text(color="black"))
  }
  
}

dt <- data.taxa.modified[,c("Lineage","Taxonomic.group","Crown.age..Mya.",
                            "No..of.species","Net.diversification.rate..ClaDS.",
                            "Biogeographic.origin","Habitat","Habitat.differentiation.among.clades.",
                            "Final.Results","DDD","TDD","Him.oro"
                            # ,"C4.Exp"
)]
#df <- df[-which(df$Taxonomic.group=="Bird" | df$Taxonomic.group=="Plant" | df$Taxonomic.group=="Mollusc"),]

colnames(dt) <- c("Lineage","Taxonomic groups","Clade Age",
                  "Species Richness","Net-diversification Rate",
                  "Biogeographic Origins","Habitat Types","Habitat Structuring",
                  "Scenarios","Diversity Dependence","Temperature Dependence","Himalayan-Orogeny Dependence"
                  # ,"C4-Plants-Expansion\nDependence"
)

# # Primary plot defining the diversification space in terms of all variables and how variables are correlated
# CairoPDF("./Graphs/FAMD_drivers.pdf",8.3,5.9,bg="transparent")
# mk_FAMD_plot(dt,list.variables = NULL,c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))
# dev.off()

# Plot to show how patterns and other variables are distributed and structured across the diversification space
CairoPDF("./Graphs/5_scenarios_drivers.pdf",8.3,11.7,bg="transparent")
## Split quantitative and qualitative variables
pl_all <- mk_FAMD_plot(dt,list.variables = NULL,c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))
pl_FAMD_SC <- mk_FAMD_plot(dt,c("Scenarios"),"jco")+guides(color=guide_legend(ncol=3,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
pl_FAMD_Tax <- mk_FAMD_plot(dt,c("Taxonomic groups"),
                            c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  guides(color=guide_legend(ncol=3,title.position = "top"))+
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
pl_FAMD_HS <- mk_FAMD_plot(dt,c("Habitat Structuring"),c("grey10","grey40"))+
  guides(color=guide_legend(ncol=2,title.position = "top"))+
  theme(legend.position="top",legend.text=element_text(size=6.5),legend.title=element_text(size=8),legend.title.align=0.5,
        plot.margin = unit(c(-0.2, 0.05,0.05,0.05), #top,right,bottom,left
                           "inches"))+labs(title = "")
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
grid.arrange(pl_all,arrangeGrob(pl_FAMD_SC,pl_FAMD_Tax,pl_FAMD_Hab,
                                pl_FAMD_Biog,pl_FAMD_DDD,pl_FAMD_TDD,nrow=2,ncol=3),nrow = 2)
dev.off()
##############################################################



###########################################################################################################################
########################################## Phylogenetic Generalised Least Squares Regression -  #############################################
##################################################### SR vs ND vs CA ################################################
mk_regression_plot_PGLS <- function(df.pgls,tree,x.lab,y.lab,xlim.min,ylim.min,plot.title){
  # df.pgls <- dt.temp
  # tree <- tree.temp
  # x.lab <- "Net-diversification rate"
  # y.lab <- "log(Species Richness)"
  # xlim.min <- min(dt.temp$`Net_Diversification_rate (ClaDS)`)
  # ylim.min <- 0
  # plot.title <- "Gondwanan"
  
  colnames(df.pgls) <- c("Predictor","Response")
  # first check whether pagel's lambda=1 or lambda=0 fits your data better
  fitPagel0 <- gls((Response) ~ Predictor, correlation = corPagel(value = 0, phy = tree,
                                                                                            fixed = TRUE), data = df.pgls) # independence
  fitPagel1 <- gls((Response) ~ Predictor, correlation = corPagel(value = 1, phy = tree,
                                                                                            fixed = TRUE), data = df.pgls) # Brownian motion
  
  compare.fits <- anova(fitPagel0,fitPagel1)$AIC
  # compare the two fits
  corLambda <- corPagel(value=c(0,1)[which(compare.fits == min(compare.fits))],phy=tree)
  
  glsControl(maxIter = 500,msMaxIter = 2000,returnObject = TRUE,optimMethod="Nelder-Mead")
  fit.lambda.r <- tryCatch(expr = {suppressWarnings(gls((Response) ~ Predictor, correlation=corLambda,data=df.pgls))},
           error = function(e){return(NA)},
           warning = function(w){},
           finally = {}
  ) 
  if(is.na(fit.lambda.r)[1]){
    fit.lambda.r <- lm((Response) ~ Predictor, data=df.pgls)
    sum.df.r <- summary(fit.lambda.r)
    lambda.r <- NA
    coeff.r <- as.data.frame(sum.df.r$coefficients)
    p.val.r <- coeff.r$`Pr(>|t|)`[2]
    intercept.r <- coeff.r[1,"Estimate"]
    slope.r <- coeff.r[2,"Estimate"]
  }else{
    fit.lambda.r <- gls((Response) ~ Predictor, correlation=corLambda,data=df.pgls)
    sum.df.r <- summary(fit.lambda.r)
    lambda.r <- as.numeric(sum.df.r$modelStruct)
    coeff.r <- as.data.frame(sum.df.r$tTable)
    p.val.r <- coeff.r$`p-value`[2]
    intercept.r <- coeff.r[1,"Value"]
    slope.r <- coeff.r[2,"Value"]
  }
  
  if(lambda.r < 0){ # this can lead to weird R2 estimates. 
    # matt pennel suggests that if lambda is already so low then one might as well just perform lm
    fit.lambda.r <- lm((Response) ~ Predictor, data=df.pgls)
    sum.df.r <- summary(fit.lambda.r)
    lambda.r <- lambda.r
    coeff.r <- as.data.frame(sum.df.r$coefficients)
    p.val.r <- coeff.r$`Pr(>|t|)`[2]
    intercept.r <- coeff.r[1,"Estimate"]
    slope.r <- coeff.r[2,"Estimate"]
  }
  #r2.r <- as.data.frame(t(R2_pred(fit.lambda.r,phy=tree)))$R2_pred
  r2.r <- R2_pred(fit.lambda.r,phy=tree)
  
  
  print.p <- if(p.val.r >= 0.05){sym="="; p.val.r}else{sym <- "<"; if(p.val.r < 0.05 & p.val.r >= 0.01){0.05}else if(p.val.r < 0.01 & p.val.r >= 0.001){0.01}else if(p.val.r < 0.001){0.001}}
  
  #  10.1093/sysbio/syy060: Rpred2 gives the most direct answer to the question of 
  # how much variance in the data is explained by a model. Rresid2 is most appropriate for comparing models fit to different 
  # datasets, because it does not depend on sample sizes. And Rlik2 is most appropriate to assess the importance of different 
  # components within the same model applied to the same data, because it is most closely associated with statistical significance tests.
  stats <- as.character(as.expression(substitute(italic("R")^2~"="~r2*","~italic("p")~s~p,
                                                 list(r2 = round(r2.r,2),
                                                      s = sym,
                                                      p = print.p
                                                      ))))
  if(is.na(lambda.r)){
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = "NA"))))
  }else if(lambda.r < 0){
    lambda.r <- 0
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"<"~l,
                                                    list(l = round(lambda.r,3)))))
  }else{
    lambda <- as.character(as.expression(substitute(italic("\u03BB")~"="~l,
                                                    list(l = round(lambda.r,3)))))
  }
  
  
  # eq <- as.character(as.expression(substitute(italic(y) == a + b %.% italic("* x"), 
  #                                             list(a = format(intercept.r, digits = 2),
  #                                                  b = format(slope.r, digits = 2)))))
  
  if(p.val.r > 0.05){
    lty<-3
  }else{lty<-1}
  
  
  get_labels <- function(x,lab) {
    if(!(-1 %in% unlist(gregexpr("log",lab)))){
      exp(x)
    }else{
      x
    }
    
  }
  get_breaks_custom <- function(vals,lab){
    if(!(-1 %in% unlist(gregexpr("log",lab)))){
      log(c(1,2,5,10,20,50,100,200))
    }else{
      seq(round(min(vals)*1.1,1),round(max(vals)*1.7,1),round(max(vals)*1.1/2,1))
    }
  }
  
  breaks.x <- get_breaks_custom(df.pgls$Predictor,x.lab)
  breaks.y <- get_breaks_custom(df.pgls$Response,y.lab)
    
  ggplot(df.pgls,aes(Predictor,Response)) + 
    geom_point() +
    scale_x_continuous(labels = get_labels(breaks.x,lab = x.lab), breaks = breaks.x,limits = c(xlim.min,(max(df.pgls$Predictor)*1.1)))+
    scale_y_continuous(labels = get_labels(breaks.y,lab = y.lab), breaks = breaks.y,limits = c(ylim.min,(max(df.pgls$Response)*1.5)))+
    geom_abline(intercept = intercept.r, slope = slope.r,linetype=lty) + 
    xlab(ifelse(!(-1 %in% unlist(gregexpr("log",x.lab))),substr(x.lab,(unlist(gregexpr("log",x.lab))+4),(nchar(x.lab)-1)),x.lab))+
    ylab(ifelse(!(-1 %in% unlist(gregexpr("log",y.lab))),substr(y.lab,(unlist(gregexpr("log",y.lab))+4),(nchar(y.lab)-1)),y.lab))+labs(title=plot.title)+
    theme_classic()+
    # geom_text(data = as.data.frame(eq), aes(label = eq,x = max(df.pgls$Predictor)/1.75,
    #                                         y = (max(df.pgls$Response))*1.45,family = 'serif'), parse = TRUE,size=3.5)+
    geom_text(data = as.data.frame(stats), aes(label = stats,x = max(df.pgls$Predictor)/1.3,
                                               y = (max(df.pgls$Response))*1.45,family = 'serif'), parse = TRUE,size=2.5)+
    geom_text(data = as.data.frame(lambda), aes(label = lambda,x = max(df.pgls$Predictor)/1.5,
                                               y = (max(df.pgls$Response))*1.2,family = 'serif'), parse = TRUE,size=2.5)+
    theme(axis.text = element_text(size = 7,colour = "black"))
}


super.tree <- read.tree("/SUPER_TREE.nwk")
#super.tree <- drop.tip(super.tree,data.taxa.modified$Lineage[which(data.taxa.modified$No..of.species < 10)])#c("Montecincla","Parreysia"))

super.tree$tip.label %in% data.taxa.modified$Lineage
data.pgls <- data.taxa.modified[,c("Lineage","Crown.age..Mya.","Net.diversification.rate..RPANDA.","Net.diversification.rate..CoMET.",
                                                                                "Net.diversification.rate..ClaDS.","No..of.species",
             "Taxonomic.group","Biogeographic.origin","Habitat","Final.Results")]
rownames(data.pgls) <- data.pgls$Lineage
data.pgls <- data.pgls[,-which(colnames(data.pgls) == "Lineage")]
data.pgls <- data.pgls[super.tree$tip.label,]

colnames(data.pgls) <- c("Clade_age","Net_Diversification_rate (RPANDA)","Net_Diversification_rate (CoMET)","Net_Diversification_rate (ClaDS)","Species_Richness","Taxonomic.group","Biogeographic.origin","Habitat","Scenarios")
data.pgls$Species_Richness <- log(data.pgls$Species_Richness)
data.pgls$Clade_age <- log(data.pgls$Clade_age)

CairoPDF("./Graphs/6_PGLS.pdf",11.7,2.75,bg="transparent")
## SPECIES RICHNESS VS. CLADE AGE
# main
pgls_main_sr_ca <- mk_regression_plot_PGLS(data.pgls[,c("Clade_age","Species_Richness")],
                                           super.tree,"log(Clade age)","log(Species Richness)",xlim.min = 0,ylim.min = 0,"All lineages")

## SPECIES RICHNESS VS. NET-DIVERSIFICATION RATES
# main
pgls_main_sr_nd_RPANDA <- mk_regression_plot_PGLS(data.pgls[,c("Net_Diversification_rate (RPANDA)","Species_Richness")],
                                           super.tree,"Net-diversification rate (RPANDA)","log(Species Richness)",xlim.min = -0.25,ylim.min = 0,"All lineages")

pgls_main_sr_nd_CoMET <- mk_regression_plot_PGLS(data.pgls[,c("Net_Diversification_rate (CoMET)","Species_Richness")],
                                                 super.tree,"Net-diversification rate (CoMET)","log(Species Richness)",xlim.min = 0,ylim.min = 0,"All lineages")

pgls_main_sr_nd_ClaDS <- mk_regression_plot_PGLS(data.pgls[,c("Net_Diversification_rate (ClaDS)","Species_Richness")],
                                                 super.tree,"Net-diversification rate (ClaDS)","log(Species Richness)",xlim.min = 0,ylim.min = 0,"All lineages")
## NET-DIVERSIFICATION VS. CLADE AGE
# main
pgls_main_nd_ca_RPANDA <- mk_regression_plot_PGLS(data.pgls[,c("Clade_age","Net_Diversification_rate (RPANDA)")],
                                           super.tree,"log(Clade age (Mya))","Net-diversification rate (RPANDA)",xlim.min = 0,ylim.min = -0.25,"All lineages")

pgls_main_nd_ca_CoMET <- mk_regression_plot_PGLS(data.pgls[,c("Clade_age","Net_Diversification_rate (CoMET)")],
                                                  super.tree,"log(Clade age (Mya))","Net-diversification rate (CoMET)",xlim.min = 0,ylim.min = 0,"All lineages")

pgls_main_nd_ca_ClaDS <- mk_regression_plot_PGLS(data.pgls[,c("Clade_age","Net_Diversification_rate (ClaDS)")],
                                                 super.tree,"log(Clade age (Mya))","Net-diversification rate (ClaDS)",xlim.min = 0,ylim.min = 0,"All lineages")

# grid.arrange(arrangeGrob(ggplot()+theme_classic(),
#                          pgls_main_sr_ca,
#                          ggplot()+theme_classic(),
#                          pgls_main_sr_nd_RPANDA,
#                          pgls_main_sr_nd_CoMET,
#                          pgls_main_sr_nd_ClaDS,
#                          pgls_main_nd_ca_RPANDA,
#                          pgls_main_nd_ca_CoMET,
#                          pgls_main_nd_ca_ClaDS,ncol=3,nrow=3))

grid.arrange(pgls_main_sr_ca,pgls_main_sr_nd_ClaDS,pgls_main_nd_ca_ClaDS,ncol=3)
dev.off()









##################### SUPPLEMENTARY FIGURES ################
# supplementary data representation
CairoPDF("./Graphs/Species_richness_taxonomic_group.pdf",6,6,bg="transparent")
endemicity <- ggplot(data.taxa.modified,aes(x=(Endemic.species/No..of.species),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  xlab("Endemicity")+
  geom_density(alpha=0.25)+theme_classic()+theme(legend.position = "inside",
                                                 legend.position.inside = c(0.3,0.75),
                                                 legend.key.size = unit(0.2, 'cm'),
                                                 legend.text = element_text(size = 5,colour = "black"),
                                                 legend.title = element_text(size = 5.5,colour = "black"))
species_richness <- ggplot(data.taxa.modified,aes(x=(No..of.species),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  xlab("Species Richness")+
  geom_density(alpha=0.25)+theme_classic()+theme(legend.position = "inside",
                                                 legend.position.inside = c(0.5,0.75),
                                                 legend.key.size = unit(0.2, 'cm'),
                                                 legend.text = element_text(size = 5,colour = "black"),
                                                 legend.title = element_text(size = 5.5,colour = "black"),
                                                 axis.title.y=element_blank())
clade_age <- ggplot(data.taxa.modified,aes(x=(Crown.age..Mya.),fill=Taxonomic.group,color=Taxonomic.group))+
  scale_fill_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
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
cm <- list( c("Arthropod", "Mollusc"), c("Arthropod", "Amphibian"), c("Arthropod", "Reptile"),c("Arthropod","Plant"),
            c("Mollusc","Amphibian"),c("Mollusc","Reptile"),c("Mollusc","Plant"),
            c("Amphibian","Reptile"),c("Amphibian","Plant"),
            c("Reptile","Plant"))
#with text at 1
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

n.Arthropod <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Arthropod")]
n.Amphibian <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Amphibian")]
n.Mollusc <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Mollusc")]
n.Bird <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Bird")]
n.Reptile <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Reptile")]
n.Plant <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Plant")]


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

# comparing across taxonomic groups

CairoPDF("./Graphs/Supplementary_Speciation_Rates.pdf",8.3,11.7,bg="transparent")
cm1 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd1<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Speciation.rate..RPANDA.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate\n(RPANDA)")+
  #ylim(min(data.taxa.modified$Speciation.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm1,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))



# comparing across habitat types
cm2 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd2 <- ggplot(data.habi.modified, aes(x=habi,y=Speciation.rate..RPANDA.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Speciation rate RPANDA.")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm2,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins

cm3 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd3 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Speciation.rate..RPANDA.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..RPANDA.)-0.001,0.8))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm3,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))




cm4 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
            )
pd4<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Speciation.rate..CoMET.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate\n(CoMET)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm4,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm5 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd5 <- ggplot(data.habi.modified, aes(x=habi,y=Speciation.rate..CoMET.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm5,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm6 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd6 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Speciation.rate..CoMET.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..CoMET.),0.7))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm6,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))


cm7 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd7 <- data.taxa.modified2 %>% 
  ggplot(., aes(x=Taxonomic.group2,y=Speciation.rate..ClaDS.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Speciation rate\n(ClaDS)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm7,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))




cm8 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd8 <- ggplot(data.habi.modified, aes(x=habi,y=Speciation.rate..ClaDS.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm8,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins
cm9 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd9 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Speciation.rate..ClaDS.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Speciation.rate..ClaDS.),0.7))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm9,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))


# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,nrow=3,ncol=3)
dev.off()

# Extinction rates
mean.rate.rpanda <- mean(data.taxa.modified$Extinction.rate..RPANDA.)
mean.rate.CoMET <- mean(data.taxa.modified$Extinction.rate..CoMET.)
mean.rate.ClaDS <- mean(data.taxa.modified$Extinction.rate..ClaDS.)

CairoPDF("./Graphs/Supplementary_Extinction_Rates.pdf",8.3,11.7,bg="transparent")
cm1 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd1<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Extinction.rate..RPANDA.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate\n(RPANDA)")+
  #ylim(min(data.taxa.modified$Extinction.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..RPANDA.)-0.001,1.2))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm1,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))



# comparing across habitat types
cm2 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd2 <- ggplot(data.habi.modified, aes(x=habi,y=Extinction.rate..RPANDA.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Extinction rate RPANDA.")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..RPANDA.)-0.001,1.2))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm2,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins

cm3 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd3 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Extinction.rate..RPANDA.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..RPANDA.)-0.001,1.2))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.rpanda, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.6,6),color="red")+
  stat_compare_means(comparisons = cm3,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))




cm4 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd4<- ggplot(data.taxa.modified2, aes(x=Taxonomic.group2,y=Extinction.rate..CoMET.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate\n(CoMET)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..CoMET.),0.45))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.43,6),color="red")+
  stat_compare_means(comparisons = cm4,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm5 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd5 <- ggplot(data.habi.modified, aes(x=habi,y=Extinction.rate..CoMET.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..CoMET.),0.45))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.43,6),color="red")+
  stat_compare_means(comparisons = cm5,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

cm6 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd6 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Extinction.rate..CoMET.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..CoMET.),0.45))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.CoMET, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.43,6),color="red")+
  stat_compare_means(comparisons = cm6,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))


cm7 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd7 <- data.taxa.modified2 %>% 
  ggplot(., aes(x=Taxonomic.group2,y=Extinction.rate..ClaDS.,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Extinction rate\n(ClaDS)")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..ClaDS.),0.11))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.10,6),color="red")+
  stat_compare_means(comparisons = cm7,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))




cm8 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd8 <- ggplot(data.habi.modified, aes(x=habi,y=Extinction.rate..ClaDS.,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..ClaDS.),0.11))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.10,6),color="red")+
  stat_compare_means(comparisons = cm8,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins
cm9 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd9 <- ggplot(data.biogeo.modified, aes(x=biogeo,y=Extinction.rate..ClaDS.,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified$Extinction.rate..ClaDS.),0.11))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, Inf, label = paste0("n=",n)), vjust=1.5, size = 2.5)+
  geom_hline(yintercept = mean.rate.ClaDS, linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(0.10,6),color="red")+
  stat_compare_means(comparisons = cm9,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))


# arrange plots in canvas
grid.arrange(pd1,pd2,pd3,pd4,pd5,pd6,pd7,pd8,pd9,nrow=3,ncol=3)
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

n.Arthropod <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Arthropod")]
n.Amphibian <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Amphibian")]
n.Mollusc <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Mollusc")]
n.Bird <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Bird")]
n.Reptile <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Reptile")]
n.Plant <- table_taxonomic_group$n[which(as.character(table_taxonomic_group$`<fct>`) == "Plant")]


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


CairoPDF("./Graphs/Pulled_Diversification_Rates.pdf",8,3.9,bg="transparent")
cm1 <- list(c(paste0("Arthropod(",n.Arthropod,")"), paste0("Mollusc(",n.Mollusc,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Arthropod(",n.Arthropod,")"), paste0("Reptile(",n.Reptile,")")),
            c(paste0("Arthropod(",n.Arthropod,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Amphibian(",n.Amphibian,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Mollusc(",n.Mollusc,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Reptile(",n.Reptile,")")),
            c(paste0("Amphibian(",n.Amphibian,")"),paste0("Plant(",n.Plant,")")),
            c(paste0("Reptile(",n.Reptile,")"),paste0("Plant(",n.Plant,")"))
)
pd1<- ggplot(data.taxa.modified3, aes(x=Taxonomic.group2,y=Pulled.Diversification.Rate,color=Taxonomic.group2))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("#bba40d","#0dc2c7","#243763","#db7093","#3C6255"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Taxonomic group")+ylab("Pulled Diversification Rate")+
  #ylim(min(data.taxa.modified$Speciation.rate..RPANDA.),0.7)+ 
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.taxa.modified3$Pulled.Diversification.Rate)-0.001,1.2))+
  #geom_text(data = give.n(quo(Taxonomic.group)),aes(Taxonomic.group, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean(data.taxa.modified3$Pulled.Diversification.Rate), linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(1.1,5),color="red")+
  stat_compare_means(comparisons = cm1,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))



# comparing across habitat types
cm2 <- list(c(paste0("FW(",n.FW,")"), paste0("Wet + Dry(",n.WetDry,")")),
            c(paste0("FW(",n.FW,")"), paste0("Wet(",n.Wet,")")),
            c(paste0("Wet(",n.Wet,")"), paste0("Wet + Dry(",n.WetDry,")")))
pd2 <- ggplot(data.habi.modified2, aes(x=habi,y=Pulled.Diversification.Rate,color=habi))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40","grey70"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Habitat types")+ylab("")+#ylab("Speciation rate RPANDA.")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.habi.modified2$Pulled.Diversification.Rate)-0.001,1.2))+
  #geom_text(data = give.n(quo(Habitat)),aes(Habitat, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean(data.taxa.modified3$Pulled.Diversification.Rate), linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(1.1,5),color="red")+
  stat_compare_means(comparisons = cm2,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position="none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

#comparing across lineages with varying biogeographic origins

cm3 <- list(c(paste0("Gondwanan(",n.Gondwanan,")"), paste0("Asian(",n.Asian,")")))
pd3 <- ggplot(data.biogeo.modified2, aes(x=biogeo,y=Pulled.Diversification.Rate,color=biogeo))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  scale_color_manual(values=c("grey10","grey40"))+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Biogeographic origins")+ylab("")+#ylab("Net-Diversification rate")+
  #ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.7)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2),limits = c(min(data.biogeo.modified2$Pulled.Diversification.Rate)-0.001,1.2))+
  #geom_text(data = give.n(quo(Biogeographic.origin)),aes(Biogeographic.origin, -Inf, label = paste0("n=",n)), vjust=-0.5, size = 2.5)+
  geom_hline(yintercept = mean(data.taxa.modified3$Pulled.Diversification.Rate), linetype = 2,size=0.5,color="#000000")+
  stat_compare_means(label="p.signif",hide.ns = TRUE,ref.group = ".all.",size=6,label.y=rep(1.1,5),color="red")+
  stat_compare_means(comparisons = cm3,hide.ns = TRUE,label="p.signif",size=5,vjust=0.7)+ # Add pairwise comparisons p-value
  theme(#plot.margin = margin(-30.0,5.5,5.5,5.5, "pt"),
    legend.position = "none",axis.text=element_text(color="black"),
    axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),axis.text.y=element_text(size=12),axis.title.x=element_blank(),axis.title.y=element_text(size=14))

grid.arrange(pd1,pd2,pd3,ncol=3)
dev.off()





# SUPPLEMENTARY #
CairoPDF("./Graphs/supplementary_scenarios.pdf",8.3,5.9,bg="transparent")
pb0 <- ggplot(data.taxa.modified, aes(x=Final.Results))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+
  theme(
        axis.text.y=element_text(size=12,angle=0,hjust=1,vjust=0.9),
        axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=20),
        axis.title.y=element_text(size=14))#+coord_flip()

pb1 <-ggplot(data.taxa.modified, aes(x=Final.Results,fill=Taxonomic.group))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Biogeographic origins")+
  scale_fill_manual(name = "Taxonomic\ngroups",values = c("#bba40d","#0dc2c7","#243763","#db7093","#0000ff","#3C6255"))+
  guides(fill=guide_legend(ncol=2))+
  theme(axis.text.y=element_text(size=12,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=14),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()

pb2 <- ggplot(data.taxa.modified, aes(x=Final.Results,fill=Habitat))+geom_bar(width=0.5)+theme_classic()+xlab("Diversification scenarios")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Habitat types",values = c("grey10","grey40","grey70"))+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=12,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=15),axis.title.y=element_text(size=14),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()

pb3 <- ggplot(data.taxa.modified, aes(x=Final.Results,fill=Biogeographic.origin))+geom_bar(width=0.5)+theme_classic()+xlab("")+ylab("No. of lineages")+#ggtitle("Diversification pattern vs. Habitat")+
  scale_fill_manual(name = "Biogeographic\norigins",values = c("grey10","grey40"))+
  #scale_y_discrete(limits=c(1,3,6,9,12))+
  #scale_fill_discrete(name="Diversification pattern")+
  theme(axis.text.y=element_text(size=12,angle=0,hjust=1,vjust=0.9),axis.text.x=element_text(size=12,angle=45,hjust=1,vjust=0.9),
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=14),
        legend.position = "inside",
        legend.position.inside = c(0.5,0.75),legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 5,colour = "black"),
        legend.title = element_text(size = 6,colour = "black"))#+coord_flip()
# SUPPLEMENTARY #
grid.arrange(pb0,pb1,pb2,pb3,
             nrow=2,ncol=2)
dev.off()

CairoPDF("./Graphs/supplementary_pgls.pdf",11.7,8.3,bg="transparent")
# SUPPLEMENTARY #
# only gondwanan
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Gondwanan')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Gondwanan'),c("Net_Diversification_rate (ClaDS)","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_nd_gond <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                           "Net-diversification rate","log(Species Richness)",min(dt.temp$`Net_Diversification_rate (ClaDS)`),0,"Gondwanan")

# only asian
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Asian')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Asian'),c("Net_Diversification_rate (ClaDS)","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_nd_sea <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "Net-diversification rate","log(Species Richness)",min(dt.temp$`Net_Diversification_rate (ClaDS)`),0,"Asian")

# only wet habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet'),c("Net_Diversification_rate (ClaDS)","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_nd_wet <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "Net-diversification rate","log(Species Richness)",min(dt.temp$`Net_Diversification_rate (ClaDS)`),0,"Wet Habitat")

# only wet+dry habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet + Dry')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet + Dry'),c("Net_Diversification_rate (ClaDS)","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_nd_wet_dry <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                              "Net-diversification rate","log(Species Richness)",min(dt.temp$`Net_Diversification_rate (ClaDS)`),0,"Wet+Dry Habitat")

# only gondwanan
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Gondwanan')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Gondwanan'),c("Clade_age","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_ca_gond <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                           "log(Clade age)","log(Species Richness)",0,0,"Gondwanan")

# only asian
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Asian')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Asian'),c("Clade_age","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_ca_sea <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "log(Clade age)","log(Species Richness)",0,0,"Asian")

# only wet habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet'),c("Clade_age","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_ca_wet <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "log(Clade age)","log(Species Richness)",0,0,"Wet Habitat")

# only wet+dry habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet + Dry')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet + Dry'),c("Clade_age","Species_Richness")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_sr_ca_wet_dry <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                              "log(Clade age)","log(Species Richness)",0,0,"Wet+Dry Habitat")

# only gondwanan
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Gondwanan')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Gondwanan'),c("Clade_age","Net_Diversification_rate (ClaDS)")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_nd_ca_gond <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                           "log(Clade age)","Net-diversification rate",0,min(dt.temp$`Net_Diversification_rate (ClaDS)`),"Gondwanan")

# only se-asian
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Biogeographic.origin != 'Asian')])
dt.temp <- data.pgls[which(data.pgls$Biogeographic.origin == 'Asian'),c("Clade_age","Net_Diversification_rate (ClaDS)")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_nd_ca_sea <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "log(Clade age)","Net-diversification rate",0,min(dt.temp$`Net_Diversification_rate (ClaDS)`),"Asian")

# only wet habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet'),c("Clade_age","Net_Diversification_rate (ClaDS)")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_nd_ca_wet <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                          "log(Clade age)","Net-diversification rate",0,min(dt.temp$`Net_Diversification_rate (ClaDS)`),"Wet Habitat")

# only wet+dry habitat
tree.temp <- drop.tip(super.tree,rownames(data.pgls)[which(data.pgls$Habitat != 'Wet + Dry')])
dt.temp <- data.pgls[which(data.pgls$Habitat == 'Wet + Dry'),c("Clade_age","Net_Diversification_rate (ClaDS)")]
dt.temp <- dt.temp[tree.temp$tip.label,]
pgls_nd_ca_wet_dry <- mk_regression_plot_PGLS(dt.temp,tree.temp,
                                              "log(Clade age)","Net-diversification rate",0,min(dt.temp$`Net_Diversification_rate (ClaDS)`),"Wet+Dry Habitat")

grid.arrange(pgls_main_sr_ca,
             pgls_main_sr_nd_ClaDS,
             pgls_main_nd_ca_ClaDS,
             arrangeGrob(pgls_sr_ca_gond,pgls_sr_ca_sea,pgls_sr_ca_wet,pgls_sr_ca_wet_dry,nrow=2,ncol=2),
             arrangeGrob(pgls_sr_nd_gond,pgls_sr_nd_sea,pgls_sr_nd_wet,pgls_sr_nd_wet_dry,nrow=2,ncol=2),
             arrangeGrob(pgls_nd_ca_gond,pgls_nd_ca_sea,pgls_nd_ca_wet,pgls_nd_ca_wet_dry,nrow=2,ncol=2),
             nrow=2,ncol=3)

# SUPPLEMENTARY #
dev.off()

# Supplementary box plots

CairoPDF("./Graphs/supplementary_scenarios_drivers.pdf",8.3,11.7,bg="transparent")
cm_TDD <- list( c("TDD","No TDD"))
cm_DDD <- list(c("DDD","No DDD"))

pTDD_SR <- ggplot(data.taxa.modified, aes(x=TDD,y=No..of.species))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Drivers")+ylab("Species Richness")+
  ylim(0,110)+
  stat_compare_means(comparisons = cm_TDD,hide.ns = TRUE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

pDDD_SR <- ggplot(data.taxa.modified, aes(x=DDD,y=No..of.species))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Drivers")+ylab("Species Richness")+
  ylim(0,110)+
  stat_compare_means(comparisons = cm_DDD,hide.ns = TRUE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

pTDD_DR <- ggplot(data.taxa.modified, aes(x=TDD,y=Net.diversification.rate..RPANDA.))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Drivers")+ylab("Net-diversification rate")+
  ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.9)+
  stat_compare_means(comparisons = cm_TDD,hide.ns = TRUE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

pDDD_DR <- ggplot(data.taxa.modified, aes(x=DDD,y=Net.diversification.rate..RPANDA.))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Drivers")+ylab("Net-diversification rate")+
  ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.9)+
  stat_compare_means(comparisons = cm_DDD,hide.ns = FALSE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

cm_DS <- list(c("SC1","SC2"),c("SC1","SC2+Shift"),c("SC1","SC3+Shift"),c("SC1","Shift"))
pDS_SR <- ggplot(data.taxa.modified, aes(x=Final.Results,y=No..of.species))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Scenarios")+ylab("Species Richness")+
  ylim(0,130)+
  stat_compare_means(comparisons = cm_DS,hide.ns = FALSE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

pDS_DR <- ggplot(data.taxa.modified, aes(x=Final.Results,y=Net.diversification.rate..RPANDA.))+
  geom_boxplot(lwd=0.35,outlier.alpha = 0)+
  geom_jitter(shape=16, position=position_jitter(0.2),size=1)+
  theme_classic()+xlab("Scenarios")+ylab("Net-diversification rate")+
  ylim(min(data.taxa.modified$Net.diversification.rate..RPANDA.),0.75)+
  stat_compare_means(comparisons = cm_DS,hide.ns = FALSE,label="p.signif",size=7)+ # Add pairwise comparisons p-value
  theme(plot.margin = margin(5.5,5.5,35.5,5.5, "pt"),axis.text.y=element_text(size=12),axis.title.y=element_text(size=14))

grid.arrange(arrangeGrob(pTDD_SR,pDDD_SR,ncol=2),
             arrangeGrob(pTDD_DR,pDDD_DR,ncol=2),
             arrangeGrob(pDS_SR,pDS_DR,ncol=2),ncol=1)
dev.off()
pdf_combine(input=c("./Graphs/supplementary_FAMD.pdf","./Graphs/supplementary_scenarios_drivers.pdf"),output = "./Graphs/supplementary_FAMD2.pdf")

# species richness and div rate box plots comparing TD, DD and Scenarios #
# supplementary tables #
# LTT and RTT supplementary #

############################################################################################################################
############################################################################################################################


t <- read.tree("/media/evol-eco/New_volume/Pragyadeep_RPANDA_analysis/Summary and Graphs/SUPER_TREE.nwk")
CairoPNG("/media/evol-eco/New_volume/Pragyadeep_RPANDA_analysis/Summary and Graphs/Graphs/For_MS_Final/supplementary/SI5_pgls_Tree.png",
         1170,830,bg="transparent")
p <- ggtree(t) +
  coord_geo(xlim = c(-1650, 0), ylim = c(-2, Ntip(t)), neg = TRUE, abbrv = TRUE) +
  geom_tiplab(as_ylab = TRUE,face="italic",)+
  scale_x_continuous(breaks = seq(-1600, 0, 250), labels = abs(seq(-1600, 0, 250))) +
  theme_minimal()
revts(p)
dev.off()
############################################################################################################################
