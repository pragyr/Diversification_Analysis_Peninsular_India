# Diversification in the Peninsular Indian Plate
## Authors
- Pragyadeep Roy
- Jahnavi Joshi

## Relevant Manuscript Title:
Idiosyncrasies unveiled: examining the pace, patterns and predictors of biotic diversification in peninsular India

## Summary
Hello,
Welcome to the this project. It deals with understanding the pace, pattern and predictors of in-situ diversification in the Peninsular Indian Plate (PIP). The data that we have used in the project includes dated tree files and paleoclimate data.

## Data used in the study
All the data would be present in the "data" directory.

### Phylogenetic trees
The "all_trees" folder as the name suggests contains all the dated phylogenetic trees of the clades of our interest. These files were used the analyses that involved using the homogeneous birth-death models of package RPANDA and to run the CoMET (CPP on Mass Extinction Time) analysis. Few of the clades used in our analysis have smaller clade size and it has been seen that diversification estimates show more uncertainty with smaller clades. Hence, in addition to using the homogenous birth-death model from RPANDA and the CoMET analysis, we analysed our data using a hetergenous birth-death model - ClaDS (Cladogenetic Diverisification rate Shift) model. And, this model was implemented on larger trees which contained those small clades of our interest and those "super trees" were kept in the "all_trees_or_super_trees" folder. For clades with fairly large clade size we have used the same phylogenetic tree as is present in the "all_trees" folder.

### Paleoclimatic data
Paleoclimate data includes global paleotemperature (for the cenozoic era), Himalayan elevation, Expansion of C4-plants through time. The "paleoclimate_data" folder doesn't however include the temperature data. It is available in the RPANDA package
```{r}
library(RPANDA)
data(InfTemp)
head(InfTemp)
```

Outside these folders the file - "Final_summary_updated.csv" summary information about the identity of each of the lineages and the corresponding inference about its diversification. More details about the analyses can be found in the supplementary information of the study. The "SUPER_TREE.tre" file is a dated phylogenetic tree connecting all the clades used in the study. It was generated in [www.timetree.org](www.timetree.org) and also using the information in the clades and respective studies and was used in one of the analysis (PGLS) in the study.

## Codes
All the relevant codes are present in the "codes" folder. Here is the list of all the codes.

Although the code names are self-explanatory, here's a brief account of the actions that the codes would perform:

"1_spline_interpolation_PR_JJ.R" - for smoothing the paleoclimatic variables used in the analyses

"2_COMET_analysis_and_CRABS_PR_JJ.R" - for performing the CoMET analysis, to extract information about diverisification rates (speciation and extinction) through time and times of strong episodic rate shifts, and to assess congruence classes using the CRABS analysis

"3_rpanda_code_and_RTT_comparisons_PR_JJ.R" - for estimating diversification rates, patterns and drivers using homogenous birth-death models from the package RPANDA and to compare the estimated rates through time with that of CoMET analysis

"4_run_ClaDS_in_Julia_PR_JJ.jl" - to peform the ClaDS analysis in Julia and to save the results in Rdata format

"5_export_ClaDS_data_and_summarise_in_R_PR_JJ.R" - to load the Rdata result files relevant to the ClaDS analysis and to extract rates through time for clades of interest

"6_pulled_diversification_rates_finding_optimal_age_gridding.R" - to fit the homogenous pulled diversification rate model on the data to assess general trends in pulled rates

"7_make_plots_PR_JJ.R" - to make plots submitted with this manuscript

## Smoothing paleoclimatic data through spline interpolation
```{r}
library(npreg)
library(RPANDA)
# MAIN CODE STARTS HERE
# Provide all your tree files in a folder with the name - 'all_trees' and then set the working directory to the path of that 'all_trees' folder
# PLEASE CHECK THE FOLLOWING TWO LINES
# if you are using the linux/mac terminal
setwd(paste0(getwd(),"/../data/all_trees"))

# if you are using RStudio IDE
setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"../data/paleoclimate_data/"))


# raw paleodata
data(InfTemp)
C_dat_potwar <- read.csv("C_dat2_Potwar.csv")
C_dat_nwi <- read.csv("C_dat2_NWI.csv")
Him.oro_dat <- read.csv("Him_Oro_data.csv")

# spline interpolation
inftemp.spline <- ss(InfTemp$Age,InfTemp$Temperature)
C_dat_potwar.spline <- ss(C_dat_potwar$Age,C_dat_potwar$Delta.13_C)
C_dat_nwi.spline <- ss(C_dat_nwi$Age,C_dat_nwi$Delta.13_C)
Him.oro_dat.spline <- ss(Him.oro_dat$Age,Him.oro_dat$Elevation)

# summary of interpolations
summary(inftemp.spline)
summary(C_dat_potwar.spline)
summary(C_dat_nwi.spline)
summary(Him.oro_dat.spline)

# plotting smoothened curves
CairoPDF("splines.pdf",11.7,8.3,bg="transparent")
par(mfrow=(c(3,1)),mar=c(5.1,14.1,4.1,2.1))
plot(InfTemp$Age,InfTemp$Temperature,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab='',ylab='',xlim=c(max(InfTemp$Age),0),cex.axis=2)
title(ylab="Temperature (\u00B0C)",line=4,cex.lab=2.5)
lines(inftemp.spline$x,inftemp.spline$y,lwd=1.5)

plot(C_dat_potwar$Age,C_dat_potwar$Delta.13_C,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab="",ylab="",xlim=c(max(InfTemp$Age),0),cex.axis=2)
title(ylab=paste0("\u03B4","13 C"),line=4,cex.lab=2.5)
lines(C_dat_potwar.spline$x,C_dat_potwar.spline$y,lwd=1.5)

#plot(C_dat_nwi.spline,xlab="Age",ylab="Delta13.C")
plot(Him.oro_dat$Age,Him.oro_dat$Elevation,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab="Age (Mya)",ylab='',xlim=c(max(InfTemp$Age),0),cex.axis=2,cex.lab=2.5)
title(ylab="Elevation (m)",line=4,cex.lab=2.5)
lines(Him.oro_dat.spline$x,Him.oro_dat.spline$y,xlab='Age',ylab="Elevation",lwd=1.5)
dev.off()



# make them individual dataframes
smooth.Temp.dat <- as.data.frame(cbind("Age"=inftemp.spline$x,"Temperature"=inftemp.spline$y))
smooth.C_dat_potwar <- as.data.frame(cbind("Age"=C_dat_potwar.spline$x,"Delta.13_C"=C_dat_potwar.spline$y))
smooth.C_dat_nwi <- as.data.frame(cbind("Age"=C_dat_nwi.spline$x,"Delta.13_C"=C_dat_nwi.spline$y))
smooth.Him.oro_dat <- as.data.frame(cbind("Age"=Him.oro_dat.spline$x,"Elevation"=Him.oro_dat.spline$y))

# write the dataframes as files
write.csv(smooth.Temp.dat,file = "Smooth_Inftemp.csv")
write.csv(smooth.C_dat_potwar,file = "Smooth_C_dat_Potwar.csv")
write.csv(smooth.C_dat_nwi,file = "Smooth_C_dat_nwi.csv")
write.csv(smooth.Him.oro_dat,file = "Smooth_Him_oro_dat.csv")

pdf("./Curve_Differentials.pdf",8.3,11.7,bg="transparent")
par(mfrow=c(3,1))
# plot the derivatives of the functions
# C4 Expansion
Y<-diff(smooth.C_dat_potwar$Delta.13_C)
X<-diff(smooth.C_dat_potwar$Age)
dy_dx.c4 <- abs(Y/X)
plot(smooth.C_dat_potwar$Age[-101],dy_dx.c4,type="l",xlab="Time-t (Mya)",ylab="d(Delta.13_C)/dt",xlim=c(0,max(smooth.Temp.dat$Age)))

# Himalayan Elevation
Y<-diff(smooth.Him.oro_dat$Elevation)
X<-diff(smooth.Him.oro_dat$Age)
dy_dx.him.oro <- abs(Y/X)
plot(smooth.Him.oro_dat$Age[-57],dy_dx.him.oro,type="l",xlab="Time-t (Mya)",ylab="d(Elevation)/dt",xlim=c(0,max(smooth.Temp.dat$Age)))

# Temperature
Y<-diff(smooth.Temp.dat$Temperature)
X<-diff(smooth.Temp.dat$Age)
dy_dx.temp <- abs(Y/X)
plot(smooth.Temp.dat$Age[-15363],dy_dx.temp,type="l",xlab="Time-t (Mya)",ylab="d(Temperature)/dt",xlim=c(0,max(smooth.Temp.dat$Age)))
dev.off()

write.csv(as.data.frame(cbind("Age"=smooth.Temp.dat$Age[-15363],dy_dx.temp)),file = "Smooth_Inftemp_derivative.csv")
write.csv(as.data.frame(cbind("Age"=smooth.C_dat_potwar$Age[-101],dy_dx.c4)),file = "Smooth_C_dat_Potwar_derivative.csv")
write.csv(as.data.frame(cbind("Age"=smooth.Him.oro_dat$Age[-57],dy_dx.him.oro)),file = "Smooth_Him_oro_dat_derivative.csv")
```
