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
