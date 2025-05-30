library(npreg)
library(RPANDA)
# MAIN CODE STARTS HERE
# Provide all your tree files in a folder with the name - 'all_trees' and then set the working directory to the path of that 'all_trees' folder


# raw paleodata
data(InfTemp) # data for global paleotemperature - available with the RPANDA package
C_dat_potwar <- read.csv("data/paleoclimate_data/C_dat2_Potwar.csv") # data for the expansion of C4 plants - collected from Clift and Webb 2017 and Quade et al. 1989
Him.oro_dat <- read.csv("data/paleoclimate_data/Him_Oro_data.csv") # data for himalayan orogeny - collected from Ding et al. 2017

# spline interpolation
inftemp.spline <- ss(InfTemp$Age,InfTemp$Temperature)
C_dat_potwar.spline <- ss(C_dat_potwar$Age,C_dat_potwar$Delta.13_C)
Him.oro_dat.spline <- ss(Him.oro_dat$Age,Him.oro_dat$Elevation)

# plotting smoothened curves
pdf("data/paleoclimate_data/splines.pdf",11.7,8.3,bg="transparent")
par(mfrow=(c(3,1)),mar=c(5.1,14.1,4.1,2.1))
plot(InfTemp$Age,InfTemp$Temperature,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab='',ylab='',xlim=c(max(InfTemp$Age),0),cex.axis=2)
title(ylab="Temperature (\u00B0C)",line=4,cex.lab=2.5)
lines(inftemp.spline$x,inftemp.spline$y,lwd=1.5)

plot(C_dat_potwar$Age,C_dat_potwar$Delta.13_C,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab="",ylab="",xlim=c(max(InfTemp$Age),0),cex.axis=2)
title(ylab=paste0("\u03B4","13 C"),line=4,cex.lab=2.5)
lines(C_dat_potwar.spline$x,C_dat_potwar.spline$y,lwd=1.5)

plot(Him.oro_dat$Age,Him.oro_dat$Elevation,type="l",lty=1,lwd=0.1,col=rgb(red=0,green=0,blue=0,alpha=0.5),
     xlab="Age (Mya)",ylab='',xlim=c(max(InfTemp$Age),0),cex.axis=2,cex.lab=2.5)
title(ylab="Elevation (m)",line=4,cex.lab=2.5)
lines(Him.oro_dat.spline$x,Him.oro_dat.spline$y,xlab='Age',ylab="Elevation",lwd=1.5)
dev.off()



# make them individual dataframes
smooth.Temp.dat <- as.data.frame(cbind("Age"=inftemp.spline$x,"Temperature"=inftemp.spline$y))
smooth.C_dat_potwar <- as.data.frame(cbind("Age"=C_dat_potwar.spline$x,"Delta.13_C"=C_dat_potwar.spline$y))
smooth.Him.oro_dat <- as.data.frame(cbind("Age"=Him.oro_dat.spline$x,"Elevation"=Him.oro_dat.spline$y))

# write the dataframes as files
write.csv(smooth.Temp.dat,file = "data/paleoclimate_data/Smooth_Inftemp.csv")
write.csv(smooth.C_dat_potwar,file = "data/paleoclimate_data/Smooth_C_dat_Potwar.csv")
write.csv(smooth.Him.oro_dat,file = "data/paleoclimate_data/Smooth_Him_oro_dat.csv")
