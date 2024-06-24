# Plotting diversity metrics along chromosome 6 and 9 to look to signatures of selective sweeps:

rm(list=ls())
setwd("~/Desktop/oysterdata")

# Reading in windowed pi (nucleotide diversity) data:
chr6 <- read.delim("chr6.windowed.pi", header=TRUE)
chr9 <- read.delim("chr9.windowed.pi", header=TRUE)

#Plotting to highlight the outlier data points:
plot(1, type="n", xlab="LG6 position", ylab="Pi", xlim=c(0, 60140001), ylim=c(0, 0.05))
points(chr6$BIN_START[1:3412], chr6$PI[1:3412], col = "lightgray")
points(chr6$BIN_START[3423:5441], chr6$PI[3423:5441], col = "lightgray")
points(chr6$BIN_START[3413:3422], chr6$PI[3413:3422], col = "red")

#Plotting to highlight the outlier data points:
plot(1, type="n", xlab="LG9 position", ylab="Pi", xlim=c(0, 37060001), ylim=c(0, 0.05))
points(chr9$BIN_START[1:2909], chr9$PI[1:2909], col = "lightgray")
points(chr9$BIN_START[2919:2921], chr9$PI[2919:2921], col = "lightgray")
points(chr9$BIN_START[2910:2918], chr9$PI[2910:2918], col = "red")

# Now doing the same for Tajima's D:
# Reading in data:
chr6D <- read.delim("chr6_TajimasD.txt", header=TRUE)
chr9D <- read.delim("chr9_TajimasD.txt", header=TRUE)

#Plotting to highlight the outlier data points:
plot(1, type="n", xlab="LG6 position", ylab="Tajimas D", xlim=c(0, 60140000), ylim=c(-3, 3))
points(chr6D$BIN_START[1:3798], chr6D$TajimaD[1:3798], col = "lightgray")
points(chr6D$BIN_START[3809:6013], chr6D$TajimaD[3809:6013], col = "lightgray")
points(chr6D$BIN_START[3799:3808], chr6D$TajimaD[3799:3808], col = "red")

#Plotting to highlight the outlier data points:
plot(1, type="n", xlab="LG9 position", ylab="Tajimas D", xlim=c(0, 37060000), ylim=c(-3, 3))
points(chr9D$BIN_START[1:3688], chr9D$TajimaD[1:3688], col = "lightgray")
points(chr9D$BIN_START[3699:3705], chr9D$TajimaD[3699:3705], col = "lightgray")
points(chr9D$BIN_START[3689:3698], chr9D$TajimaD[3689:3698], col = "red")
