# Plotting G scores as a function of sample size

rm(list=ls())
setwd("~/Desktop/oysterdata")

# Reading in data:
Gscores14 <- read.csv("TotalGscores14.csv", header=TRUE)
samplesize14 <- read.csv("samplesizes14.csv", header=TRUE)
Gscores24 <- read.csv("TotalGscores24.csv", header=TRUE)
samplesize24 <- read.csv("samplesizes24.csv", header=TRUE)

# Trying to fit a linear model, didn't work very well....
col2<-lm((as.numeric(Gscores14[,2])~as.numeric(samplesize14[,2])), na.action=na.omit)
x <- seq(0, 6133, 1)
y <- predict(col2, list(samplesize14[,2]=x),type="response")


# Testing to plot a few families in 14

plot(Gscores14[,2]~samplesize14[,2], col="red")
abline(lm(as.numeric(Gscores14[,2])~as.numeric(samplesize14[,2])))
points(Gscores14[,3]~samplesize14[,3], col="blue")
points(Gscores14[,4]~samplesize14[,4], col="green")
abline(a=0, b=1)

# And also in 24...

plot(Gscores24[,2]~samplesize24[,2], col="red")
points(Gscores24[,3]~samplesize24[,3], col="blue")
points(Gscores24[,4]~samplesize24[,4], col="green")
abline(a=0, b=1)

# Now, plotting all families in 14:
par(mfrow=c(5,7))
for (i in 2:36) {
  fam <- colnames(Gscores14[i])
  plot(Gscores14[,i]~samplesize14[,i], xlab="Nloci", ylab="Gscore", main=fam)
  abline(a=0, b=1)
}
par(mfrow=c(1,1))

# And all families in 24:
par(mfrow=c(6,7))
for (i in 2:38) {
  fam <- colnames(Gscores24[i])
  plot(Gscores24[,i]~samplesize24[,i], xlab="Nloci", ylab="Gscore", main=fam)
  abline(a=0, b=1)
}
par(mfrow=c(1,1))

# Now, let's plot each family in both salinities, so that we can compare:

# Have to remove families that exist only in one salinity
# These are: 
# In 14: F14_M13, F03_M03
# In 24: F13_M15, F01_M02, F01_M03, F06_M05

Gscores14sub = subset(Gscores14, select = -c(F14_M13, F03_M03))
Gscores24sub = subset(Gscores24, select = -c(F13_M15, F01_M02, F01_M03, F06_M05))
samplesize14sub = subset(samplesize14, select = -c(F14_M13, F03_M03))
samplesize24sub = subset(samplesize24, select = -c(F13_M15, F01_M02, F01_M03, F06_M05))

# Plotting:
par(mfrow=c(5,7))
for (i in 2:34) {
  fam <- colnames(Gscores14sub[i])
  plot(Gscores24sub[,fam]~samplesize24sub[,fam], xlab="Nloci", ylab="Gscore", main=fam, col= rgb(red = 0, green = 0, blue = 1, alpha = 0.25))
  points(Gscores14sub[,fam]~samplesize14sub[,fam], col=rgb(red = 1, green = 0, blue = 0, alpha = 0.25))
  abline(a=0, b=1)
}
par(mfrow=c(1,1))


#Let's try doing some multiregressions:

Gscores <- read.delim("TotalGscores_all.txt", header=TRUE)
samplesize <- read.delim("samplesizes_all.txt", header=TRUE)

test<-lm((as.numeric(Gscores[,3])~as.numeric(samplesize[,3])*Gscores[,2]), na.action=na.omit)
summary(test)
plot(test)

#Plotting with ggplot:

library(ggplot2)

# Testing to plot family F01_M01 as a test:
x <- ggplot(data = data.frame(x = samplesize[,3], y = Gscores[,3]), aes(x = x, y = y, col = as.factor(Gscores[,2]))) +
  labs(y="delta G", x="Nloci", title="F01_M01") +
  guides(col=guide_legend(title="Sal")) +
  scale_color_brewer(palette="Set1") +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  geom_abline(slope=1, intercept = 0, linetype = "dotted")
x

# Now, plotting all of them:

my_plots <- list()
for (i in 3:35) {
  fam <- colnames(Gscores[i])
  x <- ggplot(data = data.frame(x = samplesize[,i], y = Gscores[,i]), aes(x = x, y = y, col = as.factor(Gscores[,2]))) +
    labs(y="delta G", x="Nloci", title=fam) +
    guides(col=guide_legend(title="Sal")) +
    scale_color_brewer(palette="Set1") +
    geom_point() +
    geom_smooth(method="lm", se=FALSE) +
    geom_abline(slope=1, intercept = 0, linetype = "dotted")
  print(x)
  my_plots[[i]] = x
}

library(ggpubr)
figure <- ggarrange(my_plots[[3]], my_plots[[4]], my_plots[[5]], my_plots[[6]], my_plots[[7]], my_plots[[8]], my_plots[[9]], my_plots[[10]], my_plots[[11]], my_plots[[12]], my_plots[[13]], my_plots[[14]], my_plots[[15]], my_plots[[16]], my_plots[[17]], my_plots[[18]], my_plots[[19]], my_plots[[20]], my_plots[[21]], my_plots[[22]], my_plots[[23]], my_plots[[24]], my_plots[[25]], my_plots[[26]], my_plots[[27]], my_plots[[28]], my_plots[[29]], my_plots[[30]], my_plots[[31]], my_plots[[32]], my_plots[[33]], my_plots[[34]], my_plots[[35]],
                    ncol = 7, nrow = 5)
figure


