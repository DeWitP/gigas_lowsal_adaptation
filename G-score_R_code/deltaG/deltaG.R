#Testing to calculate standardized G-values and p-values for chromosome LR761634 in salinity 24. 

#Clearing env
rm(list = ls())
setwd("~/Desktop/oysterdata")
#Reading data:
TotG14<-read.csv("embryos14.TotGmat.csv", header=TRUE)
df14<-read.csv("embryos14.df.csv", header=TRUE)
TotG24<-read.csv("embryos24.TotGmat.csv", header=TRUE)
df24<-read.csv("embryos24.df.csv", header=TRUE)

#Creating new dataframes to standaridize the G values into:
Gstand14 <- TotG14
Gstand24 <- TotG24

#Standardizing by number of loci:
Gstand14[,3:37] = TotG14[,3:37]/df14[,3:37]
Gstand24[,3:39] = TotG24[,3:39]/df24[,3:39]

# Have to remove families that exist only in one salinity
# These are: 
# In 14: F14_M13, F03_M03
# In 24: F13_M15, F01_M02, F01_M03, F06_M05

Gstand14sub = subset(Gstand14, select = -c(F14_M13, F03_M03))
Gstand24sub = subset(Gstand24, select = -c(F13_M15, F01_M02, F01_M03, F06_M05))

#Now, calculating the difference in G value due to salinity (deltaG):
# A positive deltaG means that the deviation from Mendelian inheritance is stronger in 14 than in 24.

deltaG = Gstand14sub
deltaG[,3:35] = Gstand14sub[,3:35]-Gstand24sub[,3:35]
deltaG$LG <- as.factor(deltaG$LG)

# Probably dividing by zeroes above has introduced NaNs in deltaG, need to remove them.
deltaG = na.omit(deltaG)

# Now, calculating the average deltaG across all families and also bootstrapping the 95 % confidence intervals.

avgdeltaG=data.frame(matrix(ncol = 6, nrow = nrow(deltaG)))
colnames(avgdeltaG)<-c("LG","Block","meandeltaG","lower95","upper95","track_no")
avgdeltaG[,1:2] = deltaG[,1:2]
avgdeltaG[,6] <- c(1:nrow(deltaG))

library(simpleboot)

for(x in 1:nrow(deltaG)) {
  avgdeltaG[x,3] = mean(as.numeric(deltaG[x,3:35]))
  bootstrap <- one.boot(as.numeric(as.vector(deltaG[x,3:35])), mean, trim = 0, R=10000)
  CI <- boot.ci(bootstrap, type = "perc")
  avgdeltaG[x,4] = CI$percent[4]
  avgdeltaG[x,5] = CI$percent[5]
}

#Saving the data frame to disk:
write.csv(avgdeltaG, file="avgdeltaG.csv")
save(avgdeltaG, file="avgdeltaG.Rdata")

#For next time, can skip thebootstrap and just do
# load(file="avgdeltaG.Rdata")


# Finally, plotting:

# Plotting the distribution of mean deltaG as a histogram
hist(avgdeltaG$meandeltaG, xlim=c(0,10), breaks = 500, main ="", xlab="Mean delta G")

#  Now let's try to plot means and confidence intervals along the chromosomes with ggplot!

library(ggplot2)
ggplot(avgdeltaG,aes(x=track_no, y=meandeltaG, col=as.factor(LG))) +
  labs(y="delta G", colour="Linkage Group") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point() +
  theme(panel.background = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())



# Plotting family-wise distributions of deltaG in random subsets of 20 blocks
random = runif(20, min = 0, max = nrow(deltaG))
dgblock <- deltaG[random,]
par(mfrow=c(4,5))
for(x in 1:nrow(dgblock)) {
  hist(c(as.numeric(dgblock[x,3:35])), breaks = 8, main = "")
}
par(mfrow=c(1,1))

# Looks normal in the great majority of cases.

# Extracting outlier blocks (deltaG > 10) and plotting histograms
outliers <- subset(avgdeltaG, meandeltaG > 10)
outlierrows = c("1","1620","2398","2945","3090","3096","3523","3775","5320")
par(mfrow=c(3,3))
for(x in 1:9) {
  hist(c(as.numeric(deltaG[outlierrows[x],3:35])), breaks = 8, main = "", xlim = c(-200,1200),xlab = paste(deltaG[outlierrows[x],1],deltaG[outlierrows[x],2],sep="; "))
}
par(mfrow=c(1,1))

#Does NOT look normal in several of them!

# Plotting pairwise cross-correlations across families in the 9 outlier blocks. 
test <- deltaG[3:35]
par(mfrow=c(9,9))
for (x in 1:9){
  for (y in 1:9){
    plot(as.numeric(test[outlierrows[y],])~as.numeric(test[outlierrows[x],]), main = paste(deltaG[outlierrows[x],1],deltaG[outlierrows[x],2],deltaG[outlierrows[y],1],deltaG[outlierrows[y],2],sep="; "))
    
  }
}
par(mfrow=c(1,1))

# Interestingly, most of these are completely orthogonal, so there's high deltaGs in either one or the other, never both.
