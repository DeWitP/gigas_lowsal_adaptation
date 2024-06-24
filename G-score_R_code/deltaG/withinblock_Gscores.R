#!/usr/bin/env Rscript
setwd("~/Desktop/oysterdata")

# Loading in the data:
#Expected proportions (calculated in the "parents_LR761634.1.R" script):
expect_ref_AF <- read.delim("parent_outlier_AFs.txt", header=TRUE)

#Reading in the observed counts for salinity 14 (where "NA" has been replaced with zeroes):
OBS_REF_14 <- read.delim("outliers_14_REF.txt", header=TRUE)
OBS_ALT_14 <- read.delim("outliers_14_ALT.txt", header=TRUE)

#And the same for salinity 24):
OBS_REF_24 <- read.delim("outliers_24_REF.txt", header=TRUE)
OBS_ALT_24 <- read.delim("outliers_24_ALT.txt", header=TRUE)


#Creating a new dataframe to calculate the total observed counts, sal 14 (needed to filter based on expected counts >= 5 below - a prerequisite for the G test to work well)
TOTAL_14 <- data.frame(matrix(ncol = 39, nrow = 5980))
#provide column names
colnames(TOTAL_14) <- colnames(OBS_REF_14)
TOTAL_14$CHROM <- OBS_REF_14$CHROM
TOTAL_14$POS <- OBS_REF_14$POS
TOTAL_14$REF <- OBS_REF_14$REF
TOTAL_14$ALT <- OBS_REF_14$ALT
TOTAL_14[,5:39] = OBS_REF_14[,5:39] + OBS_ALT_14[,5:39]
head(TOTAL_14)

#Creating a new dataframe to fill with the observed Major allele frequency at sal 14, for information later.
MaAF_14 <- data.frame(matrix(ncol = 39, nrow = 5980))
#provide column names
colnames(MaAF_14) <- colnames(OBS_REF_14)
MaAF_14$CHROM <- OBS_REF_14$CHROM
MaAF_14$POS <- OBS_REF_14$POS
MaAF_14$REF <- OBS_REF_14$REF
MaAF_14$ALT <- OBS_REF_14$ALT
MaAF_14[,5:39] = (OBS_REF_14[,5:39] / TOTAL_14[,5:39])
head(MaAF_14)

#Saving the data frame to disk:
write.csv(MaAF_14, file="outlier_14_MaAFs.csv")

#Creating a new dataframe to calculate the total observed counts, sal 24 (needed to filter based on expected counts >= 5 below - a prerequisite for the G test to work well)
TOTAL_24 <- data.frame(matrix(ncol = 41, nrow = 5980))
#provide column names
colnames(TOTAL_24) <- colnames(OBS_REF_24)
TOTAL_24$CHROM <- OBS_REF_24$CHROM
TOTAL_24$POS <- OBS_REF_24$POS
TOTAL_24$REF <- OBS_REF_24$REF
TOTAL_24$ALT <- OBS_REF_24$ALT
TOTAL_24[,5:41] = OBS_REF_24[,5:41] + OBS_ALT_24[,5:41]
head(TOTAL_24)

# G-tests for each SNP, sal 14:

# Creating an empty dataframe to fill with the G values for each SNP
G_14 <- data.frame(matrix(ncol = 39, nrow = 5980))
colnames(G_14) <- colnames(OBS_REF_14)
G_14$CHROM <- OBS_REF_14$CHROM
G_14$POS <- OBS_REF_14$POS
G_14$REF <- OBS_REF_14$REF
G_14$ALT <- OBS_REF_14$ALT

# Need this package for the G tests:
library(RVAideMemoire)

# Doing G tests as a for loop

fam14 <- c("F10_M10","F10_M11","F10_M12","F11_M10","F11_M11","F11_M12","F12_M10","F12_M11","F12_M12","F13_M13","F13_M14","F14_M13","F14_M14","F14_M15","F15_M13","F15_M14","F15_M15","F01_M01","F02_M01","F02_M02","F02_M03","F03_M01","F03_M02","F03_M03","F06_M04","F06_M06","F07_M07","F07_M08","F07_M09","F08_M07","F08_M08","F08_M09","F09_M07","F09_M08","F09_M09")

for (x in 1:35) {
  fam = fam14[x]
  for (row in 1:nrow(G_14)) {
      pos <- OBS_REF_14[row, "POS"]
      line <- which(expect_ref_AF$POS == pos)
      #print(row)
      #print(pos)
      #print(line)
      if ((expect_ref_AF[line,fam] < 1) & (expect_ref_AF[line,fam] > 0) & (TOTAL_14[row,fam]*expect_ref_AF[line,fam] >= 5) & (TOTAL_14[row,fam]*(1-expect_ref_AF[line,fam]) >= 5)) {
          test <-G.test(x=c(OBS_REF_14[row,fam], OBS_ALT_14[row,fam]), p=c(expect_ref_AF[line,fam],(1-expect_ref_AF[line,fam])))
          print(as.numeric(test$statistic))
      if(as.numeric(test$statistic) != "Inf") { 
      G_14[row,fam] = as.numeric(test$statistic)
        }
      
      }
  }
}
# G-tests for each SNP, sal 24:

# Creating an empty dataframe to fill with the G values for each SNP
G_24 <- data.frame(matrix(ncol = 41, nrow = 5980))
colnames(G_24) <- colnames(OBS_REF_24)
G_24$CHROM <- OBS_REF_24$CHROM
G_24$POS <- OBS_REF_24$POS
G_24$REF <- OBS_REF_24$REF
G_24$ALT <- OBS_REF_24$ALT

# Doing G tests as a for loop

fam24 <- c("F10_M10","F10_M11","F10_M12","F11_M10","F11_M11","F11_M12","F12_M10","F12_M11","F12_M12","F13_M13","F13_M14","F13_M15","F14_M14","F14_M15","F15_M13","F15_M14","F15_M15","F01_M01","F01_M02","F01_M03","F02_M01","F02_M02","F02_M03","F03_M01","F03_M02","F06_M04","F06_M05","F06_M06","F07_M07","F07_M08","F07_M09","F08_M07","F08_M08","F08_M09","F09_M07","F09_M08","F09_M09")

for (x in 1:37) {
  fam = fam24[x]
  for (row in 1:nrow(G_24)) {
    pos <- OBS_REF_24[row, "POS"]
    line <- which(expect_ref_AF$POS == pos)
    #print(row)
    #print(pos)
    #print(line)
    if ((expect_ref_AF[line,fam] < 1) & (expect_ref_AF[line,fam] > 0) & (TOTAL_24[row,fam]*expect_ref_AF[line,fam] >= 5) & (TOTAL_24[row,fam]*(1-expect_ref_AF[line,fam]) >= 5)) {
      test <-G.test(x=c(OBS_REF_24[row,fam], OBS_ALT_24[row,fam]), p=c(expect_ref_AF[line,fam],(1-expect_ref_AF[line,fam])))
      print(as.numeric(test$statistic))
      if(as.numeric(test$statistic) != "Inf") { 
        G_24[row,fam] = as.numeric(test$statistic)
      }
      
    }
  }
}

# Have to remove families that exist only in one salinity
# These are: 
# In 14: F14_M13, F03_M03
# In 24: F13_M15, F01_M02, F01_M03, F06_M05

G14sub = subset(G_14, select = -c(F14_M13, F03_M03))
G24sub = subset(G_24, select = -c(F13_M15, F01_M02, F01_M03, F06_M05))

#Now, calculating the difference in G value due to salinity (deltaG):
# A positive deltaG means that the deviation from Mendelian inheritance is stronger in 14 than in 24.

deltaG = G14sub
deltaG[,5:37] = G14sub[,5:37]-G24sub[,5:37]

#Subsetting the outlier regions:

deltaG_LG1 <- deltaG[deltaG$CHROM == "LR761634.1", ]
deltaG_LG3 <- deltaG[deltaG$CHROM == "LR761636.1", ]
deltaG_LG5 <- deltaG[deltaG$CHROM == "LR761638.1", ]
deltaG_LG6 <- deltaG[deltaG$CHROM == "LR761639.1", ]
deltaG_LG7 <- deltaG[deltaG$CHROM == "LR761640.1", ]
deltaG_LG9 <- deltaG[deltaG$CHROM == "LR761642.1", ]

plot(F09_M09~POS, data=deltaG_LG6)

# Now, calculating the average deltaG across all families and also bootstrapping the 95 % confidence intervals.

avgdeltaG=data.frame(matrix(ncol = 6, nrow = nrow(deltaG)))
colnames(avgdeltaG)<-c("CHROM","POS","meandeltaG","lower95","upper95","track_no")
avgdeltaG[,1:2] = deltaG[,1:2]
avgdeltaG[,6] <- c(1:nrow(deltaG))

library(simpleboot)
library(boot)

for(x in 1:nrow(deltaG)) {
  nasum <- rowSums(!is.na(deltaG[x,5:37]))
  if (nasum >= 6) {
    avgdeltaG[x,3] = mean(as.numeric(deltaG[x,5:37]), na.rm = TRUE)
    bootstrap <- one.boot(as.numeric(as.vector(deltaG[x,5:37])), mean, trim = 0, R=10000, na.rm = TRUE)
    CI <- boot.ci(bootstrap, type = "perc")
    avgdeltaG[x,4] = CI$percent[4]
    avgdeltaG[x,5] = CI$percent[5]
  }
}

#Saving the data frame to disk:
write.csv(avgdeltaG, file="outlier_avgdeltaG.csv")
save(avgdeltaG, file="outlier_avgdeltaG.Rdata")

#For next time, can skip the bootstrap and just do
# load(file="outlier_avgdeltaG.Rdata")

#Subsetting the outlier regions:

avgdeltaG_LG1 <- avgdeltaG[avgdeltaG$CHROM == "LR761634.1", ]
avgdeltaG_LG3 <- avgdeltaG[avgdeltaG$CHROM == "LR761636.1", ]
avgdeltaG_LG5 <- avgdeltaG[avgdeltaG$CHROM == "LR761638.1", ]
avgdeltaG_LG6 <- avgdeltaG[avgdeltaG$CHROM == "LR761639.1", ]
avgdeltaG_LG7 <- avgdeltaG[avgdeltaG$CHROM == "LR761640.1", ]
avgdeltaG_LG9 <- avgdeltaG[avgdeltaG$CHROM == "LR761642.1", ]


# Finally, plotting:

# All LGs together:
library(ggplot2)
ggplot(avgdeltaG,aes(x=track_no, y=meandeltaG, col=as.factor(CHROM))) +
  labs(y="delta G", colour="Linkage Group") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point() +
  theme(panel.background = element_blank(), 
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

# One LG at a time:

library(ggpubr)

LG1 <- ggplot(avgdeltaG_LG1,aes(x=POS, y=meandeltaG)) +
  xlim(0,100000) +
  labs(y="delta G", x="Position", title="Linkage Group 1") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

LG3 <- ggplot(avgdeltaG_LG3,aes(x=POS, y=meandeltaG)) +
  xlim(32800000,32900000) +
  labs(y="delta G", x="Position", title="Linkage Group 3") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

LG5 <- ggplot(avgdeltaG_LG5,aes(x=POS, y=meandeltaG)) +
  xlim(68800000,68900000) +
  labs(y="delta G", x="Position", title="Linkage Group 5") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

LG6 <- ggplot(avgdeltaG_LG6,aes(x=POS, y=meandeltaG)) +
  xlim(38000000,38100000) +
  labs(y="delta G", x="Position", title="Linkage Group 6") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

LG7 <- ggplot(avgdeltaG_LG7,aes(x=POS, y=meandeltaG)) +
  xlim(3000000,3100000) +
  labs(y="delta G", x="Position", title="Linkage Group 7") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

LG9 <- ggplot(avgdeltaG_LG9,aes(x=POS, y=meandeltaG)) +
  xlim(36900000,37000000) +
  labs(y="delta G", x="Position", title="Linkage Group 9") +
  scale_color_brewer(palette="Paired") +
  geom_errorbar(aes(ymin=lower95, ymax=upper95), size=0.5) +
  geom_point()

figure <- ggarrange(LG1, LG3, LG5, LG6, LG7, LG9,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 2, nrow = 3)
figure

