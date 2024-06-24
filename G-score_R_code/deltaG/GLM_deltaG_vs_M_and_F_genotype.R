rm(list=ls())
setwd("~/Desktop/oysterdata/linear_model")

##    Pierre's DeltaG analysis

##    Load Libraries

library(lme4)
library(lmerTest)

##    Load data:

PDW<- read.table ("LG6.txt", sep="\t", header=TRUE)
str(PDW)

PDW$M<-as.character(PDW$M)
PDW$FGEN<-as.character(PDW$FGEN)
PDW$MGEN<-as.character(PDW$MGEN)
PDW$FHET<-as.character(PDW$FHET)
PDW$MHET<-as.character(PDW$MHET)


## LME (lmer) of identity (male, female), parental genotype (Mgen Fgen) on DeltaG 

lmerDeltaG <- lmer(deltaG ~ matrix + (1 | M) + (1 | F) + MHET + FHET, data = PDW)
summary(lmerDeltaG)

plot(resid(lmerDeltaG)~lmerDeltaG@frame$MHET)

## and if no effect of random factors above (male & female identities):

lmDeltaG <- lm(deltaG ~ MHET + FHET, data = PDW)
summary(lmDeltaG)

lmDeltaG <- lm(deltaG ~ M + F + MHET + FHET, data = PDW)
summary(lmDeltaG)
library(car)
Anova(lmDeltaG)





