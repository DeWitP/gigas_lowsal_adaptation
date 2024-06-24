#!/usr/bin/env Rscript
# rm(list = ls())
setwd("/scratch/05301/dewitp/dynamo/R/LR761637")
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("One argument must be supplied (family).n", call.=FALSE)
}

# Loading in the data:
#Expected proportions (calcaulated in the "parents_LR761634.1.R" script):
load("expected_props_LR761637.1.Rdata")

#Reading in the observed counts (where "NA" has been replaced with zeroes):
OBS_REF <- read.delim("embryos14_LR761637.1_REF.txt", header=TRUE)
OBS_ALT <- read.delim("embryos14_LR761637.1_ALT.txt", header=TRUE)

# dim(OBS_REF)
#Creating a new dataframe to calculate the total observed counts (needed to filter based on expected counts >= 5 below - a prerequisite for the G test to work well)
TOTAL <- data.frame(matrix(ncol = 39, nrow = 1728960))
#provide column names
colnames(TOTAL) <- colnames(OBS_REF)
TOTAL$CHROM <- OBS_REF$CHROM
TOTAL$POS <- OBS_REF$POS
TOTAL$REF <- OBS_REF$REF
TOTAL$ALT <- OBS_REF$ALT
TOTAL[,5:39] = OBS_REF[,5:39] + OBS_ALT[,5:39]
# head(TOTAL)

# G-tests for 100 KB windows:

# max(OBS_REF$POS)
## Max position is: 53118008 (53.11 MB)

# This means that the G-matrix should have 532 100-KB blocks, to cover the whole chromosome.
# Creating 3 empty dataframes, to fill with 1. the Total G values; 2. The heterogeneity G-values (not really useful though as the expected values differ across loci); 3. The degrees of freddom (Nloci) in each block.
Total_G_matrix <- data.frame(matrix(ncol = 39, nrow = 532))
Het_G_matrix <- data.frame(matrix(ncol = 39, nrow = 532))
df_matrix <- data.frame(matrix(ncol = 39, nrow = 532))

colnames(Total_G_matrix) <- colnames(OBS_REF)
colnames(Het_G_matrix) <- colnames(OBS_REF)
colnames(df_matrix) <- colnames(OBS_REF)

# Need this package for the G tests:
library(RVAideMemoire)

# Defining a function to calculate the additive (Total) G value for a part of the chromosome. Syntax: Fun.G(minPOS,maxPOS,"FAMILY",Block#)
# Outputs the sum of the G values, the heterogeneity G-value and the degrees of freedom (nloci) in the interval, in the three separate matrices (Total_G_matrix, Het_G_matrix and df_matrix).
# Also outputs the values to the screen after each iteration.

Fun.G = function(X,Y,Z,W) {
require(RVAideMemoire)

sumtest = 0
df = 0
het_test = "NA"
data <- matrix(ncol = 2, nrow = 100000)
startrow <- Position(function(A) A > X, OBS_REF[,"POS"])
endrow <- Position(function(A) A > Y, OBS_REF[,"POS"])

for (row in startrow:(endrow-1)) {
    pos <- OBS_REF[row, "POS"]
    line <- which(expect_ref_AF$POS == pos)
    #print(row)
    #print(pos)
    #print(line)
    if ((expect_ref_AF[line,Z] < 1) & (expect_ref_AF[line,Z] > 0) & (TOTAL[row,Z]*expect_ref_AF[line,Z] >= 5) & (TOTAL[row,Z]*(1-expect_ref_AF[line,Z]) >= 5)) {
        test <-G.test(x=c(OBS_REF[row,Z], OBS_ALT[row,Z]), p=c(expect_ref_AF[line,Z],(1-expect_ref_AF[line,Z])))
        # print(as.numeric(test$statistic))
    if(as.numeric(test$statistic) != "Inf") { 
    sumtest = sumtest + as.numeric(test$statistic)
      }
    df = df+1
    data[df,1] <- OBS_REF[row,Z]
    data[df,2] <- OBS_ALT[row,Z]
    }
}

hettest_data <- subset(data, data[,1] != "NA")
if(nrow(hettest_data) >1) {
   het_test <- G.test(hettest_data)
   print("Heterogeneity G-value:")
   print(as.numeric(het_test$statistic))
   Het_G_matrix[W,Z] <<- as.numeric(het_test$statistic)
} else {
  print("Heterogeneity G-value not applicable")
  Het_G_matrix[W,Z] <<- "NA"
}
print("Total G-value:")
print(sumtest)
print("Df:")
print(df)
print("Family:")
print(Z)
print("G matrix block number:")
print(W)
Total_G_matrix[W,Z] <<- sumtest
df_matrix[W,Z] <<- df
}

# Now, looping the function over all 100 KB blocks in the chromosome:
# Creating a new function, which loops over all blocks in the G matrix:
# Usage: Family_loop("F01_M01") 


Family_loop = function(X) {
  for (n in 1:nrow(Total_G_matrix)) {
  endPOS = n * 100000
  startPOS = endPOS - 99999 
  Fun.G(startPOS,endPOS,X,n)
  }
write.csv(Total_G_matrix,file=paste(X, ".TotGmat.csv", sep=""))
write.csv(Het_G_matrix,file=paste(X, ".HetGmat.csv", sep=""))
write.csv(df_matrix,file=paste(X, ".df.csv", sep=""))
}

#Now, running this for the family specified in the command-line argument:
Family_loop(args)
