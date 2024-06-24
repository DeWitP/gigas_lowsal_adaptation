#!/usr/bin/env Rscript
# rm(list = ls())
setwd("/scratch/05301/dewitp/dynamo/R/embryos24/CADCXH")
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("One argument must be supplied (family).n", call.=FALSE)
}

# Loading in the data:
#Expected proportions (calcaulated in the "parents_LR761634.1.R" script):
load("expected_props_CADCXH.Rdata")

#Reading in the observed counts (where "NA" has been replaced with zeroes):
OBS_REF <- read.delim("embryos24_CADCXH_REF.txt", header=TRUE)
OBS_ALT <- read.delim("embryos24_CADCXH_ALT.txt", header=TRUE)

# dim(OBS_REF)
#Creating a new dataframe to calculate the total observed counts (needed to filter based on expected counts >= 5 below - a prerequisite for the G test to work well)
TOTAL <- data.frame(matrix(ncol = 41, nrow = 748486))
#provide column names
colnames(TOTAL) <- colnames(OBS_REF)
TOTAL$CHROM <- OBS_REF$CHROM
TOTAL$POS <- OBS_REF$POS
TOTAL$REF <- OBS_REF$REF
TOTAL$ALT <- OBS_REF$ALT
TOTAL[,5:41] = OBS_REF[,5:41] + OBS_ALT[,5:41]
# head(TOTAL)

# G-tests for the satellite contigs:
# I'll do just one G score for each contig, I think.

CADCXHlist <- read.delim("CADCXH_contiglist.txt", header=TRUE)
# dim(CADCXHlist)
# 226 contigs.
# This means that the G-matrix should have 226 blocks, to cover the whole chromosome.
# Creating 3 empty dataframes, to fill with 1. the Total G values; 2. The heterogeneity G-values (not really useful though as the expected values differ across loci); 3. The degrees of freddom (Nloci) in each block.
Total_G_matrix <- data.frame(matrix(ncol = 41, nrow = 226))
Het_G_matrix <- data.frame(matrix(ncol = 41, nrow = 226))
df_matrix <- data.frame(matrix(ncol = 41, nrow = 226))

colnames(Total_G_matrix) <- colnames(OBS_REF)
colnames(Het_G_matrix) <- colnames(OBS_REF)
colnames(df_matrix) <- colnames(OBS_REF)

# Need this package for the G tests:
library(RVAideMemoire)

# Defining a function to calculate the additive (Total) G value for a user-defined chromosome. Syntax: Fun.G("contigname","FAMILY",Block#)
# Outputs the sum of the G values, the heterogeneity G-value and the degrees of freedom (nloci) in the interval, in the three separate matrices (Total_G_matrix, Het_G_matrix and df_matrix).
# Also outputs the values to the screen after each iteration.

Fun.G = function(X,Z,W) {
require(RVAideMemoire)
sumtest = 0
df = 0
het_test = "NA"
data <- matrix(ncol = 2, nrow = max(CADCXHlist[,"contiglength"]))
#startrow <- Position(function(A) A > X, OBS_REF[,"POS"])
#endrow <- Position(function(A) A > Y, OBS_REF[,"POS"])
contigOBSREF = subset.data.frame(OBS_REF, CHROM == X)
contigOBSALT = subset.data.frame(OBS_ALT, CHROM == X)
contigEXP = subset.data.frame(expect_ref_AF, CHROM == X)
contigTOTAL = subset.data.frame(TOTAL, CHROM == X)
for (row in 1:nrow(contigOBSREF)) {
    pos <- contigOBSREF[row, "POS"]
    line <- which(contigEXP$POS == pos)
    #print(row)
    #print(pos)
    #print(line)
    if ((contigEXP[line,Z] < 1) & (contigEXP[line,Z] > 0) & (contigTOTAL[row,Z]*contigEXP[line,Z] >= 5) & (contigTOTAL[row,Z]*(1-contigEXP[line,Z]) >= 5)) {
        test <-G.test(x=c(contigOBSREF[row,Z], contigOBSALT[row,Z]), p=c(contigEXP[line,Z],(1-contigEXP[line,Z])))
        # print(as.numeric(test$statistic))
    if(as.numeric(test$statistic) != "Inf") { 
    sumtest = sumtest + as.numeric(test$statistic)
      }
    df = df+1
    data[df,1] <- contigOBSREF[row,Z]
    data[df,2] <- contigOBSALT[row,Z]
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

# Test:
# Fun.G("CADCXH010000001.1","F01_M01",1)
# Seems to work!

# Now, looping the function over all 226 contigs:
# Creating a new function, which loops over all contigs in the G matrix:
# Usage: Family_loop("F01_M01") 


Family_loop = function(X) {
  for (n in 1:nrow(Total_G_matrix)) {
  contigNAME = CADCXHlist[n,"contigname"]
  nloci = nrow(subset.data.frame(OBS_REF, CHROM == contigNAME))
  if(nloci > 0) {
    Fun.G(contigNAME,X,n)
  }
}
write.csv(Total_G_matrix,file=paste(X, ".TotGmat.csv", sep=""))
write.csv(Het_G_matrix,file=paste(X, ".HetGmat.csv", sep=""))
write.csv(df_matrix,file=paste(X, ".df.csv", sep=""))
}

#Now, running this for the family specified in the command-line argument:
Family_loop(args)
