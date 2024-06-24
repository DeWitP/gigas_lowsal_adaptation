
# Going to try to calculate expected allele frequencies in the different crosses. 

# As the dataset is so big, it makes sens to do it for each chromosome separately (n=10), and then also one round for the satellite chromosomes.

data<-read.delim("parents_LR761641.1.txt", header=TRUE)
dim(data)
#rn<-read.delim("rownames_parentsLR761634.1.txt", header=FALSE)

#Creating an empty dataframe to fill with the cross
expect_ref_AF <- data.frame(matrix(ncol = 47, nrow = 1762985))
#provide column names
colnames(expect_ref_AF) <- c('CHROM','POS','F01_M01','F01_M02','F01_M03','F02_M01','F02_M02','F02_M03','F03_M01','F03_M02','F03_M03','F04_M04','F04_M05','F04_M06','F05_M04','F05_M05','F05_M06','F06_M04','F06_M05','F06_M06','F07_M07','F07_M08','F07_M09','F08_M07','F08_M08','F08_M09','F09_M07','F09_M08','F09_M09','F10_M10','F10_M11','F10_M12','F11_M10','F11_M11','F11_M12','F12_M10','F12_M11','F12_M12','F13_M13','F13_M14','F13_M15','F14_M13','F14_M14','F14_M15','F15_M13','F15_M14','F15_M15')
#rownames(expect_ref_AF) <- rn
head(expect_ref_AF)

#Starting with assigning position data form the parent data file:
expect_ref_AF$CHROM=data$CHROM
expect_ref_AF$POS=data$POS

# Now, calculating the expected proportion of reference reads in the crosses given the parental genotypes at all loci:
# Starting with the first cross (F01 * M01)
expect_ref_AF$F01_M01 = 1-((data$F01 + data$M01)/4)

# Sanity check, top 20 loci:
head(expect_ref_AF$F01_M01, n=20)
head(data$F01, n=20)
head(data$M01, n=20)
# Looks correct, proceeding with the rest of the crosses.

expect_ref_AF$F01_M02 = 1-((data$F01 + data$M02)/4)
expect_ref_AF$F01_M03 = 1-((data$F01 + data$M03)/4)
expect_ref_AF$F02_M01 = 1-((data$F02 + data$M01)/4)
expect_ref_AF$F02_M02 = 1-((data$F02 + data$M02)/4)
expect_ref_AF$F02_M03 = 1-((data$F02 + data$M03)/4)
expect_ref_AF$F03_M01 = 1-((data$F03 + data$M01)/4)
expect_ref_AF$F03_M02 = 1-((data$F03 + data$M02)/4)
expect_ref_AF$F03_M03 = 1-((data$F03 + data$M03)/4)
expect_ref_AF$F04_M04 = 1-((data$F04 + data$M04)/4)
expect_ref_AF$F04_M05 = 1-((data$F04 + data$M05)/4)
expect_ref_AF$F04_M06 = 1-((data$F04 + data$M06)/4)
expect_ref_AF$F05_M04 = 1-((data$F05 + data$M04)/4)
expect_ref_AF$F05_M05 = 1-((data$F05 + data$M05)/4)
expect_ref_AF$F05_M06 = 1-((data$F05 + data$M06)/4)
expect_ref_AF$F06_M04 = 1-((data$F06 + data$M04)/4)
expect_ref_AF$F06_M05 = 1-((data$F06 + data$M05)/4)
expect_ref_AF$F06_M06 = 1-((data$F06 + data$M06)/4)
expect_ref_AF$F07_M07 = 1-((data$F07 + data$M07)/4)
expect_ref_AF$F07_M08 = 1-((data$F07 + data$M08)/4)
expect_ref_AF$F07_M09 = 1-((data$F07 + data$M09)/4)
expect_ref_AF$F08_M07 = 1-((data$F08 + data$M07)/4)
expect_ref_AF$F08_M08 = 1-((data$F08 + data$M08)/4)
expect_ref_AF$F08_M09 = 1-((data$F08 + data$M09)/4)
expect_ref_AF$F09_M07 = 1-((data$F09 + data$M07)/4)
expect_ref_AF$F09_M08 = 1-((data$F09 + data$M08)/4)
expect_ref_AF$F09_M09 = 1-((data$F09 + data$M09)/4)
expect_ref_AF$F10_M10 = 1-((data$F10 + data$M10)/4)
expect_ref_AF$F10_M11 = 1-((data$F10 + data$M11)/4)
expect_ref_AF$F10_M12 = 1-((data$F10 + data$M12)/4)
expect_ref_AF$F11_M10 = 1-((data$F11 + data$M10)/4)
expect_ref_AF$F11_M11 = 1-((data$F11 + data$M11)/4)
expect_ref_AF$F11_M12 = 1-((data$F11 + data$M12)/4)
expect_ref_AF$F12_M10 = 1-((data$F12 + data$M10)/4)
expect_ref_AF$F12_M11 = 1-((data$F12 + data$M11)/4)
expect_ref_AF$F12_M12 = 1-((data$F12 + data$M12)/4)
expect_ref_AF$F13_M13 = 1-((data$F13 + data$M13)/4)
expect_ref_AF$F13_M14 = 1-((data$F13 + data$M14)/4)
expect_ref_AF$F13_M15 = 1-((data$F13 + data$M15)/4)
expect_ref_AF$F14_M13 = 1-((data$F14 + data$M13)/4)
expect_ref_AF$F14_M14 = 1-((data$F14 + data$M14)/4)
expect_ref_AF$F14_M15 = 1-((data$F14 + data$M15)/4)
expect_ref_AF$F15_M13 = 1-((data$F15 + data$M13)/4)
expect_ref_AF$F15_M14 = 1-((data$F15 + data$M14)/4)
expect_ref_AF$F15_M15 = 1-((data$F15 + data$M15)/4)

# Check that it all worked out:
head(expect_ref_AF)
# Looks good! Saving the new data frame:
save(expect_ref_AF, file="expected_props_LR761641.1.Rdata")




