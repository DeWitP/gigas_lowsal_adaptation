#Clearing env
rm(list = ls())
# Loading in the data:
setwd("~/Downloads/Dynamo_ms/bioinformatics/data/embryos24")
familylist <- c("F10_M10","F10_M11","F10_M12","F11_M10","F11_M11","F11_M12","F12_M10","F12_M11","F12_M12","F13_M13","F13_M14","F13_M15","F14_M14","F14_M15","F15_M13","F15_M14","F15_M15","F01_M01","F01_M02","F01_M03","F02_M01","F02_M02","F02_M03","F03_M01","F03_M02","F06_M04","F06_M05","F06_M06","F07_M07","F07_M08","F07_M09","F08_M07","F08_M08","F08_M09","F09_M07","F09_M08","F09_M09")

for (X in 1:37) {
  TG<-read.csv(paste("LR761634/", familylist[X], ".TotGmat.csv", sep=""))
  HG<-read.csv(paste("LR761634/", familylist[X], ".HetGmat.csv", sep=""))
  df<-read.csv(paste("LR761634/", familylist[X], ".df.csv", sep=""))
  if (X==1) {
    TotGmat <- data.frame(matrix(ncol = 37, nrow = nrow(TG)))
    HetGmat <- data.frame(matrix(ncol = 37, nrow = nrow(TG)))
    dfmat <- data.frame(matrix(ncol = 37, nrow = nrow(TG)))
    colnames(TotGmat) <- familylist
    colnames(HetGmat) <- familylist
    colnames(dfmat) <- familylist
    
  }
  TotGmat[,familylist[X]] <- TG[,familylist[X]]
  HetGmat[,familylist[X]] <- HG[,familylist[X]]
  dfmat[,familylist[X]] <- df[,familylist[X]]
}
write.csv(TotGmat, file="LR761634.TotGmat.csv")
write.csv(HetGmat, file="LR761634.HetGmat.csv")
write.csv(dfmat, file="LR761634.df.csv")
