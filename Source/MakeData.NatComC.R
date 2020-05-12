# Load and transform data for analysis
# Data: Prometheus / Convalescent positives
# This script was written to analyse the COVAM multi-viral data with both IgA and IgG
# Date: January 2020
# Authors: Rafael Joshua

#load libraries
library(tidyverse)
library(matrixStats)
library(preprocessCore)
# Set the working directory.
setwd("Y:/Rafael/COVAM/NatCom")


###################################################################################
############################### SETUP #############################################
ROOT <- getwd()
InputFiles_Dir <- "InputFolder"

# ##############################################################################
# #..............................Begining of the code .........................#
# ##############################################################################

inputData  <- read.csv(paste(InputFiles_Dir,"/natcom.wide.csv", sep = ""), stringsAsFactors = FALSE, check.names = FALSE)

pheno.df <- read.csv(paste(InputFiles_Dir,"/NatCom.pheno.csv", sep = ""), stringsAsFactors = FALSE, check.names = FALSE)
ann.df <- read.csv(paste(InputFiles_Dir,"/ann.csv", sep = ""), stringsAsFactors = FALSE, check.names = FALSE)
ann.avr <- read.csv(paste(InputFiles_Dir,"/ann.targets.avr.csv", sep = ""), stringsAsFactors = FALSE, check.names = FALSE)
order.df <- read.csv(paste(InputFiles_Dir,"/Order.csv", sep = ""), stringsAsFactors = FALSE, check.names = FALSE)


colnames(inputData) <- gsub(" ","",colnames(inputData))

rownames(inputData) <- ann.df$Unique.ID

pheno.df <- subset(pheno.df, Infection == "Positive" | Infection == "Negative")
pheno.df <- subset(pheno.df, Isotype == "IgG" | Isotype == "IgA")
pheno.df$Unique.Sample.ID <- gsub(" ", "", pheno.df$Unique.Sample.ID)

pheno.df <- subset(pheno.df, Status == "OK")

rownames(ann.df) <- ann.df$Unique.ID
rownames(order.df) <- order.df$Unique.ID
ann.df <- ann.df[order.df$Unique.ID,]

natcom.wide <- inputData[, pheno.df$Unique.Sample.ID]

natcom.pheno.IgG <- subset(pheno.df, Isotype == "IgG")
natcom.pheno.IgA <- subset(pheno.df, Isotype == "IgA")



zero.wide.df <- natcom.wide

for(s in colnames(natcom.wide)){
  for(a in rownames(natcom.wide)){
    if(natcom.wide[a,s] < 0){
      zero.wide.df[a,s] <- 0
    }
  }
}



#rownames(wide.df) <- data.df$Unique.ID
zero.wide.IgG <- zero.wide.df[ann.df$Unique.ID[which(ann.df$Spot.Type=="Antigen")], natcom.pheno.IgG$Unique.Sample.ID]
zero.wide.IgA <- zero.wide.df[ann.df$Unique.ID[which(ann.df$Spot.Type=="Antigen")], natcom.pheno.IgA$Unique.Sample.ID]


mock.IgG.vec <- rowMedians(as.matrix(zero.wide.IgG))
mock.IgA.vec <- rowMedians(as.matrix(zero.wide.IgA))


norm.wide.IgG <- as.data.frame(normalize.quantiles.use.target(as.matrix(zero.wide.IgG), mock.IgG.vec, copy=TRUE, subset=NULL))
rownames(norm.wide.IgG) <- rownames(zero.wide.IgG)
colnames(norm.wide.IgG) <- colnames(zero.wide.IgG)
norm.wide.IgG$Unique.ID <- rownames(norm.wide.IgG)

norm.wide.IgA <- as.data.frame(normalize.quantiles.use.target(as.matrix(zero.wide.IgA), mock.IgA.vec, copy=TRUE, subset=NULL))
rownames(norm.wide.IgA) <- rownames(zero.wide.IgA)
colnames(norm.wide.IgA) <- colnames(zero.wide.IgA)
norm.wide.IgA$Unique.ID <- rownames(norm.wide.IgA)


zero.wide.df$Unique.ID <- rownames(zero.wide.df)

norm.wide <- cbind(norm.wide.IgG, norm.wide.IgA)
norm.wide <- norm.wide[, pheno.df$Unique.Sample.ID]
norm.wide$Unique.ID <- rownames(norm.wide)



long.df <- gather(zero.wide.df, key = "Unique.Sample.ID", value = "MFI", -Unique.ID)
long.norm <- gather(norm.wide, key = "Unique.Sample.ID", value = "Norm.MFI", -Unique.ID)
long.df <- left_join(long.df, long.norm, by = c("Unique.Sample.ID", "Unique.ID"))
long.df <- left_join(long.df, pheno.df, by = c("Unique.Sample.ID"))
long.df <- left_join(long.df, ann.df, by = c("Unique.ID"))

rm(long.norm, mock.IgG.vec, mock.IgA.vec)

#save.image("SavedSessions/MakeData.R")


####################################################################################
############### Average
####################################################################################

targets.ann.df <- subset(ann.df, Spot.Type == "Antigen")

targets.wide.avr <- data.frame(matrix(nrow = 0, ncol = ncol(natcom.wide)))
colnames(targets.wide.avr) <- colnames(natcom.wide)



for(s in names(zero.wide.df)){
  if(s != "Unique.ID"){
    for(a in unique(targets.ann.df$ID)){
      average <- mean(zero.wide.df[ann.df$Unique.ID[which(ann.df$ID == a)], s])
      targets.wide.avr[a,s] <- average
    }
  }
}

rownames(ann.avr) <- ann.avr$ID
ann.avr <- subset(ann.avr, Keep =="Yes")

targets.wide.avr <- targets.wide.avr[ann.avr$ID,]
IgG.wide.avr <- targets.wide.avr[, pheno.df$Unique.Sample.ID[which(pheno.df$Isotype == "IgG")]]
IgA.wide.avr <- targets.wide.avr[, pheno.df$Unique.Sample.ID[which(pheno.df$Isotype == "IgA")]]

mock.IgG.vec <- rowMedians(as.matrix(IgG.wide.avr))
mock.IgA.vec <- rowMedians(as.matrix(IgA.wide.avr))


norm.IgG.avr <- as.data.frame(normalize.quantiles.use.target(as.matrix(IgG.wide.avr), mock.IgG.vec, copy=TRUE, subset=NULL))
rownames(norm.IgG.avr) <- rownames(IgG.wide.avr)
colnames(norm.IgG.avr) <- colnames(IgG.wide.avr)
norm.IgG.avr$Unique.ID <- rownames(norm.IgG.avr)

norm.IgA.avr <- as.data.frame(normalize.quantiles.use.target(as.matrix(IgA.wide.avr), mock.IgA.vec, copy=TRUE, subset=NULL))
rownames(norm.IgA.avr) <- rownames(IgA.wide.avr)
colnames(norm.IgA.avr) <- colnames(IgA.wide.avr)



# Make long data frames 



wide.avr <- cbind(IgG.wide.avr, IgA.wide.avr)
wide.avr <- wide.avr[, pheno.df$Unique.Sample.ID]
wide.avr$Unique.ID <- rownames(wide.avr)

norm.wide.avr <- cbind(norm.IgG.avr, norm.IgA.avr)
norm.wide.avr <- norm.wide.avr[, pheno.df$Unique.Sample.ID]
norm.wide.avr$Unique.ID <- rownames(norm.wide.avr)

long.avr <- gather(wide.avr, key = "Unique.Sample.ID", value = "MFI", -Unique.ID)
long.norm.avr <- gather(norm.wide.avr, key = "Unique.Sample.ID", value = "Norm.MFI", -Unique.ID)
long.avr <- left_join(long.avr, long.norm.avr, by = c("Unique.Sample.ID", "Unique.ID"))
long.avr <- left_join(long.avr, pheno.df, by = c("Unique.Sample.ID"))
long.avr <- left_join(long.avr, ann.avr, by = c("Unique.ID"))

rm(targets.ann.df, long.norm.avr, targets.wide.avr, data.df.IgA, data.df.IgG, mock.IgA.vec, mock.IgG.vec)


save.image("SavedSessions/MakeData.Nat.Com.R")
