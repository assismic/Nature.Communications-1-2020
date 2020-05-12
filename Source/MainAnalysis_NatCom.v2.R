setwd("Y:/Rafael/COVAM/NatCom")
load("SavedSessions/MakeData.Nat.Com.R")
#source("Y:/Rafael/Scripts/R Scripts/SupporFunctions/SupportFunctions.R")


library(sjmisc)
library(caret)
library(ROCR)
library(matrixStats)
library(pROC)
library(ggpubr)
library(ggplot2)


#################################################################################################
# Pos vcs nown negatives (prometheus)
  #Generate Excel dataframes
  outputFolder <- "Output"
  FiguresFolder <- 'Figures'
  
  long.avr$Unique.ID <- factor(long.avr$Unique.ID, levels = c(ann.avr$Unique.ID))
  
  
        ################################################################################################  
            # Step1
        ################################################################################################  
  #Box PLots
  #full array
  for(i in unique(long.avr$Isotype)){
    temp.long <- subset(long.avr, Isotype == i)
    print(
      ggplot(temp.long, aes(x = Unique.ID, y = Norm.MFI, fill=Infection, color = Infection)) +
        geom_boxplot(outlier.shape = NA, width = 0.8) +
        scale_fill_manual(values=c("#0361A6","#d10a0a")) + #(Blue: #17A9D8, red)
        scale_color_manual(name = "Infection", values=c("#0361A6","#d10a0a"))+
        ylab("Normalized MFI") +
        labs(title = paste("CoV-SARS-2 array reactivity -", i, sep = " ")) +
        theme_mine() +
        theme(plot.title = element_text(color = "black", size  = 12),
              axis.text.x = element_text(face = "plain", size = 8, angle = 90, hjust=1, vjust = 0.5),
              axis.title.x = element_blank(), 
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              plot.caption = element_text(size = 13, face = "italic", color = "coral2"))
    )
     #ggsave(file.path(ROOT,FiguresFolder,paste("Boxplot_FullArray.avr.NoOutline", i, ".tiff", sep="")), width = 12, height = 6,compression = "lzw")
    rm(temp.long)
    
  }

  
  ################################################################################################
  #Step 2 - Combination AUCs
    #2.1. Bootstrap individual AUCs (Select top picks on excel)
    #2.2. combinat AUC with top picks set
    #2.3. DRAW plots
  ################################################################################################   
  
  #2.1. Bootstrap individual AUCs (Select top picks on excel)
  # library(grDevices)

  temp.df <- norm.wide.avr[ann.avr$Unique.ID[which(ann.avr$KillerTrio == "Yes")],]
  temp.df <- as.data.frame(t(subset(temp.df, select = -c(Unique.ID))))
  top.antigens.vec <- colnames(temp.df)
  #reference.name <- data.frame(Reference = rownames(norm.wide.avr), Transformed = colnames(temp.df))
  classifier <- pheno.df$Infection
  temp.df$Classifier <- classifier
  
  
  #rownames(auc.bootstrap.antigens.IgG.df) <- rownames(auc.bootstrap.antigens.IgA.df) <- top.antigens.vec
  iterations <- 1000
  
  auc.bootstrap.antigens.IgG.df <- auc.bootstrap.antigens.IgA.df <- data.frame(matrix(nrow = 0, ncol = 0))
  for(iso in unique(pheno.df$Isotype)){
    iso.sub.temp.df <- temp.df[pheno.df$Unique.Sample.ID[which(pheno.df$Isotype == iso)], ]
      
      for (p in top.antigens.vec){
        protein.sub.temp.df <- data.frame(iso.sub.temp.df[,p])
        
        Classifier <- iso.sub.temp.df$Classifier
        Classifier <- ifelse(test = Classifier == "Positive", yes = 1, no = 0) 
        protein.sub.temp.df$Classifier <- factor(Classifier)
        
        for(ite in 1:iterations){
          print(paste("Ite:", ite, ": ",p , " Iso: ", iso))
          trainIndex.index <- createDataPartition(y=protein.sub.temp.df$Classifier, p=0.75, list=FALSE)
          trainingSet <- protein.sub.temp.df[trainIndex.index,]
          testingSet <- protein.sub.temp.df[-trainIndex.index,]
          
          model <- glm(Classifier ~., family = binomial(link = 'logit'), data = trainingSet)########
          pr <- predict(model, newdata = testingSet, type='response')
          pre <- prediction(pr, testingSet$Classifier)
          prf <- performance(pre, measure = "tpr", x.measure = "fpr")
            auc <- performance(pre, measure = "auc")
          auc <- auc@y.values[[1]]
          if(iso == "IgA"){
            auc.bootstrap.antigens.IgA.df[p,ite] <- auc
            
          }
          if(iso == "IgG"){
            auc.bootstrap.antigens.IgG.df[p,ite] <- auc
            
          }
        }
      }
  }
  
  bootstrap.Ind.AUC <- data.frame(matrix(nrow = nrow(auc.bootstrap.antigens.IgA.df), ncol = 0))
  rownames(bootstrap.Ind.AUC) <- rownames(auc.bootstrap.antigens.IgA.df)

  bootstrap.Ind.AUC$IgG_1K.median <- rowMedians(as.matrix(auc.bootstrap.antigens.IgG.df))
  bootstrap.Ind.AUC$n_perfect.IgG <- row_count(auc.bootstrap.antigens.IgG.df, count = 1, append = FALSE)$rowcount*100/iterations
  
  bootstrap.Ind.AUC$IgA_1K.median <- rowMedians(as.matrix(auc.bootstrap.antigens.IgA.df))
  bootstrap.Ind.AUC$n_perfect.IgA <- row_count(auc.bootstrap.antigens.IgA.df, count = 1, append = FALSE)$rowcount*100/iterations
  rm(auc.bootstrap.antigens.IgG.df, auc.bootstrap.antigens.IgA.df)
  #write.csv(data.frame(bootstrap.Ind.AUC), "Output/bootstrap.Indv.AUC.median.1000.075div.csv")
  
  
#########################################################################################
##2.2. combinat AUC with top picks set
#########################################################################################
  #2.2. combinat AUC with top picks set
  #Get bootstrap combination AUCs reference
  
  top.antigens.vec <-  c("SARS-CoV-2_S1+S2",
                         "SARS-CoV-2_NP",
                         "SARS-CoV-2_S2",
                         "SARS-CoV_NP",
                         "SARS-CoV-2_S1, (mFc Tag)",
                         "MERS-CoV_S2",
                         "SARS-CoV-2_S1-RBD"
  )
  
  temp.df <- norm.wide.avr[top.antigens.vec,]
  temp.df <- as.data.frame(t(subset(temp.df, select = -c(Unique.ID))))
  classifier <- pheno.df$Infection
  temp.df$Classifier <- classifier

  iterations <- 1000
  auc.bootstrap.antigens.IgG.df <- auc.bootstrap.antigens.IgA.df <- data.frame(matrix(nrow = 0, ncol = 0))
for(iso in unique(pheno.df$Isotype)){
    iso.sub.temp.df <- temp.df[pheno.df$Unique.Sample.ID[which(pheno.df$Isotype == iso)], ]
  for(l in 1:2){
    proteinSubsetCombination <- combn(top.antigens.vec,l, simplify = TRUE)
    
    for (p in 1:ncol(proteinSubsetCombination)){
      protein.sub.temp.df <- as.data.frame(iso.sub.temp.df[,proteinSubsetCombination[,p]])
      
      Classifier <- iso.sub.temp.df$Classifier
      Classifier <- ifelse(test = Classifier == "Positive", yes = 1, no = 0) 
      protein.sub.temp.df$Classifier <- factor(Classifier)
      
      for(ite in 1:iterations){
        print(paste(iso,"l: ",l, " col: ", p, "/", ncol(proteinSubsetCombination), " iter: ", ite, sep = ""))
        row <- paste(proteinSubsetCombination[,p], collapse = "; ")
        
        trainIndex.index <- createDataPartition(y=protein.sub.temp.df$Classifier, p=0.75, list=FALSE)
        trainingSet <- protein.sub.temp.df[trainIndex.index,]
        testingSet <- protein.sub.temp.df[-trainIndex.index,]
        
        model <- glm(Classifier ~., family = binomial(link = 'logit'), data = trainingSet)########
        pr <- predict(model, newdata = testingSet, type='response')
        pre <- prediction(pr, testingSet$Classifier)
        prf <- performance(pre, measure = "tpr", x.measure = "fpr")
        auc <- performance(pre, measure = "auc")
        auc <- auc@y.values[[1]]
        if(iso == "IgA"){
          auc.bootstrap.antigens.IgA.df[row,"Len"] <- l
          auc.bootstrap.antigens.IgA.df[row,ite] <- auc
          
        }
        if(iso == "IgG"){
          auc.bootstrap.antigens.IgG.df[row,"Len"] <- l
          auc.bootstrap.antigens.IgG.df[row,ite] <- auc
          
        }
      }
    }
  }
}
  
  bootstrap.Combinat.AUC <- data.frame(matrix(nrow = nrow(auc.bootstrap.antigens.IgA.df), ncol = 0))
  rownames(bootstrap.Combinat.AUC) <- rownames(auc.bootstrap.antigens.IgA.df)
  
  bootstrap.Combinat.AUC$Len <- auc.bootstrap.antigens.IgG.df$Len
  bootstrap.Combinat.AUC$IgG_1K.median <- rowMedians(as.matrix(auc.bootstrap.antigens.IgG.df))
  bootstrap.Combinat.AUC$n_perfect.IgG <- row_count(auc.bootstrap.antigens.IgG.df, count = 1, append = FALSE)$rowcount*100/iterations
  
  bootstrap.Combinat.AUC$IgA_1K.median <- rowMedians(as.matrix(auc.bootstrap.antigens.IgA.df))
  bootstrap.Combinat.AUC$n_perfect.IgA <- row_count(auc.bootstrap.antigens.IgA.df, count = 1, append = FALSE)$rowcount*100/iterations
  #write.csv(data.frame(bootstrap.Combinat.AUC), "Output/bootstrap.Indv.AUC.median.1000.075div.csv")

  
#########################################################################################
##2.3. DRAW plots
#########################################################################################

  top.antigens.vec <-  c("SARS-CoV-2_S1+S2",
                         "SARS-CoV-2_NP",
                         "SARS-CoV-2_S2",
                         "SARS-CoV_NP",
                         "SARS-CoV-2_S1, (mFc Tag)",
                         "MERS-CoV_S2",
                         "SARS-CoV-2_S1-RBD"
  )
  
  colorsBrew <- c("#a60303", "#54a603","#03a6a6", "#5403a6")
  
  bootstrap.Combinat.AUC <- read.csv("Output/bootstrap.Combinat.AUC.median.1000.075div.csv", stringsAsFactors = FALSE, check.names = FALSE)
  rownames(bootstrap.Combinat.AUC) <- bootstrap.Combinat.AUC$AntigenCombinat
  
  temp.df <- norm.wide.avr[top.antigens.vec,]
  temp.df <- as.data.frame(t(subset(temp.df, select = -c(Unique.ID))))
  
  classifier <- pheno.df$Infection
  temp.df$Classifier <- classifier
  
  i1 <- i2 <- i3 <- i4 <- c()
  Youden.df <- data.frame(matrix(nrow = nrow(bootstrap.Combinat.AUC)))
  Youden.df$IgG.TH <- Youden.df$IgG.Spe <- Youden.df$IgG.Sensi <- Youden.df$IgA.TH <- Youden.df$IgA.Spe <- Youden.df$IgA.Sensi <- ""
  
  rownames(Youden.df) <- rownames(bootstrap.Combinat.AUC)
  for(iso in unique(pheno.df$Isotype)){
    iso.sub.temp.df <- temp.df[pheno.df$Unique.Sample.ID[which(pheno.df$Isotype == iso)], ]
    
    for (i in 1:4){ 
      proteinSubsetCombination <- combn(top.antigens.vec,i, simplify = TRUE)
      if(i == 1){
        title.size <- 1
      } else if(i == 2){
        title.size <- 0.75
      }else if(i == 3){
        title.size <- 0.5
      }else if(i == 4){
        title.size <- 0.25
      }
      for (p in 1:ncol(proteinSubsetCombination)){
        row <- paste(proteinSubsetCombination[,p], collapse = "; ")
        #for(row in rownames(bootstrap.Combinat.AUC)){
        if(row %in% paste(proteinSubsetCombination[,p], collapse = "; ")){
          r <- 1
          if(iso == "IgG"){
            ref.auc <- round(bootstrap.Combinat.AUC[row, "IgG_1K.median"], digits =  3)
            #print(paste(row, ref.auc))
          }
          if(iso == "IgA"){
            ref.auc <- round(bootstrap.Combinat.AUC[row, "IgA_1K.median"], digits =  3)
            #print(paste(row, ref.auc))
          }
          a.sub.temp.df <- data.frame(iso.sub.temp.df[, proteinSubsetCombination[,p]])
          Classifier <- iso.sub.temp.df$Classifier
          Classifier <- ifelse(test = Classifier == "Positive", yes = 1, no = 0)
          a.sub.temp.df$Classifier <- factor(Classifier)
          trainIndex.index <- createDataPartition(y=a.sub.temp.df$Classifier, p=0.75, list=FALSE)
          trainingSet <- a.sub.temp.df[trainIndex.index,]
          
          model <- glm(Classifier ~., family = binomial(link = 'logit'), data = trainingSet)
          auc <- round(roc(trainingSet$Classifier, model$fitted.values, plot=FALSE)$auc, digits = 3)
          singles <- proteinSubsetCombination[,p]
          print(paste("first round auc:", auc))
          print(row)
          
          if(auc == ref.auc){
            antigen.title <- row
            mainModel <- model
            print(paste("MATCH", " Round", r, "auc:", auc, "@ digits ", ndigits))
            par(pty ="s")
            #png(paste("Figures/MultiAUC/",iso,".",i,".",row, ".AUC.100.075.tiff", sep=""),width = 960, height = 960,
            #pointsize = 12, res = 300)
            
            main.roc <- roc(trainingSet$Classifier, mainModel$fitted.values, plot=TRUE, print.auc=F,print.auc.x = 0.35, print.auc.y = 0.1,
                            legacy.axes = T, col = "#0361A6", lwd = 3, main = row, print.auc.cex = 0.75, cex.main=title.size)
            roc.df <- data.frame(
              tpp=(main.roc$sensitivities)*100,
              fpp=(1-main.roc$specificities)*100,
              thresholds = ifelse(test =  main.roc$thresholds == "-Inf", yes = "Inf", no = main.roc$thresholds)
            )
            #write.csv(roc.df , paste("Output/rocDFs/",iso,".",i,".",row, ".AUC.csv", sep=""))
            aucList <- c(auc)
            elements <- c("")
            coor <- round(coords(main.roc, x="best", input="threshold", best.method="youden", transpose = TRUE),3)
            specificity <- coor[2]
            Sensitivity <- coor[3]
            Youden.df[row, paste(iso, "TH", sep = ".")] <- coor[1]
            Youden.df[row, paste(iso, "Spe", sep = ".")] <- coor[2]
            Youden.df[row, paste(iso, "Sensi", sep = ".")] <- coor[3]
            if(length(singles) > 1){
              print("HERE SINGLES")
              elements <- c("Combination")
              for(s in 1:length(singles)){
                print(s)
                temp.data <- data.frame(trainingSet[,c(colnames(trainingSet)[s], "Classifier")])
                antigen.title <- colnames(trainingSet)[s]
                
                submodel <- glm(Classifier ~., family = binomial(link = 'logit'), data = temp.data)
                sub.roc <- roc(trainingSet$Classifier, submodel$fitted.values, plot=TRUE, print.auc=F, print.auc.x = 0.35, print.auc.y = 0.1+(0.1*s),
                               legacy.axes = T, col = colorsBrew[s], lwd = 2, main = antigen.title, print.auc.cex = 0.75, add = T)
                aucList  <- c(aucList, round(sub.roc$auc, 3))
                elements <- c("Combination", singles)
                
                roc.df <- data.frame(
                  tpp=(sub.roc$sensitivities)*100,
                  fpp=(1-sub.roc$specificities)*100,
                  thresholds = ifelse(test =  sub.roc$thresholds == "-Inf", yes = "Inf", no = sub.roc$thresholds)
                )
                #write.csv(roc.df , paste("Output/rocDFs/sub.roc.",iso,".",i,".",singles[s], ".AUC.csv", sep=""))

              }
            }
            
            #dev.off()
            
          }
          
          
          
          while(auc != ref.auc){
            print(r)
            if(r < 5000){
              ndigits <- 3
            }else{
              ndigits <- 2
            }
            ref.auc <- round(ref.auc, digits = ndigits)
            #print(row)
            print(paste("auc: ", auc, "ref:", ref.auc))
            trainIndex.index <- createDataPartition(y=a.sub.temp.df$Classifier, p=0.70, list=FALSE)
            trainingSet <- a.sub.temp.df[trainIndex.index,]
            testingSet <- a.sub.temp.df[-trainIndex.index,]
            

            model <- glm(Classifier ~., family = binomial(link = 'logit'), data = trainingSet)
            mainModel <- model
            
            auc <- round(roc(trainingSet$Classifier, mainModel$fitted.values, plot=FALSE, print.auc=TRUE, print.auc.y = 0.2,
                             legacy.axes = T, col = "#0361A6", lwd = 3, main = a)$auc, digits = ndigits)
            print(paste("Round", r, "auc:", auc, "@ digits ", ndigits))
            if(auc == ref.auc){
              print(paste("MATCH", " Round", r, "auc:", auc, "@ digits ", ndigits))
              par(pty ="s")
              #png(paste("Figures/MultiAUC/",iso,".",i,".",row, ".AUC.100.075.tiff", sep=""),width = 960, height = 960,
              #pointsize = 12, res = 300)
              
              main.roc <- roc(trainingSet$Classifier, mainModel$fitted.values, plot=TRUE, print.auc=F,print.auc.x = 0.35, print.auc.y = 0.1,
                              legacy.axes = T, col = "#0361A6", lwd = 3, main = row, print.auc.cex = 0.75, cex.main=title.size)
              roc.df <- data.frame(
                tpp=(main.roc$sensitivities)*100,
                fpp=(1-main.roc$specificities)*100,
                thresholds = ifelse(test =  main.roc$thresholds == "-Inf", yes = "Inf", no = main.roc$thresholds)
              )
              #write.csv(roc.df , paste("Output/rocDFs/",iso,".",i,".",row, ".AUC.csv", sep=""))
              aucList <- c(auc)
              elements <- c("")
              coor <- round(coords(main.roc, x="best", input="threshold", best.method="youden"),3)
              specificity <- coor[2]
              Sensitivity <- coor[3]
              Youden.df[row, paste(iso, "TH", sep = ".")] <- coor[1]
              Youden.df[row, paste(iso, "Spe", sep = ".")] <- coor[2]
              Youden.df[row, paste(iso, "Sensi", sep = ".")] <- coor[3]
              if(length(singles) > 1){
                print("HERE SINGLES")
                elements <- c("Combination")
                for(s in 1:length(singles)){
                  print(s)
                  temp.data <- data.frame(trainingSet[,c(colnames(trainingSet)[s], "Classifier")])
                  antigen.title <- colnames(trainingSet)[s]
                  
                  submodel <- glm(Classifier ~., family = binomial(link = 'logit'), data = temp.data)
                  sub.roc <- roc(trainingSet$Classifier, submodel$fitted.values, plot=TRUE, print.auc=F, print.auc.x = 0.35, print.auc.y = 0.1+(0.1*s),
                                 legacy.axes = T, col = colorsBrew[s], lwd = 2, main = antigen.title, print.auc.cex = 0.75, add = T)
                  aucList  <- c(aucList, round(sub.roc$auc, 3))
                  elements <- c("Combination", singles)
                  
                  roc.df <- data.frame(
                    tpp=(sub.roc$sensitivities)*100,
                    fpp=(1-sub.roc$specificities)*100,
                    thresholds = ifelse(test =  sub.roc$thresholds == "-Inf", yes = "Inf", no = sub.roc$thresholds)
                  )
                  #write.csv(roc.df , paste("Output/rocDFs/sub.roc.",iso,".",i,".",singles[s], ".AUC.csv", sep=""))
                  
                  
                }
              }
              #dev.off()
              
            }
            r <- r + 1
          }
          if(i == 1){
            i1<- c(i1, r)
          } else if(i == 2){
            i2 <- c(i2, r)
          }else if(i == 3){
            i3 <- c(i3, r)
          }else if(i == 4){
            i4 <- c(i4, r)
          }
          
          print(paste("End auc: ", auc, "ref:", ref.auc))
          
          
        }
        
        # }
        
      }
    }
  }
  
#write.csv(Youden.df, "Output/Youden.df.csv")
  
  
################################################################################################
#Step 4 - Draw Individual Antigen Boxplots
################################################################################################   

temp.long <- subset(long.avr, Unique.ID %in% top.antigens.vec)
temp.long <- subset(temp.long, Isotype == "IgG" | Isotype == "IgA")
my_comparisons <- list(c("Positive", "Negative"))

temp.long$Isotype <- factor(temp.long$Isotype, levels = c("IgG", "IgA"))
temp.long <- subset(temp.long, Infection == "Positive" | Infection == "Negative")
  
for (i in unique(temp.long$Isotype)){
  sub.temp.long <- subset(temp.long, Isotype == i)
  for(a in unique(sub.temp.long$Unique.ID)){
    sub.3.temp.long <- subset(sub.temp.long, Unique.ID == a)
    print(a)
    print(
      ggplot(sub.3.temp.long, aes(x = Infection, y = Norm.MFI, fill=Infection))+
        geom_boxplot(outlier.shape = NA, width = 0.5) +
        geom_point(size = 1, pch = 21, position = position_jitterdodge())+
        scale_fill_manual(values=c("#0361A6","#d10a0a"))+#, "#3D835D")) + #(Blue: #17A9D8, red, #3D835D= green)
        scale_color_manual(name = "Infection", values=c("#0361A6","#d10a0a"))+#, "#3D835D"))+
        stat_compare_means(fontface = "italic", size = 4, comparisons = my_comparisons, label = "p.format",label.y = 20000 )+#, comparisons = my_comparisons, ) +
        ylab("Normalized MFI") +
        xlab("Sample group") +
        ylim(0, 21000) +
        labs(title = paste(a, i , sep = " - "))+
        theme_mine() +
        theme(strip.text = element_text(size = 12),
              strip.text.x = element_text(size = 10),
              plot.title = element_text(color = "black", size = 10, face = "bold"),
              axis.text.x = element_text(face = "plain", size = 12, angle = 90, hjust=1, vjust = 0.5),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 12),
              legend.title = element_text(size = 12), legend.text = element_text(size = 11),
              plot.caption = element_text(size = 16, face = "italic", color = "coral2"))
    )
    #ggsave(file.path(ROOT,FiguresFolder,paste("Antigens/", i,".Boxplot.", a, ".tiff", sep="")), width = 5, height = 4,compression = "lzw")
    #  
    #  
  }
}
 
print("stop")
  
  
  
################################################################################################
#Step 5 - Excel Heatmaps 
################################################################################################   

  
  write.csv(norm.IgA.avr[ann.avr$Unique.ID[which(ann.avr$Spot.Type == "Antigen")],
                         c(natcom.pheno.IgA$Unique.Sample.ID[which(natcom.pheno.IgA$Infection == "Positive")],
                           natcom.pheno.IgA$Unique.Sample.ID[which(natcom.pheno.IgA$Infection == "Negative")]
                           )], "Output/norm.IgA.avr.csv")
  
  write.csv(norm.IgG.avr[ann.avr$Unique.ID[which(ann.avr$Spot.Type == "Antigen")],
                         c(natcom.pheno.IgG$Unique.Sample.ID[which(natcom.pheno.IgG$Infection == "Positive")],
                           natcom.pheno.IgG$Unique.Sample.ID[which(natcom.pheno.IgG$Infection == "Negative")]
                         )], "Output/norm.IgG.avr.csv")

 
