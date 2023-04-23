#----------------------------------------------------------------------------------------------------------
#Installing Packages
#----------------------------------------------------------------------------------------------------------
#install.packages('ranger')
#install.packages('survival')
#install.packages('pROC')
#install.packages("randomForestSRC", repos = "https://cran.us.r-project.org")
#install.packages('survivalROC')
#install.packages('ggRandomForests')
#install.packages('survminer')
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Loading in Packages
#----------------------------------------------------------------------------------------------------------
library(ranger)
library(survival) 
library(caTools)
library(pROC)
library(randomForestSRC)
library(magrittr)
library(dplyr)
library(survivalROC)
library(ggRandomForests)
library("survival")
library("survminer")
library("Rcpp")
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Load in Dataset 07-GenomicPhenotypicData
#----------------------------------------------------------------------------------------------------------
data = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL - RSF/GenomicPhenotypicData.csv") 
data <- na.omit(data)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
noriskdata = data %>% select(-c("RISK"))
df <- data %>% select(-c("PATIENT","RISK"))
featureNames <- paste(colnames(df), collapse = " + ")
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
survival_formula <-paste('Surv(MONTHS,STATUS) ~', featureNames)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
#set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
set.seed(16)

sample <- sample.split(noriskdata$PATIENT, SplitRatio = 0.3)
train  <- subset (noriskdata%>% select(-c("PATIENT")), sample == FALSE)
test   <- subset(noriskdata %>% select(-c("PATIENT")), sample == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
out.rsf.3 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                    data = train, 
                    mtry = sqrt(196), 
                    ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                    importance=TRUE, xtest = test, nodesize = 6, nsplit = 1, seed = 3)

error <- gg_error(out.rsf.3)
lowest <- min( error$error[error$error!=min(error$error)] )
row <- error[which(error == lowest),]
plot(gg_error(out.rsf.3))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

numtree = 900 #row$ntree[1]
#RANDOM SURVIVAL FOREST
survival_modelrfsrc<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = train, mtry = sqrt(196), 
                           ntree = numtree, cv.fold = 10, verbose = TRUE, block.size=1,
                           importance=TRUE, nodesize = 6, nsplit = 1, seed = 4)
set.seed(3)
nobs <- NROW(df)
surv.fit11 <- survivalROC(train$MONTHS, train$STATUS, survival_modelrfsrc$predicted.oob, 
                           predict.time = 11, method = "KM")
surv.fit16 <- survivalROC(train$MONTHS, train$STATUS, survival_modelrfsrc$predicted.oob, 
                         predict.time = 16, method = "KM")
surv.fit20 <- survivalROC(train$MONTHS, train$STATUS, survival_modelrfsrc$predicted.oob, 
                         predict.time = 20, method = "KM")

plot(surv.fit11$FP, surv.fit11$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",numtree," Trees Genomic and Phenotypic Data"))
lines(surv.fit16$FP, surv.fit16$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP")
lines(surv.fit20$FP, surv.fit20$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP")
abline(0,1)
variable11 <- paste0("11 month AUC =", round(surv.fit11$AUC,3))
variable16 <- paste0("16 month AUC =", round(surv.fit16$AUC,3))
variable20 <- paste0("20 month AUC =", round(surv.fit20$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(variable11,variable16,variable20))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
fit <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = data)
print(fit)
summary(fit)$table
ggsurvplot(fit,
           data = data,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
pred <- predict(survival_modelrfsrc,test[, (colnames(test))],block.size=1)
plot(pred$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))

legend("topright",fill=c("orange"),c("test","OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = train$MONTHS, censoring = train$STATUS, predicted = survival_modelrfsrc$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
bs.km <- get.brier.survival(survival_modelrfsrc, cens.mode = "km")$brier.score
bs.rsf <- get.brier.survival(survival_modelrfsrc, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(bs.km, type = "s", col = 2)
lines(bs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------
#Section 2: Dataset with Clinical Data
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 09-GenomicPhenotypicClinicalData
#----------------------------------------------------------------------------------------------------------
Cdata = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL - RSF/GenomicPhenotypicClinicalData.csv") #Cdata = Clinical Data
Cdata <- na.omit(Cdata)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
Cdata$Male <- ifelse(Cdata$SEX == 'Male', 1, 0)
Cdata$Female <- ifelse(Cdata$SEX == 'Female', 1, 0)
Cdata$White <- ifelse(Cdata$RACE == 'WHITE', 1, 0)
Cdata$Black <- ifelse(Cdata$RACE == 'BLACK OR AFRICAN AMERICAN', 1, 0)
Cdata$Asian <- ifelse(Cdata$RACE == 'ASIAN', 1, 0)
Cdata$T1 <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T1', 1, 0)
Cdata$T1a <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T1a', 1, 0)
Cdata$T1b <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T1b', 1, 0)
Cdata$T2 <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T2', 1, 0)
Cdata$T2a <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T2a', 1, 0)
Cdata$T2b <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T2b', 1, 0)
Cdata$T3 <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T3', 1, 0)
Cdata$T3a <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T3a', 1, 0)
Cdata$T3b <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T3b', 1, 0)
Cdata$T3c <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T3c', 1, 0)
Cdata$T4 <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'T4', 1, 0)
Cdata$TX <- ifelse(Cdata$AJCC_TUMOR_PATHOLOGIC_PT == 'TX', 1, 0)
CDatawD <- Cdata %>% select(-c("SEX","RACE","AJCC_TUMOR_PATHOLOGIC_PT")) #Clinical data with dummy variables
CDataR <- CDatawD %>% select(-c("RISK")) #Clinical Data no Risk
Cdf <- CDatawD %>% select(-c("PATIENT","RISK")) #Clinical Dataframe
CfeatureNames <- paste(colnames(Cdf), collapse = " + ") #Clinical Feature Names
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Csurvival_formula <-paste('Surv(MONTHS,STATUS) ~', CfeatureNames) #Clinical Survival Formula
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
#set.seed(12)
#set.seed(13)
set.seed(14) #run later
#set.seed(15)
#set.seed(16)

sample <- sample.split(CDataR$PATIENT, SplitRatio = 0.3)
Ctrain  <- subset (CDataR%>% select(-c("PATIENT")), sample == FALSE) #clinical Train
Ctest   <- subset(CDataR %>% select(-c("PATIENT")), sample == TRUE) #Clinical Test
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
Cout.rsf.3 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                    data = Ctrain, 
                    mtry = sqrt(196), 
                    ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                    importance=TRUE, xtest = Ctest, nodesize = 6, nsplit = 1, seed = 3)
Cerror <- gg_error(Cout.rsf.3)
lowest <- min( Cerror$error[Cerror$error!=min(Cerror$error)] )
Crow <- Cerror[which(Cerror == lowest),]
plot(gg_error(Cout.rsf.3))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

Cnumtree = 900 #Crow$ntree
num_row <- nrow(Ctrain)
#RANDOM SURVIVAL FOREST
Csurvival_modelrfsrc<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Ctrain, mtry = sqrt(196), 
                            ntree = Cnumtree, cv.fold = 10, verbose = TRUE, block.size=1,
                            importance=TRUE, nodesize = 6, nsplit = 1,seed = 4)
set.seed(3)
Csurv.fit11 <- survivalROC(Ctrain$MONTHS, Ctrain$STATUS, Csurvival_modelrfsrc$predicted.oob, 
                          predict.time = 11, method = "KM")
Csurv.fit16 <- survivalROC(Ctrain$MONTHS, Ctrain$STATUS, Csurvival_modelrfsrc$predicted.oob, 
                          predict.time = 16, method = "KM")
Csurv.fit20 <- survivalROC(Ctrain$MONTHS, Ctrain$STATUS, Csurvival_modelrfsrc$predicted.oob, 
                          predict.time = 20, method = "KM")

plot(Csurv.fit11$FP, Csurv.fit11$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Cnumtree," Trees for Phenotypic, Genmoic and Clinical data")) 
lines(Csurv.fit16$FP, Csurv.fit16$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
lines(Csurv.fit20$FP, Csurv.fit20$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
abline(0,1)
Cvariable11 <- paste0("11 month AUC =", round(Csurv.fit11$AUC,3))
Cvariable16 <- paste0("16 month AUC =", round(Csurv.fit16$AUC,3))
Cvariable20 <- paste0("20 month AUC =", round(Csurv.fit20$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Cvariable11,Cvariable16,Cvariable20))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Cfit <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = CDatawD)
print(Cfit)
summary(Cfit)$table
ggsurvplot(Cfit,
           data = CDatawD,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Cpred <- predict(Csurvival_modelrfsrc,Ctest[, (colnames(Ctest))],block.size=1)
plot(Csurvival_modelrfsrc$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))
legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Ctrain$MONTHS, censoring = Ctrain$STATUS, predicted = Csurvival_modelrfsrc$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
Cbs.km <- get.brier.survival(Csurvival_modelrfsrc, cens.mode = "km")$brier.score
Cbs.rsf <- get.brier.survival(Csurvival_modelrfsrc, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(Cbs.km, type = "s", col = 2)
lines(Cbs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------
#Section 3: Dataset with Genomic Data Alone
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 10-Genomic Data
#----------------------------------------------------------------------------------------------------------
Gdata = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL - RSF/Genomic-Data.csv") #Gdata = Genomic Data
Gdata <- na.omit(Gdata)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
GDataR <- Gdata %>% select(-c("RISK")) #Genomic Data no Risk
Gdf <- Gdata %>% select(-c("PATIENT","RISK")) #Genomic Dataframe without Risk and Patient
GfeatureNames <- paste(colnames(Gdf), collapse = " + ") #Genomic Feature Names
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Gsurvival_formula <-paste('Surv(MONTHS,STATUS) ~', GfeatureNames) #Clinical Survival Formula
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
#set.seed(12)
#set.seed(13)
#set.seed(144)#13 and 14 have similar splits
set.seed(15)
#set.seed(16)
sample <- sample.split(GDataR$PATIENT, SplitRatio = 0.3)
Gtrain  <- subset (GDataR%>% select(-c("PATIENT")), sample == FALSE)
Gtest   <- subset(GDataR %>% select(-c("PATIENT")), sample == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
Gout.rsf.3 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                    data = Gtrain, 
                    mtry = sqrt(196), 
                    ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                    importance=TRUE, xtest = Gtest, nodesize = 6, nsplit = 1, seed = 3)
Gerror <- gg_error(Gout.rsf.3)
Glowest <- min( Gerror$error[Gerror$error!=min(Gerror$error)] )
Grow <- Gerror[which(Gerror == lowest),]
plot(gg_error(Gout.rsf.3))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

Gnumtree = 900 #Grow$ntree[1]
#RANDOM SURVIVAL FOREST
Gsurvival_modelrfsrc<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Gtrain, mtry = sqrt(196), 
                            ntree = numtree, cv.fold = 10, verbose = TRUE, block.size=1,
                            importance=TRUE, nodesize = 6, nsplit = 1, seed = 4)
set.seed(3)
Gsurv.fit11 <- survivalROC(Gtrain$MONTHS, Gtrain$STATUS, Gsurvival_modelrfsrc$predicted.oob, 
                          predict.time = 11, method = "KM")
Gsurv.fit16 <- survivalROC(Gtrain$MONTHS, Gtrain$STATUS, Gsurvival_modelrfsrc$predicted.oob, 
                          predict.time = 16, method = "KM")
Gsurv.fit20 <- survivalROC(Gtrain$MONTHS, Gtrain$STATUS, Gsurvival_modelrfsrc$predicted.oob, 
                          predict.time = 20, method = "KM")

plot(Gsurv.fit11$FP, Gsurv.fit11$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Gnumtree," Trees for Genomic Data"))
lines(Gsurv.fit16$FP, Gsurv.fit16$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
lines(Gsurv.fit20$FP, Gsurv.fit20$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
abline(0,1)
Gvariable11 <- paste0("11 month AUC =", round(Gsurv.fit11$AUC,3))
Gvariable16 <- paste0("16 month AUC =", round(Gsurv.fit16$AUC,3))
Gvariable20 <- paste0("20 month AUC =", round(Gsurv.fit20$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Gvariable11,Gvariable16,Gvariable20))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Gfit <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = Gdata)
print(Gfit)
summary(Gfit)$table
ggsurvplot(Gfit,
           data = Gdata,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Gpred <- predict(Gsurvival_modelrfsrc,Gtest[, (colnames(Gtest))],block.size=1)
plot(Gsurvival_modelrfsrc$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))
legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Gtrain$MONTHS, censoring = Gtrain$STATUS, predicted = Gsurvival_modelrfsrc$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
Gbs.km <- get.brier.survival(Gsurvival_modelrfsrc, cens.mode = "km")$brier.score
Gbs.rsf <- get.brier.survival(Gsurvival_modelrfsrc, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(Gbs.km, type = "s", col = 2)
lines(Gbs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------
#Section 4: Phenotypic Data Alone
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 11-Phenotypic Data
#----------------------------------------------------------------------------------------------------------
Pdata = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL - RSF/Phenotypic-Data.csv") #Pdata = Phenotypic Data
Pdata <- na.omit(Pdata)
#----------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
PDataR <- Pdata %>% select(-c("RISK")) #Genomic Data no Risk
Pdf <- Pdata %>% select(-c("PATIENT","RISK")) #Genomic Dataframe without Risk and Patient
PfeatureNames <- paste(colnames(Pdf), collapse = " + ") #Genomic Feature Names
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Psurvival_formula <-paste('Surv(MONTHS,STATUS) ~', PfeatureNames) #Clinical Survival Formula
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
#set.seed(12)
set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)
sample <- sample.split(PDataR$PATIENT, SplitRatio = 0.3)
Ptrain  <- subset (PDataR%>% select(-c("PATIENT")), sample == FALSE)
Ptest   <- subset(PDataR %>% select(-c("PATIENT")), sample == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
Pout.rsf.3 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                     data = Ptrain, 
                     mtry = sqrt(196), 
                     ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                     importance=TRUE, xtest = Ptest, nodesize = 6, nsplit = 1, seed = 3)
Perror <- gg_error(Pout.rsf.3)
Plowest <- min( Perror$error[Perror$error!=min(Perror$error)] )
Prow <- Perror[which(Perror == lowest),]
plot(gg_error(Pout.rsf.3))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)

Pnumtree = 900 #Prow$ntree
#RANDOM SURVIVAL FOREST
Psurvival_modelrfsrc<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Ptrain, mtry = sqrt(196), 
                            ntree = numtree, cv.fold = 10, verbose = TRUE, block.size=1,
                            importance=TRUE, nodesize = 6, nsplit = 1, seed = 4)
set.seed(3)
Psurv.fit11 <- survivalROC(Ptrain$MONTHS, Ptrain$STATUS, Psurvival_modelrfsrc$predicted.oob, 
                           predict.time = 11, method = "KM")
Psurv.fit16 <- survivalROC(Ptrain$MONTHS, Ptrain$STATUS, Psurvival_modelrfsrc$predicted.oob, 
                           predict.time = 16, method = "KM")
Psurv.fit20 <- survivalROC(Ptrain$MONTHS, Ptrain$STATUS, Psurvival_modelrfsrc$predicted.oob, 
                           predict.time = 20, method = "KM")

plot(Psurv.fit11$FP, Psurv.fit11$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Pnumtree," Trees for Phenotypic Data"))
lines(Psurv.fit16$FP, Psurv.fit16$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
lines(Psurv.fit20$FP, Psurv.fit20$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
abline(0,1)
Pvariable11 <- paste0("11 month AUC =", round(Psurv.fit11$AUC,3))
Pvariable16 <- paste0("16 month AUC =", round(Psurv.fit16$AUC,3))
Pvariable20 <- paste0("20 month AUC =", round(Psurv.fit20$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Pvariable11,Pvariable16,Pvariable20))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Pfit <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = Pdata)
print(Pfit)
summary(Pfit)$table
ggsurvplot(Pfit,
           data = Pdata,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Ppred <- predict(Psurvival_modelrfsrc,Ptest[, (colnames(Ptest))],block.size=1)
plot(Psurvival_modelrfsrc$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))
legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Ptrain$MONTHS, censoring = Ptrain$STATUS, predicted = Psurvival_modelrfsrc$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
Pbs.km <- get.brier.survival(Psurvival_modelrfsrc, cens.mode = "km")$brier.score
Pbs.rsf <- get.brier.survival(Psurvival_modelrfsrc, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(Pbs.km, type = "s", col = 2)
lines(Pbs.rsf, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------
#Section 5: Egiengenes and Phenotypic data
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 16-THESIS_GENOMIC&PHENOTYPIC_DATA
#----------------------------------------------------------------------------------------------------------
data2 = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL RUN THESIS/16-THESIS_GENOMIC&PHENOTYPIC_DATA.csv") 
data2 <- na.omit(data2)
#----------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
noriskdata2 = data2 %>% select(-c("RISK"))
df2 <- data2 %>% select(-c("PATIENT","RISK"))
featureNames2 <- paste(colnames(df2), collapse = " + ")
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
survival_formula2 <-paste('Surv(MONTHS,STATUS) ~', featureNames2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)

sample2 <- sample.split(noriskdata2$PATIENT, SplitRatio = 0.3)
train2  <- subset (noriskdata2%>% select(-c("PATIENT")), sample2 == FALSE)
test2   <- subset(noriskdata2%>% select(-c("PATIENT")), sample2 == TRUE) #this was false

# Split for RSF Kaplan Meier
# Note: Use the same seed for previous split
sample2_withRISK <- sample.split(data2$PATIENT, SplitRatio = 0.3)
train2_withRISK  <- subset (data2%>% select(-c("PATIENT")), sample2_withRISK == FALSE)
test2_withRISK   <- subset(data2%>% select(-c("PATIENT")), sample2_withRISK == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#GRID SEARCH FOR OPTIMUM nodesize and mtry
#----------------------------------------------------------------------------------------------------------
#set default values
nodesize <- c(5,10,20)

#model to solve for optimum values using 1000 trees
nodemodel2 <- tune(Surv(MONTHS, STATUS) ~ ., data = train2,
                    mtryStart = 2,
                    nodesizeTry= nodesize,
                    cv.fold = 10,
                    nsplitTry = 10, ntree=1000,
                    blocksize=1, importance=TRUE, seed=15)

grid2 <- as.data.frame(nodemodel2$results)
low2 <- (grid2[grid2$err == min(grid2$err), ])

#Optimized Values
nodesize2 <- low2$nodesize
mtry2 <- low2$mtry
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
out.rsf.3.2 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                    data = train2, 
                    mtry = mtry2, 
                    ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                    importance=TRUE, nodesize = nodesize2, nsplit = 10, seed = 3)#6,7 has higher test performance)#3) current paper

error2 <- gg_error(out.rsf.3.2)
lowest2 <- min( error2$error[error2$error!=min(error2$error)] )
row2 <- error2[which(error2 == lowest2),]
plot(gg_error(out.rsf.3.2))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

numtree = row2$ntree[1]
numtree
numtree <- 23
mtry2 <- 11
nodesize2 <- 10

#RANDOM SURVIVAL FOREST
survival_modelrfsrc2<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = train2, mtry = mtry2, 
                           ntree = numtree, cv.fold = 10, verbose = TRUE, block.size=1,
                           importance=TRUE, nodesize = nodesize2, nsplit = 10, seed = 4)

set.seed(3)
surv.fit11.2 <- survivalROC(train2$MONTHS, train2$STATUS, survival_modelrfsrc2$predicted.oob, 
                          predict.time = 11, method = "KM")
surv.fit16.2 <- survivalROC(train2$MONTHS, train2$STATUS, survival_modelrfsrc2$predicted.oob, 
                          predict.time = 16, method = "KM")
surv.fit20.2 <- survivalROC(train2$MONTHS, train2$STATUS, survival_modelrfsrc2$predicted.oob, 
                          predict.time = 20, method = "KM")

plot(surv.fit11.2$FP, surv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",numtree," Trees Genomic and Phenotypic Data \n Parameters: mtry (", mtry2,
                             ") nodesize (", nodesize2,")"))
lines(surv.fit16.2$FP, surv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
lines(surv.fit20.2$FP, surv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
variable11.2 <- paste0("11 month AUC =", round(surv.fit11.2$AUC,3))
variable16.2 <- paste0("16 month AUC =", round(surv.fit16.2$AUC,3))
variable20.2 <- paste0("20 month AUC =", round(surv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(variable11.2,variable16.2,variable20.2))
#----------------------------------------------------------------------------------------------------------
#RANDOM SURVIVAL FOREST TEST

pred2 <- predict(survival_modelrfsrc2,test2)

set.seed(3)
surv.fit11.2 <- survivalROC(test2$MONTHS, test2$STATUS, pred2$predicted, 
                            predict.time = 11, method = "KM")
surv.fit16.2 <- survivalROC(test2$MONTHS, test2$STATUS, pred2$predicted, 
                            predict.time = 16, method = "KM")
surv.fit20.2 <- survivalROC(test2$MONTHS, test2$STATUS, pred2$predicted, 
                            predict.time = 20, method = "KM")

plot(surv.fit11.2$FP, surv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",numtree," Trees Genomic and Phenotypic Data \n Parameters: mtry (", mtry2,
                             ") nodesize (", nodesize2,")"))
lines(surv.fit16.2$FP, surv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
lines(surv.fit20.2$FP, surv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
variable11.2 <- paste0("11 month AUC =", round(surv.fit11.2$AUC ,3))
variable16.2 <- paste0("16 month AUC =", round(surv.fit16.2$AUC,3))
variable20.2 <- paste0("20 month AUC =", round(surv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(variable11.2,variable16.2,variable20.2))

#----------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
fit2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = data2)
print(fit2)
summary(fit2)$table
ggsurvplot(fit2,
           data = data2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
pred2 <- predict(survival_modelrfsrc2,test2[, (colnames(test2))],block.size=1)
plot(pred2$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))

legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = train2$MONTHS, censoring = train2$STATUS, predicted = survival_modelrfsrc2$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
bs.km.2 <- get.brier.survival(survival_modelrfsrc2, cens.mode = "km")$brier.score
bs.rsf.2 <- get.brier.survival(survival_modelrfsrc2, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(bs.km.2, type = "s", col = 2)
lines(bs.rsf.2, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#KAPLAN MEIER of RSF Results
#----------------------------------------------------------------------------------------------------------
Run2 <- test2
Run2_withRISK <- test2_withRISK 

surv_prob2 <- predict(survival_modelrfsrc2, newdata = Run2, type = "survival" ) #type is the same thing as "response" 
surv_times2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = Run2_withRISK, weights = surv_prob2$survival[,11] )
#plot(surv_times, main = "Kaplan-Meier Plot", xlab = "Time", ylab = "Survival Probability")

ggsurvplot(surv_times2,
           data = Run2_withRISK,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#----------------------------------------------------------------------------------------------------------





#----------------------------------------------------------------------------------------------------------
#Section 6: Eigengenes, Phenotypic and Clinical
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 12.1-THESIS_OVERALL_DATA_RISK
#----------------------------------------------------------------------------------------------------------
Cdata2 = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL RUN THESIS/12.1-THESIS_OVERALL_DATA_RISK.csv") 
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
Cdata2$Male <- ifelse(Cdata2$SEX == 'Male', 1, 0)
Cdata2$Female <- ifelse(Cdata2$SEX == 'Female', 1, 0)
Cdata2$White <- ifelse(Cdata2$RACE == 'WHITE', 1, 0)
Cdata2$Black <- ifelse(Cdata2$RACE == 'BLACK OR AFRICAN AMERICAN', 1, 0)
Cdata2$Asian <- ifelse(Cdata2$RACE == 'ASIAN', 1, 0)
Cdata2$Indian <- ifelse(Cdata2$RACE == 'AMERICAN INDIAN OR ALASKA NATIVE', 1, 0)
Cdata2$T1 <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1', 1, 0)
Cdata2$T1a <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1a', 1, 0)
Cdata2$T1b <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1b', 1, 0)
Cdata2$T2 <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2', 1, 0)
Cdata2$T2a <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2a', 1, 0)
Cdata2$T2b <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2b', 1, 0)
Cdata2$T3 <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3', 1, 0)
Cdata2$T3a <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3a', 1, 0)
Cdata2$T3b <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3b', 1, 0)
Cdata2$T3c <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3c', 1, 0)
Cdata2$T4 <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T4', 1, 0)
Cdata2$TX <- ifelse(Cdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'TX', 1, 0)
CDatawD2 <- Cdata2 %>% select(-c("SEX","RACE","AJCC_TUMOR_PATHOLOGIC_PT")) #Clinical data with dummy variables
CDataR2 <- CDatawD2 %>% select(-c("RISK")) #Clinical Data no Risk
Cdf2 <- CDatawD2 %>% select(-c("PATIENT","RISK")) #Clinical Dataframe
CfeatureNames2 <- paste(colnames(Cdf2), collapse = " + ") #Clinical Feature Names
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Csurvival_formula2 <-paste('Surv(MONTHS,STATUS) ~', CfeatureNames2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)
CDataR2 <- na.omit(CDataR2)
Csample2 <- sample.split(CDataR2$PATIENT, SplitRatio = 0.3)
Ctrain2  <- subset (CDataR2%>% select(-c("PATIENT")), Csample2 == FALSE)
Ctest2   <- subset(CDataR2%>% select(-c("PATIENT")), Csample2 == TRUE) #this was false

# Split for RSF Kaplan Meier
# Note: Use the same seed for previous split
Csample2_withRISK <- sample.split(CDatawD2$PATIENT, SplitRatio = 0.3)
Ctrain2_withRISK  <- subset (CDatawD2%>% select(-c("PATIENT")), Csample2_withRISK == FALSE)
Ctest2_withRISK   <- subset(CDatawD2%>% select(-c("PATIENT")), Csample2_withRISK == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#GRID SEARCH FOR OPTIMUM nodesize and mtry
#----------------------------------------------------------------------------------------------------------
#set default values  start again 
nodesize <- c(5,10,20)
nsplit <- c(5,10)

#model to solve for optimum values using 1000 trees
Cnodemodel2 <- tune(Surv(MONTHS, STATUS) ~ ., data = Ctrain2,
                   mtryStart = 2,
                   nodesizeTry= nodesize,
                   cv.fold = 10,
                   nsplitTry = 10, ntree=1000,
                   blocksize=1, importance=TRUE, seed=8)#8 

Cgrid2 <- as.data.frame(Cnodemodel2$results)
Clow2 <- (Cgrid2[Cgrid2$err == min(Cgrid2$err), ])

#Optimized Values
Cnodesize2 <- Clow2$nodesize
Cnodesize2_ <- as.data.frame(Cnodesize2)
Cnodesize2 <- Cnodesize2_[2,1]

Cmtry2 <- Clow2$mtry
Cmtry2_ <- as.data.frame(Cmtry2)
Cmtry2 <- Cmtry2_[2,1]
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
C.out.rsf.3.2 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                      data = Ctrain2, 
                      mtry = Cmtry2, #sqrt(196), 
                      ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                      importance=TRUE, nodesize = Cnodesize2, #6, 
                      nsplit = 10, seed = 3)#3  13 16

Cerror2 <- gg_error(C.out.rsf.3.2)
Clowest2 <- min(Cerror2$error[Cerror2$error!=min(Cerror2$error)] )
Crow2 <- Cerror2[which(Cerror2 == Clowest2),]
plot(gg_error(C.out.rsf.3.2))
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)


Cnumtree2 = Crow2$ntree[1]

#Use this if you want to manually test specific parameters
Cnumtree2 = 23
Cmtry2 = 10
Cnodesize2 = 11

#RANDOM SURVIVAL FOREST
Csurvival_modelrfsrc2<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Ctrain2, mtry = Cmtry2,#sqrt(196), 
                            ntree = Cnumtree2, 
                            cv.fold = 10, verbose = TRUE, block.size=1,
                            importance=TRUE, nodesize = Cnodesize2, 
                            nsplit = 10, seed = 4)#4  11, 16 GOOD 

set.seed(3)
Csurv.fit11.2 <- survivalROC(Ctrain2$MONTHS, Ctrain2$STATUS, Csurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 11, method = "KM")
Csurv.fit16.2 <- survivalROC(Ctrain2$MONTHS, Ctrain2$STATUS, Csurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 16, method = "KM")
Csurv.fit20.2 <- survivalROC(Ctrain2$MONTHS, Ctrain2$STATUS, Csurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 20, method = "KM")

plot(Csurv.fit11.2$FP, Csurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Cnumtree2," Trees Clinical, Genomic, and Phenotypic Data \n Parameters: mtry (", Cmtry2,
                             ") nodesize (", Cnodesize2,")"))
lines(Csurv.fit16.2$FP, Csurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP"),
      ylab="TP")
lines(Csurv.fit20.2$FP, Csurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP"),
      ylab="TP")
abline(0,1)
Cvariable11.2 <- paste0("11 month AUC =", round(Csurv.fit11.2$AUC,3))
Cvariable16.2 <- paste0("16 month AUC =", round(Csurv.fit16.2$AUC,3))
Cvariable20.2 <- paste0("20 month AUC =", round(Csurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Cvariable11.2,Cvariable16.2,Cvariable20.2))
#----------------------------------------------------------------------------------------------------------
#RANDOM SURVIVAL FOREST TEST

Cpred2 <- predict(Csurvival_modelrfsrc2,Ctest2)

set.seed(3)
Csurv.fit11.2 <- survivalROC(Ctest2$MONTHS, Ctest2$STATUS, Cpred2$predicted, 
                            predict.time = 11, method = "KM")
Csurv.fit16.2 <- survivalROC(Ctest2$MONTHS, Ctest2$STATUS, Cpred2$predicted, 
                            predict.time = 16, method = "KM")
Csurv.fit20.2 <- survivalROC(Ctest2$MONTHS, Ctest2$STATUS, Cpred2$predicted, 
                            predict.time = 20, method = "KM")

plot(Csurv.fit11.2$FP, Csurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Cnumtree2," Trees Clinical, Genomic, and Phenotypic Data \n Parameters: mtry (", Cmtry2,
                             ") nodesize (", Cnodesize2,")"))
lines(Csurv.fit16.2$FP, Csurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP"),
      ylab="TP")
lines(Csurv.fit20.2$FP, Csurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP"),
      ylab="TP")
abline(0,1)
Cvariable11.2 <- paste0("11 month AUC =", round(Csurv.fit11.2$AUC ,3))
Cvariable16.2 <- paste0("16 month AUC =", round(Csurv.fit16.2$AUC,3))
Cvariable20.2 <- paste0("20 month AUC =", round(Csurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Cvariable11.2,Cvariable16.2,Cvariable20.2))

#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Cfit2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = CDatawD2)
print(Cfit2)
summary(Cfit2)$table
ggsurvplot(Cfit2,
           data = CDatawD2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Cpred2 <- predict(Csurvival_modelrfsrc2,Ctest2[, (colnames(Ctest2))],block.size=1)
plot(Cpred2$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))

legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Ctrain2$MONTHS, censoring = Ctrain2$STATUS, predicted = Csurvival_modelrfsrc2$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
C.bs.km.2 <- get.brier.survival(Csurvival_modelrfsrc2, cens.mode = "km")$brier.score
C.bs.rsf.2 <- get.brier.survival(Csurvival_modelrfsrc2, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(C.bs.km.2, type = "s", col = 2)
lines(C.bs.rsf.2, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#KAPLAN MEIER of RSF Results
#----------------------------------------------------------------------------------------------------------
CRun2 <- Ctest2
CRun2_withRISK <- Ctest2_withRISK 

Csurv_prob2 <- predict(Csurvival_modelrfsrc2, newdata = CRun2, type = "survival" ) #type is the same thing as "response" 
Csurv_times2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = CRun2_withRISK, weights = Csurv_prob2$survival[,11] )
#plot(surv_times, main = "Kaplan-Meier Plot", xlab = "Time", ylab = "Survival Probability")

ggsurvplot(Csurv_times2,
           data = CRun2_withRISK,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#----------------------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------------------
#Section 7: Egiengenes data
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 16-THESIS_GENOMIC_DATA
#----------------------------------------------------------------------------------------------------------
Gdata2 = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL RUN THESIS/14-THESIS_GENOMIC_DATA.csv") 
Gdata2 <- na.omit(Gdata2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
Gnoriskdata2 = Gdata2 %>% select(-c("RISK"))
Gdf2 <- Gdata2 %>% select(-c("PATIENT","RISK"))
GfeatureNames2 <- paste(colnames(Gdf2), collapse = " + ")
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Gsurvival_formula2 <-paste('Surv(MONTHS,STATUS) ~', GfeatureNames2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)

Gsample2 <- sample.split(Gnoriskdata2$PATIENT, SplitRatio = 0.3)
Gtrain2  <- subset (Gnoriskdata2%>% select(-c("PATIENT")), Gsample2 == FALSE)
Gtest2   <- subset(Gnoriskdata2%>% select(-c("PATIENT")), Gsample2 == TRUE) #this was false

# Split for RSF Kaplan Meier
# Note: Use the same seed for previous split
Gsample2_withRISK <- sample.split(Gdata2$PATIENT, SplitRatio = 0.3)
Gtrain2_withRISK  <- subset (Gdata2%>% select(-c("PATIENT")), Gsample2_withRISK == FALSE)
Gtest2_withRISK   <- subset(Gdata2%>% select(-c("PATIENT")), Gsample2_withRISK == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#GRID SEARCH FOR OPTIMUM nodesize and mtry
#----------------------------------------------------------------------------------------------------------
#set default values
nodesize <- c(5,10,20)

#model to solve for optimum values using 1000 trees
Gnodemodel2 <- tune(Surv(MONTHS, STATUS) ~ ., data = Gtrain2,
                    mtryStart = 2,
                    nodesizeTry= nodesize,
                    cv.fold = 10,
                    nsplitTry = 10, ntree=1000,
                    blocksize=1, importance=TRUE, seed=3)

Ggrid2 <- as.data.frame(Gnodemodel2$results)
Glow2 <- (Ggrid2[Ggrid2$err == min(Ggrid2$err), ])

#Optimized Values
Gnodesize2 <- Glow2$nodesize
Gmtry2 <- Glow2$mtry
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
Gout.rsf.3.2 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                      data = Gtrain2, 
                      mtry = Gmtry2, 
                      ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                      importance=TRUE, nodesize = Gnodesize2, nsplit = 13, seed = 3)

Gerror2 <- gg_error(Gout.rsf.3.2)
Glowest2 <- min(Gerror2$error[Gerror2$error!=min(Gerror2$error)] )
Grow2 <- Gerror2[which(Gerror2 == Glowest2),]
plot(gg_error(Gout.rsf.3.2))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

Gnumtree2 = Grow2$ntree[1]

Gnumtree2 <-186
Gmtry2 <-6
Gnodesize2 <- 5

#RANDOM SURVIVAL FOREST
Gsurvival_modelrfsrc2<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Gtrain2, mtry = Gmtry2, 
                            ntree = Gnumtree2, cv.fold = 10, verbose = TRUE, block.size=1,
                            importance=TRUE, nodesize = Gnodesize2, nsplit = 10, seed = 4)
set.seed(3)
Gsurv.fit11.2 <- survivalROC(Gtrain2$MONTHS, Gtrain2$STATUS, Gsurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 11, method = "KM")
Gsurv.fit16.2 <- survivalROC(Gtrain2$MONTHS, Gtrain2$STATUS, Gsurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 16, method = "KM")
Gsurv.fit20.2 <- survivalROC(Gtrain2$MONTHS, Gtrain2$STATUS, Gsurvival_modelrfsrc2$predicted.oob, 
                            predict.time = 20, method = "KM")

plot(Gsurv.fit11.2$FP, Gsurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n" ),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Gnumtree2," Trees Genomic Data \n Parameters: mtry (", Gmtry2,
                             ") nodesize (", Gnodesize2,")"))
lines(Gsurv.fit16.2$FP, Gsurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n" ),
      ylab="TP")
lines(Gsurv.fit20.2$FP, Gsurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
Gvariable11.2 <- paste0("11 month AUC =", round(Gsurv.fit11.2$AUC,3))
Gvariable16.2 <- paste0("16 month AUC =", round(Gsurv.fit16.2$AUC,3))
Gvariable20.2 <- paste0("20 month AUC =", round(Gsurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Gvariable11.2,Gvariable16.2,Gvariable20.2))
#----------------------------------------------------------------------------------------------------------
#RANDOM FOREST TEST
Gpred2 <- predict(Gsurvival_modelrfsrc2,Gtest2)

set.seed(3)
Gsurv.fit11.2 <- survivalROC(Gtest2$MONTHS, Gtest2$STATUS, Gpred2$predicted, 
                             predict.time = 11, method = "KM")
Gsurv.fit16.2 <- survivalROC(Gtest2$MONTHS, Gtest2$STATUS, Gpred2$predicted, 
                             predict.time = 16, method = "KM")
Gsurv.fit20.2 <- survivalROC(Gtest2$MONTHS, Gtest2$STATUS, Gpred2$predicted, 
                             predict.time = 20, method = "KM")

plot(Gsurv.fit11.2$FP, Gsurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n" ),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Gnumtree2," Trees Genomic Data \n Parameters: mtry (", Gmtry2,
                             ") nodesize (", Gnodesize2,")"))
lines(Gsurv.fit16.2$FP, Gsurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n" ),
      ylab="TP")
lines(Gsurv.fit20.2$FP, Gsurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n" ),
      ylab="TP")
abline(0,1)
Gvariable11.2 <- paste0("11 month AUC =", round(Gsurv.fit11.2$AUC,3))
Gvariable16.2 <- paste0("16 month AUC =", round(Gsurv.fit16.2$AUC,3))
Gvariable20.2 <- paste0("20 month AUC =", round(Gsurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Gvariable11.2,Gvariable16.2,Gvariable20.2))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Gfit2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = Gdata2)
print(Gfit2)
summary(Gfit2)$table
ggsurvplot(Gfit2,
           data = Gdata2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Gpred2 <- predict(Gsurvival_modelrfsrc2,Gtest2[, (colnames(Gtest2))],block.size=1)
plot(Gpred2$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))

legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Gtrain2$MONTHS, censoring = Gtrain2$STATUS, predicted = Gsurvival_modelrfsrc2$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
Gbs.km.2 <- get.brier.survival(Gsurvival_modelrfsrc2, cens.mode = "km")$brier.score
Gbs.rsf.2 <- get.brier.survival(Gsurvival_modelrfsrc2, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(Gbs.km.2, type = "s", col = 2)
lines(Gbs.rsf.2, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#KAPLAN MEIER of RSF Results
#----------------------------------------------------------------------------------------------------------
GRun2 <- Gtest2
GRun2_withRISK <- Gtest2_withRISK 

Gsurv_prob2 <- predict(Gsurvival_modelrfsrc2, newdata = GRun2, type = "survival" ) #type is the same thing as "response" 
Gsurv_times2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = GRun2_withRISK, weights = Gsurv_prob2$survival[,11] )
#plot(surv_times, main = "Kaplan-Meier Plot", xlab = "Time", ylab = "Survival Probability")

ggsurvplot(Gsurv_times2,
           data = GRun2_withRISK,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#----------------------------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------------------------
#Section 8: Phenotypic data
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 15-THESIS_PHENOTYPIC_DATA
#----------------------------------------------------------------------------------------------------------
Pdata2 = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL RUN THESIS/15-THESIS_PHENOTYPIC_DATA.csv") 
Pdata2 <- na.omit(Pdata2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
Pnoriskdata2 = Pdata2 %>% select(-c("RISK"))
Pdf2 <- Pdata2 %>% select(-c("PATIENT","RISK"))
PfeatureNames2 <- paste(colnames(Pdf2), collapse = " + ")
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
Psurvival_formula2 <-paste('Surv(MONTHS,STATUS) ~', PfeatureNames2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)

Psample2 <- sample.split(Pnoriskdata2$PATIENT, SplitRatio = 0.3)
Ptrain2  <- subset (Pnoriskdata2%>% select(-c("PATIENT")), Psample2 == FALSE)
Ptest2   <- subset(Pnoriskdata2%>% select(-c("PATIENT")), Psample2 == TRUE) #this was false

# Split for RSF Kaplan Meier
# Note: Use the same seed for previous split
Psample2_withRISK <- sample.split(Pdata2$PATIENT, SplitRatio = 0.3)
Ptrain2_withRISK  <- subset (Pdata2%>% select(-c("PATIENT")), Psample2_withRISK == FALSE)
Ptest2_withRISK   <- subset(Pdata2%>% select(-c("PATIENT")), Psample2_withRISK == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#GRID SEARCH FOR OPTIMUM nodesize and mtry
#----------------------------------------------------------------------------------------------------------
#set default values
nodesize <- c(5,10,20)

#model to solve for optimum values using 1000 trees
Pnodemodel2 <- tune(Surv(MONTHS, STATUS) ~ ., data = Ptrain2,
                    mtryStart = 2,
                    nodesizeTry= nodesize,
                    cv.fold = 10,
                    nsplitTry = 10, ntree=1000,
                    blocksize=1, importance=TRUE, seed=5)

Pgrid2 <- as.data.frame(Pnodemodel2$results)
Plow2 <- (Pgrid2[Pgrid2$err == min(Pgrid2$err), ])

#Optimized Values
Pnodesize2 <- Plow2$nodesize
Pmtry2 <- Plow2$mtry
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
Pout.rsf.3.2 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                       data = Ptrain2, 
                       mtry = Pmtry2, 
                       ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                       importance=TRUE, nodesize = Pnodesize2, nsplit = 10, seed = 3)

Perror2 <- gg_error(Pout.rsf.3.2)
Plowest2 <- min(Perror2$error[Perror2$error!=min(Perror2$error)] )
Prow2 <- Perror2[which(Perror2 == Plowest2),]
plot(gg_error(Pout.rsf.3.2))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)

Pnumtree2 = Prow2$ntree[1]

Pnumtree2 <- 73
Pmtry2 <- 1
Pnodesize2 <- 10
#RANDOM SURVIVAL FOREST
Psurvival_modelrfsrc2<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = Ptrain2, mtry = Pmtry2, 
                             ntree = Pnumtree2, cv.fold = 10, verbose = TRUE, block.size=1,
                             importance=TRUE, nodesize = Pnodesize2, nsplit = 10, seed = 4)
set.seed(3)
Psurv.fit11.2 <- survivalROC(Ptrain2$MONTHS, Ptrain2$STATUS, Psurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 11, method = "KM")
Psurv.fit16.2 <- survivalROC(Ptrain2$MONTHS, Ptrain2$STATUS, Psurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 16, method = "KM")
Psurv.fit20.2 <- survivalROC(Ptrain2$MONTHS, Ptrain2$STATUS, Psurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 20, method = "KM")

plot(Psurv.fit11.2$FP, Psurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Pnumtree2," Trees Phenotypic Data \n Parameters: mtry (", Pmtry2,
                             ") nodesize (", Pnodesize2,")"))
lines(Psurv.fit16.2$FP, Psurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
lines(Psurv.fit20.2$FP, Psurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
Pvariable11.2 <- paste0("11 month AUC =", round(Psurv.fit11.2$AUC,3))
Pvariable16.2 <- paste0("16 month AUC =", round(Psurv.fit16.2$AUC,3))
Pvariable20.2 <- paste0("20 month AUC =", round(Psurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Pvariable11.2,Pvariable16.2,Pvariable20.2))
#----------------------------------------------------------------------------------------------------------
#RANDOM FOREST TEST
Ppred2 <- predict(Psurvival_modelrfsrc2,Ptest2)

set.seed(3)
Psurv.fit11.2 <- survivalROC(Ptest2$MONTHS, Ptest2$STATUS, Ppred2$predicted, 
                             predict.time = 11, method = "KM")
Psurv.fit16.2 <- survivalROC(Ptest2$MONTHS, Ptest2$STATUS, Ppred2$predicted, 
                             predict.time = 16, method = "KM")
Psurv.fit20.2 <- survivalROC(Ptest2$MONTHS, Ptest2$STATUS, Ppred2$predicted, 
                             predict.time = 20, method = "KM")

plot(Psurv.fit11.2$FP, Psurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",Pnumtree2," Trees Phenotypic Data \n Parameters: mtry (", Pmtry2,
                             ") nodesize (", Pnodesize2,")"))
lines(Psurv.fit16.2$FP, Psurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
lines(Psurv.fit20.2$FP, Psurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
Pvariable11.2 <- paste0("11 month AUC =", round(Psurv.fit11.2$AUC,3))
Pvariable16.2 <- paste0("16 month AUC =", round(Psurv.fit16.2$AUC,3))
Pvariable20.2 <- paste0("20 month AUC =", round(Psurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(Pvariable11.2,Pvariable16.2,Pvariable20.2))
#----------------------------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
Pfit2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = Pdata2)
print(Pfit2)
summary(Pfit2)$table
ggsurvplot(Pfit2,
           data = Pdata2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
Ppred2 <- predict(Psurvival_modelrfsrc2,Ptest2[, (colnames(Ptest2))],block.size=1)
plot(Ppred2$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,1))

legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = Ptrain2$MONTHS, censoring = Ptrain2$STATUS, predicted = Psurvival_modelrfsrc2$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
Pbs.km.2 <- get.brier.survival(Psurvival_modelrfsrc2, cens.mode = "km")$brier.score
Pbs.rsf.2 <- get.brier.survival(Psurvival_modelrfsrc2, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(Pbs.km.2, type = "s", col = 2)
lines(Pbs.rsf.2, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#KAPLAN MEIER of RSF Results
#----------------------------------------------------------------------------------------------------------
PRun2 <- Ptest2
PRun2_withRISK <- Ptest2_withRISK 

Psurv_prob2 <- predict(Psurvival_modelrfsrc2, newdata = PRun2, type = "survival" ) #type is the same thing as "response" 
Psurv_times2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = PRun2_withRISK, weights = Psurv_prob2$survival[,11] )
#plot(surv_times, main = "Kaplan-Meier Plot", xlab = "Time", ylab = "Survival Probability")

ggsurvplot(Psurv_times2,
           data = PRun2_withRISK,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Section 9: Clinical
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Load in the Dataset 12.1-THESIS_OVERALL_DATA_RISK
#----------------------------------------------------------------------------------------------------------
CCdata2 = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FINAL RUN THESIS/17-THESIS_CLINICAL_DATA_RISK.csv") 
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Dataset Preparation
#----------------------------------------------------------------------------------------------------------
CCdata2$Male <- ifelse(CCdata2$SEX == 'Male', 1, 0)
CCdata2$Female <- ifelse(CCdata2$SEX == 'Female', 1, 0)
CCdata2$White <- ifelse(CCdata2$RACE == 'WHITE', 1, 0)
CCdata2$Black <- ifelse(CCdata2$RACE == 'BLACK OR AFRICAN AMERICAN', 1, 0)
CCdata2$Asian <- ifelse(CCdata2$RACE == 'ASIAN', 1, 0)
CCdata2$Indian <- ifelse(CCdata2$RACE == 'AMERICAN INDIAN OR ALASKA NATIVE', 1, 0)
CCdata2$T1 <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1', 1, 0)
CCdata2$T1a <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1a', 1, 0)
CCdata2$T1b <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T1b', 1, 0)
CCdata2$T2 <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2', 1, 0)
CCdata2$T2a <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2a', 1, 0)
CCdata2$T2b <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T2b', 1, 0)
CCdata2$T3 <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3', 1, 0)
CCdata2$T3a <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3a', 1, 0)
CCdata2$T3b <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3b', 1, 0)
CCdata2$T3c <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T3c', 1, 0)
CCdata2$T4 <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'T4', 1, 0)
CCdata2$TX <- ifelse(CCdata2$AJCC_TUMOR_PATHOLOGIC_PT == 'TX', 1, 0)
CCDatawD2 <- CCdata2 %>% select(-c("SEX","RACE","AJCC_TUMOR_PATHOLOGIC_PT")) #Clinical data with dummy variables
CCDataR2 <- CCDatawD2 %>% select(-c("RISK")) #Clinical Data no Risk
CCdf2 <- CCDatawD2 %>% select(-c("PATIENT","RISK")) #Clinical Dataframe
CCfeatureNames2 <- paste(colnames(CCdf2), collapse = " + ") #Clinical Feature Names
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Survival Formula
#----------------------------------------------------------------------------------------------------------
CCsurvival_formula2 <-paste('Surv(MONTHS,STATUS) ~', CCfeatureNames2)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Splitting Data
#----------------------------------------------------------------------------------------------------------
set.seed(12)
#set.seed(13)
#set.seed(14)
#set.seed(15)
#set.seed(16)
CCDataR2 <- na.omit(CCDataR2)
CCsample2 <- sample.split(CCDataR2$PATIENT, SplitRatio = 0.3)
CCtrain2  <- subset (CCDataR2%>% select(-c("PATIENT")), CCsample2 == FALSE)
CCtest2   <- subset(CCDataR2%>% select(-c("PATIENT")), CCsample2 == TRUE) #this was false

# Split for RSF Kaplan Meier
# Note: Use the same seed for previous split
CCsample2_withRISK <- sample.split(CCDatawD2$PATIENT, SplitRatio = 0.3)
CCtrain2_withRISK  <- subset (CCDatawD2%>% select(-c("PATIENT")), CCsample2_withRISK == FALSE)
CCtest2_withRISK   <- subset(CCDatawD2%>% select(-c("PATIENT")), CCsample2_withRISK == TRUE) #this was false
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#GRID SEARCH FOR OPTIMUM nodesize and mtry
#----------------------------------------------------------------------------------------------------------
#set default values  start again 
nodesize <- c(5,10,20)
nsplit <- c(5,10)

#model to solve for optimum values using 1000 trees
CCnodemodel2 <- tune(Surv(MONTHS, STATUS) ~ ., data = CCtrain2,
                    mtryStart = 2,
                    nodesizeTry= nodesize,
                    cv.fold = 10,
                    nsplitTry = 10, ntree=1000,
                    blocksize=1, importance=TRUE, seed=8)#8 

CCgrid2 <- as.data.frame(CCnodemodel2$results)
CClow2 <- (CCgrid2[CCgrid2$err == min(CCgrid2$err), ])

#Optimized Values
CCnodesize2 <- CClow2$nodesize
CCnodesize2_ <- as.data.frame(CCnodesize2)
CCnodesize2 <- CCnodesize2_[2,1]

CCmtry2 <- CClow2$mtry
CCmtry2_ <- as.data.frame(CCmtry2)
CCmtry2 <- CCmtry2_[2,1]
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## Towards an optimal number of trees RFS
# https://github.com/pedroconcejero/survival/blob/master/survival_random_forests_churn.rmd
#----------------------------------------------------------------------------------------------------------
set.seed(9)
CC.out.rsf.3.2 <- rfsrc( Surv(MONTHS, STATUS) ~ . , 
                        data = CCtrain2, 
                        mtry = CCmtry2, #sqrt(196), 
                        ntree = 1000, cv.fold = 10, verbose = TRUE, block.size=1,
                        importance=TRUE, nodesize = CCnodesize2, #6, 
                        nsplit = 10, seed = 3)#3  13 16

CCerror2 <- gg_error(CC.out.rsf.3.2)
CClowest2 <- min(CCerror2$error[CCerror2$error!=min(CCerror2$error)] )
CCrow2 <- CCerror2[which(CCerror2 == CClowest2),]
plot(gg_error(CC.out.rsf.3.2))
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
#Training RF and Plotting AUC
#----------------------------------------------------------------------------------------------------------
set.seed(3)
#set.seed(4)
#set.seed(5)
#set.seed(6)
#set.seed(7)


CCnumtree2 = CCrow2$ntree[1]

#Use this if you want to manually test specific parameters
CCnumtree2 = 265#6
CCmtry2 = 4
CCnodesize2 = 20

#RANDOM SURVIVAL FOREST
CCsurvival_modelrfsrc2<-rfsrc(Surv(MONTHS, STATUS) ~ ., data = CCtrain2, mtry = CCmtry2,#sqrt(196), 
                             ntree = CCnumtree2, 
                             cv.fold = 10, verbose = TRUE, block.size=1,
                             importance=TRUE, nodesize = CCnodesize2, 
                             nsplit = 10, seed = 4)#4  11, 16 GOOD 

set.seed(3)
CCsurv.fit11.2 <- survivalROC(CCtrain2$MONTHS, CCtrain2$STATUS, CCsurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 11, method = "KM")
CCsurv.fit16.2 <- survivalROC(CCtrain2$MONTHS, CCtrain2$STATUS, CCsurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 16, method = "KM")
CCsurv.fit20.2 <- survivalROC(CCtrain2$MONTHS, CCtrain2$STATUS, CCsurvival_modelrfsrc2$predicted.oob, 
                             predict.time = 20, method = "KM")

plot(CCsurv.fit11.2$FP, CCsurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n", "AUC = "),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",CCnumtree2," Trees Clinical Data \n Parameters: mtry (", CCmtry2,
                             ") nodesize (", CCnodesize2,")"))
lines(CCsurv.fit16.2$FP, CCsurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
lines(CCsurv.fit20.2$FP, CCsurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n", "AUC = "),
      ylab="TP")
abline(0,1)
CCvariable11.2 <- paste0("11 month AUC =", round(CCsurv.fit11.2$AUC,3))
CCvariable16.2 <- paste0("16 month AUC =", round(CCsurv.fit16.2$AUC,3))
CCvariable20.2 <- paste0("20 month AUC =", round(CCsurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(CCvariable11.2,CCvariable16.2,CCvariable20.2))
#----------------------------------------------------------------------------------------------------------
#RANDOM SURVIVAL FOREST TEST

CCpred2 <- predict(CCsurvival_modelrfsrc2,CCtest2)

set.seed(3)
CCsurv.fit11.2 <- survivalROC(CCtest2$MONTHS, CCtest2$STATUS, CCpred2$predicted, 
                             predict.time = 11, method = "KM")
CCsurv.fit16.2 <- survivalROC(CCtest2$MONTHS, CCest2$STATUS, CCpred2$predicted, 
                             predict.time = 16, method = "KM")
CCsurv.fit20.2 <- survivalROC(CCtest2$MONTHS, CCtest2$STATUS, CCpred2$predicted, 
                             predict.time = 20, method = "KM")

plot(CCsurv.fit11.2$FP, CCsurv.fit11.2$TP, type = "l", col = "darkgoldenrod1", xlim = c(0,1), ylim = c(0,1),
     xlab = paste("FP",  "\n"),
     ylab="TP",main = paste0("Area under the curve (AUC) of a time-dependent receiver operating curve (ROC) \n",CCnumtree2," Trees Clinical Data \n Parameters: mtry (", CCmtry2,
                             ") nodesize (", CCnodesize2,")"))
lines(CCsurv.fit16.2$FP, CCsurv.fit16.2$TP, type = "l", col = "deeppink", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
lines(CCsurv.fit20.2$FP, CCsurv.fit20.2$TP, type = "l", col = "deepskyblue", xlim = c(0,1), ylim = c(0,1),
      xlab = paste("FP",  "\n"),
      ylab="TP")
abline(0,1)
CCvariable11.2 <- paste0("11 month AUC =", round(CCsurv.fit11.2$AUC ,3))
CCvariable16.2 <- paste0("16 month AUC =", round(CCsurv.fit16.2$AUC,3))
CCvariable20.2 <- paste0("20 month AUC =", round(CCsurv.fit20.2$AUC,3))
legend("bottomright",fill=c("darkgoldenrod1","deeppink", "deepskyblue"),c(CCvariable11.2,CCvariable16.2,CCvariable20.2))

#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Kaplan Meier Curve
#----------------------------------------------------------------------------------------------------------
CCfit2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = CCDatawD2)
print(CCfit2)
summary(CCfit2)$table
ggsurvplot(CCfit2,
           data = CCDatawD2,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Predicting Error Rate
#https://stackoverflow.com/questions/61478162/depth-and-oob-error-of-a-randomforest-and-randomforestsrc
#----------------------------------------------------------------------------------------------------------
set.seed(5)
CCpred2 <- predict(CCsurvival_modelrfsrc2,CCtest2[, (colnames(CCtest2))],block.size=1)
plot(CCpred2$err.rate,type="l",col="orange",xlab="Number of trees",ylab="Error rate",
     ylim=c(0,0.5))

legend("topright",fill=c("orange"),c("OOB.train"))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#Getting Cindex
#----------------------------------------------------------------------------------------------------------
get.cindex(time = CCtrain2$MONTHS, censoring = CCtrain2$STATUS, predicted = CCsurvival_modelrfsrc2$predicted.oob)
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
## obtain Brier score using KM and RSF censoring distribution estimators
#----------------------------------------------------------------------------------------------------------
CC.bs.km.2 <- get.brier.survival(CCsurvival_modelrfsrc2, cens.mode = "km")$brier.score
CC.bs.rsf.2 <- get.brier.survival(CCsurvival_modelrfsrc2, cens.mode = "rfsrc")$brier.score
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
## plot the brier score
#----------------------------------------------------------------------------------------------------------
plot(CC.bs.km.2, type = "s", col = 2)
lines(CC.bs.rsf.2, type ="s", col = 4)
legend("bottomright", legend = c("cens.model = km", "cens.model = rfsrc"), fill = c(2,4))
#----------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------
#KAPLAN MEIER of RSF Results
#----------------------------------------------------------------------------------------------------------
CCRun2 <- CCtest2
CCRun2_withRISK <- CCtest2_withRISK 

CCsurv_prob2 <- predict(CCsurvival_modelrfsrc2, newdata = CCRun2, type = "survival" ) #type is the same thing as "response" 
CCsurv_times2 <- survfit(Surv(MONTHS, STATUS) ~ RISK, data = CCRun2_withRISK, weights = CCsurv_prob2$survival[,11] )
#plot(surv_times, main = "Kaplan-Meier Plot", xlab = "Time", ylab = "Survival Probability")

ggsurvplot(CCsurv_times2,
           data = CCRun2_withRISK,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#----------------------------------------------------------------------------------------------------------

