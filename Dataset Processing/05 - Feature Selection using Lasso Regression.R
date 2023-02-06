#code adapted from https://www.youtube.com/watch?v=5GZ5BHOugBQ
#github version https://github.com/davidcaughlin/R-Tutorial-Data-Files

#loading data
df=read.csv(file= 'FinalDataset_NEW1.csv')
df$Target = as.factor(df$Target)

str(df)

#statistical assumptions

####Partitioning the Data
#install.packages ("caret")
library(caret)

######Fix error (Confirmation if caret is loaded)
loaded_packages  <- library()$results[,1]

# confirm caret is loaded
"caret" %in% tolower(loaded_packages)
######

#Set seed
set.seed(1985)

#Partition (split) and create index matrix of selected values
index <- createDataPartition(df$Target, p=.7, list=FALSE, times=1)

#Create test and training data frames
train_df <- df[index,]
test_df <- df[-index,]

#k-fold cross-validation (5-fold cross-validation) framework

#Specify 5-fold cross-validation as training method (framework)
ctrlspecs <- trainControl(method="cv", number=5,
                          savePredictions = "all")

####Specify & Train LASSO Regression Model

#Create vector of potential lambda values
lambda_vector <- 10^seq(5, -5, length=500)

#Set seed
set.seed(1985)

#Specify LASSO regression model to be estimated using training data
#and 5-fold cross-validation framework/process
model1 <- train(Target ~ .,
                data=train_df,
                preProcess=c("center","scale"),
                method="glmnet",
                tuneGrid=expand.grid(alpha=1, lambda=lambda_vector),
                trControl=ctrlspecs,
                na.action=na.omit)

#Best (optimal) tuning parameter (alpha, lambda)
model1$bestTune
model1$bestTune$lambda

#LASSO regression model coefficients (parameter estimates)
round(coef(model1$finalModel, model1$bestTune$lambda),3)



#Plot log(lambda) & RMSE
plot(log(model1$results$lambda),
     model1$results$Accuracy,
     xlab="log(lambda)",
     ylab="Accuracy"
     #xlim=c(-10,-2)
     )

log(model1$bestTune$lambda)

#Variable importance
library(caret)
varImp(model1)

########

#Data visualization of variable importance
#install.packages("ggplot2")
library(ggplot2)
varImp(model1)
ggplot(varImp(model1))
###


#####Model Prediction
predictions1 <- predict(model1, newdata=test_df)

#Model performance/accuracy
tbl <-confusionMatrix(predictions1, test_df$Target)
tbl

##saving varImp
str(varImp(model1)$importance)
imps <- as.matrix(varImp(model1)$importance)
imps

data_zero <- imps[apply(imps, 1, function(row) all(row !=0 )), ]  # Remove zero-rows
View(data_zero) 


LASSO_data <-as.data.frame(data_zero)
LASSO_data
library(MASS)
write.csv(LASSO_data,file="LASSO_results.csv")
