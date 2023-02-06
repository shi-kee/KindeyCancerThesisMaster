# Copyright (C) 2011  John Colby
# http://github.com/johncolby/SVM-RFE

# Start up an R session in the SVM-RFE directory. Then work through these commands.

# Set up R environment
set.seed(7777777) #Three features intersect with Lasso results
#set.seed(147258)
library(e1071)
source('msvmRFE.R')
load('demo/input.Rdata')
DF=read.csv(file= 'FinalDataset_NEW1.csv')
DF

# Take a look at the expected input structure
dim(DF)
DF[1:5,1:5]

# Basic usage
model <- svmRFE(DF, k=5, halve.above=100)#changed

# Set up cross validation
nfold = 5
nrows = nrow(DF)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds
folds = lapply(1:nfold, function(x) which(folds == x))
folds

# Perform feature ranking    on all training sets
results = lapply(folds, svmRFE.wrap, DF, k=5, halve.above=100)
length(results)
results

# Obtain top features across ALL folds
top.features = WriteFeatures(results, DF, save=F)
head(top.features)

########


########

# Estimate generalization error using a varying number of top features
featsweep = lapply(1:5, FeatSweep.wrap, results, DF)
featsweep

# Make plot
no.info = min(prop.table(table(DF[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=4, height=4, bg='white')
PlotErrors(errors, no.info=no.info)
dev.off()


