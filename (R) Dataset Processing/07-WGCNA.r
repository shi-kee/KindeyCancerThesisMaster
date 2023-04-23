#install.packages("CorLevelPlot")
#install.packages("remotes")
#remotes::install_github("kevinblighe/CorLevelPlot")
install.packages("rasterVis")




library(WGCNA)
library(DESeq2)
#library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(rasterVis)

allowWGCNAThreads()   #allow multi-threading (optional)

# PREPARE DATA
options(stringsAsFactors = FALSE);

data0=read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/FILTERED_data_mrna_seq_v2_rsem_.csv")
data1 <- data0[,-1]
#rownames(data1) <- data0[,1]
datExpr0 = as.matrix(data0)
datExpr1 = as.data.frame(t(datExpr0))#rows to columns
datExpr2 = as.matrix(datExpr1)

#install.packages("janitor")
library(janitor)
duplicated(datExpr2[1,])
datExpr3 <- row_to_names(datExpr2,row_number=1)

# Normalization with log2(FPKM+1)
sample_metadata = read.csv(file = "C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/PATIENT_FILTERED.csv")


datExpr3 = as.data.frame(datExpr3)#rows to columns; change to dataframe to use goodsamplegenes
#################################################################

#QUALITY CONTROL - Detecting outlier samples
gsg <- goodSamplesGenes(datExpr3) 
#gsg variable holds logical vectors and shows whether there are any genes or samples
#that are detected to be outliers, so if you want to evaluate whether there are any outliers 
#detected in your data, then you need to look at the allOK vector
summary(gsg)
gsg$allOK 
#Since this is false, then it indicates either or both the genes or the samples
#are detected as outliers and we need to find and exclude them from our analysis

table(gsg$goodGenes) #this many genes are detected as outliers
table(gsg$goodSamples) #all samples have passed the filters and are good to go

#remove genes that are detected as outliers
data_goodgenes <- datExpr3[,gsg$goodGenes == TRUE] #dataframe with good genes

#detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(data_goodgenes), method = "average")
plot(htree)

#PCA - method 2
pca <- prcomp(data_goodgenes)
pca.dat <- pca$x #principal components computed for all samples

pca.var <- pca$sdev^2#calculate the variance explained by each principal component
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2) #percentage variance

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))


#Exclude outlier samples
samples.to.be.excluded <- c('TCGA.J7.A8I2.01','TCGA.Y8.A8S1.01','TCGA.BQ.7049.01','TCGA.2Z.A9JJ.01')
datExpr4 <- t(datExpr3)
data.subset <- datExpr4[,!(colnames(datExpr4) %in% samples.to.be.excluded)]
data.subset <- as.data.frame(t(data.subset)) #removed samples

#PCA - show results of removal
pca <- prcomp(data.subset)
pca.dat <- pca$x #principal components computed for all samples

pca.var <- pca$sdev^2#calculate the variance explained by each principal component
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2) #percentage variance

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], '%'),
       y = paste0('PC2: ', pca.var.percent[2], '%'))

#Checking for batch effects
#install.packages("rgl")
library(rgl) #to view samples in 3D pca
plot3d(pca.dat[,1:3], 
       size=5,
       col = seq(nrow(pca.dat)))

text3d(pca.dat[,1:3],
       texts=c(rownames(pca.dat)), 
       cex= 0.7, pos=3)

#############
#normalization
data.normalize <- as.data.frame(t(data.subset)) 
data.log = log(data.normalize+1) 
head(data.log[1:5,1:5])
rownames(data.log) <- substr(rownames(data.log),2,nchar(rownames(data.log)))#remove X in beginning

sample_metadata0 <-column_to_rownames(sample_metadata,var = 'PATIENT') #change column Patient as row index
colData <- sample_metadata0 %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

#checking rownames and column names identical
all(rownames(colData) %in% colnames(data.log))
all(rownames(colData) == colnames(data.log))
data.log <- as.data.frame(t(data.log))
################################################################################

#NETWORK CONSTRUCTION
#Choose a set of soft-thresholding powers
'''power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

###test
test <- t(data.log) #[1:5000,]


###

#Call the network topology analysis function
sft <- pickSoftThreshold(test, #changed this from test
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices #use 2 metrics to pick power - rsquared values(max) and mean connectivity (min)

#visualization to pick power
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power' , y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2) #12
'''#I did not use this part since it gave wonky results
#---------
# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))

sft = pickSoftThreshold(data.log, powerVector = powers, verbose = 5) #this had test before
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


# Turn data expression into topological overlap matrix
power=sft$powerEstimate #4

################################################################################

#convert matrix to numeric 
test <- data.log #copied this from commented section

test[] <-sapply(test, as.numeric)

soft_power <- 12 #did not use this
temp_cor <- cor
cor <-WGCNA::cor

#memory estimate w.r.t. blocksize
bwnet <- blockwiseModules(data.log,
                 maxBlockSize = 21000,
                 TOMType = "signed",
                 power = power, #soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

###############################################################################
#MODULE EIGENGENES

module_eigengenes <- bwnet$MEs

#Print out a preview
head(module_eigengenes)

#get number of genes for each module
table(bwnet$colors)
mergedColors = labels2colors(bwnet$colors)
unmergedColors = labels2colors(bwnet$unmergedColors)

#Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(unmergedColors[bwnet$blockGenes[[1]]], mergedColors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

sizeGrWindow(12, 9)
mergedColors = labels2colors(bwnet$colors)
unmergedColors = labels2colors(bwnet$unmergedColors)
pdf(file = "module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(unmergedColors[bwnet$blockGenes[[1]]],mergedColors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

###############################################################################
#grey module = all genes that does not fall into other modules were assigned to the grey module

#RELATE MODULES TO TRAITS
#module trait associations

#Define numbers of genes and samples
nSamples <- nrow(test)
nGenes <- ncol(test)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(test, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
####write.csv(MEs,"C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/MEs.csv",row.names=TRUE)
#commented out the previous line to not overwrite current file

# Read patient-CellProfiler data as traits
traits = read.csv("C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/PATIENT_FILTERED.csv")
#line below is temporary

rownames(traits) = traits[, 1]
traits = traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
traits = traits[match(rownames(MEs), rownames(traits)), ]
table(rownames(MEs) == rownames(traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, traits, use = 'p');
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
#write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");

#  Plot heatmap of module-traits relationship
heatmap.data <- merge(MEs, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

heatmap.data.new <- heatmap.data %>% select(c("MEsalmon","MEdarkturquoise", "MEblue", "MEturquoise", "MEgreen", "MEbrown", "MEdarkred", "MEpurple", "MEyellow", "AreaShape_Zernike_7_3", "Granularity_3_ErodeImage"))

CorLevelPlot(#heatmap.data,
             #x = names(heatmap.data)[44:45],
             #y = names(heatmap.data)[1:43],
              heatmap.data.new,
              x = names(heatmap.data.new)[10:11],
              y = names(heatmap.data.new)[1:9],
             col = c("royalblue", "skyblue", "white", "pink", "deeppink"),
             rotLabX = 45,
             #rotLabY = 45,
             corFUN = "pearson",
             corUSE = "pairwise.complete.obs", fontLabX = 1, fontLabY = 1, titleY = "Gene Modules",
             titleX = "Image Features", rotTitleY = 90, 
             )



#determine which genes are in specific modules
module.gene.mapping <- as.data.frame(bwnet$colors)

#TEMPLATE SAVING
MEpurple <- module.gene.mapping %>%
  filter(`bwnet$colors`=='purple') %>%
  rownames()

MEpurple <- as.data.frame(MEpurple)
write.csv(MEpurple,"C:/Users/keesh/OneDrive/Documents/Thesis_Life Depends on It/SVM-RFE-master/New WGNCA/MEpurple.csv",row.names=FALSE)

#Saving Turquoise genes
MEturquoise <- module.gene.mapping %>%
  filter(`bwnet$colors`=='turquoise') %>%
  rownames()

MEturquoise <- as.data.frame(MEturquoise)
write.csv(MEturquoise,"MEturquoise.csv",row.names=FALSE)

#Saving Red genes
MEred <- module.gene.mapping %>%
  filter(`bwnet$colors`=='red') %>%
  rownames()

MEred <- as.data.frame(MEred)
write.csv(MEred,"MEred.csv",row.names=FALSE)

#Saving Yellow genes
MEyellow <- module.gene.mapping %>%
  filter(`bwnet$colors`=='yellow') %>%
  rownames()

MEyellow <- as.data.frame(MEyellow)
write.csv(MEyellow,"MEyellow.csv",row.names=FALSE)

#Saving Lightgreen genes
MElightgreen <- module.gene.mapping %>%
  filter(`bwnet$colors`=='lightgreen') %>%
  rownames()

lightgreen <- as.data.frame(lightgreen)
write.csv(lightgreen,"lightgreen.csv",row.names=FALSE)

#Saving Salmon genes
MEsalmon <- module.gene.mapping %>%
  filter(`bwnet$colors`=='salmon') %>%
  rownames()

salmon <- as.data.frame(salmon)
write.csv(salmon,"salmon.csv",row.names=FALSE)




###############################################################################
#INTRAMODULAR ANALYSIS: Identifying driver genes
#identify highly intramodular hub genes in the modules of interest by calculating the correlation of the module eigengenes and gene expression profile 

#Calculate the module membership and the associated p-values
#The module membership/itnramodular connectivity is calculated as the correlation of the eignenges and the gene expression profile
#This quantifies the similarity of all genes on the array to every module

module.membership.measure <- cor(module_eigengenes, test, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

#Calculate the gene significance and associated p-values
gene.signf.corr <- cor(test, traits$AreaShape_Zernike_7_3, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


#see top 25 genes that are sig associated with AreaShape_Zernike_7_3
gene.signf.corr.pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25) 
