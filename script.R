# a function that install and load all the required libraries to the workspace
LoadLibraries= function (){
  install.packages("WGCNA")
  install.packages("BiocManager")
  BiocManager::install("Biobase")
  BiocManager::install("GO.db")
  BiocManager::install("htmlTable")
  BiocManager::install("impute")
  BiocManager::install("preprocessCore")
  BiocManager::valid()
  install.packages("reshape2")
  library(WGCNA)
  library(preprocessCore)
  library("reshape2")
  library(MASS)
  library(ggplot2)
  library(dplyr)
  install.packages("ggbeeswarm")
  library(ggbeeswarm)
  library(ISLR)
  install.packages('devtools')
  install.packages('remotes')
  library(devtools)
  remotes::install_github('vqv/ggbiplot')
  library(ggbiplot)
  library(smbinning)
  install.packages("smbinning")
  library(randomForest)
  install.packages("caret")
  library(caret)


  print("The libraries have been loaded .")
}
LoadLibraries()
#----------------------------------------------------------------------------------------------------------------------------
# read/load the data files:

load(paste0("../GTEx6/","mat.f.coding.RData"), verbose = "TRUE")
load(paste0("../GTEx6/","pheno.f.RData"), verbose = "TRUE")
load(paste0("../GTEx6/","gene.f.RData"), verbose = "TRUE")
load(paste0("../GTEx6/","gene.f.with.entrez.RData"), verbose = "TRUE")
#----------------------------------------------------------------------------------------------------------------------------
# function to get the data for the relevant tissue
ts.list = as.character(unique(pheno.f$SMTSD))

get.raw.tissue.edata<-function(tissue.name, mat.f.coding, pheno.f){
  tiss.cols.1 = which(pheno.f$SMTSD %in% tissue.name)
  mat.1 = mat.f.coding[, tiss.cols.1]
  return(mat.1)
}
#----------------------------------------------------------------------------------------------------------------------------
# at index 36 are the "Lung" tissue that we picked to work on
tmp.tissue = ts.list[36]
print(paste0("loading ", tmp.tissue, " edata"))
reads.src1 = get.raw.tissue.edata(tmp.tissue, mat.f.coding, pheno.f)
# transpose the dataframe
t.reads.src = t(reads.src1)
length(t.reads.src)
#//6340480 length
#----------------------------------------------------------------------------------------------------------------------------
# get only values whose greater than log 0.1
vec.1 = apply(reads.src1 , 1, function(x) length(which( x > log(0.1+1, 2) )))

row.index = which(vec.1 > (0.8*(ncol(reads.src1 ))))
src.reads = reads.src1 [row.index, ]
#----------------------------------------------------------------------------------------------------------------------------

var.data <- apply(src.reads, 1, var) #generate variance of each row - gene
low.var.indxs = which(var.data == 0)
if(length(low.var.indxs) > 0)
{
  data.free = src.reads
  #now we get smaller matrix, with no genes with variance 0
  src.reads <- data.free[-low.var.indxs,]
}
length(src.reads)
src.reads
#//4705920
#//deleted : 1634560
#//14703 rows
#----------------------------------------------------------------------------------------------------------------------------

remove.outliers.with.SD<-function(t.reads.src)
{
  #remove outliers
  #cluster the samples and not the genes to find outliers
  A = adjacency(t(t.reads.src), type = "distance")
  #the connectivity of each human. -1 is to remove the diagonal, the cor to itself
  k = as.numeric(apply(A,2,sum))-1
  Z.k = scale(k) #standardized k
  thresholdZ.k = -3 #standard deviation
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")#the red is the outlier
  my.outliers = which(outlierColor == "red")
  #printing the outlier samples
  my.outliers.samp = (rownames(t.reads.src))[my.outliers]
  print("outlier samples to remove")
  print(my.outliers.samp)
  my.good.samples = which(outlierColor == "black")
  my.good.samples.names = (rownames(t.reads.src))[my.good.samples]
  #printing the outlier samples
  #print(my.good.samples.names)
  #this is the final mat after outliers removal
  t.reads.src = t.reads.src[my.good.samples.names, ]
  return(t.reads.src)
}

#----------------------------------------------------------------------------------------------------------------------------

t.reads.src = remove.outliers.with.SD(t(src.reads))
tissue.edata = t(t.reads.src)
tissue.edata
sampleTree = hclust(dist(t(tissue.edata)), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.3);#change this to change the size of the text
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

#----------------------------------------------------------------------------------------------------------------------------
#rows are genes, columns are samples
quantile.normalize.raw.gtex <- function(edata.mat)
{
  norm_edata = normalize.quantiles(as.matrix(edata.mat))
  rownames(norm_edata) = rownames(edata.mat)
  colnames(norm_edata) = colnames(edata.mat)
  return(norm_edata)
}
tissue.edata.qn = quantile.normalize.raw.gtex(tissue.edata)
#----------------------------------------------------------------------------------------------------------------------------
# plotting normalized data
boxplot(tissue.edata.qn)
#----------------------------------------------------------------------------------------------------------------------------
# plotting raw data
boxplot(tissue.edata)
#----------------------------------------------------------------------------------------------------------------------------
boxplot(src.reads)
#----------------------------------------------------------------------------------------------------------------------------
# a function to get thew relevant samples data that belongs to the lung tissue
get.pheno.mat.from.samplIDs<-function(pheno.mat,samples.vec){
  indexes=which(pheno.f$SAMPID %in% samples.vec)
  tissue.pheno=pheno.f[indexes,]
  ordering=match(tissue.pheno$SAMPID,samples.vec)
  tissue.pheno=tissue.pheno[ordering,]
  if(!identical(as.character(tissue.pheno$SAMPID),samples.vec)) {print("ERROR 2: samples in norm.mat and pheno do not match")}
  return(tissue.pheno)
}

samples.id=colnames(tissue.edata)
tmp.pheno=get.pheno.mat.from.samplIDs(pheno.f,samples.id)
#----------------------------------------------------------------------------------------------------------------------------
# transposing the matrix (switching between columns and rows)
tissue.edata.qn.t=t(tissue.edata.qn)
#----------------------------------------------------------------------------------------------------------------------------
# to get rid of big dataframes only work with the relevant columns
relevant.features.pheno=select(tmp.pheno,SUBJID,SAMPID,SMRIN,SMTSISCH,SMGEBTCH,AGE,DTHHRDY)
#----------------------------------------------------------------------------------------------------------------------------
# getting the relevant genes data from the genes.f data (maybe extra)
genes.id=colnames(tissue.edata.qn.t)
gene.f=as.data.frame(gene.f)
indexs=which(gene.f$Name %in% genes.id)
releveant.genes=gene.f[indexs,]
#----------------------------------------------------------------------------------------------------------------------------
#plotting the ages as a bar chart
ggplot(relevant.features.pheno, aes(x=AGE)) +
  geom_bar()
#----------------------------------------------------------------------------------------------------------------------------
#plotting the death cause as a bar chart

ggplot(relevant.features.pheno, aes(x=DTHHRDY)) +
  geom_bar()

#----------------------------------------------------------------------------------------------------------------------------0
ages.lbl=c("20-29","30-39","40-49","50-59","60-69","70-79")
# ages pie chart for the ventilator cases
inx=which(relevant.features.pheno$DTHHRDY %in% 0)
agg.ages=(relevant.features.pheno[inx,] %>% group_by(AGE) %>% summarise(counter = n()))[2]

agg.ages=c(21,14,45,68,40,2)

piepercent<- round(100*agg.ages/sum(agg.ages), 1)

pie(agg.ages, labels =piepercent, main = "Ventilator Cases Ages",col = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))
legend("topright", c("20-29","30-39","40-49","50-59","60-69","70-79"), cex = 0.8,fill = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))

#----------------------------------------------------------------------------------------------------------------------------1
# ages pie chart for the Violent death cases

inx=which(relevant.features.pheno$DTHHRDY %in% 1)

agg.ages=(relevant.features.pheno[inx,] %>% group_by(AGE) %>% summarise(counter = n()))[2]

agg.ages=c(2,3,4,5,3,0)

piepercent<- round(100*agg.ages/sum(agg.ages), 1)
pie(agg.ages, labels =piepercent , main = "Violent Death Ages",col = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))

#----------------------------------------------------------------------------------------------------------------------------2
# ages pie chart for the fast death cases

inx=which(relevant.features.pheno$DTHHRDY %in% 2)
agg.ages=(relevant.features.pheno[inx,] %>% group_by(AGE) %>% summarise(counter = n()))[2]
agg.ages=c(2,2,6,23,26,3)
piepercent<- round(100*agg.ages/sum(agg.ages), 1)


pie(agg.ages, labels =piepercent , main = "Fast Death Ages",col = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))

#----------------------------------------------------------------------------------------------------------------------------3
# ages pie chart for the intermediate death cases

inx=which(relevant.features.pheno$DTHHRDY %in% 3)

agg.ages=(relevant.features.pheno[inx,] %>% group_by(AGE) %>% summarise(counter = n()))[2]


agg.ages=c(0,0,0,4,13,0)
piepercent<- round(100*agg.ages/sum(agg.ages), 1)

pie(agg.ages, labels =piepercent , main = "Intermidiate Death Ages",col = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))

#----------------------------------------------------------------------------------------------------------------------------4
# ages pie chart for the slow death cases

inx=which(relevant.features.pheno$DTHHRDY %in% 4)

agg.ages=(relevant.features.pheno[inx,] %>% group_by(AGE) %>% summarise(counter = n()))[2]

agg.ages=c(0,2,3,9,13,0)

piepercent<- round(100*agg.ages/sum(agg.ages), 1)
pie(agg.ages, labels =piepercent , main = "Slow Death Ages",col = c("azure","bisque1","cadetblue1","chartreuse1","chocolate1","brown1"))

#----------------------------------------------------------------------------------------------------------------------------
# plotting age groups the RNA integrity feature
ggplot(data=relevant.features.pheno)+
  aes(y=SMRIN,x=AGE)+
  geom_beeswarm(col="#72CB5C")
#----------------------------------------------------------------------------------------------------------------------------
# RNA integrity to time delay for each sample , while coloring according to the death cause
ggplot(data=relevant.features.pheno, aes(y=SMRIN, x=SMTSISCH)) +
  geom_point(cex = 1.5,aes(colour=DTHHRDY))
#----------------------------------------------------------------------------------------------------------------------------
# RNA integrity to age groups and coloring by age
ggplot(data=relevant.features.pheno, aes(y=SMRIN, x=AGE)) +
  geom_point(cex =2,aes(colour=AGE))
#----------------------------------------------------------------------------------------------------------------------------
# generate a pca model
genes.pca <- prcomp(tissue.edata.qn.t, center = TRUE, scale. = TRUE)
summary(genes.pca)
genes.pca.df=as.data.frame(genes.pca$x)
#----------------------------------------------------------------------------------------------------------------------------
# merging between the models (genes data the samples for each person)
merged.raw=merge(x=tissue.edata.qn.t,y=relevant.features.pheno,by.x=c('row.names'),by.y=c("SAMPID"))
merged.pca=merge(x=genes.pca.df,y=relevant.features.pheno,by.x=c('row.names'),by.y=c("SAMPID"))
#----------------------------------------------------------------------------------------------------------------------------
# casting the column for plotting purposes
merged.pca$DTHHRDY=as.character(merged.pca$DTHHRDY)
#----------------------------------------------------------------------------------------------------------------------------
# plotting the PC1 and PC2 for each feature
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 1.5,aes(colour=AGE))+scale_colour_brewer(palette = "Set1")
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 1.5,aes(colour=DTHHRDY))+scale_colour_brewer(palette = "Set1")
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 1.5,aes(shape=DTHHRDY,colour=AGE))+scale_colour_brewer(palette = "Set1")
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 1.5,aes(shape=AGE,colour=DTHHRDY))+scale_colour_brewer(palette = "Set1")
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 2,aes(colour=SMRIN))
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 2,aes(colour=SMTSISCH))
ggplot(merged.pca, aes(x=PC1, y=PC2)) + geom_point(cex = 1.5,aes(colour=SMGEBTCH))
#----------------------------------------------------------------------------------------------------------------------------
# linear regression
l =lm(PC1~PC2,data=merged.pca)
ggplot(l,aes(x=PC1,y=PC2))+geom_point()+geom_abline(col="red")
#----------------------------------------------------------------------------------------------------------------------------
cor.test(as.numeric(merged.pca$PC1),as.numeric(merged.pca$PC2),exact=False)
#----------------------------------------------------------------------------------------------------------------------------
# data split by 80-20
smp_size <- floor(0.80 * nrow(merged.pca))

## set the seed to make your partition reproducible
set.seed(123)
train_ind <- sample(seq_len(nrow(merged.pca)), size = smp_size)

train <- merged.pca[train_ind, ]
test <- merged.pca[-train_ind, ]
# fill null values with zeros (constraint to get the algorithm running)
train[is.na(train)] <- 0
test[is.na(test)] <- 0

#----------------------------------------------------------------------------------------------------------------------------
# Logistic Regression to predict the SMRIN feature
lm.fit.result =glm(SMRIN~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=train)
model=data.frame(name=test$Row.names,prediction=predict(lm.fit.result,newdata=test),actual=test$SMRIN)
ggplot(model,aes(x=prediction,y=actual))+geom_point()+geom_abline(col="red")

#----------------------------------------------------------------------------------------------------------------------------
# Logistic Regression to predict the SMTSISCH feature

lm.fit =glm(SMTSISCH~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20,data=train)
model2=data.frame(name=test$Row.names,prediction=predict(lm.fit,newdata=test),actual=test$SMTSISCH)
ggplot(model2,aes(x=prediction,y=actual))+geom_point()+geom_abline(col="red")

#----------------------------------------------------------------------------------------------------------------------------
# function to normalize the two different model results (for comparsion)
normalized = function(x){
  return (x-min(x))/(max(x)-min(x))
}

d=data.frame(normalized=normalized(model$prediction)/max(model$prediction),actual=normalized(model$actual)/max(model$actual))
d1=data.frame(normalized=normalized(model2$prediction)/max(model2$prediction),actual=normalized(model2$actual)/max(model2$actual))

# calculate the error/mean distance for each model
sum(abs(d$normalized-d$actual))/nrow(d)
sum(abs(d1$normalized-d1$actual))/nrow(d1)

#----------------------------------------------------------------------------------------------------------------------------
# Random Forest to predict DTHHRDY
set.seed(1)
# converting the columns/feature to categorical
train$DTHHRDY=as.factor(train$DTHHRDY)
test$DTHHRDY=as.factor(test$DTHHRDY)

forest=randomForest(y=train$DTHHRDY,x=train[2:317],subset = train,ntree=1000,mtry=sqrt(316),importance=TRUE)

# [2:317] are all the PC# values in the dataframe
y.pred=predict(forest,newdata = test[2:317])

# generating a confusion matrix to see the predicted vs actual results
confusion.matrix=table(y.pred,test$DTHHRDY)
confusionMatrix(confusion.matrix)
#----------------------------------------------------------------------------------------------------------------------------
# Random Forest to predict AGE

train$AGE=as.factor(train$AGE)
test$AGE=as.factor(test$AGE)

forest2=randomForest(y=train$AGE,x=train[2:317],subset = train,ntree=1000,mtry=sqrt(316),importance=TRUE)

# [2:317] are all the PC# values in the dataframe
y.pred2=predict(forest2,newdata = test[2:317])

# generating a confusion matrix to see the predicted vs actual results
confusion.matrix2=table(y.pred2,test$AGE)
confusionMatrix(confusion.matrix2)
#----------------------------------------------------------------------------------------------------------------------------
