#Name : Abdallah Osama Mohamed 
#ID : 20198053
#Group:B1
######################
#Name :Amr Hosney Eid 
#ID: 20198059
#Group:B1
#########################
#TA : Mohamed Ramadan
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("antiProfilesData")
library(Biobase)
library("Hmisc")
library(corrplot)
#suitable normalization.
library(MASS)
#assign = each of feature, phenotype and expression data to different variables
apColonData<-antiProfilesData::apColonData
pdata=pData(apColonData)
edata=exprs(apColonData)
fdata = fData(apColonData)
##############################################################################
#1.a Show the type of each column

type_column<-sapply(pdata,class)
print(type_column)
###############################################################################
#1.b Show column names and rows name
column_names<-colnames(pdata)
print(column_names)
row_names<-rownames(pdata)
print(row_names)
################################################################################
#1.c Calculate summary of each column
summary_column<-summary(edata)
print(summary_column)
################################################################################
#1.d Show frequency of categorical data, taking into the consideration, NA values frequency if any.
frequency_categorical<-table(pdata$Status,useNA="ifany")
##############################################################################
#1.e Calculate the correlation and covariance between the first 10 columns only of our data set and draw full correlation matrix.
my_data <- edata[, c(1:10)]
mydata.cov<-cov(my_data, y = NULL)
print(mydata.cov)
mydata.cor<-cor(my_data)
print(mydata.cor)
corrplot(mydata.cor,tl.col = "black")
################################################################################
#1.f For both genes: GSM95478,GSM95473 show the plot with a line of their relation.
qplot(edata[,"GSM95478"],edata[,"GSM95473"],geom="point") + geom_smooth(method = "lm",se=FALSE)
################################################################################
#2 Using PCA and SVD, Prove by plotting and values that both can return the same result by
#suitable normalization.
library(MASS)
PCA = prcomp(edata)
pcvar<-PCA$sdev^2
pcavarper <- round(pcvar/sum(pcvar)*100, 1)
barplot(pcavarper, main="PCA PLOT", xlab="Principal Component", ylab="Percent Variation")
table(pcavarper)
###############################################################################
svd1 = svd(edata)
svdvarper <- round(svd1$d^2/sum(svd1$d^2)/sum(svd1$d^2/sum(svd1$d^2))*100, 1)
table(svdvarper)
barplot(svdvarper,ylab="Percent Variance Explained", main="SVD PLOT",col=2)
plot(PCA$rotation[,1],svd1$v[,1])
#############################################################################
#3
observed <-c (29, 24, 22, 19, 21, 18, 19, 20, 23, 18, 20, 23)
#p_value <- pchisq(q = sum((observed - (256/12))^2 / (256/12)), df = 12 - 1,lower.tail = FALSE)
chisq.test(observed)

################################################################################
#Hypothesis H0 is uniformly births distribution where H1 is not 
#4
col10=edata[,1:10]
dist1 = dist(t(col10))
hclust1 = hclust(dist1)
plot(hclust1,hang = -1) # hang = - 1 to make labels written on the same level
#******************#
kmeanss = kmeans(edata,8)
kmeanss$centers