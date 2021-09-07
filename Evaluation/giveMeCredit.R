
#################################################################################################################################
## Evaluation on the dataset GIVE ME SOME CREDIT KAGGLE 
## https://www.kaggle.com/c/GiveMeSomeCredit/data
#################################################################################################################################


######################################
library(robustbase)
library(pROC)
library(glm)
library(rrcovHD)
library(MASS) #cov.mcd
#####################################

source("glmrobBy_UPDATE.R")
source("helpFunctions.R")

givemecredit <- read.csv("GiveMeSomeCredit/cs-training.csv")
str(givemecredit)

sum(complete.cases(givemecredit)) ## 120269 / 150000
givemecredit <- givemecredit[complete.cases(givemecredit), ]
summary(as.factor(givemecredit$SeriousDlqin2yrs))/nrow(givemecredit)

results <- calculate_estimatorsGMC()
# 0.3932158 0.5994681 0.1383382 0.1139128 0.5667799 0.6281824 
names(results) <- c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY", "WBY cost-sensitive") 


##################################################################################################################################
## introduce 10% of outliers in the majority class 
##################################################################################################################################

set.seed(100)
missclasified <- sample(which(givemecredit$SeriousDlqin2yrs == 0), size = round(nrow(givemecredit)/10))
givemecredit[missclasified, "SeriousDlqin2yrs"] <- 1

resultsM <- calculate_estimatorsGMC()
resultsM

resulting.table <- data.frame(cbind(results, resultsM), row.names = c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY", "WBY cost-sensitive"))
colnames(resulting.table) <- c("original", "10% missclasified")
