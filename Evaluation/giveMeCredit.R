
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


######################################################################################################################################
## Univariate distribution shows extreme skeweness of the majority of inputs ==> log transformation

#######################################################################################################################################
## variables to transform 
#######################################################################################################################################

vars.to.transform <- names(givemecredit)[-c(1, 3, 7, 11)]
LOGgivemecredit <- givemecredit

for (var in vars.to.transform) {
  LOGgivemecredit[, var] <- log(givemecredit[, var] +1)
}

str(LOGgivemecredit)
names(LOGgivemecredit)[-c(1, 3, 7, 11)] <- paste0("LOG", names(LOGgivemecredit)[-c(1, 3, 7, 11)])

########################################################################################################################################
## results
#########################################################################################################################################

results <- calculate_estimatorsGMC(givemecredit)
names(results) <- c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY", "WBY cost-sensitive") 

resultsLOG <- calculate_estimatorsGMC(LOGgivemecredit)
names(results) <- c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY", "WBY cost-sensitive") 




###########################################################################################################################################
## score plot 
###########################################################################################################################################

get_betas <- function(data) {
  model <- glm(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ], family=binomial)
  out1 <- model$coefficients
  model <- glm(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ], 
               weights = get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")] , family=binomial)
  out2 <- model$coefficients
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(data, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  out3 <- model$coefficients
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights =  get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  out4 <- model$coefficients
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(data, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  out5 <- model$coefficients
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights =  get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  out6 <- model$coefficients
  return(list("glm" = out1, "glm.c.s" = out2, "BY" = out3, "BY.c.s" = out4, "WBY" = out5, "WBY.c.s" = out6))
}

betas <- get_betas(givemecredit)
betasLOG <- get_betas(LOGgivemecredit)


###################################################################################################################################

data.matrix <- cbind(rep(1, length(get_train(LOGgivemecredit, "SeriousDlqin2yrs"))), as.matrix(LOGgivemecredit[get_train(LOGgivemecredit, "SeriousDlqin2yrs"), -1]))

scores.original <- data.frame(Sglm = as.numeric(data.matrix %*% betas[[1]]), 
                              Sglm.c.s. = as.numeric(data.matrix %*% betas[[2]]), 
                              SBY = as.numeric(data.matrix %*% betas[[3]]), 
                              SBY.c.s = as.numeric(data.matrix %*% betas[[4]]), 
                              SWBY = as.numeric(data.matrix %*% betas[[5]]), 
                              SWBY.c.s. = as.numeric(data.matrix %*% betas[[6]]))

outpcd <- OutlierPCDist(x = data.matrix[, -1], grouping = as.factor(givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), "SeriousDlqin2yrs"]))
wrd <- as.logical(outpcd@flag)

scores.LOG <- data.frame(Sglm = as.numeric(data.matrix %*% betasLOG[[1]]), 
                         Sglm.c.s. = as.numeric(data.matrix %*% betasLOG[[2]]), 
                         SBY = as.numeric(data.matrix %*% betasLOG[[3]]), 
                         SBY.c.s = as.numeric(data.matrix %*% betasLOG[[4]]), 
                         SWBY = as.numeric(data.matrix %*% betasLOG[[5]]), 
                         SWBY.c.s. = as.numeric(data.matrix %*% betasLOG[[6]]))

outpcd <- OutlierPCDist(x = data.matrix[, -1], grouping = as.factor(givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), "SeriousDlqin2yrs"]))
wrdLOG <- as.logical(outpcd@flag)

scores.plot <- data.frame(variable = factor(c(rep(c("glm", "glm c.s", "BY", "BY c.s"), each=nrow(scores.original)), rep(c("WBY", "WBY c.s."), each = sum(wrd))), 
                                            levels = c("glm", "glm c.s", "BY", "BY c.s", "WBY", "WBY c.s.")), 
                          value = c(scores.original[, 1], scores.original[, 2], scores.original[, 3], scores.original[, 4], scores.original[wrd, 5], scores.original[wrd, 6]), 
                          class = as.factor(c(rep(givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), "SeriousDlqin2yrs"], 4), rep(givemecredit$SeriousDlqin2yrs[get_train(givemecredit, "SeriousDlqin2yrs")][wrd], 2))))

scoresLOG.plot <- data.frame(variable = factor(c(rep(c("glm", "glm c.s", "BY", "BY c.s"), each=nrow(scores.LOG)), rep(c("WBY", "WBY c.s."), each = sum(wrd))), 
                                               levels = c("glm", "glm c.s", "BY", "BY c.s", "WBY", "WBY c.s.")), 
                             value = c(scores.LOG[, 1], scores.LOG[, 2], scores.LOG[, 3], scores.LOG[, 4], scores.LOG[wrd, 5], scores.LOG[wrd, 6]), 
                             class = as.factor(c(rep(LOGgivemecredit[get_train(LOGgivemecredit, "SeriousDlqin2yrs"), "SeriousDlqin2yrs"], 4), rep(LOGgivemecredit$SeriousDlqin2yrs[get_train(LOGgivemecredit, "SeriousDlqin2yrs")][wrd], 2))))


## both classes
ggplot(data = scores.plot, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=class)) + 
  ylim(c(-10, 10))

ggplot(data = scoresLOG.plot, aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=class)) + 
  ylim(c(-10, 10))
