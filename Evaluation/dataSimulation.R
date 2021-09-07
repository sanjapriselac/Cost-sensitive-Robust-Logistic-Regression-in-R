################################################################################################################################
## Code for the data simulation with various settings - imbalance proportaions, leverage detection, 
## dimension of the space of explanatory veriables 
################################################################################################################################

source("glmrobBy_UPDATE.R")
source("additionalFunctions.R")

library(pROC)
library(xtable)
library(robustbase)
library(glmnet)
library(tidyverse)
library(dplyr)
library(qwraps2)   ## masks auc

#################################################################################################################################

simulate_results <- function(p, gamma, outmethod = "mcd") {
  n = 5000
  set.seed(24)
  numrep = 100
  
  ## define tables
  models <- c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY", "WBY cost-sensitive") 
  caseI <- data.frame(matrix(ncol=24, nrow = numrep))
  colnames(caseI) <- c(paste(models, 20), paste(models, 10), paste(models, 5), paste(models, 1))

  caseII <- caseI
  caseIII <- caseI
  caseMiss <- caseI

  
  for (i in 1:numrep) {
    cat("Iteration", i, "\n")
    x <- replicate(p, rnorm(n))
    eps <- rlogis(n, scale = 2)
  
    ## find c - binary search 
    c <- rep(NA, 4)
    yraw <- x%*%gamma + eps
    c[1] <- find_c(yraw, 0.2)
    cat("done 20% \n")
    c[2] <- find_c(yraw, 0.1)
    cat("done 10% \n")
    c[3] <- find_c(yraw, 0.05)
    cat("done 5% \n")
    c[4] <- find_c(yraw, 0.01)
    cat("done 1% \n")
    
    x1 <- replicate((p-1), runif(0.1*n,  min=(-1), max = 4))
  
    for (j in 1:length(c)) {
      ### !!!! put c/2 for beta
      caseIIpoints <- cbind(x1, apply(x1, 1, funy, const = c[j], case = 2.5, p = p) )
      caseIIIpoints <- cbind(x1, apply(x1, 1, funy, const = c[j], case = 5, p = p) )
      y <- as.numeric(yraw>c[j])
      cat("c = ", c[j], "\n")
      cat("Outliers: ", sum(y)/length(y), "\n")
      
      data <- data.frame(x,Y =y)
      dataII <- data.frame(rbind(cbind(x, y), cbind(caseIIpoints, y =rep(0, nrow(caseIIpoints)))))
      dataIII <- data.frame(rbind(cbind(x, y), cbind(caseIIIpoints, y =rep(0, nrow(caseIIIpoints)))))
      names(data) <- c(paste0("X", 1:p), "Y")
      names(dataII) <- c(paste0("X", 1:p), "Y")
      names(dataIII) <- c(paste0("X", 1:p), "Y")
      
      missclasified <- sample(which(y == 0), size = round(length(y)/10))
      y[missclasified] <- 1
      dataMiss <- data.frame(x,Y =y)
      names(dataMiss) <- c(paste0("X", 1:p), "Y")
      
      caseI[i, ((j-1)*6+1):((j-1)*6+6)] <- calculate_estimators(data, outmethod)
      caseII[i, ((j-1)*6+1):((j-1)*6+6)] <- calculate_estimators(dataII, outmethod)
      caseIII[i, ((j-1)*6+1):((j-1)*6+6)] <- calculate_estimators(dataIII, outmethod)
      caseMiss[i, ((j-1)*6+1):((j-1)*6+6)] <- calculate_estimators(dataMiss, outmethod)
    }
  }

  write.csv(caseI, paste0("SimulationResults/caseIp", p, outmethod, ".csv"), row.names = FALSE)  
  write.csv(caseII, paste0("SimulationResults/caseIIp", p, outmethod, ".csv"), row.names = FALSE) 
  write.csv(caseIII, paste0("SimulationResults/caseIIIp", p, outmethod,".csv"), row.names = FALSE) 
  write.csv(caseMiss, paste0("SimulationResults/casep", p, outmethod, "miss.csv"), row.names = FALSE)  
}





################################################################################################################################
## Simulation
################################################################################################################################

simulate_results(2, c(2,2),outmethod = "mcd")
simulate_results(2, c(2,2),outmethod = "pcout")

## change the constant!! 
simulate_results(10, rep(1, 10),outmethod = "mcd")
simulate_results(10, rep(1, 10),outmethod = "pcout")

