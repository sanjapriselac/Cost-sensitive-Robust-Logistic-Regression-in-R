################################################################################################################################
## Code for the data simulation with various settings - imbalance proportaions, leverage detection, 
## dimension of the space of explanatory veriables 
################################################################################################################################

source("glmrobBy_UPDATE.R")
source("additionalFunctions.R")

library(pROC)
library(robustbase)
library(glmnet)

#################################################################################################################################

simulate_results <- function(p, gamma) {
  n = 5000
  set.seed(243)
  numrep = 500
  
  ##############################
  ## define tables
  models <- c("glm", "glm cost-sensitive", "BY", "BY cost-sensitive", "WBY (mcd)", "WBY cost-sensitive (mcd)", "WBY (pcdist)", "WBY cost-sensitive (pcdist)") 
  
  betasI <- data.frame(matrix(ncol=32, nrow = numrep*(p+1)))
  colnames(betasI) <- c(paste(models, 20), paste(models, 10), paste(models, 5), paste(models, 1))
  betasII <- betasI
  betasIII <- betasI
  betasIV <- betasI
  ###############################
  
  for (i in 1:numrep) {
    cat("Iteration", i, "\n")
    x <- replicate(p, rnorm(n))
    eps <- rlogis(n, scale = 1)
    
    ## find c - binary search 
    c <- rep(NA, 4)
    yraw <- x%*%gamma + eps
    skip_to_next <- FALSE
    tryCatch(c[1] <- find_c(yraw, 0.2),
             error = function(e) { skip_to_next <<- TRUE})
    tryCatch(c[2] <- find_c(yraw, 0.1),
             error = function(e) { skip_to_next <<- TRUE})
    tryCatch(c[3] <- find_c(yraw, 0.05),
             error = function(e) { skip_to_next <<- TRUE})
    tryCatch(c[4] <- find_c(yraw, 0.01),
             error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { next }    

    for (j in 1:length(c)) {
      y <- as.numeric(yraw>c[j])
      
      repl.obs <- sample(which(y == 0), size = round(length(y)/10))
      repl.obsIV <- repl.obs[1:round(length(repl.obs)/2)]
      caseIIpoints <- cbind(x[repl.obs, -p], apply(as.matrix(x[repl.obs, -p]), 1, funy, const = c[j]/gamma[1], case = 5, p = p))
      caseIVpoints <- cbind(x[repl.obsIV, -p], apply(as.matrix(x[repl.obsIV, -p]), 1, funy, const = c[j]/gamma[1], case = 5, p = p) )
      yIII <- y
      yIII[repl.obs] <- 1
      yIV <- y
      yIV[repl.obs[(round(length(repl.obs)/2)+1):length(repl.obs)]] <- 1
      
      dataI <- data.frame(x,Y =y)
      dataII <- data.frame(rbind(cbind(x[-repl.obs, ], y[-repl.obs]), cbind(caseIIpoints, y =rep(0, nrow(caseIIpoints)))))
      dataIII <- data.frame(x,Y =yIII)
      dataIV <- data.frame(rbind(cbind(x[-repl.obsIV, ], yIV[-repl.obsIV]), cbind(caseIVpoints, y =rep(0, nrow(caseIVpoints)))))
      
      names(dataI) <- c(paste0("X", 1:p), "Y")
      names(dataII) <- c(paste0("X", 1:p), "Y")
      names(dataIII) <- c(paste0("X", 1:p), "Y")
      names(dataIV) <- c(paste0("X", 1:p), "Y")
      
      betasI[((i-1)*(p+1)+1):((p+1)*i), ((j-1)*8+1):(8*j)] <- calculate_estimators(dataI)
      betasII[((i-1)*(p+1)+1):((p+1)*i), ((j-1)*8+1):(8*j)] <-  calculate_estimators(dataII)
      betasIII[((i-1)*(p+1)+1):((p+1)*i), ((j-1)*8+1):(8*j)] <-  calculate_estimators(dataIII)
      betasIV[((i-1)*(p+1)+1):((p+1)*i), ((j-1)*8+1):(8*j)] <- calculate_estimators(dataIV)
    }
  }
  write.csv(betasI, paste0("SimulationResults/betasIp", p, ".csv"), row.names = FALSE)
  write.csv(betasII, paste0("SimulationResults/betasIIp", p,".csv"), row.names = FALSE)
  write.csv(betasIII, paste0("SimulationResults/betasIIIp", p,".csv"), row.names = FALSE)
  write.csv(betasIV, paste0("SimulationResults/betasIVp", p, ".csv"), row.names = FALSE)
}





################################################################################################################################
## Simulation
################################################################################################################################

simulate_results(2, c(2,2))
simulate_results(2, c(2,2))



