
##########################################################################################################
## Additional functions - simulation 
##########################################################################################################

find_c <- function(ybar, proc, const = seq(1, 17, by = 0.000001), tol= 0.0005) {
  index <- round(length(const)/2)
  y <- as.numeric(ybar>const[index])
  p = sum(y)/length(y)
  #cat("index:", index, "c = ", const[index], "perc = ", p, "\n")
  
  if(p < (proc+tol) && p > (proc-tol)) {
    return(const[index])
  }  else {
    if (p > (proc+tol)) {
      find_c(ybar, proc, const = const[(index + 1):length(const)])
    } else {
      find_c(ybar, proc, const = const[1:(index-1)])
    }
  }
}

## put const = const/2 for p=2; case = 2.5/5
funy <- function(x, const, case, p) {
  return(const + case*sqrt(p) - sum(x))
}



##########################################################################################################
## Additional functions - general
##########################################################################################################

gini <- function(pred.prob, class, data) {
  ## train: indicates if the index vector is train or test
  
  mydata <- data
  mydata$pred <- pred.prob 
  roc <- roc(mydata[, class]~pred, data=mydata, quiet = TRUE)
  gini <- 2*pROC::auc(roc) - 1
  #plot(roc)
  return(gini)
}

calculate_estimators <- function (data) {
  betas <- matrix(NA, nrow = ncol(data), ncol = 8)
  
  model <- glm(Y ~ ., data = data, family=binomial)
  betas[, 1] <- model$coefficients
  
  model <- glm(Y ~ ., data = data, weights = get_weights(data, "Y"), family=binomial)
  betas[, 2] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "BY", trace.lev =FALSE,
                    weights = rep(1, nrow(data)), outmethod = "mcd",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 3] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "BY", trace.lev = FALSE,
                    weights = get_weights(data, "Y"), outmethod = "mcd",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 4] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "WBY", trace.lev = FALSE,
                    weights = rep(1, nrow(data)), outmethod = "mcd",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 5] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "WBY", trace.lev = FALSE,
                    weights = get_weights(data, "Y"), outmethod = "mcd",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 6] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "WBY", trace.lev = FALSE,
                    weights = rep(1, nrow(data)), outmethod = "pcdist",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 7] <- model$coefficients
  
  model <- glmrobSP(Y ~ ., data = data,
                    family=binomial, method = "WBY", trace.lev = FALSE,
                    weights = get_weights(data, "Y"), outmethod = "pcdist",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  betas[, 8] <- model$coefficients
  
  return(betas)
}


calculate_estimatorsGMC <- function (data) {
  out <- rep(NA, 6)
  
  model <- glm(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ], family=binomial)
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"), ], type = "response")
  out[1] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  
  
  model <- glm(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ], weights = get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")] , family=binomial)
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"), ], type = "response")
  out[2] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(data, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[3] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights =  get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[4] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(data, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[5] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = data[get_train(data, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights =  get_weights(data, "SeriousDlqin2yrs")[get_train(data, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[6] <- gini(pred, index = get_train(data, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = data)
  return(out)
}

###################
### train
###################

## 75% in train
get_train <- function(data, class){
  set.seed(100)
  index_grp0 <- which(data[, class]==0)
  index_grp1 <- which(data[, class]==1)
  train_grp0 <- sample(index_grp0,size=round(3/4*length(index_grp0)))
  train_grp1 <- sample(index_grp1,size=round(3/4*length(index_grp1)))
  train_index <- c(train_grp0,train_grp1)
  return(train_index)
}

###################
### weights
###################
get_weights <- function(data, class) {
  tb <- table(data[, class])
  wght <- rep(NA,nrow(data))
  wght[data[, class] == 0] <- tb[2]/sum(tb)
  wght[data[, class] == 1] <- tb[1]/sum(tb)
  return(wght)
}

