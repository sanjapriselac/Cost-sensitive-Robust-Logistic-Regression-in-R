
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

gini <- function(pred.prob, index, class, data, train=TRUE) {
  ## train: indicates if the index vector is train or test
  
  if (train) {
    mydata <- data[-index,]
  } else {
    mydata <- data[index,]
  }
  
  mydata$pred <- pred.prob 
  roc <- roc(mydata[, class]~pred, data=mydata)
  gini <- 2*pROC::auc(roc) - 1
  plot(roc)
  return(gini)
}

calculate_estimators <- function (data, outmethod) {
  out <- rep(NA, 6)
  
  model <- glm(Y ~ ., data = data[get_train(data, "Y"), ], family=binomial)
  pred <- predict(model, data[-get_train(data, "Y"), ], type = "response")
  out[1] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)
  
  model <- glm(Y ~ ., data = data[get_train(data, "Y"), ], weights = get_weights(data, "Y")[get_train(data, "Y")], family=binomial)
  pred <- predict(model, data[-get_train(data, "Y"), ], type = "response")
  out[2] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)
  
  model <- glmrobSP(Y ~ ., data = data[get_train(data, "Y"), ],
                    family=binomial, method = "BY", trace.lev =FALSE, 
                    weights = rep(1, length(get_train(data, "Y"))), outmethod = outmethod, 
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "Y"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[3] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)
  
  model <- glmrobSP(Y ~ ., data = data[get_train(data, "Y"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights = get_weights(data, "Y")[get_train(data, "Y")], outmethod = outmethod, 
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "Y"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[4] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)

  model <- glmrobSP(Y ~ ., data = data[get_train(data, "Y"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(data, "Y"))), outmethod = outmethod, 
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "Y"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[5] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)

  model <- glmrobSP(Y ~ ., data = data[get_train(data, "Y"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights = get_weights(data, "Y")[get_train(data, "Y")], outmethod = outmethod, 
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, data[-get_train(data, "Y"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[6] <- gini(pred, index = get_train(data, "Y"), class= "Y", data = data)
  return(out)
}

calculate_estimatorsGMC <- function () {
  out <- rep(NA, 5)
  
  model <- glm(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ], family=binomial)
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"), ], type = "response")
  out[1] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
  
  
  model <- glm(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ], weights = get_weights(givemecredit, "SeriousDlqin2yrs")[get_train(givemecredit, "SeriousDlqin2yrs")] , family=binomial)
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"), ], type = "response")
  out[2] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
  
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(givemecredit, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[3] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "BY", trace.lev = FALSE, 
                    weights =  get_weights(givemecredit, "SeriousDlqin2yrs")[get_train(givemecredit, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[4] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights = rep(1, length(get_train(givemecredit, "SeriousDlqin2yrs"))), outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[5] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
  
  model <- glmrobSP(SeriousDlqin2yrs ~ ., data = givemecredit[get_train(givemecredit, "SeriousDlqin2yrs"), ],
                    family=binomial, method = "WBY", trace.lev = FALSE, 
                    weights =  get_weights(givemecredit, "SeriousDlqin2yrs")[get_train(givemecredit, "SeriousDlqin2yrs")], outmethod = "pcout",
                    control = list(maxit = 300, const = 0.5, maxhalf = 10))
  
  pred <- predict(model, givemecredit[-get_train(givemecredit, "SeriousDlqin2yrs"),], type = "link") #response not supported
  pred <- exp(pred)/(1+exp(pred))
  out[6] <- gini(pred, index = get_train(givemecredit, "SeriousDlqin2yrs"), class= "SeriousDlqin2yrs", data = givemecredit)
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

