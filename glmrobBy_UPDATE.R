#### http://www.econ.kuleuven.be/public/NDBAE06/programs/roblog/ :

##  Computation of the estimator of Bianco and Yohai (1996) in logistic regression
##  -------------
##  Christophe Croux, Gentiane Haesbroeck
##  (thanks to Kristel Joossens and Valentin Todorov for improving the code) -
##  ==> Now "contains" both the *weighted* and regular, unweighted  BY-estimator
##
##  This program computes the estimator of Bianco and Yohai in
##  logistic regression. By default, an intercept term is included
##  and p parameters are estimated.

## The adptation of the existing functions for the imbalance learning

##################################################################################################
## glmrob
##################################################################################################

glmrobSP <- function (formula, family, data, weights, subset, na.action, 
            start = NULL, offset, method = c("Mqle", "BY", "WBY", "MT"), outmethod = c("mcd", "pcout"), 
            weights.on.x = c("none", "hat", "robCov", "covMcd"), control = NULL, 
            model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, trace.lev = 0, 
            ...) 
  {
    call <- match.call()
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    fami <- family$family
    if (is.null(fami)) 
      stop(gettextf("'%s' is not a valid family (see ?family)", 
                    as.character(call[["family"]])), domain = NA)
    if (!(fami %in% c("binomial", "poisson", "Gamma", "gaussian"))) {
      stop(gettextf("Robust GLM fitting not yet implemented for family %s", 
                    fami), domain = NA)
    }
    if (missing(data)) 
      data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) 
      return(mf)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm)) 
        names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
      model.matrix(mt, mf, contrasts)
    else matrix(NA_real_, NROW(Y), 0)
    weights <- model.weights(mf)
    offset <- model.offset(mf)
    if (!is.null(weights) && any(weights < 0)) 
      stop("'weights' must be non-negative")
    if (!is.null(offset) && length(offset) != NROW(Y)) 
      stop(gettextf("Number of offsets is %d, should rather equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
    method <- match.arg(method)
    meth. <- if (method == "WBY")  #SP commented
      "BY"
    else method
    if (is.null(control)) 
      control <- get(paste0("glmrob", meth., "SP.control"))(...)  ## SP other contol function added  
    if (missing(weights.on.x) || is.character(weights.on.x)) 
      weights.on.x <- match.arg(weights.on.x)
    else if (!(is.function(weights.on.x) || is.list(weights.on.x) || 
               (is.numeric(weights.on.x) && length(weights.on.x) == 
                NROW(Y)))) 
      stop("'weights.on.x' must be a string, function, list or numeric n-vector")
    if (!is.null(start) && !is.numeric(start)) {
      if (!is.character(start)) 
        stop("'start' must be a numeric vector, NULL, or a character string")
      start <- switch(start, lmrob = , lmrobMM = {
        if (!is.null(weights)) warnings("weights are not yet used in computing start estimate")
        lmrob.fit(x = X, y = family$linkinv(Y), control = lmrob.control())$coefficients
      }, stop("invalid 'start' string"))
    }
    ## SP fit is only for BY/WBY
    fit <- if (fami != "binomial") {
              stop(gettextf("method='%s' is only applicable for binomial family, but family=\"\"", method, fami), domain = NA)
    } else {  ##SP glmrobBYSP
      glmrobBYSP(X = X, y = Y, weights = weights, start = start, 
               method = method, outmethod = outmethod, weights.on.x = weights.on.x, 
               control = control, intercept = attr(mt, "intercept") > 
                 0, trace.lev = trace.lev)
    }
    
    fit$na.action <- attr(mf, "na.action")
    if (model) 
      fit$model <- mf
    if (x) 
      fit$x <- X
    if (!y) 
      warning("setting 'y = FALSE' has no longer any effect")
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
                       data = data, offset = offset, control = control, method = method, 
                       prior.weights = if (is.null(weights)) rep.int(1, nrow(X)) else weights, 
                       contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
                                                                               mf)))
    class(fit) <- c("glmrob", "glm")
    fit
  }


##################################################################################################
## BYlogregSP
##################################################################################################

BYlogregSP <- function(x0, y, method = "BY", initwml=TRUE, weights=rep(1, nrow(x0)), outmethod = outmethod,   ## SP chnage default weights to rep.int(1, nobs); outmethod
                     addIntercept=TRUE,
                     const=0.5, kmax = 1000, maxhalf = 10, 
                     sigma.min = 1e-4, trace.lev=0)
{
  if(!is.numeric(y))
    y <- as.numeric(y)
  #if(!is.null(w.x))
  #     warning("x weights  'w.x'  are not yet made use of")
  if(!is.null(dim(y))) {
    if(ncol(y) != 1)
      stop("y is not onedimensional")
    y <- as.vector(y)
  }
  n <- length(y)
  
  if(is.data.frame(x0)) {
    x0 <- data.matrix(x0)
  } else if (!is.matrix(x0)) {
    x0 <- matrix(x0, length(x0), 1,
                 dimnames = list(names(x0), deparse(substitute(x0))))
  }
  if(nrow(x0) != n)
    stop("Number of observations in x and y not equal")
  
  na.x <- !is.finite(rowSums(x0))
  na.y <- !is.finite(y)
  ok <- !(na.x | na.y)
  if(!all(ok)) {
    x0 <- x0[ok, , drop = FALSE]
    y  <- y [ok] # y[ok, , drop = FALSE]
  }
  
  if(addIntercept) {
    x <- cbind("Intercept" = 1, x0)
  } else { # x0 := x without the  "intercept column"
    x <- x0
    all1 <- apply(x == 1, 2, all)
    if(any(all1))
      x0 <- x[,!all1, drop = FALSE]
    else message("no intercept in the model")
  }
  
  dx <- dim(x)
  n <- dx[1]   
  if(n == 0)
    stop("All observations have missing values!")
  p <- dx[2] # == ncol(x)
  
  family <- binomial()
  ## Computation of the initial value of the optimization process
  gstart <-
    if(initwml) {  
      if (outmethod == "mcd") {
        mcd <- covMcd(x0, alpha=0.75, tolSolve = 1e-20) 
        D <- mahalanobis(mcd$X, mcd$center, mcd$cov) 
        vc  <- qchisq(0.975, p-1) 
        wrd <- D <= vc 
        
        ## SP added weighted Bianco Yohai
        if (method == "WBY") {
          wby <- as.numeric(wrd)
        }
      } else {
        outpcd <- OutlierPCDist(x0, grouping = as.factor(y))
        wrd <- as.logical(outpcd@flag)
        #cat("weights: ", wrd, "\n", "sum(y)", sum(y[wrd]))
        
        ## SP added weighted Bianco Yohai
        if (method == "WBY") {
          wby <- outpcd@flag
        }
      }
      
      
      glm.fit(x[wrd,], y[wrd], weights = weights[wrd], family=family)$coef   ## SP "prior" weights added in the glm.fit formula
    } else {
      
      ## SP added weighted Bianco Yohai
      if (method == "WBY") { 
        if (outmethod == "mcd") {
          mcd <- covMcd(x0, alpha=0.75, tolSolve = 1e-40) 
          D <-  mahalanobis(mcd$X, mcd$center, mcd$cov) 
          vc  <- qchisq(0.975, p-1)
          wby <- as.numeric(D <= vc)
        } else {
          outpcd <- OutlierPCDist(x0, grouping = as.factor(y))
          wby <- outpcd@flag
        }
        
      }
      glm.fit(x, y, weights = weights, family=family)$coef   ## SP "prior" weights added in the glm.fit formula
    }
  
  if (method == "WBY") {
    weights <- weights * wby
  }
  
  sigma1 <- 1/sqrt(sum(gstart^2))
  xistart <- gstart*sigma1
  stscores <- x %*% xistart
  
  ## Initial value for the objective function
  oldobj <- mean(phiBY3(stscores/sigma1, y, const) * weights)  ## SP weights added in the objective (str checked)
  
  converged <- FALSE
  kstep <- 1L
  while(kstep < kmax && !converged)
  {
    unisig <- function(sigma) mean(phiBY3(stscores/sigma, y, const) * weights)  ## SP objective function as a fnction of sigma, weights added
    ##                             ------
    #optimsig <- nlminb(sigma1, unisig, lower=0)# "FIXME" arguments to nlminb()    ## SP nlminb optimisation algorithm!  
    ##                  ======
    optimsig <- optimize(unisig, interval = c(0, 10^10))   ## SP chnaged optimize optimsig$minimum
    if(trace.lev) cat(sprintf("k=%2d, s1=%12.8g: => new s1= %12.8g",
                              kstep, sigma1, optimsig$minimum)) # MM: jhalf =!?= 1 here ?? ## SP change optimsig$par because of optimize method
    sigma1 <- optimsig$minimum ## SP changed to minimum because of the optimize function
    
    if(sigma1 < sigma.min) {
      if(trace.lev) cat("\n")
      warning(gettextf("Implosion: sigma1=%g became too small", sigma1))
      kstep <- kmax #-> *no* convergence
    } else {
      ## gamma1 <- xistart/sigma1
      scores <- stscores/sigma1
      newobj <- mean(phiBY3(scores, y,const) * weights)   ## SP weights added
      oldobj <- newobj
      grad.BY <- colMeans(((derphiBY3(scores,y,const)*weights) %*% matrix(1,ncol=p))*x) ## SP weights added
      h <- -grad.BY + as.numeric(grad.BY %*% xistart) *xistart ## SP as.numeric added
      finalstep <- h/sqrt(sum(h^2))   ## SP h/||h||
      
      if(trace.lev) {
        if(trace.lev >= 2) cat(sprintf(", obj=%12.9g: ", oldobj))
        cat("\n")
      }
      
      ## FIXME repeat { ... }   {{next 4 lines are also inside while(..) below}}
      xi1 <- xistart+finalstep
      xi1 <- xi1/sum(xi1^2)
      scores1 <- (x %*% xi1)/sigma1
      newobj <- mean(phiBY3(scores1,y,const) * weights) ## SP weights added 
      
      ## If 'newobj' is not better, try taking a smaller step size:
      hstep <- 1.
      jhalf <- 1L
      ## SP maxhalf changed from 10 to 100
      while(jhalf <= maxhalf & newobj > oldobj)
      {
        hstep <- hstep/2
        xi1 <- xistart+finalstep*hstep
        xi1 <- xi1/sqrt(sum(xi1^2))
        scores1 <- x %*% xi1/sigma1
        newobj <- mean(phiBY3(scores1,y,const) * weights)   ## SP weights added
        if(trace.lev >= 2)
          cat(sprintf("  jh=%2d, hstep=%13.8g => new obj=%13.9g\n",
                      jhalf, hstep, newobj))
        jhalf <- jhalf+1L
      }
      
      converged <-
        not.improved <- (jhalf > maxhalf && newobj > oldobj)
      if(not.improved) {
        ## newobj is "worse" and step halving did not improve
        message("Convergence Achieved")
      } else {
        jhalf <- 1L
        xistart <- xi1
        oldobj <- newobj
        stscores <- x %*% xi1
        kstep <- kstep+1L
      }
    }
  } ## while( kstep )
  
  if(kstep == kmax) {
    warning("No convergence in ", kstep, " steps.")
    #list(convergence=FALSE, objective=0, coefficients= rep(NA,p)) ## SP delete
  } #else {
  gammaest <- xistart/sigma1  # SP the estimator
  V <- vcovBY3(x, y, const, estim=gammaest, addIntercept=FALSE)  ## SP variance covariance matrix, does it change???
  list(convergence=TRUE, objective=oldobj, coefficients=gammaest,
        cov = V, sterror = sqrt(diag(V)),
        iter = kstep)
  #}
}


##################################################################################################
## glmrobBYSP.control
##################################################################################################

glmrobBYSP.control <-
  function(maxit = 1000, const = 0.5, maxhalf = 10)  {
    if(!is.numeric(maxit) || maxit <= 0)
      stop("maximum number of \"kstep\" iterations must be > 0")
    if(!is.numeric(maxhalf) || maxhalf <= 0)
      stop("maximal number of *inner* step halvings must be > 0")
    ## if (!is.numeric(tcc) || tcc <= 0)
    ##     stop("value of the tuning constant c (tcc) must be > 0")
    if(!is.numeric(const) || const <= 0)
      stop("value of the tuning constant c ('const') must be > 0")
    
    list(## acc = acc, consttest.acc = test.acc,
      const=const,
      maxhalf=maxhalf,
      maxit=maxit #, tcc = tcc
    )
  }


##################################################################################################
## glmrobBYSP
##################################################################################################

##' @param intercept logical, if true, X[,] has an intercept column which should
##'                  not be used for rob.wts
glmrobBYSP <- function(X, y,
                     weights = NULL, start = NULL, offset = NULL,
                     method = c("WBY","BY"), weights.on.x = "none",
                     outmethod = outmethod,                           ## SP added
                     control = glmrobBYSP.control(...), intercept = TRUE,
                     trace.lev = 0, ...)
{
  ### THIS is *NOT* exported
  
  method <- match.arg(method)
  ## SP comment error prior weights not implemented
  #if(!is.null(weights) || any(weights != 1)) ## FIXME (?)
  #  stop("non-trivial prior 'weights' are not yet implemented for \"BY\"")
  if(!is.null(start))
    stop(" 'start' cannot yet be passed to glmrobBYSP()")
  if(!is.null(offset))
    stop(" 'offset' is not yet implemented for \"BY\"")
  const   <- if(is.null(cc <- control$const  )) 0.5 else cc
  kmax    <- if(is.null(cc <- control$maxit  )) 1e3 else cc 
  maxhalf <- if(is.null(cc <- control$maxhalf))  10 else cc
  if(!identical(weights.on.x, "none"))
    stop("'weights.on.x' = ", format(weights.on.x)," is not implemented")

    r <- BYlogregSP(x0=X, y=y, initwml = TRUE, method = method, outmethod = outmethod, weights = weights,
                addIntercept = !intercept, ## add intercept if there is none
                const=const, kmax=kmax, maxhalf=maxhalf,
                ## FIXME sigma.min  (is currently x-scale dependent !????)
                trace.lev=trace.lev)   # SP changed initwml = (method == "WBY")
  ## FIXME: make result more "compatible" with other glmrob() methods
  r
}


##################################################################################################
### Functions needed for the computation of estimator of Bianco and Yohai ----------------------
##################################################################################################

## From their paper:

## A last remark is worth mentioning: when huge outliers occur in
## the logistic regression setting, often numerical imprecision occurs in the computation
## of the deviances given by
##    d(s;y_i)= -y_i log F(s) - (1-y_i) log{1-F(s)} .
##
## Instead of directly computing this expression, it can be seen that a
## numerically more stable and accurate formula is given by
##    log(1 + exp(-abs(s))) + abs(s)* ((y-0.5)*s < 0)
## in which the second term equals abs(s) if the observation is misclassified, 0 otherwise.
dev1 <- function(s,y) log(1+exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev2 <- function(s,y) log1p(exp(-abs(s))) + abs(s)*((y-0.5)*s<0)
dev3 <- function(s,y) -( y  * plogis(s, log.p=TRUE) +
                           (1-y)*plogis(s, lower.tail=FALSE, log.p=TRUE))
## for now,  100% back-compatibility:
devBY <- dev1
rm(dev1, dev2, dev3)


phiBY3 <- function(s,y,c3)
{
  s <- as.double(s)
  dev. <- devBY(s,y)
  ## FIXME: GBY3Fs()  computes the 'dev' above *again*, and
  ##        GBY3Fsm() does with 's>0' instead of 's<0'
  sc3 <- sqrt(c3)
  ## SP G1 for the exact value of rho function
  #G1 <- exp(0.25)*sqrt(pi)*(pnorm(sqrt(2)*(0.5+sc3))-1) + exp(-sc3)
  rhoBY3(dev.,c3) + GBY3Fs(s,c3) + GBY3Fsm(s,c3) #- G1
}

rhoBY3 <- function(t,c3)   
  {
  ec3 <- exp(-sqrt(c3))
  t*ec3* (t <= c3) +
    (ec3*(2+(2*sqrt(c3))+c3) - 2*exp(-sqrt(t))*(1+sqrt(t)))* (t > c3)
}


psiBY3 <- function(t,c3)
{
  exp(-sqrt(c3)) *(t <= c3) +
    exp(-sqrt( t)) *(t > c3)
}

derpsiBY3 <- function(t, c3)
{
  r <- t
  r[in. <- (t <= c3)] <- 0
  if(any(out <- !in.)) {
    t <- t[out]
    st <- sqrt(t)
    r[out] <- -exp(-st)/(2*st)
  }
  r
}

sigmaBY3 <- function(sigma,s,y,c3)
{
  mean(phiBY3(s/sigma,y,c3))
}

derphiBY3 <- function(s,y,c3)
{
  Fs <- exp(-devBY(s,1))  ## SP pi(1) = exp(s)/(1+exp(s))
  ds <- Fs*(1-Fs) ## MM FIXME: use expm1()
  dev. <- devBY(s,y)
  Gprim1 <- devBY(s,1)
  Gprim2 <- devBY(-s,1)
  
  -psiBY3(dev.,c3)*(y-Fs) + ds*(psiBY3(Gprim1,c3) - psiBY3(Gprim2,c3))
}

## SP TBC
der2phiBY3 <- function(s, y, c3)
{
  s <- as.double(s)
  Fs <- exp(-devBY(s,1))
  ds <- Fs*(1-Fs) ## MM FIXME: use expm1()
  dev. <- devBY(s,y)
  Gprim1 <- devBY(s,1)
  Gprim2 <- devBY(-s,1)
  der2 <- derpsiBY3(dev.,c3)*(Fs-y)^2  + ds*psiBY3(dev.,c3)
  der2 <- der2+ ds*(1-2*Fs)*(psiBY3(Gprim1,c3) - psiBY3(Gprim2,c3))
  der2 - ds*(derpsiBY3(Gprim1,c3)*(1-Fs) +
               derpsiBY3(Gprim2,c3)*  Fs )
}


GBY3Fs <- function(s,c3)
{
  e.f <- exp(0.25)*sqrt(pi)  
  ## MM FIXME: Fs = exp(..) and below use  log(Fs) !!
  Fs <- exp(-devBY(s,1))  ## SP pi(s) 
  resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fs))))-1) 
  ## MM FIXME: use expm1():
  ## SP indicators come from F(s) = exp(s)/(1+exp(s)) and basic log calculation 
  resGinf <- (resGinf+(Fs*exp(-sqrt(-log(Fs)))))*as.numeric(s <= -log(exp(c3)-1))
  resGsup <- ((Fs*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
    as.numeric(s > -log(exp(c3)-1))
  resGinf + resGsup
}

GBY3Fsm <- function(s,c3)
{
  e.f <- exp(0.25)*sqrt(pi)
  ## MM FIXME: Fsm = exp(..) and below use  log(Fsm) !!
  Fsm <- exp(-devBY(-s,1)) ## SP pi(-s) = pi(0) = 1/(1+exp(s))
  resGinf <- e.f*(pnorm(sqrt(2)*(0.5+sqrt(-log(Fsm))))-1)
  ## MM FIXME: use expm1():
  resGinf <- (resGinf+(Fsm*exp(-sqrt(-log(Fsm))))) * as.numeric(s >= log(exp(c3)-1))
  resGsup <- ((Fsm*exp(-sqrt(c3)))+(e.f*(pnorm(sqrt(2)*(0.5+sqrt(c3)))-1))) *
    as.numeric(s < log(exp(c3)-1))
  resGinf + resGsup
}

## Compute the standard erros of the estimates -
##  this is done by estimating the asymptotic variance of the normal
##  limiting distribution of the BY estimator - as derived in Bianco
##  and Yohai (1996)
##
sterby3 <- function(x0, y, const, estim, addIntercept)
{
  sqrt(diag(vcovBY3(x0, y, const=const, estim=estim, addIntercept=addIntercept)))
}

vcovBY3 <- function(z, y, const, estim, addIntercept)
{
  stopifnot(length(dim(z)) == 2)
  if(addIntercept) z <- cbind(1, z)
  d <- dim(z)
  n <- d[1]
  p <- d[2]
  argum <- z %*% estim
  matM <- IFsqr <- matrix(0, p, p)
  for(i in 1:n)
  {
    myscalar <- as.numeric(der2phiBY3(argum[i],y[i], c3=const))
    zzt <- tcrossprod(z[i,])
    matM <- matM + myscalar * zzt
    IFsqr <- IFsqr + derphiBY3(argum[i],y[i], c3=const)^2 * zzt
  }
  
  matM    <- matM/n
  matMinv <- solve(matM)
  IFsqr <- IFsqr/n
  ## Now,  asymp.cov  =  matMinv %*% IFsqr %*% t(matMinv)
  
  ## provide  vcov(): the full matrix
  (matMinv %*% IFsqr %*% t(matMinv))/n
}