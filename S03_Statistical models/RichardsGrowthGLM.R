# R code implementing the method described in
# 
# Alaimo Di Loro P., Divino F., Farcomeni A. Jona Lasinio G., 
#   Lovison G., Maruotti A., Mingione M. (2021) Nowcasting COVID-19 
#   incidence indicators during the Italian first outbreak. 
#   Statistics in Medicine, 40: 3843-3864. 
#   https://doi.org/10.1002/sim.9004
#
# Code kindly provided by the authors.

# Packages -------------------------------------------------------------

load_install_package <- function(pkg)
{
# Load a package and install if not already installed
  if(!suppressWarnings(require(pkg, character.only = TRUE)) )
  { 
    install.packages(pkg)
    require(pkg, character.only = TRUE)
  }
}

load_install_package("MASS")
load_install_package("BB")
load_install_package("optimr")
load_install_package("GA")
load_install_package("nplr")
load_install_package("mvtnorm")
load_install_package("corpcor")
load_install_package("Matrix")


# Growth GLM --------------------------------------------------------------

growthGLM <- function(di, ti,  X, family = "Poisson",
                      tPred = NA, alpha = 0.05,
                      maxiter = 1e4, runs = 500, nBoot = 1000)
{
  
  di <- as.vector(di)
  ti <- as.vector(ti)
  tDiff <- diff(ti)[1]
  # Covariates and parameters
  X <- as.matrix(X)
  ncovs <- ncol(X)
  # prediction horizon
  if (is.na(tPred))
  {
    tPred <- max(ti)
  }
  
  # Initial values for optim
  smooth <- loess(cumsum(di)~ti, degree = 2, span=0.4)
  
  np <- tryCatch(suppressWarnings(suppressMessages(
                 as.vector(unlist(nplr(smooth$x, smooth$fitted,
                                       useLog=FALSE, npars=5)@pars))
                 )), 
                 error=function(e) {
                   a <- as.vector(unlist(nplr(ti, convertToProp(di), 
                                              useLog=FALSE,npars=5)@pars))
                   return(a)
                 })
  
  inits <- c(log(np[4]), np[3], np[5], rep(0, ncovs), log(np[2]))
  lb <- c(-20, -20, 0, rep(-Inf, ncovs+1))
  ub <- c(20, 20, sum(di), rep(Inf, ncovs+1))
  
  # Picking the selected distribution
  if(family=="Poisson") 
  {
    npars <- length(inits)
    nrich <- npars-ncovs-1
    
    lik <- function(x) likPois(x, ti = ti, di=di, X=X)
    neglik <- function(x) -lik(x)
    
    neggr <- function(x) -PoisRichGradient(x, ti=ti, di=di, X=X)
    neghessfun <- function(x) -PoisRichHessian(x, ti=ti, di=di, X=X)
  }
  if(family=="Negative Binomial") 
  {
    inits <- c(inits, log(mean(di)))
    lb <- c(lb, 0)
    ub <- c(ub, 15)
    npars <- length(inits)
    nrich <- npars-ncovs-2
    
    lik <- function(x) likNB(x, ti = ti, di=di, X=X)
    neglik <- function(x) -lik(x)
    
    neggr <- function(x) -NBRichGradient(x, ti=ti, di=di, X=X)
    neghessfun <- function(x) -NBRichHessian(x, ti=ti, di=di, X=X)
  }
  if(family!="Negative Binomial" & family!="Poisson") 
  {
    stop("Family must be Negative Binomial or Poisson")
  }
  
  # First optimization
  inits2 <- optim(inits, neglik, gr=neggr, method="BFGS")$par
  inits3 <- suppressWarnings(
              BB::multiStart(rbind(inits, inits2), 
                             fn=neglik, gr=neggr, 
                             quiet = TRUE,
                             control=list(M=20, trace=FALSE))$par)

  # Second optimization using Genetic algorithm
  rg <- ga("real-valued", function(x) lik(x), 
           lower=lb, upper=ub,
           maxiter=maxiter/2, run=runs, optim=TRUE, 
           suggestions=rbind(inits, inits3), popSize=100, 
           monitor = FALSE)
  rg2 <- ga("real-valued", function(x) lik(x),
            lower=lb, upper=ub, suggestions=rbind(inits2, inits3),
            maxiter=maxiter/2, run=runs, optim=TRUE, popSize=100,
            monitor = FALSE)
  rgsols <- rbind(rg@solution, rg2@solution)
  rg3 <- ga("real-valued", function(x)  lik(x),
            lower=apply(rgsols, 2, min)-3,
            upper=apply(rgsols, 2, max)+3,
            maxiter=maxiter, run=runs, optim=TRUE, popSize=100,
            suggestions=rgsols,
            optimArgs=list(control=list(maxit=1000)),
            monitor = FALSE)
  
  # Third optimization using multistart
  notshown <- capture.output(
    propPars <- optimr::multistart(parmat = rbind(rgsols,
                                                  rg3@solution),
                                   fn = neglik, gr = neggr,
                                   control=list(trace=0)))
  
  # Picking the best solution
  best <- which.min(propPars$value)
  pars <- as.numeric(propPars[best, 1:npars])
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(nrich+ncovs)]
  logr <- pars[nrich+ncovs+1]
  
  # Point prediction
  newt <- seq(min(ti), tPred, by=tDiff)
  newX <- rbind(X, as.matrix(X[rep(nrow(X), length(newt)-length(ti)), ]))
  linPredDiff <- exp(newX%*%betapars)+exp(logr+ldiffRich(richpars, newt))
  linPredCum <- cumsum(linPredDiff)
  
  # Computing Hessian
  hess <- neghessfun(pars)
  if (any(is.nan(hess)))
  {
    hess <- jacobian(neggr, pars)
    if (any(is.nan(hess)))
    {
      hess <- hessian(neglik, pars)
    }
  }
  
  # Verifying positive definiteness
  PD <-  tryCatch(corpcor::is.positive.definite(hess),
                  error=function(e) {
                    return(corpcor::is.positive.definite(hess, 
                                                         method = "chol"))
                  })
  
  if(!PD)
  {
    hess <- as.matrix(Matrix::nearPD(hess)$mat)
    warning("Information matrix is not positive definite, 
            possible local optimum")
  }
  
  # Symmetrize
  invFish <- solve(hess)
  invFish <- (invFish+t(invFish))/2
  # Standard errors
  se <- sqrt(abs(diag(invFish)))
  
  # Prediction interval
  nPred <- length(newt)
  
  # Simulating parameters
  rpars <- rmvnorm(nBoot, pars, invFish)
  rpars <- rpars[rpars[,3]>=0, ]
  if (family=="Negative Binomial")
  {
    rpars[, npars] <- pars[npars]
  }
  
  # Simulating trajectories
  rlpmat <- apply(rpars, 1, function(x)
    {
      richpars <- x[1:nrich]
      betapars <- x[(nrich+1):(nrich+ncovs)]
      logr <- x[nrich+ncovs+1]
      exp(newX%*%betapars)+exp(logr + ldiffRich(pars=richpars, ti = newt))
  })
  
  # Storing objects for intervals
  # First differences
  diffcurves <- matrix(NA, nrow=nBoot, ncol=length(newt))
  
  # Simulating
  for(i in 1:length(newt))
  {
    if (family=="Poisson")
    {
      diffcurves[,i] <- rpois(nBoot, rlpmat[i, ])
    }
    if (family=="Negative Binomial")
    {
      diffcurves[,i] <- rnbinom(nBoot, 
                                size=exp(rpars[, npars]),
                                mu=rlpmat[i, ])
    }
  }
  
  # Cumulative
  curves <- apply(diffcurves, 1, cumsum)
  
  # Quantiles
  qtCurves <- apply(curves, 1, quantile, 
                    probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)
  qtU <- qtCurves[2,]
  qtL <- qtCurves[1,]
  diffqtCurves <- apply(diffcurves, 2, quantile, 
                        probs=c(alpha/2, 1-alpha/2), na.rm=TRUE)
  diffqtU <- diffqtCurves[2,]
  diffqtL <- diffqtCurves[1,]
  
  return(list(pars=pars, 
              optim = propPars[[best]],
              linPredDiff = linPredDiff[,1], 
              linPredCum = linPredCum,
              hessian = hess, 
              asyVar = invFish, 
              se = se,
              rpars=rpars, 
              lowdiff = diffqtL, updiff = diffqtU, 
              lowcum = qtL, upcum = qtU, 
              loglik = -propPars$value[best], 
              R2diff = cor(di, linPredDiff[1:length(ti)])^2,
              R2cum = cor(cumsum(di), linPredCum[1:length(ti)])^2,
              NoConv=(!PD)))
}


# Functions - Richards -----------------------------------------------------

# Linear Predictor
lRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # Denominator
  l1p10hptInf <- log(1+10^(h*(p-ti)))
  l1p10hpt <- ifelse(is.infinite(l1p10hptInf), h*(p-ti)*log(10), l1p10hptInf)
  logden <- -s*l1p10hpt
  
  # Out
  lout <- logden
  
  return(lout)
}

# Linear Predictor monotone
ldiffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # Auxiliary
  l1p10hptInf <- log(1+10^(h*(p-tAll)))
  l1p10hpt <- ifelse(is.infinite(l1p10hptInf), h*(p-tAll)*log(10), l1p10hptInf)
  
  # Terms
  # 1
  laux1 <- log(diff(exp( -s * l1p10hpt)))
  
  # Out  
  lout <- laux1
  
  return(lout)
}

# Derivate prime
d1diffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # Auxiliary
  l1p10hptInf <- log(1+10^(h*(p-tAll)))
  l1p10hpt <- ifelse(is.infinite(l1p10hptInf), h*(p-tAll)*log(10), l1p10hptInf)
  l10hpt <- h*(p-tAll)*log(10)
  
  # Terms
  #1
  laux1 <- log(diff(exp( -s * l1p10hpt)))
  
  # First derivatives
  
  # h
  dlh <- -s * exp(log(log(10)) + 
                    logh) *               # Jacobian
    diff((p-tAll) * exp(l10hpt - (s+1) * l1p10hpt))
  
  # p
  dlp <- -s * exp(log(log(10)) + logh) *
    diff(exp(l10hpt - (s+1)*l1p10hpt))
  
  # s  
  dls <- -diff(exp(log(l1p10hpt) - s*l1p10hpt))
  
  out <- matrix(c(dlh, dlp, dls), nrow=length(ti), ncol=3)
  
  return(out)
}

# Derivate seconde
d2diffRich <- function(pars, ti)
{
  # Parameters
  logh <- pars[1]
  h <- exp(logh)
  p <- pars[2]
  s <- pars[3]
  
  # All times
  tAll <- c(ti[1]-1, ti)
  
  # First derivatives
  d1lt <- d1diffRich(ti = ti, pars = pars)
  dlh <- d1lt[, 1]
  dlp <- d1lt[, 2]
  dls <- d1lt[, 3]
  
  # Auxiliary
  l1p10hptInf <- log(1+10^(h*(p-tAll)))
  l1p10hpt <- ifelse(is.infinite(l1p10hptInf), h*(p-tAll)*log(10), l1p10hptInf)
  l10hpt <- h*(p-tAll)*log(10)
  
  # Terms
  #1
  laux1 <- log(diff(exp( -s * l1p10hpt)))
  
  # Second derivatives
  # h
  dlhh <- dlh + -exp(log(s)+2*log(log(10)) + 
                       2*logh) *                      # Jacobian
    diff( exp(log((p-tAll)^2) + l10hpt + (-s-1)*l1p10hpt) * 
            (1 - exp(log(s+1) - l1p10hpt + l10hpt)) )
  dlhp <- -exp(log(s)+log(log(10)) + 
                 logh) *                      # Jacobian
    diff( exp(l10hpt + (-s-1)*l1p10hpt) * 
            (1 + (p-tAll)*h*log(10)*
               (1 - exp(log(s+1) - l1p10hpt + l10hpt))) )
  dlhs <- - exp(log(log(10)) + 
                  logh) *                       # Jacobian
    diff( (p-tAll)*exp(l10hpt + (-s-1) * l1p10hpt) * (1-exp(log(s)+log(l1p10hpt))))
  
  dlpp <- - exp(log(s) + 2*log(log(10)) + 2*log(h)) *
    diff( exp((-s-1)*l1p10hpt + l10hpt) * (1-exp(log(s+1)+l10hpt-l1p10hpt)))
  dlps <- - exp(log(log(10)) + logh) *
    diff( exp(l10hpt + (-s-1) * l1p10hpt) * (1-exp(log(s)+log(l1p10hpt))))
  
  dlss <- diff(exp(2*log(l1p10hpt) - s * l1p10hpt))
  
  # Output
  out <- list(matrix(c(dlhh, dlhp, dlhs), nrow=length(ti), ncol=3),
              matrix(c(dlpp, dlps), nrow=length(ti), ncol=2),
              matrix(c(dlss), nrow=length(ti), ncol=1))
  
  return(out)
}

# Functions - Poisson -----------------------------------------------------

# Likelihood
likPois <- function(pars, ti, di, X) 
{
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- (npars-ncovs-1)
  richpars <- pars[1:(nrich)]
  betapars <- pars[(nrich+1):(npars-1)]
  logr <- pars[npars]
  
  linPred <- exp(X%*%betapars) + exp(logr+ldiffRich(richpars, ti))
  
  out <- sum(dpois(di, linPred, log=TRUE))
  
  return(out)
}

# Gradient Richard-Poisson
PoisRichGradient <- function(pars, ti, di, X)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- (npars-ncovs-1)
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-1)]
  logr <- pars[npars]
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- eXb+exp(logr+llambdat)
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(pars, ti)
  d1eXb <- c(eXb)*X
  d1logr <- rep(exp(logr), length(ti))
  d1mut <- cbind(exp(logr)*d1lt, d1eXb, lambdat*d1logr)
  
  # Auxiliary
  aux1 <- exp(log(di)-lmut)
  
  # Gradient computation
  out <- colSums(c(aux1-1)*d1mut)
  
  return(out)
}

# Hessian Richard-Poisson
PoisRichHessian <- function(pars, ti, di, X)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- (npars-ncovs-1)
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-1)]
  logr <- pars[npars]
  r <- exp(logr)
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- eXb+exp(logr+llambdat)
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(pars, ti)
  d1eXb <- c(eXb)*X
  d1logr <- rep(exp(logr), length(ti))
  d1mut <- cbind(exp(logr)*d1lt, d1eXb, lambdat*d1logr)
  
  # Second derivatives
  d2lt <- d2diffRich(richpars, ti)
  d2eXb <- list()
  for(i in 1:ncovs)
  {
    d2eXb[[i]] <- d1eXb[,i]*as.matrix(X[,i:ncovs])
  }
  d2logr <- rep(exp(logr), length(ti))
  
  # Auxiliary
  aux1 <- (di-mut)/mut
  aux2 <- exp(log(di)-2*lmut)
  
  # Hessian computation
  out <- matrix(NA, npars, npars)
  for (i in 1:npars)
  {
    if (i<=nrich)
    {
      out[i, i:npars] <- c(colSums(c(aux1*r)*d2lt[[i]]-c(aux2*r^2)*d1lt[,i]*d1lt[,i:nrich]),
                           colSums(-c(aux2)*(r*d1lt[,i])*d1eXb), 
                           sum(c(aux1)*d1lt[,i]*d1logr-c(aux2*exp(llambdat+logr))*d1lt[,i]*d1logr))
      out[i:npars, i] <- out[i, i:npars]
    }
    
    if (i>nrich & i<npars)
    {
      k <- i-nrich
      
      out[i, i:npars] <- c(colSums(c(aux1)*d2eXb[[k]]-c(aux2)*d1eXb[,k]*d1eXb[,k:ncovs]),
                           sum(-c(aux2)*(lambdat*d1logr)*d1eXb[,k]))
      out[i:npars, i] <- out[i, i:npars]
    }
    if (i==npars)
    {
      out[npars, npars] <- sum(c(aux1*lambdat)*d2logr-c(aux2*lambdat^2)*d1logr^2)
    }
  }
  
  return(out)
}

# Functions Negative Binomial ---------------------------------------------

# Likelihood
likNB <- function(pars, ti, di, X) 
{
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- npars-ncovs-2
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-2)]
  logr <- pars[npars-1]
  n <- exp(pars[npars])
  
  linPred <- exp(X%*%betapars)+exp(logr + ldiffRich(richpars, ti))
  
  out <- sum(dnbinom(di, size=n, mu=linPred, 
                     log=TRUE))
  
  return(out)
}

# Functions Negative Binomial - Richards ---------------------------------------------

NBRichGradient <- function(pars, ti, di, X)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- npars-ncovs-2
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-2)]
  logr <- pars[npars-1]
  n <- exp(pars[npars])
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- eXb+exp(logr+llambdat)
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(richpars, ti)
  d1eXb <- c(eXb)*X
  d1logr <- rep(exp(logr), length(ti))
  d1mut <- cbind(exp(logr)*d1lt, d1eXb, lambdat*d1logr)
  d1n <- (log(n) - digamma(n) + digamma(n+di) - log(mut+n) + (mut-di)/(mut+n))*n
  
  # Auxiliary
  aux1 <- exp(log(n+di)-log(mut+n))
  aux2 <- exp(log(di)-lmut)
  
  # Gradient computation
  out <- c(colSums(-c(aux1)*d1mut+c(aux2)*d1mut), sum(d1n))
  
  return(out)
}


NBRichHessian <- function(pars, ti, di, X)
{
  # Parameters
  npars <- length(pars)
  ncovs <- ncol(X)
  nrich <- npars-ncovs-2
  richpars <- pars[1:nrich]
  betapars <- pars[(nrich+1):(npars-2)]
  logr <- pars[npars-1]
  r <- exp(logr)
  n <- exp(pars[npars])
  
  # Linear predictor
  lambdat <- exp(ldiffRich(richpars, ti))
  llambdat <- log(lambdat)
  Xb <- X%*%betapars
  eXb <- exp(Xb)
  mut <- eXb+exp(logr+llambdat)
  lmut <- log(mut)
  
  # First derivatives
  d1lt <- d1diffRich(richpars, ti)
  d1eXb <- c(eXb)*X
  d1logr <- rep(exp(logr), length(ti))
  d1mut <- cbind(exp(logr)*d1lt, d1eXb, lambdat*d1logr)
  d1n <- (log(n) - digamma(n) + digamma(n+di) - log(mut+n) + (mut-di)/(mut+n))*n
  
  # Second derivatives
  d2lt <- d2diffRich(ti = ti, pars = richpars)
  d2eXb <- list()
  for(i in 1:ncovs)
  {
    d2eXb[[i]] <- d1eXb[,i]*as.matrix(X[,i:ncovs])
  }
  d2logr <- rep(exp(logr), length(ti))
  d2n <- d1n + (1/n-trigamma(n)+trigamma(n+di)-1/(mut+n)-(mut-di)/(mut+n)^2)*n^2
  
  # Auxiliary
  aux1 <- exp(log(di)-lmut)-exp(log(di+n)-log(mut+n))
  aux2 <- exp(log(di+n)-2*log(mut+n))-exp(log(di)-2*lmut)
  auxnb <- (di-mut)*exp(-2*log(n+mut)+log(n))
  
  # Hessian computation
  out <- matrix(NA, npars, npars)
  
  # Hessian computation
  for (i in 1:npars)
  {
    if (i<=nrich)
    {
      out[i, i:(npars-1)] <- c(colSums(c(aux1*r)*d2lt[[i]]+c(aux2*r^2)*d1lt[,i]*d1lt[,i:nrich]),
                               colSums(c(aux2*r)*d1lt[,i]*d1eXb), 
                               sum(c(aux1)*d1lt[,i]*d1logr+c(aux2*exp(logr+llambdat))*d1lt[,i]*d1logr))
      out[i:(npars-1), i] <- out[i, i:(npars-1)]
    }
    if (i>nrich & i<(npars-1))
    {
      k <- i-nrich
      
      out[i, i:(npars-1)] <- c(colSums(c(aux1)*d2eXb[[k]]+c(aux2)*d1eXb[,k]*d1eXb[,k:ncovs]),
                               sum(c(aux2*lambdat)*d1eXb[,k]*d1logr))
      out[i:(npars-1), i] <- out[i, i:(npars-1)]
    }
    if (i==(npars-1))
    {
      out[npars-1, npars-1] <- sum(c(aux1*lambdat)*d2logr+c(aux2*lambdat^2)*d1logr^2)
    }
    if (i==npars)
    {
      out[i, -npars] <- colSums(c(auxnb)*d1mut)
      out[-npars, i] <- colSums(c(auxnb)*d1mut)
      out[npars, npars] <- sum(d2n)
    }
  }
  
  return(out)
}

