#
# COVINDEX based on a GAM Beta regression model ----
#


load_install_package <- function(pkg)
{
# Load a package and install if not already installed
  if(!suppressWarnings(require(pkg, character.only = TRUE)) )
  { 
    install.packages(pkg)
    require(pkg, character.only = TRUE)
  }
}

load_install_package("mgcv")
load_install_package("data.table")

covindex <- function(x) x/data.table::shift(x, n = 7)

covindex_gam_betareg <- function(y, t = NULL, 
                                 x = NULL, wts = NULL, 
                                 k = 3:51, bs = "tp",
                                 criterion = c("AIC", "BIC"))
{
# Beta GAM model for covindex
# y = a vector of daily test positive rate values
# t = a vector of times 
# x = a data.frame of named further covariates
# wts = a vector of swabs to be used as weights
# k = a vector of integers to search for the optimal dimension of the 
#     spline basis used to represent the smooth term

  y <- as.vector(y)
  n <- length(y)
  if(is.null(t)) t <- seq.int(1,n)
  t <- as.vector(t)
  stopifnot(length(t) == n)
  if(!is.null(x)) 
    x <- as.data.frame(x)
  k <- as.integer(k)
  if(is.null(wts)) wts <- rep(1, n)
  if(any(y == 0 | y == 1, na.rm = TRUE))  
  { # if needed apply transformation for boundary correction
    y <- (y*(n-1) + 0.5)/n  
  }
  data <- data.frame(y, t, wts)
  if(!is.null(x)) 
    data <- cbind(data, x)
  
  # model selection
  criterion <- match.arg(criterion, choices = c("AIC", "BIC"))
  crit <- matrix(as.double(NA), nrow = length(k), ncol = 2, 
                 dimnames = list(NULL, c("AIC", "BIC")))
  if(length(k) > 1)
  {
    # select smoothing by criterion (AIC or BIC)
    for(i in seq(k))
    {
      fmla <- as.formula(paste0("y ~ ", 
                                paste0("s(t, k = ", k[i], 
                                       " , bs = '", bs, "')")))
      if(!is.null(x)) 
         fmla <- update(fmla, paste("~ . +", colnames(x), collapse = " + "))
      mod <- gam(fmla, 
                 data = data,
                 family = betar(link="logit"), 
                 weights = wts,
                 method = "REML")
      crit[i,] <- c(AIC(mod), BIC(mod))
    }
    k_opt <- if(criterion == "AIC") k[which.min(crit[,1])] 
             else                   k[which.min(crit[,2])] 
    names(k_opt) <- criterion
  } else
  {
    k_opt <- k
  }
  
  # fit best model
  fmla <- as.formula(paste0("y ~ ", 
                            paste0("s(t, k = ", k_opt, 
                                   " , bs = '", bs, "')")))
  if(!is.null(x)) 
    fmla <- update(fmla, paste("~ . +", colnames(x), collapse = " + "))
  mod <- gam(fmla, 
             data = data,
             family = betar(link="logit"), 
             weights = wts,
             method = "REML")
  mod$k <- k
  mod$AIC <- crit[,1]
  mod$BIC <- crit[,2]
  mod$k_opt <- k_opt 
  class(mod) <- c("covindex_gam_betareg", class(mod))
  return(mod)
}

# Simulation-based credible and prediction intervals
# References:
# Gelman Hill (2007) Data Analysis Using Regression and Multilevel 
#   Hierarchical Models, Sec. 7.2
# Gelman Hill Vehtari (2020) Regression and Other Stories, Sec. 14.3

simpred_covindex_gam_betareg <- function(object, newdata, 
                                         nsim = 1e4, level = 0.95, ...)
{

  if(missing(newdata) || is.null(newdata)) 
    newdata <- object$model

  conf_levels <- c((1-level)/2, (1+level)/2)
  
  linkinv <- object$family$linkinv
  mu  <- predict(object, newdata = newdata, type = "response")
  phi <- exp(object$family$getTheta())
  rndbeta <- object$family$rd
  
  # prediction intervals
  beta_sim <- rmvn(nsim, coef(object), vcov(object, 
		freq=FALSE, unconditional=TRUE))
  Xlp <- predict(object, newdata = newdata, type = "lpmatrix")
  mu_sim  <- linkinv(Xlp %*% t(beta_sim))
  credint <- apply(mu_sim, 1, quantile, 
                   prob = conf_levels, na.rm = TRUE)
  # prediction intervals
  y_sim <- apply(mu_sim, 2, function(m) rndbeta(m))
  predint <- apply(y_sim, 1, quantile, prob = conf_levels, na.rm = TRUE)
  
  out <- list(mu = mu, phi = phi, 
              level = level, nsim = nsim,
              mu_sim = mu_sim,
              y_sim = y_sim,
              lower_credint = credint[1,], 
              upper_credint = credint[2,],
              lower_predint = predint[1,], 
              upper_predint = predint[2,])
  return(out)
}

minor_breaks = function(x, digits = 0)
{
  x = round(x[is.finite(x)], digits)
  seq(min(x), max(x), by = 1/10^digits)
}

monthyear_labels <- function(x, ...) 
{
  fmt <- rep("%b", length(x))
  fmt[which(!is.na(x))[1]] <- "%b\n%Y"
  fmt[which(diff(year(x)) > 0)+1] <- "%b\n%Y"
  format(x, fmt)
}
