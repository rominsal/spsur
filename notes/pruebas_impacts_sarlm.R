impacts.sarlm <- function (obj, ..., tr = NULL, R = NULL, listw = NULL, evalues = NULL, 
          useHESS = NULL, tol = 1e-06, empirical = FALSE, Q = NULL) 
{
  if (obj$type == "error") {
    if (obj$etype == "emixed") {
      return(impactSDEM(obj))
    }
    else {
      stop("impact measures not for error models")
    }
  }
  if (is.null(listw) && !is.null(obj$listw_style) && obj$listw_style != 
      "W") 
    stop("Only row-standardised weights supported")
  rho <- obj$rho
  beta <- obj$coefficients
  s2 <- obj$s2
  if (obj$type == "sac" || obj$type == "sacmixed") 
    lambda <- obj$lambda
  usingHESS <- NULL
  iNsert <- obj$insert
  if (!is.null(R)) {
    resvar <- obj$resvar
    usingHESS <- FALSE
    irho <- 2
    drop2beta <- 1:2
    if (obj$type == "sac" || obj$type == "sacmixed") 
      drop2beta <- c(drop2beta, 3)
    if (is.logical(resvar)) {
      fdHess <- obj$fdHess
      if (is.logical(fdHess)) 
        stop("coefficient covariance matrix not available")
      usingHESS <- TRUE
      if (!iNsert) {
        irho <- 1
        drop2beta <- 1
        if (obj$type == "sac" || obj$type == "sacmixed") 
          drop2beta <- c(drop2beta, 2)
      }
    }
    if (!is.null(useHESS) && useHESS) {
      fdHess <- obj$fdHess
      if (is.logical(fdHess)) 
        stop("Hessian matrix not available")
      usingHESS <- TRUE
      if (!iNsert) {
        irho <- 1
        drop2beta <- 1
        if (obj$type == "sac" || obj$type == "sacmixed") 
          drop2beta <- c(drop2beta, 2)
      }
    }
    interval <- obj$interval
    if (is.null(interval)) 
      interval <- c(-1, 0.999)
    cat("drop2beta: ",drop2beta,"\n")
  }
  icept <- grep("(Intercept)", names(beta))
  iicept <- length(icept) > 0L
  zero_fill <- NULL
  dvars <- obj$dvars
  if (obj$type == "lag" || obj$type == "sac") {
    if (iicept) {
      P <- matrix(beta[-icept], ncol = 1)
      bnames <- names(beta[-icept])
    }
    else {
      P <- matrix(beta, ncol = 1)
      bnames <- names(beta)
    }
    p <- length(beta)
  } else if (obj$type == "mixed" || obj$type == "sacmixed") {
    if (!is.null(dvars)) 
      zero_fill <- attr(dvars, "zero_fill")
    if (iicept) {
      b1 <- beta[-icept]
    }
    else {
      b1 <- beta
    }
    if (!is.null(zero_fill)) {
      if (length(zero_fill) > 0L) {
        inds <- attr(dvars, "inds")
        b1_long <- rep(0, 2 * (dvars[1] - 1))
        b1_long[1:(dvars[1] - 1L)] <- b1[1:(dvars[1] - 
                                              1)]
        names(b1_long)[1:(dvars[1] - 1L)] <- names(b1)[1:(dvars[1] - 
                                                            1)]
        for (i in seq(along = inds)) {
          b1_long[(dvars[1] - 1L) + (inds[i] - 1L)] <- b1[(dvars[1] - 
                                                             1L) + i]
        }
        b1 <- b1_long
      }
    }
    p <- length(b1)
    if (p%%2 != 0) 
      stop("non-matched coefficient pairs")
    P <- cbind(b1[1:(p/2)], b1[((p/2) + 1):p])
    bnames <- names(b1[1:(p/2)])
  }
  n <- length(obj$residuals)
  mu <- NULL
  Sigma <- NULL
  if (!is.null(R)) {
    if (usingHESS && !iNsert) {
      mu <- c(rho, beta)
      if (obj$type == "sac" || obj$type == "sacmixed") 
        mu <- c(rho, lambda, beta)
      Sigma <- fdHess
    } else {
      mu <- c(s2, rho, beta)
      if (obj$type == "sac" || obj$type == "sacmixed") 
        mu <- c(s2, rho, lambda, beta)
      if (usingHESS) {
        Sigma <- fdHess
      }
      else {
        Sigma <- resvar
      }
    }
  }
  res <- spatialreg::intImpacts(rho = rho, beta = beta, P = P, n = n, 
                    mu = mu, Sigma = Sigma, irho = irho, drop2beta = drop2beta, 
                    bnames = bnames, interval = interval, type = obj$type, 
                    tr = tr, R = R, listw = listw, evalues = evalues, tol = tol, 
                    empirical = empirical, Q = Q, icept = icept, iicept = iicept, 
                    p = p, zero_fill = zero_fill, dvars = dvars)
  attr(res, "useHESS") <- usingHESS
  attr(res, "insert") <- iNsert
  attr(res, "iClass") <- class(obj)
  res
}
