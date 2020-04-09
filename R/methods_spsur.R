#' Methods for class spsur
#'
#' @name methods_spsur
#'   
#' @param object An \code{spsur} object created by \code{\link{spsurml}},
#'            \code{\link{spsur3sls}} or \code{\link{spsurtime}}.
#' @param x Same than *object* for `print` method.
#' @param ... further arguments passed to or from other methods.          
NULL

#' @name anova
#' @rdname methods_spsur
#' @export
anova.spsur <- function(object, ...) {
  if (length(list(object, ...)) > 1L) {
    getRespFormula <- function(object) {
      form <- formula(object$call)
      if (!(inherits(form, "formula") && (length(form) == 
                                          3L))) {
        stop("\"Form\" must be a two sided formula")
      }
      eval(parse(text = paste("~", deparse(form[[2]]))))
    }
    object <- list(object, ...)
    ancall <- sys.call()
    nmodels <- length(object)
    if (nmodels == 1) 
      return(anova(object))
    termsClass <- unlist(lapply(object, data.class))
    if (!all(match(termsClass, c("spsur"), 
                   0))) {
      stop(paste("Objects must inherit from class spsur"))
    }
    resp <- unlist(lapply(object, 
                          function(el) deparse(getRespFormula(el)[[2]])))
    subs <- as.logical(match(resp, resp[1], FALSE))
    if (!all(subs)) 
      warning(paste("Some fitted objects deleted because", 
                    "response differs from the first model"))
    if (sum(subs) == 1) 
      stop("First model has a different response from the rest")
    object <- object[subs]
    aux <- lapply(object, logLik)
    if (length(unique(unlist(lapply(object, 
                                    function(el) length(residuals(el)))))) > 1L) {
      stop("All fitted objects must use the same number of observations")
    }
    dfModel <- unlist(lapply(aux, function(el) attr(el, "df")))
    logLik <- unlist(lapply(aux, function(el) c(el)))
    AIC <- unlist(lapply(aux, AIC))
    aod <- data.frame(Model = (1:nmodels), df = dfModel, 
                      AIC = AIC, logLik = logLik, check.names = FALSE)
    ddf <- diff(dfModel)
    if (sum(abs(ddf)) > 0) {
      effects <- rep("", nmodels)
      for (i in 2:nmodels) {
        if (ddf[i - 1] != 0) {
          effects[i] <- paste(i - 1, i, sep = " vs ")
        }
      }
      pval <- rep(NA, nmodels - 1)
      ldf <- as.logical(ddf)
      lratio <- 2 * abs(diff(logLik))
      lratio[!ldf] <- NA
      pval[ldf] <- 1 - pchisq(lratio[ldf], abs(ddf[ldf]))
      aod <- data.frame(aod, Test = effects, L.Ratio = c(NA, 
                                                         lratio),
                        `p-value` = c(NA, pval), check.names = FALSE)
    }
    row.names(aod) <- unlist(lapply(as.list(ancall[-1]), 
                                    deparse))
    attr(aod, "nmodels") <- nmodels
    class(aod) <- c("anova", "data.frame")
    return(aod)
  } else {
    if (!inherits(object, "spsur")) 
      stop("object not a fitted spsur model")
    LL <- logLik(object)
    AIC <- AIC(LL)
    res <- data.frame(AIC = AIC, `Log likelihood` = LL, 
                      df = attr(LL, "df"), 
                      row.names = deparse(substitute(object)))
    class(res) <- c("anova", "data.frame")
    return(res)
  }
}

#' @name coef
#' @rdname methods_spsur
#' @export
coef.spsur <- function(object, ...) 
{
  ret <- NULL
  if (!(object$type == "sim" || object$type == "slx"))
    ret <- c(ret, object$deltas)
  ret <- c(ret, object$coefficients)
  ret
}

#' @name fitted
#' @rdname methods_spsur
#' @export
fitted.spsur <- function(object, ...) 
{
  if (is.null(object$na.action)) 
    fits <- object$fitted.values
  else fits <- napredict(object$na.action, object$fitted.values)
  G <- object$G
  if (!(length(fits) %% G == 0)) 
    stop("Length of fitted values incompatible with number of equations")
  N <- length(fits) %/% G
  res <- vector("list", G)
  indx <- 1
  for (i in 1:G) {
    res[[i]] <- fits[indx:(indx + N - 1)]
    indx <- indx + N
  }
  names(res) <- paste0("fitted_values_eq",1:G)
  res
}


#' @name logLik
#' @rdname methods_spsur
#' @export
logLik.spsur <- function(object, ...) 
{
  LL <- c(object$LL)
  class(LL) <- "logLik"
  N <- length(residuals(object))
  attr(LL, "nall") <- N
  attr(LL, "nobs") <- N
  attr(LL, "df") <- object$parameters
  LL
}

#' @name print
#' @rdname methods_spsur
#' @export
print.spsur <- function(x, ...) 
{
  if (!(x$type == "sim" || x$type == "slx" 
        || x$type == "sarar")) {
    G <- x$G
    for (i in 1:G) {
      delta <- x$deltas[i]
      if (isTRUE(all.equal(unname(delta), x$interval[1])) 
          || isTRUE(all.equal(unname(delta), x$interval[2]))) {
        cat("Warning in equation ",i,"\n")
        warning("spatial parameter on interval bound - 
                      results should not be used")
      }
    }
  } 
  cat("\nCall:\n")
  print(x$call)
  cat("Type:", x$type, "\n")
  cat("\nCoefficients:\n")
  print(coef(x))
  cat("\nLog likelihood:", logLik(x), "\n")
  invisible(x)
}

#' @name residuals
#' @rdname methods_spsur
#' @export
residuals.spsur <- function(object, ...) 
{
  if (is.null(object$na.action)) 
    resids <- object$residuals
  else resids <- napredict(object$na.action, object$residuals)
  G <- object$G
  if (!(length(resids) %% G == 0)) 
    stop("Length of residuals incompatible with number of equations")
  N <- length(resids) %/% G
  res <- vector("list", G)
  indx <- 1
  for (i in 1:G) {
    res[[i]] <- resids[indx:(indx + N - 1)]
    indx <- indx + N
  }
  names(res) <- paste0("residuals_eq",1:G)
  res
}

#' @name vcov
#' @rdname methods_spsur
#' @export
vcov.spsur <- function(object, ...) 
{
  res <- object$resvar
  res
}