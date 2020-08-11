#' @name methods_spsur
#' @title Methods for class spsur
#' @description The \code{anova()} function provides tables of fitted 
#'  spsur models including information criteria (AIC and BIC), 
#'  log-likelihood and degrees of freedom of each fitted model. The 
#'  argument \code{lrtest} allows to perform LR tests between nested models.
#'  The \code{plot()} function allows the user to plot both beta and spatial
#'  coefficients for all equations of the spsur model. The argument 
#'  \code{viewplot} is used to choose between interactive or non-interactive
#'  plots. The \code{print()} function is used to print short tables including the values of beta and
#'  spatial coefficients as well as p-values of significance test for each 
#'  coefficient. This can be used as an alternative to 
#'  \code{\link{summary.spsur}} when a brief output is needed. 
#'  The rest of methods works in the usual way. 
#'       
#' @param object a \code{spsur} object created by \code{\link{spsurml}},
#'            \code{\link{spsur3sls}} or \code{\link{spsurtime}}.
#' @param x similar to \code{object} argument for \code{print()} 
#'  and \code{plot} functions.
#' @param digits number of digits to show in printed tables.
#'   Default: max(3L, getOption("digits") - 3L).
#' @param lrtest logical value to compute likelihood ratio
#'   test for nested models in `anova` method. Default = \code{TRUE}
#' @param ci confidence level for the intervals in `plot` method. 
#'   Default \code{ci = 0.95}
#' @param viewplot logical value to show interactively the plots. 
#'   Default = \code{TRUE}
#' @param ... further arguments passed to or from other methods.  
#' @examples
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' spcsur.sim <-spsurml(formula = Tformula, data = spc, type = "sim")
#' ## Print Table       
#' print(spcsur.sim)
#' \donttest{
#' spcsur.slm <-spsurml(formula = Tformula, data = spc, type = "slm", 
#'                      listw = Wspc)
#' # ANOVA table and LR test for nested models:
#' anova(spcsur.sim, spcsur.slm)
#' ## Plot spatial and beta coefficients
#' # Interactive plot
#' plot(spcsur.slm)
#' # Non-interactive plot
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.slm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#'}
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }

NULL

#' @name anova
#' @rdname methods_spsur
#' @export
anova.spsur <- function(object, ..., lrtest = TRUE) {
  object <- list(object, ...)
  ancall <- sys.call()
  nmodels <- length(object)
  vtypes <- character(length = nmodels)
  vLL <- vdf <- vector(mode = "numeric", length=nmodels)
  vAIC <- vBIC <- vector(mode = "numeric", length=nmodels)
  if (lrtest) {
    vlrtest <- vpval <- vector(mode = "numeric", length=nmodels)
    vlrtest[1] <- vpval[1] <- NA
  } 
  for (i in 1:nmodels) {
    obj_i <- object[[i]]
    if (!inherits(obj_i, "spsur"))
      stop("object not a fitted spsur model")
    if (is.null(obj_i$LL))
      stop("object not fitted by maximum likelihood")
    LL <- logLik(obj_i)
    vAIC[i] <- stats::AIC(LL)
    vBIC[i] <- stats::BIC(LL)
    vLL[i] <- LL
    vdf[i] <- attr(LL, "df")
    vtypes[i] <- paste("model ",i,": ",obj_i$type, sep="")
    if (lrtest) {
      if(i>1) {
        vlrtest[i] <- 2*(vLL[i] - vLL[i-1])
        vpval[i] <- pchisq(vlrtest[i],
                           df = vdf[i] - vdf[i-1],
                           lower.tail = FALSE)
      }
    }
  }
  if (lrtest) {
    res <- data.frame(logLik = vLL,
                      df = vdf,
                      AIC = vAIC, 
                      BIC = vBIC,
                      LRtest = vlrtest,
                      `p-val` = vpval,
                      row.names = vtypes)
  } else {
    res <- data.frame(`logLik` = vLL,
                      df = vdf,
                      AIC = vAIC, 
                      BIC = vBIC,
                      row.names = vtypes)
    
  }
  class(res) <- c("anova", "data.frame")
  return(res)  
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
# #' @name print
# #' @rdname methods_spsur
# #' @export
# print.spsur <- function(x, ...) 
# {
#   if (!(x$type == "sim" || x$type == "slx" 
#         || x$type == "sarar")) {
#     G <- x$G
#     for (i in 1:G) {
#       delta <- x$deltas[i]
#       if (isTRUE(all.equal(unname(delta), x$interval[1])) 
#           || isTRUE(all.equal(unname(delta), x$interval[2]))) {
#         cat("Warning in equation ",i,"\n")
#         warning("spatial parameter on interval bound - 
#                       results should not be used")
#       }
#     }
#   } 
#   cat("\nCall:\n")
#   print(x$call)
#   cat("Type:", x$type, "\n")
#   cat("\nCoefficients:\n")
#   print(coef(x))
#   cat("\nLog likelihood:", logLik(x), "\n")
#   invisible(x)
# }

# #' @name linearHypothesis
# #' @rdname methods_spsur
# #' @importFrom car linearHypothesis linearHypothesis.default
# #' @method linearHypothesis spsur
# #' @examples
# #'   ## Linear Hypothesis (Wald test)
# #'   linearHypothesis(spcsur.slm, c("rho_1 = rho_2"))
# #'  ## Linear Hypothesis in matrix form (first the spatial coefficients)
# #'   coef(spcsur.slm)
# #'   R1 <- matrix(c(1,-1,0,0,0,0,0,0,0,0), nrow=1)
# #'   b1 <- c(0)
# #' linearHypothesis(spcsur.slm, hypothesis.matrix = R1, rhs = b1)
# #' @export
# linearHypothesis.spsur <- function(object, ...) {
#   car::linearHypothesis.default(object, ..., test = "Chisq")
# }

# #' @name lrtest
# #' @rdname methods_spsur
# #' @importFrom lmtest lrtest lrtest.default
# #' @method lrtest spsur
# #' @export
# lrtest.spsur <- function(object, ...) {
#   typespsur <- function(x) {
#     if (!is.null(x$type)) {
#       return(x$type)
#     } else stop ("object must be of spsur class")
#   }
#   lmtest::lrtest.default(object, ...,
#                          name = typespsur)
# }

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
  # For compatibility with linearHypothesis method
  idxsigma <- grepl("sigma", rownames(res))
  res <- res[!idxsigma, !idxsigma]
  res
}

#' @name print
#' @rdname methods_spsur
#' @export
print.spsur <- function(x, digits = max(3L, getOption("digits") - 3L),
                        ...) {
  if (!inherits(x, "spsur")) stop("Argument must be a spsur object")
  summx <- summary(x)
  fulltable <- summx$coef_table
  G <- length(fulltable)  
  fullrownames <- NULL
  for (i in 1:G) {
    rownames_i <- rownames(fulltable[[i]]) 
    rownames_i <- gsub(paste("_",i,sep=""), "", 
                       rownames_i)
    newrownames <- which (!(rownames_i %in% fullrownames))
    rownames_i <- rownames_i[newrownames]
    fullrownames <- c(fullrownames, rownames_i)
  }
  mtable <- matrix(NA, nrow = length(fullrownames),
                   ncol = 2*G)
  rownames(mtable) <- fullrownames
  colnames(mtable) <- paste(rep(c("coeff_","pval_"), G),
                            rep(1:G, each = 2), sep="")
  idxcol <- 1
  for (i in 1:G) {
    table_i <- fulltable[[i]][,c("Estimate", "Pr(>|t|)")]
    table_i[, 2] <- round(table_i[, 2], 
                          digits = min(digits, 3))
    rownames(table_i) <- gsub(paste("_",i,sep=""), "", 
                              rownames(table_i))
    mtable[rownames(table_i), c(idxcol,idxcol+1)] <- table_i
    idxcol <- idxcol + 2
  }
  if (any(grepl("rho",rownames(mtable)))) {
    idxrowrho <- which(grepl("rho",
                             rownames(mtable)))
    rowrho <- matrix(mtable[idxrowrho, ], nrow=1)
    rownames(rowrho) <- c("rho")
    mtable <- mtable[-idxrowrho, ]
    mtable <- rbind(mtable, rowrho)
  }  
  if (any(grepl("lambda",rownames(mtable)))) {
    idxrowlambda <- which(grepl("lambda",
                                rownames(mtable)))
    rowlambda <- matrix(mtable[idxrowlambda, ], nrow=1)
    rownames(rowlambda) <- c("lambda")
    mtable <- mtable[-idxrowlambda, ]
    mtable <- rbind(mtable, rowlambda)
  }
  print(round(mtable, digits = digits))
  invisible(x)
}

#' @name plot
#' @rdname methods_spsur
#' @export
plot.spsur <- function(x, ci = 0.95, viewplot = TRUE, ...) {
  if (!inherits(x, "spsur")) stop("Argument must be a spsur object")
  critval <- qnorm((1-ci)/2, lower.tail = FALSE)
  G <- x$G
  p <- x$p
  eq <- NULL
  for (i in 1:G) {
    eq <- c(eq, rep(i,p[i]))
  }
  eq <- as.factor(eq)
  dfx <- data.frame(eq = eq, 
                    nbetas = names(x$coefficients),
                    betas = x$coefficients,
                    sebetas = x$rest.se,
                    lwbetas = x$coefficients - critval*x$rest.se,
                    upbetas = x$coefficients + critval*x$rest.se,
                    row.names = NULL)
  lplbetas <- vector("list", G)
  for (i in 1:G) {
    plbetasi <- ggplot2::ggplot(dfx[dfx$eq==i, ], 
                                ggplot2::aes(.data$nbetas, 
                                             .data$betas, 
                                             colour = i)) + 
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lwbetas, 
                                            ymax = .data$upbetas),
                               show.legend = FALSE) +
      ggplot2::labs(title = paste("Equation ",i), 
                    x = "", y = "betas")
    lplbetas[[i]] <- plbetasi
    if (viewplot) {
      gridExtra::grid.arrange(lplbetas[[i]], newpage = TRUE)
      readline(prompt="Press [enter] to continue")
    }
  }
  if(!is.null(x$deltas)) {
    eq <- rep(0, length(x$deltas))
    for (i in 1:G) {
      ni <- i*as.integer(grepl(paste("_", i, sep=""), names(x$deltas)))
      eq <- eq + ni
    }
    eq <- as.factor(eq)
    dfxsp <- data.frame(eq = eq,
                        ndeltas = names(x$deltas),
                        deltas = x$deltas,
                        sedeltas = x$deltas.se,
                        lwdeltas = x$deltas - critval*x$deltas.se,
                        updeltas = x$deltas + critval*x$deltas.se,
                        row.names = NULL)
    pldeltas <- ggplot2::ggplot(dfxsp, ggplot2::aes(.data$ndeltas,
                                                    .data$deltas, 
                                                    colour = eq)) + 
      ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$lwdeltas, 
                                            ymax = .data$updeltas),
                               show.legend = FALSE) +
      ggplot2::labs(title = "", x = "", 
                    y = "Spatial Coefficients")  
    if (viewplot) gridExtra::grid.arrange(pldeltas, nrow = 1)
  } else pldeltas <- NULL
  plots <- list(lplbetas = lplbetas, 
                pldeltas = pldeltas)
}
