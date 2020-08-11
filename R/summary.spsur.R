#' @name summary.spsur
#' @rdname summary.spsur
#'
#' @title Summary of estimated objects of class \emph{spsur}.
#'
#' @description  This function summarizes estimated \emph{spsur} objects. The tables in the output
#'  include basic information for each equation. The report also shows other complementary results
#'  corresponding to the SUR model like the \emph{(GxG)} covariance matrix of the residuals of the
#'  equations of the SUR, the estimated log-likelihood, the Breusch-Pagan diagonality test or the Marginal
#'  Lagrange Multiplier, LMM, tests of spatial dependence.
#'  
#' @param object An \emph{spsur} object estimated using \code{\link{spsurml}},
#'  \code{\link{spsur3sls}} or \code{\link{spsurtime}} functions.
#' @param ... further arguments passed to or from other methods.
#'
#' @return An object of class \emph{summary.spsur}
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @seealso
#' \code{\link{print.summary.spsur}}; \code{\link{spsurml}}; \code{\link{spsur3sls}}.
#'
#'
#' @examples
#'  # See examples for \code{\link{spsurml}} or
#'  # \code{\link{spsur3sls}} functions.
#'
#' @export
summary.spsur <- function(object, ...) {
  z <- object
  stopifnot(inherits(z, "spsur"))
  G <- z$G; N <- z$N; Tm <- z$Tm; p <- z$p 
  rdf <- z$df.residual
  r <- z$residuals
  f <- z$fitted.values
  rss <- sum(r^2)
  resvar <- rss/rdf
  betas <- z$coefficients
  se_betas <- z$rest.se
  t_betas <- betas / se_betas
  deltas <- z$deltas
  se_deltas <- z$deltas.se
  t_deltas <- deltas / se_deltas
  Sigma <- z$Sigma
  allSigmas <- get_Sigma(r, N, G, Tm)
  z$Sigma_corr <- allSigmas$Sigma_corr
  z$Sigma_inv <- allSigmas$Sigma_inv
  rm(allSigmas)
  z$coef_table <- list(NULL)
 # Build coefficients table by Equation
  for (i in 1:G)
  {
    if (i == 1) {
      z$coef_table[[i]] <- cbind(betas[1:p[i]], se_betas[1:p[i]],
                                         t_betas[1:p[i]],
                                         2 * pt(abs(t_betas[1:p[i]]), rdf,
                                                lower.tail = FALSE))
      colnames(z$coef_table[[i]] ) <- c("Estimate", "Std. Error",
                                        "t value", "Pr(>|t|)")
    } else {
      z$coef_table[[i]] <-  cbind(betas[
                                   (cumsum(p)[i-1]+1):cumsum(p)[i]],
                                   se_betas[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                   t_betas[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                   2 * pt(abs(t_betas[(cumsum(p)[i-1]+1):
                                                        cumsum(p)[i]]),rdf,
                                                        lower.tail = FALSE))
    }
    if (!is.null(deltas)) {
      z$coef_table[[i]] <-  rbind(z$coef_table[[i]],
                                  cbind(deltas[i],se_deltas[i],t_deltas[i],
                                        2 * pt(abs(t_deltas[i]),rdf,
                                               lower.tail = FALSE)) )
      if (any(z$typ == c("sarar","gnm"))) {
        z$coef_table[[i]] <-  rbind(z$coef_table[[i]],
                                    cbind(deltas[i+G],
                                          se_deltas[i+G],t_deltas[i+G],
                                          2 * pt(abs(t_deltas[i+G]),rdf,
                                                 lower.tail = FALSE)) )
      }
    }
    colnames(z$coef_table[[i]]) <- c("Estimate", "Std. Error",
                                     "t value", "Pr(>|t|)")
  }
  class(z) <- c("summary.spsur", class(z))
  z
}
