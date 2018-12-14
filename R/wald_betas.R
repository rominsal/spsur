#' @name wald_betas
#' @rdname wald_betas
#'
#' @title Wald tests on beta coefficients
#'
#' @description
#' Wald tests on beta coefficients
#'
#' @param    results : An object create with \code{\link{spsurml}} or \code{\link{spsur3sls}}.
#' @param    R       : Coefficient matrix for betas.
#' @param    r       : Vector of independent terms.
#'
#' @details
#' The hypothesis H_0 is that the beta coefficients is specify in R matrix
#'
#' @return
#' statistics and p-value
#' The LR tests
#'   \tabular{ll}{
#'   \code{stat} \tab The value of LR test. \cr
#'   \code{p_val}    \tab The p-value of LR test. \cr
#'   \code{q}    \tab The degrees of freedom. \cr
#'   }
#'
#' @references
#' Mur, J., López, F., & Herrera, M. (2010). Testing for spatial effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#' \cr
#' \cr
#' López, F.A., Mur, J., & Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' \cr
#' \cr
#' López, F.A., Martínez-Ortiz, P.J., & Cegarra-Navarro, J.G. (2017). Spatial spillovers in public expenditure on a municipal level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{lr_betas_spsur}}
#'
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (nG>1; nT=1) ########
#' #################################################
#'
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## Estimate SUR-SAR model
#' spcsur.sar <-spsurml(Form = Tformula, data = spc,
#'                      type = "sar", W = Wspc)
#' summary(spcsur.sar)
#' ## H_0: equality between SMSA coefficients in both equations.
#' R1 <- matrix(c(0,0,0,1,0,0,0,-1), nrow=1)
#' r1 <- matrix(0, ncol=1)
#' Wald_beta <- wald_betas(results = spcsur.sar, R = R1, r = r1)
#'
#' ## H_0: equality between intercepts and SMSA coefficients in both equations.
#' R2 <- matrix(c(1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,-1),
#'              nrow = 2, ncol = 8, byrow = TRUE)
#' r2 <- matrix(c(0,0),ncol=1)
#' res2 <- wald_betas(results = spcsur.sar, R = R2, r = r2)
#' @export
 wald_betas <- function(results , R , r){
  z <- results # OBJETO QUE INCLUYE ESTIMACIÓN EN Rbetas <- z$betas
  betas <- Matrix::Matrix(matrix(z$betas,ncol=1))
  rownames(betas) <- names(z$betas)
  cov_betas <- Matrix::Matrix(z$cov[rownames(betas),rownames(betas)])
  R <- Matrix::Matrix(R)
  colnames(R) <- rownames(betas)
  r <- Matrix::Matrix(matrix(r,ncol=1))
  holg <- R %*% betas - r
  q <- nrow(as.matrix(R))
  Wald <- as.numeric( Matrix::t(holg) %*%
                Matrix::solve(R %*% cov_betas %*% Matrix::t(R),holg) )
  p_val <- pchisq(Wald,df=q,lower.tail=FALSE)
  # cat("\n R: "); print(as.matrix(R))
  # cat("\n r: "); print(as.matrix(r))
  # cat("\n statistical discrepancies: "); print(as.matrix(holg))
  cat("\n Wald stat.: ",round(Wald,3)," (",round(p_val,3),")",sep="")
  res <- list(stat = Wald,
              p_val = p_val,
              q = q,
              R = as.matrix(R),
              r = as.matrix(r),
              discr = as.matrix(holg) )
}
