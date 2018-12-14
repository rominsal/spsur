#' @name wald_deltas
#' @rdname wald_deltas
#'
#' @title Wald tests for spatial parameters coefficients.
#'
#' @description Wald tests for linear hypothesis about delta coefficients
#'   (spatial parameters).
#'
#' @param results An object create with \code{\link[spsur]{spsurml}}.
#' @param R Coefficient matrix for linear hypothesis about deltas.
#' @param r Vector of independent terms.
#'
#' @return statistics and p-value of test.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{spsurml}}
#' @examples
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' #################################
#' ## Estimate SUR-SAR model
#' spcsur.sar <-spsurml(Form = Tformula, data = spc, type = "sar", W= Wspc)
#' summary(spcsur.sar)
#' ## H_0: equality between lambda parameters in both equations.
#' R1 <- matrix(c(1,-1), nrow=1)
#' r1 <- matrix(0, ncol=1)
#' res1 <- wald_deltas(results = spcsur.sar, R = R1, r = r1)
#' res1$stat; res1$p_val
#' #################################
#' ## Estimate SUR-SEM model
#' spcsur.sem <-spsurml(Form = Tformula, data = spc, type = "sem", W = Wspc)
#' summary(spcsur.sem)
#' ## H_0: equality between delta parameters in both equations.
#' R2 <- matrix(c(1,-1), nrow=1)
#' r2 <- matrix(0, ncol=1)
#' res2 <- wald_deltas(results = spcsur.sem, R = R2, r = r2)
#' res2$stat; res2$p_val
#' #################################
#' ## Estimate SUR-SARAR model
#' spcsur.sarar <-spsurml(Form = Tformula, data = spc,
#'                        type = "sarar", W = Wspc)
#' summary(spcsur.sarar)
#' ## H_0: equality between lambda and delta parameters in both equations.
#' R3 <- matrix(c(1,-1,0,0,0,0,1,-1),nrow=2,ncol=4,byrow=TRUE)
#' r3 <- matrix(c(0,0), ncol=1)
#' res3 <- wald_deltas(results = spcsur.sarar, R = R3, r = r3)
#' res3$stat; res3$p_val
#' @export
wald_deltas <- function(results , R , r){
  z <- results # OBJETO QUE INCLUYE ESTIMACIÓN EN Rbetas <- z$betas
  deltas <- Matrix::Matrix(matrix(z$deltas,ncol=1))
  rownames(deltas) <- names(z$deltas)
  cov_deltas <- Matrix::Matrix(z$cov[rownames(deltas),rownames(deltas)])
  R <- Matrix::Matrix(R)
  colnames(R) <- rownames(deltas)
  r <- Matrix::Matrix(matrix(r,ncol=1))
  holg <- R %*% deltas - r
  q <- nrow(as.matrix(R))
  Wald <- as.numeric(Matrix::t(holg) %*%
                    Matrix::solve(R %*% cov_deltas %*% Matrix::t(R),holg)  )
  p_val <- pchisq(Wald,df=q,lower.tail=FALSE)
  #cat("\n R: "); print(as.matrix(R))
  #cat("\n r: "); print(as.matrix(r))
  #cat("\n statistical discrepancies: "); print(as.matrix(holg))
  cat("\n Wald stat.: ",round(Wald,3)," (",round(p_val,3),")",sep="")
  res <- list(stat = Wald,
              p_val = p_val,
              q = q,
              R = as.matrix(R),
              r = as.matrix(r),
              discr = as.matrix(holg) )
}
