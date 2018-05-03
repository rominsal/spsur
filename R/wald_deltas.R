#' @name wald_deltas
#' @rdname wald_deltas
#'
#' @title Wald tests for spatial parameters coefficients
#'
#' @description
#' Wald tests for linear hypothesis about delta coefficients (spatial parameters)
#'
#' @param    results : An object create with \code{\link[spSUR2]{spsurml}}.
#' @param    R       : Coefficient matrix for deltas. 
#' @param    r       : Vector of independent terms.
#' 
#' @details
#' ¿?¿?
#'
#' @return
#' statistics and p-value 
#' 
#' @references
#' 
#' CAMBIAR
#' J. LeSage and R.K. Pace. \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10.1.6, 2009
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
#'   Fernando López  \tab \email{fernando.lopez@upct.es} \cr
#'   Román Mínguez  \tab \email{Roman.Minguez@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{spsurml}}
#' @examples
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ################################# 
#' ## Estimate SUR-SAR model 
#' spcSUR.sar <-spsurml(Form=Tformula,data=spc,type="sar",W=Wspc)
#' summary(spcSUR.sar)
#' ## H_0: equality between lambda parameters in both equations. 
#' R1 <- matrix(c(1,-1),nrow=1)
#' r1 <- matrix(0,ncol=1)
#' res1 <- wald_deltas(results=spcSUR.sar,R=R1,r=r1)
#' res1$stat; res1$p_val    
#' #################################
#' ## Estimate SUR-SEM model 
#' spcSUR.sem <-spsurml(Form=Tformula,data=spc,type="sem",W=Wspc)
#' summary(spcSUR.sem)
#' ## H_0: equality between delta parameters in both equations. 
#' R2 <- matrix(c(1,-1),nrow=1)
#' r2 <- matrix(0,ncol=1)
#' res2 <- wald_deltas(results=spcSUR.sem,R=R2,r=r2)
#' res2$stat; res2$p_val    
#' ################################# 
#' ## Estimate SUR-SARAR model 
#' spcSUR.sarar <-spsurml(Form=Tformula,data=spc,type="sarar",W=Wspc)
#' summary(spcSUR.sarar)
#' ## H_0: equality between lambda and delta parameters in both equations. 
#' R3 <- matrix(c(1,-1,0,0,0,0,1,-1),nrow=2,ncol=4,byrow=TRUE)
#' r3 <- matrix(c(0,0),ncol=1)
#' res3 <- wald_deltas(results=spcSUR.sarar,R=R3,r=r3)
#' res3$stat; res3$p_val    


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
  cat("\n R: "); print(as.matrix(R))
  cat("\n r: "); print(as.matrix(r))
  cat("\n statistical discrepancies: "); print(as.matrix(holg))
  cat("\n Wald stat.: ",round(Wald,3)," p-value (",p_val,")")
  res <- list(stat = Wald,
              p_val = p_val,
              q = q,
              R = as.matrix(R),
              r = as.matrix(r),
              discr = as.matrix(holg) )
}  