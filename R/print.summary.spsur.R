#' @name print.summary.spsur
#' @rdname print.summary.spsur
#'
#' @title Print method for objects of class summary.spsur.
#'
#' @param x object of class \emph{summary.spsur}.
#' @param digits number of digits to show in printed tables.
#'   Default: max(3L, getOption("digits") - 3L).
#' @param ... further arguments passed to or from other methods.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @seealso
#' \code{\link{summary.spsur}}.
#'
#'
#' @examples
#'  # See examples for \code{\link{spsurml}} or
#'  # \code{\link{spsur3sls}} functions.
#' @export
print.summary.spsur <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  G <- x$G
  cat("Call:\n")
  print(x$call)
  cat("\n","\n")
  cat("Spatial SUR model type: ",x$type,"\n\n")
  for (i in 1:G){
    cat("Equation ",i,"\n")
    printCoefmat(x$coef_table[[i]], P.values = TRUE, has.Pvalue = TRUE)
    cat("R-squared: ", formatC(x$R2[i+1], digits = 4), " ", sep = "")
    cat("\n"," ")
  }
  cat("\n")
  
  if (x$Tm>1 | x$G>1){
  cat("Variance-Covariance Matrix of inter-equation residuals:")
  prmatrix(x$Sigma,digits=4,
        rowlab=rep("", nrow(x$Sigma)), collab = rep("", ncol(x$Sigma)))
  } else {
  cat("Residual standard error:",formatC(sqrt(x$Sigma), digits = 4, width = 6))
  }
    
  if (x$Tm>1 | x$G>1){
  cat("Correlation Matrix of inter-equation residuals:")
  prmatrix(x$Sigma_corr,digits=3,
           rowlab=rep("",nrow(x$Sigma)), collab = rep("",ncol(x$Sigma)))
  }
  if (x$Tm>1 | x$G>1){
  cat("\n R-sq. pooled: ", formatC(x$R2[1], digits = 4)," ", sep = "")
  if(!is.null(x$llsur)){
    cat("\n Log-Likelihood: ",formatC(x$llsur, digits = 6, width = 6))
  }
  }
  # Se ajusta a una Chi con G*(G-1)/2 gl
  if (x$Tm>1 | x$G>1){
  # Only report the BP test for multiequations
  if(!is.null(x$BP)){
    cat("\n Breusch-Pagan: ",formatC(x$BP, digits = 4),
        "  p-value: (",formatC(pchisq(x$BP, df = G*(G - 1)/2, lower.tail = FALSE),
                               digits = 3, width = 4),") ", sep = "")

  }
  }
  if(!is.null(x$LMM)){
    # If Tm=G=1 the LMM test have df=1
    if (x$Tm==1 & x$G==1){
      df= 1}
        else {df= G*(G - 1)/2
    }
    
    cat("\n LMM: ",formatC(x$LMM, digits = 5),
        "  p-value: (",formatC(pchisq(x$LMM, df = df, 
                                      lower.tail = FALSE),
                               digits = 3, width = 4),")\n", sep = "")

  }
  invisible(x)
}
