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
                                ...){
  G <- x$G
  cat("Call:\n")
  print(x$call)
  cat("\n","\n")
  cat("Spatial SUR model type: ",x$type,"\n\n")
  for (i in 1:G){
    cat("Equation ",i,"\n")
    printCoefmat(x$coef_table[[i]], P.values = TRUE, has.Pvalue = TRUE)
    cat("R-squared: ",formatC(x$R2[i+1], digits=4), " ", sep="")
    cat("\n"," ")
  }

  cat("Variance-Covariance Matrix of inter-equation residuals:")
  prmatrix(x$Sigma,digits=5,
        rowlab=rep("",nrow(x$Sigma)), collab=rep("",ncol(x$Sigma)))

  cat("Correlation Matrix of inter-equation residuals:")
  prmatrix(x$Sigma_corr,digits=3,
           rowlab=rep("",nrow(x$Sigma)), collab=rep("",ncol(x$Sigma)))


  cat("\n R-sq. pooled: ", formatC(x$R2[1],digits=4)," ",sep="")
  if(!is.null(x$llsur)){
    cat("\n Log-Likelihood: ",formatC(x$llsur, digits=6, width=6))
  }

  # Se ajusta a una Chi con G*(G-1)/2 gl
  if(!is.null(x$BP)){
    cat("\n Breusch-Pagan: ",formatC(x$BP, digits = 4),
        "  p-value: (",formatC(pchisq(x$BP,df=G*(G-1)/2,lower.tail=FALSE),
                               digits=3, width=4),") ",sep="")

  }
  if(!is.null(x$LMM)){
    cat("\n LMM: ",formatC(x$LMM, digits=5),
        "  p-value: (",formatC(pchisq(x$LMM,df=G*(G-1)/2,lower.tail=FALSE),
                               digits=3, width=4),")\n",sep="")

  }
  invisible(x)
}
