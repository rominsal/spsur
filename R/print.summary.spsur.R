print.summary.spsur <- function(z, digits = max(3L, getOption("digits") - 3L), ...)
{
  nG <- z$nG
  cat("Call:\n")
  print(z$call)
  cat("\n","\n")
  cat("Spatial SUR model type: ",z$type,"\n\n")
  for (i in 1:nG){
    cat("Equation ",i,"\n")
    printCoefmat(z$coef_table[[i]], P.value=TRUE, has.Pvalue=TRUE)
    cat("\n R-squared: ",formatC(z$R2[i+1],digits=4)," ",sep="")
    cat("\n","\n")
  }
  
  cat("\n Variance-Covariance Matrix of inter-equation residuals: \n")
  prmatrix(z$Sigma,digits=5,
        rowlab=rep("",nrow(z$Sigma)), collab=rep("",ncol(z$Sigma)))
  
  cat("\n Correlation Matrix of inter-equation residuals: \n")
  prmatrix(z$Sigma_corr,digits=3,
           rowlab=rep("",nrow(z$Sigma)), collab=rep("",ncol(z$Sigma)))

  
  # for (i in 1:nG)
  #   cat(" sigma",i,"_sq.: ",formatC(z$Sigma[i,i],digits=5,width=5),sep="")
  # cat("\nCovariances: ")
  # for (i in 1:nG)
  #   for (j in 1:nG){
  #     if(j>i) { cat(" sigma",i,j,":",formatC(z$Sigma[i,j],digits=5,width=5),sep="") }
  #   }
  # 
  # cat("\nCorrelations: ")
  # for (i in 1:nG)
  #   for (j in 1:nG){
  #     if(j>i) { cat(" phi",i,j,":",formatC(z$Sigma_corr[i,j],digits=5,width=5),sep="") }
  #   }
  
  cat("\n R-sq. pooled: ", formatC(z$R2[1],digits=4)," ",sep="")
  #for(i in 2:(nG+1))  cat("\n R-sq.(eq.",i-1,"): ",formatC(z$R2[i],digits=4)," ",sep="")
  
  if(!is.null(z$llsur))
    cat("\n Log-Likelihood: ",formatC(z$llsur, digits=6, width=6))
    
  # Se ajusta a una Chi con G*(G-1)/2 gl
  if(!is.null(z$BP))
    cat("\n Breusch-Pagan: ",formatC(z$BP, digits=5),"  p-value: (",formatC(pchisq(z$BP,df=nG*(nG-1)/2,lower.tail=FALSE),digits=3, width=4),") ",sep="")
  if(!is.null(z$LMM)) 
    cat("\n LMM: ",formatC(z$LMM, digits=5),"  p-value: (",formatC(pchisq(z$LMM,df=nG*(nG-1)/2,lower.tail=FALSE),digits=3, width=4),")\n",sep="")
  invisible(z)
}