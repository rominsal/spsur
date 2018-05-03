summary.spsur <- function(object,...)
{

  z <- object
  stopifnot(inherits(z, "spsur"))
  nG <- z$nG; nR <- z$nR; nT <- z$nT; p <- z$p
  rdf <- z$df.residual
  r <- z$residuals
  f <- z$fitted.values
  rss <- sum(r^2)
  resvar <- rss/rdf
  betas <- z$betas
  se_betas <- z$se_betas
  t_betas <- betas / se_betas
  deltas <- z$deltas
  se_deltas <- z$se_deltas
  t_deltas <- deltas / se_deltas
  z$coef_table <- list(NULL)
 # Build coefficients table by Equation
  for (i in 1:nG)
  {
    if(i==1){
      z$coef_table[[i]] <- cbind(betas[1:p[i]], se_betas[1:p[i]],
                                         t_betas[1:p[i]],
                                         2 * pt(abs(t_betas[1:p[i]]),rdf,
                                                lower.tail = FALSE))
      colnames(z$coef_table[[i]] ) <- c("Estimate", "Std. Error", 
                                        "t value", "Pr(>|t|)")
    } else {
      z$coef_table[[i]] <-  cbind(betas[(cumsum(p)[i-1]+1):cumsum(p)[i]], 
                                         se_betas[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                         t_betas[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                         2 * pt(abs(t_betas[(cumsum(p)[i-1]+1):cumsum(p)[i]]),rdf,
                                                lower.tail = FALSE))

      }
    if(!is.null(deltas)){
      z$coef_table[[i]] <-  rbind(z$coef_table[[i]],
                                  cbind(deltas[i],se_deltas[i],t_deltas[i],
                                        2 * pt(abs(t_deltas[i]),rdf,
                                               lower.tail = FALSE)) )
      if(z$type=="sarar") {
        z$coef_table[[i]] <-  rbind(z$coef_table[[i]],
                                    cbind(deltas[i+nG],
                                          se_deltas[i+nG],t_deltas[i+nG],
                                          2 * pt(abs(t_deltas[i+nG]),rdf,
                                                 lower.tail = FALSE)) )
      }
    } 
    colnames(z$coef_table[[i]]) <- c("Estimate", "Std. Error", 
                                     "t value", "Pr(>|t|)")    
  }
  class(z) <- c("summary.spsur",class(z))
  z
}
