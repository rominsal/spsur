### HAY QUE EXPORTARLO....

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