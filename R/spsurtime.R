#' @name spsurtime
#' @rdname spsurtime
#'
#' @title Estimation of spatio-temporal SUR model either using maximum
#'        likelihood (ml) or three stages least squares (3sls).
#'
#' @description Estimation of spatio-temporal SUR-LAG; SUR-ERR; SUR-SLX;
#'   SUR-SDM and SUR-SARAR models. Assumption: \emph{nR} and \emph{nT}>1
#'   and \emph{nG} = 1. The correlation matrix of the noise only gathers
#'   temporal covariances. The database is a spatio-temporal panel.
#'
#' @inheritParams spsurml
#' @param time Variable including temporal periods.
#' @param method Method of estimation either "ml" or "3sls". Default = "ml".
#' @param maxlagW Maximum order of the instruments in 3sls. Default = 2.
#' @param trace Logical value to get intermediate results.
#'   Default = \code{TRUE}.
#'
#' @details The function \code{\link{spsurtime}} estimate a spatial SUR model
#'   with multiple temporal equations for the same cross-section.
#'   The estimation allows for temporal correlations.
#'
#' @return Results of ML estimation of a SUR model with spatial effects.
#'   A list including:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of specified model. \cr
#'     \code{betas} \tab Fitted coefficients for the covariates. \cr
#'     \code{deltas} \tab Fitted spatial autocorrelation coefficients. \cr
#'     \code{se_betas} \tab Standard errors of beta. \cr
#'     \code{se_deltas} \tab Standard errors of delta. \cr
#'     \code{cov} \tab Fitted covariance matrix of betas and deltas.\cr
#'     \code{llsur} \tab Value of likelihood function at maximum. \cr
#'     \code{R2} \tab Determination coefficient between the explained
#'       variable and its predictions. Global and for each equation. \cr
#'     \code{Sigma} \tab Fitted covariance matrix of residuals between
#'       equations. \cr
#'     \code{Sigma_corr} \tab Fitted correlation matrix of residuals
#'       between equations. \cr
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom of the residuals. \cr
#'     \code{fitted.values} \tab Fitted values of the dependent variables. \cr
#'     \code{BP} \tab Breusch-Pagan test. Null hypothesis: diagonality of
#'       \code{Sigma} matrix. \cr
#'     \code{LMM} \tab Marginal test LM(\eqn{\rho}|\eqn{\lambda}) or
#'       LM(\eqn{\lambda}|\eqn{\rho}) of spatial autocorrelation in the
#'       residuals. \cr
#'     \code{nG} \tab Number of equations. \cr
#'     \code{nR} \tab Number of cross-section or spatial observations. \cr
#'     \code{nT} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation
#'       (including independent terms). \cr
#'     \code{demean} \tab Logical value used for demeaning. \cr
#'     \code{Y} \tab Vector \emph{Y} of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the SUR model. \cr
#'     \code{W} \tab Spatial weight matrix. \cr
#'     \code{Sigma_inv} \tab Inverse of \code{Sigma}. \cr#'
#'   }
#' @references
#'   \itemize{
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'      \item López, F.A., Martínez-Ortiz, P.J., & Cegarra-Navarro, J.G. (2017).
#'        Spatial spillovers in public expenditure on a municipal level in
#'        Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsur3sls}}
#' @examples
#'
#' ####################################
#' ######## PANEL DATA (nG=1; nT>1) ###
#' ####################################
#'
#' data("unemp_it")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#'
#' ## SUR-SAR by Maximum Likelihood
#' unempitml_sar <- spsurtime(Form = form_un, data = unemp_it,
#'                            time = unemp_it$year, W= W_italy,
#'                            type = "sar", method = "ml")
#' summary(unempitml_sar)
#'
#' ## SUR-SEM by Maximum Likelihood
#' unempitml_sem <- spsurtime(Form = form_un, data = unemp_it,
#'                            time = unemp_it$year, W = W_italy,
#'                            type = "sem", method = "ml")
#' summary(unempitml_sem)
#'
#' ## SUR-SAR by Three-Stages Least Squares Estimation
#' unempit3sls_sar <- spsurtime(Form = form_un, data = unemp_it,
#'                              time = unemp_it$year, W = W_italy,
#'                              type = "sar", method = "3sls")
#' summary(unempit3sls_sar)
#'
#' ## SUR-SDDM by Three-Stages Least Squares Estimation
#' unempit3sls_sdm <- spsurtime(Form = form_un, data = unemp_it,
#'                              time = unemp_it$year, W = W_italy,
#'                              type = "sdm", method = "3sls")
#' summary(unempit3sls_sdm)

#' ############################ LR tests about betas in spatio-temporal models
#' ## H0: equal emprgrowth beta in equations 1, 3, 4 and 5
#' R <- matrix(0,nrow=2,ncol=30)
#' R[1,2] <- 1; R[1,14] <- -1
#' R[2,2] <- 1; R[2,20] <- -1
#' r <- matrix(0,nrow=2,ncol=1)
#' lr_partrate <-  lr_betas_spsur(Form = form_un, data = unemp_it,
#'                                time = unemp_it$year, W = W_italy,
#'                                type = "sar", R = R, r = r, trace = TRUE,
#'                                printmodels = TRUE)
#'
#' ############################ Wald tests about betas in spatio-temporal models
#' wald_betas(unempitml_sar , R = R , r = r)

#' ############################ Wald tests about spatial-parameters in
#' ############################ spatio-temporal models
#' ## H0: equal lambdas in sar model for equations 1, 2, 3.
#' R2 <- matrix(0,nrow=2,ncol=5)
#' R2[1,1] <- 1; R2[1,2] <- -1
#' R2[2,1] <- 1; R2[2,3] <- -1
#' r2 <- matrix(0,nrow=2,ncol=1)
#' wald_deltas(unempitml_sar, R = R2, r = r2)
#' @export
spsurtime <- function(Form, data, time, type = "sim",  method = "ml",
                      maxlagW = 2, W = NULL, cov = TRUE, trace = TRUE,
                      R = NULL, r = NULL) {
  if (class(time) != "factor") time <- as.factor(time)
  time <- droplevels(time)
  if (length(time) != nrow(data)) stop("time must have same length than the
                                       number of rows in data")
  mt <- terms(Form)
  nG <- length(levels(time))
  Ylist <- vector("list",nG)
  Xlist <- vector("list",nG)
  p <- NULL
  namesX <- NULL
  levels_time <- levels(time)
  for (i in 1:nG) {
    data_i <- model.frame(mt,data=data[time==levels_time[i],])
    Ylist[[i]] <- data_i[,1]
    Xlist[[i]] <- model.matrix(mt,data=data[time==levels_time[i],])
    p <- c(p,ncol(Xlist[[i]]))
    namesX <- c(namesX,paste(colnames(Xlist[[i]]),i,sep="_"))
  }
  Y <- matrix(unlist(Ylist),ncol=1)
  X <- as.matrix(Matrix::bdiag(Xlist))
  colnames(X) <- namesX
  nR <- length(Ylist[[1]]); nT <- 1
  if (method == "ml"){
    res <- spsurml(X = X, Y = Y, W = W, nG = nG, nR = nR, nT = nT,
                   p = p, R = R, r = r, type = type, cov = cov,
                   control = list(tol = 1e-3, maxit = 200,
                                  trace = trace) )
  }
  if (method == "3sls"){
    res <- spsur3sls(X = X, Y = Y, W = W, nG = nG, nR = nR, nT = nT,
                     p = p, type = type, maxlagW = maxlagW)
  }
  res
}




