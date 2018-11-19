#' @name spsurtime
#' @rdname spsurtime
#'
#' @title Estimation of spatio-temporal SUR model either using maximum
#'        likelihood (ml) or three stages least squares (3sls).
#'
#' @description
#' Estimation of spatio-temporal SUR-LAG; SUR-ERR; SUR-SLX; SUR-SDM & SUR-SARAR
#'  models. Assumption: nR, nT>1 and nG = 1. The correlation matrix of the noise
#'  models temporal covariances. Database is a spatio-temporal panel.
#' @description
#' Estimation of spatio-temporal SUR-LAG; SUR-ERR; SUR-SLX; SUR-SDM & SUR-SARAR models
#'
#' @param    Form   : A object create with \code{\link{formula}}. Only allows
#'                    the same dependent and independent variables
#'                    for every year.
#' @param    data   : An object of class data.frame.
#' @param    time   : A variable including the temporal index in data frame.
#' @param    W      : A nRxnR spatial weight matrix
#' @param    R      : Default=NULL. A matrix including coefficients of linear hypothesis on beta
#' @param    r      : Default=NULL. A vector including independent vector of linear hypothesis on beta
#' @param    type   : type of model "sim" "slx" "sar" "sem" "sdm" "sdem" "sarar"
#' @param    method : Method of estimation either "ml" or "3sls". Default = "ml".
#' @param    maxlagW: Maximum number of spatial lags of independent variables to get instruments in 3sls. Default = 2.
#' @param    cov    : get the covariance matrix. Default = TRUE.
#' @param    trace  : TRUE to see intermediate results. Default = TRUE.
#'
#' @details
#' The \code{\link{spsurtime}} estimate a spatial SUR model with multiple temporal equations
#' for the same cross-section. The estimation consider temporal correlations.
#'
#' @return
#' Results of ML estimation of a SUR model with spatial effects. A list with:
#'   \tabular{ll}{
#'   \code{type}   \tab the matched call. \cr
#'   \code{method}   \tab Maximum Likelihood (by now!). \cr
#'   \code{call}   \tab the model estimate. \cr
#'   \code{coefficientes} \tab Coefficients estimated by maximum likelihood. \cr
#'   \code{p.value}    \tab the p-value of coefficientes. \cr
#'   \code{Cov}    \tab estimated covariance matrix of residuals between equations. \cr
#'   \code{Corr}    \tab estimated correlation matrix of residuals between equations. \cr
#'   \code{residuals}    \tab residual of model. \cr
#'   \code{LogLik}    \tab the maximum value of likelihood function. \cr
#'   \code{BP}    \tab Breush-Pagan test of diagonaltiy of covariance matrix. \cr
#'   \code{LM}    \tab Marginal test of spatial autocorrelation in the residuals. \cr
#'   }
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
#'   Fernando López  \tab \email{fernando.lopez@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsur3sls}}
#' @examples
#'
#' ####################################
#' ######## PANEL DATA (nG=1; nT>1) ###
#' ####################################
#'
#' data("unemp_it_short")
#' data("W_italy")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#'
#' ## SUR-SAR by Maximum Likelihood
#' unempitml_sar <- spsurtime(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy,type="sar",method="ml")
#' summary(unempitml_sar)
#'
#' ## SUR-SEM by Maximum Likelihood
#' unempitml_sem <- spsurtime(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy,type="sem",method="ml")
#' summary(unempitml_sem)
#'
#' ## SUR-SAR by Three-Stages Least Squares Estimation
#' unempit3sls_sar <- spsurtime(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy,type="sar",method="3sls")
#' summary(unempit3sls_sar)
#'
#' ## SUR-SDDM by Three-Stages Least Squares Estimation
#' unempit3sls_sdm <- spsurtime(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy,type="sdm",method="3sls")
#' summary(unempit3sls_sdm)

#' ############################ LR tests about betas in spatio-temporal models
#' H0: equal emprgrowth beta in equations 1, 3, 4 and 5
#' R <- matrix(0,nrow=2,ncol=30)
#' R[1,2] <- 1; R[1,14] <- -1
#' R[2,2] <- 1; R[2,20] <- -1
#' r <- matrix(0,nrow=2,ncol=1)
#' lr_partrate <-  lr_betas_spsur(Form=form_un,data=unemp_it,time=unemp_it$year,
#'                           W=W_italy,type="sar",R=R,r=r,
#'                          trace=TRUE,printmodels = TRUE)
#' ############################ Wald tests about betas in spatio-temporal models
#' wald_betas(unempitml_sar , R = R , r = r)

#' ############################ Wald tests about spatial-parameters in spatio-temporal models
#' H0: equal lambdas in sar model for equations 1, 2, 3.
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
    res <- spsurml(X=X,Y=Y,W=W,nG=nG,nR=nR,nT=nT,p=p,R=R,r=r,
                   type=type,cov=cov,trace=trace)
  }
  if (method == "3sls"){
    res <- spsur3sls(X=X,Y=Y,W=W,nG=nG,nR=nR,nT=nT,p=p,
                     type=type,maxlagW=maxlagW,trace=trace)
  }
  res
}




