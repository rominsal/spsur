#' @name lrtestspsur
#' @rdname lrtestspsur
#'
#' @title Likelihood ratio test
#'
#'
#' @description
#' The function report Likelihood ratios tests for nested SUR models
#'
#' @param    Form   : An object create with \code{\link[Formula]{Formula}} package allowing for multiple responses and multiple parts of regressors
#' @param    data   : An object of class data.frame or a matrix
#' @param    W      : A nRxnR spatial weight matrix
#' @param    Y      : Default NULL. Data vector nRxnTx1 (firs: space dimension | second: time periods)
#' @param    X      : Default NULL. Data matrix nRxnTxp (p=sum(\eqn{p_{g}}) where \eqn{p_{g}} is the number of independent variables for g-th equation, g=1,...,nG)
#' @param    nG     : Default NULL. number of equations
#' @param    nR     : Default NULL. number of spatial observations
#' @param    nT     : Default NULL. Number of time periods
#'
#' @return
#' The code \code{\link{lrtestspsur}} print the LR tests\cr
#' \cr
#' A list with the values of LR test and p-values
#'
#' @details
#' The code \code{\link{lrtestspsur}} obtain the LogLik of the models SUR-SIM; SUR-SAR; SUR-SEM and SUR-SARAR
#' and print the Likelihood ratio tests:\cr
#' \cr
#' SUR-SIM vs SUR-SAR\cr
#' SUR-SIM vs SUR-SEM\cr
#' SUR-SIM vs SUR-SARAR\cr
#' SUR-SAR vs SUR-SARAR\cr
#' SUR-SEM vs SUR-SARAR\cr
#'
#' @references
#' Mur, J., Lopez, F., and Herrera, M. (2010). Testing for spatial effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#' \cr
#' \cr
#' Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' \cr
#' \cr
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsur3sls}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (nG>1; nT=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' LRtests <- lrtestspsur(Form=Tformula,data=spc,W=Wspc)
#'
#' #################################################
#' ######## PANEL DATA (nG>1; nT>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' Homicides and selected socio-economic characteristics for continental U.S. counties.
#' Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#'
#' data(NAT)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' LRtests <- lrtestspsur(Form=Tformula,data=NAT,W=W)
#'
#' ##################################################
#' #### PANEL DATA. TEMP. CORRELAT. (nG=1;nT>1) #####
#' ##################################################
#'
#' #### Example 3: Italian unemployment
#' data("unemp_it_short")
#' data("W_italy")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' LR_time <- lrtestspsur(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy)
#' @export
lrtestspsur <- function(Form=NULL,data=NULL,W=NULL,
                         X=NULL,Y=NULL,time=NULL,
                         nG=NULL,nR=NULL,nT=NULL)
{
  ## PURPOSE:
  # Realiza todos los test LR.
  if (is.null(W) && !type=="sim") stop("W matrix is needed")
  if (!is.null(W)) W <- Matrix::Matrix(W)
  if (is.null(time)){  # nG > 1 (no temporal correlations are modelled)
    if(!is.null(Form) && !is.null(data)){
      # Lectura datos
      if (!class(Form) == "Formula") Form <- Formula::Formula(Form)
      get_XY <- get_data_spsur(formula=Form,data=data,W=W)
      Y <- get_XY$Y
      X <- get_XY$X
      nG <- get_XY$nG
      nR <- get_XY$nR
      nT <- get_XY$nT
      p <- get_XY$p
      rm(get_XY)
    }
  } else { #nG = 1 and nT > 1 (temporal correlations are modelled)
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
  }

  ## Lik modelo SIM
  model_sim <- fit_spsursim(nT=nT,nG=nG,nR=nR,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sim <- model_sim$llsur
  cat("LogLik SUR-SIM: ", round(llsur_sim,3),"\n")
  ## Lik modelo SAR
  model_sar <- fit_spsursar(nT=nT,nG=nG,nR=nR,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sar <- model_sar$llsur
  cat("LogLik SUR-SAR: ", round(llsur_sar,3),"\n")
  ## Lik modelo SEM
  model_sem <- fit_spsursem(nT=nT,nG=nG,nR=nR,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sem <- model_sem$llsur
  cat("LogLik SUR-SEM: ", round(llsur_sem,3),"\n")
  ## Lik modelo SARAR
  model_sarar <- fit_spsursarar(nT=nT,nG=nG,nR=nR,
                                Y=Y,X=X,W=W,trace=FALSE)
  llsur_sarar <- model_sarar$llsur
  cat("LogLik SUR-SARAR: ", round(llsur_sarar,3),"\n\n")
  lr_sim_sar <- -2*(llsur_sim - llsur_sar)
  pval_lr_sim_sar <- pchisq(lr_sim_sar,df=nG,
                            lower.tail=FALSE)
  lr_sim_sem <- -2*(llsur_sim - llsur_sem)
  pval_lr_sim_sem <- pchisq(lr_sim_sem,df=nG,
                            lower.tail=FALSE)

  lr_sim_sarar <- -2*(llsur_sim - llsur_sarar)
  pval_lr_sim_sarar <- pchisq(lr_sim_sarar,df=2*nG,
                            lower.tail=FALSE)

  lr_sar_sarar <- -2*(llsur_sar - llsur_sarar)
  pval_lr_sar_sarar <- pchisq(lr_sar_sarar,df=nG,
                              lower.tail=FALSE)
  lr_sem_sarar <- -2*(llsur_sem - llsur_sarar)
  pval_lr_sem_sarar <- pchisq(lr_sem_sarar,df=nG,
                              lower.tail=FALSE)

  cat("LR test SIM-SIM versus SUR-SAR: \n")
  cat("statistic: ",round(lr_sim_sar,3),
      " p-value: ",round(pval_lr_sim_sar,3),"\n\n")
  cat("LR test SIM-SIM versus SUR-SEM: \n")
  cat("statistic: ",round(lr_sim_sem,3),
      " p-value: ",round(pval_lr_sim_sem,3),"\n\n")
  cat("LR test SIM-SIM versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_sim_sarar,3),
      " p-value: ",round(pval_lr_sim_sarar,3),"\n\n")
  cat("LR test SAR-SAR versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_sar_sarar,3),
      " p-value: ",round(pval_lr_sar_sarar,3),"\n\n")
  cat("LR test SEM-SEM versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_sem_sarar,3),
      " p-value: ",round(pval_lr_sem_sarar,3),"\n\n")

  res <- list(lr_sim_sar = lr_sim_sar,
              pval_lr_sim_sar=pval_lr_sim_sar,
              lr_sim_sem = lr_sim_sem,
              pval_lr_sim_sem = pval_lr_sim_sem,
              lr_sim_sarar  = lr_sim_sarar,
              pval_lr_sim_sarar = pval_lr_sim_sarar,
              lr_sar_sarar = lr_sar_sarar,
              pval_lr_sar_sarar = pval_lr_sar_sarar,
              lr_sem_sarar = lr_sem_sarar,
              pval_lr_sem_sarar = pval_lr_sem_sarar
              )
 res
}



