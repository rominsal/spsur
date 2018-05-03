#' @name lrtestspsur
#' @rdname lrtestspsur
#'
#' @title Likelihood ratio test
#'
#'
#' @description
#' The function report serveral LR tests
#'
#' @param    nT     : Number of time periods
#' @param    nG     : number of equations
#' @param    nR     : number of spatial observations
#' @param    Y       : Data vector R*Tx1 (firs space, second time)
#' @param    X       : Data matrix R*TxK (K is a vector ...)
#' @param    W       : A RxR spatial weight matrix.
#'
#' @return
#' The LR tests
#'
#' @details
#' The SUR models
#' Various tests of spatial autocorrelation are obtain based on the principle of the Lagrange Multiplier in a maximum-likelihood framework.
#'
#' W is a spatial weights matrix
#'
#' @references
#' López, F.A., Mur, J., & Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. The Annals of Regional Science, 53(1), 197-220.
#' @seealso
#' \code{\link{spsurml}}
#' 
#' @examples
#'
#' data(Sar)
#' nT <- 4 # Number of periods
#' nG <- 3 # Number equations
#' nR <- ncol(Ws)
#' lrtestspsur(W=Ws,X=XXsar,Y=Ysar,nG=nG,nR=nR,nT=nT)

#' ########################################################################
#' ##########     TEMPORAL CORRELATIONS (nG = 1 and nT > 1)     ###########
#' ########################################################################
#' data("unemp_it_short") 
#' data("W_italy")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' lr_time <- lrtestspsur(Form=form_un,data=unemp_it,time=unemp_it$year,W=W_italy)
#' 
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

  ## Lik modelo SAR
  model_sar <- fit_spsursar(nT=nT,nG=nG,nR=nR,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sar <- model_sar$llsur

  ## Lik modelo SEM
  model_sem <- fit_spsursem(nT=nT,nG=nG,nR=nR,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sem <- model_sem$llsur

  ## Lik modelo SARAR
  model_sarar <- fit_spsursarar(nT=nT,nG=nG,nR=nR,
                                Y=Y,X=X,W=W,trace=FALSE)
  llsur_sarar <- model_sarar$llsur

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

  cat("LR test SIM-SAR: \n")
  cat("statistic: ",round(lr_sim_sar,3),
      " p-value: ",round(pval_lr_sim_sar,3),"\n\n")
  cat("LR test SIM-SEM: \n")
  cat("statistic: ",round(lr_sim_sem,3),
      " p-value: ",round(pval_lr_sim_sem,3),"\n\n")
  cat("LR test SIM-SARAR: \n")
  cat("statistic: ",round(lr_sim_sarar,3),
      " p-value: ",round(pval_lr_sim_sarar,3),"\n\n")
  cat("LR test SAR-SARAR: \n")
  cat("statistic: ",round(lr_sar_sarar,3),
      " p-value: ",round(pval_lr_sar_sarar,3),"\n\n")
  cat("LR test SEM-SARAR: \n")
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



