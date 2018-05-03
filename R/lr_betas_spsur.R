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
#' examples
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' H0: equal  beta for SMSA in both equations.
#' R1 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
#' r1 <- matrix(0,ncol=1)
#' lr_smsa <-  lr_betas_spsur(Form=Tformula,data=spc,
#'                               W=Wspc,type="sar",R=R1,r=r1,
#'                               trace=TRUE,printmodels = TRUE)


#' ####### Example with spatio-temporal SUR. 
#' ####### Assumption: nR,nT > 1 and nG = 1. Database is a spatio-temporal panel
#' data("unemp_it_short") 
#' data("W_italy")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' H0: equal emprgrowth beta in equations 1, 3, and 4 
#' R <- matrix(0,nrow=2,ncol=30)
#' R[1,2] <- 1; R[1,14] <- -1
#' R[2,2] <- 1; R[2,20] <- -1
#' r <- matrix(0,nrow=2,ncol=1)
#' lr_partrate <-  lr_betas_spsur(Form=form_un,data=unemp_it,time=unemp_it$year,
#'                           W=W_italy,type="sar",R2=R,r=r2,
#'                          trace=TRUE,printmodels = TRUE)
#' lr_betas_spcSUR_sar <- lr_betas_spsur(Form=Tformula,data=spc,
#'                                  R=R1,r=r1,type="sar",W=Wspc,
#'                                  printmodels=TRUE)


lr_betas_spsur<- function(Form=NULL,data=NULL,R=NULL,r=NULL,W=NULL,time=NULL,
                        X=NULL,Y=NULL,
                        nG=NULL,nR=NULL,nT=NULL,p=NULL,
                        type="sim",
                        printmodels=FALSE,
                        cov=FALSE,trace=FALSE) {
  
  if (is.null(R) || is.null(r)) stop("R and r must be specified as arguments")
  if (printmodels) cov <- TRUE
  
  start_fit <- proc.time()[3]
  cat("\n Fitting unrestricted model ...")
  if (is.null(time)){
    spcSUR_unr <-spsurml(Form=Form,data=data,type=type,W=W,
                         cov=cov,trace=trace)
  } else {
    spcSUR_unr <-spsurtime(Form=Form,data=data,type=type,W=W,
                         time=time,method="ml",cov=cov,trace=trace)
  }
  ll_unr <- spcSUR_unr$llsur
  #betas_unr <-  spcSUR_unr$betas
  #holg <- R %*% matrix(betas_unr,ncol=1) - r  
  end_fit <- proc.time()[3]
  cat("\n Time to fit unrestricted model: ",end_fit-start_fit," seconds\n\n")
  
  start_fit <- proc.time()[3]
  cat("\n Fitting restricted model ...")
  if (is.null(time)){
    spcSUR_res <-spsurml(Form=Form,data=data,R=R,r=r,type=type,W=W,
                           cov=cov,trace=trace)
  } else {
    spcSUR_res <-spsurtime(Form=Form,data=data,R=R,r=r,type=type,W=W,
                           time=time,method="ml",cov=cov,trace=trace)
  }  
  ll_res <- spcSUR_res$llsur
  end_fit <- proc.time()[3]
  cat("\n Time to fit restricted model: ",end_fit-start_fit," seconds\n\n")
  
  lr_stat <- -2*(ll_res-ll_unr) 
  lr_df <- nrow(R)
  lr_pval <- pchisq(lr_stat,df=lr_df,lower.tail=FALSE)

  cat("\n LR-Test \n")
  cat("\n Log-likelihood unrestricted model: ",round(ll_unr,4))
  cat("\n Log-likelihood restricted model: ",round(ll_res,4))
  cat("\n LR statistic: ",round(lr_stat,3)," degrees of freedom: ",lr_df,
      " p-value: (",lr_pval,")")
  if (printmodels) {
    cat("\n\n UNRESTRICTED MODEL \n\n")
    print(summary(spcSUR_unr))
    cat("\n\n RESTRICTED MODEL \n\n")
    print(summary(spcSUR_res))
  }
  res <- list(statistic=lr_stat,
              p_val=lr_pval,
              df=lr_df,
              llik_unr=ll_unr,
              llik_res=ll_res,
              R=R,
              r=r,
              type=type)
}  
