#' @name lr_betas_spsur
#' @rdname lr_betas_spsur
#'
#' @title Likelihood ratio test on coefficient homogeneity across equations.
#'
#' @description This function estimate both unrestricted and restricted models
#'   for testing equality between beta parameters of different equations.
#'
#' @param printmodels Logical value to print unrestricted and restricted
#'   models. Default = \code{FALSE}.
#' @param time Variable including temporal periods.
#' @param trace Logical value to get intermediate results.
#'   Default = \code{FALSE}.
#' @inheritParams spsurml
#'
#' @return The LR tests
#'   \tabular{ll}{
#'   \code{statistic} \tab The Values of LR tests. \cr
#'   \code{p_val}    \tab The p-value of LR tests. \cr
#'   \code{df}    \tab The degrees of freedom. \cr
#'   \code{llik_unr}    \tab The Log-likelihood of unrestricted model. \cr
#'   \code{llik_res} \tab The Log-likelihood of restricted model. \cr
#'   }
#'
#' @references
#'   \itemize{
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'   }
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
#' ## H0: equal beta for SMSA in both equations.
#' R <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
#' r <- matrix(0,ncol=1)
#' LR_SMSA <-  lr_betas_spsur(Form = Tformula, data = spc, W = Wspc,
#'                            type = "sar", R = R, r = r, trace = TRUE,
#'                            printmodels = TRUE)
#'
#' #################################################
#' ######## PANEL DATA (nG>1; nT>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#'
#' data(NAT)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## H0: equal beta for PS80 and PS90 in both equations.
#' R <- matrix(c(0,1,0,0,-1,0),nrow=1)
#' r <- matrix(0,ncol=1)
#' LR_PS <-  lr_betas_spsur(Form = Tformula, data = NAT, W = W,
#'                          type = 'sar', R = R, r = r, printmodels = TRUE)
#'
#' ##################################################
#' #### PANEL DATA. TEMP. CORRELAT. (nG=1;nT>1) #####
#' ##################################################
#'
#' ## Example 3: with spatio-temporal SUR.
#' ## Assumption: nR,nT > 1 and nG = 1. Database is a spatio-temporal panel
#' data("unemp_it")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' # H0: equal emprgrowth beta in equations 1, 3, and 4
#' R <- matrix(0,nrow=2,ncol=30)
#' R[1,2] <- 1; R[1,14] <- -1
#' R[2,2] <- 1; R[2,20] <- -1
#' r <- matrix(0,nrow=2,ncol=1)
#' lr_partrate <-  lr_betas_spsur(Form = form_un, data = unemp_it,
#'                                time = unemp_it$year, W = W_italy,
#'                                type = "sar", R = R, r = r, trace = TRUE,
#'                                printmodels = FALSE)
#' @export
lr_betas_spsur <- function(Form = NULL, data = NULL, R = NULL, r = NULL,
                          W = NULL, time = NULL, X = NULL, Y = NULL,
                          nG = NULL, nR = NULL, nT = NULL, p = NULL,
                          type = "sim", printmodels = FALSE,
                          cov = FALSE, trace = FALSE) {

  if (is.null(R) || is.null(r)) stop("R and r must be specified as arguments")
  if (printmodels) cov <- TRUE

  start_fit <- proc.time()[3]
  cat("\n Fitting unrestricted model ... \n")
  if (is.null(time)){
    spcsur_unr <-spsurml(Form = Form, data = data, type = type, W = W,
                         cov = cov, control = list(tol = 1e-3,
                                                   maxit = 200,
                                                   trace = trace) )
  } else {
    spcsur_unr <-spsurtime(Form = Form, data =data, type = type, W = W,
                           time = time, method = "ml", cov = cov,
                           trace = trace)
  }
  ll_unr <- spcsur_unr$llsur
  # betas_unr <-  spcsur_unr$betas
  # holg <- R %*% matrix(betas_unr,ncol=1) - r
  end_fit <- proc.time()[3]
  cat("\n Time to fit unrestricted model: ",end_fit-start_fit," seconds\n")

  start_fit <- proc.time()[3]
  cat("\n Fitting restricted model ...")
  if (is.null(time)){
    spcsur_res <-spsurml(Form = Form, data = data, R = R, r = r,
                         type = type, W = W, cov = cov,
                         control = list(tol = 1e-3, maxit = 200,
                                        trace = trace) )
  } else {
    spcsur_res <-spsurtime(Form = Form, data = data, R = R, r = r,
                           type = type, W = W, time = time,
                           method = "ml", cov = cov, trace = trace)
  }
  ll_res <- spcsur_res$llsur
  end_fit <- proc.time()[3]
  cat("Time to fit restricted model: ",end_fit-start_fit," seconds\n")

  lr_stat <- -2*(ll_res-ll_unr)
  lr_df <- nrow(R)
  lr_pval <- pchisq(lr_stat,df=lr_df,lower.tail=FALSE)

  cat("\n LR-Test \n")
  cat("\n Log-likelihood unrestricted model: ",round(ll_unr,4))
  cat("\n Log-likelihood restricted model: ",round(ll_res,4))
  cat("\n LR statistic: ",round(lr_stat,3)," degrees of freedom: ",lr_df,
      " p-value: (",lr_pval,")")
  if (printmodels) {
    cat("\n\n UNRESTRICTED MODEL \n")
    print(summary(spcsur_unr))
    cat("\n\n RESTRICTED MODEL \n")
    print(summary(spcsur_res))
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
