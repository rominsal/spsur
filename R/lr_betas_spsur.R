#' @name lr_betas_spsur
#' @rdname lr_betas_spsur
#'
#' @title Likelihood ratio for testing homogeneity constraints on beta coefficients of the SUR equations.
#'
#' @description Function \code{\link{lr_betas_spsur}} obtains a Likelihood Ratio test, LR in what follows,
#'   with the purpose of testing if some of the \eqn{\beta} coefficients in the \emph{G} equations of the
#'   \emph{SUR} model are equal. This function has a straightforward application, especially when \eqn{G=1},
#'   to the case of testing for the existence of structural breaks in the \eqn{\beta} parameters.
#'
#'   The function can test for the homogeneity of only one coefficient, of a few of them
#'   or even the homogeneity of all the slope terms. The testing procedure implies, \emph{first},
#'   the estimation of both a constrained and a unconstrained model and, \emph{second}, the comparison
#'   of the log-likelihoods to compute the LR statistics.
#'
#' @param printmodels Logical value: prints constrained and unconstrained models. Default = \code{FALSE}.
#'   models. Default = \code{FALSE}.
#' @param time Time variable.
#' @param trace Logical value to show intermediate results during estimation. Default = \code{FALSE}.
#' @inheritParams spsurml
#'

#' @return
#'  The function estimates two variants of the same SUR model, namely, a restricted and a unrestricted
#'  version. The purpose, as indicated above, is testing homogeneity constraints between the
#'  \eqn{\beta} parameters of the different equations of the SUR.
#'  The output of \code{\link{lr_betas_spsur}} shows, in first place, the iteration sequence of the
#'  maximum-likelihood algorithm. Then appears the Likelihood  Ratio, LR, test.
#'  The output also includes the maximum-likelihood estimation of the two models.
#'
#'   \tabular{ll}{
#'   \code{statistic} \tab The Values of the LR test. \cr
#'   \code{p_val}    \tab The p-value of the LR test. \cr
#'   \code{df}    \tab Degrees of freedom. \cr
#'   \code{llik_unr}    \tab The Log-likelihood of the unrestricted model. \cr
#'   \code{llik_res} \tab The Log-likelihood of the restricted model. \cr
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @references
#'   \itemize{
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsurtime}}, \code{\link{wald_betas}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1)  ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## H0: equal beta for SMSA in both equations.
#' R <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
#' b <- matrix(0,ncol=1)
#' LR_SMSA <-  lr_betas_spsur(Form = Tformula, data = spc, W = Wspc,
#'                            type = "sim", R = R, b = b, trace = TRUE,
#'                            printmodels = TRUE)
#'
#' ################################################
#' ####### PANEL DATA (G>1; Tm>1)          ########
#' ################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#' \donttest{
#' ## It usually requires 1-2 minutes maximum
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## H0: equal beta for PS80 and PS90 in both equations.
#' R <- matrix(c(0,1,0,0,-1,0),nrow=1)
#' b <- matrix(0,ncol=1)
#' LR_PS <-  lr_betas_spsur(Form = Tformula, data = NCOVR, W = W,
#'                         type = 'slm', R = R, b = b, printmodels = FALSE)
#'
#'}
#' ################################################################
#' ######## PANEL DATA: TEMPORAL CORRELATIONS (nG=1; nT>1) ########
#' ################################################################
#' ## Example 3: with classical panel data set. Database is
#' ##            a spatio-temporal panel
#' \donttest{
#' ### Only execute if you have enough memory...
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' N <- nrow(NCOVR)
#' Tm <- 4
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' pHR <- c(NCOVR$HR60, NCOVR$HR70, NCOVR$HR80, NCOVR$HR90)
#' pPS <- c(NCOVR$PS60, NCOVR$PS70, NCOVR$PS80, NCOVR$PS90)
#' pUE <- c(NCOVR$UE60, NCOVR$UE70, NCOVR$UE80, NCOVR$UE90)
#' pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
#'                      HR = pHR, PS = pPS, UE = pUE)
#' rm(NCOVR,pHR,pPS,pUE,index_time,index_indiv)
#' form_pHR <- HR ~ PS + UE
#' # H0: equal PS beta coefficient in equations 1, 3, and 4
#' R <- matrix(0,nrow=2,ncol=12) # nrow = number of restrictions ; ncol = number of beta parameters
#' R[1,2] <- 1; R[1,8] <- -1 # PS beta coefficient in equations 1 equal to 3
#' R[2,2] <- 1; R[2,11] <- -1 # PS beta coefficient in equations 1 equal to 4
#' b <- matrix(0,nrow=2,ncol=1)
#' lr_partrate <-  lr_betas_spsur(Form = form_pHR, data = pNCOVR,
#'                                time = pNCOVR$time, W = W,
#'                                type = "slm", R = R, b = b, trace = TRUE,
#'                                printmodels = FALSE)
#' }
#' @export
lr_betas <- function(object, R, b) {
  modelu <- object 
  lliku <- logLik(modelu)
  cat("\n Fitting restricted model...\n")
  modelr <- update(modelu, R = R, b = b) 
  llikr <- logLik(modelr)
  statistic <- -2*(llikr-lliku)
  attr(statistic, "names") <- "Likelihood ratio"
  parameter <- abs(attr(llikr, "df") - attr(lliku, "df"))
  if (parameter < 1) 
    stop("non-positive degrees of freedom: no test possible")
  attr(parameter, "names") <- "df"
  p.value <- 1 - pchisq(abs(statistic), parameter)
  estimate <- c(llikr, lliku)
  attr(estimate, "names") <- c("Log likelihood of restricted model", 
                               "Log likelihood of unrestricted model")
  method <- "Likelihood ratio on beta parameters"
  data.name <- modelu$call[[3]]
  res <- list(statistic = statistic, parameter = parameter, 
              p.value = p.value, estimate = estimate, 
              method = method, data.name = data.name)
  class(res) <- "htest"
  res
}
