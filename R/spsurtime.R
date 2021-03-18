#' @name spsurtime
#' @rdname spsurtime
#'
#' @title Estimation of SUR models for simple spatial panels (G=1).
#'
#' @description This function estimates SUR models for simple spatial 
#'   panel datasets.
#'   \code{\link{spsurtime}} is restricted, specifically, to cases where 
#'   there is only one equation, \emph{G=1}, and a varying number of spatial 
#'   units, \emph{N}, and time periods, \emph{Tm}. The SUR structure appears
#'   in form of serial dependence among the error terms corresponding to the 
#'   same spatial unit.
#'   Note that it is assumed that all spatial units share a common pattern 
#'   of serial dependence.
#'
#'   The user can choose between different types of spatial specifications, 
#'   as described below, and the estimation algorithms allow for the 
#'   introduction of linear restrictions on the \eqn{\beta} parameters
#'   associated to the regressors. The spatial panels with SUR structure 
#'   can be estimated by maximum-likelihood methods or three-stages least 
#'   squares procedures, using spatial instrumental variables.
#'   
#' @usage spsurtime (formula, data, time, na.action,
#'                   listw = NULL, type = "sim", Durbin = NULL, 
#'                   method = "eigen", fit_method = "ml", maxlagW = NULL,
#'                   zero.policy = NULL, interval = NULL, trs = NULL,
#'                   R = NULL, b = NULL, demean = FALSE, control = list() )
#' @param time Time variable.
#' @param fit_method Method of estimation for the spatial panel SUR model, 
#'  either \emph{ml} or \emph{3sls}. Default = \emph{ml}.
#' @param demean  Logical value to allow for the demeaning of panel data. In this case,
#'   \code{\link{spsurml}} subtracts the individual mean to each spatial 
#'   or cross-sectional unit. Default = \code{FALSE}.
#' @param maxlagW Maximum spatial lag order of the regressors employed to 
#'   produce spatial instruments for the spatial lags of the explained 
#'   variables. Default = 2. Note that in case of \code{type} = "sdm", the 
#'   default value for \code{maxlagW} is set to 3 because the first lag of 
#'   the regressors, \eqn{WX_{tg}}, can not be used as spatial instruments.
#' @param trs Either \code{NULL} or vector of powered spatial weights 
#'       matrix traces output by \code{trW}. Default = \code{NULL}.
#' @inheritParams spsurml
#'
#' @details
#'  Function \code{\link{spsurtime}} only admits a formula, created with 
#'  \code{\link[Formula]{Formula}} and a dataset of class \code{data.frame} 
#'  or \code{matrix}. That is, the data cannot be uploaded using data 
#'  matrices \eqn{Y} and \eqn{X} provided for other functions in 
#'  this package. \cr
#'  The argument \code{time} selects the variable, in the \code{data.frame}, 
#'  associated to the time dimension in the panel dataset. Then 
#'  \code{\link{spsurtime}} operates as in Anselin (1988), that is,
#'  each cross-section is treated as if it were an equation in a SUR model, 
#'  which now has \emph{Tm} 'equations' and \emph{N} individuals. \cr
#'  The SUR structure appears because there is serial dependence in the errors 
#'  of each individual in the panel. The serial dependence in the errors is 
#'  not parameterized, but estimated non-parametrically in the \eqn{Sigma} 
#'  covariance matrix returned by the function. An important constraint to 
#'  mention is that the serial dependence assumed to be the same for all 
#'  individuals in the sample. Serial dependence among individuals is 
#'  excluded from Anselin approach.
#'
#'
#' @return An \code{spsur} object with the output of the maximum-likelihood or 
#' three-stages least-squares estimation of the spatial panel SUR model. 
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
#'     \item Anselin, L. (1988). Spatial econometrics: methods and models. Dordrecht,
#'      Kluwer Academic Publishers.
#'     \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'     \item López, F.A., Martínez-Ortiz, P.J., and Cegarra-Navarro, J.G. (2017).
#'        Spatial spillovers in public expenditure on a municipal level in
#'        Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 
#'       5(4), 399-440. \url{https://doi.org/10.1080/17421772.2010.516443}
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsur3sls}}, 
#' \code{\link{wald_betas}},  \code{\link{wald_deltas}},
#' \code{\link{lmtestspsur}}, \code{\link{lrtest}}
#'
#' @examples
#'
#' ####################################
#' ######## PANEL DATA (G=1; Tm>1) ###
#' ####################################
#'
#' ## Example 1:
#' rm(list = ls()) # Clean memory
#' data(spc)
#' lwspc <- spdep::mat2listw(Wspc, style = "W")
#' N <- nrow(spc)
#' Tm <- 2
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' WAGE <- c(spc$WAGE83, spc$WAGE81)
#' UN <- c(spc$UN83, spc$UN80)
#' NMR <- c(spc$NMR83, spc$NMR80)
#' SMSA <- c(spc$SMSA, spc$SMSA)
#' pspc <- data.frame(index_indiv, index_time, WAGE, UN,
#'                     NMR, SMSA)
#' form_pspc <- WAGE ~ UN + NMR + SMSA
#' form2_pspc <- WAGE | NMR ~ UN  | UN + SMSA
#'
#' # SLM 
#' pspc_slm <- spsurtime(formula = form_pspc, data = pspc, 
#'                       listw = lwspc,
#'                       time = pspc$index_time, 
#'                       type = "slm", fit_method = "ml")
#'summary(pspc_slm)
#' pspc_slm2 <- spsurtime(formula = form2_pspc, data = pspc, 
#'                       listw = lwspc,
#'                       time = pspc$index_time, 
#'                       type = "slm", fit_method = "ml")
#'summary(pspc_slm2)
#'
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_spsurtime, package="spsur")
#' 
#' \donttest{ 
#' ### Example 2:
#' rm(list = ls()) # Clean memory
#' ### Read NCOVR.sf object
#' data(NCOVR, package="spsur")
#' nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ### Some regions with no links...
#' lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#' N <- nrow(NCOVR.sf)
#' Tm <- 4
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' pHR <- c(NCOVR.sf$HR60, NCOVR.sf$HR70, NCOVR.sf$HR80, NCOVR.sf$HR90)
#' pPS <- c(NCOVR.sf$PS60, NCOVR.sf$PS70, NCOVR.sf$PS80, NCOVR.sf$PS90)
#' pUE <- c(NCOVR.sf$UE60, NCOVR.sf$UE70, NCOVR.sf$UE80, NCOVR.sf$UE90)
#' pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
#'                      HR = pHR, PS = pPS, UE = pUE)
#' form_pHR <- HR ~ PS + UE
#'
#' ## SIM
#'
#' pHR_sim <- spsurtime(formula = form_pHR, data = pNCOVR, 
#'                     time = pNCOVR$time, type = "sim", fit_method = "ml")
#' summary(pHR_sim)
#'
#' ## SLM by 3SLS. 
#'
#' pHR_slm <- spsurtime(formula = form_pHR, data = pNCOVR, listw = lwncovr,
#'                      time = pNCOVR$time, type = "slm", 
#'                      fit_method = "3sls")
#' summary(pHR_slm)
#'
#' ############################# Wald tests about betas in spatio-temporal models
#' ### H0: equal betas for PS in equations 1, 3 and 4.
#' R <- matrix(0, nrow = 2, ncol = 12) 
#' ## nrow = number of restrictions 
#' ## ncol = number of beta parameters
#' R[1, 2] <- 1; R[1, 8] <- -1 # PS beta coefficient in equations 1 equal to 3
#' R[2, 2] <- 1; R[2, 11] <- -1 # PS beta coefficient in equations 1 equal to 4
#' b <- matrix(0, nrow=2, ncol=1)
#' wald_betas(pHR_sim , R = R , b = b) # SIM model
#' wald_betas(pHR_slm , R = R , b = b) # SLM model

#' ############################# Wald tests about spatial-parameters in
#' ############################# spatio-temporal models
#' ### H0: equal rhos in slm model for equations 1 and 2.
#' R2 <- matrix(0, nrow = 1, ncol = 4)
#' R2[1, 1] <- 1; R2[1, 2] <- -1
#' b2 <- matrix(0, nrow = 1, ncol = 1)
#' wald_deltas(pHR_slm, R = R2, b = b2)
#' }
#' @export
spsurtime <- function(formula, data, time, na.action,
                      listw = NULL, type = "sim", 
                      Durbin = NULL, method = "eigen",
                      fit_method = "ml", maxlagW = NULL,
                      zero.policy = NULL, 
                      interval = NULL, 
                      trs = NULL,
                      R = NULL, b = NULL, demean = FALSE,
                      control = list() ) {
  if (!inherits(formula, "Formula")) 
    formula <- Formula::Formula(formula)
  if (!inherits(time, "factor")) time <- as.factor(time)
  time <- droplevels(time)
  if (length(time) != nrow(data)) 
    stop("time must have same length than the number of rows in data")  
  mt <- terms(formula)
  G <- length(levels(time))
  Ylist <- vector("list", G)
  Xlist <- vector("list", G)
  p <- NULL
  namesX <- NULL
  levels_time <- levels(time)
  for (i in 1:G) {
    data_i <- model.frame(mt, 
                data = data[time == levels_time[i],])
    na.act <- attr(data_i, "na.action")
    if (G > length(formula)[1]) {
      Ylist[[i]] <- data_i[, 1]
    } else {
      Ylist[[i]] <- Formula::model.part(formula,
                                      data = data_i, 
                                      lhs = i)
    }
    if (G > length(formula)[2]) {
      Xlist[[i]] <- model.matrix(mt, 
                                 data = data[time == levels_time[i],])
    } else {
      Xlist[[i]] <- model.matrix(formula, 
                                 data = data_i, 
                                 rhs = i)
    }
    p <- c(p, ncol(Xlist[[i]]))
    namesX <- c(namesX, 
                paste(colnames(Xlist[[i]]), 
                      i, sep = "_"))
  }
  Y <- matrix(unlist(Ylist), ncol=1)
  X <- as.matrix(Matrix::bdiag(Xlist))
  colnames(X) <- namesX
  if (inherits(Ylist[[1]], 
               "data.frame"))
    N <- nrow(Ylist[[1]])
  if (is.vector(Ylist[[1]])) 
    N <- length(Ylist[[1]]) 
  Tm <- 1
  if (demean) { # demeaning for pure panel data
    # First reorder X matrix
    X <- NULL
    for (i in 1:G) X <- rbind(X, Xlist[[i]])
    #Tm == G
    demeanXY <- demeaning(X = X, Y = Y, G = G, N = N, Tm = G,
                          p = p, pdemean = TRUE)
      X <- demeanXY$Xdem
      Y <- demeanXY$Ydem
      p <- demeanXY$pdem
  }
  if (fit_method == "ml") {
    res <- spsurml(X = X, Y = Y, listw = listw, method = method,
                   G = G, N = N, Tm = Tm, na.action = na.act,
                   p = p, R = R, b = b, type = type, 
                   control = control)
  }
  if (fit_method == "3sls") {
    res <- spsur3sls(X = X, Y = Y, listw = listw, na.action = na.act,
                     G = G, N = N, Tm = Tm,
                     p = p, type = type, maxlagW = maxlagW,
                     R = R, b = b)
  }
  res$call <- match.call()
  res
}




