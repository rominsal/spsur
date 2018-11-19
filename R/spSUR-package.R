#'spSUR
#'
#' @docType package
#' @name spSUR-package
#' @rdname spSUR-package
#'
#' @title Spatial Seemingly Unrelated Regression Models
#'
#' @description
#' A collection of functions to test and estimate Spatial Seemingly Unrelated Regression (SUR)
#' models by maximum likelihood and three-stage least square the usual
#' SUR especifications with spatial effects
#' (SUR-SLX; SUR-SAR; SUR-SEM; SUR-SDM; SUR-SDEM and SUR-SARAR) and non spatial SUR model (SUR-SIM)\cr
#'
#' @details
#' Three functionalities have the \link[pkg]{spSUR} package:
#'
#' @section 1. Testing for spatial effects:
#'
#' The function \code{\link{lmtestspsur}} give five statistis for testing spatial dependence in Semmingly Unrelated Regression
#' models (SUR). The five statistics are based on Lagrangian Multipliers principe (SUR-LM-).
#' The statistics are similar to the LM tests of spatial autocorrelation for a single equation
#' (see \code{\link[spdep]{lm.LMtests}} in package \link[pkg]{spdep}).
#' The tests can be applied to a data set with nG>1 equations for one cross sectional data set (nR>1) in the same
#' time period (nT=1) or in a panel data framework. Two alternatives in panel: nT>1 time periods for the same cross
#' section (nR>1) assuming temporal correlation in residuals or nG>1 equations in nT>1 time
#' periods assuming correlation only in residuals of nG equations.\cr
#'
#' @section 2. Estimation of the Spatial SUR models:
#' Two alternative estimation methods are implemented
#' \cr
#' @section 2.1. Maximum likelihood estimation:
#' Maximum likelihood estimation for differents spatial SUR models using \code{\link{spsurml}}
#' . The models are:\cr
#' SUR-SIM: SUR model without spatial effects\cr
#' \deqn{y_{gt} = y_{gt} + X_{gt} \beta_{g} + \epsilon_{gt}}
#' SUR-SLX: \cr
#' \deqn{y_{gt} = y_{gt} + X_{gt} \beta_{g} + WX_{gt} \theta_{g} + \epsilon_{gt}}
#' SUR-SAR: Spatial autorregresive model\cr
#' \deqn{y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} + \epsilon_{gt}}
#' SUR-SEM: Spatial error model\cr
#' \deqn{y_{gt} = X_{gt} \beta_{g} + u_{gt} ; u_{gt} = \rho_{g} Wu_{gt} + \epsilon_{gt}}
#' SUR-SDM: Spatial Durbin model\cr
#' \deqn{y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} + WX_{gt} \theta_{g} + \epsilon_{gt}}
#' SUR-SDEM: Spatial Durbin error model\cr
#' \deqn{y_{gt} = X_{gt} \beta_{g} + WX_{gt} \theta_{g} + u_{gt} ; u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt}}
#' SUR-SARAR: Spatial autoregressive model with spatial autoregressive error term\cr
#' \deqn{y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} + u_{gt} ; u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt}}
#' where \eqn{\y_{gt}}, u_{gt} and \eqn{\epsilon_{gt}} are (nRx1) vectors; \eqn{\X_{gt}} is a matrix of exogenous variables of
#' order (nRxsum(p)); \eqn{\lambda_g} and \eqn{\rho_g} are parametres of spatial dependence; W is a nRxnR matrix of spatial interactions.
#' \cr
#' @section 2.2 Three-Stage Least Squate estimation:
#' The SUR-SAR and SUR-SDM could be estimate by 3sls using \code{\link{spsur3sls}}\cr
#'
#' @section 3. Diagnostic tests:
#' Several diagnostic tests could be calculate to looking for the corect specification:\cr
#'
#' @section 3.1 Marginal test:
#' The Marginal Multiplier tests are used to test for no correlation in one part of the model allowing for spatial correlation in the other.
#' The LM(\eqn{\lambda}|\eqn{\rho}) is the test for sustantive spatial autocorrelation in a model with spatial autocorrelation
#' in error term (SUR-SEM; SUR-SDEM).
#' The LM(\eqn{\rho}|\eqn{\lambda}) is the test for spatial error correlation in a model
#' with sustanative spatial correlation (SUR-SAR; SUR-SDM).\cr
#'
#' @section 3.2 Coefficients stability:
#' The Walt tests of coefficent stability (or equity) could be calculate for beta coefficients \code{\link{wald_betas}} or
#' for spatial autocorrelation coefficients \code{\link{wald_deltas}}
#'
#' @section 4. Marginal effects:
#' The marginal effects \code{\link{impacts}} of spatial autoregressive models (SAR; SDM; SARAR) has been calculated
#' following the propose of LeSage and Pace (2009)
#'
#' @section 5. Additional funtionalities:
#' A Data Generating Process of a spatial SUR models is avalible using the function \code{\link{dgp_spSUR}}.
#'
#' @section Dataset:
#' This package comes with four different datasets: spc; NAT; Italian Unemployment\cr
#' The spc (Spatial Phillips-Curve) is a classical data set from Anselin (1988, p. 203)\cr
#' The NAT data set is come from
#' \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#' Homicides and selected socio-economic characteristics for continental U.S. counties.
#' Data for four decennial census years: 1960, 1970, 1980 and 1990.\cr
#' A data set 'Sar' obtain with the function dgp_spSUR to imput data as matrices (alternative to use the \code{\link[Formula]{Formula}} package)
#' @references
#'
#' Breusch T, Pagan A (1980). The Lagrange multiplier test and its applications to model specification in
#'  econometrics. \emph{Review of Economic Studies} 47: 239-254
#' \cr
#' \cr
#' LeSage, J., and Pace, R. K. (2009). \emph{Introduction to spatial econometrics}. Chapman and Hall/CRC.
#' \cr
#' \cr
#' Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model selection strategies in a SUR framework.
#' The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' \cr
#' \cr
#' López, F.A., Martínez-Ortiz, P.J., and Cegarra-Navarro, J.G. (2017). Spatial spillovers in public
#' expenditure on a municipal level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#' \cr
#' \cr
#' Mur, J., Lopez, F., and Herrera, M. (2010). Testing for spatial effects in seemingly
#' unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#' \cr
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@unizar.es} \cr
#'   }
#'
#' @importFrom Formula Formula model.part
#' @importFrom MASS ginv
#' @importFrom Matrix bdiag crossprod Diagonal Matrix solve t
#' @importFrom methods as
#' @importFrom minqa bobyqa
#' @importFrom numDeriv hessian
#' @importFrom sparseMVN rmvn.sparse
#' @importFrom spdep knearneigh knn2nb nb2mat
#' @importFrom stats cor cov optim pchisq pnorm pt qnorm rnorm runif
#' @importFrom stats coefficients fitted lm residuals printCoefmat
#' @importFrom stats model.frame model.matrix terms
#'
#' @exportPattern "^[[:alpha:]]+"
NULL
