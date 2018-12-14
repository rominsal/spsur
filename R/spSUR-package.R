#' @docType package
#' @name spsur-package
#' @rdname spsur-package
#'
#' @title Spatial Seemingly Unrelated Regression Models.
#'
#' @description A collection of functions to test and estimate Spatial
#'   Seemingly Unrelated Regression (SUR) models by maximum likelihood
#'   and three-stage least square the usual SUR especifications with spatial
#'   effects (SUR-SLX; SUR-SAR; SUR-SEM; SUR-SDM; SUR-SDEM and SUR-SARAR)
#'   and non spatial SUR model (SUR-SIM).
#'
#' @details Some functionalities have been included in \pkg{spsur} package:
#'
#' @section 1. Testing for spatial effects:
#'   The function \code{\link{lmtestspsur}} provides five statistics
#'   for testing spatial dependence in Semmingly Unrelated Regression(SUR)
#'   models. The five statistics are based on Lagrangian Multipliers
#'   principles (SUR-LM-). The statistics are similar to the LM tests of
#'   spatial autocorrelation for a single equation (see
#'   \code{\link[spdep]{lm.LMtests}} in package \pkg{spdep}). \cr
#'   These tests can be applied to a data set with \emph{nG}>1 equations for
#'   one cross sectional data set (\emph{nR}>1) in the same time period
#'   (\emph{nT}=1) or in a panel data framework. Two alternatives in
#'   panel data are implemented: \emph{nT}>1 time periods for the same
#'   cross section (\emph{nR}>1) assuming temporal correlation in
#'   residuals or \emph{nG}>1 equations in \emph{nT}>1 time periods
#'   assuming correlation only in residuals of \emph{nG} equations.
#'
#' @section 2. Estimation of the Spatial SUR models:
#'   Two alternative estimation methods are implemented:
#'      \enumerate{
#'        \item Maximum likelihood estimation. \cr
#'          Using \code{\link{spsurml}} function it is possible to estimate
#'          any of the next spatial SUR models:
#'          \itemize{
#'            \item \emph{SUR-SIM}: SUR model without spatial effects
#'              \deqn{ y_{gt} = y_{gt} + X_{gt} \beta_{g} + \epsilon_{gt} }
#'            \item \emph{SUR-SLX}:
#'              \deqn{ y_{gt} = y_{gt} + X_{gt} \beta_{g} +
#'                     WX_{gt} \theta_{g} + \epsilon_{gt} }
#'            \item \emph{SUR-SAR}: Spatial autorregresive model
#'              \deqn{y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} +
#'                    \epsilon_{gt} }
#'            \item \emph{SUR-SEM}: Spatial error model
#'              \deqn{ y_{gt} = X_{gt} \beta_{g} + u_{gt} }
#'              \deqn{ u_{gt} = \rho_{g} Wu_{gt} + \epsilon_{gt} }
#'            \item \emph{SUR-SDM}: Spatial Durbin model
#'              \deqn{ y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} +
#'                     WX_{gt} \theta_{g} + \epsilon_{gt} }
#'            \item \emph{SUR-SDEM}: Spatial Durbin error model
#'              \deqn{ y_{gt} = X_{gt} \beta_{g} +
#'                     WX_{gt} \theta_{g} + u_{gt} }
#'              \deqn{ u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt} }
#'            \item \emph{SUR-SARAR}: Spatial autoregressive model with
#'              spatial autoregressive error term
#'              \deqn{ y_{gt} = \lambda_{g} Wy_{gt} +
#'                     X_{gt} \beta_{g} + u_{gt} }
#'              \deqn{ u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt} }
#'          }
#'          where \eqn{y_{gt}}, \eqn{u_{gt}} and \eqn{\epsilon_{gt}} are
#'          (\emph{nR}x1) vectors; \eqn{X_{gt}} is a matrix of exogenous
#'          variables of order (\emph{nR}x\code{sum(p)}); \eqn{\lambda_g}
#'          and \eqn{\rho_g} are parameters of spatial dependence; \emph{W} is
#'          an \emph{nR}x\emph{nR} matrix of spatial interactions.

#'         \item Three-Stage Least Squares estimation. \cr
#'           \emph{SUR-SAR} and \emph{SUR-SDM} models can be estimated
#'           by 3sls using \code{\link{spsur3sls}} function.
#'       }
#'
#' @section 3. Diagnostic tests:
#'   Several diagnostic tests could be calculate to looking for
#'   the correct specification:
#'  \enumerate{
#'    \item \emph{Marginal tests} \cr
#'      The Marginal Multiplier tests are used to test for no correlation
#'      in one part of the model allowing for spatial correlation in the other.
#'      The LM(\eqn{\lambda}|\eqn{\rho}) is the test for sustantive
#'      spatial autocorrelation in a model with spatial autocorrelation
#'      in error term (SUR-SEM; SUR-SDEM). \cr
#'      The LM(\eqn{\rho}|\eqn{\lambda}) is the test for spatial error
#'      correlation in a model with sustanative spatial correlation
#'      (SUR-SAR; SUR-SDM).
#'
#'    \item \emph{Coefficients stability} \cr
#'      The Wald tests of coefficent stability (that is equal coefficients
#'      for different equations) can be calculate for beta coefficients
#'      using function \code{\link{wald_betas}} or for spatial
#'      autocorrelation coefficients using function
#'      \code{\link{wald_deltas}}.
#'  }
#'
#' @section 4. Marginal effects:
#'   The marginal effects of spatial autoregressive models (SAR; SDM; SARAR)
#'   can be computed with function \code{\link{impacts}}. The computations
#'   have been implemented following the propose of LeSage and Pace (2009).
#'
#' @section 5. Additional functionalities:
#'   To carry on simulations, data generating process of spatial SUR models
#'   are available  using the function \code{\link{dgp_spsur}}.
#'
#' @section Datasets:
#'   This package comes with four different datasets: spc, NAT,
#'   Italian Unemployment and simulated data. \cr
#'   \itemize{
#'     \item The \emph{spc} dataset (Spatial Phillips-Curve) is a classical
#'       dataset from Anselin (1988, p. 203).
#'      \item The \emph{NAT} dataset comprises Homicides and selected
#'        socio-economic characteristics for continental U.S. counties in
#'        four decennial census years: 1960, 1970, 1980 and 1990.
#'        It is available from
#'        \url{https://geodacenter.github.io/data-and-lab/ncovr/}.
#'      \item \emph{Italian Unemployment} includes a panel data of
#'        unemployment in Italy in period 1996-2014 at province level.
#'      \item A data set \emph{Sar} obtained with the function
#'        \code{\link{dgp_spsur}} to input data as matrices. The other input
#'         alternative is to use the \code{\link[Formula]{Formula}} package and
#'         a data frame.
#'   }
#'
#' @references
#'   \itemize{
#'     \item Breusch T, Pagan A (1980). The Lagrange multiplier test
#'       and its applications to model specification in econometrics.
#'       \emph{Review of Economic Studies} 47: 239-254.
#'
#'      \item LeSage, J., and Pace, R. K. (2009). \emph{Introduction to
#'        spatial econometrics}. Chapman and Hall/CRC.
#'
#'      \item Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science},
#'        53(1), 197-220.
#'
#'      \item López, F.A., Martínez-Ortiz, P.J., and Cegarra-Navarro, J.G.
#'        (2017). Spatial spillovers in public expenditure on a municipal
#'        level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#'
#'      \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'        effects in seemingly unrelated regressions. \emph{Spatial Economic
#'        Analysis}, 5(4), 399-440.
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
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
NULL
