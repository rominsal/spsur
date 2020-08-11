#' @name spsurml
#' @rdname spsurml
#' @title Maximum likelihood estimation of spatial SUR model.
#' @description This function estimates spatial SUR models using 
#'   maximum-likelihood methods.The number of equations, time periods 
#'   and cross-sectional units is not restricted.The user can choose 
#'   between different spatial specifications as described below.
#'   The estimation procedure allows for the introduction of linear 
#'   restrictions on the \eqn{\beta} parameters associated to the 
#'   regressors.
#' @usage spsurml(formula = NULL, data = NULL, na.action,
#'                listw = NULL, type = "sim", Durbin = NULL,
#'                method = "eigen", zero.policy = NULL, interval = NULL,
#'                trs = NULL, R = NULL, b = NULL, X = NULL, Y = NULL, 
#'                G = NULL, N = NULL, Tm = NULL,p = NULL, 
#'                control = list() )
#' @param formula An object type \code{\link[Formula]{Formula}} 
#'   similar to objects created with the package \pkg{Formula} 
#'   describing the equations to be estimated in the model. 
#'   This model may contain several responses (explained 
#'   variables) and a varying number of regressors in each equation.
#' @param data An object of class data.frame or a matrix.
#' @param na.action	A function (default \code{options("na.action")}),
#'  can also be \code{na.omit} or \code{na.exclude} with consequences 
#'   for residuals and fitted values. It may be necessary to set 
#'   \code{zero.policy} to \code{TRUE} because this subsetting may 
#'   create no-neighbour observations. 
#' @param listw A \code{listw} object created for example by 
#'   \code{\link[spdep]{nb2listw}} from \pkg{spatialreg} package; if 
#'   \code{\link[spdep]{nb2listw}} not given, set to 
#'   the same spatial weights as the \code{listw} argument. It can
#'   also be a spatial weighting matrix of order \emph{(NxN)} instead of
#'   a \code{listw} object. Default = \code{NULL}.
#' @param method Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package. 
#'   "eigen" (default) - the Jacobian is computed as the product of 
#'   (1 - rho*eigenvalue) using \code{\link[spatialreg]{eigenw}}, and 
#'   "spam" or "Matrix_J" for strictly symmetric weights lists of 
#'   styles "B" and "C", or made symmetric by similarity 
#'   (Ord, 1975, Appendix C) if possible for styles "W" and "S", 
#'   using code from the spam or Matrix packages to calculate the 
#'   determinant; "Matrix" and "spam_update" provide updating Cholesky 
#'   decomposition methods; "LU" provides an alternative sparse matrix 
#'   decomposition approach. In addition, there are "Chebyshev" and 
#'   Monte Carlo "MC" approximate log-determinant methods; 
#'   the Smirnov/Anselin (2009) trace approximation is available 
#'   as "moments". Three methods: "SE_classic", "SE_whichMin", 
#'   and "SE_interp" are provided experimentally, the first to 
#'   attempt to emulate the behaviour of Spatial Econometrics 
#'   toolbox ML fitting functions. All use grids of log determinant 
#'   values, and the latter two attempt to ameliorate some features 
#'   of "SE_classic".
#' @param interval Search interval for autoregressive parameter.
#'   Default = \code{NULL}.
#' @param trs Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package.
#'   Default \code{NULL}, if given, a vector of powered spatial weights 
#'   matrix traces output by \code{\link[spdep]{trW}}.
#' @param zero.policy Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package. 
#'   If \code{TRUE} assign zero to the lagged value of zones without 
#'   neighbours, if \code{FALSE} assign \code{NA} - causing 
#'   \code{spsurml()} to terminate with an error. Default = \code{NULL}.  
#' @param Durbin If a formula object and model is type "sdm", "sdem" 
#'   or "slx" the subset of explanatory variables to lag for each equation.        
#' @param Y A column vector of order \emph{(NTmGx1)}, with the 
#'   observations of the explained variables. The ordering of the data 
#'   must be (first) equation, (second) time dimension and (third) 
#'   cross-sectional/spatial units. The specification of \emph{Y} is 
#'   only necessary if not available a \code{\link[Formula]{Formula}}
#'   and a data frame. Default = \code{NULL}.
#' @param X A data matrix of order \emph{(NTmGxp)} with the observations
#'   of the regressors. The number of covariates in the SUR model is 
#'   \emph{p} = \eqn{sum(p_{g})} where \emph{\eqn{p_{g}}} is the number 
#'   of regressors (including the intercept) in the g-th equation, 
#'   \emph{g = 1,...,G}). The specification of "X" is only 
#'   necessary if not available a \code{\link[Formula]{Formula}} and a 
#'   data frame. Default = \code{NULL}.
#' @param p Number of regressors by equation, including the intercept. 
#'   \emph{p} can be a row vector of order \emph{(1xG)}, if the number 
#'    of regressors is not the same for all the equations, or a scalar, 
#'    if the \emph{G} equations have the same number of regressors. The 
#'   specification of \emph{p} is only necessary if not available a 
#'   \code{\link[Formula]{Formula}} and a data frame.
#' @param G Number of equations.
#' @param N Number of cross-section or spatial units
#' @param Tm Number of time periods.
#' @param type Type of spatial model specification: "sim",
#'   "slx", "slm", "sem", "sdm", 
#'   "sdem", "sarar" or "gnm". Default = "sim".
#' @param R A row vector of order \emph{(1xpr)} with  the set of 
#'   \emph{r} linear constraints on the \emph{beta} parameters. The 
#'   \emph{first} restriction appears in the first \emph{p} terms,
#'   the second restriction in the next \emph{p} terms and so on. 
#'   Default = \code{NULL}.
#' @param b A column vector of order \emph{(rx1)} with the values of 
#'   the linear restrictions on the \emph{beta} parameters. 
#'   Default = \code{NULL}.
#' @param control List of additional control arguments.
#' @details
#'  The list of (spatial) models that can be estimated with the \emph{spsurml} function are:
#'  \itemize{
#'     \item "sim": SUR model with no spatial effects
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item "slx": SUR model with spatial lags of the regressors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item "slm": SUR model with spatial lags of the explained variables
#'       \deqn{y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item "sem": SUR model with spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} Wu_{tg} + \epsilon_{tg} }
#'     \item "sdm": SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tt} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item "sdem": SUR model with spatial lags of the regressors and spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + WX_{tg} \theta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
#'     \item "sarar": SUR model with spatial lags of the explained variables and spatial
#'       errors
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
#'     \item "gnm": SUR model with spatial lags of the explained variables, 
#'       regressors and spatial errors
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + 
#'                         WX_{tg} \theta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
#'   }
#'   
#' @return Object of \code{spsur} class with the output of the  
#'  maximum-likelihood estimation of the specified spatial SUR model. 
#'  A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of model specified. \cr
#'     \code{method} \tab Value of \code{method} argument to compute the 
#'       Jacobian \cr
#'     \code{Durbin} \tab Value of \code{Durbin} argument. \cr     
#'     \code{coefficients} \tab Estimated coefficients for the regressors. \cr
#'     \code{deltas} \tab Estimated spatial coefficients. \cr
#'     \code{rest.se} \tab Estimated standard errors for the 
#'       estimates of \emph{beta}. \cr
#'     \code{deltas.se} \tab Estimated standard errors for the estimates of 
#'       the spatial coefficients (\code{deltas}). \cr
#'     \code{resvar} \tab Estimated covariance matrix for the estimates of 
#'       \emph{beta's} and spatial coefficients (\code{deltas}).\cr
#'     \code{LL} \tab Value of the likelihood function at the 
#'       maximum-likelihood estimates. \cr
#'     \code{R2} \tab Coefficient of determination for each equation, 
#'       obtained as the squared of the correlation coefficient between the 
#'       corresponding explained variable and its estimate. 
#'       \code{spsurml} also shows a \emph{global} coefficient of
#'        determination obtained, in the same manner, for the set of 
#'          the \emph{G} equations. \cr
#'     \code{Sigma} \tab Estimated covariance matrix for the residuals of 
#'       the \emph{G} equations. \cr
#'     \code{fdHess} \tab Logical value of \code{fdHess} argument when 
#'       computing numerical covariances.  \cr  
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom for the residuals. \cr
#'     \code{fitted.values} \tab Estimated values for the dependent 
#'       variables. \cr
#'     \code{BP} \tab Value of the Breusch-Pagan statistic to test the 
#'       null hypothesis of diagonality among the errors of the \emph{G} 
#'       equations. \cr
#'     \code{LMM} \tab Marginal Lagrange Multipliers, 
#'       LM(\eqn{\rho}|\eqn{\lambda}) and
#'       LM(\eqn{\lambda}|\eqn{\rho}), to test for omitted spatial effects 
#'       in the specification. \cr
#'     \code{G} \tab Number of equations. \cr
#'     \code{N} \tab Number of cross-sections or spatial units. \cr
#'     \code{Tm} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation (including intercepts). \cr
#'     \code{Y} \tab Vector \emph{Y} of the explained variables of the 
#'       SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the regressors of the SUR model. \cr
#'     \code{W} \tab Spatial weighting matrix. \cr
#'     \code{zero.policy} \tab Logical value of \code{zero.policy} . \cr
#'     \code{interval} \tab	Search interval for spatial parameter. \cr
#'     \code{listw_style} \tab	Style of neighborhood matrix \code{W}. \cr
#'     \code{trs} \tab Either \code{NULL} or vector of powered spatial weights 
#'       matrix traces output by \code{trW}. \cr
#'     \code{insert} \tab Logical value to check if \code{is.null(trs)}. \cr 
#'  }
#'
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{tol} \tab Numerical value for the tolerance for the estimation 
#'       algorithm until convergence. Default = 1e-3. \cr
#'     \code{maxit} \tab Maximum number of iterations until convergence; 
#'       it must be an integer value. Default = 200. \cr
#'     \code{trace} \tab A logical value to show intermediate results during 
#'       the estimation process. Default = \code{TRUE}. \cr
#'     \code{fdHess} \tab Compute variance-covariance matrix using the numerical 
#'       hessian. Suited for large samples. Default = \code{FALSE} \cr
#'     \code{Imult} \tab default 2; used for preparing the Cholesky 
#'       decompositions for updating in the Jacobian function \cr
#'     \code{super} \tab  if \code{NULL} (default), set to \code{FALSE} to use 
#'       a simplicial decomposition for the sparse Cholesky decomposition and 
#'       method "Matrix_J", set to as.logical(NA) for method "Matrix", if 
#'       \code{TRUE}, use a supernodal decomposition \cr
#'     \code{cheb_q} \tab default 5; highest power of the approximating 
#'       polynomial for the Chebyshev approximation \cr
#'     \code{MC_p} \tab default 16; number of random variates \cr
#'     \code{MC_m} \tab default 30; number of products of random variates 
#'       matrix and spatial weights matrix \cr
#'     \code{spamPivot} \tab  default "MMD", alternative "RCM" \cr
#'     \code{in_coef} \tab default 0.1, coefficient value for initial Cholesky 
#'       decomposition in "spam_update" \cr
#'     \code{type} \tab default "MC", used with method "moments"; alternatives 
#'       "mult" and "moments", for use if trs is missing \cr 
#'     \code{correct} \tab default \code{TRUE}, used with method "moments" to 
#'       compute the Smirnov/Anselin correction term \cr
#'     \code{trunc} \tab default \code{TRUE}, used with method "moments" to 
#'       truncate the Smirnov/Anselin correction term \cr
#'     \code{SE_method} \tab default "LU", may be "MC" \cr
#'     \code{nrho} \tab default 200, as in SE toolbox; the size of the first 
#'       stage lndet grid; it may be reduced to for example 40 \cr
#'     \code{interpn} \tab default 2000, as in SE toolbox; the size of the 
#'       second stage lndet grid \cr
#'     \code{SElndet} \tab default \code{NULL}, may be used to pass a 
#'       pre-computed SE toolbox style matrix of coefficients and their lndet 
#'       values to the "SE_classic" and "SE_whichMin" methods \cr
#'     \code{LU_order} \tab default \code{FALSE}; used in "LU_prepermutate", 
#'       note warnings given for lu method \cr
#'     \code{pre_eig} \tab default \code{NULL}; may be used to pass a 
#'       pre-computed vector of eigenvalues \cr
#' }
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
#'     \item Anselin, L. (1988). \emph{Spatial econometrics: methods and models.} 
#'       Dordrecht: Kluwer
#'     \item Bivand, R.S. and  Piras G. (2015). Comparing Implementations of 
#'      Estimation Methods for Spatial Econometrics. \emph{Journal of 
#'      Statistical Software}, 63(18), 1-36. 
#'      https://www.jstatsoft.org/v63/i18/.
#'     \item Bivand, R. S., Hauke, J., and Kossowski, T. (2013). 
#'       Computing the Jacobian in Gaussian spatial autoregressive models: An 
#'       illustrated comparison of available methods. \emph{ Geographical 
#'       Analysis}, 45(2), 150-179.
#'     \item Breusch T, Pagan A (1980). The Lagrange multiplier test and its
#'        applications to model specification in econometrics. 
#'        \emph{Rev Econ Stud} 47: 239-254
#'     \item Cliff, A. D. and Ord, J. K. (1981). \emph{Spatial processes: Models 
#'       and applications}, Pion. 
#'     \item LeSage J and Pace, R.K. (2009). \emph{Introduction to Spatial 
#'       Econometrics.} CRC Press, Boca Raton.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'     \item Ord, J. K. (1975). Estimation methods for models of spatial 
#'       interaction, \emph{Journal of the American Statistical Association}, 
#'       70, 120-126; 
#'   }
#'
#' @seealso
#' \code{\link{spsur3sls}}, \code{\link[spatialreg]{lagsarlm}}, 
#' \code{\link{lmtestspsur}}, \code{\link{wald_betas}}, 
#' \code{\link{lrtest}}
#'
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' spcsur.sim <-spsurml(formula = Tformula, data = spc, type = "sim")
#' summary(spcsur.sim)
#' # All the coefficients in a single table.
#' print(spcsur.sim)
#' # Plot of the coefficients of each equation in different graphs
#' plot(spcsur.sim) 
#'
#' ## A SUR-SLX model 
#' ## (listw argument can be either a matrix or a listw object )
#' spcsur.slx <-spsurml(formula = Tformula, data = spc, type = "slx", 
#'   listw = Wspc)
#' summary(spcsur.slx)
#' # All the coefficients in a single table.
#' print(spcsur.slx)
#' # Plot of the coefficients in a single graph
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.slx, viewplot = FALSE)
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                         nrow = 2)
#' } 
#' 
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_spsurml, package="spsur")
#'   
#' \donttest{
#' ### A SUR-SLM model
#' spcsur.slm <-spsurml(formula = Tformula, data = spc, type = "slm", 
#'                      listw = Wspc)
#' summary(spcsur.slm)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.slm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' #' ### A SUR-SEM model
#' spcsur.sem <-spsurml(formula = Tformula, data = spc, type = "sem", 
#'                      listw = Wspc)
#' summary(spcsur.sem)                      
#' print(spcsur.sem)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.sem, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ### A SUR-SDM model
#' spcsur.sdm <-spsurml(formula = Tformula, data = spc, type = "sdm", 
#'                      listw = Wspc)
#' summary(spcsur.sdm)
#' print(spcsur.sdm)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.sdm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ## A SUR-SDM model with different spatial lags in each equation
#' TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
#' spcsur.sdm2 <-spsurml(formula = Tformula, data = spc, type = "sdm", 
#'                       listw = Wspc, Durbin = TformulaD)
#' summary(spcsur.sdm2)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.sdm2, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' ### A SUR-SDEM model
#' spcsur.sdem <-spsurml(formula = Tformula, data = spc, type = "sdem", 
#'                       listw = Wspc)
#' print(spcsur.sdem)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.sdem, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ### A SUR-SARAR model
#' spcsur.sarar <-spsurml(formula = Tformula, data = spc, type = "sarar", 
#'                        listw = Wspc, control = list(tol = 0.1))
#' print(spcsur.sarar)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.sarar, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ### A SUR-GNM model
#' spcsur.gnm <-spsurml(formula = Tformula, data = spc, type = "gnm", 
#'                        listw = Wspc, control = list(tol = 0.1))
#' print(spcsur.gnm)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.gnm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ## A A SUR-GNM model model with different spatial lags in each equation
#' TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
#' spcsur.gnm2 <-spsurml(formula = Tformula, data = spc, type = "gnm", 
#'                        listw = Wspc, Durbin = TformulaD,
#'                        control = list(tol = 0.1))
#' print(spcsur.gnm2)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur.gnm2, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#'
#' }
#' 
#' ##################################################
#' #########  G=1; Tm>1         ########
#' ##################################################
#' #
#' ##### Example 2: Homicides + Socio-Economics (1960-90)
#' ## Homicides and selected socio-economic characteristics for continental
#' ## U.S. counties.
#' ## Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' ## \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#' 
#'\donttest{
#' ### It usually requires 1-2 minutes maximum...
#' rm(list = ls()) # Clean memory
#' ### Read NCOVR.sf object
#' data(NCOVR, package = "spsur")
#' nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ### Some regions with no links...
#' lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ### A SUR-SIM model
#' NCOVRSUR.sim <-spsurml(formula = Tformula, data = NCOVR.sf, type = "sim")
#' summary(NCOVRSUR.sim)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.sim, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], nrow = 3)
#' }
#' ### A SUR-SLX model
#' NCOVRSUR.slx <-spsurml(formula = Tformula, data = NCOVR.sf, type = "slx", 
#'                        listw = lwncovr, zero.policy = TRUE)
#' print(NCOVRSUR.slx)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.slx, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], nrow = 2)
#' }
#' 
#' ### A SUR-SLM model
#' ### method = "Matrix" (Cholesky) instead of "eigen"
#' ### (fdHess = TRUE to compute numerical covariances )
#' NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                        type = "slm", listw = lwncovr, method = "Matrix", 
#'                        zero.policy = TRUE, control = list(fdHess = TRUE))
#' summary(NCOVRSUR.slm)
#' 
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.slm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' # LR test for nested models
#' #lrtest(NCOVRSUR.sim, NCOVRSUR.slm)
#' 
#' ### A SUR-SDM model with different spatial lags in each equation
#' ### Analytical covariances (default)
#' TformulaD <- ~ PS80 + UE80 | PS90 
#' NCOVRSUR.sdm <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                        type = "sdm", listw = lwncovr, method = "Matrix",
#'                        Durbin = TformulaD, zero.policy = TRUE)
#' print(NCOVRSUR.sdm)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.sdm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' ### A SUR-SEM model
#' NCOVRSUR.sem <- spsurml(formula = Tformula, data = NCOVR.sf, 
#'                         type = "sem", listw = lwncovr, method = "Matrix",
#'                         zero.policy = TRUE, control = list(fdHess = TRUE))
#' print(NCOVRSUR.sem)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.sem, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ### A SUR-SDEM model
#' NCOVRSUR.sdem <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                         type = "sdem", listw = lwncovr, method = "Matrix",
#'                         zero.policy = TRUE, control = list(fdHess = TRUE))
#' print(NCOVRSUR.sdem)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR.sdem, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#'}
#' 
#' ###############################################
#' ######### SUR with G>1; Tm>1  ###
#' ###############################################
#' ##### Reshape NCOVR in panel format
#' \donttest{
#' N <- nrow(NCOVR.sf)
#' Tm <- 4
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' pHR <- c(NCOVR.sf$HR60, NCOVR.sf$HR70, NCOVR.sf$HR80, NCOVR.sf$HR90)
#' pPS <- c(NCOVR.sf$PS60, NCOVR.sf$PS70, NCOVR.sf$PS80, NCOVR.sf$PS90)
#' pUE <- c(NCOVR.sf$UE60, NCOVR.sf$UE70, NCOVR.sf$UE80, NCOVR.sf$UE90)
#' pDV <- c(NCOVR.sf$DV60, NCOVR.sf$DV70, NCOVR.sf$DV80, NCOVR.sf$DV90)
#' pFP <- c(NCOVR.sf$FP59, NCOVR.sf$FP70, NCOVR.sf$FP80, NCOVR.sf$FP90)
#' pSOUTH <- rep(NCOVR.sf$SOUTH, Tm)
#' pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
#'                      HR = pHR, PS = pPS, UE = pUE, DV = pDV,
#'                      FP = pFP, SOUTH = pSOUTH)
#' pform <- HR | DV | FP ~ PS + UE | PS + UE + SOUTH | PS
#' ### SIM 
#' ### Remark: It is necessary to provide Tm value as argument 
#' ### when G>1 && Tm>1
#' pNCOVRSUR.sim <- spsurml(formula = pform, data = pNCOVR, 
#'                          type = "sim", Tm = Tm)
#' print(pNCOVRSUR.sim)
#' if (require(gridExtra)) {
#'   pl <- plot(pNCOVRSUR.sim, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$lplbetas[[3]], nrow = 3)
#' }
#' # SLM
#' pNCOVRSUR.slm <- spsurml(formula = pform, data = pNCOVR, 
#'                          listw = lwncovr, type = "slm", method = "Matrix", Tm = Tm,
#'                          zero.policy = TRUE, control= list(fdHess = TRUE))
#' print(pNCOVRSUR.slm)
#' if (require(gridExtra)) {
#'   pl <- plot(pNCOVRSUR.slm, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$lplbetas[[3]], pl$pldeltas, nrow = 4)
#' }
#' 
#' pNCOVRSUR.sem <- spsurml(formula = pform, data = pNCOVR, 
#'                          listw = lwncovr, type = "sem", method = "Matrix", Tm = Tm,
#'                          zero.policy = TRUE, control= list(fdHess = TRUE))
#' print(pNCOVRSUR.sem) 
#' if (require(gridExtra)) {
#'   pl <- plot(pNCOVRSUR.sem, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$lplbetas[[3]], pl$pldeltas, nrow = 4)
#' }
#' }
#' @export
spsurml <- function(formula = NULL, data = NULL, na.action, 
                    listw = NULL, type = "sim", Durbin = NULL, 
                    method = "eigen",
                    zero.policy = NULL, 
                    interval = NULL, 
                    trs = NULL,
                    R = NULL, b = NULL, 
                    X = NULL, Y = NULL, 
                    G = NULL, N = NULL, Tm = NULL,
                    p = NULL, control = list() ) {
  con <- list(tol = 0.001, maxit = 200, trace = TRUE, 
              fdHess = NULL,
              Imult = 2, cheb_q = 5, MC_p = 16L, MC_m = 30L, super = NULL,
              spamPivot = "MMD", in_coef = 0.1, type = "MC", correct = TRUE,
              trunc = TRUE, SE_method = "LU", nrho = 200, interpn = 2000,
              SElndet = NULL, LU_order = FALSE,
              pre_eig = NULL) 
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
   if (!(type == "sim")) {
    if (is.null(listw) || !inherits(listw,c("listw","Matrix","matrix")))
      stop("listw format unknown or NULL")
    if (inherits(listw, "listw")) {
      if (is.null(formula) || is.null(data)) {
        W <- Matrix::Matrix(spdep::listw2mat(listw))
      }  
    }
    if (inherits(listw, "matrix")) {
      W <- Matrix::Matrix(listw)
      listw <- spdep::mat2listw(W)
    }  
    if (inherits(listw, "Matrix")) {
      W <- listw
      listw <- spdep::mat2listw(as.matrix(W))
    } 
  } else W <- NULL  
  
  if (is.null(zero.policy))
    zero.policy <- spatialreg::get.ZeroPolicyOption()
  can.sim <- FALSE
  if (!(is.null(listw)) && listw$style %in% c("W", "S")) {
    can.sim <- spatialreg::can.be.simmed(listw)
  } 
  
  if (!is.null(Tm) && !is.null(G) && Tm > 1 && G == 1){
    # Change dimensions (assumption: matrix as data)
    G <- Tm
    Tm <- 1
  }
  
  if (!is.null(formula) && (!inherits(formula, "Formula"))) 
    formula <- Formula::Formula(formula)
  cl <- match.call()
  if (!is.null(formula) && !is.null(data)) {
    mt <- terms(formula, data = data)
    mf <- lm(formula, data = data, na.action = na.action, 
             method = "model.frame")
    mf$drop.unused.levels <- TRUE
    na.act <- attr(mf, "na.action")
    if (!(type == "sim")) {
      if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy = zero.policy)
      }
      W <- Matrix::Matrix(spdep::listw2mat(listw))
    }
    if (any(type == c("gnm", "sdm", "sdem", "slx"))) {
      if(!inherits(Durbin, "formula")) Durbin <- TRUE
    } else { Durbin <- FALSE }
    get_XY <- get_data_spsur(formula = formula, mf = mf, 
                             Durbin = Durbin,
                             listw = listw, 
                             zero.policy = zero.policy, 
                             N = N, Tm = Tm)
    Y <- get_XY$Y
    X <- get_XY$X
    G <- get_XY$G
    N <- get_XY$N
    Tm <- get_XY$Tm
    p <- get_XY$p
    dvars <- get_XY$dvars
    if (Tm > 1 && G == 1) {
       # Change dimensions in this case with Matrix Data
       G <- Tm
       Tm <- 1
     }
    rm(get_XY)
    if (length(p) == 1) p <- rep(p,G)
  } else {# Input data in matrix form...
    dvars <- vector("list", G)
    for (i in 1:G) {
      dvars[[i]] <- c(p[i], 0L)
    } 
  }
  if (length(p) == 1) p <- rep(p,G)
  names(p) <- NULL
  if (!is.null(R) && !is.null(b)) {
    Xorig <- X
    porig <- p
    restr <- X_restr(X = X, R = R, b = b, p = p)
    X <- restr$Xstar
    p <- restr$pstar
  }

  #### ASIGNACIONES DE DATOS A NEW ENVIRONMENT ##############
  similar <- FALSE
  env <- new.env()
  assign("Y", Y, envir = env)
  assign("X", X, envir = env)
  assign("N", N, envir = env)
  assign("G", G, envir = env)
  assign("Tm", Tm, envir = env)
  assign("p", p, envir = env)
  assign("dvars", dvars, envir = env)
  # CÓDIGO EJEMPLO PARA DETERMINANTE JACOBIANO
  if (!(is.null(listw))) {
    assign("listw", listw, envir = env)
    assign("n", length(listw$neighbours), envir = env)
    assign("similar", FALSE, envir = env)
    assign("can.sim", can.sim, envir = env)
    assign("verbose", con$trace, envir = env)
    assign("family", "SAR", envir = env) # CHEQUEAR OTRAS OPCIONES
    if (!(type == "sim" || type == "slx")) {
      interval <- spatialreg::jacobianSetup(method, env, con, 
                                            pre_eig = con$pre_eig,
                                            trs = trs, 
                                            interval = interval)
      assign("interval", interval, envir = env)
      
    }
  }
  if (any(type == c("sim","slm","sem","sarar")))
    name_fit <- paste("fit_spsur", type, sep = "")
  if (type == "sdm") name_fit <- "fit_spsurslm"
  if (type == "sdem") name_fit <- "fit_spsursem"
  if (type == "slx") name_fit <- "fit_spsursim"
  if (type == "gnm") name_fit <- "fit_spsursarar"
  fit <- get(name_fit)
  if (con$trace)  start_fit <- proc.time()[3]
  # Maximize concentrate likelihood
  z <- fit(env = env, con = con)
  if (con$trace) {
    end_fit <- proc.time()[3]
    cat("Time to fit the model: ",
        end_fit-start_fit," seconds\n")
  }
  coefficients <- z$coefficients
  deltas <- z$deltas
  Sigma <- z$Sigma
  names_sigma <- NULL
  for (i in 1:G){
    new_name <- paste0("sigma",i,sep="")
    names_sigma <- c(names_sigma,new_name)
  }
  colnames(Sigma) <- rownames(Sigma) <- names_sigma
  LL <- z$LL
  parameters <- length(coefficients) + 
                  length(deltas) + G*(G + 1)/2
  df.residual <- G*N*Tm - parameters
  dn <- colnames(X)
  if (is.null(dn)) dn <- paste0("x", 1L:(G*sum(p)), sep = "")
  names(coefficients) <- dn
  names_deltas <- NULL
  for (i in 1:G) {
    if (any(type == c("slm","sdm")))
      names_deltas[i] <- paste("rho", i, sep = "_")
    if (any(type == c("sem","sdem")))
      names_deltas[i] <- paste("lambda", i, sep = "_")
    if (any(type == c("sarar", "gnm"))) {
      names_deltas[i] <- paste("rho", i, sep = "_")
      names_deltas[G + i] <- paste("lambda", i, sep = "_")
    }
  }
  names(deltas) <- names_deltas
  assign("Sigma", Matrix::Matrix(z$Sigma), envir = env)
  if (!(type == "sim" || type == "slx")) {
    assign("deltas",Matrix::Diagonal(length(deltas),deltas),
           envir = env)
  }
  if (con$trace) start_cov <- proc.time()[3]
  fdHess <- con$fdHess
  if (is.null(fdHess) || !(fdHess) || 
      any(type == c("sim","slx")))  {
    # ANALYTICAL VARIANCE-COVARIANCE MATRIX
    if (any(type == c("sim","slx")))
      name_cov_fit <- "cov_spsursim_f"
    if (any(type == c("slm","sem","sarar")))
      name_cov_fit <- paste("cov_spsur", type, "_f", sep = "")
    if (type == "sdm") name_cov_fit <- "cov_spsurslm_f"
    if (type == "sdem") name_cov_fit <- "cov_spsursem_f"
    if (type == "gnm") name_cov_fit <- "cov_spsursarar_f"
    cov_fit <- get(name_cov_fit)
    allcov <- try( cov_fit(env = env) )
    if (inherits(allcov, "try-error")) {
      cat("Impossible to compute analytical covariances ","\n")
      fdHess <- TRUE
    } else {
      fdHess <- FALSE
      rest.se <- allcov$rest.se
      names(rest.se) <- names(coefficients)
      deltas.se <- allcov$deltas.se
      names(deltas.se) <- names(deltas)
      if (!is.null(allcov$LMM)) LMM <- allcov$LMM else LMM <- NULL
      if (!is.null(allcov$BP)) BP <- allcov$BP else BP <- NULL
      if (any(type == c("sim","slx"))) {
        resvar <- allcov$vcov
        colnames(resvar) <- rownames(resvar) <- 
          names(coefficients)
      } else {
        resvar <- allcov$vcov
        names_var_Sigma <- NULL
        for (k in 1:G) {
          for (l in k:G) {
            new_name <- paste0("sigma",k,l,sep="")
            names_var_Sigma <- c(names_var_Sigma,new_name)
          }
        }
        colnames(resvar) <- rownames(resvar) <-
                c(names(coefficients),
                  names(deltas),names_var_Sigma)
        ## VIP: CAMBIO ORDEN MATRIZ COVARIANZAS
        ## IGUAL ORDEN QUE SPDEP Y SPATIALREG PACKAGES...
        resvar <- resvar[c(names_var_Sigma,names(deltas),names(coefficients)),
                         c(names_var_Sigma,names(deltas),names(coefficients))]
      }
    }
  }
  if (fdHess) {
    if (con$trace){
      cat("Computing numerical covariances...","\n")
    }
    if (any(type == c("slm","sdm"))) name_cov_fit <- "f_sur_lag"
    if (any(type == c("sem","sdem"))) name_cov_fit <- "f_sur_sem"
    if (any(type == c("sarar","gnm"))) name_cov_fit <- "f_sur_sarar"
    cov_fit <- get(name_cov_fit)
    vardeltas <- solve(numDeriv::hessian(func = cov_fit, 
                                         x = deltas, env = env))
    deltas.se <- as.vector(sqrt(diag(vardeltas)))
    names(deltas.se) <- names(deltas)
    IT <- Matrix::Diagonal(Tm)
    IR <- Matrix::Diagonal(N)
    Sigmainv <- Matrix::solve(z$Sigma)
    OMEinv <- kronecker(kronecker(IT, Sigmainv), IR)
    varbetas <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X))
    rest.se <- sqrt(diag(as.matrix(varbetas)))
    names(rest.se) <- names(coefficients)
    rm(IT,IR,OMEinv)
    resvar <- as.matrix(Matrix::bdiag(list(vardeltas, varbetas)))
    colnames(resvar) <- c(names(deltas), names(coefficients))
    rownames(resvar) <- colnames(resvar)
    
    ## Añadido por Fernando 16/03/2020: 
    # Incluyo BP
    Sigma_corr <- stats::cov2cor(Sigma)
    index_ltri <- lower.tri(Sigma_corr)
    BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
    # Para calcular los marginales es necesario la matriz de varianzas y covarainzas de todos cos parñametros del modelo
    # Cuando se calculan la cov numéricas no de calculas las covarianzas de los sigma (solo de los deltas y los betas)
    # Propongo poner un aviso que diga que para obtener los test marginales se seleccionen las cov-no-numericas
    LMM <- NULL
  }
  if (con$trace) {
    end_cov <- proc.time()[3]
    cat("Time to compute covariances: ",
        end_cov - start_cov," seconds \n")
  }
  # Compute R^2 general and for each equation
  Yhat <- z$fitted.values
  R2_pool <- as.numeric((cor(Y,Yhat))^2)
  names(R2_pool) <- c("R2_pool")
  arrYhat <- array(Yhat,c(N,G,Tm))
  arrY <- array(Y,c(N,G,Tm))
  R2_eq <- rep(0,G)
  for (i in 1:G) {
    R2_eq[i] <- cor(matrix(arrY[,i,], ncol = 1),
                    matrix(arrYhat[,i,], ncol = 1))^2
  }
  names(R2_eq) <- paste0("R2_eq", 1:G)
  if (!is.null(R) && !is.null(b)) {
    namesXorig <- colnames(Xorig)
    coefforig <- seorig <- rep(0, ncol(Xorig))
    names(coefforig) <- names(seorig) <- colnames(Xorig)
    coefforig[names(coefficients)] <- coefficients
    seorig[names(rest.se)] <- rest.se
    b <- as.numeric(b)
    for (i in 1:nrow(R)) {
      lidxRi <- R[i,] != 0
      widxRi <- which(lidxRi)
      vidxRi <- R[i,lidxRi]
      ## Check if the constraint is the equality between coefficients
      if ((length(widxRi) == 2) && (sum(vidxRi) == 0) 
          && (b[i] == 0)) {
        coefforig[widxRi[2]] <- coefforig[widxRi[1]]
        seorig[widxRi[2]] <- seorig[widxRi[1]]
        # Updates covariance matrix to include constrained coefficient
        name1 <- names(coefforig)[widxRi[1]]
        name2 <- names(coefforig)[widxRi[2]]
        pr1 <- rbind(resvar, resvar[name1, ])
        rownames(pr1) <- c(rownames(resvar), name2)
        pr2 <- cbind(pr1, c(resvar[, name1], resvar[name1, name1]))
        colnames(pr2) <- rownames(pr2)
        resvar <- pr2
        rm(pr1, pr2)
      }
      # ## Check if the constraint is individual coefficient = 0
      # if ((length(widxRi) == 1) && (vidxRi == 0)&& (b[i] == 0)) {
      #   coefforig[widxRi] <- seorig[widxRi] <- 0
      # }
    }
    X <- Xorig
    p <- porig
    coefficients <- coefforig
    rest.se <- seorig
  }  
  
  ret <- new_spsur(list(call = cl, type = type, 
                       method = method, Durbin = Durbin, 
                       G = G, N = N, Tm = Tm, 
                       deltas = deltas, deltas.se = deltas.se,  
                       coefficients = coefficients, rest.se = rest.se,
                       resvar = resvar, fdHess = fdHess,
                       p = p, dvars = dvars,
                       parameters = parameters,
                       LL = LL, R2 = c(R2_pool,R2_eq),
                       Sigma = Sigma,
                       BP = BP, LMM = LMM,
                       residuals = z$residuals, df.residual = df.residual,
                       fitted.values = z$fitted.values, se.fit = NULL,
                       Y = Y, X = X, W = W, 
                       similar = similar, can.sim = can.sim, 
                       zero.policy = zero.policy, listw_style = listw$style, 
                       interval = interval,  
                       insert = !is.null(trs)))
  if (zero.policy) {
    zero.regs <- attr(listw$neighbours, "region.id")[which(
                      spdep::card(listw$neighbours) == 0)]
    if (length(zero.regs) > 0L) 
      attr(ret, "zero.regs") <- zero.regs
  }
  if(exists("na.act")) { # It could don't exist with data matrices
    if (!is.null(na.act)) ret$na.action <- na.act
  } 
  ret
}

