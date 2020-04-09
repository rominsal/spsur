#' @name spsurml
#' @rdname spsurml
#'
#' @title Maximum likelihood estimation of spatial SUR model.
#'
#' @description This function estimates spatial SUR models using maximum-likelihood methods.
#'  The number of equations, time periods and cross-sectional units is not restricted.
#'  The user can choose between different spatial specifications as described below.
#'  The estimation procedure allows for the introduction of linear restrictions
#'  on the \eqn{\beta} parameters associated to the regressors.
#'
#'
#' @param formula An object created with the package \code{\link[Formula]{Formula}}
#' that describes the model to be estimated. This model may contain several
#' responses (explained variables) and a varying number of regressors in each equation.
#' @param data An object of class data.frame or a matrix.
#' @param W A spatial weighting matrix of order \emph{(NxN)}, assumed to be the
#' same for all equations and time periods.
#' @param Y A column vector of order \emph{(NTmGx1)}, with the observations of the
#' explained variables. The ordering of the data must be (first) equation,
#' (second) time dimension and (third) Cross-sectional/spatial units.
#' The specification of \emph{Y} is only necessary if not available a \code{\link[Formula]{Formula}}
#' and a data frame. Default = \code{NULL}.
#' @param X A data matrix of order \emph{(NTmGxp)} with the observations of the regressors
#' The number of covariates in the SUR model is p = \eqn{sum(p_{g})} where \emph{\eqn{p_{g}}}
#' is the number of regressors (including the intercept) in the g-th equation, \emph{g = 1,...,G}).
#' The specification of \emph{X} is only necessary if not available a \code{\link[Formula]{Formula}}
#' and a data frame. Default = \code{NULL}.
#' @param p Number of regressors by equation, including the intercept. \emph{p} can be a row
#' vector of order \emph{(1xG)}, if the number of regressors is not the same for all the
#' equations, or a scalar, if the \emph{G} equations have the same number of regressors.
#' The specification of \emph{p} is only necessary if not available a \code{\link[Formula]{Formula}}
#' and a data frame.
#' @param G Number of equations.
#' @param N Number of cross-section or spatial units
#' @param Tm Number of time periods.
#' @param type Type of spatial model specification: \strong{"sim"},\strong{"slx"}, \strong{"slm"},
#'  \strong{"sem"}, \strong{"sdm"}, \strong{"sdem"} or \strong{"sarar"}. Default = \code{"sim"}.
#' @param demean  Logical value to allow for the demeaning of panel data. In this case,
#'  \code{\link{spsurml}} substracts the individual mean to each spatial or cross-sectional
#'  unit. Default = \code{FALSE}.
#' @param cov Logical value to show the covariance matrix of the \emph{beta} coefficients.
#'   Default = \code{TRUE}.
#' @param R A row vector of order \emph{(1xpr)} with  the set of \emph{r} linear constraints
#'  on the \emph{beta} parameters. The \emph{first} restriction appears in the first \emph{p} terms,
#'  the second restriction in the next \emph{p} terms and so on. Default = \code{NULL}.
#' @param b A column vector of order \emph{(rx1)} with the values of the linear restrictions on the
#' \emph{beta} parameters. Default = \code{NULL}.
#' @param control List of additional control arguments.
#'
#' @details
#'  The list of (spatial) models that can be estimated with the \emph{spsurml} function are:
#'  \itemize{
#'     \item \strong{"sim"}: SUR model with no spatial effects
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"slx"}: SUR model with spatial lags of the regressors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item \strong{"slm"}: SUR model with spatial lags of the explained variables
#'       \deqn{y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"sem"}: SUR model with spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} Wu_{tg} + \epsilon_{tg} }
#'     \item \strong{"sdm"}: SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tt} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item \strong{"sdem"}: SUR model with spatial lags of the regressors and spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + WX_{tg} \theta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
#'     \item \strong{"sarar"}: SUR model with spatial lags of the explained variables and spatial
#'       errors
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
#'   }
#'
#' @return Output of the  maximum-likelihood estimation of the specified spatial SUR model. A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of model specified. \cr
#'     \code{coefficients} \tab Estimated coefficients for the regressors. \cr
#'     \code{deltas} \tab Estimated spatial coefficients. \cr
#'     \code{rest.se} \tab Estimated standard errors for the estimates of \emph{beta}. \cr
#'     \code{deltas.se} \tab Estimated standard errors for the estimates of the spatial coefficients. \cr
#'     \code{cov} \tab Estimated covariance matrix for the estimates of \emph{beta's} and spatial coefficients.\cr
#'     \code{LL} \tab Value of the likelihood function at the maximum-likelihood estimates. \cr
#'     \code{R2} \tab Coefficient of determination for each equation, obtained as the squared
#'      of the correlation coefficient between the corresponding explained variable and
#'       its estimate. \emph{spsurml} also shows a \emph{global} coefficient of
#'        determination obtained, in the same manner, for the set of \emph{G} equations. \cr
#'     \code{Sigma} \tab Estimated covariance matrix for the residuals of the \emph{G} equations. \cr
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom for the residuals. \cr
#'     \code{fitted.values} \tab Estimated values for the dependent variables. \cr
#'     \code{BP} \tab Value of the Breusch-Pagan statistic to test the null hypothesis
#'      of diagonality among the errors of the \emph{G} equations. \cr
#'     \code{LMM} \tab Marginal Lagrange Multipliers, LM(\eqn{\rho}|\eqn{\lambda}) and
#'       LM(\eqn{\lambda}|\eqn{\rho}), to test for omitted spatial effects in the specification. \cr
#'     \code{G} \tab Number of equations. \cr
#'     \code{N} \tab Number of cross-sections or spatial units. \cr
#'     \code{Tm} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation (including intercepts). \cr
#'     \code{demean} \tab Logical value used for demeaning. \cr
#'     \code{y} \tab Vector \emph{Y} of the explained variables of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the regressors of the SUR model. \cr
#'      \code{W} \tab Spatial weighting matrix. \cr
#'  }
#'
#' @section Control arguments:
#'   \tabular{ll}{
#'     \code{tol} \tab Numerical value for the tolerance for the estimation algorithm until
#'      convergence. Default = 1e-3. \cr
#'     \code{maxit} \tab Maximum number of iterations until convergence; it must be an integer
#'      value. Default = 200. \cr
#'     \code{trace} \tab A logical value to show intermediate results during the estimation process.
#'       Default = \emph{TRUE}. \cr
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
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'      \item Breusch T, Pagan A (1980) The Lagrange multiplier test and its
#'        applications to model specification in econometrics. \emph{Rev Econ Stud} 47: 239-254
#'   }
#'
#' @seealso
#'
#' \code{\link{spsur3sls}}, \code{\link{lmtestspsur}}, \code{\link{wald_betas}}, \code{\link{lrtestspsur}}
#'
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' ## It usually requires 2-3 minutes maximum...
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' spcsur.sim <-spsurml(formula = Tformula, data = spc, type = "sim", W = Wspc)
#' summary(spcsur.sim)
#'
#' ## A SUR-SLX model
#' spcsur.slx <-spsurml(formula = Tformula, data = spc, type = "slx", W = Wspc)
#' summary(spcsur.slx)
#' \donttest{
#' ## A SUR-SLM model
#' spcsur.slm <-spsurml(formula = Tformula, data = spc, type = "slm", W = Wspc)
#' summary(spcsur.slm)
#' rm(spcsur.slm) # remove
#'
#' ## A SUR-SEM model
#' spcsur.sem <-spsurml(formula = Tformula, data = spc, type = "sem", W = Wspc)
#' summary(spcsur.sem)
#' rm(spcsur.sem) # remove
#'
#' ## A SUR-SDM model
#' spcsur.sdm <-spsurml(formula = Tformula, data = spc, type = "sdm", W = Wspc)
#' summary(spcsur.sdm)
#' rm(spcsur.sdm) # remove
#'
#' ## A SUR-SDEM model
#' spcsur.sdem <-spsurml(formula = Tformula, data = spc, type = "sdem", W = Wspc)
#' summary(spcsur.sdem)
#' rm(spcsur.sdem) # remove
#'
#' ## A SUR-SARAR model
#' spcsur.sarar <-spsurml(formula = Tformula, data = spc, type = "sarar", W = Wspc)
#' summary(spcsur.sarar)
#' rm(spcsur.sarar) # remove
#' }
#'
#' #################################################
#' ########  G=1; Tm>1         ########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # \url{https://geodacenter.github.io/data-and-lab/ncovr/}
#' \donttest{
#' ## It usually requires 1-2 minutes maximum...
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SIM model
#' NCOVRSUR.sim <-spsurml(formula = Tformula, data = NCOVR, type = "sim", W = W)
#' summary(NCOVRSUR.sim)
#' rm(NCOVRSUR.sim)
#'
#' ## A SUR-SLX model
#' NCOVRSUR.slx <-spsurml(formula = Tformula, data = NCOVR, type = "slx", W = W)
#' summary(NCOVRSUR.slx)
#' rm(NCOVRSUR.slx)
#' }
#'\donttest{
#' ## It usually requires 1-2 minutes maximum...
#' ## A SUR-SLM model
#' NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR, type = "slm", W = W)
#' summary(NCOVRSUR.slm)
#' rm(NCOVRSUR.slm)
#'
#' ## A SUR-SDM model
#' NCOVRSUR.sdm <-spsurml(formula = Tformula, data = NCOVR, type = "sdm", W = W)
#' summary(NCOVRSUR.sdm)
#' rm(NCOVRSUR.sdm)
#'
#' ## A SUR-SEM model
#' NCOVRSUR.sem <-spsurml(formula = Tformula, data = NCOVR, type = "sem", W = W)
#' summary(NCOVRSUR.sem)
#' ## A SUR-SDEM model
#' NCOVRSUR.sdem <-spsurml(formula = Tformula, data = NCOVR, type = "sdem",W = W)
#' summary(NCOVRSUR.sdem)
#'
#' ## A SUR-SARAR model
#' NCOVRSUR.sarar <-spsurml(formula = Tformula, data = NCOVR,
#'                           type = "sarar", W = W)
#'  summary(NCOVRSUR.sarar)
#' }
#'
#' ##############################################
#' ######## SUR with G>1; Tm>1  ###
#' ##############################################
#' \donttest{
#' ## It usually requires 2-3 minutes maximum...
#' rm(list = ls())  # Clean memory
#' #### Reshape NCOVR in panel format
#' data(NCOVR,package="spsur")
#' N <- nrow(NCOVR)
#' Tm <- 4
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' pHR <- c(NCOVR$HR60, NCOVR$HR70, NCOVR$HR80, NCOVR$HR90)
#' pPS <- c(NCOVR$PS60, NCOVR$PS70, NCOVR$PS80, NCOVR$PS90)
#' pUE <- c(NCOVR$UE60, NCOVR$UE70, NCOVR$UE80, NCOVR$UE90)
#' pDV <- c(NCOVR $DV60, NCOVR$DV70, NCOVR$DV80, NCOVR$DV90)
#' pFP <- c(NCOVR$FP59, NCOVR$FP70, NCOVR$FP80, NCOVR$FP90)
#' pSOUTH <- rep(NCOVR$SOUTH, Tm)
#' pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
#'                      HR = pHR, PS = pPS, UE = pUE, DV = pDV,
#'                      FP = pFP, SOUTH = pSOUTH)
#' rm(NCOVR) # Free memory...
#' pform <- HR | DV | FP ~ PS + UE | PS + UE + SOUTH | PS
#'
#' ## SIM (easy to compute...)
#'
#' psur_sim <- spsurml(formula = pform, data = pNCOVR, W = W, type = "sim")
#' summary(psur_sim)
#'
#' ## SLM (cov = FALSE to prevent overflows of memory)
#' psur_slm <- spsurml(formula = pform, data = pNCOVR, W = W,
#'                         type = "slm", cov = FALSE)
#' psur_slm$deltas
#' psur_slm$betas
#' psur_slm$Sigma_corr
#' rm(psur_slm)
#' ## SEM (cov = FALSE to prevent overflows of memory)
#' ### Only execute if you have enough memory...
#' psur_sem <- spsurml(formula = pform, data = pNCOVR, W = W,
#'                         type = "sem", cov = FALSE)
#' psur_sem$deltas
#' psur_sem$betas
#' psur_sem$Sigma_corr
#' psur_sem
#' }
#'
#' ##############################################
#' ######## Demeaning Examples with G>1; Tm>>1  ###
#' ##############################################
#' \donttest{
#' #' rm(list = ls())  # Clean memory
#' set.seed(123456)
#' Tm <- 10 # Number of time periods
#' G <- 3 # Number of equations
#' N <- 100 # Number of spatial elements
#' p <- 3 # Number of independent variables
#' Sigma <- matrix(0.5, ncol = G, nrow = G)
#' diag(Sigma) <- 1
#' Betas <- rep(1:3, G)
#' rho <- 0.5
#' lambda <- 0.0 # spatial autocorrelation error term = 0
#' #  random coordinates
#' co <- cbind(runif(N,0,1),runif(N,0,1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co, k = 5, longlat = FALSE)))
#' DGPsim <- dgp_spsur(Sigma = Sigma, Betas = Betas, rho = rho, lambda = lambda,
#'                  Tm = Tm, G = G, N = N, p = p, W = W)
#'
#' ## SLM without demeaning
#'
#' SUR_slm  <-spsurml(Y= DGPsim$Y, X = DGPsim$X, G = G, N = N, Tm = Tm,
#'                   p = p, W = W, type = "slm")
#' summary(SUR_slm)
#'
#' # SLM with demeaning
#'
#' SUR_slm_dem  <-spsurml(Y= DGPsim$Y, X = DGPsim$X, G = G, N = N, Tm = Tm,
#'                   p = p, W = W, type = "slm", demean = TRUE)
#' summary(SUR_slm_dem)
#' }
#' @export
spsurml <- function(formula = NULL, data = NULL, na.action, 
                    R = NULL, b = NULL, listw = NULL, 
                    quiet = NULL, zero.policy = NULL, 
                    interval = NULL, trs = NULL,
                    X = NULL, Y = NULL, 
                    G = NULL, N = NULL, Tm = NULL,
                    p = NULL, demean = FALSE,
                    type = "sim", method = "eigen",
                    control = list() ) {
  # Spatial Models: sim, slx, slm, sdm, sem, sdem, sarar
  
  con <- list(tol = 0.05, maxit = 200, trace = TRUE, 
              tol.opt = .Machine$double.eps^0.5, fdHess = NULL,
              optimHess = FALSE, optimHessMethod = "optimHess",
              compiled_sse = FALSE,
              Imult = 2, cheb_q = 5, MC_p = 16L, MC_m = 30L, super = NULL,
              spamPivot = "MMD", in_coef = 0.1, type = "MC", correct = TRUE,
              trunc = TRUE, SE_method = "LU", nrho = 200, interpn = 2000,
              small_asy = TRUE, small = 1500, SElndet = NULL, LU_order = FALSE,
              pre_eig = NULL, OrdVsign = 1)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]))
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  if (is.null(quiet)) quiet <- TRUE
  stopifnot(is.logical(quiet))
  
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
    # Change dimensions in this case with matrix Data
    G <- Tm
    Tm <- 1
  }
  
  if (!is.null(formula) && !any(class(formula) == "Formula")) 
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
    if (any(type == c("sdm","sdem","slx"))) {
      Durbin <- TRUE } else { Durbin <- FALSE }
    get_XY <- get_data_spsur(formula = formula, mf = mf, Durbin = Durbin,
                             listw = listw, zero.policy = zero.policy, 
                             N = N)
    Y <- get_XY$Y
    X <- get_XY$X
    G <- get_XY$G
    N <- get_XY$N
    Tm <- get_XY$Tm
    p <- get_XY$p
    dvars <- get_XY$dvars
    if (Tm > 1 && G == 1){
       # Change dimensions in this case with Matrix Data
       G <- Tm
       Tm <- 1
     }
    rm(get_XY)
    if (length(p) == 1) p <- rep(p,G)
  }
  if (length(p) == 1) p <- rep(p,G)
  names(p) <- NULL

  if (!is.null(R) & !is.null(b)) {
    demean <- FALSE   # Demeaning isn't allowed in restricted case
    restr <- X_restr(X = X, R = R, b = b, p = p)
    X <- restr$Xstar
    p <- restr$pstar
  }
  
  if (Tm < 2) demean <- FALSE # Demeaning is not allowed for Tm < 2
  if (demean) {
    demeanXY <- demeaning(X = X, Y = Y, G = G, N = N, Tm = Tm, p = p)
    X <- demeanXY$Xdem
    Y <- demeanXY$Ydem
    p <- demeanXY$pdem
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
    assign("verbose", !quiet, envir = env)
    assign("family", "SAR", envir = env) # CHEQUEAR OTRAS OPCIONES
    interval <- spatialreg::jacobianSetup(method, env, con, 
                                          pre_eig = con$pre_eig,
                                          trs = trs, interval = interval)
    assign("interval", interval, envir = env)
  }
  if (any(type == c("sim","slm","sem","sarar")))
    name_fit <- paste("fit_spsur", type, sep = "")
  if (type == "sdm") name_fit <- "fit_spsurslm"
  if (type == "sdem") name_fit <- "fit_spsursem"
  if (type == "slx") name_fit <- "fit_spsursim"
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
    if (type == "sarar") {
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
      name_cov_fit <- "cov_spsursim"
    if (any(type == c("slm","sem","sarar")))
      name_cov_fit <- paste("cov_spsur", type, sep = "")
    if (type == "sdm") name_cov_fit <- "cov_spsurslm"
    if (type == "sdem") name_cov_fit <- "cov_spsursem"
    cov_fit <- get(name_cov_fit)
    #z_cov <- cov_fit(Tm=Tm,G=G,N=N, 
    #      Y=Y,X=X,W=W,
    #      deltas=deltas,Sigma=z$Sigma,trace=trace)
    allcov <- try( cov_fit(env = env) ) 
    ## *****VIP****: VÁLIDO PARA types == any("sim","slx",slm"). 
    ## HAY QUE REHACER CÓDIGO COVARIANZAS PARA RESTO TIPOS (SEM Y SARAR)
    if (class(allcov) == "try-error") {
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
        ## OJO: NO SE TIENEN EN CUENTA COVARIANZAS
        ## ENTRE SIGMA Y COEFFICIENTS...
      } else {
        resvar <- allcov$vcov
        colnames(resvar) <- rownames(resvar) <-
                c(names(coefficients),
                  names(deltas),names_sigma)
        ## VIP: CAMBIO ORDEN MATRIZ COVARIANZAS
        ## IGUAL ORDEN QUE SPDEP Y SPATIALREG PACKAGES...
        resvar <- resvar[c(names_sigma,names(deltas),names(coefficients)),
                         c(names_sigma,names(deltas),names(coefficients))]
      }
    }
  }
  if (fdHess) {
    cat("Computing numerical covariances...","\n")
    if (any(type == c("slm","sdm"))) name_cov_fit <- "f_sur_lag"
    if (any(type == c("sem","sdem"))) name_cov_fit <- "f_sur_sem"
    if (type == "sarar") name_cov_fit <- "f_sur_sarar"
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
  names(R2_eq) <- paste0("R2_eq",1:G)
  ret <- structure(list(call = cl, type = type, method = method,
                       G = G, N = N, Tm = Tm, 
                       deltas = deltas, deltas.se = deltas.se,  
                       coefficients = coefficients, rest.se = rest.se,
                       resvar = resvar, fdHess = fdHess,
                       p = p, dvars = dvars,
                       parameters = parameters,
                       LL = LL, R2 = c(R2_pool,R2_eq),
                       Sigma = Sigma, 
                       residuals = z$residuals, df.residual = df.residual,
                       fitted.values = z$fitted.values, se.fit = NULL,
                       y = Y, X = X, W = W, 
                       similar = similar, can.sim = can.sim, 
                       demean = demean,   
                       zero.policy = zero.policy, listw_style = listw$style, 
                       interval = interval,  
                       optimHess = con$optimHess, 
                       insert = !is.null(trs), trs = trs), 
                       class = c("spsur"))
  if (zero.policy) {
    zero.regs <- attr(listw$neighbours, "region.id")[which(
                      spdep::card(listw$neighbours) == 0)]
    if (length(zero.regs) > 0L) 
      attr(ret, "zero.regs") <- zero.regs
  }
  if (!is.null(na.act)) ret$na.action <- na.act
  ret  
}

