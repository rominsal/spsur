#' @name spsurgs3sls
#' @rdname spsurgs3sls
#'
#' @title General Spatial 3SLS for systems of spatial equations. 
#' 

#' @description The function estimates spatial SUR models using 
#' general spatial three stages least squares. This is a 
#' system instrumental variable procedure which also include 
#' GMM estimation when there is spatial correlations in the errors.
#' The procedure allows for additional endogenous regressors in
#' addition to spatial lags of the dependent variable. It could 
#' be applied to "slm", "sdm", "sem" and "sarar" spatial models. 
#' Furthermore, for non-spatial models including endogenous 
#' regressors ("iv"), it could be used to estimate using 
#' instrumental variables and Feasible Generalized Least Squares.
#'  
#' @usage spsurgs3ls (formula = NULL, data = NULL, na.action,
#'                   listw = NULL, zero.policy = NULL, 
#'                   type = "slm", Durbin = FALSE,
#'                   endog = NULL, instruments = NULL,
#'                   lag.instr = FALSE, initial.value = 0.2, 
#'                   het = FALSE, trace = TRUE)
#'     
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
#' @param zero.policy Similar to the corresponding parameter of 
#'   \code{\link[spatialreg]{lagsarlm}} function in \pkg{spatialreg} package. 
#' (NxN) instead of a listw object. Default = NULL.
#' @param type Type of spatial model specification: "sim", "iv",
#'   "slm", "sem", "sdm" or "sarar" . Default = "slm".
#' @param Durbin If a formula object and model is type "sdm"  
#'   the subset of explanatory variables to lag for each equation.   
#' @param endog Additional endogenous variables. Default \emph{NULL}. 
#'   If not \emph{NULL} should be specified as a 
#'   \code{\link[Formula]{Formula}} with no dependent variable. 
#'   Examples: ~ x1 | x2 (x1 endogeous regressor for the first 
#'   equation and x2 endogeneous regressor for the second equation) 
#'   or ~ x1 | . (x1 endogenous regressor for the first equation and 
#'   none endogenous regressors for the second equation)
#' @param instruments external instruments. Default \emph{NULL}. If not 
#'   \emph{NULL} should be specified as a formula with no 
#'   dependent variable in the same way than 
#'   previous \emph{endog} argument.
#' @param lag.instr should the external instruments be spatially lagged?      
#' @param initial.value he initial value for ρ. It can be either numeric
#'  (default is 0.2) or set to 'SAR', in which case the optimization 
#'  will start from the estimated coefficient of a regression of 
#'  the 2SLS residuals over their spatial lag (i.e. a spatial 
#'  AR model)
#' @param het default {FALSE}: if {TRUE} uses the methods 
#'   developed for heteroskedasticity for each equation. 
#'   Wrapper using \code{\link[sphet]{spreg}} function.
#' @param trace A logical value to show intermediate results during 
#'       the estimation process. Default = \code{TRUE}.
#'              
#' @details
#'  \emph{spsurg3sls} generalize the \code{\link[sphet]{spreg}} function
#'  to multiequational spatial SUR models. The methodology to estimate
#'  spatial SUR models by Generalized 3SLS follows the steps outlined in
#'  Kelejian and Piras (pp. 304-305). The summary of the algorithm is
#'  the next one:
#'   \itemize{ 
#'    \item Estimate each equation by 2SLS and obtain the estimated 
#'     residuals \eqn{\hat{u}_j} for each equation. 
#'    \item If the model includes a spatial lag for the errors. 
#'     (that is, it is a SEM/SARAR model), apply GMM to obtain 
#'     the spatial parameters \eqn{\hat{\lambda}_j} for the residuals
#'     in each equation. In this case the \code{\link[sphet]{spreg}} 
#'     function is used as a wrapper for the GMM estimation. If the
#'     model does not include a spatial lag for the errors (that is, 
#'     it is a "sim", "iv", "slm" or "sdm" model), then \eqn{\hat{\lambda}_j = 0}
#'     \item Compute 
#'      \deqn{\hat{v}_j = \hat{u}_j-\hat{\lambda}_j W \hat{u}_j} 
#'      and the covariances 
#'      \deqn{\hat{\sigma}_{i,j} = N^{-1}\hat{v}_i\hat{v}_j}. 
#'      Build \eqn{\hat{Sigma}=\lbrace \hat{\sigma_{i,j}} \rbrace}
#'    \item Compute 
#'      \deqn{y_j^* = y_j - \hat{\lambda}_j W y_j}
#'     and
#'      \deqn{X_j^* = X_j- \hat{\lambda}_j W X_j}
#'     Compute 
#'      \deqn{ \hat{X}_j^* = H_j(H_j^T H_j)^{-1} H_j^T X_j^*}      
#'     where \eqn{H_j} is the matrix including all the instruments and
#'     the exogenous regressors for each equation. That is,
#'     \eqn{\hat{X}_j^*} 
#'     is the projection of \eqn{X_j^*} using
#'     the instruments matrix \eqn{H_j}. 
#'    \item Compute, in a multiequational way, 
#'     the Feasible Generalized Least Squares estimation using the 
#'     new variables \eqn{\hat{y}_j^*}, \eqn{\hat{X}_j^*} and 
#'     \eqn{\hat{Sigma}}. This is the 3sls step.
#'   }
#' @return Object of \code{spsur} class with the output of the  three-stages 
#'   least-squares estimation of the specified spatial model.
#'   A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of model specified. \cr
#'     \code{Durbin} \tab Value of \code{Durbin} argument. \cr
#'     \code{coefficients} \tab Estimated coefficients for the regressors. \cr
#'     \code{deltas} \tab Estimated spatial coefficients. \cr
#'     \code{rest.se} \tab Estimated standard errors for the estimates of 
#'       \emph{\eqn{\beta}} coefficients. \cr
#'     \code{deltas.se} \tab Estimated standard errors for the estimates of 
#'       the spatial coefficients. \cr
#'     \code{resvar} \tab Estimated covariance matrix for the estimates of 
#'       \emph{beta's} and spatial coefficients.\cr
#'     \code{R2} \tab Coefficient of determination for each equation, 
#'       obtained as the squared of the correlation coefficient between 
#'       the corresponding explained variable and its estimates. 
#'       \emph{spsur3sls} also shows a \emph{global} coefficient of
#'        determination obtained, in the same manner, for the set of 
#'        \emph{G} equations. \cr
#'     \code{Sigma} \tab Estimated covariance matrix for the residuals of the 
#'       \emph{G} equations. \cr
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom for the residuals. \cr
#'     \code{fitted.values} \tab Estimated values for the dependent variables. \cr
#'     \code{G} \tab Number of equations. \cr
#'     \code{N} \tab Number of cross-sections or spatial units. \cr
#'     \code{Tm} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation (including intercepts). \cr
#'     \code{Y} \tab Vector \emph{Y} of the explained variables of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the regressors of the SUR model. \cr
#'     \code{W} \tab Spatial weighting matrix. \cr
#'     \code{zero.policy} \tab Logical value of \code{zero.policy} . \cr
#'     \code{listw_style} \tab	Style of neighborhood matrix \code{W}. \cr
#'  }
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
#'    \item Kelejian, H. H. and Piras, G. (2017). 
#'      \emph{ Spatial Econometrics}. Academic Press.
#'    \item Kelejian, H.H. and Prucha, I.R. (2010).
#'     Specification and Estimation of Spatial 
#'     Autoregressive Models with Autoregressive 
#'     and Heteroskedastic Disturbances. 
#'     \emph{Journal of Econometrics}, 157, 
#'     pp. 53-67.
#'     \item Kelejian, H.H. and Prucha, I.R. (1999). 
#'       A Generalized Moments Estimator for the 
#'       Autoregressive Parameter in a Spatial 
#'       Model. \emph{International Economic Review},
#'       40, pp. 509-533.
#'     \item Kelejian, H.H. and Prucha, I.R. (1998).
#'       A Generalized Spatial Two Stage Least 
#'       Square Procedure for Estimating a Spatial 
#'       Autoregressive Model with Autoregressive 
#'       Disturbances. \emph{Journal of Real 
#'       Estate Finance and Economics}, 17, 
#'       pp. 99–121.
#'     \item Piras, G. (2010). sphet: Spatial
#'       Models with Heteroskedastic Innovations 
#'       in R. \emph{Journal of Statistical 
#'       Software}, 35(1), pp. 1-21. 
#'       \url{http://www.jstatsoft.org/v35/i01/}. -          
#'   }
#'
#' @seealso
#' \code{\link[sphet]{spreg}},
#' \code{\link{spsur3sls}},
#' \code{\link[spatialreg]{stsls}},
#' \code{\link{spsurml}} 
#' 
#' @examples
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data(spc)
#' lwspc <- spdep::mat2listw(Wspc, style = "W")
#' ## No endogenous regressors  
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## Endogenous regressors and Instruments
#' Tformula2 <- WAGE83 | WAGE81 ~  NMR83 | NMR80 
#' ## Endogenous regressors: UN83 , UN80
#' ## Instrumental variable: SMSA
#'
#' ## A IV model with endogenous regressors only in first equation
#' spciv <- spsurgs3sls(formula = Tformula2, data = spc,
#'                       type = "iv", listw = lwspc,
#'                       endog = ~ UN83 | ., 
#'                       instruments = ~ SMSA | .)
#' summary(spciv)
#' print(spciv)#
#' #########################################################################
#' ## A SLM model with endogenous regressors 
#' spcslm <- spsurgs3sls(formula = Tformula2, data = spc,
#'                               endog = ~ UN83 | ., 
#'                               instruments = ~ SMSA |.,
#'                               type = "slm", 
#'                               listw = lwspc)
#' summary(spcslm)
#' print(spcslm)                           
#' impacts_spcslm <- impactspsur(spcslm, listw = lwspc, R = 1000)
#' summary(impacts_spcslm[[1]], zstats = TRUE, short = TRUE)
#' summary(impacts_spcslm[[2]], zstats = TRUE, short = TRUE)
#' #########################################################################
#' ## A SDM model with endogenous regressors 
#' spcsdm <- spsurgs3sls(formula = Tformula2, data = spc,
#'                    endog = ~ UN83 | UN80, 
#'                    instruments = ~ SMSA | SMSA,
#'                    type = "sdm", listw = lwspc,
#'                    Durbin =  ~  NMR83 | NMR80)
#' summary(spcsdm)
#' ## Durbin only in one equation
#' spcsdm2 <- spsurgs3sls(formula = Tformula2, data = spc,
#'                     endog = ~ UN83 | UN80, 
#'                     instruments = ~ SMSA | SMSA,
#'                     type = "sdm", listw = lwspc,
#'                     Durbin =  ~  NMR83 | .)
#' summary(spcsdm2)
#' #########################################################################
#' ## A SEM model with endogenous regressors 
#' spcsem <- spsurgs3sls(formula = Tformula2, data = spc,
#'                       endog = ~ UN83 | UN80, 
#'                       instruments = ~ SMSA | SMSA,
#'                       type = "sem", listw = lwspc)
#' summary(spcsem)
#' print(spcsem)                           
#' #########################################################################
#' ## A SARAR model with endogenous regressors 
#' spcsarar <- spsurgs3sls(formula = Tformula2, data = spc,
#'                         endog = ~ UN83 | UN80, 
#'                         instruments = ~ SMSA | SMSA,
#'                         type = "sarar", listw = lwspc)
#' summary(spcsarar)
#' print(spcsarar)                           
#' impacts_spcsarar <- impactspsur(spcsarar, listw = lwspc, R = 1000)
#' summary(impacts_spcsarar[[1]], zstats = TRUE, short = TRUE)
#' summary(impacts_spcsarar[[2]], zstats = TRUE, short = TRUE)
#' @export
spsurgs3sls <- function(formula = NULL, data = NULL, na.action, 
                      listw = NULL, zero.policy = NULL, 
                      type = "slm", Durbin = FALSE,  
                      endog = NULL, instruments = NULL,
                      lag.instr = FALSE, initial.value = 0.2,
                      het = FALSE, trace = TRUE) {
  if (is.null(formula) || is.null(data))
    stop("formula and data arguments need to be specified")
  if (type == "slm" || type == "sdm") model_spreg = "lag" 
  else if (type == "sem") model_spreg = "error" 
  else if (type == "sarar") model_spreg = "sarar" 
  else if (type == "sim") model_spreg = "ols"
  else if (type == "iv") model_spreg = "ivhac"
  else stop("spsur_spreg can only be used with slm, sdm, sem, sarar, 
            iv or sim models")
  if (is.null(listw) || 
      !inherits(listw, c("listw", "Matrix", "matrix")))
    stop("listw format unknown or NULL")
  if (inherits(listw, "matrix")) {
    W <- Matrix::Matrix(listw)
    listw <- spdep::mat2listw(W)
  }  
  if (inherits(listw, "Matrix")) {
    W <- listw
    listw <- spdep::mat2listw(as.matrix(W))
  }
  if (is.null(zero.policy))
    zero.policy <- spatialreg::get.ZeroPolicyOption()
  cl <- match.call()
  if (!inherits(formula, "Formula")) 
    formula <- Formula::Formula(formula)
  if (!is.null(endog)) {
    if (!inherits(endog, "Formula")) 
      endog <- Formula::Formula(endog)
    if (is.null(instruments)) {
      stop("Instruments needed to be supplied for endogenous variables") 
    } else {
      if (!inherits(instruments, "Formula")) 
        instruments <- Formula::Formula(instruments)
    }
  }
  ####### VIP: Assumption (G>1 && Tm==1) || (G==1 && Tm>1)
  if (trace)  start_fit <- proc.time()[3]
  G <- length(formula)[1] # number of lhs's 
  Tm <- 1
  lspreg <- vector(mode = "list", length = G)
  yg <- Xg <- Zg <- Phig <- Hg <- vector("list", length = G)
  uhatg <- vector("list", length = G)
  betashatg <- coeffg <- vector(mode = "list", length = G)
  VChatg <- vector(mode = "list", length = G)
  rhog <- lambdag <- vector(mode = "list", length = G)
  Ng <- pg <- vector(mode = "integer", length = G)
  ## Notation: Coefficients similar to "Spatial Econometrics" (Chapter 14)
  for (i in 1:G) {
    if (trace) cat("Fitting equation ",i,"...","\n")
    form_i <- formula(formula, lhs = i, rhs = i)
    if (!is.null(endog) && !is.null(instruments)) {
      endog_i <- formula(endog, rhs = i)
      if (endog_i == ~.) endog_i <- NULL
      instruments_i <- formula(instruments, rhs = i)
      if (instruments_i == ~.) instruments_i <- NULL
    } else {
      endog_i <- instruments_i <- NULL
    }
    if (!is.null(endog_i) && is.null(instruments_i))
      stop("Endogenous regressors without instruments")
    if (inherits(Durbin, "formula")) {
      if (!inherits(Durbin, "Formula")) 
        Durbin <- Formula::Formula(Durbin)
      Durbin_i <- formula(Durbin, rhs = i)
      if (Durbin_i == ~.) Durbin_i <- FALSE
    } else Durbin_i <- FALSE
  ## Step 1,2,3 textbook "Spatial Econometrics", pp. 305
    lspreg[[i]] <- sphet::spreg(formula = form_i, data = data, 
                        listw = listw, na.action = na.action,
                        endog = endog_i, 
                        instruments = instruments_i,
                        lag.instr = lag.instr,
                        initial.value = initial.value,
                        model = model_spreg, het = het, 
                        verbose = FALSE, HAC = FALSE, 
                        distance = NULL, type = "Epanechnikov",
                        bandwidth = "variable",
                        step1.c = step1.c,
                        Durbin = Durbin_i)
    coeffg[[i]] <- lspreg[[i]]$coefficients
    VChatg[[i]] <- lspreg[[i]]$var
    namescoeffg <- rownames(coeffg[[i]])
    ## Change lambda rho and lag names to be consistent with 
    ## the rest of spsur package...
    if (any(grepl("lambda", namescoeffg)))
      namescoeffg[grepl("lambda", namescoeffg)] <- "phi1"
    if (any(grepl("rho", namescoeffg)))
      namescoeffg[grepl("rho", namescoeffg)] <- "phi2"
    if (any(grepl("phi1", namescoeffg)))
      namescoeffg[grepl("phi1", namescoeffg)] <- "rho"
    if (any(grepl("phi2", namescoeffg)))
      namescoeffg[grepl("phi2", namescoeffg)] <- "lambda"
    if (any(grepl("lag_", namescoeffg)))
      namescoeffg <- gsub("lag_","lag.", namescoeffg)
    namescoeffg <- paste(namescoeffg, i, sep = "_")
    rownames(coeffg[[i]]) <- namescoeffg
    rownames(VChatg[[i]]) <- namescoeffg
    colnames(VChatg[[i]]) <- namescoeffg
    idxbetashag <- (!grepl("rho", 
                 rownames(coeffg[[i]])) & 
      !grepl("lambda", rownames(coeffg[[i]])))
    betashatg[[i]] <- coeffg[[i]][idxbetashag]
    names(betashatg[[i]]) <- namescoeffg[idxbetashag]
    rhog[[i]] <- coeffg[[i]][grepl("rho", namescoeffg)]
    if (length(rhog[[i]]) < 1) rhog[[i]] <- 0
    lambdag[[i]] <- coeffg[[i]][grepl("lambda", namescoeffg)]
    if (length(lambdag[[i]]) < 1) lambdag[[i]] <- 0
    Ng[[i]] <- nrow(lspreg[[i]]$model) 
    yg[[i]] <- lspreg[[i]]$model[, c("y")]
    if (Durbin == FALSE) {
      Xg[[i]] <- model.matrix(form_i, data = data)
    } else { # Adding lags of X's variables...
      Xg[[i]] <- as.matrix(lspreg[[i]]$model[, -1]) # Exclude y column
      Xg[[i]] <- cbind(1, Xg[[i]]) # Add intercept column
    }
    if (!is.null(endog_i)) {
    ## Exclude intercept from model.matrix for endog. regressors
      Xg[[i]] <- cbind(Xg[[i]], 
                       model.matrix(endog_i, 
                                    data = data)[, -1, 
                                                 drop = FALSE])
    }
    #Add colnames to Xg[[i]]
    colnames(Xg[[i]]) <- names(betashatg[[i]])
    pg[[i]] <- ncol(Xg[[i]])
    Zg[[i]] <- Xg[[i]]
    W <- sphet::listw2dgCMatrix(listw, 
                                zero.policy = zero.policy)
    if (model_spreg == "lag" ||  model_spreg == "sarar") {
      ## Add spatial lag of y to Zg[[i]]
      Wyg <-  W %*% yg[[i]]
      colnames(Wyg) <- paste("Wy", i, sep = "_")
      Zg[[i]] <- cbind(Zg[[i]], Wyg)
    }
    ## Add spatial lag of Xg[[i]] to Hg[[i]]
    if (all(Xg[[i]][, 1] == 1)) { 
    ## Exclude intercept from the spatial lag
      WXg <-  W %*% Xg[[i]][, -1, drop = FALSE]
    } else { WXg <-  W %*% Xg[[i]] }
    colnames(WXg) <- paste("W", colnames(WXg), sep = "")
    Hg[[i]] <- cbind(Xg[[i]], WXg)
    if (!is.null(instruments_i)) {
      Phig[[i]] <- model.matrix(instruments_i, 
                                data = data)[, -1, drop = FALSE]
      colnames(Phig[[i]]) <- paste(colnames(Phig[[i]]), i, 
                                   sep = "_")
      Hg[[i]] <- cbind(Hg[[i]], Phig[[i]])
      if (lag.instr) { 
        WPhig <- W %*% Phig[[i]]
        colnames(WPhig) <- paste("W", colnames(WPhig), sep = "")
        Hg[[i]] <- cbind(Hg[[i]], WPhig[[i]])
      } 
    } else Phig[[i]] <- NULL
    uhatg[[i]] <- lspreg[[i]]$residuals
  } ## end for (i in 1:G)
  # Check equals Ni (assuming G>1 && Tm == 1)
  if (!length(unique(unlist(Ng))) == 1) {
    stop ("Some equations have different N")
  } else N <- unique(Ng)
  yf <- unlist(yg)
  uhatf <- unlist(uhatg)    
  Xf <- Matrix::bdiag(Xg)
  Zf <- Matrix::bdiag(Zg)
  Hf <- Matrix::bdiag(Hg)
  rhof <- Matrix::bdiag(rhog)
  lambdaf <- Matrix::bdiag(lambdag)
  names_colX <- names_colH <- names_colZ <- NULL
  for (i in 1:G) {
    names_colX <- c(names_colX, colnames(Xg[[i]]))
    names_colZ <- c(names_colZ, colnames(Zg[[i]]))
    names_colH <- c(names_colH, colnames(Hg[[i]]))
  }
  colnames(Xf) <- names_colX
  colnames(Zf) <- names_colZ
  colnames(Hf) <- names_colH
  N <- Ng[[1]]
  ## Step 4 textbook "Spatial Econometrics", pp. 305
  vhatf <- (Matrix::Diagonal(N*G) - 
              Matrix::kronecker(lambdaf, W)) %*% uhatf
  Sigma <- get_Sigma(resids = vhatf, N = N, G = G, Tm = Tm)
  Sigmahat <- Matrix::Matrix(Sigma$Sigma)
  Sigmahatinv <- Matrix::Matrix(Sigma$Sigma_inv)
  ## Step 5 textbook "Spatial Econometrics", pp. 305
  yfrho2 <-  (Matrix::Diagonal(N*G) - 
                Matrix::kronecker(lambdaf, W)) %*% yf
  Zfrho2 <-  (Matrix::Diagonal(N*G) - 
                Matrix::kronecker(lambdaf, W)) %*% Zf
  ## Step 6 textbook "Spatial Econometrics", pp. 305
  Zhatrho2 <- Matrix::Matrix(lm.fit(as.matrix(Hf), 
                                    as.matrix(Zfrho2), 
                  singular.ok = TRUE)$fitted.values)
  gammahat <- Matrix::solve(Matrix::t(Zhatrho2) %*% 
                              Matrix::kronecker(Sigmahatinv, 
                                        Matrix::Diagonal(N)) %*% 
                              Zhatrho2, 
                            Matrix::t(Zhatrho2) %*% 
                              Matrix::kronecker(Sigmahatinv, 
                                        Matrix::Diagonal(N)) %*% 
                              yfrho2)
  ## Change names Wy_i by rho_i
  if (any(grepl("Wy", rownames(gammahat)))) {
    rownames(gammahat)[grepl("Wy", 
                        rownames(gammahat))] <- paste("rho", 
                                                   1:G, sep = "_")
  }
  VCgammahat <- Matrix::solve(Matrix::t(Zhatrho2) %*% 
                                Matrix::kronecker(Sigmahatinv, 
                                          Matrix::Diagonal(N)) %*% 
                              Zhatrho2)
  rownames(VCgammahat) <- rownames(gammahat)
  colnames(VCgammahat) <- rownames(gammahat)
  rhohat <- gammahat[grepl("rho", rownames(gammahat))]
  VCrhohat <- VCgammahat[grepl("rho", rownames(VCgammahat)),
                         grepl("rho", rownames(VCgammahat))]
  if (length(rhohat) < 1) rhohat <- rep(0, G)
  if (length(VCrhohat) < 1) VCrhohat <- matrix(0, G, G)
  rhohat.se <- sqrt(Matrix::diag(VCrhohat))
  names(rhohat) <- names(rhohat.se) <- paste("rho", 1:G, 
                                             sep = "_")
  lambdahat <- Matrix::diag(lambdaf)
  VClambdahat <- Matrix(0, G, G)
  for (i in 1:G) {
    lambi <- paste("lambda", i, sep="_")
    if (any(grepl(lambi, rownames(VChatg[[i]]))))
      VClambdahat[i, i] <- VChatg[[i]][lambi, lambi]
  }
  lambdahat.se <- sqrt(Matrix::diag(VClambdahat))
  names(lambdahat) <- names(lambdahat.se) <- paste("lambda", 1:G, 
                                             sep = "_")
  deltas <- c(rhohat[abs(rhohat) > 0], 
              lambdahat[abs(lambdahat) > 0])
  deltas.se <- c(rhohat.se[rhohat.se > 0], 
                 lambdahat.se[lambdahat.se > 0])
  coefficients <- gammahat[!grepl("rho", rownames(gammahat))]
  idxcoeff <- (!grepl("rho", rownames(gammahat)))
  names(coefficients) <- rownames(gammahat)[idxcoeff]
  idxrest.se <- (!grepl("rho", rownames(VCgammahat)))
  rest.se <- sqrt(Matrix::diag(VCgammahat[idxrest.se, 
                                          idxrest.se]))
  names(rest.se) <- names(coefficients)
  Sigma <- as.matrix(Sigmahat)
  names_sigma <- NULL
  for (i in 1:G) {
    new_name <- paste0("sigma", i, sep = "")
    names_sigma <- c(names_sigma, new_name)
  }
  colnames(Sigma) <- rownames(Sigma) <- names_sigma
  resvar <- VCgammahat
  if (any(grepl("lambda", names(deltas)))) {
    ## Add variances of lambda (assumption covariances = 0)
    fullnames <- rownames(resvar)
    resvar <- rbind(resvar, matrix(0, nrow = G, 
                                   ncol = ncol(resvar)))
    resvar <- cbind(resvar, matrix(0, nrow = nrow(resvar), 
                                   ncol = G))
    rownames(resvar) <- c(fullnames, 
                          paste("lambda", 1:G, sep = "_"))
    colnames(resvar) <- rownames(resvar)
    for (i in 1:G) {
      lambi <- paste("lambda", i, sep = "_")
      resvar[lambi, lambi] <- VClambdahat[i, i]
    }
  }
  parameters <- length(coefficients) + 
    length(deltas) + G*(G + 1)/2
  df.residual <- G*N*Tm - parameters
  # Compute R^2 general and for each equation
  fitted.values <-  solve(Matrix::Diagonal(N*G) - 
                            Matrix::kronecker(diag(rhohat), W), 
                          Xf %*% coefficients)
  yhat <- as.numeric(fitted.values)
  residuals <- yf - yhat
  R2_pool <- as.numeric((cor(yf, yhat))^2)
  names(R2_pool) <- c("R2_pool")
  arryhat <- array(yhat, c(N, G, Tm))
  arry <- array(yf, c(N, G, Tm))
  R2_eq <- rep(0,G)
  for (i in 1:G) {
    R2_eq[i] <- cor(matrix(arry[, i, ], ncol = 1),
                    matrix(arryhat[, i, ], ncol = 1))^2
  }
  names(R2_eq) <- paste0("R2_eq", 1:G)
  R2 <- c(R2_pool, R2_eq)
  Sigma_corr <- stats::cov2cor(as.matrix(Sigmahat))
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
  if (trace) {
    end_fit <- proc.time()[3]
    cat("Time to fit the model: ",
        end_fit-start_fit," seconds\n")
  }  
  ret <- new_spsur(list(call = cl, formula = formula,
                        type = type, data = data, W = W,
                        Durbin = Durbin, 
                        G = G, N = N, Tm = Tm, 
                        deltas = deltas, 
                        deltas.se = deltas.se,  
                        coefficients = coefficients, 
                        rest.se = rest.se,
                        resvar = resvar, 
                        p = pg, 
                        parameters = parameters,
                        BP = BP,
                        R2 = c(R2_pool,R2_eq),
                        Sigma = Sigmahat, 
                        residuals = residuals, 
                        df.residual = df.residual,
                        fitted.values = fitted.values,
                        Y = if(is.null(data)) Y else yf, 
                        X = if(is.null(data)) X else Xf,  
                        zero.policy = zero.policy, 
                        listw_style = listw$style))
  
  if (zero.policy) {
    zero.regs <- attr(listw$neighbours, "region.id")[which(
      spdep::card(listw$neighbours) == 0)]
    if (length(zero.regs) > 0L) 
      attr(ret, "zero.regs") <- zero.regs
  }
  if(exists("na.act")) { # It could not exist with data matrices
    if (!is.null(na.act)) ret$na.action <- na.act
  } 
  ret  
}
