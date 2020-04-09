#' @name spsur3sls
#' @rdname spsur3sls
#'
#' @title Three Stages Least Squares estimation,3sls, of spatial SUR models.
#'
#' @description The function estimates spatial SUR models using three stages
#'  least squares, where the instruments are obtained from the spatial lags
#'  of the \emph{X} variables, assumed to be exogenous. The number of equations, time periods
#'  and spatial units is not restricted. The user can choose between a Spatial Durbin Model
#'  or a Spatial Lag Model, as described below. The estimation procedure allows for the introduction
#'  of linear restrictions on the \eqn{\beta} parameters associated to the regressors.
#'
#' @param type Type of spatial model, restricted to cases where lags of the explainded variable appear
#' in the righ hand side of the equations. There are two possibilities: \strong{"slm"} or
#' \strong{"sdm"}. Default = "slm".
#' @param maxlagW Maximum spatial lag order of the regressors employed to produce spatial
#'  instruments for the spatial lags of the explained variables. Default = 2. Note that in case of
#'  \emph{type="sdm"}, the default value for maxlagW is set to 3 because the first lag of the
#'  regressors, \eqn{WX_{tg}}, can not be used as spatial instruments.
#' @inheritParams spsurml
#'
#' @details
#'  \emph{spsur3sls} can be used to estimate two groups of spatial models:
#'   \itemize{
#'     \item \strong{"slm"}: SUR model with spatial lags of the endogenous in the right hand
#'     side of the equations
#'       \deqn{y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"sdm"}: SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} +
#'              WX_{tg} \theta_{g} + \epsilon_{tg} }
#'              }
#'
#'   where \eqn{y_{tg}} and \eqn{\epsilon_{tg}} are \emph{(Nx1)} vectors,
#'   corresponding to the g-th equation and time period t; \eqn{X_{tg}} is the matrix
#'   of regressors, of order \emph{(Nxp_{g})}. Moreover, \eqn{\lambda_{g}} is a
#'   spatial coefficient and \emph{W} is a \emph{(NxN)} spatial weighting matrix.
#'
#'  By default, the input of this function is an object created with \code{\link[Formula]{Formula}} and
#'  a data frame. However, \emph{spsur3sls} also allows for the direct especification of vector
#'  \emph{Y} and matrix \emph{X}, with the explained variables and  regressors respectively, as
#'  inputs (these terms may be the result, for example, of \code{\link{dgp_spsur}}). \cr
#'
#'
#'  \emph{spsur3sls} is a Least-Squares procedure in three-stages designed to circumvent
#'   the endogeneity problems due to the presence of spatial lags of the explained variable
#'   in the right hand side of the equations do the SUR. The instruments are produced internally
#'   by \emph{spsur3sls} using a sequence of spatial lags of the \emph{X} variables, which are assumed
#'   to be exogenous. The user must define the number of (spatial) instruments to be used in the
#'   procedure, through the argument \emph{maxlagW} (i.e. maxlagW=3). Then, the collection of
#'   instruments generated is \eqn{[WX_{tg}; W*WX_{tg}; W*W*WX_{tg}]}. In the case of a \emph{SDM},
#'   the first lag of the \emph{X} matrix already is in the equation and cannot be used as instrument.
#'   In the example above, the list of instruments for a \emph{SDM} model would be
#'   \eqn{[W^{2}X_{tg}; W^{3}X_{tg}]}.
#'
#'   The \emph{first} stage of the procedure consists in the least squares of the \emph{Y} variables
#'   on the set of instruments. From this estimation, the procedure retains the estimates of \emph{Y}
#'   in the so-called \emph{Yls} variables. In the \emph{second} stage, the \emph{Y} variables that
#'   appear in the right hand side of the equation are substituted by \emph{Yls} and the SUR model
#'   is estimated by Least Squares. The \emph{third} stage improves the estimates of the second stage
#'   through a Feasible Generalized Least Squares estimation of the parameters of the model,
#'   using the residuals of the \emph{second} stage to estimate the \emph{Sigma} matrix.
#'
#'
#'  The arguments \emph{R} and \emph{b} allows to introduce linear restrictions on the \emph{beta}
#'  coefficients of the \emph{G} equations. \code{\link{spsur3sls}}, first, introduces the
#'  linear restrictions in the SUR model and builds, internally, the corresponding constrained
#'  SUR model. Then, the function estimates the restricted model which is shown in the output.
#'  The function does not compute the unconstrained model nor test for the linear restrictions.
#'  The user may ask for the unconstrained estimation using another \code{\link{spsurml}}
#'  estimation. Moreover, the function \code{\link{wald_betas}} obtains the Wald test
#'  of a set of linear restrictions for an object created previously
#'  by \code{\link{spsurml}} or \code{\link{spsur3sls}}.
#'
#' @return
#'  Output of the  three-stages Least-Squares estimation of the specified spatial model.
#'   A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of model specified. \cr
#'     \code{betas} \tab Estimated coefficients for the regressors. \cr
#'     \code{deltas} \tab Estimated spatial coefficients. \cr
#'     \code{se_betas} \tab Estimated standard errors for the estimates of \emph{\eqn{\beta}} coefficients. \cr
#'     \code{se_deltas} \tab Estimated standard errors for the estimates of the spatial coefficients. \cr
#'     \code{cov} \tab Estimated covariance matrix for the estimates of \emph{beta's} and spatial coefficients.\cr
#'     \code{R2} \tab Coefficient of determination for each equation, obtained as the squared
#'      of the correlation coefficient between the corresponding explained variable and
#'       its estimates. \emph{spsur3sls} also shows a \emph{global} coefficient of
#'        determination obtained, in the same manner, for the set of \emph{G} equations. \cr
#'     \code{Sigma} \tab Estimated covariance matrix for the residuals of the \emph{G} equations. \cr
#'     \code{Sigma_corr} \tab Estimated correlation matrix for the residuals of the \emph{G} equations. \cr
#'     \code{Sigma_inv} \tab Inverse of \code{Sigma}, the \emph{(GxG)} covariance matrix of
#'      the residuals of the SUR model. \cr
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom for the residuals. \cr
#'     \code{fitted.values} \tab Estimated values for the dependent variables. \cr
#'     \code{G} \tab Number of equations. \cr
#'     \code{N} \tab Number of cross-sections or spatial units. \cr
#'     \code{Tm} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation (including intercepts). \cr
#'     \code{demean} \tab Logical value used for demeaning. \cr
#'     \code{Y} \tab Vector \emph{Y} of the explained variables of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the regressors of the SUR model. \cr
#'     \code{W} \tab Spatial weighting matrix. \cr
#'  }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#
#'
#' @references
#'   \itemize{
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{wald_betas}}
#'
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (G=1; Tm>1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' ## A SUR model without spatial effects
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#'
#' ## A SUR-SLM model (3SLS Estimation)
#' spcsur.slm.3sls <-spsur3sls(Form = Tformula, data = spc,
#'                             type = "slm", W = Wspc)
#' summary(spcsur.slm.3sls)
#'
#' ## A SUR-SDM model (3SLS Estimation)
#' spcsur.sdm.3sls <-spsur3sls(Form = Tformula, data = spc,
#'                             type = "sdm", W = Wspc)
#' summary(spcsur.sdm.3sls)
#'
#' #################################################
#' ######## PANEL DATA (G>1; Tm>1)         #########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#'
#' ## A SUR-SLR model
#' NCOVRSUR.slm.3sls <-spsur3sls(Form = Tformula, data = NCOVR, type = "slm",
#'                             W = W, maxlagW = 2)
#' summary(NCOVRSUR.slm.3sls)
#'
#' @export
spsur3sls <- function(formula = NULL, data = NULL, na.action, 
                      R = NULL, b = NULL, listw = NULL,
                      quiet = NULL, zero.policy = NULL, 
                      X = NULL, Y = NULL, G = NULL, N = NULL,
                      Tm = NULL, p = NULL, demean = FALSE, 
                      type = "slm",
                      maxlagW = 2) {
  # Función para estimar models SUR-SLM o SUR-SDM espaciales.
  # a través de Mínimos Cuadrados Tres Etapas (3SLS)
  # Spatial Models:  slm, sdm
  
  if (!((type=="slm") || (type=="sdm")))
    stop("3sls can only be used with slm or sdm models")
  
  if (is.null(listw) || 
      !inherits(listw,c("listw","Matrix","matrix")))
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
  
  if (is.null(zero.policy))
    zero.policy <- spatialreg::get.ZeroPolicyOption()
  
  if (!is.null(Tm) && !is.null(G) && Tm > 1 && G == 1){
    #Change dimensions in this case with Matrix Data
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
    if (!is.null(na.act)) {
      subset <- !(1:length(listw$neighbours) %in% na.act)
      listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    W <- Matrix::Matrix(spdep::listw2mat(listw))
    if (type == "sdm") {
      Durbin <- TRUE 
      } else { 
      Durbin <- FALSE 
    }
    get_XY <- get_data_spsur(formula = formula, mf = mf, 
                             Durbin = Durbin,
                             listw = listw, 
                             zero.policy = zero.policy, 
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
  
  #VIP: CAMBIA MATRIZ R SI HAY RESTRICCIONES DADAS POR R Y b
  # REDUCE EL VALOR DE p EN LA ÚLTIMA ECUACIÓN
  if (!is.null(R) & !is.null(b)) {
    demean <- FALSE   # Demeaning isn't allowed in restricted case
    restr <- X_restr(X=X,R=R,b=b,p=p)
    X <- restr$Xstar
    p <- restr$pstar
  }
  if (Tm < 2) demean <- FALSE # Demeaning is not allowed for Tm < 2
  if (demean) {
    demeanXY <- demeaning(X=X,Y=Y,G=G,N=N,Tm=Tm,p=p)
    X <- demeanXY$Xdem
    Y <- demeanXY$Ydem
    p <- demeanXY$pdem
  }
  start_fit <- proc.time()[3]
  # Fits using 3sls
  z <- fit_spsurslm_3sls(Tm = Tm, G = G, N = N, Y = Y,
                         X = X, W = W, p = p,
                         type = type, maxlagW = maxlagW)
  end_fit <- proc.time()[3]
  cat("Time to fit the model: ",
       end_fit-start_fit," seconds\n\n")
  coefficients <- z$coefficients
  deltas <- z$deltas
  Sigma <- z$Sigma
  names_sigma <- NULL
  for (i in 1:G){
    new_name <- paste0("sigma",i,sep="")
    names_sigma <- c(names_sigma,new_name)
  }
  colnames(Sigma) <- rownames(Sigma) <- names_sigma
  parameters <- length(coefficients) + 
    length(deltas) + G*(G + 1)/2
  df.residual <- G*N*Tm - parameters
  dn <- colnames(X)
  if (is.null(dn)) dn <- paste0("x", 1L:(G*sum(p)), sep = "")
  names(coefficients) <- dn
  names_deltas <- NULL
  for (i in 1:G) {
    names_deltas[i] <- paste("rho", i, sep = "_")
  }  
  names(deltas) <- names_deltas
  rest.se <- z$rest.se
  names(rest.se) <- names(coefficients)
  deltas.se <- z$deltas.se
  names(deltas.se) <- names(deltas)
  resvar <- z$resvar
  # Compute R^2 general and for each equation
  Yhat <- z$fitted.values
  R2_pool <- as.numeric((cor(Y,Yhat))^2)
  names(R2_pool) <- c("R2_pool")
  arrYhat <- array(Yhat,c(N,G,Tm))
  arrY <- array(Y,c(N,G,Tm))
  R2_eq <- rep(0,G)
  for(i in 1:G){
    R2_eq[i] <- cor(matrix(arrY[,i,],ncol=1),
                    matrix(arrYhat[,i,],ncol=1))^2
  }
  names(R2_eq) <- paste0("R2_eq",1:G)
  z$R2 <- c(R2_pool,R2_eq)
  ret <- structure(list(call = cl, type = type, 
                        G = G, N = N, Tm = Tm, 
                        deltas = deltas, 
                        deltas.se = deltas.se,  
                        coefficients = coefficients, 
                        rest.se = rest.se,
                        resvar = resvar, 
                        p = p, dvars = dvars,
                        parameters = parameters,
                        R2 = c(R2_pool,R2_eq),
                        Sigma = Sigma, 
                        residuals = z$residuals, 
                        df.residual = df.residual,
                        fitted.values = z$fitted.values, 
                        se.fit = NULL,
                        y = Y, X = X, W = W, 
                        demean = demean,   
                        zero.policy = zero.policy, 
                        listw_style = listw$style), 
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
