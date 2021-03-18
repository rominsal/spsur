#' @name spsur3sls
#' @rdname spsur3sls
#'
#' @title Three Stages Least Squares estimation,3sls, of spatial SUR models.
#'
#' @description The function estimates spatial SUR models using three stages
#'  least squares, where the instruments are obtained from the spatial lags
#'  of the \emph{X} variables, assumed to be exogenous. The number of 
#'  equations, time periods and spatial units is not restricted. The user can 
#'  choose between a Spatial Durbin Model or a Spatial Lag Model, as described 
#'  below. The estimation procedure allows for the introduction of linear 
#'  restrictions on the \eqn{\beta} parameters associated to the regressors.
#'  
#' @usage spsur3sls (formula = NULL, data = NULL, na.action,
#'                   R = NULL, b = NULL, listw = NULL, 
#'                   zero.policy = NULL, X= NULL, Y = NULL, G = NULL, 
#'                   N = NULL, Tm = NULL, p = NULL,  
#'                   type = "slm", Durbin = NULL, maxlagW = NULL,
#'                   trace = TRUE)
#'     
#' @param type Type of spatial model, restricted to cases where lags of the 
#' explained variable appear in the rigth hand side of the equations. There 
#' are two possibilities: "slm" or "sdm". Default = "slm".
#' @param maxlagW Maximum spatial lag order of the regressors employed to 
#' produce spatial instruments for the spatial lags of the explained variables. 
#' Default = 2. Note that in case of \code{type}="sdm", the default value for 
#' \code{maxlagW} is set to 3 because the first lag of the regressors, 
#' \eqn{WX_{tg}}, can not be used as spatial instruments.
#' @param Durbin If a formula object and model is type "sdm" the subset 
#'   of explanatory variables to lag for each equation.   
#' @param trace A logical value to show intermediate results during 
#'       the estimation process. Default = \code{TRUE}. 
#' @inheritParams spsurml
#'
#' @details
#'  \emph{spsur3sls} can be used to estimate two groups of spatial models:
#'   \itemize{
#'     \item "slm": SUR model with spatial lags of the endogenous in 
#'     the right hand side of the equations
#'       \deqn{y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item "sdm": SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} +
#'              WX_{tg} \theta_{g} + \epsilon_{tg} }
#'              }
#'
#'   where \eqn{y_{tg}} and \eqn{\epsilon_{tg}} are \emph{(Nx1)} vectors,
#'   corresponding to the g-th equation and time period t; \eqn{X_{tg}} is 
#'   the matrix of regressors, of order \emph{(Nxp_{g})}. Moreover, 
#'   \eqn{\rho_{g}} is a spatial coefficient and \emph{W} is a 
#'   \emph{(NxN)} spatial weighting matrix.
#'
#'  By default, the input of this function is an object created with 
#'  \code{\link[Formula]{Formula}} and a data frame. However, 
#'  \emph{spsur3sls} also allows for the direct specification of vector
#'  \emph{Y} and matrix \emph{X}, with the explained variables and  regressors 
#'  respectively, as inputs (these terms may be the result, for example, 
#'  of \code{\link{dgp_spsur}}). \cr
#'
#'
#'  \code{spsur3sls} is a Least-Squares procedure in three-stages designed 
#'  to circumvent the endogeneity problems due to the presence of spatial 
#'  lags of the explained variable in the right hand side of the equations 
#'  do the SUR. The instruments are produced internally by \code{spsur3sls} 
#'  using a sequence of spatial lags of the \emph{X} variables, which are 
#'  assumed to be exogenous. The user must define the number of (spatial) 
#'  instruments to be used in the procedure, through the argument 
#'  \code{maxlagW} (i.e. maxlagW = 3). Then, the collection of instruments 
#'  generated is \eqn{[WX_{tg}; W*WX_{tg}; W*W*WX_{tg}]}. In the case of 
#'  a \emph{SDM}, the first lag of the \emph{X} matrix already is in the 
#'  equation and cannot be used as instrument. In the example above, the 
#'  list of instruments for a \emph{SDM} model would be 
#'  \eqn{[W^{2}X_{tg}; W^{3}X_{tg}]}.
#'
#'   The \emph{first} stage of the procedure consists in the least squares 
#'   of the \emph{Y} variables on the set of instruments. From this 
#'   estimation, the procedure retains the estimates of \emph{Y}
#'   in the so-called \emph{Yls} variables. In the \emph{second} stage, 
#'   the \emph{Y} variables that appear in the right hand side of the equation 
#'   are substituted by \emph{Yls} and the SUR model
#'   is estimated by Least Squares. The \emph{third} stage improves the 
#'   estimates of the second stage through a Feasible Generalized Least Squares 
#'   estimation of the parameters of the model,
#'   using the residuals of the \emph{second} stage to estimate the 
#'   \emph{Sigma} matrix.
#'
#'
#'  The arguments \emph{R} and \emph{b} allows to introduce linear 
#'  restrictions on the \emph{beta} coefficients of the \emph{G} equations. 
#'  \code{\link{spsur3sls}}, first, introduces the linear restrictions in 
#'  the SUR model and builds, internally, the corresponding constrained
#'  SUR model. Then, the function estimates the restricted model which is 
#'  shown in the output. The function does not compute the unconstrained 
#'  model nor test for the linear restrictions. The user may ask for the 
#'  unconstrained estimation using another \code{\link{spsurml}}
#'  estimation. Moreover, the function \code{\link{wald_betas}} obtains 
#'  the Wald test of a set of linear restrictions for an object created 
#'  previously by \code{\link{spsurml}} or \code{\link{spsur3sls}}.
#'
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
#
#'
#' @references
#'   \itemize{
#'       \item López, F. A., Mínguez, R., Mur, J. (2020). ML versus IV estimates 
#'       of spatial SUR models: evidence from the case of Airbnb in Madrid urban 
#'       area. \emph{The Annals of Regional Science}, 64(2), 313-347.
#'       \url{https://doi.org/10.1007/s00168-019-00914-1}
#'       \item Anselin, L. (2016) Estimation and Testing in the Spatial Seemingly 
#'       Unrelated Regression (SUR). \emph{Geoda Center for Geospatial Analysis 
#'       and Computation, Arizona State University}. Working Paper 2016-01.
#'       \url{http://dx.doi.org/10.13140/RG.2.2.15925.40163}
#'       \item, Anselin, L. (1988). \emph{Spatial Econometrics: Methods and Models}. 
#'       Kluwer Academic Publishers, Dordrecht, The Netherlands (p. 146).
#'       \item Anselin, L., Le Gallo, J., Hubert J. (2008) Spatial Panel Econometrics. 
#'       In \emph{The econometrics of panel data. 
#'       Fundamentals and recent developments in theory and practice}. (Chap 19, p. 653)
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link[spatialreg]{stsls}}, 
#' \code{\link{wald_betas}}
#'
#' @examples
#'
#' #################################################
#' ######## CLASSIC PANEL DATA (G=1; Tm>1)  ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' ## A SUR model without spatial effects
#' rm(list = ls()) # Clean memory
#' data(spc)
#' lwspc <- spdep::mat2listw(Wspc, style = "W")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#'
#' ## A SUR-SLM model (3SLS Estimation)
#' spcsur_slm_3sls <-spsur3sls(formula = Tformula, data = spc,
#'                             type = "slm", listw = lwspc)
#' summary(spcsur_slm_3sls)
#' print(spcsur_slm_3sls)
#' 
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur_slm_3sls, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' 
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_spsur3sls, package="spsur")
#' 
#' \donttest{
#' ## A SUR-SDM model (3SLS Estimation)
#' spcsur_sdm_3sls <- spsur3sls(formula = Tformula, data = spc,
#'                              type = "sdm", listw = lwspc)
#' summary(spcsur_sdm_3sls)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur_sdm_3sls, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' rm(spcsur_sdm_3sls)
#' 
#' ## A SUR-SDM model with different spatial lags in each equation
#'  TformulaD <-  ~ UN83 + NMR83 + SMSA | UN80 + NMR80  
#'  spcsur_sdm2_3sls <-spsur3sls(formula = Tformula, data = spc,
#'                              type = "sdm", listw = lwspc,
#'                              Durbin = TformulaD)
#'  summary(spcsur_sdm2_3sls)
#' if (require(gridExtra)) {
#'   pl <- plot(spcsur_sdm2_3sls, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' rm(spcsur_sdm2_3sls)
#' }
#'
#' #################################################
#' ###  MULTI-DIMENSIONAL PANEL DATA (G>1; Tm>1) ###
#' #################################################
#'
#' #### Example 3: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' 
#' \donttest{
#' rm(list = ls()) # Clean memory
#' data(NCOVR, package = "spsur")
#' nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ## Some regions with no links...
#' lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SLM model
#' NCOVRSUR_slm_3sls <- spsur3sls(formula = Tformula, data = NCOVR.sf, 
#'                                type = "slm", zero.policy = TRUE,
#'                                listw = lwncovr, trace = FALSE)
#' summary(NCOVRSUR_slm_3sls)
#' if (require(gridExtra)) {
#'   pl <- plot(NCOVRSUR_slm_3sls, viewplot = FALSE) 
#'   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
#'                pl$pldeltas, nrow = 3)
#' }
#' rm(NCOVRSUR_slm_3sls)
#' }
#' @export
spsur3sls <- function(formula = NULL, data = NULL, na.action, 
                      R = NULL, b = NULL, listw = NULL,
                      zero.policy = NULL, 
                      X = NULL, Y = NULL, G = NULL, N = NULL,
                      Tm = NULL, p = NULL,  
                      type = "slm", Durbin = NULL,
                      maxlagW = NULL, trace = TRUE) {
  # Función para estimar models SUR-SLM o SUR-SDM espaciales.
  # a través de Mínimos Cuadrados Tres Etapas (3SLS)
  # Spatial Models:  slm, sdm
  if (!((type == "slm") || (type == "sdm")))
    stop("3sls can only be used with slm or sdm models")
  if (is.null(maxlagW)) {
    if (type == "slm") {
      maxlagW <- 2 
    } else {
      maxlagW <- 3
    }
  }  
  
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
    # Change dimensions in this case with Matrix Data
    G <- Tm
    Tm <- 1
  }
  if (!is.null(formula) && !inherits(formula, "Formula")) 
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
    Xorig <- X
    porig <- p
    restr <- X_restr(X = X, R = R, b = b, p = p)
    X <- restr$Xstar
    p <- restr$pstar
  }
  start_fit <- proc.time()[3]
  # Fits using 3sls
  z <- fit_spsurslm_3sls(Tm = Tm, G = G, N = N, Y = Y,
                         X = X, W = W, p = p,
                         type = type, maxlagW = maxlagW)
  end_fit <- proc.time()[3]
  if (trace) 
    cat("Time to fit the model: ", end_fit-start_fit," seconds\n\n")
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
  for (i in 1:G) {
    R2_eq[i] <- cor(matrix(arrY[,i,],ncol=1),
                    matrix(arrYhat[,i,],ncol=1))^2
  }
  names(R2_eq) <- paste0("R2_eq",1:G)
  z$R2 <- c(R2_pool,R2_eq)
  # Añadido por Fernando: Hay que incluir BP y LMM con las covarianzas numéricas
  Sigma_corr <- stats::cov2cor(Sigma)
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
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
                        Durbin = Durbin, 
                        G = G, N = N, Tm = Tm, 
                        deltas = deltas, 
                        deltas.se = deltas.se,  
                        coefficients = coefficients, 
                        rest.se = rest.se,
                        resvar = resvar, 
                        p = p, 
                        dvars = if(exists("dvars")) dvars else NULL,
                        parameters = parameters,
                        BP=BP,
                        R2 = c(R2_pool,R2_eq),
                        Sigma = Sigma, 
                        residuals = z$residuals, 
                        df.residual = df.residual,
                        fitted.values = z$fitted.values, 
                        Y = Y, X = X, W = W, 
                        zero.policy = zero.policy, 
                        listw_style = listw$style, 
                        maxlagW = maxlagW))
  
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
