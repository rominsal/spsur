#' @name spsur3sls
#' @rdname spsur3sls
#'
#' @title 3sls estimation of spatial SUR model including restrictions
#'   of beta coefficients.
#'
#' @description The function fits by three-stage least squares
#'   estimation of spatial a spatial SUR-lag spatial or SUR-SDM model.
#'
#' @param type Type of chosen model between "sar" or "sdm". Default = "sar".
#' @param maxlagW Maximum order of the instruments. Default = 2.
#' @inheritParams spsurml
#'
#' @details  Three-stage least squate (3sls) estimation of the SUR
#'   models by using spatially lagged X variables as instruments for
#'   the spatially lagged dependent variable: \cr \cr
#'   SUR-SAR: Spatial autorregresive model \cr
#'  \deqn{ y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} +
#'    \epsilon_{gt} }
#'   SUR-SDM: Spatial Durbin model \cr
#'   \deqn{ y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} + WX_{gt}
#'     \theta_{g} + \epsilon_{gt} }
#'   where \eqn{y_{gt}}, \eqn{u_{gt}} and \eqn{\epsilon_{gt}}
#'   are (\emph{nR}x1) vectors; \emph{X_{gt}} is a matrix of exogenous
#'   variables of order (\emph{nR}x\code{sum(p)});
#'   \eqn{\lambda_g} and \eqn{\rho_g} are parametres of spatial
#'   dependence; \emph{W} is a \emph{nR}x\emph{nR} matrix of spatial
#'   interactions.
#'
#'   The instruments used are \emph{WX} and \eqn{W^{2}X} by default. \cr
#'   By default an object create with \code{\link[Formula]{Formula}}
#'   is the input of this function. Alternatively a \emph{Y} vector of
#'   dependent variables and \emph{X} matrix with independent
#'   variables could be include as inputs obtained by example
#'   with \code{\link{dgp_spSUR}} \cr \cr
#'
#'
#' @return Results of 3sls estimation of a SUR model with spatial effects.
#'   A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of specified model. \cr
#'     \code{betas} \tab Fitted coefficients for the covariates. \cr
#'     \code{deltas} \tab Fitted spatial autocorrelation coefficients. \cr
#'     \code{se_betas} \tab Standard errors of beta. \cr
#'     \code{se_deltas} \tab Standard errors of delta. \cr
#'     \code{cov} \tab Fitted covariance matrix of betas and deltas.\cr
#'     \code{R2} \tab Determination coefficient between the explained
#'       variable and its predictions. Global and for each equation. \cr
#'     \code{Sigma} \tab Fitted covariance matrix of residuals between
#'       equations. \cr
#'     \code{Sigma_corr} \tab Fitted correlation matrix of residuals
#'       between equations. \cr
#'     \code{residuals} \tab Residuals of the model. \cr
#'     \code{df.residuals} \tab Degrees of freedom of the residuals. \cr
#'     \code{fitted.values} \tab Fitted values of the dependent variables. \cr
#'     \code{nG} \tab Number of equations. \cr
#'     \code{nR} \tab Number of cross-section or spatial observations. \cr
#'     \code{nT} \tab Number of time periods. \cr
#'     \code{p} \tab Number of regressors by equation
#'       (including independent terms). \cr
#'     \code{demean} \tab Logical value used for demeaning. \cr
#'     \code{Y} \tab Vector \emph{Y} of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the SUR model. \cr
#'     \code{W} \tab Spatial weight matrix. \cr
#'     \code{Sigma_inv} \tab Inverse of \code{Sigma}. \cr
#'  }
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
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @seealso
#' \code{\link{lmtestspsur}}, \code{\link{spsurml}}
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (nG>1; nT=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' ## A SUR model without spatial effects
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#'
#' ## A SUR-SAR model (3SLS Estimation)
#' spcSUR.sar.3sls <-spsur3sls(Form = Tformula, data = spc,
#'                             type = "sar", W = Wspc)
#' summary(spcSUR.sar.3sls)
#'
#' ## A SUR-SDM model (3SLS Estimation)
#' spcSUR.sdm.3sls <-spsur3sls(Form = Tformula, data = spc,
#'                             type = "sdm", W = Wspc)
#' summary(spcSUR.sdm.3sls)
#'
#' #################################################
#' ######## PANEL DATA (nG>1; nT>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for continental
#' # U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' data(NAT)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#'
#' ## A SUR-SAR model
#' NATSUR.sar.3sls <-spsur3sls(Form = Tformula, data = NAT, type = "sar",
#'                             W = W, maxlagW = 2)
#' summary(NATSUR.sar.3sls)
#' @export
spsur3sls <- function(Form = NULL, data = NULL, R = NULL, r = NULL,
                      W = NULL, X = NULL, Y = NULL, nG = NULL, nR = NULL,
                      nT = NULL, p = NULL, demean = FALSE, type = "sar",
                      maxlagW = 2) {

  # Función para estimar models SUR-SAR o SUR-SDM espaciales.
  # a través de Mínimos Cuadrados Tres Etapas (3SLS)
  # Spatial Models:  sar, sdm

  if (!((type=="sar") || (type=="sdm")))
  {
    stop("3sls can only be used with sar or sdm models")
  }
  if(is.null(W) && !type=="sim") stop("W matrix is needed")
  if(!is.null(W)) W <- Matrix::Matrix(W)
  cl <- match.call()
  if(!is.null(Form) && !is.null(data)){
    if (!class(Form) == "Formula") Form <- Formula::Formula(Form)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("Formula", "data", "subset",
                 "weights", "na.action",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    # Obtener Datos
    if(!is.null(Form) & !any(class(Form)=="Formula")) {
      Form <- Formula::Formula(Form)
    }
    get_XY <- get_data_spsur(formula=Form,data=data,W=W)
    Y <- get_XY$Y
    X <- get_XY$X
    nG <- get_XY$nG
    nR <- get_XY$nR
    nT <- get_XY$nT
    p <- get_XY$p
    rm(get_XY)
    if (length(p)==1) p <- rep(p,nG)
  }
  if (length(p)==1) p <- rep(p,nG)
  names(p) <- NULL
  # CAMBIA MATRIZ X EN EL CASO DURBIN
  if (type == "sdm") {
    IT <- Matrix::Diagonal(nT)
    IG <- Matrix::Diagonal(nG)
    WX <- (IT %x% IG %x% W) %*% X
    colnames(WX) <- paste0("W_",colnames(X))
    Xdurbin <- NULL
    for (i in 1:length(p))
    {
      if(i==1){
        Xdurbin <- cbind(X[,1:p[i]],WX[,2:p[i]])
      } else {
        Xdurbin <- cbind(Xdurbin,
                         X[,(cumsum(p)[i-1]+1):cumsum(p)[i]],
                         WX[,(cumsum(p)[i-1]+2):cumsum(p)[i]]) # Sin intercepto
      }
    }
    X <- as.matrix(Xdurbin); rm(Xdurbin)
    p <- p + (p-1)  # Para el caso sdm cambia el p (ojo Intercepto)
  }
  #VIP: CAMBIA MATRIZ R SI HAY RESTRICCIONES DADAS POR R Y r
  # REDUCE EL VALOR DE p EN LA ÚLTIMA ECUACIÓN
  if (!is.null(R) & !is.null(r)) {
    demean <- FALSE   # Demeaning isn't allowed in restricted case
    restr <- X_restr(X=X,R=R,r=r,p=p)
    X <- restr$Xstar
    p <- restr$pstar
  }
  if (nT < 2) demean <- FALSE # Demeaning is not allowed for nT < 2
  if (demean) {
    demeanXY <- demeaning(X=X,Y=Y,nG=nG,nR=nR,nT=nT,p=p)
    X <- demeanXY$Xdem
    Y <- demeanXY$Ydem
    p <- demeanXY$pdem
  }
  start_fit <- proc.time()[3]
  # Fits using 3sls
  z <- fit_spsursar_3sls(nT=nT,nG=nG,nR=nR,Y=Y,X=X,W=W,p=p,
            type=type,maxlagW=maxlagW,trace=trace)
  end_fit <- proc.time()[3]
  cat("Time to fit the model: ",
       end_fit-start_fit," seconds\n\n")
  z$call <- cl
  if(!is.null(Form) && !is.null(data)){
    z$model <- mf
  }
  z$type <- type
  z$df.residual <- nG*nR*nT -
    (length(z$betas) + length(z$deltas) + nG*(nG+1)/2)
  dn <- colnames(X)
  if (is.null(dn)) dn <- paste0("x", 1L:(nG*sum(p)),sep="")
  names(z$betas) <- dn
  z$LMM <- NULL
  z$BP <- NULL
  z$llsur <- NULL
  # Compute R^2 general and for each equation
  Yhat <- z$fitted.values
  R2_pool <- as.numeric((cor(Y,Yhat))^2)
  names(R2_pool) <- c("R2_pool")
  arrYhat <- array(Yhat,c(nR,nG,nT))
  arrY <- array(Y,c(nR,nG,nT))
  R2_eq <- rep(0,nG)
  for(i in 1:nG){
    R2_eq[i] <- cor(matrix(arrY[,i,],ncol=1),
                    matrix(arrYhat[,i,],ncol=1))^2
  }
  names(R2_eq) <- paste0("R2_eq",1:nG)
  z$R2 <- c(R2_pool,R2_eq)
  z$p <- p
  z$nG <- nG
  z$nR <- nR
  z$nT <- nT
  z$Y <- Y
  z$X <- X
  z$W <- W
  z$demean <- demean
  class(z) <- c("spsur")
  z
}
