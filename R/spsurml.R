#' @name spsurml
#' @rdname spsurml
#'
#' @title Maximum likelihood estimation of spatial SUR model
#'
#' @description
#' Maximum likelihood estimation of the SUR model without spatial effects
#'   SUR-SIM, and  the spatial SUR-SLX; SUR-SAR; SUR-SEM; SUR-SDM; SUR-SEDM
#'   and SUR-SARAR models (including possible restrictions of beta coefficients).
#'
#' @param    Form   : An object create with \code{\link[Formula]{Formula}} package allowing for multiple responses and multiple parts of regressors
#' @param    data   : An object of class data.frame or a matrix
#' @param    W      : A RxR spatial weight matrix
#' @param    Y      : Default NULL. Data vector nRxnTx1 (firs: space dimension | second: time periods)
#' @param    X      : Default NULL. Data matrix nRxnTxp (p=sum(p_g) where p_g is the number of independent variables for g-th equation, g=1,...,nG)
#' @param    nG     : number of equations
#' @param    nR     : number of spatial observations
#' @param    nT     : Number of time periods
#' @param    p      : Number of regressors by equation (including independent term). A number (resp. vector) if equal (resp. no equal) regressors by equation
#' @param    type   : type of model "sim" "slx" "sar" "sem" "sdm" "sdem" "sarar" . Defalult is "sim"
#' @param    R      : Default NULL. A matrix including coefficients of linear hypothesis on beta
#' @param    r      : Default NULL. A vector including independent vector of linear hypothesis on beta
#' @param    demean : Default FALSE.
#' @param    cov    : get the covariance matrix (TRUE or FALSE). Defalult is TRUE
#' @param    trace  : To get intermediate results (TRUE or FALSE). Defalult is TRUE
#'
#' @details
#' Maximum likelihood estimation of the models: \cr
#' \cr
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
#' where y_{gt}, u_{gt} and \eqn{\epsilon_{gt}} are (nRx1) vectors; X_{gt} is a matrix of exogenous variables of
#' order (nRxsum(p)); \eqn{\lambda_g} and \eqn{\rho_g} are parametres of spatial dependence; W is a nRxnR matrix of spatial interactions.
#'
#' By default an object create with \code{\link[Formula]{Formula}} is the input of this function. Alternatively a Y vector of dependent
#' variables and X matrix with independent varaibles could be include as imputs obtain by exaple with \code{\link{dgp_spSUR}}\cr
#' \cr
#' The marginal Multiplier tests are used to test for no correlation in one part of the model allowing for spatial correlation in the other. \cr
#' The LM(LM(\eqn{\rho}|\eqn{\lambda})) is the test for sustantive spatial autocorrelation in a model with spatial autocorrelation in error term.
#' The LM(\eqn{\lambda}|\eqn{\rho}) is the test for spatial error correlation in a model with sustanative spatial correlation.
#'
#' @return
#' Results of ML estimation of a SUR model with spatial effects. A list with:
#'   \tabular{ll}{
#'   \code{betas} \tab Coefficients estimated by maximum likelihood. \cr
#'   \code{llsur}    \tab the maximum value of likelihood function. \cr
#'   \code{Sigma}    \tab estimated covariance matrix of residuals between equations. \cr
#'   \code{Sigma_corr}    \tab estimated correlation matrix of residuals between equations. \cr
#'   \code{residuals}    \tab residual of model. \cr
#'   \code{fitted.values}    \tab fitted values of model. \cr
#'   \code{call}   \tab the model estimate. \cr
#'   \code{type}   \tab the matched call. \cr
#'   \code{se_betas}   \tab Standard desviation of beta coeeffients. \cr
#'   \code{BP}   \tab Breusch-Pagan test of diagonality of the SUR matrix of covariance residuals. \cr
#'   \code{R2}    \tab Determination coefficient between the explained varaible and its predictionsthe p-value of coefficientes. Global and for each equation. \cr
#'   \code{BP}    \tab Breush-Pagan test (1980) of diagonaltiy of covariance matrix. \cr
#'   \code{LMM}    \tab Marginal test (LM(\eqn{\rho}|\eqn{\lambda}) or LM(\eqn{\lambda}|\eqn{\rho})) of spatial autocorrelation in the residuals. \cr
#'   }
#' @references
#'
#' Mur, J., Lopez, F., & Herrera, M. (2010). Testing for spatial effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#' \cr
#' \cr
#' Lopez, F.A., Mur, J., & Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' \cr
#' \cr
#' Breusch T, Pagan A (1980) The Lagrange multiplier test and its applications to model specification in
#'  econometrics. \emph{Rev Econ Stud} 47: 239-254
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@unizar.es} \cr
#'   }
#' @seealso
#'
#' \code{\link{lmtestspsur}}, \code{\link{spsurml}}, \code{\link{spsur3sls}}
#'
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
#' spcSUR.sim <-spsurml(Form=Tformula,data=spc,type="sim",W=Wspc)
#' summary(spcSUR.sim)
#'
#' ## A SUR-SLX model
#' spcSUR.slx <-spsurml(Form=Tformula,data=spc,type="slx",W=Wspc)
#' summary(spcSUR.slx)
#'
#' ## A SUR-SAR model
#' spcSUR.sar <-spsurml(Form=Tformula,data=spc,type="sar",W=Wspc)
#' summary(spcSUR.sar)
#'
#' ## A SUR-SEM model
#' spcSUR.sem <-spsurml(Form=Tformula,data=spc,type="sem",W=Wspc)
#' summary(spcSUR.sem)
#'
#' ## A SUR-SDM model
#' spcSUR.sdm <-spsurml(Form=Tformula,data=spc,type="sdm",W=Wspc)
#' summary(spcSUR.sdm)
#'
#' ## A SUR-SDEM model
#' spcSUR.sdem <-spsurml(Form=Tformula,data=spc,type="sdem",W=Wspc)
#' summary(spcSUR.sdem)
#'
#' ## A SUR-SARAR model
#' spcSUR.sarar <-spsurml(Form=Tformula,data=spc,type="sarar",W=Wspc)
#' summary(spcSUR.sarar)
#'
#' #################################################
#' ######## PANEL DATA (nG>1; nT>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' Homicides and selected socio-economic characteristics for continental U.S. counties. Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' https://geodacenter.github.io/data-and-lab/ncovr/
#'
#' data(NAT)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SIM model
#' NATSUR.sim <-spsurml(Form=Tformula,data=NAT,type="sim",W=W)
#' summary(NATSUR.sim)
#'
#' ## A SUR-SLX model
#' NATSUR.slx <-spsurml(Form=Tformula,data=NAT,type="slx",W=W)
#' summary(NATSUR.slx)
#'
#' ## A SUR-SAR model
#' NATSUR.sar <-spsurml(Form=Tformula,data=NAT,type="sar",W=W)
#' summary(NATSUR.sar)
#'
#' ## A SUR-SEM model
#' NATSUR.sem <-spsurml(Form=Tformula,data=NAT,type="sem",W=W)
#' summary(NATSUR.sem)
#'
#' ## A SUR-SDM model
#' NATSUR.sdm <-spsurml(Form=Tformula,data=NAT,type="sdm",W=W)
#' summary(NATSUR.sdm)
#'
#' ## A SUR-SDEM model
#' NATSUR.sdem <-spsurml(Form=Tformula,data=NAT,type="sdem",W=W)
#' summary(NATSUR.sdem)
#'
#' #' ## A SUR-SARAR model
#' NATSUR.sarar <-spsurml(Form=Tformula,data=NAT,type="sarar",W=W)
#' summary(NATSUR.sarar)
#'
#'
#' ##############################################
#' ######## PANEL DATA (nG>1; nT>1) DEMEANING ###
#' ##############################################
#'
#' #### Example 3: Matrix as input argument
#' data(Sar)
#' nT <- 4 # Number of periods
#' nG <- 3 # Number equations
#' nR <- 49 # Number of spatial units
#' SUR.sar.ml <-spsurml(Y=Ysar,X=XXsar,nG=nG,nR=nR,nT=nT,p=5,W=Ws,type="sar")
#' summary(SUR.sar.ml)
#'
#' #### Estimation with demeaning in nT (intercept of each equation dissapears)
#' SUR.sar.ml.dem <-spsurml(Y=Ysar,X=XXsar,nG=nG,nR=nR,nT=nT,p=5,W=Ws,type="sar",demean=TRUE)
#' summary(SUR.sar.ml.dem)
#'
#' # Durbin case with demeaning in nT
#' SUR.sdm.ml.dem <-spsurml(Y=Ysar,X=XXsar,nG=nG,nR=nR,nT=nT,p=5,W=Ws,type="sdm",demean=TRUE)
#' summary(SUR.sdm.ml.dem)

spsurml <- function(Form=NULL,data=NULL,R=NULL,r=NULL,W=NULL,
                  X=NULL,Y=NULL,
                  nG=NULL,nR=NULL,nT=NULL,p=NULL,demean=FALSE,
                  type="sim",
                  cov=TRUE,trace=TRUE)
{
  # Función para estimar cualquier modelo SUR espacial.
  # imponiendo restricciones
  # Spatial Models: sim, slx, sar, sdm,
  #                 sem, sdem, sarar
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
    if(!is.null(Form) & !any(class(Form)=="Formula")) Form <- Formula::Formula(Form)
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

  ##p <- diff(c(which(colSums(X)/nT==nR), length(colSums(X))))
  names(p) <- NULL
  # CAMBIA MATRIZ X EN EL CASO DURBIN
  if (any(type == c("sdm","sdem","slx"))) {
    IT <- Matrix::Diagonal(nT)
    IG <- Matrix::Diagonal(nG)

    # index_interc_X <- c(1,1+cumsum(p[1:(length(p)-1)]))
    # X_nointerc <- X[,-c(index_interc_X)]
    # WX <- (IT %x% IG %x% W) %*% X_nointerc
    # colnames(WX) <- paste0("W_",colnames(X_nointerc))
    # Xdurbin <- as.matrix(cbind(X,WX))
    # Construye Matriz uniendo X y WX (sin intercepto) para cada ecuación
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
    p <- p + (p-1)  # Para el caso sdm, slx y sdem se cambia el p (ojo Intercepto)
  }
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

  if (any(type == c("sim","sar","sem","sarar")))
    name_fit <- paste("fit_spsur",type,sep="")
  if(type=="sdm") name_fit <-"fit_spsursar"
  if(type=="sdem") name_fit <-"fit_spsursem"
  if(type=="slx") name_fit <-"fit_spsursim"
  fit <- get(name_fit)
  if (trace) start_fit <- proc.time()[3]
  # Maximize concentrate likelihood
  z <- fit(nT=nT,nG=nG,nR=nR,
                 Y=Y,X=X,W=W,trace=trace)
  if (trace){
    end_fit <- proc.time()[3]
    cat("Time to fit the model: ",
        end_fit-start_fit," seconds\n\n")
  }
  z$call <- cl
  if(!is.null(Form) && !is.null(data)){
    #z$terms <- mt
    z$model <- mf
  }
  z$type <- type
  z$df.residual <- nG*nR*nT -
    (length(z$betas) + length(z$deltas) + nG*(nG+1)/2)
  dn <- colnames(X)
  if (is.null(dn)) dn <- paste0("x", 1L:(nG*sum(p)),sep="")
  names(z$betas) <- dn
  if(type=="sarar"){
    names_delta <- rep("",2*nG)
  } else names_deltas <- rep("",nG)
  names_deltas <- NULL
  for(i in 1:nG){
    if(any(type==c("sar","sdm")))
      names_deltas[i] <- paste("lambda",i,sep="_")
    if(any(type==c("sem","sdem")))
      names_deltas[i] <- paste("rho",i,sep="_")
    if(type=="sarar"){
      names_deltas[i] <- paste("lambda",i,sep="_")
      names_deltas[nG+i] <- paste("rho",i,sep="_")
    }
  }
  names(z$deltas) <- names_deltas
  if(cov==TRUE){
    if (any(type == c("sim","sar","sem","sarar")))
      name_cov_fit <- paste("cov_spsur",type,sep="")
    if(type=="sdm") name_cov_fit <-"cov_spsursar"
    if(type=="sdem") name_cov_fit <-"cov_spsursem"
    if(type=="slx") name_cov_fit <-"cov_spsursim"
    cov_fit <- get(name_cov_fit)
    if (trace) start_cov <- proc.time()[3]
    z_cov <- cov_fit(nT=nT,nG=nG,nR=nR,
                           Y=Y,X=X,W=W,
                           deltas=z$deltas,Sigma=z$Sigma)
    if (trace){
      end_cov <- proc.time()[3]
      cat("Time to compute covariances: ",
          end_cov-start_cov," seconds \n\n")
    }
    z$se_betas <- z_cov$se_betas
    z$se_deltas <- z_cov$se_deltas
    z$LMM <- z_cov$LMM
    z$BP <- z_cov$BP
    z$cov <- z_cov$cov
    # CAMBIAR NOMBRES MATRIZ COVARIANZAS
     names_sigma <- NULL
     for (i in 1:nG){
       for (j in i:nG){
         new_name <- paste0("sigma",i,j,sep="")
         names_sigma <- c(names_sigma,new_name)
       }
     }
     if(type=="sim" | type=="slx"){
       colnames(z$cov) <- names(z$betas)
       rownames(z$cov) <- names(z$betas)
     } else {
       colnames(z$cov) <- c(names(z$betas),names(z$deltas),names_sigma)
       rownames(z$cov) <- c(names(z$betas),names(z$deltas),names_sigma)
     }
  }
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

