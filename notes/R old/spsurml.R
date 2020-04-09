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
#' @param Form An object created with the package \code{\link[Formula]{Formula}}
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
#'       \deqn{y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"sem"}: SUR model with spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \rho_{g} Wu_{tg} + \epsilon_{tg} }
#'     \item \strong{"sdm"}: SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \lambda_{g} Wy_{tg} + X_{tt} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item \strong{"sdem"}: SUR model with spatial lags of the regressors and spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + WX_{tg} \theta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \rho_{g} W u_{tg} + \epsilon_{tg} }
#'     \item \strong{"sarar"}: SUR model with spatial lags of the explained variables and spatial
#'       errors
#'       \deqn{ y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \rho_{g} W u_{tg} + \epsilon_{tg} }
#'   }
#'
#' @return Output of the  maximum-likelihood estimation of the specified spatial SUR model. A list with:
#'   \tabular{ll}{
#'     \code{call} \tab Matched call. \cr
#'     \code{type} \tab  Type of model specified. \cr
#'     \code{betas} \tab Estimated coefficients for the regressors. \cr
#'     \code{deltas} \tab Estimated spatial coefficients. \cr
#'     \code{se_betas} \tab Estimated standard errors for the estimates of \emph{beta}. \cr
#'     \code{se_deltas} \tab Estimated standard errors for the estimates of the spatial coefficients. \cr
#'     \code{cov} \tab Estimated covariance matrix for the estimates of \emph{beta's} and spatial coefficients.\cr
#'     \code{llsur} \tab Value of the likelihood function at the maximum-likelihood estimates. \cr
#'     \code{R2} \tab Coefficient of determination for each equation, obtained as the squared
#'      of the correlation coefficient between the corresponding explained variable and
#'       its estimate. \emph{spsurml} also shows a \emph{global} coefficient of
#'        determination obtained, in the same manner, for the set of \emph{G} equations. \cr
#'     \code{Sigma} \tab Estimated covariance matrix for the residuals of the \emph{G} equations. \cr
#'     \code{Sigma_corr} \tab Estimated correlation matrix for the residuals of the \emph{G} equations. \cr
#'     \code{Sigma_inv} \tab Inverse of \code{Sigma}, the \emph{(GxG)} covariance matrix of
#'      the residuals of the SUR model. \cr
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
#'     \code{Y} \tab Vector \emph{Y} of the explained variables of the SUR model. \cr
#'     \code{X} \tab Matrix \emph{X} of the regressors of the SUR model. \cr
#      \code{W} \tab Spatial weighting matrix. \cr
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
#' ## A SUR model without spatial effects
#' \donttest{
#' ## It usually requires 2-3 minutes maximum...
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' spcsur.sim <-spsurml(Form = Tformula, data = spc, type = "sim", W = Wspc)
#' summary(spcsur.sim)
#' rm(spcsur.sim) # remove
#'
#' ## A SUR-SLX model
#' spcsur.slx <-spsurml(Form = Tformula, data = spc, type = "slx", W = Wspc)
#' summary(spcsur.slx)
#' rm(spcsur.slx) # remove
#'
#' ## A SUR-SLM model
#' spcsur.slm <-spsurml(Form = Tformula, data = spc, type = "slm", W = Wspc)
#' summary(spcsur.slm)
#' rm(spcsur.slm) # remove
#'
#' ## A SUR-SEM model
#' spcsur.sem <-spsurml(Form = Tformula, data = spc, type = "sem", W = Wspc)
#' summary(spcsur.sem)
#' rm(spcsur.sem) # remove
#'
#' ## A SUR-SDM model
#' spcsur.sdm <-spsurml(Form = Tformula, data = spc, type = "sdm", W = Wspc)
#' summary(spcsur.sdm)
#' rm(spcsur.sdm) # remove
#'
#' ## A SUR-SDEM model
#' spcsur.sdem <-spsurml(Form = Tformula, data = spc, type = "sdem", W = Wspc)
#' summary(spcsur.sdem)
#' rm(spcsur.sdem) # remove
#'
#' ## A SUR-SARAR model
#' spcsur.sarar <-spsurml(Form = Tformula, data = spc, type = "sarar", W = Wspc)
#' summary(spcsur.sarar)
#' rm(spcsur.sarar) # remove
#' }
#'
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
#'
#' \donttest{
#' ## It usually requires 1-2 minutes maximum...
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SIM model
#' NCOVRSUR.sim <-spsurml(Form = Tformula, data = NCOVR, type = "sim", W = W)
#' summary(NCOVRSUR.sim)
#' rm(NCOVRSUR.sim)
#'
#' ## A SUR-SLX model
#' NCOVRSUR.slx <-spsurml(Form = Tformula, data = NCOVR, type = "slx", W = W)
#' summary(NCOVRSUR.slx)
#' rm(NCOVRSUR.slx)
#' }
#'
#'\donttest{
#' ## It usually requires 1-2 minutes maximum...
#' ## A SUR-SLM model
#' NCOVRSUR.slm <-spsurml(Form = Tformula, data = NCOVR, type = "slm", W = W)
#' summary(NCOVRSUR.slm)
#' rm(NCOVRSUR.slm)
#'
#' ## A SUR-SDM model
#' NCOVRSUR.sdm <-spsurml(Form = Tformula, data = NCOVR, type = "sdm", W = W)
#' summary(NCOVRSUR.sdm)
#' rm(NCOVRSUR.sdm)
#'
#'
#' ## A SUR-SEM model
#' NCOVRSUR.sem <-spsurml(Form = Tformula, data = NCOVR, type = "sem", W = W)
#' summary(NCOVRSUR.sem)
#' ## A SUR-SDEM model
#' NCOVRSUR.sdem <-spsurml(Form = Tformula, data = NCOVR, type = "sdem",W = W)
#' summary(NCOVRSUR.sdem)
#'
#' ## A SUR-SARAR model
#'  NCOVRSUR.sarar <-spsurml(Form = Tformula, data = NCOVR,
#'                           type = "sarar", W = W)
#'  summary(NCOVRSUR.sarar)
#'}
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
#' psur_sim <- spsurml(Form = pform, data = pNCOVR, W = W, type = "sim")
#' summary(psur_sim)
#'
#' ## SLM (cov = FALSE to prevent overflows of memory)
#'  psur_slm <- spsurml(Form = pform, data = pNCOVR, W = W,
#'                         type = "slm", cov = FALSE)
#'  psur_slm$deltas
#'  psur_slm$betas
#'  psur_slm$Sigma_corr
#'  rm(psur_slm)
#' ## SEM (cov = FALSE to prevent overflows of memory)
#' ### Only execute if you have enough memory...
#'  psur_sem <- spsurml(Form = pform, data = pNCOVR, W = W,
#'                         type = "sem", cov = FALSE)
#'  psur_sem$deltas
#'  psur_sem$betas
#'  psur_sem$Sigma_corr
#'  psur_sem
#' }
#'
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
#' lambda <- 0.5
#' rho <- 0.0 # spatial autocorrelation error term = 0
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
spsurml <- function(Form = NULL, data = NULL, R = NULL,
                    b = NULL, W = NULL,
                    X = NULL,Y = NULL,
                    G = NULL, N = NULL, Tm = NULL,
                    p = NULL, demean = FALSE,
                    type = "sim", cov = TRUE,
                    control = list(tol = 0.05, maxit = 200,
                                   trace = TRUE)) {
  # Spatial Models: sim, slx, slm, sdm, sem, sdem, sarar
  #check for row-standardization of W
  if (!is.null(W)){
    if (class(W) != "matrix") W <- as.matrix(W)
    rsumW <- rowSums(W)
    rsumW[rsumW == 0] <- 1
    nW <- dim(W)[1]
    W <- W / matrix(rep(rsumW, each = nW),
                    nrow = nW, ncol = nW, byrow = TRUE)
    W <- Matrix::Matrix(W)
  }
  if (!is.null(Tm) && !is.null(G) && Tm > 1 && G == 1){
    # Change dimensions in this case with Matrix Data
    G <- Tm
    Tm <- 1
  }

  if(is.null(W) && !type=="sim") stop("W matrix is needed")
  cl <- match.call()
  if(!is.null(Form) && !is.null(data)){
    if (!any(class(Form) == "Formula")) Form <- Formula::Formula(Form)
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
    G <- get_XY$G
    N <- get_XY$N
    Tm <- get_XY$Tm
    p <- get_XY$p

    if (Tm > 1 && G == 1){
       # Change dimensions in this case with Matrix Data
       G <- Tm
       Tm <- 1
     }

    rm(get_XY)
    if (length(p)==1) p <- rep(p,G)
  }
  if (length(p) == 1) p <- rep(p,G)
  names(p) <- NULL
  # CAMBIA MATRIZ X EN EL CASO DURBIN
  if (any(type == c("sdm","sdem","slx"))) {
    IT <- Matrix::Diagonal(Tm)
    IG <- Matrix::Diagonal(G)
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
    p <- p + (p-1)  # Para el caso sdm, slx y sdem se cambia el p
  }
  if (!is.null(R) & !is.null(b)) {
    demean <- FALSE   # Demeaning isn't allowed in restricted case
    restr <- X_restr(X=X,R=R,b=b,p=p)
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

  if (any(type == c("sim","slm","sem","sarar")))
    name_fit <- paste("fit_spsur",type,sep="")
  if(type == "sdm") name_fit <-"fit_spsurslm"
  if(type == "sdem") name_fit <-"fit_spsursem"
  if(type == "slx") name_fit <-"fit_spsursim"
  fit <- get(name_fit)
  tol <- control$tol
  maxit <- control$maxit
  trace <- control$trace
  if (trace) start_fit <- proc.time()[3]
  # Maximize concentrate likelihood
  z <- fit(Tm = Tm, G = G, N = N,
           Y = Y, X = X, W = W, trace = trace,
           tol = tol, maxit = maxit)
  if (trace){
    end_fit <- proc.time()[3]
    cat("Time to fit the model: ",
        end_fit-start_fit," seconds\n")
  }
  z$call <- cl
  if(!is.null(Form) && !is.null(data)){
    z$model <- mf
  }
  z$type <- type
  z$df.residual <- G*N*Tm -
    (length(z$betas) + length(z$deltas) + G*(G+1)/2)
  dn <- colnames(X)
  if (is.null(dn)) dn <- paste0("x", 1L:(G*sum(p)),sep="")
  names(z$betas) <- dn
  if(type=="sarar"){
    names_delta <- rep("",2*G)
  } else names_deltas <- rep("",G)
  names_deltas <- NULL
  for(i in 1:G){
    if(any(type==c("slm","sdm")))
      names_deltas[i] <- paste("lambda",i,sep="_")
    if(any(type==c("sem","sdem")))
      names_deltas[i] <- paste("rho",i,sep="_")
    if(type=="sarar"){
      names_deltas[i] <- paste("lambda",i,sep="_")
      names_deltas[G+i] <- paste("rho",i,sep="_")
    }
  }
  names(z$deltas) <- names_deltas
  if(cov){
    if (any(type == c("sim","slm","sem","sarar")))
      name_cov_fit <- paste("cov_spsur",type,sep="")
    if(type=="sdm") name_cov_fit <-"cov_spsurslm"
    if(type=="sdem") name_cov_fit <-"cov_spsursem"
    if(type=="slx") name_cov_fit <-"cov_spsursim"
    cov_fit <- get(name_cov_fit)
    if (trace) start_cov <- proc.time()[3]
    z_cov <- cov_fit(Tm=Tm,G=G,N=N,
                           Y=Y,X=X,W=W,
                           deltas=z$deltas,Sigma=z$Sigma,trace=trace)
    if (trace){
      end_cov <- proc.time()[3]
      cat("Time to compute covariances: ",
          end_cov-start_cov," seconds \n")
    }
    z$se_betas <- z_cov$se_betas
    names(z$se_betas) <- names(z$betas)
    z$se_deltas <- z_cov$se_deltas
    names(z$se_deltas) <- names(z$deltas)
    z$LMM <- z_cov$LMM
    z$BP <- z_cov$BP
    z$cov <- z_cov$cov
    # CAMBIAR NOMBRES MATRIZ COVARIANZAS
     names_sigma <- NULL
     for (i in 1:G){
       for (j in i:G){
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
  arrYhat <- array(Yhat,c(N,G,Tm))
  arrY <- array(Y,c(N,G,Tm))
  R2_eq <- rep(0,G)
  for(i in 1:G){
    R2_eq[i] <- cor(matrix(arrY[,i,],ncol=1),
                    matrix(arrYhat[,i,],ncol=1))^2
  }
  names(R2_eq) <- paste0("R2_eq",1:G)
  z$R2 <- c(R2_pool,R2_eq)
  z$p <- p
  z$G <- G
  z$N <- N
  z$Tm <- Tm
  z$Y <- Y
  z$X <- X
  z$W <- W
  z$demean <- demean
  class(z) <- c("spsur")
  z
}

