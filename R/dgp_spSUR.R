#' @name dgp_spSUR
#' @rdname dgp_spSUR
#'
#' @title Data Generating Process spSUR
#'
#' @description
#' Data Generating Process of a spatial SUR. nRxnGxnT with \deqn{p_{g}} (g=1,...,nG) independent variables (all N(0,1))
#' with constant term and equal coefficients beta for all independent variables
#'
#' @param  Sigma  : Covariance matrix between equations.
#' @param  Betas  : Vector of beta coefficients including intercept term
#' @param  Thetas : Vector of delta coefficients for Durbin processes. Default: NULL
#' @param  durbin : Logical for durbin type processes. Default: FALSE.
#' @param  lambda : Level of spatial autocorrelation in lag term. A nGx1 vector rho[g]=level of spatial autocorrelation in equation 'g'. If rho=number, then the same level for all equations.
#' @param  rho    : Level of spatial autocorrelation in error term. A nGx1 vector lambda[g]=level of spatial autocorrelation in equation 'g'. If rho=number, then the same level for all equations.
#' @param  p      : A nGx1 vector \eqn{p_{g}}=number of covariates (including intercept) in equation 'g'. If \eqn{p_{g}}=number the same number of variables in all equations.
#' @param  nR     : Number of spatial units.
#' @param  nT     : Number of time periods.
#' @param  nG     : Number of equations in SUR model.
#' @param  W      : A nRxnR spatial weight matrix.
#'
#' @return
#' A list with a matrix Y (nRxnGxnT)x1 with dependent variable and matrix X (nRxnGxnT)xsum(p)
#'
#' @details
#' This function allows to simulate data samples of SUR-SAR; SUR-SEM and SUR-SARAR models.
#' It can be used to make Monte Carlo experiments.
#'
#' The model specification for nR individuals, nT cross sections and nG equations, with spatial interaction mechanisms,
#' as follows:\cr
#' \cr
#' SUR-SAR\cr
#' \deqn{y_{gt}=\lambda_{g} W y_{gt} + X_{gt} \beta_{g} + \epsilon_{gt}} \cr
#' SUR-SEM\cr
#' \deqn{y_{gt}= X_{gt} \beta_{g} + u_{gt}}
#' \deqn{u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt}}\cr
#' SUR-SARAR\cr
#' \deqn{y_{gt}=\lambda_{g} W y_{gt} + X_{gt} \beta_{g} + u_{gt}}
#' \deqn{u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt}}
#' where \eqn{y_{gt}}, \eqn{u_{gt}} and \eqn{\epsilon_{gt}} are (nRx1) vectors; \eqn{X_{gt}} is a matrix
#' of exogenous variables of order (nRxsum(p)); \eqn{\lambda_{g}} and \eqn{\rho_{g}} are
#'  parametres of spatial dependence; W is the nRxnR matrix.

#' @references
#' López, F.A., Mur, J., and Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#' \cr
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@upct.es} \cr
#'   Román Mínguez  \tab \email{Roman.Minguez@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@unizar.es}\cr
#'   }
#' @seealso
#' \code{\link{lmtestspsur}},\code{\link{spsurml}}, \code{\link{spsur3sls}}
#' @examples
#'
#' ####################################
#' ######## CROSS SECTION DATA ########
#' ####################################
#'
#' #### Example 1: DGP SDM model
#' nT <- 1 # Number of time periods
#' nG <- 3 # Number of equations
#' nR <- 500 # Number of spatial elements
#' p <- 3 # Number of independent variables
#' Sigma <- matrix(0.3,ncol=nG,nrow=nG)
#' diag(Sigma) <- 1
#' Betas <- c(1,2,3,1,-1,0.5,1,-0.5,2)
#' Thetas <- c(1,-1,0.5,-0.5,1,0)
#' lambda <- 0.5 # level of spatial dependence
#' rho <- 0.0 # spatial autocorrelation error term = 0
#' #  random coordinates
#' co <- cbind(runif(nR,0,1),runif(nR,0,1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co,k=5,longlat=F)))
#' DGP <- dgp_spSUR(Sigma=Sigma,Betas=Betas,Thetas=Thetas,durbin=TRUE,
#'                  rho=rho,lambda=lambda,nT=nT,nG=nG,nR=nR,p=p,W=W)
#' SDM <- spsurml(W=W,X=DGP$X,Y=DGP$Y,nT=nT,nR=nR,nG=nG,
#'                p=c(3,3,3),type="sdm")
#' SAR <- spsurml(W=W,X=DGP$X,Y=DGP$Y,nT=nT,nR=nR,nG=nG,
#'                p=c(3,3,3),type="sar")

#' ### Example 2: DGP SEM model
#' nT <- 1 # Number of time periods
#' nG <- 3 # Number of equations
#' nR <- 500 # Number of spatial elements
#' p <- c(2,3,4) # Number of independent variables
#' Sigma <- matrix(0.8,ncol=nG,nrow=nG)
#' diag(Sigma) <- 1
#' Betas <- c(1,2,1,2,3,1,2,3,4)
#' lambda <- 0 # level of spatial dependence = 0
#' rho <- c(0.2,0.5,0.8) # spatial autocorrelation error terms for each equation
#' # random coordinates
#' co <- cbind(runif(nR, 0, 1),runif(nR, 0, 1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co,k=5,longlat=F)))
#' DGP <- dgp_spSUR(Sigma=Sigma,Betas=Betas,rho=rho,lambda=lambda,nT=nT,nG=nG,nR=nR,p=p,W=W)
#' SAR <- spsurml(W=W,X=DGP$X,Y=DGP$Y,nT=nT,nR=nR,nG=nG,
#'                p=c(2,3,4),type="sar")
#' SEM <- spsurml(W=W,X=DGP$X,Y=DGP$Y,nT=nT,nR=nR,nG=nG,
#'                p=c(2,3,4),type="sem")
#'
#' ####################################
#' ######## PANEL DATA (nG>1;nT>1) ####
#' ####################################
#'
#' ######## Example 3: Equal number of regressors
#' rho12 <- 0.8; rho13 <- 0; rho23 <- 0.5
#' sigma11 <- 1; sigma22 <- 2; sigma33 <- 1.5
#' sigma12 <- rho12*sqrt(sigma11)*sqrt(sigma22)
#' sigma13 <- rho13*sqrt(sigma11)*sqrt(sigma33)
#' sigma23 <- rho23*sqrt(sigma22)*sqrt(sigma33)
#' Sigma <- matrix(c(sigma11,sigma12,sigma13,sigma12,sigma22,sigma23,sigma13,sigma23,sigma33),nrow=3,ncol=3,byrow=TRUE)
#' # Check Sigma is positive definite
#' eigen(Sigma)
#' Betas <- c(1:3,-1:1,3:5)
#'
#' # Para nR=3000, nG=3 y nT=5 da problemas al estimar matriz covarianzas en SAR.
#' SAR <- spsurml(W=W,Betas=Betas,X=DGP$X,Y=DGP$Y,nT=nT,nR=nR,nG=nG,
#'                p=c(3,3,3),type="sar")
#' summary(SAR)

#' ######## Example 4: With different number of regressors an levels of spatial dependence
#' rho <- c(0.3,0.8) # spatial dependence
#' lambda <- c(-0.5,0.1) # spatial autocorrelation error term
#' nT <- 5 # Number of time periods
#' nG <- 2 # Number of equations
#' nR <- 2000 # Number of spatial elements
#' p <- c(4,2) # Number of independent variables
#' Betas <- c(1,-1,2,0.5,2,3)
#' # ramdom coordinates
#' co <- cbind(runif(nR, 0, 1),runif(nR, 0, 1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co,k=5,longlat=F)))
#'
#' DGP <- dgp_spSUR(Sigma=Sigma,Betas=Betas,rho=rho,lambda=lambda,nT=nT,nG=nG,nR=nR,p=p,W=W)
#'
#' @export
dgp_spSUR <- function(Sigma, nT = 1, nG, nR, Betas,
                      Thetas = NULL, durbin = FALSE, rho = NULL,
                      lambda = NULL, p = NULL, W = NULL, X = NULL)
{
  if(!is.null(W)) W <- Matrix::Matrix(W) else W <- Matrix::Diagonal(nR)
  if (!is.null(p) & length(p)==1) p <- matrix(p,nrow = nG,ncol = 1)
  if (is.null(lambda)) lambda <- rep(0,nG)
  if (is.null(rho)) rho <- rep(0,nG)
  if (length(lambda)==1) lambda <- as.numeric(matrix(lambda,nrow = nG,ncol = 1))
  if (length(rho)==1) rho <- as.numeric(matrix(rho,nrow = nG,ncol = 1))
  if(is.null(X)) {
    if (is.null(p)) stop("Arguments X and p can not be NULL simultaneously")
    X <- cbind(matrix(1,nR,1),
               matrix(runif(nR*(p[1]-1)),nR,(p[1]-1)))
    for (i in 1:(nG-1)){
      X <- Matrix::bdiag(X,cbind(matrix(1,nR,1),
                                 matrix(runif(nR*(p[i+1]-1)),nR,(p[i+1]-1))))
    }
    # nam <- matrix(0,ncol = p[1]+1)
    # for(i in 1:(p[1]+1)){
    #   assign(paste("C", i, sep = ""), i)
    # }
    if (nT>1){
      for (i in 1:(nT-1)){
        X2 <- cbind(matrix(1,nR,1),
                    matrix(runif(nR*(p[1]-1)),
                           nR,(p[1]-1)))
        for (i in 1:(nG-1)){
          X2 <- Matrix::bdiag(X2,cbind(matrix(1,nR,1),
                                       matrix(runif(nR*(p[i+1]-1)),
                                              nR,(p[i+1]-1))))
        }
        X <- rbind(X,X2)
      }
    }

    # Nombro las columnas de X
    nam <- c(paste0("Intercep_",1),
             paste(paste0("X",1,"_"), 1:(p[1]-1), sep=""))
    if (length(p>1)) {
      for (i in 2:(length(p))){
        nam <- c(nam,c(paste0("Intercep_",i),
                       paste(paste0("X",i,"_"), 1:(p[i]-1),
                             sep="")))
      }
    }
    dimnames(X)[[2]] <- nam
  }

  if(is.null(p)) {
    if((ncol(X) %% nG) != 0) stop("Argument p need to be set")
    p <- rep(ncol(X) / nG, nG)
  }

  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)

  # CAMBIA MATRIZ X Y COEFICIENTES EN EL CASO DURBIN
  if (durbin) {
    WX <- (IT %x% IG %x% W) %*% X
    dimnames(WX)[[2]] <- paste0("W_",colnames(X))
    Xdurbin <- NULL
    pdurbin <- p - 1 # Sin intercepto
    for (i in 1:length(p))
    {
      if (i == 1) {
        Xdurbin <- cbind(X[,1:p[i]],WX[,2:p[i]])
        Coeff <- c(Betas[1:p[1]],Thetas[1:pdurbin[1]])
      } else {
        Xdurbin <- cbind(Xdurbin,
                         X[,(cumsum(p)[i-1]+1):cumsum(p)[i]],
                         WX[,(cumsum(p)[i-1]+2):cumsum(p)[i]])
        # Sin intercepto
        Coeff <- c(Coeff,
                   Betas[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                   Thetas[(cumsum(pdurbin)[i-1]+1):cumsum(pdurbin)[i]])
      }
    }
    #p <- p + (p-1)  # Para el caso sdm cambia el p (ojo Intercepto)
  }

  S <- Sigma
  #S <- (Matrix::Matrix(1,nrow = nG,ncol=nG)-IG)*cor + IG
  OME <- (IT %x% S) %x% IR
  # Factor Cholesky covarianzas
  chol_OME <- Matrix::Cholesky(OME)
  #factors_chol_OME <- Matrix::expand(chol_OME)
  #Lchol_OME <- Matrix::t(factors_chol_OME$P) %*% factors_chol_OME$L
  # Uchol_OME <- Matrix::t(factors_chol_OME$L) %*% factors_chol_OME$P
  # Comprobación factor Choleski
  #pr <- as(Matrix::tcrossprod(Lchol_OME),"dgTMatrix")
  #all.equal(OME,pr); range(OME-pr)

  M <- Matrix::Matrix(0,ncol=1,nrow=nT*nG*nR)
  U <- matrix(sparseMVN::rmvn.sparse(n=1, mu=M, CH=chol_OME,prec=FALSE),ncol=1)
  # Alternativa también rápida pero no utiliza sparse y tiene problemas para gran tamaño
  #U <- try(matrix(mvnfast::rmvn(1,mu=M,sigma=OME,isChol=FALSE),ncol=1))
  U <- Matrix::Matrix(U)

  # Si nT*nG*nR es muy grande (>30000 ó 40000) hay problemas
  # Solución: Evitar cálculo inversa matrices
  #IA <- IT %x% Matrix::solve(IGR - diag(rho) %x% W)
  #IB <- IT %x% Matrix::solve(IGR - diag(lambda) %x% W)

  # beta <- as.matrix(c(Betas[1], Betas[2]*matrix(1,ncol=p[1],nrow=1)))
  # for (i in 1:(nG-1)){
  #   beta <- rbind(beta,as.matrix(c(Betas[1], Betas[2]*matrix(1,ncol=p[i+1],nrow=1))))
  # }
  #Y <- IA %*% (X %*% beta + IB %*% U)
  # Con este código evitamos inversas
  IBU <- Matrix::solve(IT %x% (IGR - diag(rho) %x% W),U)
  if (durbin) {
    Y <- Matrix::solve(IT %x% (IGR - diag(lambda) %x% W) ,
                       (Xdurbin %*% Coeff + IBU))
  } else {
    Y <- Matrix::solve(IT %x% (IGR - diag(lambda) %x% W) ,
                       (X %*% Betas + IBU))
  }

  results <- list(X = as.matrix(X), Y = as.matrix(Y))
}

