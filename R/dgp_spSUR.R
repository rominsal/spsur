#' @name dgp_spsur
#' @rdname dgp_spsur
#'
#' @title Generation of a random dataset with a spatial SUR structure.
#'
#' @description
#'  The purpose of the function \code{dgp_spsur} is to generate a random dataset
#'  with the dimensions and spatial structure decided by the user. This function may be useful in
#'  pure simulation experiments or with the aim of showing specific properties and characteristics
#'  of a spatial SUR dataset and inferential procedures related to them.
#'
#'  The user of \code{dgp_spsur} should think in terms of a Monte Carlo experiment.
#'  The arguments of the funtion specify the dimensions of the dataset to be generated, the spatial
#'  mechanism underlying the data, the intensity of the SUR structure among the equations
#'  and the values of the parameters to be used to obtain the simulated data, which includes the error terms,
#'  the regressors and the explained variables.
#'
#'
#'dgp_spsur <- function(Sigma, Tm = 1, G, N, Betas,
#'Thetas = NULL, durbin = FALSE, rho = NULL,
#'lambda = NULL, p = NULL, W = NULL, X = NULL,
#'pdfU = "nvrnorm", pdfX = "nvrnorm")
#'
#'
#' @param Tm Number of time periods.Default = \code{1}
#' @param p Number of regressors by equation, including the intercept. \emph{p} can be a row
#' vector of order \emph{(1xG)}, if the number of regressors is not the same for all the
#' equations, or a scalar, if the \emph{G} equations have the same number of regressors.
#' @param Sigma Covariance matrix between the \emph{G} equations of the SUR model.
#'   This matrix should be definite positive and the user must check for that.
#' @param Betas A row vector of order \eqn{(1xP)} showing the values for the \emph{beta} coefficients.
#'  The first \eqn{P_{1}} terms correspond to the first equation (where the first element is the intercept),
#'  the second \eqn{P_{2}} terms to the coefficients of the second equation and so on.
#' @param  Thetas Values for the \eqn{\theta} coefficients in the \emph{G} equations of the model,
#'  when the type of spatial SUR model to be simulated is a \strong{"slx"}, \strong{"sdm"}  or \strong{"sdem"}.
#'   \emph{Thetas} is a row vector of order \emph{\eqn{1xPTheta}}, where \emph{\eqn{PThetas=p-G}}; let us note
#'   that the intercept cannot appear among the spatial lags of the regressors. The first \emph{\eqn{1xKTheta_{1}}}
#'   terms correspond to the first equation, the second \emph{\eqn{1xPTheta_{2}}} terms correspond to the
#'   second equation, and so on. Default = \code{NULL}.
#' @param lambda Values of the coefficients \eqn{\lambda_{g}; g=1,2,..., G} related to the spatial lag of
#'   the explained variable of the g-th equation. If \eqn{lambda} is an scalar and there are \emph{G} equations
#'   in the model, the same value will be used for all the equations. If \eqn{lambda} is a row vector,
#'   of order \emph{(1xG)}, the function \code{dgp_spsur} will use these values,
#'   one for each equation. Default = \code{NULL}.
#' @param rho Values of the coefficients \eqn{\rho_{g}; g=1,2,..., G} related to the spatial lag of
#'   the errors in the \emph{G} equations. If \eqn{rho} is an scalar and there are \emph{G} equations
#'   in the model, the same value will be used for all the equations. If \eqn{rho} is a row vector,
#'   of order \emph{(1xG)}, the function \code{dgp_spsur} will use these values,
#'   one for each equation of the spatial errors. Default = \code{NULL}.
#' @param X This argument tells the function \code{dgp_spsur} which \strong{X} matrix should be used to
#'   generate the SUR dataset. If \strong{X} is different from \code{NULL}, \code{{dgp_spsur}}
#'   will upload the \strong{X} matrix selected in this argument. Note that the \strong{X} must be consistent
#'   with the dimensions of the model. If \strong{X} is \code{NULL}, \code{dgp_spsur} will
#'   generate the desired matrix of regressors from a multivariate Normal distribution with mean value zero and
#'   identity \eqn{(PxP)} covariance matrix. As an alternative, the user may change this probability distribution
#'   function to the uniform case, \eqn{U(0,1)}, through the argument \strong{pdfX}. Default = \code{NULL}.
#' @param pdfX  Multivariate probability distribution function, mdf, from which the values of the regressors
#'   will be drawn. The regressors are assumed to be independent. \code{dgp_spsur} provides
#'   two mdf, the multivariate Normal, which is the default, and the uniform in the interval \eqn{U[0,1]}, using the
#'   dunif function, \code{\link[stats]{dunif}}, from the \pkg{stats} package. Default = \code{"nvrnorm"}.
#' @param pdfU Multivariate probability distribution function, mdf, from which the values of the error terms will
#'  be drawn. The covariance matrix is the \eqn{\Sigma} matrix specificied by the user in the argument \strong{Sigma}.
#'  The funtion \code{dgp_spsur} provides two mdf, the multivariate Normal, which is the default,
#'  and the log-Normal distribution funtion which means just exponenciate the sampling drawn form a \eqn{N(0,\Sigma)}
#'  distribution. Default = \code{"nvrnorm"}.
#' @inheritParams spsurml
#'
#'
#' @details
#'  The purpose of the function \code{dgp_spsur} is to generate random datasets, of a SUR
#'  nature, with the spatial structure decided by the user. The function requires certain information to be
#'  supplied externally because, in fact, \code{dgp_spsur} constitutes a Data Generation
#'  Process, DGP. The following aspects should be addressed:
#'  \itemize{
#'     \item The user must define the dimensions of the dataset, that is, number of equations,
#'     \emph{G}, number of time periods, \emph{Tm}, and number of cross-sectional units, \emph{N}.
#'
#'     \item Then, the user must choose the type of spatial structure desired for the model from among
#'      the list of candidates of \strong{"sim"}, \strong{"slx"}, \strong{"slm"}, \strong{"sem"},
#'      \strong{"sdm"}, \strong{"sdem"} or  \strong{"sarar"}; the default is the \strong{"sim"} specification
#'      which does not have spatial structure. The decision is made implicitly, just omiting the
#'      specification of the spatial parameters which are not involved  in the model (i.e., in a \strong{"slm"}
#'      there are no \eqn{\rho} parameters but appear \eqn{\lambda} parameters; in a \strong{"sdem"} model there
#'      are \eqn{\rho} and \eqn{\theta} parameters but no \eqn{\lambda} coefficients). Of course, if the user
#'      needs a model with spatial structure, a \emph{(nxN)}  weighting matrix, \emph{W}, should be chosen.
#'     \item The next step builds the equations of the SUR model. In this case, the user must specify
#'      the number of regressors that intervene in each equation and the coefficients, \eqn{\beta} parameters,
#'      associated with each regressor. The \emph{first} question is solved through the argument \emph{p}
#'      which, if a scalar, indicates that the same number of regressors should appear in all the equations
#'      of the model; if the user seeks for a model with different number of regressors in the
#'      \emph{G} equations, the argument \emph{p} must be a \emph{(1xG)} row vector with the required
#'      information. It must be remembered that \code{dgp_spsur} assumes that an
#'      intercept appears in all equations of the model.
#'
#'      The \emph{second} part of the problem posited above is solved through the argument \emph{Betas},
#'      which is a row vector of order \emph{(1xp)} with the information requiered for this set of
#'      coeficcients.
#'    \item The user must specify, also, the values of the spatial parameters corresponding to the chosen
#'      specification; we are refering to the \eqn{\lambda_{g}}, \eqn{\rho_{g}} and  \eqn{\theta_{g}},
#'      for \eqn{g=1, ..., G and k=1,..., K_{g}} parameters. This is done throught the arguments \emph{lambda},
#'      \emph{rho} and \emph{theta}. The firs two, \emph{lambda} and \emph{rho}, work as \emph{K}: if they are
#'      scalar, the same value will be used in the \emph{G} equations of the SUR model; if they are
#'      \emph{(1xG)} row vectors, a different value will be assigned for each equation.
#'
#'      Moreover, \emph{theta} works like the argument \emph{beta}. The user must define a row vector of
#'      order \eqn{1xPTheta} showing these values. It is worth to remember that in no case the intercept
#'      will appear among the lagged regressors.
#'    \item Finally, the user must decide which values of the regressors and of the error terms are to be used
#'     in the simulation. The regressors can be uploaded from an external matrix generated previously by the
#'     user. This is the argument \emph{X}. It is the responsability of the user to check that the dimensions
#'     of the external matrix are consistent with the dataset required for the SUR model. A second possibility
#'     implies  the regressors to be generated randomly by the function \code{\link{dgp_spsur}}.
#'     In this case, the user must select the probability distribution function from which the corresponding
#'     data (of the regressors and the error terms) are to be drawn.\cr
#'
#'      \code{dgp_spsur} provides two multivariate distibution functions, namely, the Normal and the log-Normal for the
#'     errors (the second should be taken as a clear departure from the standard assumption of
#'     normality). In both cases, random matrices of order \emph{(TmNxG)} are obtained from a multivariate
#'     normal distribution, with a mean value of zero and the covariance matrix specified in the argument
#'     \emph{Sigma}; then, this matrix is exponentiated for the log-Normal case. Roughly, the same procedure
#'     applies for drawing the values of the regressor. There are two distribution functions available, the
#'     normal and the uniform in the interval \eqn{U[0,1]}; the regressors are always independent.
#'
#'   }
#'
#' @return
#' A list with a vector \eqn{Y} of order \emph{(TmNGx1)} with the values generated for the explained
#' variable in the G equations of the SUR and a matrix \eqn{XX} of order (\emph{(TmNGxsum(p))}, with the values
#' generated for the regressors of the SUR, including an intercept for each equation.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{spsur3sls}}, \code{\link{spsurtime}}
#' @examples
#'
#' ####################################
#' ######## CROSS SECTION DATA ########
#' ####################################
#'
#' ####################################
#' #### Example 1: DGP SDM model
#' ####################################
#' rm(list = ls()) # Clean memory
#' Tm <- 1 # Number of time periods
#' G <- 3 # Number of equations
#' N <- 500 # Number of spatial elements
#' p <- 3 # Number of independent variables
#' Sigma <- matrix(0.3, ncol = G, nrow = G)
#' diag(Sigma) <- 1
#' Betas <- c(1,2,3,1,-1,0.5,1,-0.5,2)
#' Thetas <- c(1,-1,0.5,-0.5,1,0)
#' lambda <- 0.5 # level of spatial dependence
#' rho <- 0.0 # spatial autocorrelation error term = 0
#' #  random coordinates
#' co <- cbind(runif(N,0,1),runif(N,0,1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co, k = 5,
#'                                                    longlat = FALSE)))
#'
#' DGP <- dgp_spsur(Sigma = Sigma, Betas = Betas, Thetas = Thetas,
#'                  rho = rho, lambda = lambda, Tm = Tm,
#'                  G = G, N = N, p = p, W = W)
#' SDM <- spsurml(W = W, X = DGP$X, Y = DGP$Y, Tm = Tm, N = N, G = G,
#'                p = c(3,3,3), type = "sdm")
#' summary(SDM)
#' SLM <- spsurml(W = W, X = DGP$X, Y = DGP$Y, Tm = Tm, N = N, G = G,
#'                p = c(3,3,3), type = "slm")
#' summary(SLM)
#'
#' ####################################
#' ### Example 2: DGP SEM model with Tm>1; G=1 and
#' ### different p for each equation
#' ####################################
#' rm(list = ls()) # Clean memory
#' Tm <- 3 # Number of time periods
#' G <- 1 # Number of equations
#' N <- 500 # Number of spatial elements
#' p <- c(2,3,4) # Number of independent variables
#' Sigma <- matrix(0.8, ncol = Tm, nrow = Tm)
#' diag(Sigma) <- 1
#' Betas <- c(1,2,1,2,3,1,2,3,4)
#' lambda <- 0 # level of spatial dependence = 0
#' rho <- c(0.2,0.5,0.8) # spatial autocorrelation error terms for each equation
#' # random coordinates
#' co <- cbind(runif(N, 0, 1),runif(N, 0, 1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co, k = 5,
#'                                                    longlat = FALSE)))
#' DGP2 <- dgp_spsur(Sigma = Sigma, Betas = Betas, rho = rho, lambda = lambda,
#'                  Tm = Tm, G = G, N = N, p = p, W = W)
#' SLM2 <- spsurml(W = W, X = DGP2$X, Y = DGP2$Y, Tm = Tm, N = N, G = G,
#'                p = c(2,3,4), type = "slm")
#' summary(SLM2)
#' SEM2 <- spsurml(W = W, X = DGP2$X, Y = DGP2$Y, Tm = Tm, N = N, G = G,
#'                p = c(2,3,4), type = "sem")
#' summary(SEM2)
#'
#' ####################################
#' #### Example with G>1 and Tm>>1
#' ####################################
#' rm(list = ls()) # Clean memory
#' Tm <- 10 # Number of time periods
#' G <- 3 # Number of equations
#' N <- 100 # Number of spatial elements
#' p <- 3 # Number of independent variables
#' Sigma <- matrix(0.5, ncol = G, nrow = G)
#' diag(Sigma) <- 1
#' Betas <- rep(1:3, G)
#' lambda <- c(0.5, 0.1, 0.8)
#' rho <- 0.0 # spatial autocorrelation error term = 0
#' #  random coordinates
#' co <- cbind(runif(N,0,1),runif(N,0,1))
#' W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co, k = 5,
#'                                                    longlat = FALSE)))
#' DGP3 <- dgp_spsur(Sigma = Sigma, Betas = Betas, rho = rho, lambda = lambda,
#'                  Tm = Tm, G = G, N = N, p = p, W = W)
#' SLM3  <-spsurml(Y= DGP3$Y, X = DGP3$X, G = G, N = N, Tm = Tm,
#'                   p = p, W = W, type = "slm")
#' summary(SLM3)
#'
#' ## slm with demeaning
#'
#' SLM3_dem  <-spsurml(Y= DGP3$Y, X = DGP3$X, G = G, N = N, Tm = Tm,
#'                   p = p, W = W, type = "slm", demean = TRUE)
#' summary(SLM3_dem)
#'
#' @export
dgp_spsur <- function(Sigma, Tm = 1, G, N, Betas,
                      Thetas = NULL, rho = NULL,
                      lambda = NULL, p = NULL, W = NULL, X = NULL,
                      pdfU = "nvrnorm", pdfX = "nvrnorm")
{
  #check for row-standardization of W
  if (!is.null(W)){
    if (class(W) != "matrix") W <- as.matrix(W)
    rsumW <- rowSums(W)
    rsumW[rsumW == 0] <- 1
    nW <- dim(W)[1]
    W <- W / matrix(rep(rsumW, each = nW),
                    nrow = nW, ncol = nW, byrow = TRUE)
    W <- Matrix::Matrix(W)
  } else W <- Matrix::Diagonal(N)

  if (Tm > 1 && G == 1){ #Change dimensions
    G <- Tm
    Tm <- 1
  }
  if (!is.null(Thetas)) durbin <- TRUE else durbin <- FALSE
  if (!is.null(p) & length(p)==1) p <- matrix(p,nrow = G,ncol = 1)
  if (is.null(lambda)) lambda <- rep(0,G)
  if (is.null(rho)) rho <- rep(0,G)
  if (length(lambda) == 1) lambda <- as.numeric(matrix(lambda,nrow = G,ncol = 1))
  if (length(rho) == 1) rho <- as.numeric(matrix(rho,nrow = G,ncol = 1))
  if (is.null(X)){
    if (is.null(p)) stop("Arguments X and p can not be NULL simultaneously")
    if (pdfX == "nvrunif"){
      X <- cbind(matrix(1,N,1),
                 matrix(runif(N*(p[1]-1)),N,(p[1]-1)))
      for (i in 1:(G-1)){
        X <- Matrix::bdiag(X,cbind(matrix(1,N,1),
                                   matrix(runif(N*(p[i+1]-1)),N,(p[i+1]-1))))
      }
      if (Tm > 1){
        for (i in 1:(Tm-1)){
          X2 <- cbind(matrix(1,N,1),
                      matrix(runif(N*(p[1]-1)),
                             N,(p[1]-1)))
          for (i in 1:(G-1)){
            X2 <- Matrix::bdiag(X2,cbind(matrix(1,N,1),
                                         matrix(runif(N*(p[i+1]-1)),
                                                N,(p[i+1]-1))))
          }
          X <- rbind(X,X2)
        }
      }
    } else if (pdfX == "nvrnorm"){
      X <- cbind(matrix(1,N,1),
                 matrix(rnorm(N*(p[1]-1),0,1),N,(p[1]-1)))
      for (i in 1:(G-1)){
        X <- Matrix::bdiag(X,cbind(matrix(1,N,1),
                               matrix(rnorm(N*(p[i+1]-1),0,1),N,(p[i+1]-1))))
      }
      if (Tm > 1){
        for (i in 1:(Tm-1)){
          X2 <- cbind(matrix(1,N,1),
                      matrix(rnorm(N*(p[1]-1),0,1),N,(p[1]-1)))
          for (i in 1:(G-1)){
            X2 <- Matrix::bdiag(X2,cbind(matrix(1,N,1),
                                         matrix(rnorm(N*(p[i+1]-1),0,1),
                                                N,(p[i+1]-1))))
          }
          X <- rbind(X,X2)
        }
      }
    } else stop("pdfX only can be nvrnorm or nvrunif")

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
    if((ncol(X) %% G) != 0) stop("Argument p need to be set")
    p <- rep(ncol(X) / G, G)
  }

  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)

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
  #S <- (Matrix::Matrix(1,nrow = G,ncol=G)-IG)*cor + IG
  OME <- (IT %x% S) %x% IR
  # Factor Cholesky covarianzas
  chol_OME <- Matrix::Cholesky(OME)
  #factors_chol_OME <- Matrix::expand(chol_OME)
  #Lchol_OME <- Matrix::t(factors_chol_OME$P) %*% factors_chol_OME$L
  # Uchol_OME <- Matrix::t(factors_chol_OME$L) %*% factors_chol_OME$P
  # Comprobación factor Choleski
  #pr <- as(Matrix::tcrossprod(Lchol_OME),"dgTMatrix")
  #all.equal(OME,pr); range(OME-pr)

  M <- Matrix::Matrix(0,ncol=1,nrow=Tm*G*N)
  U <- matrix(sparseMVN::rmvn.sparse(n=1, mu=M, CH=chol_OME,prec=FALSE),ncol=1)
  U <- Matrix::Matrix(U)
  if (pdfU == "lognvrnorm") U <- exp(U)
  if (pdfU != "lognvrnorm" && pdfU != "nvrnorm")
    print(" Improper pdf. The errors will be drawn from a multivariate Normal ")

  # Si Tm*G*N es muy grande (>30000 ó 40000) hay problemas
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

