
fit_spsurslm_3sls <- function(Tm, G, N, Y, X, W, p,
                              type, maxlagW){
  ## Crear matriz instrumentos para retardo espacial Wy
  IT <- Matrix::Diagonal(Tm)
  IG <- Matrix::Diagonal(G)
  IR <- Matrix::Diagonal(N)
  if (inherits(X, "matrix")) X <- Matrix::Matrix(X)
  if (inherits(Y, "matrix")) Y <- Matrix::Matrix(Y)
  if (inherits(W, "matrix")) W <- Matrix::Matrix(W)

  ########## FULL REGRESSOR MATRIX #######################

  #### Fill with 0 the columns of spatial lags Wy_i
  Z <- NULL
  for (i in 1:G)
  {
    if (i==1) {
      Z <- cbind(0,X[,c(1:p[i])])
      colnames(Z)[i] <- "Wy_1"
    } else {
      ncolz <- ncol(Z)
      Z <- cbind(Z,0,X[,c((sum(p[1:(i-1)])+1):sum(p[1:i]))])
      colnames(Z)[ncolz+1] <- paste("Wy_",i,sep="")
    }
  }

  #### Add the values of spatial lags Wy_i
  WY<- (IT %x% IG %x% W) %*%  Y
  m <- nrow(WY) / G
  if ( (nrow(WY) %% G) != 0) stop("full sample size not divisible by number of equations ")
  for (i in 1:G) {
    col_Wy <- paste("Wy_",i,sep="")
    Z[((i-1)*m+1):(i*m),c(col_Wy)] <- WY[((i-1)*m+1):(i*m)]
   }


  ######################### FULL INSTRUMENTS MATRIX ################

  # Build instrument matrix for spatial lag regressor
  # slm model: WX as instruments of Wy
  # sdm model: W^2X as instruments of Wy
  # Drop intercept for the spatial lag of X
  X_noIn <- X[,!grepl("Interc",colnames(X))]
  H <- X
  # Build instruments for slm and sdm
  Wi_X_noIn <- X_noIn
  for (i in 1:maxlagW){
    Wi_X_noIn <- (IT %x% IG %x% W) %*% Wi_X_noIn
    colnames(Wi_X_noIn) <- paste("W",i,"_", colnames(X_noIn), sep = "")
    H <- cbind(Wi_X_noIn,H)
  }

  ############# 2SLS ESTIMATION #####################

  Z_hat <- fitted(lm(as.matrix(Z) ~ as.matrix(H) - 1))
  tsls <- lm(as.matrix(Y) ~ Z_hat-1)
  # summary(tsls)
  resids.tsls <- residuals(tsls)
  Sigmas <- get_Sigma(resids=resids.tsls,N=N,G=G,Tm=Tm)
  Sigma <- Matrix::Matrix(Sigmas$Sigma)
  Sigmainv <- Matrix::Matrix(Sigmas$Sigma_inv)

  #################### 3SLS ESTIMATION
  Z_hat <- Matrix::Matrix(Z_hat)
  OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
  full_3SLS <- Matrix::solve(Matrix::crossprod(Z_hat,OMEinv %*% Z_hat),
                        Matrix::crossprod(Z_hat,OMEinv %*% Y))
  rownames(full_3SLS) <- sub("Wy","rho",rownames(full_3SLS))
  full_3SLS <- as.matrix(full_3SLS)
  cov_3SLS <- Matrix::solve( Matrix::t(Z_hat) %*%
                                     OMEinv %*% Z_hat )
  colnames(cov_3SLS) <- rownames(full_3SLS)
  rownames(cov_3SLS) <- rownames(full_3SLS)
  cov_3SLS <- as.matrix(cov_3SLS)

  se_full_3SLS <- as.vector(sqrt(diag(cov_3SLS)))
  names(se_full_3SLS) <- rownames(full_3SLS)

  Yhat <- Z %*% full_3SLS ## OJO: USAMOS REGRESORES ORIGINALES
  Res <- Y - Yhat

  betas_3SLS <- full_3SLS[!grepl("rho",rownames(full_3SLS)),]
  se_betas_3SLS <- se_full_3SLS[!grepl("rho",names(se_full_3SLS))]

  deltas_3SLS <- full_3SLS[grepl("rho",rownames(full_3SLS)),]
  se_deltas_3SLS <- se_full_3SLS[grepl("rho",names(se_full_3SLS))]
   res <- list(deltas = deltas_3SLS,
              deltas.se = se_deltas_3SLS,
              coefficients = betas_3SLS,
              rest.se = se_betas_3SLS,
              resvar = as.matrix(cov_3SLS),
              Sigma = as.matrix(Sigma),
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )

 }
