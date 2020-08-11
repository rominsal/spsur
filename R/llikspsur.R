f_sur_sim <- function(env) {
  # Log-lik SUR-SIM Spatio-Temporal Model (No Spatial Effects)
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X; Sigma <- env$Sigma
  Sigmainv <- Matrix::solve(Sigma)
  # Log-lik SUR-SLM Spatio-Temporal Model
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  OMEinv <- Matrix::kronecker(Matrix::kronecker(IT, Sigmainv),IR)
  B_sim <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                         Matrix::crossprod(X, OMEinv %*% Y))
  Res <- Y - X %*% B_sim
  ldet_Sigma <- Matrix::determinant(Sigma, logarithm = TRUE)$modulus
  llike <- as.numeric( -(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma
                       - (1/2)*Matrix::crossprod(Res, OMEinv %*% Res) )
  return((-1)*llike)
}

######################################################

f_sur_lag <- function(deltag, env){
  # Log-lik SUR-SLM Spatio-Temporal Model
  if (!is.null(env$W)) {
    W <- env$W
  }  else {
    W <- as(env$listw, "CsparseMatrix")
  }  
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X; Sigma <- env$Sigma
  Sigmainv <- Matrix::solve(Sigma)
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(N*G)
  lDAg <- Matrix::Matrix(0, nrow = length(deltag), ncol = 1)
  for (i in 1:G)
  {
    lDAg[i] <- spatialreg::do_ldet(deltag[i], env)
  }
  delta <- Matrix::Diagonal(length(deltag),deltag)
  AY <- Matrix::kronecker(IT, 
          (IGR - Matrix::kronecker(delta, W))) %*% matrix(Y, ncol = 1)
  OMEinv <- Matrix::kronecker(Matrix::kronecker(IT, Sigmainv), IR)
  B_slm <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                         Matrix::crossprod(X, OMEinv %*% AY))
  Res <- AY - X %*% B_slm
  ldet_Sigma <- Matrix::determinant(Sigma,logarithm = TRUE)$modulus
  llike <- as.numeric( -(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDAg)
                     - (1/2)*Matrix::crossprod(Res, OMEinv %*% Res) )
  # Minimize function
  return((-1)*llike)
}

######################################################

f_sur_sem <- function(deltag, env){
  # Log-lik SUR-SEM Spatio-Temporal Model
  if (!is.null(env$W)) {
    W <- env$W
  }  else {
    W <- as(env$listw, "CsparseMatrix")
  }  
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X; Sigma <- env$Sigma
  Sigmainv <- Matrix::solve(Sigma)
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(N*G)
  lDBg <- Matrix::Matrix(0, nrow = length(deltag), ncol = 1)
  for (i in 1:G)
  {
    lDBg[i] <- spatialreg::do_ldet(deltag[i], env)
  }
  delta <- Matrix::Diagonal(length(deltag), deltag)
  BY <- Matrix::kronecker(IT, 
              (IGR - Matrix::kronecker(delta, W))) %*% Y
  BX <- Matrix::kronecker(IT, 
              (IGR - Matrix::kronecker(delta, W))) %*% X
  OMEinv <- Matrix::kronecker(Matrix::kronecker(IT, Sigmainv), IR)
  B_sem <- Matrix::solve(Matrix::crossprod(BX, OMEinv %*% BX),
                        Matrix::crossprod(BX, OMEinv %*% BY))
  Res <- Y - (X %*% B_sem)
  BRes <- Matrix::kronecker(IT, 
              (IGR - Matrix::kronecker(delta, W))) %*% Res
  ldet_Sigma <- Matrix::determinant(Sigma,logarithm = TRUE)$modulus
  llike <- as.numeric( -(Tm*G*N/2)*log(2*pi)
                       - (N*Tm/2)*ldet_Sigma
                       + Tm*sum(lDBg)
                       - (1/2)*Matrix::crossprod(BRes, 
                                                 OMEinv %*% BRes) )
  return((-1)*llike)
}

######################################################

f_sur_sarar <- function(DELTA, env){
  # Log-lik SUR-SARAR Spatio-Temporal Model
  if (!is.null(env$W)) {
    W <- env$W
  }  else {
    W <- as(env$listw, "CsparseMatrix")
  }  
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X; Sigma <- env$Sigma
  Sigmainv <- Matrix::solve(Sigma)
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  deltag <- DELTA[1:G]
  deltah <- DELTA[(G+1):(2*G)]
  lDAg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    #lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDAg[i] <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  lDBg <- Matrix::Matrix(0,nrow=length(deltah),ncol=1)
  for (i in 1:G)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    lDBg[i] <- Matrix::determinant(IR - deltah[i]*W, 
                                   logarithm = TRUE)$modulus
  }
  deltaG <- Matrix::Diagonal(length(deltag), deltag)
  deltaH <- Matrix::Diagonal(length(deltah), deltah)
  AY <- Matrix::kronecker(IT, 
            IGR - Matrix::kronecker(deltaG,W)) %*% matrix(Y, ncol = 1)
  OME <- Matrix::kronecker(IT, Matrix::kronecker(Sigma, IR))
  BX <- Matrix::kronecker(IT, IGR - Matrix::kronecker(deltaH,W)) %*% X
  BAY <- Matrix::kronecker(IT, 
                           IGR - Matrix::kronecker(deltaH,W)) %*% AY
  B_sarar <- Matrix::solve(Matrix::crossprod(BX, Matrix::solve(OME, BX)),
                           Matrix::crossprod(BX, Matrix::solve(OME, BAY)))
  rm(BX)
  Res <- AY - (X %*% B_sarar)
  BRes <- Matrix::kronecker(IT, 
                            IGR - Matrix::kronecker(deltaH,W)) %*% Res
  ldet_Sigma <- Matrix::determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi)
                      - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDAg) + Tm*sum(lDBg)
                      - (1/2)*Matrix::crossprod(BRes,
                                            Matrix::solve(OME,BRes)))

  return((-1)*llike)
}

######################################################


