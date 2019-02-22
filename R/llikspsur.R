f_sur_sim <- function(Tm,G,N,Y,X,Sigma)
{
  # Log-lik SUR-SIM Spatio-Temporal Model
  # (No Spatial Effects)
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  # OME <- kronecker(IT, kronecker(Sigma,IR))
  Sigma_Inv <- solve(Sigma)
  YY <- matrix(Y,nrow = N)
  ITKSI <- kronecker(IT,Sigma_Inv)
  IOME <- kronecker(ITKSI,IR)
  YYOME <- matrix(YY%*%ITKSI,ncol = 1)
  B_sim <- Matrix::solve(Matrix::crossprod(X,IOME%*%X),
                         Matrix::crossprod(X,YYOME))
  # B_sim <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
  #                        Matrix::crossprod(X,Matrix::solve(OME,Y)))
  Res <- matrix(Y - X%*%B_sim,ncol=1)
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  RES <- matrix(Res,nrow = N)
  IOMER<-matrix(RES%*%kronecker(IT,Sigma_Inv),ncol = 1)
  # chol_Sigma <- chol(Sigma)
  # ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  # llike <- as.numeric(-(Tm*G*N/2)*log(2*pi) -
  #                     (N*Tm/2)*ldet_Sigma  -
  #                    (1/2)*Matrix::crossprod(Res,Matrix::solve(OME,Res)))
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi) -
                        (N*Tm/2)*ldet_Sigma  -
                        (1/2)*Matrix::crossprod(Res,IOMER))
  return((-1)*llike)

}
######################################################

f_sur_lag <- function(deltag,Tm,G,N,Y,X,W,Sigma){
  W <- as(W,"dgCMatrix")
  # Log-lik SUR-SLM Spatio-Temporal Model
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  lDAg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    #lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDAg[i]  <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  delta <- Matrix::Matrix(diag(as.vector(deltag)))
  # A <- as(kronecker(IT,(IGR - kronecker(delta,W))),"dgCMatrix")
  YY <- matrix(Y,nrow = G*N)
  AY <- matrix((IGR - kronecker(delta,W))%*%YY,ncol=1)
  OME <- as(kronecker(kronecker(IT,Sigma),IR),"dgCMatrix")
  B_slm <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
                         Matrix::crossprod(X,Matrix::solve(OME,AY)))
  Res <- matrix(AY - X%*%B_slm,nrow=nrow(Y))
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  #chol_Sigma <- chol(Sigma)
  #ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDAg) - (1/2)*Matrix::crossprod(Res,Matrix::solve(OME,Res)))
  # Minimize function
  return((-1)*llike)
}

######################################################

f_sur_sem <- function(deltag,Tm,G,N,Y,X,W,Sigma){
  W <- as(W,"dgCMatrix")
  # Log-lik SUR-SEM Spatio-Temporal Model
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  lDBg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDBg[i]  <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  delta <- Matrix::Matrix(diag(as.vector(deltag)))
  B <- kronecker(IT,(IGR - kronecker(delta,W)))
  YY <- matrix(Y,nrow = G*N)
  BY <- matrix((IGR - kronecker(delta,W))%*%YY,ncol=1)
  OME <- as(kronecker(kronecker(IT,Sigma),IR),"dgCMatrix")
  BX <- B%*%X
  B_sem <-Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                        Matrix::crossprod(BX,Matrix::solve(OME,BY)))
  rm(BX)
  Res <- Y - X%*%B_sem
  BRes <- B%*%Res
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi)
                      - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDBg)
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))
  return((-1)*llike)
}

######################################################

f_sur_sarar <- function(DELTA,Tm,G,N,Y,X,W,Sigma){
  W <- as(W,"dgCMatrix")
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
  lDBg <- matrix(0,nrow=length(deltah),ncol=1)
  for (i in 1:G)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    lDBg[i] <- Matrix::determinant(IR-deltah[i]*W,logarithm=TRUE)$modulus
  }
  deltaG <- Matrix::Matrix(diag(as.vector(deltag)))
  deltaH <- Matrix::Matrix(diag(as.vector(deltah)))
  # A <- as(kronecker(IT,(IGR - kronecker(deltaG,W))),"dgCMatrix")
  YY <- matrix(Y,nrow = G*N)
  AY <- matrix((IGR - kronecker(deltaG,W))%*%YY,ncol=1)
  B <- kronecker(IT,IGR - kronecker(deltaH,W))
  OME <- as(kronecker(IT,kronecker(Sigma,IR)),"dgCMatrix")
  BX <- B%*%X
  # Sigma_inv <- solve(Sigma)
  B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                           Matrix::crossprod(BX,Matrix::solve(OME,B%*%AY)))
  rm(BX)
  Res <- AY - X%*%B_sarar
  BRes <- B%*%Res
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi)
                      - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDAg) + Tm*sum(lDBg)
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))

  return((-1)*llike)
}

######################################################

llsur_sarar_full <- function(param,p,Tm,G,N,Y,X,W)
{
  W <- as(W,"dgCMatrix")
  pfull <- sum(p)
  B_sarar <- param[1:pfull]
  DELTA <- param[(pfull+1):((pfull+1)+2*G-1)]
  sigmas <- param[((pfull+1)+2*G):length(param)]
  Sigma <- Matrix::Matrix(0,nrow=G,ncol=G)
  for (i in 1:G) {
    for (j in 1:G) {
      if (i<=j) { #Upper Triangular
        Sigma[i,j] <- sigmas[paste("sigma_",i,j,sep="")]
      }
      if (j<i) Sigma[i,j] <- Sigma[j,i] # Lower Triangular
    }
  }
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
  lDBg <- matrix(0,nrow=length(deltah),ncol=1)
  for (i in 1:G)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    lDBg[i] <- Matrix::determinant(IR-deltah[i]*W,logarithm=TRUE)$modulus
  }
  deltaG <- Matrix::Matrix(diag(deltag))
  deltaH <- Matrix::Matrix(diag(deltah))
  A <- as(kronecker(IT,(IGR - kronecker(deltaG,W))),"dgCMatrix")
  B <- as(kronecker(IT,(IGR - kronecker(deltaH,W))),"dgCMatrix")
  OME <- kronecker(IT,kronecker(Sigma,IR))
  # BX <- B %*% X
  # Sigma_inv <- Matrix::solve(Sigma)
  # B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
  #                          Matrix::crossprod(BX,
  #                                    Matrix::solve(OME,B%*%(A%*%Y))))
  # rm(BX)
  Res <- A %*% Y - X %*% B_sarar
  BRes <- B %*% Res
  ldet_Sigma <- Matrix::determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi)
                      - (N*Tm/2)*ldet_Sigma
                      + Tm*sum(lDAg) + Tm*sum(lDBg)
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))
  return(-1*llike)
}
