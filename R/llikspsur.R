f_sur_sim <- function(nT,nG,nR,Y,X,Sigma)
{
  # Log-lik SUR-SIM Spatio-Temporal Model
  # (No Spatial Effects)
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  OME <- kronecker(IT, kronecker(Sigma,IR))
  B_sim <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
                         Matrix::crossprod(X,Matrix::solve(OME,Y)))
  Res <- matrix(Y - X%*%B_sim,ncol=1)
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  # chol_Sigma <- chol(Sigma)
  # ldet_Sigma <- sum(2*log(diag(chol_Sigma)))

  llike <- as.numeric(-(nT*nG*nR/2)*log(2*pi) -
                      (nR*nT/2)*ldet_Sigma  -
                     (1/2)*Matrix::crossprod(Res,Matrix::solve(OME,Res)))
  return((-1)*llike)                   
}
######################################################


f_sur_lag <- function(deltag,nT,nG,nR,Y,X,W,Sigma)
{
  # Log-lik SUR-SLM Spatio-Temporal Model
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  lDAg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    #lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDAg[i]  <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  delta <- Matrix::Matrix(diag(as.vector(deltag)))
  A <- kronecker(IT,(IGR - kronecker(delta,W)))
  OME <- kronecker(kronecker(IT,Sigma),IR)
  B_sar <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
                          Matrix::crossprod(X,Matrix::solve(OME,A%*%Y)))
  Res <- matrix(A%*%Y - X%*%B_sar,nrow=nrow(Y))
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
   #chol_Sigma <- chol(Sigma)
  #ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  llike <- as.numeric(-(nT*nG*nR/2)*log(2*pi) - (nR*nT/2)*ldet_Sigma
                      + nT*sum(lDAg) - (1/2)*Matrix::crossprod(Res,Matrix::solve(OME,Res)))
  # Minimize function
  return((-1)*llike)
}
######################################################

f_sur_sem <- function(deltag,nT,nG,nR,Y,X,W,Sigma)
{
  # Log-lik SUR-SEM Spatio-Temporal Model
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  lDBg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDBg[i]  <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  delta <- Matrix::Matrix(diag(as.vector(deltag)))
  B <- kronecker(IT,(IGR - kronecker(delta,W)))
  OME <- kronecker(kronecker(IT,Sigma),IR)
  BX <- B%*%X
  B_sem <-Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                                Matrix::crossprod(BX,Matrix::solve(OME,B%*%Y)))  
  rm(BX)
  Res <- Y - X%*%B_sem
  BRes <- B%*%Res
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(nT*nG*nR/2)*log(2*pi) 
                      - (nR*nT/2)*ldet_Sigma
                      + nT*sum(lDBg) 
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))
  return((-1)*llike)                    
}
######################################################

f_sur_sarar <- function(DELTA,nT,nG,nR,Y,X,W,Sigma)
{
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  deltag <- DELTA[1:nG]
  deltah <- DELTA[(nG+1):(2*nG)]
  lDAg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    #lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDAg[i] <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  lDBg <- matrix(0,nrow=length(deltah),ncol=1)
  for (i in 1:nG)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    lDBg[i] <- Matrix::determinant(IR-deltah[i]*W,logarithm=TRUE)$modulus
  }
  deltaG <- Matrix::Matrix(diag(as.vector(deltag)))
  deltaH <- Matrix::Matrix(diag(as.vector(deltah)))
  A <- kronecker(IT,(IGR - kronecker(deltaG,W)))
  B <- kronecker(IT,(IGR - kronecker(deltaH,W)))
  OME <- kronecker(IT,kronecker(Sigma,IR))
  BX <- B%*%X
  Sigma_inv <- solve(Sigma)
  B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                           Matrix::crossprod(BX,Matrix::solve(OME,B%*%(A%*%Y)))) 
  rm(BX)
  Res <- A%*%Y - X%*%B_sarar
  BRes <- B%*%Res
  ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(nT*nG*nR/2)*log(2*pi) 
                      - (nR*nT/2)*ldet_Sigma
                      + nT*sum(lDAg) + nT*sum(lDBg) 
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))
  return((-1)*llike)
}

######################################################

llsur_sarar_full <- function(param,p,nT,nG,nR,Y,X,W)
{
  pfull <- sum(p)
  B_sarar <- param[1:pfull]
  DELTA <- param[(pfull+1):((pfull+1)+2*nG-1)]
  sigmas <- param[((pfull+1)+2*nG):length(param)]
  Sigma <- Matrix::Matrix(0,nrow=nG,ncol=nG)
  for (i in 1:nG) {
    for (j in 1:nG) {
      if (i<=j) { #Upper Triangular
        Sigma[i,j] <- sigmas[paste("sigma_",i,j,sep="")]
      }
      if (j<i) Sigma[i,j] <- Sigma[j,i] # Lower Triangular
    }
  }
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  deltag <- DELTA[1:nG]
  deltah <- DELTA[(nG+1):(2*nG)]
  lDAg <- Matrix::Matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    #lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    lDAg[i] <- Matrix::determinant(IR-deltag[i]*W,logarithm=TRUE)$modulus
  }
  lDBg <- matrix(0,nrow=length(deltah),ncol=1)
  for (i in 1:nG)
  {
    #lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    lDBg[i] <- Matrix::determinant(IR-deltah[i]*W,logarithm=TRUE)$modulus
  }
  deltaG <- Matrix::Matrix(diag(deltag))
  deltaH <- Matrix::Matrix(diag(deltah))
  A <- kronecker(IT,(IGR - kronecker(deltaG,W)))
  B <- kronecker(IT,(IGR - kronecker(deltaH,W)))
  OME <- kronecker(IT,kronecker(Sigma,IR))
  # BX <- B %*% X
  # Sigma_inv <- Matrix::solve(Sigma)
  # B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
  #                          Matrix::crossprod(BX,Matrix::solve(OME,B%*%(A%*%Y))))
  # rm(BX)
  Res <- A %*% Y - X %*% B_sarar
  BRes <- B %*% Res
  ldet_Sigma <- Matrix::determinant(Sigma,logarithm=TRUE)$modulus
  llike <- as.numeric(-(nT*nG*nR/2)*log(2*pi)
                      - (nR*nT/2)*ldet_Sigma
                      + nT*sum(lDAg) + nT*sum(lDBg)
                      - (1/2)*Matrix::crossprod(BRes,Matrix::solve(OME,BRes)))
  return(-1*llike)
}
