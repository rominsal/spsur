f_surGRT_beta <- function(beta,Tm,G,N,Y,X,Sigma){
  chol_Sigma <- chol(Sigma)
  ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  Sigma_inv <- try(chol2inv(chol_Sigma))
  IT <- diag(Tm) # as.matrix(Matrix::Diagonal(Tm))
  IR <- diag(N) # as.matrix(Matrix::Diagonal(N))
  #momegap1 <- Matrix::Matrix(kronecker(IT,S0_inv))
  # k0 <- kronecker(Sigma_inv,IR)
  Res <-  (Y - X%*%beta)
  #IOME <- kronecker(IT,kronecker(Sigma_inv,IR))
  #IOMER<- IOME%*%Res
  RES=matrix(as.matrix(Res),nrow = N)
  IOMER<-matrix(RES%*%kronecker(IT,Sigma_inv),ncol = 1)
  # llike <- as.numeric( -(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma - (1/2)*(t(Res) %*% IOME %*% Res))
  llike <- as.numeric( -(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma - (1/2)*Matrix::crossprod(Res,IOMER))
  return((-1)*llike)
}
