f_surGRT_beta <- function(beta,nT,nG,nR,Y,X,Sigma)
{
  chol_Sigma <- chol(Sigma)
  ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  Sigma_inv <- try(chol2inv(chol_Sigma))
  IT <- diag(nT) # as.matrix(Matrix::Diagonal(nT))
  IR <- diag(nR) # as.matrix(Matrix::Diagonal(nR))
  
  #momegap1 <- Matrix::Matrix(kronecker(IT,S0_inv))
  k0 <- kronecker(Sigma_inv,IR)
  IOME <- kronecker(IT,k0)
  Res <-  Y - X%*%beta
  llike <- as.numeric( -(nT*nG*nR/2)*log(2*pi) - (nR*nT/2)*ldet_Sigma
                       - (1/2)*(t(Res) %*% IOME %*% Res)) 
  llike
}