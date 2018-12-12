f_sur_sarar <- function(DELTA, nT, nG, nR, Y, X, W, Sigma)
{
  IT <- Matrix::Diagonal(nT)
  IG <- Matrix::Diagonal(nG)
  IG <- Matrix::Diagonal(nG*nR)
  deltag <- DELTA[1:nG]
  deltah <- DELTA[(nG+1):(2*nG)]
  lDAg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    #det_ldAg <- determinant(as.matrix(IR-deltag[i]*W),logarithm=TRUE)
    #lDag[i] <- det_lDAg$modulus*det_lDAg$sign
  }
  deltaG <- diag(deltag)
  lDBg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:nG)
  {
    lDBg[i] <- log(det(as.matrix(IR-deltah[i]*W)))
    #det_lDBg <- determinant(as.matrix(IR-deltah[i]*W),logarithm=TRUE)
    #lDBg[i] <- det_lDBg$modulus*det_lDBg$sign
  }
  deltaH <- diag(deltah)
  A <- IT %x% (IGR - (deltaG %x% W))
  B <- IT %x% (IGR - (deltaH %x% W))
  Sigma_inv <- solve(Sigma)
  IOME <- (IT %x% Sigma_inv %x% IR)
  BX <- as.matrix(B %*% X)
  Bsarar <- solve(t(BX) %*% IOME %*% BX,
            as.matrix(t(BX)%*% IOME %*% (B %*% (A%*%Y))))
  Res <- A%*%Y - X%*%Bsarar
  BRes <- as.matrix(B%*%Res)
  llike <- as.numeric(- (nT*nG*nR/2)*log(2*pi) -
           (nR*nT/2)*log(det(Sigma)) +
           nT*sum(lDAg) + nT*sum(lDBg) -
           (1/2)*t(BRes)%*%IOME%*%BRes)
}

