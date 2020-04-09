f_sur_sarar <- function(DELTA, Tm, G, N, Y, X, W, Sigma)
{
  IT <- Matrix::Diagonal(Tm)
  IG <- Matrix::Diagonal(G)
  IG <- Matrix::Diagonal(G*N)
  deltag <- DELTA[1:G]
  deltah <- DELTA[(G+1):(2*G)]
  lDAg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    #det_ldAg <- determinant(as.matrix(IR-deltag[i]*W),logarithm=TRUE)
    #lDag[i] <- det_lDAg$modulus*det_lDAg$sign
  }
  deltaG <- diag(deltag)
  lDBg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
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
  llike <- as.numeric(- (Tm*G*N/2)*log(2*pi) -
           (N*Tm/2)*log(det(Sigma)) +
           Tm*sum(lDAg) + Tm*sum(lDBg) -
           (1/2)*t(BRes)%*%IOME%*%BRes)
}

