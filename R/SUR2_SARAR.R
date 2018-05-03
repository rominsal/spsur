# setwd("~/Dropbox/SpSUR/spSUR2/R") 
# #setwd("E:/HANDYDRIVE/Artículos en curso/spSUR/SUR2/R")
# load("XXsar.RData")
# load("Ysar.RData")
# load("Ws.RData")
# # 
# X<- as.matrix(XXsar); Y <- as.matrix(Ysar);
# W <- as.matrix(Ws)
# rm(XXsar,Ysar,Ws)
# 
# nT <- 4 # Número de periodos temporales
# nG <- 3 # Número de ecuaciones
# nR <- 49 # Sample size
# k <- (ncol(X)/nG)-1

f_sur_sarar <- function(DELTA,nT,nG,nR,Y,X,W,Sigma)
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

SUR2_SARAR <- function(T,G,R,Y,XX,W)
{
  ## PURPOSE: 
  #  Estima el modelo SARAR
  #---------------------------------------------------
  # where: 
  #        T   = # of ecuaciones
  #        Y   = Un vector de orden R*T*Gx1 con las y apiladas
  #        XX  = Un vector R*T*GxK XX=blkdiag([ones(R,1) Xi]);

  IT <- Matrix::Diagonal(nT) 
  IR <- Matrix::Diagonal(nR)
  IG <- Matrix::Diagonal(nG)
  IGR <- Matrix::Diagonal(nG*nR)
  
  E <- array(0,dim=c(nG,nG,nG,nG))  
  for(i in 1:nG){
    for(j in 1:nG){
      E[i,j,i,j] <- 1
      E[j,i,i,j] <- 1
    }
  }
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0,nrow=nG,ncol=1)
  deltah <- matrix(0,nrow=nG,ncol=1)
  DELTA <- c(deltag, deltah)
  # Introducir como punto de partida las covarianzas de los residuos del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  RR <- array(Res,dim=c(nR,nG,nT))
  Sigma <- diag(rep(1,nG))
  for (i in 1:nG){
    for (j in 1:nG){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }
  # Proceso iterativo para la obtencion de los estimadores:
  maxit <- 20;  tol <- 0.5
  opt_sur_sarar <- optim(DELTA,f_sur_sarar,
                         method="BFGS",
                         hessian=FALSE,
                         control=list(fnscale=-1,trace=TRUE),
                         nT=nT,nG=nG,nR=nR,Y=Y,X=X,
                         W=W,Sigma=Sigma)
  DELTA_t <- opt_sur_sarar$par
  llsur_sarar0 <- opt_sur_sarar$value
  # Control para que no se salga de rango
  DELTA_t <- ifelse(DELTA_t>=1,0.98,DELTA_t)
  for (i in 1:maxit){
    DELTA <- diag(DELTA_t)
    A <- IT %x% (IGR - (DELTA[1:nG,1:nG] %x% W))
    B <- IT %x% (IGR - (DELTA[(nG+1):(2*nG),(nG+1):(2*nG)] %x% W))
    OME <- IT %x% Sigma %x% IR
    BX <- as.matrix(B%*%X)
    Bsarar <- Matrix::Matrix(
      solve(t(BX)%*%solve(OME,BX),
         t(BX)%*%solve(OME,as.vector(A%*%Y))))
    Res <- matrix(B%*%(A%*%Y - X%*%Bsarar),nrow=nrow(Y))
    RR <- array(Res,dim=c(nR,nG,nT))
    for (i in 1:nG){
      for (j in 1:nG){
        Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
      }
    } 
    DELTA <- DELTA_t
    opt_sur_sarar <- optim(DELTA,f_sur_sarar,
                           method="BFGS",
                           hessian=FALSE,
                           control=list(fnscale=-1,trace=TRUE),
                           nT=nT,nG=nG,nR=nR,Y=Y,X=X,
                           W=W,Sigma=Sigma)
    DELTA_t <- opt_sur_sarar$par
    llsur_sarar <- opt_sur_sarar$value
    DELTA_t <- ifelse(DELTA_t>=1,0.98,DELTA_t)
    if (abs(llsur_sarar0 - llsur_sarar) < tol) break
    llsur_sarar0 <- llsur_sarar
  }  
  # Coeficientes finales
  DELTAfin <- DELTA_t
  llsur_sararfin <- llsur_sarar
  DELTA <- diag(DELTAfin)
  A <- IT %x% (IGR - (DELTA[1:nG,1:nG] %x% W))
  B <- IT %x% (IGR - (DELTA[(nG+1):(2*nG),(nG+1):(2*nG)] %x% W))
  OME <- IT %x% Sigma %x% IR
  DELTA_sarar <- DELTAfin
  BX <- as.matrix(B%*%X)
  B_sarar <- Matrix::Matrix(
    solve(t(BX)%*%solve(OME,BX),
          t(BX)%*%solve(OME,as.vector(A%*%Y))))
  Res_sarar <- matrix(B%*%(A%*%Y - X%*%B_sarar),
                      nrow=nrow(Y))
  Yhat_sarar <- Y - Res_sarar
  RR <- array(Res_sarar,dim=c(nR,nG,nT))
  for (i in 1:nG){
    for (j in 1:nG){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }   
  Sigma_corr <- diag(rep(1,nG))
  for (i in 1:nG){
    for (j in 1:nG){
      Sigma_corr[i,j] <- cor(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }  

 
  
## MI SARAR
  IOME <- IT %x% solve(Sigma) %x% IR
  BX <- as.matrix(B%*%X)
  XXTBTIOMEB <- t(BX) %*% IOME %*% B
  J11A <- Matrix::Matrix(XXTBTIOMEB%*%X)
  #IA=inv(A);
  #IB=inv(B);
  #iSH1=inv(SH1);
  delta_sar <- DELTA_sarar[1:nG]
  delta_sem <- DELTA_sarar[(nG+1):(2*nG)]
  J12A <- matrix(0,ncol(X),nG) 
  for (i in 1:nG){
    J12A[,i] <- matrix(XXTBTIOMEB %*% (IT %x% E[,,i,i] %x% W) %*% 
               solve(A,as.matrix(X %*% B_sarar)),ncol=1)    
  }
  rm(XXTBTIOMEB)
  fx <- nrow(J12A)
  cx <- ncol(J12A)
  J13A <- matrix(0,nrow=fx,ncol=nG)
  J14A <- matrix(0,nrow=fx,ncol=nG*(nG+1)/2)
  J21A <- t(J12A)
  J22A <- matrix(0,nrow=nG,ncol=nG)
  XXBsarar <- as.matrix(X%*%B_sarar)
  Sigma_inv <- solve(Sigma)
  BF <- t(as.matrix(B)) %*% ((IT %x% Sigma_inv %x% IR) %*% B)
  for (i in 1:nG){
    for (j in i:nG){
      HH <- t(as.matrix(IT %x% E[,,j,j] %x% W)) %*% 
        BF %*% (IT %x% E[,,i,i] %x% W)  
      t1 <- sum(diag(solve(t(as.matrix(A))%*%B,as.matrix(HH)) %*% 
        solve(t(as.matrix(B))%*%A,as.matrix(OME))))
      #HH =IA'* kronF(kronF(sIT,E(j,j).eq),W)'*
      #  BF*kronF(kronF(sIT,E(i,i).eq),W)*IA
      #t1 <- sum(sum(HH.*BG'))
      if (i==j){
        IAg <- solve(IR - delta_sar[i]*W)
        BF2 <- IAg %*% W
        J22A[i,j] <- nT*sum(BF2*t(BF2)) +
          t(XXBsarar) %*% 
          solve(t(as.matrix(A)),as.matrix(HH)) %*% 
          solve(A,XXBsarar) + t1
      } else {
        if (i>j){
          J22A[i,j] <-J22A[j,i]
        } else {
          J22A[i,j] <- t(XXBsarar) %*% 
            solve(t(as.matrix(A)),as.matrix(HH)) %*% 
            solve(A,XXBsarar) + t1
        }  
      }
    }
  }
  rm(HH,IAg,BF,XXBsarar)
  J23A <- matrix(0,nrow=nG,ncol=nG)
  #IAIB <- IA*IB;
  # OMEIBT=OME*IB';
  WW <- W%*%W
  for (i in 1:nG){
    for (j in 1:nG){
      HH <- OME %*% solve(t(as.matrix(B)),
        as.matrix((IT %x% (E[,,i,i]%*%Sigma_inv) %x% t(W)) %*%
        B %*% (IT %x% E[,,j,j] %x% W) +
        (IT %x% (E[,,i,i] %*% E[,,j,j]) %x% W)))
      J23A[i,j] = sum(diag(solve(B%*%A,as.matrix(HH))))
    }
  }
  rm(HH)
# Indexar los elementos diferentes de Sigma
  c1 <- c2 <- c3 <- 0
  ff <- rep(0,nG*(nG+1)/2)
  cf <- rep(0,nG*(nG+1)/2)  
  for (k1 in 1:nG){
    c2 <- c2+1
    for (i in 1:(nG-k1+1)){
      c1 <- c1+1; c3 <- c3+1
      ff[c1] <- c2
      cf[c1] <- c3
    }
    c3 <- c2
  }
  J24A <- matrix(0,nrow=nG,ncol=nG*(nG+1)/2)
  for (i in 1:nG){
    for (j in 1:(nG*(nG+1)/2)) {
      J24A[i,j] <- sum(diag(
        solve(B%*%A,
        as.matrix((IT %x% 
        (E[,,ff[j],cf[j]] %*% Sigma_inv) %x% IR) %*%
        B  %*% (IT %x% E[,,i,i] %x% W)))))
     }
  }
  J33A <- matrix(0,nrow=nG,ncol=nG)
  WtW = crossprod(W)
  #FFFF=IB*OME*IB';
  for (i in 1:nG){
    for (j in 1:nG){
      if (i==j){
        Bg <- IR - delta_sem[i]*W
        BG <- solve(Bg,W)
        BFF <- solve(t(as.matrix(B)),
            as.matrix(IT %x% 
          (E[,,j,j]%*%Sigma_inv%*%E[,,i,i]) %x% WtW)) 
        J33A[i,j] <- nT*sum(BG*t(BG)) +
           sum(BFF*t(as.matrix(solve(B,as.matrix(OME)))))
      } else {
        if (i>j) {
          J33A[i,j] <- J33A[j,i] 
        } else {
          BG <- IT %x% (E[,,j,j]%*%Sigma_inv%*%E[,,i,i]) %x%
                WtW
          J33A[i,j] <- sum(as.matrix(
            solve(t(as.matrix(B)),as.matrix(BG)) *
            t(solve(B,as.matrix(OME)))))  
        }
      }
    }
  }
  rm(Bg,BG,BFF) 

  J34A <- matrix(0,nrow=nG,ncol=nG*(nG+1)/2)
  for (i in 1:nG){
    for (j in 1:(nG*(nG+1)/2)){
      J34A[i,j] <- sum(diag(solve(t(as.matrix(B)),
        as.matrix((IT %x% 
        (E[,,i,i] %*% Sigma_inv %*% E[,,ff[j],cf[j]]) %x% 
          W)))))
    }
  }
  J44A <- matrix(0,nrow=nG*(nG+1)/2,ncol=nG*(nG+1)/2)
  for (i in 1:(nG*(nG+1)/2)){
    for (j in 1:(nG*(nG+1)/2)){
      if (i>j){
        J44A[i,j] <- J44A[j,i]
      } else {
        J44A[i,j] <-
          sum(((nT*nR/2)*Sigma_inv %*% E[,,ff[i],cf[i]]) *
          t(Sigma_inv%*%E[,,ff[j],cf[j]]))
      }
    }
  }
  misml <- rbind(cbind(J11A, J12A, J13A, J14A),
           cbind(t(J12A), J22A, J23A, J24A),
           cbind(t(J13A), t(J23A), J33A, J34A),
           cbind(t(J14A), t(J24A), t(J34A), J44A))
  imislm <- solve(misml)
  tmp <- sqrt(diag(imislm))
  tstatsarar <- B_sarar / tmp[1:fx]  
  tstatTsarar <- c(matrix(B_sarar,ncol=1),
                   delta_sar,delta_sem) /tmp[1:(fx+2*nG)]
  # CAMBIAN LAS ESTIMACIONES DE LOS INTERCEPTOS (EN UN
  # PAR DE DÉCIMAS) RESPECTO A CÓDIGO MATLAB
  index_ltri <- lower.tri(Sigma_corr)
  BP <- nR*nT*sum(Sigma_corr[index_ltri]^2)

}
