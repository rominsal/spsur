#
# f_surGRT_beta <- function(beta,nT,nG,nR,Y,X,Sigma)
# {
#   #ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
#   ldet_Sigma <- determinant(Sigma,logarithm=TRUE)$modulus
#   chol_Sigma <- chol(Sigma)
#   Sigma_inv <- try(chol2inv(chol_Sigma))
#   IT <- Matrix::Diagonal(nT)
#   IR <- Matrix::Diagonal(nR)
#   #momegap1 <- Matrix::Matrix(kronecker(IT,S0_inv))
#   IOME <- kronecker(kronecker(IT,Sigma_inv),IR)
#   Res <-  matrix(Y - X %*% beta,ncol=1)
#   llike <- as.numeric( -(nT*nG*nR/2)*log(2*pi) - (nR*nT/2)*ldet_Sigma
#            - (1/2)*(Matrix::t(Res) %*% IOME %*% Res))
#   llike
# }

sur2_spdiag <- function(nT,nG,nR,Y,X,W,info)
{
  # Estimación del modelo SUR2 sin Efectos Espaciales
  IT <- Matrix::Diagonal(nT)
  IR <- Matrix::Diagonal(nR)
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
  ##
  cat("Computing model whitout spatial effects... \n")
  # Condiciones para la convergencia
  criteria <- 0.0001; itermax <- 100

  for (iter in 1:itermax)
  {
    # opt_sur_sp <- optim(betaoud,f_surGRT_beta,
    #                       method="BFGS",hessian=FALSE,
    #                       control=list(fnscale=-1,trace=FALSE),
    #                       nT=nT,nG=nG,nR=nR,Y=Y,X=X,Sigma=Sigma)
    # beta <- opt_sur_sp$par
    # llsur <- opt_sur_sp$value
    # conv_llsur <- opt_sur_sp$convergence

    opt_sur_sp <- fit_spsursim(nT=nT,nG=nG,nR=nR,Y=Y,X=X,W=W,trace=FALSE)
    beta <- opt_sur_sp$betas
    llsur <- opt_sur_sp$llsur

    Res <- Y - X%*%beta
    RR <- array(Res,dim=c(nR,nG,nT))
    # Covariance residual matrix between G equations
    for (i in 1:nG){
        for (j in 1:nG){
            Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
        }
    }
    dbeta <- sum(abs(beta-betaoud))
    betaoud <- beta
    if (dbeta < criteria) break
  }
  # Test de diagnostico de la dependencia espacial
  # LM_SUR_SLM2
  cat("Computing LM-SAR test... \n")
 E <- array(0,dim=c(nG,nG,nG,nG))
 for(i in 1:nG){
     for(j in 1:nG){
         E[i,j,i,j] <- 1
         E[j,i,i,j] <- 1
     }
 }
 # Gradiente
  g_sar <- rep(0,nG)
  Res <- Y - X%*%beta # Residuos del SUR sin ee
  Sigma_inv <- solve(Sigma)

W<-as(W,"dgCMatrix")

for (i in 1:nG){
      g_sar[i] <- Matrix::t(Res) %*%(IT %x% (Sigma_inv%*%E[,,i,i]) %x% W)%*% Y
  }
  # Matriz de Informacion
  I11 <- Matrix::t(X) %*% (IT %x% Sigma_inv %x% IR) %*% X
  I12 <- matrix(0,nrow=nG,ncol=length(beta))
  for (i in 1:nG){
      I12[i,] <- matrix(t(X) %*% (IT %x% (Sigma_inv%*%E[,,i,i]) %x% W) %*% (X %*% beta),nrow=1)
  }
  I22 <- matrix(0,nrow=nG,ncol=nG)
  OO <- kronecker(IT,kronecker(Sigma,IR)) # Matrix::Matrix(IT %x% Sigma %x% IR)
  tr2 <- sum(W*Matrix::t(W)) # tr <- sum(diag(W%*%W))
  WtW <- Matrix::crossprod(W)
  for (i in 1:nG){
    for (j in 1:nG){
    H <- kronecker(IT,kronecker(E[,,j,j]%*%Sigma_inv%*%E[,,i,i],WtW)) # H <- Matrix::Matrix(IT %x% (E[,,j,j]%*%Sigma_inv%*%E[,,i,i]) %x% WtW)
      tr3 <- sum(H*OO)
      # tr3 <- nT*sum(diag(E[,,j,j]*Sigma_inv*E[,,i,i]%*%Sigma))*sum(diag(WtW))
      if (i==j){
        I22[i,j] <- as.numeric(t(X%*%beta) %*% H %*% (X%*%beta)) + tr3 + nT*tr2
      } else {
        I22[i,j] <- as.numeric(t(X%*%beta) %*% H %*% (X%*%beta)) + tr3
      }
     }
  }

  LMSURSLM2 <- as.numeric(matrix(g_sar,nrow=1) %*% solve(I22-I12%*%solve(I11)%*%t(I12)) %*% matrix(g_sar,ncol=1))



# LM_SUR_SEM2
  cat("Computing LM-SEM test... \n")
    g_sem <- matrix(0,nrow=nG)
    for (i in 1:nG){
        g_sem[i] <- as.numeric(t(Res) %*% (IT %x% (Sigma_inv%*%E[,,i,i]) %x% W) %*% Res)
    }

    J22 <- matrix(0,nrow=nG,ncol=nG)
    tr1 <- sum(Matrix::t(W)*Matrix::t(W))

    for (i in 1:nG){
        for (j in 1:nG){
            J22[i,j]=nT*Sigma_inv[i,j]*Sigma[i,j]*tr1
        }
    }
    J22 <- J22 + nT*tr1*Matrix::Diagonal(nG)
    LMSURSEM2 <- as.numeric(matrix(g_sem,nrow=1) %*% solve(J22) %*% matrix(g_sem,ncol=1))

# LM_SUR_SARMA2
    cat("Computing LM-SARAR test... \n")
    g_sarar1 <- matrix(0,nrow=nG)
    g_sarar2 <- matrix(0,nrow=nG)
    for (i in 1:nG){
        g_sarar1[i] <- t(Res) %*% (IT %x% (Sigma_inv%*%E[,,i,i]) %x% W) %*% Y
        g_sarar2[i] <- t(Res) %*% (IT %x% (Sigma_inv%*%E[,,i,i]) %x% W) %*% Res
    }
    g_sarar <- rbind(g_sarar1, g_sarar2)

# Matriz Informacion
    K11 <- I11
    K12 <- I12
    K22 <- I22
    K33 <- J22
    K23 <- matrix(0,nG,nG)

    for (i in 1:nG){
        for (j in 1:nG){
            K23[i,j] <- nT*Sigma_inv[i,j]*Sigma[i,j]*(tr1+tr2)
        }
    }
    Ird <- matrix(0,nG,nG)
    WW <- W%*%W
    Sigma<-as(Sigma,"dgCMatrix")
    for (g in 1:nG){
        for (s in 1:nG){
            P1s1 <- Sigma %*% E[,,s,s] %*% Sigma_inv%*%E[,,g,g]
            P1s2 <- IT %x% P1s1
            P1 <- P1s2 %x% WtW
            P2s1 <- E[,,g,g] %*% E[,,s,s]
            P2s2 <- IT %x% P2s1
            P2 <- (P2s2 %x% WW)
            Ird[g,s] <- sum(Matrix::diag(P1+P2))
        }
    }
    K23 <- Ird
    LMSURSARAR <- as.numeric(matrix(g_sarar,nrow=1) %*%
               solve(rbind(cbind(K22-K12%*%solve(K11)%*%t(K12), K23),
                            cbind(t(K23), K33))) %*% matrix(g_sarar,ncol=1))

# Test LM*-SUR-Lag
cat("Computing Robust LM*-SUR-SAR test... \n")
# Indexar los elementos diferentes de Sigma
 ff <- rep(0,nG*(nG+1)/2)
 cf <- rep(0,nG*(nG+1)/2)
 c1 <- 0; c2<- 0; c3 <-0
    for (k1 in 1:nG){
        c2<- c2+1
        for (i in 1:(nG-k1+1)){
            c1 <- c1+1; c3 <- c3+1
            ff[c1] <- c2
            cf[c1] <- c3
        }
        c3 <- c2
    }
    RI44 <- matrix(0,nrow=(nG*(nG+1)/2),ncol=nG*(nG+1)/2)
    for (i in 1:(nG*(nG+1)/2)){
        for (j in 1:(nG*(nG+1)/2)){
            RI44[i,j] <- (nT*nR/2)*sum(diag(Sigma_inv%*%E[,,ff[i],cf[i]]
                          %*%Sigma_inv%*%E[,,ff[j],cf[j]]))
        }
    }

    ##F Iff <- rbind(cbind(I11,matrix(0,nrow=(k+1)*nG,ncol=nG*(nG+1)/2)),
    ##F         cbind(matrix(0,nrow=nG*(nG+1)/2,ncol=(k+1)*nG),RI44))
    Iff <- as.matrix(rbind(cbind(I11,matrix(0,nrow=ncol(X),ncol=nG*(nG+1)/2)),
                 cbind(matrix(0,nrow=nG*(nG+1)/2,ncol=ncol(X)),RI44)))
    Ilf <- as.matrix(cbind(I12,matrix(0,nrow=nG,ncol=nG*(nG+1)/2)))
    Ilpf <- I22 - Ilf%*%solve(Iff)%*%t(Ilf)
    Irr <- nT*tr1*(Sigma_inv*Sigma + Matrix::Diagonal(nG))
    Irl <- nT*(tr1+tr2)*(Sigma_inv*Sigma)

LMRSURlag <- as.numeric(Matrix::t(matrix(g_sar,ncol=1) - Irl%*%solve(Irr)%*%matrix(g_sem,ncol=1))
             %*% solve(Ilpf-Irl%*%solve(Irr)%*%Irl)
             %*% (matrix(g_sar,ncol=1)-Irl%*%solve(Irr)%*%matrix(g_sem,ncol=1)))

# Test LM*-SUR-Err
cat("Computing Robust LM*-SUR-SEM test... \n")
LMRSURerr <- as.numeric(Matrix::t(matrix(g_sem,ncol=1)-Irl%*%solve(Ilpf)%*%matrix(g_sar,ncol=1))
             %*% solve(Irr-Irl%*%solve(Ilpf)%*%Irl)
             %*% (matrix(g_sem,ncol=1)-Irl%*%solve(Ilpf)%*%matrix(g_sar,ncol=1)))
res <- list( stat_names= c("LM-SUR-SAR","LM-SUR-SEM",
                          "LM*-SUR-SAR", "LM*-SUR-SEM","LM-SUR-SARAR"),
            stat = c(LMSURSLM2,LMSURSEM2,LMRSURlag, LMRSURerr,LMSURSARAR),
            df = c(2,2,2,2,4) )
}
