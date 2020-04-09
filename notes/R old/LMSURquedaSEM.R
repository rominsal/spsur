f_sur_lag <- function(deltag,Tm,G,N,Y,X,W,Sigma)
{
  # Log-lik SUR-SLM Spatio-Temporal Model
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  lDAg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    #det_lDAg <- determinant(as.matrix(IR-deltag[i]*W),logarithm=TRUE)
    #lDAG[i] <- det_lDAg$modulus*det_lDAg$sign
  }
  delta <- diag(as.vector(deltag))
  A <- Matrix::Matrix(IT %x% (IGR - delta %x% W))
  OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
  # bottleneck here. Try to prevent the computation of inverses...
  Bslm <- Matrix::Matrix(solve(t(X)%*%solve(OME,X),
                t(X)%*%solve(OME,as.vector(A%*%Y))))
  Res <- matrix(A%*%Y - X%*%Bslm,nrow=nrow(Y))
  chol_Sigma <- chol(Sigma)
  ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma
           + Tm*sum(lDAg) - (1/2)*t(Res)%*%solve(OME,Res))
}

LMSURquedaSEM <- function(Tm,G,N,Y,X,W)
{
  ## PURPOSE:
  #  Estima el modelo SUR-SLM y
  #  Test LM para contrastar si queda estructura SEM despues de estimar SLM
  #---------------------------------------------------
  # where:
  #        Tm   = numero de periodos temporales
  #        G   = numero de ecuaciones que correlacionan (sectores)
  #        N   = numero de elementos del corte transversal.
  #        Y   = Un vector de orden R*T*Gx1 con las y apiladas
  #        X  = Un vector R*T*GxK X=blkdiag([ones(R,1) Xi]);
  #
  #  SEE ALSO: LMSURquedaSLM; SUR_SARAR; LRCOMFAC LRSUR2_SARAR
  # Estimación SURE2-SLM para obtener los marginales

  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)

  E <- array(0,dim=c(G,G,G,G))
  for(i in 1:G){
    for(j in 1:G){
      E[i,j,i,j] <- 1
      E[j,i,i,j] <- 1
    }
  }
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0,nrow=G,ncol=1)
  # Introducir como punto de partida las correlaciones del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  RR <- array(Res,dim=c(N,G,Tm))
  Sigma <- diag(rep(1,G))
  for (i in 1:G){
    for (j in 1:G){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }
  # Proceso iterativo para la obtención de los estimadores:
  maxit <- 20;  tol <- 0.5;
  # Obtención del minimo bajo la hip alternativa.
  opt_sur_slm <- optim(deltag,f_sur_lag,
                      method="BFGS",
                      hessian=FALSE,
                      control=list(fnscale=-1,trace=TRUE),
                      Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
  deltag_t <- opt_sur_slm$par
  llsur_slm0 <- opt_sur_slm$value
  #conv_llsur_slm <- opt_sur_slm$convergence
  deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)

  # Proceso de estimacion iterativo
  for(i in 1:maxit)
  {
    delta <- diag(as.vector(deltag_t))
    A <- Matrix::Matrix(IT %x% (IGR - delta %x% W))
    OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
    Bslm <- Matrix::Matrix(solve(t(X)%*%solve(OME,X),
                                 t(X)%*%solve(OME,as.vector(A%*%Y))))
    Res <- matrix(A%*%Y - X%*%Bslm,nrow=nrow(Y))
    RR <- array(Res,dim=c(N,G,Tm))
    for (i in 1:G){
      for (j in 1:G){
        Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
      }
    }
    deltag <- deltag_t
    opt_sur_slm <- optim(deltag,f_sur_lag,
                         #method="L-BFGS-B",lower=-1,upper=1,
                         method="BFGS",
                         hessian=FALSE,
                         control=list(fnscale=-1,trace=TRUE),
                         Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
    deltag_t <- opt_sur_slm$par
    llsur_slm <- opt_sur_slm$value
    #conv_llsur_slm <- opt_sur_slm$convergence
    deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)
    if (abs(llsur_slm0 - llsur_slm) < tol) break
    llsur_slm0 <- llsur_slm
  }
  deltafin <- deltag_t
  delta_slm <- deltafin
  llsur_slmfin <- llsur_slm
  delta <- diag(as.vector(deltafin))
  A <- Matrix::Matrix(IT %x% (IGR - delta %x% W))
  OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
  #coeficientes finales
  B_slm <- Matrix::Matrix(solve(t(X)%*%solve(OME,X),
                   t(X)%*%solve(OME,as.vector(A%*%Y))))
  Res_slm <- matrix(A%*%Y - X%*%B_slm,nrow=nrow(Y))
  Yhat <- Y - Res_slm
  RR <- array(Res,dim=c(N,G,Tm))
  for (i in 1:G){
    for (j in 1:G){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }
  Sigma_corr <- diag(rep(1,G))
  for (i in 1:G){
    for (j in 1:G){
      Sigma_corr[i,j] <- cor(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }
  res <- list(delta_slm = delta_slm,
              betas_slm = B_slm,
              llsur_slm = llsur_slmfin,
              Sigma_slm = Sigma,
              Sigma_corr_slm = Sigma_corr,
              Res_slm = Res_slm,
              Yhat_slm = Yhat
  )
  # SLM MISLM
  # elementos de la matriz
  I11A <- t(X)%*%solve(OME,X)
  #IA=inv(A);
  I12A <- matrix(0,nrow=ncol(X),ncol=G)
  for (i in 1:G){
    I12A[,i] <- as.matrix(t(X) %*% solve(OME,
    as.vector(((IT %x% E[,,i,i]) %x% W)%*%solve(A,as.vector(X%*%Bslm)))))
  }
  fx <- nrow(I12A); cx <- ncol(I12A)
  I13A <- matrix(0,nrow=fx,ncol=(G*(G+1)/2))
  I22A <- matrix(0,G,G)
  chol_Sigma <- chol(Sigma)
  Sigma_inv <- try(chol2inv(chol_Sigma))

  X_Bslm <- as.matrix(X %*% Bslm)
  for (i in 1:G)
  {
    for (j in 1:G)
    {
      HH <- t(as.matrix(IT %x% E[,,j,j] %x% W)) %*%
             (IT %x% Sigma_inv %x% IR) %*%
             (IT %x% E[,,i,i] %x% W)
      if (i==j){
        IAg <- solve(IR - delta_slm[i]*W)
        I22A[i,j] <- as.numeric(Tm*sum(diag(IAg%*%W%*%IAg%*%W)) +
                     t(X_Bslm) %*% (solve(t(as.matrix(A)),
                                          as.matrix(HH)) %*%
                    solve(as.matrix(A),X_Bslm)) +
            sum(diag(solve(t(as.matrix(A)),as.matrix(HH)) %*%
                       solve(as.matrix(A),as.matrix(OME)))))
      } else {
        I22A[i,j] <- as.numeric(t(X_Bslm) %*%
                    (solve(t(as.matrix(A)),as.matrix(HH)) %*%
                       solve(t(as.matrix(A)),X_Bslm)) +
                    sum(diag(solve(t(as.matrix(A)),as.matrix(HH)) %*%
                             solve(t(as.matrix(A)),as.matrix(OME)))))
      }
    }
  }
  rm(HH,X_Bslm)
  # Indexar los elementos diferentes de Sigma
  ff <- rep(0,G*(G+1)/2)
  cf <- rep(0,G*(G+1)/2)
  c1 <- 0;c2 <- 0;c3 <- 0;
  for (k1 in 1:G){
    c2 <- c2+1
    for (i in 1:(G-k1+1)){
      c1 <- c1+1
      c3 <- c3+1;
      ff[c1] <- c2
      cf[c1] <- c3
    }
    c3 <- c2
  }

  I23A <- matrix(0,nrow=G,ncol=G*(G+1)/2)
  #IA <- solve(A)
  for (i in 1:G){
    for (j in 1:G*(G+1)/2){
      #I23A[i,j] <- sum(diag( (IT %x% (E[,,ff[j],cf[j]]%*%Sigma_inv) %x% IR) %*%
      #  (IT %x% E[,,i,i] %x% W) %*% IA ))
      I23A[i,j] <- sum(diag( solve(A,
                as.matrix((IT %x% (E[,,ff[j],cf[j]]%*%Sigma_inv) %x% IR) %*%
                  (IT %x% E[,,i,i] %x% W)))))
    }
  }
  I33A <- matrix(0,nrow=G*(G+1)/2,ncol=G*(G+1)/2)
  for (i in 1:(G*(G+1)/2)){
    for (j in 1:(G*(G+1)/2)){
      I33A[i,j]=(Tm*N/2)* sum(diag( (Sigma_inv%*%E[,,ff[i],cf[i]]) %*%
               (Sigma_inv%*%E[,,ff[j],cf[j]]) ))
    }
  }

  mislm <- rbind(cbind(I11A,I12A,I13A),
           cbind(t(I12A),I22A,I23A),
           cbind(t(I13A),t(I23A),I33A))
  mislm_inv <- solve(mislm)
  tmp <- sqrt(diag(mislm_inv))
  tstatslm <- as.vector(Bslm / tmp[1:fx])
  tstatTslm <- as.numeric(c(as.vector(Bslm), delta_slm)) / tmp[1:(fx+G)]
  res$sd_beta_slm <- tmp[1:ncol(X)]
  res$sd_delta_slm <- tmp[(ncol(X)+1):(ncol(X)+G)]


# MARGINALES: miro si en el slm queda estructura sem. LM_M_SEM
# grad rho
  RT <- matrix(Res,nrow=N*G,ncol=Tm)
  Cgradrest <- matrix(0,nrow=Tm,ncol=G)
  for (i in 1:G){
    for (t in 1:Tm){
      Cgradrest[t,i] <- t(RT[,t]) %*%
                    ((Sigma_inv %*% E[,,i,i]) %x% W) %*% RT[,t]
    }
  }


  if (Tm>1){
  gradrho <- colSums(Cgradrest)
  } else {
  gradrho <- Cgradrest
  }

  #matriz de informacion
  WtW <- crossprod(W)
  WW <- W %*% W
  P1 <- sum(diag(WW)) * IG
  P2 <- sum(diag(WtW)) * (Sigma_inv * Sigma)
  Irhorho <- Tm*(P1+P2)

  Ird <- matrix(0,nrow=G,ncol=G)
  #iA=inv(A);
  for (g in 1:G){
    for (s in 1:G){
      P1s1 <- (Sigma%*%E[,,s,s]) %*% (Sigma_inv %*%E[,,g,g])
      P1 <- ((IT %x% P1s1) %x% WtW)
      P2s1 <- E[,,g,g] %*% E[,,s,s]
      P2 <- ((IT %x% P2s1) %x% WW)
      Ird[g,s] <- sum(diag(solve(A,as.matrix(P1+P2))))
    }
  }
  Irhopslm <- cbind(matrix(0,nrow=G,ncol=fx),
                    Ird, matrix(0,nrow=G,ncol=G*(G+1)/2))

  imi <- solve(Irhorho-Irhopslm%*%mislm_inv%*%t(Irhopslm))
  LMMqsemenslm <- as.numeric(gradrho %*% imi  %*% gradrho)
  # Breusch_pagan Test de diagonalidad (Breusch-Pagan 1980)
  # ver: http://www.stata.com/manuals13/rsureg.pdf
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
  # Se ajusta a una Chi con G*(G-1)/2 gl
  res$LMM_sem_slm <- LMMqsemenslm
  res$BP_slm <- BP
  print(c(BP,LMMqsemenslm))
}



