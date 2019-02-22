f_sur_sarar <- function(DELTA,Tm,G,N,Y,X,W,Sigma)
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

LRSUR2_SARAR <- function(Tm, G, N, Y, X, W)
{
  ## PURPOSE:
  # Obtiene la verosimilitud del modelo SUR-SIM y SUR-SARAR y obtiene el LR
  #  ANTES LOS LLAME LMSUR_tres.m
  #---------------------------------------------------
  # where:
  #        Tm   = # of ecuaciones
  #        Y   = Un vector de orden N*Tm*Gx1 con las y apiladas
  #        X  = Un vector N*Tm*GxK X=blkdiag([ones(N,1) Xi]);
  #         info = an (optional) structure variable with input options
  #         info.print = 'yes' / ['no'] print the result

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
  deltah <- matrix(0,nrow=G,ncol=1)
  DELTA <- c(deltag, deltah)
  # Introducir como punto de partida las covarianzas de los residuos del modelo
  # OLS
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
  # Proceso iterativo para la obtencion de los estimadores:
  maxit <- 20;  tol <- 0.5
  opt_sur_sarar <- optim(DELTA,f_sur_sarar,
                         method="BFGS",
                         hessian=FALSE,
                         control=list(fnscale=-1,trace=TRUE),
                         Tm=Tm,G=G,N=N,Y=Y,X=X,
                         W=W,Sigma=Sigma)
  DELTA_t <- opt_sur_sarar$par
  llsur_sarar0 <- opt_sur_sarar$value
  # Control para que no se salga de rango
  DELTA_t <- ifelse(DELTA_t>=1,0.98,DELTA_t)
  for (i in 1:maxit){
    DELTA <- diag(DELTA_t)
    A <- IT %x% (IGR - (DELTA[1:G,1:G] %x% W))
    B <- IT %x% (IGR - (DELTA[(G+1):(2*G),(G+1):(2*G)] %x% W))
    OME <- IT %x% Sigma %x% IR
    BX <- as.matrix(B%*%X)
    Bsarar <- Matrix::Matrix(solve(t(BX)%*%solve(OME,BX),
                                   t(BX)%*%solve(OME,as.vector(A%*%Y))))
    Res <- matrix(B%*%(A%*%Y - X%*%Bsarar),nrow=nrow(Y))
    RR <- array(Res,dim=c(N,G,Tm))
    for (i in 1:G){
      for (j in 1:G){
        Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
      }
    }
    DELTA <- DELTA_t
    opt_sur_sarar <- optim(DELTA,f_sur_sarar,
                           method="BFGS",
                           hessian=FALSE,
                           control=list(fnscale=-1,trace=TRUE),
                           Tm=Tm,G=G,N=N,Y=Y,X=X,
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
  A <- IT %x% (IGR - (DELTA[1:G,1:G] %x% W))
  B <- IT %x% (IGR - (DELTA[(G+1):(2*G),(G+1):(2*G)] %x% W))
  OME <- IT %x% Sigma %x% IR
  BX <- as.matrix(B%*%X)
  B_sarar <- Matrix::Matrix(solve(t(BX)%*%solve(OME,BX),
                                 t(BX)%*%solve(OME,as.vector(A%*%Y))))
  DELTA_sarar <- DELTAfin
  Res_sarar <- matrix(B%*%(A%*%Y - X%*%B_sarar),
                      nrow=nrow(Y))
  Yhat_sarar <- Y - Res_sarar
  RR <- array(Res_sarar,dim=c(N,G,Tm))
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
  ## Obtención de la Lik del modelo SUR-LS
  Sigma_sim <- Matrix::Diagonal(G)
  for (i in 1:2){
    OME_sim <- IT %x% Sigma_sim %x% IR
    B_sim <- solve(t(X)%*%solve(OME_sim,X),
                   t(X)%*%solve(OME_sim,Y))
    Res_sim <- as.matrix(Y - X%*%B_sim)
    RR_sim <- array(Res_sim,dim=c(N,G,Tm))
    for (i in 1:G){
      for (j in 1:G){
        Sigma_sim[i,j] <- cov(matrix(RR_sim[,i,],ncol=1),
                              matrix(RR_sim[,j,],ncol=1))
      }
    }
  }
  IOME_sim <- IT %x% solve(Sigma_sim) %x% IR
  llsur_simfin <- as.numeric( -(Tm*G*N/2)*log(2*pi) -
             (N*Tm/2)*log(det(as.matrix(Sigma_sim))) -
             (1/2)*t(Res_sim) %*% IOME_sim %*% Res_sim)
  test_lr <- -2*(llsur_simfin - llsur_sararfin)

  res <- list(delta_sarar = DELTA_sarar,
              llsur_sarar=llsur_sararfin,
              testlr_sim_sarar = test_lr,
              betas_sarar = B_sarar,
              Sigma_sarar = Sigma,
              Sigma_corr_sarar = Sigma_corr,
              Res_sarar = Res_sarar,
              Yhat_sarar = Yhat_sarar,
              betas_sim = B_sim,
              llsur_sim = llsur_simfin
              )
  #class(res) <- sur_sarar
  return(res)
}

