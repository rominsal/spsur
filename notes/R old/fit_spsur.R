fit_spsursim <- function(Tm, G, N, Y, X, W,
                         maxit = 100, tol = 0.05,
                         trace = TRUE){
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  #valores iniciales sugerir punto partida a fminsearch
  #deltag <- matrix(0,nrow=G,ncol=1)
  # Introducir como punto de partida las correlaciones del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma;
  rm(Sigmas)
  Sigma_inv <- try(chol2inv(chol(Sigma)))
  if(class(Sigma_inv)=="try-error")
    Sigma_inv <- MASS::ginv(as.matrix(Sigma),tol=1e-40)
    llsur_sim0 <- f_sur_sim(Tm=Tm,G=G,N=N,Y=Y,X=X,
                            Sigma=Sigma)
   if(trace){
    cat("Initial point: ","\n")
    #cat("Sigma: ",round(Sigma,3),"\n")
    cat("log_lik: ",round(-llsur_sim0,3),"\n")
   }
  # Proceso de estimacion iterativo
  for(i in 1:maxit)
  {
    OME <- as(kronecker(IT,kronecker(Sigma,IR)),"dgCMatrix")
    OMEinv <- as(kronecker(IT,kronecker(Sigma_inv,IR)),"dgCMatrix")

    #B_sim <- solve(crossprod(X,solve(OME,X)),
    #                crossprod(X,solve(OME,Y)))
    B_sim <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X),
                          Matrix::crossprod(X,OMEinv %*% Y))
    Res <- matrix(Y - X%*%B_sim,ncol=1)
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Sigmas$Sigma; rm(Sigmas)
    Sigma_inv <- try(chol2inv(chol(Sigma)))
    if(class(Sigma_inv)=="try-error")
      Sigma_inv <- MASS::ginv(as.matrix(Sigma),tol=1e-40)
    llsur_sim <- f_sur_sim(Tm=Tm,G=G,N=N,Y=Y,X=X,Sigma=Sigma)
    if(trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sim,3),"\n")
      #cat("Sigma: ",round(Sigma,3),"\n")
    }
    if (abs(llsur_sim0 - llsur_sim) < tol) break
    llsur_sim0 <- llsur_sim
  }
  llsur_simfin <- llsur_sim
  OME <- as(kronecker(IT,kronecker(Sigma,IR)),"dgCMatrix")
  OMEinv <- as(kronecker(IT,kronecker(Sigma_inv,IR)),"dgCMatrix")
  #coeficientes finales
  #B_sim <- solve(crossprod(X,solve(OME,X)),
  #               crossprod(X,solve(OME,Y)))
  B_sim <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X),
                        Matrix::crossprod(X,OMEinv %*% Y))
  Res <- matrix(Y - X%*%B_sim,nrow=nrow(Y))
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  Sigma_inv <- Sigmas$Sigma_inv
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)
  res <- list(deltas = NULL,
              betas = as.vector(B_sim),
              llsur = -llsur_simfin,
              Sigma = Sigma,
              Sigma_inv = Sigma_inv,
              Sigma_corr = Sigma_corr,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################


fit_spsurslm <- function(Tm, G, N, Y, X, W,
                         maxit = 100, tol = 0.05,
                         trace=TRUE)
{
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0,nrow=G,ncol=1)
  # Introducir como punto de partida las correlaciones del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma; rm(Sigmas)
  # Proceso iterativo para la obtención de los estimadores:
  # Obtención del minimo bajo la hip alternativa.
  # opt_sur_slm <- optim(deltag,f_sur_lag,
  #                      method="BFGS",
  #                      hessian=FALSE,
  #                      #control=list(fnscale=-1),,
  #                      Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
  #
  opt_sur_slm <- minqa::bobyqa(par=deltag,fn=f_sur_lag,
                               lower=rep(-1,length(deltag)),
                               upper=rep(1,length(deltag)),
                               control=list(rhobeg=0.5,iprint=0),
                       Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)

  deltag_t <- opt_sur_slm$par
  llsur_slm0 <- opt_sur_slm$fval
  #deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)
  if(trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_slm0,3)," ")
    cat("lambdas: ",round(deltag_t,3),"\n")

  }
  # Proceso de estimacion iterativo
  for(i in 1:maxit)
  {
    delta <- Matrix::Matrix(diag(as.vector(deltag_t)))
    A <- kronecker(IT,(IGR - kronecker(delta,W)))
    OME <- kronecker(IT,kronecker(Sigma,IR))
    B_slm <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
                           Matrix::crossprod(X,Matrix::solve(OME,A%*%Y)))
    Res <- matrix(A%*%Y - X%*%B_slm,nrow=nrow(Y))
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Sigmas$Sigma; rm(Sigmas)
    deltag <- deltag_t
    # opt_sur_slm <- optim(deltag,f_sur_lag,
    #        method="BFGS",hessian=FALSE,
    #        #control=list(fnscale=-1),,
    #        Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)

    opt_sur_slm <- minqa::bobyqa(par=deltag,fn=f_sur_lag,
                               lower=rep(-1,length(deltag)),
                               upper=rep(1,length(deltag)),
                               control=list(rhobeg=0.5,iprint=0),
                            Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)

    deltag_t <- opt_sur_slm$par
    llsur_slm <- opt_sur_slm$fval
    #deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)
    if(trace){
      cat("Iteration: ",i,"  ")
      cat("log_lik: ",round(-llsur_slm,3)," ")
      cat("lambdas: ",round(deltag_t,3),"\n")

    }
    if (abs(llsur_slm0 - llsur_slm) < tol) break
    llsur_slm0 <- llsur_slm
  }
  deltafin <- deltag_t
  delta_slm <- deltafin
  llsur_slmfin <- llsur_slm
  delta <- Matrix::Matrix(diag(as.vector(deltag_t)))
  A <- kronecker(IT,(IGR - kronecker(delta,W)))

  OME <- kronecker(IT,kronecker(Sigma,IR))
  #coeficientes finales
  B_slm <- Matrix::solve(Matrix::crossprod(X,Matrix::solve(OME,X)),
                         Matrix::crossprod(X,Matrix::solve(OME,A%*%Y)))
  Res <- matrix(A%*%Y - X%*%B_slm,nrow=nrow(Y))
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  Sigma_inv <- Sigmas$Sigma_inv
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)
  res <- list(deltas = as.vector(delta_slm),
              betas = as.vector(B_slm),
              llsur = -llsur_slmfin,
              Sigma = Sigma,
              Sigma_inv = Sigma_inv,
              Sigma_corr = Sigma_corr,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################

fit_spsursem <- function(Tm, G, N, Y, X, W,
                         maxit = 100, tol = 0.05,
                         trace=TRUE)
{
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  ff_cf <- get_ff_cf(G=G)
  ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0,nrow=G,ncol=1)
  # Introducir como punto de partida las correlaciones del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma; rm(Sigmas)
  # Proceso iterativo para la obtención de los estimadores:
  # Obtención del minimo bajo la hip alternativa.
  opt_sur_sem <- minqa::bobyqa(par=deltag,fn=f_sur_sem,
                               lower=rep(-1,length(deltag)),
                               upper=rep(1,length(deltag)),
                               control=list(rhobeg=0.5,iprint=0),
                               Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
  deltag_t <- opt_sur_sem$par
  llsur_sem0 <- opt_sur_sem$fval
   if(trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_sem0,3)," ")
    cat("deltas: ",round(deltag_t,3),"\n")

  }
  # Proceso de estimación iterativo
  for (i in 1:maxit){

    delta <- Matrix::Matrix(diag(as.vector(deltag_t)))
    B <- kronecker(IT,(IGR - kronecker(delta,W)))
    OME <- kronecker(IT,kronecker(Sigma,IR))
    BX <- B%*%X
    B_sem <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                           Matrix::crossprod(BX,Matrix::solve(OME,B%*%Y)))
    Res <- matrix(B%*%(Y - X%*%B_sem),nrow=nrow(Y))
    RR <- array(Res,dim=c(N,G,Tm))
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Sigmas$Sigma; rm(Sigmas)
    deltag <- deltag_t
    opt_sur_sem <- minqa::bobyqa(par=deltag,fn=f_sur_sem,
                                 lower=rep(-1,length(deltag)),
                                 upper=rep(1,length(deltag)),
                                 control=list(rhobeg=0.5,iprint=0),
                                 Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
    deltag_t <- opt_sur_sem$par
    llsur_sem <- opt_sur_sem$fval
    if(trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sem,3)," ")
      cat("rhos: ",round(deltag_t,3),"\n")

    }
    if (abs(llsur_sem0 - llsur_sem) < tol) break
    llsur_sem0 <- llsur_sem
  }
  deltafin <- deltag_t
  llsur_semfin <- llsur_sem
  delta <- Matrix::Matrix(diag(as.vector(deltafin)))
  B <- kronecker(IT,(IGR - kronecker(delta,W)))
  OME <- kronecker(IT,kronecker(Sigma,IR))
  BX <- B%*%X
  #coeficientes finales
  B_sem <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                         Matrix::crossprod(BX,Matrix::solve(OME,B%*%Y)))
  delta_sem <- deltafin
  Res <- matrix(B%*%(Y - X%*%B_sem),nrow=nrow(Y))
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  Sigma_inv<-Sigmas$Sigma_inv
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)
  res <- list(deltas = as.vector(delta_sem),
              betas = as.vector(B_sem),
              llsur = -llsur_semfin,
              Sigma = Sigma,
              Sigma_inv = Sigma_inv,
              Sigma_corr = Sigma_corr,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################

fit_spsursarar <- function(Tm, G, N, Y, X, W,
                           maxit = 100, tol = 0.05,
                           trace=TRUE){
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  ff_cf <- get_ff_cf(G=G)
  ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0,nrow=G,ncol=1)
  deltah <- matrix(0,nrow=G,ncol=1)
  DELTA <- c(deltag, deltah)
  # Introducir como punto de partida las covarianzas de los residuos del
  # modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma; rm(Sigmas)
  # Proceso iterativo para la obtencion de los estimadores:
  opt_sur_sarar <- minqa::bobyqa(par=DELTA,fn=f_sur_sarar,
                               lower=rep(-1,length(DELTA)),
                               upper=rep(1,length(DELTA)),
                               control=list(rhobeg=0.5,iprint=0),
                               Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
  DELTA_t <- opt_sur_sarar$par
  llsur_sarar0 <- opt_sur_sarar$fval
  if(trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_sarar0,3)," ")
    cat("lambdas: ",round(DELTA_t[1:G],3)," ")
    cat("rhos: ",round(DELTA_t[(G+1):(2*G)],3),"\n")
  }
  for (i in 1:maxit){
    DELTA <- Matrix::Matrix(diag(DELTA_t))
    A <- kronecker(IT,(IGR - kronecker(DELTA[1:G,1:G],W)))
    B <- kronecker(IT,(IGR - kronecker(DELTA[(G+1):(2*G),(G+1):(2*G)],W)))
    OME <- kronecker(IT,kronecker(Sigma,IR))
    BX <- B%*%X
    B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                             Matrix::crossprod(BX,Matrix::solve(OME,B %*% (A%*%Y))))
    Res <- matrix(B%*%(A%*%Y - X%*%B_sarar),nrow=nrow(Y))
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Sigmas$Sigma; rm(Sigmas)
    DELTA <- DELTA_t
    opt_sur_sarar <- minqa::bobyqa(par=DELTA,fn=f_sur_sarar,
                                   lower=rep(-1,length(DELTA)),
                                   upper=rep(1,length(DELTA)),
                                   control=list(rhobeg=0.5,iprint=0),
                                   Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
    DELTA_t <- opt_sur_sarar$par
    llsur_sarar <- opt_sur_sarar$fval
    if(trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sarar,3)," ")
      cat("lambdas: ",round(DELTA_t[1:G],3)," ")
      cat("rhos: ",round(DELTA_t[(G+1):(2*G)],3),"\n")
    }
    if (abs(llsur_sarar0 - llsur_sarar) < tol) break
    llsur_sarar0 <- llsur_sarar
  }
  # Coeficientes finales
   DELTA_sarar <- DELTA_t
   llsur_sararfin <- llsur_sarar
   DELTA <- Matrix::Matrix(diag(DELTA_sarar))
   A <- kronecker(IT,(IGR - kronecker(DELTA[1:G,1:G],W)))
   B <- kronecker(IT,(IGR - kronecker(DELTA[(G+1):(2*G),(G+1):(2*G)],W)))
   OME <- kronecker(IT,kronecker(Sigma,IR))
   BX <- B %*% X
   B_sarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                     Matrix::crossprod(BX,Matrix::solve(OME,B %*% (A%*%Y))))
  Res <- matrix(B %*% (A %*% Y - X %*% B_sarar),nrow=nrow(Y))
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  Sigma_inv <- Sigmas$Sigma_inv
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)

  ##### PRUEBAS PARA COVARIANZAS NUMÉRICAS ###########################
  # hessian_optim <- Matrix::Matrix(numDeriv::hessian(
  #                                  func=llsur_sarar_deltas_betas,
  #                                  x=param,#c(B_sarar,DELTA_sarar),
  #                                  p=length(B_sarar),
  #                                  Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma))
  # var_num_sarar <- Matrix::solve(hessian_optim)
  # cat("\n var_numer_sarar",var_num_sarar)
  # cat("\n sd_numer_sarar",sqrt(diag(var_num_sarar)))
###############################################################################

  res <- list(deltas = as.vector(DELTA_sarar),
              betas = as.vector(B_sarar),
              llsur = -llsur_sararfin,
              Sigma = Sigma,
              Sigma_inv = Sigma_inv,
              Sigma_corr = Sigma_corr,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}
