fit_spsursim <- function(env, con){
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  #valores iniciales sugerir punto partida a fminsearch
  # Introducir como punto de partida las correlaciones del modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Matrix::Matrix(Sigmas$Sigma);
  rm(Sigmas)
  Sigmainv <- Matrix::solve(Sigma)
  assign("Sigma", Sigma, envir = env)
  llsur_sim0 <- f_sur_sim(env = env)
   if (con$trace){
    cat("Initial point: ","\n")
    cat("log_lik: ",round(-llsur_sim0,3),"\n")
   }
  # Proceso de estimacion iterativo
  for(i in 1:con$maxit)
  {
    OMEinv <- Matrix::kronecker(IT,Matrix::kronecker(Sigmainv,IR))
    B_sim <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                          Matrix::crossprod(X, OMEinv %*% Y))
    Res <- matrix(Y - X%*%B_sim, ncol=1)
    Sigmas <- get_Sigma(resids=Res, N=N, G=G, Tm=Tm)
    Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
    assign("Sigma", Sigma, envir = env)
    llsur_sim <- f_sur_sim(env = env)
    if(con$trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sim,3),"\n")
    }
    if (abs(llsur_sim0 - llsur_sim) < con$tol) break
    llsur_sim0 <- llsur_sim
  }
  llsur_simfin <- llsur_sim
  Sigmainv <- Matrix::solve(Sigma)
  OMEinv <- Matrix::kronecker(IT, Matrix::kronecker(Sigmainv,IR))
  B_sim <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                        Matrix::crossprod(X, OMEinv %*% Y))
  Res <- matrix(Y - X%*%B_sim, ncol=1)
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  rm(Sigmas)
  res <- list(deltas = NULL,
              coefficients = as.vector(B_sim),
              LL = -llsur_simfin,
              Sigma = Sigma,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################


fit_spsurslm <- function(env, con) {
  if(!is.null(env$W)) W <- env$W else W <- as(env$listw, "CsparseMatrix")
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- rep(0, G)
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res, N=N, G=G, Tm=Tm)
  Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
  assign("Sigma", Sigma, envir = env)
  Sigmainv <- Matrix::solve(Sigma)
  # Proceso iterativo para la obtención de los estimadores:
  # Obtención del minimo bajo la hip alternativa.
  opt_sur_slm <- minqa::bobyqa(par = deltag, fn = f_sur_lag,
                               lower = rep(env$interval[1],length(deltag)),
                               upper = rep(env$interval[2],length(deltag)),
                               control = list(rhobeg=0.5,iprint=0),
                               env = env)
  deltag_t <- opt_sur_slm$par
  llsur_slm0 <- opt_sur_slm$fval
  if (con$trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_slm0,3)," ")
    cat("rhos: ",round(deltag_t,3),"\n")

  }
  # Proceso de estimacion iterativo
  for(i in 1:con$maxit) {
    delta <- Matrix::Diagonal(length(deltag_t),deltag_t)
    AY <- (IGR - kronecker(delta, W)) %*% Y
    OMEinv <- kronecker(kronecker(IT,Sigmainv),IR)
    B_slm <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                          Matrix::crossprod(X, OMEinv %*% AY))
    Res <- matrix(AY - X %*% B_slm, ncol = 1)
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
    Sigmainv <- Matrix::solve(Sigma)
    assign("Sigma", Sigma, envir = env)
    deltag <- deltag_t
    opt_sur_slm <- minqa::bobyqa(par=deltag,fn=f_sur_lag,
                               lower=rep(env$interval[1],length(deltag)),
                               upper=rep(env$interval[2],length(deltag)),
                               control=list(rhobeg=0.5,iprint=0),
                               env = env)
    deltag_t <- opt_sur_slm$par
    llsur_slm <- opt_sur_slm$fval
    if(con$trace){
      cat("Iteration: ",i,"  ")
      cat("log_lik: ",round(-llsur_slm,3)," ")
      cat("rhos: ",round(deltag_t,3),"\n")

    }
    if (abs(llsur_slm0 - llsur_slm) < con$tol) break
    llsur_slm0 <- llsur_slm
  }
  deltafin <- deltag_t
  delta_slm <- deltafin
  llsur_slmfin <- llsur_slm
  # #coeficientes finales
  delta <- Matrix::Matrix(diag(x=deltag_t, nrow = length(deltag_t),
                       ncol = length(deltag_t)))
  AY <- (IGR - kronecker(delta, W)) %*% Y
  OMEinv <- kronecker(kronecker(IT,Sigmainv),IR)
  B_slm <- Matrix::solve(Matrix::crossprod(X, OMEinv %*% X),
                         Matrix::crossprod(X, OMEinv %*% AY))
  Res <- matrix(AY - X %*% B_slm, ncol = 1 )
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res, N=N, G=G, Tm=Tm)
  Sigma <- Sigmas$Sigma
  rm(Sigmas)
  res <- list(deltas = as.vector(delta_slm),
              coefficients = as.vector(B_slm),
              LL = -llsur_slmfin,
              Sigma = Sigma,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################

fit_spsursem <- function(env, con) {
  if(!is.null(env$W)) W <- env$W else W <- as(env$listw, "CsparseMatrix")
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
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
  Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
  assign("Sigma", Sigma, envir = env)
  Sigmainv <- Matrix::solve(Sigma)
  # Proceso iterativo para la obtención de los estimadores:
  # Obtención del minimo bajo la hip alternativa.
  opt_sur_sem <- minqa::bobyqa(par = deltag,fn=f_sur_sem,
                               lower = rep(env$interval[1],length(deltag)),
                               upper = rep(env$interval[2],length(deltag)),
                               control = list(rhobeg=0.5,iprint=0),
                               env = env)
  deltag_t <- opt_sur_sem$par
  llsur_sem0 <- opt_sur_sem$fval
   if(con$trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_sem0,3)," ")
    cat("lambdas: ",round(deltag_t,3),"\n")
  }
  # Proceso de estimación iterativo
  for (i in 1:con$maxit){
    delta <- Matrix::Diagonal(length(deltag_t),deltag_t)
    B <- kronecker(IT,(IGR - kronecker(delta,W)))
    OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
    BX <- kronecker(IT,(IGR - kronecker(delta,W))) %*% X
    BY <- kronecker(IT,(IGR - kronecker(delta,W))) %*% Y
    B_sem <- Matrix::solve(Matrix::crossprod(BX, OMEinv %*% BX),
                           Matrix::crossprod(BX, OMEinv %*% BY))
    Res <- matrix(BY - BX %*% B_sem, ncol=1)
    RR <- array(Res, dim=c(N,G,Tm))
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
    assign("Sigma", Sigma, envir = env)
    Sigmainv <- Matrix::solve(Sigma)
    deltag <- deltag_t
    opt_sur_sem <- minqa::bobyqa(par=deltag,fn=f_sur_sem,
                                 lower = rep(env$interval[1],length(deltag)),
                                 upper = rep(env$interval[2],length(deltag)),
                                 control = list(rhobeg=0.5,iprint=0),
                                 env = env)
    deltag_t <- opt_sur_sem$par
    llsur_sem <- opt_sur_sem$fval
    if(con$trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sem,3)," ")
      cat("lambdas: ",round(deltag_t,3),"\n")

    }
    if (abs(llsur_sem0 - llsur_sem) < con$tol) break
    llsur_sem0 <- llsur_sem
  }
  deltafin <- deltag_t
  llsur_semfin <- llsur_sem
  delta <- Matrix::Diagonal(length(deltafin),deltafin)
  OMEinv <- Matrix::kronecker(IT,Matrix::kronecker(Sigmainv,IR))
  BX <- Matrix::kronecker(IT,(IGR - Matrix::kronecker(delta,W))) %*% X
  BY <- Matrix::kronecker(IT,(IGR - Matrix::kronecker(delta,W))) %*% Y
  #coeficientes finales
  B_sem <- Matrix::solve(Matrix::crossprod(BX, OMEinv %*% BX),
                         Matrix::crossprod(BX, OMEinv %*% BY))
  delta_sem <- deltafin
  #Res <- matrix(B%*%(Y - X%*%B_sem),nrow=nrow(Y))
  Res <- matrix(BY - BX %*% B_sem, ncol=1)
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Sigmas$Sigma
  rm(Sigmas)
  res <- list(deltas = as.vector(delta_sem),
              coefficients = as.vector(B_sem),
              LL = -llsur_semfin,
              Sigma = Sigma,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}

###############################################

fit_spsursarar <- function(env, con) {
  if(!is.null(env$W)) W <- env$W else W <- as(env$listw, "CsparseMatrix")
  G <- env$G; N <- env$N; Tm <- env$Tm
  Y <- env$Y; X <- env$X
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(G*N)
  E <- get_array_E(G)
  ff_cf <- get_ff_cf(G=G)
  ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
  #valores iniciales sugerir punto partida a fminsearch
  deltag <- matrix(0, nrow=G, ncol=1)
  deltah <- matrix(0, nrow=G, ncol=1)
  DELTA <- c(deltag, deltah)
  # Introducir como punto de partida las covarianzas de los residuos del
  # modelo OLS
  ols_init <- lm(Y ~ X - 1)
  betaoud  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
  assign("Sigma", Sigma, envir = env)
  Sigmainv <- Matrix::solve(Sigma)
  # Proceso iterativo para la obtencion de los estimadores:
  opt_sur_sarar <- minqa::bobyqa(par=DELTA,fn=f_sur_sarar,
                                 lower = rep(env$interval[1],length(DELTA)),
                                 upper = rep(env$interval[2],length(DELTA)),
                                 control = list(rhobeg=0.5,iprint=0),
                                 env = env)
  DELTA_t <- opt_sur_sarar$par
  llsur_sarar0 <- opt_sur_sarar$fval
  if(con$trace){
    cat("Initial point: "," ")
    cat("log_lik: ",round(-llsur_sarar0,3)," ")
    cat("rhos: ",round(DELTA_t[1:G],3)," ")
    cat("lambdas: ",round(DELTA_t[(G+1):(2*G)],3),"\n")
  }
  for (i in 1:con$maxit){
    DELTA <- Matrix::Diagonal(length(DELTA_t), DELTA_t)
    #DELTA <- Matrix::Matrix(diag(DELTA_t))
    AY <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[1:G,1:G],W))) %*% Y
    BX <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[(G+1):(2*G),
                                        (G+1):(2*G)],W))) %*% X
    BAY <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[(G+1):(2*G),
                                        (G+1):(2*G)],W))) %*% AY
    OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
    B_sarar <- Matrix::solve(Matrix::crossprod(BX, OMEinv %*% BX),
                             Matrix::crossprod(BX, OMEinv %*% BAY))
    Res <- matrix(BAY - (BX %*% B_sarar), ncol=1)
    Sigmas <- get_Sigma(resids=Res, N=N, G=G, Tm=Tm)
    Sigma <- Matrix::Matrix(Sigmas$Sigma); rm(Sigmas)
    assign("Sigma", Sigma, envir = env)
    Sigmainv <- Matrix::solve(Sigma)
    DELTA <- DELTA_t
    opt_sur_sarar <- minqa::bobyqa(par=DELTA, fn=f_sur_sarar,
                                   lower = rep(env$interval[1],length(DELTA)),
                                   upper = rep(env$interval[2],length(DELTA)),
                                   control = list(rhobeg=0.5,iprint=0),
                                   env = env)
    DELTA_t <- opt_sur_sarar$par
    llsur_sarar <- opt_sur_sarar$fval
    if(con$trace){
      cat("Iteration: ",i," ")
      cat("log_lik: ",round(-llsur_sarar,3)," ")
      cat("rhos: ",round(DELTA_t[1:G],3)," ")
      cat("lambdas: ",round(DELTA_t[(G+1):(2*G)],3),"\n")
    }
    if (abs(llsur_sarar0 - llsur_sarar) < con$tol) break
    llsur_sarar0 <- llsur_sarar
  }
  # Coeficientes finales
  DELTA_sarar <- DELTA_t
  llsur_sararfin <- llsur_sarar
  DELTA <- Matrix::Diagonal(length(DELTA_sarar), DELTA_sarar)
  AY <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[1:G,1:G],W))) %*% Y
  BX <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[(G+1):(2*G),
                                      (G+1):(2*G)],W))) %*% X
  BAY <- Matrix::kronecker(IT, (IGR - Matrix::kronecker(DELTA[(G+1):(2*G),
                                      (G+1):(2*G)],W))) %*% AY
  OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
  B_sarar <- Matrix::solve(Matrix::crossprod(BX, OMEinv %*% BX),
                           Matrix::crossprod(BX, OMEinv %*% BAY))
  Res <- matrix(BAY - (BX %*% B_sarar), ncol=1)
  Yhat <- Y - Res
  Sigmas <- get_Sigma(resids=Res, N=N, G=G, Tm=Tm)
  Sigma <- Sigmas$Sigma
  rm(Sigmas)
  res <- list(deltas = as.vector(DELTA_sarar),
              coefficients = as.vector(B_sarar),
              LL = -llsur_sararfin,
              Sigma = Sigma,
              residuals = as.vector(Res),
              fitted.values = as.vector(Yhat)
  )
}
