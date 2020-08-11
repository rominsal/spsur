f_sur_sem <- function(deltag,Tm,G,N,Y,X,W,Sigma)
{
  # Log-lik SUR-SEM Spatio-Temporal Model
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(G*N)
  lDAg <- matrix(0,nrow=length(deltag),ncol=1)
  for (i in 1:G)
  {
    lDAg[i] <- log(det(as.matrix(IR-deltag[i]*W)))
    #det_lDAg <- determinant(as.matrix(IR-deltag[i]*W),logarithm=TRUE)
    #lDAg[i] <- det_lDAg$modulus*det_ldag$sign
  }
  delta <- diag(as.vector(deltag))
  B <-  Matrix::Matrix(IT %x% (IGR - delta %x% W))
  OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
  BX <- as.matrix(B%*%X)
  Bsem <- solve(t(BX)%*%solve(OME,as.matrix(BX)),
                t(BX)%*%solve(OME,as.vector(B%*%Y)))
  rm(BX)
  Res <- Y - X%*%Bsem
  chol_Sigma <- chol(Sigma)
  ldet_Sigma <- sum(2*log(diag(chol_Sigma)))
  BRes <- as.matrix(B%*%Res)
  llike <- as.numeric(-(Tm*G*N/2)*log(2*pi) - (N*Tm/2)*ldet_Sigma
      + Tm*sum(lDAg) - (1/2)*t(BRes)%*%solve(OME,BRes))
}



LMSURquedaSLM <- function(Tm,G,N,Y,X,W)
{
  #### PURPOSE:
  ##  Test LM para contrastar si queda estructura SLM despues de estimar SEM
  ##---------------------------------------------------
  ## where:
  ##        T   = # of ecuaciones
  ##        Y   = Un vector de orden R*T*Gx1 con las y apiladas
  ##        X  = Un vector R*T*GxK X=blkdiag([ones(R,1) Xi]);?
  ##         info = an (optional) structure variable with input options
  ##         info.print = 'yes' / ['no'] print the result
  ##
  ##  SEE ALSO: sur2...
  ## ---------------------------------------------------
  ## REFERENCES:
  ## ---------------------------------------------------
  ## written by:
  ## Fernando A. Lopez Hernandez, 01/05/2010
  ## Dpto. Metodos Cuantitativos e Informaticos
  ## Universidad Politecnica de Cartagena
  ## E-mail: Fernando.Lopez@@upct.es
  ##
  ## Ana Angulo Garijo
  ## Dpto. Economia Aplicada
  ## Universidad de Zaragoza

  #### Test LM queda estructura SLM despues de estimar SEM
  ## Estimacion SURE-SLM para obtener los marginales
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
  maxit <- 20;  tol <- 0.5
  # Obtención del minimo bajo la hip alternativa.
  opt_sur_sem <- optim(deltag,f_sur_sem,
                       #method="L-BFGS-B",lower=-1,upper=1,
                       method="BFGS",
                       hessian=FALSE,
                       control=list(fnscale=-1,trace=TRUE),
                       Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
  deltag_t <- opt_sur_sem$par
  llsur_sem0 <- opt_sur_sem$value
  #conv_llsur_slm <- opt_sur_slm$convergence
  deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)
  # Proceso de estimación iterativo
  for (i in 1:maxit){
    delta <- diag(as.vector(deltag_t))
    B <- Matrix::Matrix(IT %x% (IGR - delta %x% W))
    OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
    BX <- as.matrix(B%*%X)
    Bsem <- solve(t(BX)%*%solve(OME,as.matrix(BX)),
                  t(BX)%*%solve(OME,as.vector(B%*%Y)))
    Res <- matrix(B%*%(Y - X%*%Bsem),nrow=nrow(Y))
    RR <- array(Res,dim=c(N,G,Tm))
    Sigma <- diag(rep(1,G))
    for (i in 1:G){
      for (j in 1:G){
        Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
      }
    }
    deltag <- deltag_t
    opt_sur_sem <- optim(deltag,f_sur_sem,
                         method="BFGS",
                         hessian=FALSE,
                         control=list(fnscale=-1,trace=TRUE),
                         Tm=Tm,G=G,N=N,Y=Y,X=X,W=W,Sigma=Sigma)
    deltag_t <- opt_sur_sem$par
    llsur_sem <- opt_sur_sem$value
    #conv_llsur_slm <- opt_sur_slm$convergence
    deltag_t <- ifelse(deltag_t>=1,0.98,deltag_t)
    if (abs(llsur_sem0 - llsur_sem) < tol) break
    llsur_sem0 <- llsur_sem
  }
  deltafin <- deltag_t
  llsur_semfin <- llsur_sem
  delta <- diag(as.vector(deltafin))
  B <- Matrix::Matrix(IT %x% (IGR - delta %x% W))
  OME <- Matrix::Matrix((IT %x% Sigma) %x% IR)
  #coeficientes finales
  B_sem <- solve(t(BX)%*%solve(OME,as.matrix(BX)),
                t(BX)%*%solve(OME,as.vector(B%*%Y)))
  delta_sem <- deltafin
  Res_sem <- matrix(B%*%(Y - X%*%B_sem),nrow=nrow(Y))
  Yhat_sem <- Y - Res_sem
  RR <- array(Res,dim=c(N,G,Tm))
  for (i in 1:G){
    for (j in 1:G){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),
                        matrix(RR[,j,],ncol=1))
    }
  }
  Sigma_corr <- diag(rep(1,G))
  for (i in 1:G){
    for (j in 1:G){
      Sigma_corr[i,j] <- cor(matrix(RR[,i,],ncol=1),matrix(RR[,j,],ncol=1))
    }
  }
  res <- list(delta_sem = delta_sem,
              betas_sem = B_sem,
              llsur_sem = llsur_semfin,
              Sigma_sem = Sigma,
              Sigma_corr_sem = Sigma_corr,
              Res_sem = Res_sem,
              Yhat_sem = Yhat_sem
  )
  #### Matriz Informacion SUR SEM
  Sigma_inv <- solve(Sigma)
  J12A <- matrix(0,nrow=ncol(X),ncol=G)
  J13A <- matrix(0,nrow=ncol(X),ncol=G*(G+1)/2)
  #IB=sparse(inv(B));
  J22A <- matrix(0,nrow=G,ncol=G)
  J11A <- t(BX) %*% solve(OME,BX)
  WtW <- crossprod(W)
  for (i in 1:G){
    for (j in 1:G){
      if (i==j){
        IBg <- solve(IR - delta_sem[i]*W)
        J22A[i,j] <- Tm*sum(diag(IBg%*%W%*%IBg%*%W)) +
          sum(diag(solve(t(as.matrix(B)),
          as.matrix(IT %x%
          (E[,,j,j]%*%Sigma_inv%*%E[,,i,i] %x% WtW) %*%
          solve(B,as.matrix(OME))))))
      } else {
        J22A[i,j] <-sum(diag(solve(t(as.matrix(B)),
          as.matrix(IT %x%
          (E[,,j,j]%*%Sigma_inv%*%E[,,i,i] %x% WtW) %*%
          solve(B,as.matrix(OME))))))
      }
    }
  }
  J23A <- matrix(0,nrow=G,ncol=G*(G+1)/2)
  for (i in 1:G){
    for (j in 1:(G*(G+1)/2)){
      J23A[i,j] <- sum(diag(solve(t(as.matrix(B)),
          as.matrix(IT %x%
        ((E[,,i,i]%*%Sigma_inv%*%E[,,ff[j],cf[j]]) %x%
              W)))))
    }
  }
  J33A <- matrix(0,nrow=G*(G+1)/2,ncol=G*(G+1)/2)
  for (i in 1:(G*(G+1)/2)){
    for (j in 1:(G*(G+1)/2)){
      J33A[i,j] <- (Tm*N/2)*
        sum(diag((Sigma_inv %*% E[,,ff[i],cf[i]]) %*%
                  (Sigma_inv %*% E[,,ff[j],cf[j]])))
    }
  }
  misem <- rbind(cbind(J11A, J12A, J13A),
           cbind(t(J12A), J22A, J23A),
           cbind(t(J13A), t(J23A), J33A))
  imisem <- solve(misem)
  tmp <- sqrt(diag(imisem))
  tstatsem <- B_sem / tmp[1:ncol(X)]
  tstatTsem <- c(as.matrix(B_sem),as.matrix(delta_sem)) / tmp[1:(ncol(X)+G)]
  res$sd_beta_sem <- tmp[1:ncol(X)]
  res$sd_delta_sem <- tmp[(ncol(X)+1):(ncol(X)+G)]
#### MARGINALES: miro si en el SEM queda estructura SLM. LM_M_SLM
  RT <- matrix(Res,nrow=N*G,ncol=Tm)
  Yt <- matrix(Y,nrow=N*G,ncol=Tm)
  WW <- W%*%W
  Cgradrest <- matrix(0,nrow=Tm,ncol=G)
  for (i in 1:G){
    for (t in 1:Tm){
      Cgradrest[t,i] <- t(RT[,t]) %*%
    (((Sigma_inv %*% E[,,i,i]) %x% W) -
     ((Sigma_inv %*% delta %*% E[,,i,i]) %x% WW)) %*%
     Yt[,t]
    }
  }
  graddel <- colSums(Cgradrest)
  ##matriz de informacion
  PP2 <- matrix(0,nrow=G,ncol=G)
  X_Bsem <- as.matrix(X%*%Bsem)
  for (i in 1:G){
    for (j in 1:G){
      HH <- (IT %x% E[,,j,j] %x% t(W)) %*%
        t(as.matrix(B)) %*%
        (IT %x% Sigma_inv %x% IR) %*% B %*%
        (IT %x% E[,,i,i] %x% W)
      # PP2[i,j] <- as.numeric( t(X%*%Bsem) %*%
      #                           (HH%*%(X%*%Bsem)) +
    #        sum(diag(as.matrix(Binv%*%HH%*%t(Binv)%*%OME))))
      PP2[i,j] <- as.numeric( t(X_Bsem) %*%
                    (HH%*%X_Bsem) +
        sum(diag(solve(B,as.matrix(HH %*%
                    solve(t(as.matrix(B)),
                                as.matrix(OME)))))))
      }
  }
  rm(X_Bsem)
  # NO COINCIDE P22[1,1] CON CÓDIGO MATLAB. REPASAR
  ####
  Ideldel <- Tm*sum(diag(WW))*IG + PP2
  ####
  Ideltabe <- matrix(0,nrow=G,ncol=ncol(X))
  for (i in 1:G)
  {
    Ideltabe[i,] <- as.numeric((t(X%*%Bsem) %*%
      (IT %x% E[,,i,i] %x% t(W))) %*%
      (t(as.matrix(B)) %*% solve(OME,as.matrix(B%*%X))))
  }
  Irhodel <- matrix(0,nrow=G,ncol=G)
  for (i in 1:G){
    for (j in 1:G){
      if(i==j){
        Irhodel[i,j] <-
          sum(diag(solve(B,
as.matrix(OME %*% solve(t(as.matrix(B)),
 as.matrix((IT %x% (E[,,j,j]%*%Sigma_inv) %x% t(W)) %*%
            B %*% (IT %x% E[,,i,i] %x% W))) +
          (IT %x% E[,,i,i] %x% WW)))))
      } else {
        Irhodel[i,j] <-
          sum(diag(solve(B,
  as.matrix(OME %*% solve(t(as.matrix(B)),
    as.matrix((IT %x% (E[,,j,j]%*%Sigma_inv) %x% t(W)) %*%
                   B %*% (IT %x% E[,,i,i] %x% W) ))))))

      }
    }
  }
  Ideltasig <- matrix(0,nrow=G,ncol=G*(G+1)/2)
  for (i in 1:G){
    for (j in 1:(G*(G+1)/2)){
      Ideltasig[i,j] <- sum(diag(solve(B,
        as.matrix((IT %x%
          (E[,,ff[j],cf[j]]%*%Sigma_inv) %x% IR) %*%
        (B %*% (IT %x% E[,,i,i] %x% W))))))
    }
  }
  # OJO: SALE NUMÉRICAMENTE LA MATRIZ NULA... REPASAR
  Idelpsem <- cbind(Ideltabe, t(Irhodel), Ideltasig)
  imi <- solve(Ideldel - Idelpsem %*% solve(misem,t(Idelpsem)))
  LMMqslmensem <- t(graddel) %*% (imi%*%graddel)
  #### Breusch_pagan Test de diagonalidad (Breusch-Pagan 1980)
  ## ver: http://www.stata.com/manuals13/rsureg.pdf
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
  ## Se ajuta a una Chi con G*(G-1)/2 gl
  res$LMM_sem_slm <- LMMqslmensem
  res$BP_sem <- BP
  print(c(BP,LMMqslmensem))
  return(res)

}

