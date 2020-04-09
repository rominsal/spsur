cov_spsursim <- function(Tm,G,N,Y,X,W,
                         deltas=NULL,Sigma,cov=FALSE,trace=FALSE){
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)

  E <- get_array_E(G)
  Sigmainv <- try(chol2inv(chol(Sigma)))
  if(class(Sigmainv)=="try-error")
    Sigmainv <- MASS::ginv(as.matrix(Sigma), tol=1e-40)
  #delta <- Matrix::Matrix(diag(as.vector(deltas)))
  #A <- kronecker(IT,(IGR - kronecker(delta,W)))
  #OME <- kronecker(IT,kronecker(Sigma,IR))
  OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
  Bsim <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X),
                         Matrix::crossprod(X,OMEinv %*% Y))
  Res <- matrix(Y - X %*% Bsim, nrow=nrow(Y))
  misim_inv <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X))
  tmp <- sqrt(Matrix::diag(misim_inv))
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
  # Se ajusta a una Chi con G*(G-1)/2 gl
  res <- list(
    se_betas = as.vector(tmp[1:ncol(X)]),
    cov = as.matrix(misim_inv),
    BP = as.numeric(BP)
  )

}


###############################################
  cov_spsurslm <- function(Tm,G,N,Y,X,W,
                         deltas,Sigma,trace=FALSE){
    W <- as(W,"dgCMatrix")
    IT <- Matrix::Diagonal(Tm)
    IR <- Matrix::Diagonal(N)
    IG <- Matrix::Diagonal(G)
    IGR <- Matrix::Diagonal(G*N)
    WtW <- as(Matrix::crossprod(W),"dgCMatrix")
    WW <- as(W %*% W,"dgCMatrix")
    Wt <- as(Matrix::t(W),"dgCMatrix")

    E <- get_array_E(G)
    Sigma <- Matrix::Matrix(Sigma)
    Sigmainv <- try(Matrix::solve(Sigma))
    #Sigmainv <- try(chol2inv(chol(Sigma)))
    if(class(Sigmainv)=="try-error")
      Sigmainv <- MASS::ginv(as.matrix(Sigma),tol=1e-40)
    delta <- Matrix::Matrix(diag(as.vector(deltas)))
    # Auxiliar Matrices
    Aux <- IGR-kronecker(delta,W)
    Auxt <- Matrix::t(Aux)
    Auxi <- Matrix::solve(Aux)
    Auxit<-Matrix::t(Auxi)
    A <- kronecker(IT,Aux)
    Ainv <- kronecker(IT,Auxi)
    rm(Auxi,Aux)
    OME <- kronecker(IT,kronecker(Sigma,IR))
    OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
    AY <- A %*% Y
    rm(A)
    Bslm <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X),
                          Matrix::crossprod(X,OMEinv %*% AY) )
    XBslm <- X %*% Bslm
    Res <- matrix(AY - XBslm,nrow=nrow(Y))
    ################################
    #### Information Matrix SUR-SLM
    ################################

    ## I11A
    tX_OMEinv <- Matrix::crossprod(X,OMEinv)
    I11A <- tX_OMEinv %*% X
    ## I12A
    I12A <- Matrix::Matrix(0,nrow=ncol(X),ncol=G)
    Ainv_XBslm <-  Ainv %*% XBslm
    rm(Ainv)
    for (i in 1:G){
      I12A[,i] <- tX_OMEinv %*%
        (kronecker(kronecker(IT,E[,,i,i]),W) %*% Ainv_XBslm)
    }
    rm(tX_OMEinv)
    ## I13A
    fx <- nrow(I12A); cx <- ncol(I12A)
    I13A <- Matrix::Matrix(0,nrow=fx,ncol=(G*(G+1)/2))
    # ## I22A
    # I22A <- Matrix::Matrix(0,nrow=G,ncol=G)
    # for (i in 1:G){
    #   IAgW <- Matrix::solve((IR - deltas[i]*W),W)
    #   for (j in 1:G){
    #     kk1 <- kronecker(IT,kronecker(E[,,i,i],W))%*%Ainv_XBslm
    #     kk1 <- kronecker(IT,kronecker(Sigmainv,IR))%*%kk1
    #     kk1 <- kronecker(IT,kronecker(E[,,j,j],Matrix::t(W)))%*%kk1
    #     kk1 <- Matrix::solve(kronecker(IT,Auxt),kk1)
    #     tXBslm_HH_XBslm <- Matrix::crossprod(XBslm,kk1)
    #     rm(kk1)
    #     k1 <- Auxit%*%kronecker(E[,,j,j]%*%Sigmainv%*%E[,,i,i],WtW)
    #     k2 <- kronecker(Sigma,IR)%*%Auxit
    #     trace_HH_OME <- Tm*sum(k1*k2)
    #     if (i==j){
    #       I22A[i,j] <- Tm*sum(IAgW * Matrix::t(IAgW)) + tXBslm_HH_XBslm +
    #         trace_HH_OME
    #     } else {
    #       I22A[i,j] <- tXBslm_HH_XBslm + trace_HH_OME
    #     }
    #   }
    # }
    # rm(k1,k2)

    ## I22A
    I22A <- Matrix::Matrix(0,nrow=G,ncol=G)
    k2 <- kronecker(Sigma,IR)%*%Auxit
    for (i in 1:G){
      IAgW <- Matrix::solve((IR - deltas[i]*W),W)
      for (j in 1:G){
        if (i>j){I22A[i,j]<-I22A[j,i]}
        else {
          # kk1 <- kronecker(IT,Auxit%*%kronecker(E[,,j,j]%*%Sigmainv%*%E[,,i,i],WtW))%*%Ainv_XBslm
              kk1 <- kronecker(IT,kronecker(E[,,i,i],W))%*%Ainv_XBslm
              kk1 <- kronecker(IT,kronecker(Sigmainv,IR))%*%kk1
              kk1 <- kronecker(IT,kronecker(E[,,j,j],Matrix::t(W)))%*%kk1
              kk1 <- Matrix::solve(kronecker(IT,Auxt),kk1)
          tXBslm_HH_XBslm <- Matrix::crossprod(XBslm,kk1)
          rm(kk1)
          k1 <- Auxit%*%kronecker(E[,,j,j]%*%Sigmainv%*%E[,,i,i],WtW)
          trace_HH_OME <- Tm*sum(k1*k2)
          if (i==j){
            I22A[i,j] <- Tm*sum(IAgW * Matrix::t(IAgW)) + tXBslm_HH_XBslm +
              trace_HH_OME
          } else {
            I22A[i,j] <- tXBslm_HH_XBslm + trace_HH_OME
          }
        }
      }
    }
    rm(k1,k2)
    ff_cf <- get_ff_cf(G=G)
    ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
    ## I23A
    # I23A <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
    # for (i in 1:G){
    #   for (j in 1:(G*(G+1)/2)){
    #     if ((ff[j]!=i)&(cf[j]!=i)){
    #       I23A[i,j]<-0
    #     }else{
    #       k0 <- kronecker(E[,,ff[j],cf[j]]%*%Sigmainv%*%E[,,i,i],W)
    #       I23A[i,j] <-Tm*sum(k0*Auxit)
    #     }
    #   }
    # }
    # rm(k0)
    I23A <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
    for (i in 1:G){
      Wti <- Matrix::solve(IR-deltas[i]*Wt)
      for (j in 1:(G*(G+1)/2)){
        if ((ff[j]!=i)&(cf[j]!=i)){
          I23A[i,j]<-0
        }else{
          I23A[i,j] <-Tm*sum((Sigmainv[ff[j],cf[j]]*W)*Wti)
        }
      }
    }

    I33A <- Matrix::Matrix(0,nrow=G*(G+1)/2,ncol=G*(G+1)/2)
    for (i in 1:(G*(G+1)/2)){
      for (j in 1:(G*(G+1)/2)){
        I33A[i,j] <- (Tm*N/2)*sum( (Sigmainv%*%E[,,ff[i],cf[i]]) *
                                     Matrix::t((Sigmainv%*%E[,,ff[j],cf[j]])) )
      }
    }

    mislm <- rbind(cbind(I11A,I12A,I13A),
                   cbind(Matrix::t(I12A),I22A,I23A),
                   cbind(Matrix::t(I13A),Matrix::t(I23A),I33A))
    mislm_inv <- try(Matrix::solve(mislm,tol=1e-50))
    if(class(mislm_inv)=="try-error")
      mislm_inv <- MASS::ginv(as.matrix(mislm),tol=1e-70)
    tmp <- sqrt(Matrix::diag(mislm_inv))

    ################################################
    #### MARGINAL LM TEST: Test for SEM in SLM
    ################################################

    if (trace) cat("Computing marginal test... \n")

    # grad rho
    RT <- Matrix::Matrix(matrix(Res,nrow=N*G,ncol=Tm))
    Cgradrest <- Matrix::Matrix(0,nrow=Tm,ncol=G)
    for (i in 1:G){
      SEW <- kronecker(Sigmainv %*% E[,,i,i],W)
      for (t in 1:Tm){
        Cgradrest[t,i] <- Matrix::crossprod(RT[,t],SEW %*% RT[,t])
      }
    }
    gradrho <- Matrix::colSums(Cgradrest)

    ## Information Matrix
    # P1 <- Matrix::Matrix(sum(Matrix::diag(WW)) * IG)
    # P2 <- Matrix::Matrix(sum(Matrix::diag(WtW)) * (Sigmainv * Sigma))
    # Irhorho <- Tm*(P1+P2)
    # Ird <- Matrix::Matrix(0,nrow=G,ncol=G)
    # for (i in 1:G){
    #   for (j in 1:G){
    #     P1s1 <- (Sigma%*%E[,,j,j]) %*% (Sigmainv %*%E[,,i,i])
    #     P1 <- kronecker(kronecker(IT,P1s1),WtW)
    #     P2s1 <- E[,,i,i] %*% E[,,j,j]
    #     P2 <- kronecker(kronecker(IT,P2s1),WW)
    #     Ird[i,j] <- sum( (P1+P2) * Matrix::t(Ainv) )
    #   }
    # }
    P1 <- Matrix::Matrix(sum(Matrix::diag(WW)) * IG)
    P2 <- Matrix::Matrix(sum(Matrix::diag(WtW)) * (Sigmainv * Sigma))
    Irhorho <- Tm*(P1+P2)
    Ird <- Matrix::Matrix(0,nrow=G,ncol=G)
    for (i in 1:G){
      for (j in 1:G){
        P1s1 <- (Sigma%*%E[,,j,j]) %*% (Sigmainv %*%E[,,i,i])
        if (i==j){
          P2s1 <- E[,,i,i] %*% E[,,j,j]
        } else {
          P2s1<-matrix(0,ncol = G,nrow = G)}
        gg <- kronecker(P1s1,WtW)+kronecker(P2s1,WW)
        Ird[i,j] <- Tm*sum(gg*Auxit)
      }
    }

    Irhopslm <- cbind(Matrix::Matrix(0,nrow=G,ncol=fx),
                      Ird, Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2))
    imi <- Matrix::solve((Irhorho - Irhopslm %*% mislm_inv %*%
                            Matrix::t(Irhopslm)))
    LMMqsemenslm <- as.numeric(gradrho %*% imi  %*% gradrho)
    #########################################################
    # Breusch_pagan Test de diagonalidad (Breusch-Pagan 1980)
    # ver: http://www.stata.com/manuals13/rsureg.pdf
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma_corr <- Sigmas$Sigma_corr
    rm(Sigmas)
    index_ltri <- lower.tri(Sigma_corr)
    BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
    # Se ajusta a una Chi con G*(G-1)/2 gl
  res <- list(
    se_betas = as.vector(tmp[1:ncol(X)]),
    se_deltas = as.vector(tmp[(ncol(X)+1):(ncol(X)+G)]),
    cov = as.matrix(mislm_inv),
    LMM = as.numeric(LMMqsemenslm),
    BP = as.numeric(BP)
  )
}

#####################################################
  cov_spsursem <- function(Tm,G,N,Y,X,W,
                           deltas,Sigma,trace=FALSE){
    IT <- Matrix::Diagonal(Tm)
    IR <- Matrix::Diagonal(N)
    IG <- Matrix::Diagonal(G)
    IGR <- Matrix::Diagonal(G*N)
    WW <- W%*%W
    WtW <- Matrix::crossprod(W)
    Wt <- Matrix::t(W)
    E <- get_array_E(G=G)
    ff_cf <- get_ff_cf(G=G)
    ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
    Sigmainv <- try(chol2inv(chol(Sigma)))
    if(class(Sigmainv)=="try-error")
      Sigmainv <- MASS::ginv(as.matrix(Sigma),tol=1e-40)
    delta <- Matrix::Matrix(diag(as.vector(deltas)))
    # Auxiliar matrices
    Aux<-IGR-kronecker(delta,W)
    Auxt<-Matrix::t(Aux)
    Auxi<-Matrix::solve(Aux)
    Auxit<-Matrix::t(Auxi)
    #
    B <- kronecker(IT,Aux)
    Binv <- kronecker(IT,Auxi)
    OME <- kronecker(IT,kronecker(Sigma,IR))
    OMEinv <- kronecker(IT,kronecker(Sigmainv,IR))
    BX <- as(B %*% X,"dgCMatrix")
    BY <- as(B %*% Y,"dgCMatrix")
    Bsem <-Matrix::solve(Matrix::crossprod(BX,OMEinv %*% BX),
                         Matrix::crossprod(BX,OMEinv %*% BY))
    # Bsem <-Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
    #                       Matrix::crossprod(BX,Matrix::solve(OME,B%*%Y)))
    Res <- matrix(BY - (BX %*% Bsem),nrow=nrow(Y))
    ################################
    #### Information Matrix SUR-SEM
    ################################
    J12A <- Matrix::Matrix(0,nrow=ncol(X),ncol=G)
    J13A <- Matrix::Matrix(0,nrow=ncol(X),ncol=G*(G+1)/2)
    J22A <- Matrix::Matrix(0,nrow=G,ncol=G)
    J11A <- Matrix::crossprod(BX,Matrix::solve(OME,BX))
    ## J22A

    for (i in 1:G){
      k0 <- Matrix::solve(IR - deltas[i]*W)%*%W
      for (j in 1:G){
        if (i==j){
          Auxiliar <- Matrix::Matrix(0,ncol=G,nrow=G)
          Auxiliar[j,i] <- Sigmainv[j,i]
          k1<-Matrix::crossprod(Auxi,kronecker(Auxiliar,WtW))
          k2<-kronecker(Matrix::t(Sigma),IR)%*%Auxit
          J22A[i,j] <- Tm*sum(k0*Matrix::t(k0)) + Tm*sum(k1*k2)
        }
        if (j<i){J22A[i,j]=J22A[j,i]}
        if (j>i) {
          Auxiliar <- Matrix::Matrix(0,ncol=G,nrow=G)
          Auxiliar[j,i] <- Sigmainv[j,i]
          k1<-kronecker(Auxiliar,WtW)%*%Auxi
          k2=Auxi%*%kronecker(Sigma,IR)
          J22A[i,j] <- Tm*sum(k1*k2)
        }
      }
    }
rm(k0,k1,k2,Auxiliar)
    ## J23A
    J23A <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
    for (i in 1:G){
      for (j in 1:(G*(G+1)/2)){
        k0 <- kronecker(E[,,i,i] %*% Sigmainv %*% E[,,ff[j],cf[j]],W)
        J23A[i,j] <- Tm*sum(k0*Auxi)
      }
    }
    ## J33A
    J33A <- Matrix::Matrix(0,nrow=G*(G+1)/2,ncol=G*(G+1)/2)
    for (i in 1:(G*(G+1)/2)){
      for (j in 1:(G*(G+1)/2)){
        J33A[i,j] <- (Tm*N/2)*
          sum(diag((Sigmainv %*% E[,,ff[i],cf[i]]) %*%
                     (Sigmainv %*% E[,,ff[j],cf[j]])))
        }
    }

    misem <- rbind(cbind(J11A, J12A, J13A),
                           cbind(Matrix::t(J12A), J22A, J23A),
                           cbind(Matrix::t(J13A), Matrix::t(J23A), J33A))
    misem_inv <- try(Matrix::solve(misem,tol=1e-40))
    if(class(misem_inv)=="try-error")
      misem_inv <- MASS::ginv(as.matrix(misem),tol=1e-40)
    tmp <- sqrt(Matrix::diag(misem_inv))
    ################################################
    #### MARGINAL TEST: Test for SLM in SEM
    ################################################
    if (trace) cat("Computing marginal test... \n")

    RT <- Matrix::Matrix(matrix(Res,nrow=N*G,ncol=Tm))
    Yt <- Matrix::Matrix(matrix(Y,nrow=N*G,ncol=Tm))
    Cgradrest <- Matrix::Matrix(0,nrow=Tm,ncol=G)
    for (i in 1:G){
      for (t in 1:Tm){
        Cgradrest[t,i] <- Matrix::crossprod(RT[,t],
                                            ( kronecker(Sigmainv %*% E[,,i,i],W) -
                                                kronecker(Sigmainv %*% (delta %*% E[,,i,i]), WW))) %*% Yt[,t]
      }
    }
    graddel <- Matrix::colSums(Cgradrest)
    ## matriz de informacion
    PP2 <- Matrix::Matrix(0,nrow=G,ncol=G)
    X_Bsem <- X%*%Bsem
    k3 <- Auxi%*%kronecker(Sigma,IR)
    for (i in 1:G){
      for (j in 1:G){
        if (j<i){PP2[i,j]=PP2[j,i]} # Matriz simetrica
        else {
          k0 <- kronecker(E[,,j,j],Wt)%*%Auxt
          k1 <- kronecker(Sigmainv,IR)%*%Aux%*%kronecker(E[,,i,i],W)
          JJ <- k0%*%k1
          HH <- kronecker(IT,JJ)
          k2 <- Auxi%*%JJ;
          PP2[i,j] <-  Matrix::crossprod(X_Bsem,HH%*%X_Bsem) +Tm*sum(k2*k3)
        }
      }
    }
    rm(X_Bsem,k0,k1,k2,k3,HH)
    ####
    Ideldel <- Tm*sum(Matrix::diag(WW))*IG + PP2
    ####
    Ideltabe <- Matrix::Matrix(0,nrow=G,ncol=ncol(X))
    for (i in 1:G)
    {
      Ideltabe[i,] <- Matrix::crossprod(X%*%Bsem,
                                        kronecker(IT, kronecker(E[,,i,i], Wt))) %*%
        Matrix::crossprod(B,Matrix::solve(OME,B%*%X))
    }
    ## Irhodel
    Irhodel <- Matrix::Matrix(0,nrow=G,ncol=G)
    k0<-kronecker(Sigma,IR)%*%Auxit
    for (i in 1:G){
      for (j in 1:G){
        if(i==j){
          k1<-kronecker(E[,,j,j]%*%Sigmainv,Wt)%*%Aux
          k2<-kronecker(E[,,i,i],W)%*%Auxi
          k12 <- Matrix::t(k1%*%k2)
          k3<-kronecker(E[,,i,i],WW)
          Irhodel[i,j] <-Tm*sum(k0*k12)+Tm*sum(k3*Auxit)
        }
        if (i!=j) {
          k1<-kronecker(E[,,j,j]%*%Sigmainv,Wt)%*%Aux
          k2<-kronecker(E[,,i,i],W)%*%Auxi
          k12 <- Matrix::t(k1%*%k2)
          Irhodel[i,j] <-Tm*sum(k0*k12)
        }
      }
    }
    rm(k0,k1,k2,k12,k3)
  #
    Ideltasig <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
    # OJO: SALE NUMÉRICAMENTE LA MATRIZ NULA... REPASAR: Fdo --> En Matlab tb
    for (i in 1:G){
      for (j in 1:(G*(G+1)/2)){
        k0<-(kronecker(E[,,ff[j],cf[j]]%*%Sigmainv,IR)%*%Aux)%*%kronecker(E[,,i,i],W)
        Ideltasig[i,j] <- Tm*sum(k0*Auxit)
      }
    }
    rm(k0)
    Idelpsem <- cbind(Ideltabe, Matrix::t(Irhodel), Ideltasig)
    imi <- Matrix::solve(Ideldel - Idelpsem %*% Matrix::solve(misem,Matrix::t(Idelpsem)))
    LMMqslmensem <- Matrix::crossprod(graddel,imi %*% graddel)
    ###########################################################
    #### Breusch_pagan Test de diagonalidad (Breusch-Pagan 1980)
    ## ver: http://www.stata.com/manuals13/rsureg.pdf
    Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
    Sigma_corr <- Sigmas$Sigma_corr
    rm(Sigmas)
    index_ltri <- lower.tri(Sigma_corr)
    BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
    # Se ajusta a una Chi con G*(G-1)/2 gl
    res <- list(
      se_betas = as.vector(tmp[1:ncol(X)]),
      se_deltas = as.vector(tmp[(ncol(X)+1):(ncol(X)+G)]),
      cov = as.matrix(misem_inv),
      LMM = as.numeric(LMMqslmensem),
      BP = as.numeric(BP)
    )
  }


#####################################################

cov_spsursarar <- function(Tm,G,N,Y,X,W,
                         deltas,Sigma,trace=FALSE){
  W <- as(W,"dgCMatrix")
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IG <- Matrix::Diagonal(G)
  IGR <- Matrix::Diagonal(G*N)
  WW <- as(W %*% W,"dgCMatrix")
  WtW <- as(Matrix::crossprod(W),"dgCMatrix")
  E <- get_array_E(G=G)
  ff_cf <- get_ff_cf(G=G)
  ff <- ff_cf$ff; cf <- ff_cf$cf; rm(ff_cf)
  Sigmainv <- try(chol2inv(chol(Sigma)))
  if(class(Sigmainv)=="try-error")
    Sigmainv <- MASS::ginv(as.matrix(Sigma),tol=1e-40)
  delta_slm <- deltas[1:G]
  delta_sem <- deltas[(G+1):(2*G)]
  DELTA <- Matrix::Matrix(diag(deltas))
  A <- kronecker(IT,(IGR - kronecker(DELTA[1:G,1:G],W)))
  B <- kronecker(IT,(IGR - kronecker(DELTA[(G+1):(2*G),(G+1):(2*G)],W)))
  OME <- kronecker(IT,kronecker(Sigma,IR))
  BX <- B %*% X
  Bsarar <- Matrix::solve(Matrix::crossprod(BX,Matrix::solve(OME,BX)),
                           Matrix::crossprod(BX,Matrix::solve(OME,B%*%(A%*%Y))))
  Res <- matrix(B %*% (A %*% Y - X %*% Bsarar),nrow=nrow(Y))

  IOME <- kronecker(IT,kronecker(Sigmainv,IR))
  XXTBTIOMEB <- Matrix::crossprod(BX,IOME %*% B)
  ## Auxiliar Matrices
  AAux <- (IGR - kronecker(DELTA[1:G,1:G],W))
  BAux <- (IGR - kronecker(DELTA[(G+1):(2*G),(G+1):(2*G)],W))
  AAuxt <- Matrix::t(AAux)
  BAuxt <- Matrix::t(BAux)
  AAuxi <- Matrix::solve(AAux)
  BAuxi <- Matrix::solve(BAux)
  Wt <- Matrix::t(W)

  ###############################
  ## Information Matrix SUR-SARAR
  ###############################

  ## J11A
  J11A <- XXTBTIOMEB %*% X
  ## J12A
  J12A <- Matrix::Matrix(0,ncol(X),G)
  k0 <- Matrix::solve(A,X %*% Bsarar)
  for (i in 1:G){
    J12A[,i] <- XXTBTIOMEB %*%
                 kronecker(IT,kronecker(E[,,i,i],W)) %*% k0
  }
  rm(XXTBTIOMEB,k0)
  fx <- nrow(J12A)
  cx <- ncol(J12A)
  J13A <- Matrix::Matrix(0,nrow=fx,ncol=G)
  J14A <- Matrix::Matrix(0,nrow=fx,ncol=G*(G+1)/2)
  J21A <- Matrix::t(J12A)
  ## J22A
  J22A <- Matrix::Matrix(0,nrow=G,ncol=G)
  XXBsarar <- X %*% Bsarar
  XXBSARAR <- matrix(XXBsarar,nrow = G*N)
  BF <- Matrix::crossprod(B,kronecker(IT,kronecker(Sigmainv,IR)) %*% B)
  k0 <- BAuxt%*%kronecker(Sigmainv,IR)%*%BAux
  k3 <- Matrix::solve(Matrix::t(BAux),kronecker(Sigma,IR)%*%BAuxi)
  for (i in 1:G){
    for (j in 1:G){
      if (j<i) {J22A[i,j]<-J22A[j,i]}
      else {
      k1 <-kronecker(E[,,j,j],Wt)%*%k0%*%kronecker(E[,,i,i],W)%*%AAuxi
      k4 <- Matrix::solve(Matrix::t(AAux),k1)
      # HH <- kronecker(IT,k4)
      HHXXBsarar <- matrix(k4%*%XXBSARAR,ncol = 1)
      t1 <-Tm*sum(k4*k3)
      F0 <- Matrix::t(XXBsarar)%*%HHXXBsarar+t1
      if (i==j){
        k50<- (IR - delta_slm[i]*W)
        k5 <- Matrix::solve(k50)
        k6 <- k5%*%W
        J22A[i,j] <- Tm*sum(k6*Matrix::t(k6))+ F0
      }
      else {
        J22A[i,j] <- F0
      }
    }
      }
  }
  rm(BF,XXBsarar,XXBSARAR,k0,k1,k3,k4,k5,k6)
  ## J23A
  J23A <- Matrix::Matrix(0,nrow=G,ncol=G)
  k0 <- Matrix::solve(AAux,BAuxi)
  k0t <- Matrix::t(k0)
  kkk <- (kronecker(Sigma,IR) %*% Matrix::t(BAuxi))
  for (i in 1:G){
    k1 <-  kkk%*% kronecker(E[,,i,i]%*% Sigmainv,Wt)
    for (j in 1:G){
      if (j<i) {J23A[i,j]<-J23A[j,i]}
      else {
     k2 <- BAux%*%kronecker(E[,,j,j],W)
     k3 <- kronecker(E[,,i,i]%*%E[,,j,j],WW)
     J23A[i,j] <- Tm*sum((k1%*%k2+k3)*k0t)
      }
      }
  }
  rm(k1,k2,k3,kkk)
  ## J24A
  J24A <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
  for (i in 1:G){
    for (j in 1:(G*(G+1)/2)) {
      if ((ff[j]!=i) & (cf[j]!=i)) {J24A[i,j] <- 0}
      else {
        k1 <- (kronecker(E[,,ff[j],cf[j]] %*% Sigmainv, IR)%*%BAux)%*%kronecker(E[,,i,i],W)
        J24A[i,j] <- Tm*sum(k1*Matrix::t(k0))}
      }
  }
  ## J33A
  J33A <- Matrix::Matrix(0,nrow=G,ncol=G)
  k0 <- Matrix::t(BAuxi%*%kronecker(Sigma,IR))
  k3 <- Matrix::solve(BAux,kronecker(Sigma,IR)%*%Matrix::t(BAuxi))
  for (i in 1:G){
    for (j in 1:G){
      if (i==j){
        BG <- Matrix::solve(IR - delta_sem[i]*W,W)
        k1 <- Matrix::t(BAuxi)%*%kronecker(E[,,j,j]%*%Sigmainv%*%E[,,i,i],WtW)
        J33A[i,j] <- Tm*sum(BG*Matrix::t(BG)) + Tm*sum(k1*k0)
      } else {
        if (i>j) {
          J33A[i,j] <- J33A[j,i]
        } else {
          k1 <- kronecker(E[,,j,j]%*%Sigmainv%*%E[,,i,i],WtW)
          J33A[i,j] <- Tm*sum(k1*k3) # Tm*sum(k1*Matrix::t(k3))
        }
      }
    }
  }
  rm(BG,k0,k1)
  ## J34A
  J34A <- Matrix::Matrix(0,nrow=G,ncol=G*(G+1)/2)
  for (i in 1:G){
    for (j in 1:(G*(G+1)/2)){
      k1 <- kronecker(E[,,i,i] %*% Sigmainv %*% E[,,ff[j],cf[j]],W)
      J34A[i,j] <- Tm*sum(k1*BAuxi)
    }
  }
  rm(k1)
  ## J44A
  J44A <- Matrix::Matrix(0,nrow=G*(G+1)/2,ncol=G*(G+1)/2)
  for (i in 1:(G*(G+1)/2)){
    for (j in 1:(G*(G+1)/2)){
      if (i>j){
        J44A[i,j] <- J44A[j,i]
      } else {
        J44A[i,j] <- sum( (Tm*N/2)*Sigmainv %*%
                     E[,,ff[i],cf[i]] * t(Sigmainv%*%E[,,ff[j],cf[j]]) )
      }
    }
  }
  # Information Matrix
  misarar <- rbind(cbind(J11A, J12A, J13A, J14A),
                           cbind(Matrix::t(J12A), J22A, J23A, J24A),
                           cbind(Matrix::t(J13A), Matrix::t(J23A),
                                         J33A, J34A),
                           cbind(Matrix::t(J14A), Matrix::t(J24A),
                                         Matrix::t(J34A), J44A))
  misarar_inv <- try(Matrix::solve(misarar,tol=1e-40))

  if(class(misarar_inv)=="try-error")
    misarar_inv <- MASS::ginv(as.matrix(misarar),tol=1e-40)
  tmp <- sqrt(Matrix::diag(misarar_inv))
  Sigmas <- get_Sigma(resids=Res,N=N,G=G,Tm=Tm)
  Sigma_corr <- Sigmas$Sigma_corr
  rm(Sigmas)
  index_ltri <- lower.tri(Sigma_corr)
  BP <- N*Tm*sum(Sigma_corr[index_ltri]^2)
  # Se ajusta a una Chi con G*(G-1)/2 gl
  res <- list(
    se_betas = as.vector(tmp[1:ncol(X)]),
    se_deltas = as.vector(tmp[(ncol(X)+1):(ncol(X)+2*G)]),
    cov = as.matrix(misarar_inv),
    LMM = NULL,
    BP = as.numeric(BP)
  )
}
