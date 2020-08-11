sur3_spdiag <- function(Tm, G, N, Y, X, W)
{
  # Estimación del modelo SUR2 sin Efectos Espaciales
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  ols_init <- lm(Y ~ X - 1)
  beta  <- coefficients(ols_init)
  Res <- residuals(ols_init)
  RR <- array(Res, dim=c(N, G, Tm))
  Sigma <- diag(rep(1, G))
  for (i in 1:G){
    for (j in 1:G){
      Sigma[i,j] <- cov(matrix(RR[,i,], ncol = 1),
                          matrix(RR[,j,], ncol = 1))
    }
  }
  LLSUR <- -f_surGRT_beta(beta = beta, Tm = Tm, G = G,
                          N = N, Y = Y, X = X, Sigma = Sigma)
  Sigma_inv <- try(chol2inv(chol(Sigma)))
  # The (if) solve a problem with matriz OME and OMEinv in case of Tm=G=1
  if (Tm==1 & G==1){
  OME <- as(as.matrix(kronecker(IT,kronecker(Sigma,IR))),"dgCMatrix")
  OMEinv <- as(as.matrix(kronecker(IT,kronecker(Sigma_inv,IR))),"dgCMatrix")
  }
  else {
  OME <- as(kronecker(IT,kronecker(Sigma,IR)),"dgCMatrix")
  OMEinv <- as(kronecker(IT,kronecker(Sigma_inv,IR)),"dgCMatrix")
  }
  
  
  beta <- Matrix::solve(Matrix::crossprod(X,OMEinv %*% X),
                        Matrix::crossprod(X,OMEinv %*% Y))
  beta <- as.matrix(beta)
    # Test de diagnostico de la dependencia espacial
  # LM_SUR_SLM2
  #cat("Computing LM-SLM test... \n")
 E <- array(0,dim=c(G,G,G,G))
 for(i in 1:G){
     for(j in 1:G){
         E[i,j,i,j] <- 1
         E[j,i,i,j] <- 1
     }
 }
 # Gradiente
  g_slm <- rep(0,G)
  Res <- Y - X%*%beta # Residuos del SUR sin ee
  #Sigma_inv <- solve(Sigma)
  #W<-as(W,"dgCMatrix")
  for (i in 1:G){
      g_slm[i] <- Matrix::t(Res) %*%(IT %x% 
                                      (Sigma_inv %*% E[,,i,i])
                                       %x% W) %*% Y
  }
  # Matriz de Informacion
  I11 <- Matrix::t(X) %*% (IT %x% Sigma_inv %x% IR) %*% X
  I12 <- matrix(0, nrow = G, ncol = length(beta))
  for (i in 1:G){
      I12[i,] <- matrix(t(X) %*% (IT %x% 
                                   (Sigma_inv%*%E[,,i,i]) 
                                     %x% W) %*% (X %*% beta), 
                                                  nrow=1)
  }
  I22 <- matrix(0, nrow = G, ncol = G)
  OO <- kronecker(IT, kronecker(Sigma, IR)) # Matrix::Matrix(IT %x% Sigma %x% IR)
  tr2 <- sum(W*Matrix::t(W)) # tr <- sum(diag(W%*%W))
  WtW <- Matrix::crossprod(W)
  for (i in 1:G) {
    for (j in 1:G) {
      H <- kronecker(IT, kronecker(E[,,j,j] %*% 
                                     Sigma_inv %*% E[,,i,i],
                                      WtW)) 
      # H <- Matrix::Matrix(IT %x% (E[,,j,j]%*%Sigma_inv%*%E[,,i,i]) %x% WtW)
      tr3 <- sum(H*OO)
      # tr3 <- Tm*sum(diag(E[,,j,j]*Sigma_inv*E[,,i,i]%*%Sigma))*sum(diag(WtW))
      if (i == j) {
        I22[i,j] <- as.numeric(t(X%*%beta) %*% H %*% 
                                 (X%*%beta)) + tr3 + Tm*tr2
      } else {
        I22[i,j] <- as.numeric(t(X%*%beta) %*% H %*% 
                                 (X%*%beta)) + tr3
      }
    }
  }
  LMSURSLM2 <- as.numeric(matrix(g_slm, nrow=1) %*% 
                            solve(I22-I12 %*%
                             solve(I11) %*% t(I12)) %*% 
                              matrix(g_slm, ncol=1))
# LM_SUR_SEM2
    g_sem <- matrix(0,nrow=G)
    for (i in 1:G) {
        g_sem[i] <- as.numeric( t(Res) %*% 
                                  (IT %x% (Sigma_inv %*% 
                                            E[,,i,i]) %x% W) 
                                              %*% Res )
    }

    J22 <- matrix(0, nrow = G, ncol = G)
    tr1 <- sum(Matrix::t(W) * Matrix::t(W))
    tr2 <- sum(Matrix::t(W) * W)
    for (i in 1:G) {
        for (j in 1:G) {
            J22[i,j]=Tm * Sigma_inv[i,j] * Sigma[i,j] * tr1
        }
    }
    J22 <- J22 + Tm * tr2 * Matrix::Diagonal(G)
    LMSURSEM2 <- as.numeric( matrix(g_sem, nrow = 1) %*% 
                               solve(J22) %*% 
                                matrix(g_sem, ncol = 1))
# LM_SUR_SARMA2
    g_sarar1 <- matrix(0, nrow = G)
    g_sarar2 <- matrix(0, nrow = G)
    for (i in 1:G) {
        g_sarar1[i] <- t(Res) %*% (IT %x% 
                                     (Sigma_inv %*% E[,,i,i]) 
                                        %x% W) %*% Y
        g_sarar2[i] <- t(Res) %*% (IT %x% 
                                     (Sigma_inv %*% E[,,i,i]) 
                                        %x% W) %*% Res
    }
    g_sarar <- rbind(g_sarar1, g_sarar2)
# Matriz Informacion
    K11 <- I11
    K12 <- I12
    K22 <- I22
    K33 <- J22
    K23 <- matrix(0, G, G)
    for (i in 1:G) {
        for (j in 1:G) {
            K23[i,j] <- Tm * Sigma_inv[i,j] * Sigma[i,j] *
                         (tr1 + tr2)
        }
    }
    Ird <- matrix(0, G, G)
    WW <- W %*% W
    Sigma <- as(Sigma, "dgCMatrix")
    for (g in 1:G) {
        for (s in 1:G){
            P1s1 <- Sigma %*% E[,,s,s] %*% 
                      Sigma_inv %*% E[,,g,g]
            P1s2 <- IT %x% P1s1
            P1 <- P1s2 %x% WtW
            P2s1 <- E[,,g,g] %*% E[,,s,s]
            P2s2 <- IT %x% P2s1
            P2 <- (P2s2 %x% WW)
            Ird[g,s] <- sum(Matrix::diag(P1 + P2))
        }
    }
    K23 <- Ird
    LMSURSARAR <- as.numeric( matrix(g_sarar, nrow = 1) %*%
               solve(rbind(cbind(K22-K12 %*% solve(K11)
                                 %*% t(K12), K23),
                            cbind(t(K23), K33))) %*% 
                             matrix(g_sarar, ncol = 1) )

# Test LM*-SUR-Lag
#cat("Computing Robust LM*-SUR-SLM test... \n")
# Indexar los elementos diferentes de Sigma
 ff <- rep(0, G*(G+1)/2)
 cf <- rep(0, G*(G+1)/2)
 c1 <- 0; c2<- 0; c3 <-0
    for (k1 in 1:G) {
        c2<- c2+1
        for (i in 1:(G-k1+1)){
            c1 <- c1+1; c3 <- c3+1
            ff[c1] <- c2
            cf[c1] <- c3
        }
        c3 <- c2
    }

    RI44 <- matrix(0, nrow = (G*(G+1)/2), ncol = G*(G+1)/2)
    for (i in 1:(G*(G+1)/2)) {
        for (j in 1:(G*(G+1)/2)) {
            RI44[i,j] <- (Tm*N/2) * sum(diag(Sigma_inv %*%
                                           E[,,ff[i],cf[i]]
                          %*% Sigma_inv %*% E[,,ff[j],cf[j]]))
        }
    }
    ##F Iff <- rbind(cbind(I11,matrix(0,nrow=(k+1)*G,ncol=G*(G+1)/2)),
    ##F         cbind(matrix(0,nrow=G*(G+1)/2,ncol=(k+1)*G),RI44))
    Iff <- as.matrix( rbind(
                      cbind(I11, matrix(0, nrow = ncol(X),
                                           ncol = G*(G+1)/2)),
                      cbind(matrix(0, nrow = G*(G+1)/2,
                                      ncol= ncol(X)), RI44)) )
    
    Ilf <- as.matrix( cbind(I12, matrix(0, nrow = G,
                                        ncol = G*(G+1)/2)) )
    Ilpf <- I22 - Ilf %*% solve(Iff) %*% t(Ilf)
    # Irr <- Tm * tr1 * (Sigma_inv * Sigma + 
    #                      Matrix::Diagonal(G))
    # Irl <- Tm * (tr1+tr2) * (Sigma_inv * Sigma)
    Irr <- K33
    Irl <- K23

LMRSURlag <- as.numeric( Matrix::t(matrix(g_slm, ncol = 1) 
                                   - Irl %*% solve(Irr)
                                     %*% matrix(g_sem, 
                                                ncol = 1))
             %*% solve(Ilpf - Irl %*% solve(Irr) %*% Irl)
             %*% (matrix(g_slm, ncol = 1) 
                  - Irl %*% solve(Irr) %*% matrix(g_sem, 
                                                  ncol = 1)) )

# Test LM*-SUR-Err
#cat("Computing Robust LM*-SUR-SEM test... \n")
LMRSURerr <- as.numeric( Matrix::t(matrix(g_sem, ncol = 1) - 
                                     Irl %*% solve(Ilpf)
                                   %*% matrix(g_slm,ncol=1))
             %*% solve(Irr - Irl %*% solve(Ilpf) %*% Irl)
             %*% (matrix(g_sem, ncol = 1) - 
                    Irl %*% solve(Ilpf) %*%
                      matrix(g_slm, ncol = 1)) )

 statistic <-  c(LMSURSLM2, LMSURSEM2, LMRSURlag, 
                   LMRSURerr, LMSURSARAR)
 names(statistic) <- rep("LM-stat", 5)
 method <- c("LM-SUR-SLM","LM-SUR-SEM",
             "LM*-SUR-SLM", "LM*-SUR-SEM","LM-SUR-SARAR")
 parameter <- c(G, G, G, G, 2*G)
 names(parameter) <- rep("df", 5)
 p.value <- rep(0, 5)
 for (i in 1:length(p.value)) {
   p.value[i] <- pchisq(statistic[i], df = parameter[i], 
                     lower.tail = FALSE)
 }
 lres <- list(length(statistic))
 for (i in 1:length(statistic)) {
   res_i <- list(statistic = statistic[i],
                 parameter = parameter[i],
                 p.value = p.value[i],
                 method = method[i])
   class(res_i) <- "htest"
   lres[[i]] <- res_i
 }
 lres
}
