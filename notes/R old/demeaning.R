demeaning <- function(Y, X, N, G, Tm, p, pdemean = FALSE){
  if(!is.matrix(Y)) Y <-matrix(Y,ncol=1)
  rownames(X) <- rownames(Y) <- as.factor(1:nrow(X))
  Xdem <- X
  Ydem <- Y
  dataXY <- data.frame(Y = Y, X = X)
  if (!pdemean){ # Assumed Tm > 1 and G > 1
    namesX <- paste("X",1:ncol(X),sep="_")
    colnames(dataXY) <- c("Y",namesX)
    dataXY$index_space <- rep(1:N,times=G*Tm)
    dataXY$index_time <- rep(1:Tm,each=N*G)
    dataXY$index_equation <- rep(rep(1:G,each=N),times=Tm)
    for (i in 1:G) {
      for (j in 1:N) {
        dataXY_ij <- dataXY[dataXY$index_equation==i & dataXY$index_space==j,]
        Y_ij <- dataXY_ij$Y
        X_ij <- dataXY_ij[,namesX]
        J_ij <- matrix(rep(colMeans(X_ij),each=nrow(X_ij)),
                       nrow=nrow(X_ij),ncol=ncol(X_ij))
        Ydem[rownames(X_ij),] <- Ydem[rownames(X_ij),] - mean(Y_ij)
        Xdem[rownames(X_ij),] <- Xdem[rownames(X_ij),] - J_ij
      }
    }
    index_interc <- which(colMeans(Xdem) == 0)
    Xdem <- Xdem[,-c(index_interc)]
    pdem <- p
    if (length(index_interc) > 0) pdem <- p - 1 # decrease 1 parameter each equ.
  } else { # Assumed Tm > 1 and G = 1 (pure panel data). Used in spsurtime()
    dataXY$index_space <- rep(1:N, times=Tm)
    dataXY$index_time <- rep(1:Tm, each=N)
    dataXY$index_equation <- rep(1:G, each=N)
    for (j in 1:N){
      dataXY_j <- dataXY[dataXY$index_space==j,]
      Y_j <- dataXY_j[,1]
      X_j <- dataXY_j[,c(2:(ncol(X)+1))]
      J_j <- matrix(rep(colMeans(X_j), each = Tm),
                    nrow = Tm, ncol = ncol(X), byrow = FALSE)
      Ydem[rownames(X_j),] <- Ydem[rownames(X_j),] - mean(Y_j)
      Xdem[rownames(X_j),] <- Xdem[rownames(X_j),] - J_j
    }
    pdem <- p
    if(any(grepl("Interc",colnames(Xdem)))){
      Xdem <- Xdem[,!grepl("Interc",colnames(Xdem))]
      pdem <- pdem - 1
    }
    #Build block-diagonal Matrix for Xdem
    Xdemlist <- list()
    namesXdem <- NULL
    for (i in 1:G){
      namesXdem <- c(namesXdem,paste(colnames(Xdem),"_",i,sep=""))
      Xdemlist[[i]] <- Xdem[dataXY$index_equation==i,]
    }
    Xdem <- as.matrix(Matrix::bdiag(Xdemlist))
    colnames(Xdem) <- namesXdem
  }
  res <- list(Ydem = Ydem, Xdem = Xdem, pdem = pdem)
  res
}
