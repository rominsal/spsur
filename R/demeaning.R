demeaning <- function(Y,X,nR,nG,nT,p){
  if(!is.matrix(Y)) Y <-matrix(Y,ncol=1)
  rownames(X) <- rownames(Y) <- as.factor(1:nrow(X))
  Xdem <- X
  Ydem <- Y
  dataXY <- data.frame(Y=Y,X=X)
  namesX <- paste("X",1:ncol(X),sep="_")
  colnames(dataXY) <- c("Y",namesX)
  dataXY$index_space <- rep(1:nR,times=nG*nT)   
  dataXY$index_time <- rep(1:nT,each=nR*nG)
  dataXY$index_equation <- rep(rep(1:nG,each=nR),times=nT)
  for (i in 1:nG) {
    for (j in 1:nR) {
      dataXY_ij <- dataXY[dataXY$index_equation==i & dataXY$index_space==j,]
      Y_ij <- dataXY_ij$Y
      X_ij <- dataXY_ij[,namesX]
      J_ij <- matrix(rep(colMeans(X_ij),each=nrow(X_ij)),
                     nrow=nrow(X_ij),ncol=ncol(X_ij))
      Ydem[rownames(X_ij),] <- Ydem[rownames(X_ij),] - mean(Y_ij)
      Xdem[rownames(X_ij),] <- Xdem[rownames(X_ij),] - J_ij 
    }
  }
  index_interc <- which(colMeans(Xdem)==0)
  Xdem <- Xdem[,-c(index_interc)]
  pdem <- p
  if (length(index_interc) > 0) pdem <- p - 1 # decrease 1 parameter each equ.
  res <- list(Ydem=Ydem, Xdem=Xdem, pdem=pdem)
  res
}
