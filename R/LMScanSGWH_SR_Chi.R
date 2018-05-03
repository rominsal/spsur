LMScanSGWH_SAR_Chi <- function(y=y,X=XX,W=W,rho=rho,mC=mC,NNN=NNN){
  # ---------------------------------------------------------------------
  #  Inputs
  #    y    : dependent variable
  #    XX   : independent variables (with ones in the first column)
  #    W    : Matrix W
  #    mC   : minimum cluster size
  #  Outputs
  #    p_value : Monte Carlo simulated (adjusted) p-value
  # ----------------------------------------------------------------------------
  #   
  # Fernando A. Lopez Hernandez, 19/03/2018
  # Dpto. Metodos Cuantitativos e Informaticos
  # Universidad Politecnica de Cartagena
  # E-mail: fernando.lopez@upct.es
  #
  # Ws <- normw(W)
  
  R <- length(y)
  Ws<-W/matrix(rowSums(W),nrow=R,ncol=R)
  A <- diag(R)-rho*Ws
  Beta <- Matrix::solve(t(XX)%*%XX,t(XX)%*%A%*%y)
  u <- (A%*%y-XX%*%Beta)
  sigma2 <- as.numeric(Matrix::crossprod(u)/R)
  k <- dim(XX)[2]
  # Gradiente
  IA <- Matrix::solve(A)
  tIAWs <- sum(IA*t(Ws))
  I11 <- crossprod(XX)
  I12 <- crossprod(XX,Ws%*%IA%*%XX%*%Beta)
  I13 <- matrix(0,k,1)
  WsIA <- Ws%*%IA
  XXBeta <- XX%*%Beta
  I22 <- sigma2*sum(WsIA*(t(WsIA)+WsIA))+crossprod(XXBeta,t(WsIA)%*%WsIA%*%XXBeta)
  I23 <- tIAWs
  I33 <- R/(2*sigma2)
  II11 <- rbind(cbind(I11, I12, I13),cbind(t(I12),I22,I23),cbind(t(I13),I23,I33))
  a44 <- det(as.numeric(1/sigma2)*II11)
  WsIA <- diag(WsIA)
  dA33 <- det(II11)
  dA22 <- det(rbind(cbind(I11,I12),cbind(t(I12),I22)))
  dXX <- det(crossprod(XX))
  Rr <-round(R/2, digits  <-  0)
  LM <- matrix(0,R,Rr-mC+1)
  F<-matrix(WsIA[NNN[,1:Rr]],ncol = Rr)
  ff <- sigma2*t(apply(F,1,cumsum))
  uu <- matrix(u[NNN[,1:Rr]]^2,ncol=Rr)
  cuu <- t(apply(uu,1,cumsum))/as.numeric(2*sigma2)
  mRZ <-t( matrix(1:Rr,nrow=dim(cuu)[2],ncol=dim(cuu)[1]))
  cuu2 <- -mRZ/2+cuu
  ksigma2 <- (1/sigma2)^(k+3)
  DI <- as.numeric(ksigma2)*(as.numeric(sigma2)*mRZ*dA33/2-(mRZ/2)*(mRZ*dA22/2-I23*dXX*ff)+(ff*dXX)*(I23*mRZ/2-ff*R/(2*sigma2)))
  MLM <- ((cuu2*cuu2)*a44)/DI
  LM <- MLM[,mC:dim(MLM)[2]]
  # results.LMp <- LM
  # results.mLM <- max(LM)
  results <- list(LM =LM, testLM=max(LM))
}
