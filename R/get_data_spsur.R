get_data_spsur <- function(formula=NULL,data=NULL,
                           W=NULL,nR=NULL,nT=NULL)
{
  # Function to get data ordered for spatio-temporal SUR
  # Assumption: Data Frame ordered by space (nR)
  # Reordenar los datos nR*nG*nT
  nG <- length(attr(formula,"lhs"))
  if (length(attr(formula,"rhs")) < nG) {
    for(i in 2:nG){
      attr(formula,"rhs")[[i]] <- attr(formula,"rhs")[[1]]
    }
  }
  if (!is.null(W)) nR <- nrow(W)
  if (is.null(W) && is.null(nR))
    stop("Dimension of spatial sample is needed")
  if (is.null(nT)) nT <- nrow(data) / nR
  Y <- vector("list",nG)
  X <- vector("list",nG)
  p <- rep(0,nG)
  for (i in 1:nG){
    nameY <- paste("Y",i,sep="")
    nameX <- paste("X",i,sep="")
    Yi <- as.matrix(Formula::model.part(formula,data=data,lhs=i))
    Xi <- model.matrix(formula,data=data,rhs=i)
    p[i] <- ncol(Xi)
    Y[[i]] <- Yi
    X[[i]] <- Xi
    names(Y)[i] <- nameY
    names(X)[i] <- nameX
  }
  Yt <- vector("list",nT)
  Xt <- vector("list",nT)
  for (i in 1:nT){
    Yg <- vector("list",nG)
    Xg <- vector("list",nG)
    for (j in 1:nG){
      # Lee R filas de cada vector Yi y matriz Xi
      Yj <- matrix(Y[[j]][((i-1)*nR+1):(i*nR)],ncol=1)
      Xj <- X[[j]][((i-1)*nR+1):(i*nR),]
      nameYg <- paste("Y",j,sep="")
      nameXg <- paste("X",j,sep="")
      Yg[[j]] <- Yj
      Xg[[j]] <- Xj
      names(Yg)[j] <- nameYg
      names(Xg)[j] <- nameXg
    }
    nameYt <- paste("Y",i,sep="")
    nameXt <- paste("X",i,sep="")
    Yt[[i]] <- unlist(Yg)
    Xt[[i]] <- Matrix::bdiag(Xg)
    names(Yt)[i] <- nameYt
    names(Xt)[i] <- nameXt
  }
  # Matrices Finales de Datos

  Yf <- Yt[[1]]
  Xf <- Xt[[1]]
  if (nT>1){
  for(i in 2:nT){
    Yf <- c(Yf,Yt[[i]])
    Xf <- rbind(Xf,Xt[[i]])
  }
  }
  Yf <- as.matrix(Yf,ncol=1)
  Xf <- as.matrix(Xf)
  names_colX <- NULL
  for(i in 1:nG){
    names_colXi <- paste(colnames(model.matrix(formula,
                     data=data,rhs=i)),i,sep="_")
    names_colX <- c(names_colX,names_colXi)
  }
  colnames(Xf) <- names_colX
  res <- list(Y=Yf,X=Xf,nG=nG,nR=nR,nT=nT,p=p)
}
