get_data_spsur <- function(formula=NULL,data=NULL,W=NULL,N=NULL,Tm=NULL)
{
  # Function to get data ordered for spatio-temporal SUR
  # Assumption: Data Frame ordered by space (N)
  # Reordenar los datos N*G*Tm
  G <- length(attr(formula,"lhs"))
  if (length(attr(formula,"rhs")) < G) {
    for(i in 2:G){
      attr(formula,"rhs")[[i]] <- attr(formula,"rhs")[[1]]
    }
  }
  if (!is.null(W)) N <- nrow(W)
  if (is.null(W) && is.null(N))
    stop("Dimension of spatial sample is needed")
  if (is.null(Tm)) Tm <- nrow(data) / N
  Y <- vector("list",G)
  X <- vector("list",G)
  p <- rep(0,G)
  for (i in 1:G){
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
  Yt <- vector("list",Tm)
  Xt <- vector("list",Tm)
  for (i in 1:Tm){
    Yg <- vector("list",G)
    Xg <- vector("list",G)
    for (j in 1:G){
      # Lee R filas de cada vector Yi y matriz Xi
      Yj <- matrix(Y[[j]][((i-1)*N+1):(i*N)],ncol=1)
      Xj <- X[[j]][((i-1)*N+1):(i*N),]
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
  if (Tm>1){
  for(i in 2:Tm){
    Yf <- c(Yf,Yt[[i]])
    Xf <- rbind(Xf,Xt[[i]])
  }
  }
  Yf <- as.matrix(Yf,ncol=1)
  Xf <- as.matrix(Xf)
  names_colX <- NULL
  for(i in 1:G){
    names_colXi <- paste(colnames(model.matrix(formula,
                     data=data,rhs=i)),i,sep="_")
    names_colX <- c(names_colX,names_colXi)
  }
  colnames(Xf) <- names_colX
  res <- list(Y=Yf,X=Xf,G=G,N=N,Tm=Tm,p=p)
}
