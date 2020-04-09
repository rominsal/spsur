get_data_spsur <- function(formula, mf, Durbin = FALSE,
                           listw = NULL, zero.policy = NULL,
                           N = NULL, Tm = NULL)
{
  # Function to get data ordered for spatio-temporal SUR
  # Assumption: Data Frame ordered by space (N)
  # Reordenar los datos N*G*Tm
  if (is.null(zero.policy))
    zero.policy <- spatialreg::get.ZeroPolicyOption()
  
  G <- length(attr(formula,"lhs"))
  if (length(attr(formula,"rhs")) < G) {
    for (i in 2:G) {
      attr(formula,"rhs")[[i]] <- attr(formula,"rhs")[[1]]
    }
  }
  if (!is.null(listw)) { 
    N <- length(listw$neighbours)
    W <- Matrix::Matrix(spdep::listw2mat(listw))
    } else {
    if (is.null(Tm)) N <- nrow(mf) else N <- nrow(mf) / Tm
  }  
  if (is.null(W) && is.null(N))
    stop("Dimension of spatial sample is needed")
  if (is.null(Tm)) {
    if (G > 1) { Tm <- nrow(mf) / N }  else { Tm <- 1 }## REPASAR Tm
  }  
  Y <- vector("list",G)
  X <- vector("list",G)
  WX <- vector("list",G)
  p <- rep(0,G)
  dvars <- vector("list",G)
  for (i in 1:G) {
    #nameY <- paste("Y",i,sep = "")
    #nameX <- paste("X",i, sep = "")
    Yi <- as.matrix(Formula::model.part(formula, data = mf, lhs = i))
    Xi <- model.matrix(formula, data = mf, rhs = i)
    names_Yi <- paste(colnames(model.part(formula,
                      data = mf, lhs = i)), i, sep = "_")
    names_colXi <- paste(colnames(model.matrix(formula,
                        data = mf, rhs = i)), i, sep = "_")
    colnames(Yi) <- names_Yi
    colnames(Xi) <- names_colXi
    dvars[[i]] <- c(ncol(Xi), 0L)
    icept <- grep("(Intercept)", colnames(Xi))
    iicept <- length(icept) > 0L
    Y[[i]] <- Yi
    if (Durbin) {
      prefix <- "lag"
      if (iicept) {
        WXi <- spatialreg::create_WX(Xi[,-c(icept)], listw, 
                                         zero.policy = zero.policy, 
                                         prefix = prefix)
      } else {
        WXi <- spatialreg::create_WX(Xi, listw, 
                                         zero.policy = zero.policy, 
                                         prefix = prefix)
      }
      dvars[[i]] <- c(ncol(Xi),ncol(WXi))
      ### CHEQUEAR SIGUIENTES LÍNEAS DE CÓDIGO COMENTADAS....
      #attr(dvars[[i]], "f") <- Durbin
      #attr(dvars[[i]], "inds") <- inds
      #attr(dvars[[i]], "zero_fill") <- zero_fill 
      X[[i]] <- cbind(Xi,WXi)
    } else X[[i]] <- Xi
    p[i] <- ncol(X[[i]])
    #names(Y)[i] <- nameY
    #names(X)[i] <- nameX
  }
  Yt <- vector("list",Tm)
  Xt <- vector("list",Tm)
  for (i in 1:Tm) {
    Yg <- vector("list",G)
    Xg <- vector("list",G)
    for (j in 1:G){
      # Lee R filas de cada vector Yi y matriz Xi
      Yj <- matrix(Y[[j]][((i-1)*N+1):(i*N)], ncol = 1)
      Xj <- X[[j]][((i-1)*N+1):(i*N),]
      colnames(Yj) <- colnames(Y[[j]])
      colnames(Xj) <- colnames(X[[j]])
      #nameYg <- paste("Y",j, sep = "")
      #nameXg <- paste("X",j, sep = "")
      Yg[[j]] <- Yj
      Xg[[j]] <- Xj
      #names(Yg)[j] <- nameYg
      #names(Xg)[j] <- nameXg
    }
    #nameYt <- paste("Y",i, sep="")
    #nameXt <- paste("X",i,sep="")
    ## SEGUIR AQUÍ, CUIDADO CON LOS NOMBRES DE COLUMNAS...
    Yt[[i]] <- unlist(Yg)
    Xt[[i]] <- Matrix::bdiag(Xg)
    #names(Yt)[i] <- nameYt
    #names(Xt)[i] <- nameXt
  }
 
  # Matrices Finales de Datos

  Yf <- Yt[[1]]
  Xf <- Xt[[1]]
  if (Tm > 1) {
  for (i in 2:Tm) {
    Yf <- c(Yf,Yt[[i]])
    Xf <- rbind(Xf,Xt[[i]])
  }
  }
  Yf <- as.matrix(Yf, ncol = 1)
  Xf <- as.matrix(Xf)
  names_colX <- NULL
  for (i in 1:G) {
    names_colX <- c(names_colX,colnames(X[[i]]))
    #names_colXi <- paste(colnames(model.matrix(formula,
    #                 data = mf, rhs = i)), i, sep = "_")
    #names_colX <- c(names_colX,names_colXi)
  }
  colnames(Xf) <- names_colX
  res <- list(Y = Yf, X = Xf, G = G, N = N, Tm = Tm, p = p,
              dvars = dvars)
}
