get_data_spsur <- function(formula, mf, Durbin = FALSE,
                           listw = NULL, zero.policy = NULL,
                           N = NULL, Tm = NULL) {
  # Function to get data ordered for spatio-temporal SUR
  # Assumption: Data Frame ordered by space (N)
  # Order of data N*G*Tm
  if (is.null(zero.policy))
    zero.policy <- spatialreg::get.ZeroPolicyOption()
  G <- length(attr(formula, "lhs"))
  if (length(attr(formula, "rhs")) < G) { ## Repeat RHS...
    for (i in 2:G) {
      attr(formula, "rhs")[[i]] <- attr(formula, "rhs")[[1]]
    }
  }
  if (inherits(Durbin, "formula")) {
    if (!inherits(Durbin, "Formula")) 
      Durbin <- Formula::Formula(Durbin) 
    if (length(attr(Durbin, "rhs")) < G) { ## Repeat RHS...
      for (i in 2:G) {
        attr(Durbin, "rhs")[[i]] <- attr(Durbin, "rhs")[[1]]
      }
    }
  }
  if (is.null(Tm)) {
    N <- nrow(mf) 
  } else {
    N <- nrow(mf) / Tm
  }  
  if (!is.null(listw)) { 
    W <- Matrix::Matrix(spdep::listw2mat(listw))
  } 
  if (is.null(listw) && is.null(N))
    stop("Dimension of spatial sample is needed")
  if (is.null(Tm)) {
    if (G > 1) { Tm <- nrow(mf) / N }  else { Tm <- 1 }## REPASAR Tm
  }  
  Y <- vector("list", G)
  X <- vector("list", G)
  WX <- vector("list", G)
  p <- rep(0, G)
  dvars <- vector("list", G)
  for (i in 1:G) {
    Yi <- as.matrix(Formula::model.part(formula, data = mf, 
                                        lhs = i))
    Xi <- model.matrix(formula, data = mf, rhs = i)
    colnames(Xi) <- paste(colnames(Xi), i, sep = "_")
    colnames(Yi) <- paste(colnames(Yi), i, sep = "_")
    dvars[[i]] <- c(ncol(Xi), 0L)
    icept <- grep("(Intercept)", colnames(Xi))
    iicept <- length(icept) > 0L
    Y[[i]] <- Yi
    if (isTRUE(Durbin) || inherits(Durbin, "formula")) {
      prefix <- "lag"
      if (isTRUE(Durbin)) {
        if (iicept) {
          WXi <- spatialreg::create_WX(Xi[,-c(icept), drop = FALSE], 
                                       listw, 
                                       zero.policy = zero.policy, 
                                       prefix = prefix)
        } else {
          WXi <- spatialreg::create_WX(Xi, listw, 
                                       zero.policy = zero.policy, 
                                       prefix = prefix)
        }
      } else { ##  Durbin is formula ...
        fXi <- try(model.matrix(Durbin, data = mf, rhs = i),
                   silent = TRUE)
        if (inherits(fXi, "try-error"))
          stop("Durbin variable mist-match")
        if (iicept) {
          WXi <- spatialreg::create_WX(fXi[,-c(icept), drop = FALSE], 
                                       listw, 
                                       zero.policy = zero.policy, 
                                       prefix = prefix)
        } else {
          WXi <- spatialreg::create_WX(fXi, listw, 
                                       zero.policy = zero.policy, 
                                       prefix = prefix)
        }
        #WXi <- spatialreg::create_WX(fXi, listw, 
        #                             zero.policy = zero.policy, 
        #                             prefix = prefix)
        colnames(WXi) <- paste(colnames(WXi), i, sep = "_")
        inds <- match(substring(colnames(WXi), 5, 
                                nchar(colnames(WXi))), 
                      colnames(Xi))
        if (anyNA(inds)) 
          stop("WX variables not in X: ", 
               paste(substring(colnames(WXi), 5, 
                          nchar(colnames(WXi)))[is.na(inds)], 
                     collapse = " "))
        if (iicept) {
          xni <- colnames(Xi)[-1]
        } else {
          xni <- colnames(Xi)
        }
        wxni <- substring(colnames(WXi), nchar(prefix) + 2, 
                          nchar(colnames(WXi)))
        zero_fill <- NULL
        if (length((which(!(xni %in% wxni)))) > 0L) 
          zero_fill <- length(xni) + (which(!(xni %in% wxni)))
      }
      dvars[[i]] <- c(ncol(Xi), ncol(WXi))
      if (inherits(Durbin, "formula")) {
        attr(dvars[[i]], "f") <- attr(Durbin, "rhs")[[i]]
        attr(dvars[[i]], "inds") <- inds
        attr(dvars[[i]], "zero_fill") <- zero_fill 
      }
      X[[i]] <- cbind(Xi, WXi)
    } else X[[i]] <- Xi
    p[i] <- ncol(X[[i]])
  }
  Yt <- vector("list", Tm)
  Xt <- vector("list", Tm)
  for (i in 1:Tm) {
    Yg <- vector("list", G)
    Xg <- vector("list", G)
    for (j in 1:G){
      # Lee R filas de cada vector Yi y matriz Xi
      Yj <- matrix(Y[[j]][((i-1)*N+1):(i*N)], ncol = 1)
      Xj <- X[[j]][((i-1)*N+1):(i*N),]
      colnames(Yj) <- colnames(Y[[j]])
      colnames(Xj) <- colnames(X[[j]])
      Yg[[j]] <- Yj
      Xg[[j]] <- Xj
    }
    Yt[[i]] <- unlist(Yg)
    Xt[[i]] <- Matrix::bdiag(Xg)
  }
  # Final Matrices
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
  }
  colnames(Xf) <- names_colX
  res <- list(Y = Yf, X = Xf, G = G, N = N, Tm = Tm, p = p,
              dvars = dvars)
}
