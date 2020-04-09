function (formula, data = list(), listw, zero.policy = NULL, 
          na.action = na.fail, robust = FALSE, HC = NULL, legacy = FALSE, 
          W2X = TRUE) 
{
  if (!inherits(listw, "listw")) 
    stop("No neighbourhood list")
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
  stopifnot(is.logical(zero.policy))
  if (class(formula) != "formula") 
    formula <- as.formula(formula)
  mt <- terms(formula, data = data)
  mf <- lm(formula, data, na.action = na.action, method = "model.frame")
  na.act <- attr(mf, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
  }
  y <- model.extract(mf, "response")
  if (any(is.na(y))) 
    stop("NAs in dependent variable")
  X <- model.matrix(mt, mf)
  if (any(is.na(X))) 
    stop("NAs in independent variable")
  if (robust) {
    if (is.null(HC)) 
      HC <- "HC0"
    if (!any(HC %in% c("HC0", "HC1"))) 
      stop("HC must be one of HC0, HC1")
  }
  Wy <- lag.listw(listw, y, zero.policy = zero.policy)
  dim(Wy) <- c(nrow(X), 1)
  colnames(Wy) <- c("Rho")
  n <- NROW(X)
  m <- NCOL(X)
  xcolnames <- colnames(X)
  K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
  if (m > 1) {
    WX <- matrix(nrow = n, ncol = (m - (K - 1)))
    if (W2X) 
      WWX <- matrix(nrow = n, ncol = ncol(WX))
    for (k in K:m) {
      wx <- lag.listw(listw, X[, k], zero.policy = zero.policy)
      if (W2X) 
        wwx <- lag.listw(listw, wx, zero.policy = zero.policy)
      if (any(is.na(wx))) 
        stop("NAs in lagged independent variable")
      WX[, (k - (K - 1))] <- wx
      if (W2X) 
        WWX[, (k - (K - 1))] <- wwx
    }
    if (W2X) 
      inst <- cbind(WX, WWX)
    else inst <- WX
  }
  if (K == 2 && listw$style != "W") {
    wx1 <- as.double(rep(1, n))
    wx <- lag.listw(listw, wx1, zero.policy = zero.policy)
    if (W2X) 
      wwx <- lag.listw(listw, wx, zero.policy = zero.policy)
    if (m > 1) {
      inst <- cbind(wx, inst)
      if (W2X) 
        inst <- cbind(wwx, inst)
    }
    else {
      inst <- matrix(wx, nrow = n, ncol = 1)
      if (W2X) 
        inst <- cbind(inst, wwx)
    }
  }
  result <- tsls(y = y, yend = Wy, X = X, Zinst = inst, robust = robust, 
                 HC = HC, legacy = legacy)
  result$zero.policy <- zero.policy
  result$robust <- robust
  if (robust) 
    result$HC <- HC
  result$legacy <- legacy
  result$listw_style <- listw$style
  result$call <- match.call()
  class(result) <- "stsls"
  result
}
