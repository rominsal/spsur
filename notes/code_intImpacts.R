function (rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames, 
          interval, type, tr, R, listw, evalues, tol, empirical, Q, 
          icept, iicept, p, mess = FALSE, samples = NULL, zero_fill = NULL, 
          dvars = NULL) 
{
  if (is.null(evalues)) {
    if (is.null(listw) && is.null(tr)) 
      stop("either tr or listw must be given")
  }
  else {
    if (!is.null(listw)) {
      warning("evalues given: listw will be ignored")
      listw <- NULL
    }
    if (!is.null(tr)) {
      warning("evalues given: listw will be ignored")
      tr <- NULL
    }
  }
  timings <- list()
  .ptime_start <- proc.time()
  if (is.null(listw)) {
    q <- length(tr) - 1L
    g <- rho^(0:q)
    T <- matrix(c(1, tr[-(q + 1)]/n), nrow = 1)
    if (type == "mixed" || type == "sacmixed") {
      T <- rbind(T, tr/n)
    }
    if (is.null(evalues)) {
      res <- lagImpacts(T, g, P)
      cmethod <- "trace"
    }
    else {
      if (type == "mixed" || type == "sacmixed") 
        stop("eigenvalue mixed impacts not available")
      if (length(evalues) != n) 
        stop("wrong eigenvalue vector length")
      res <- lagImpacts_e(rho, P, n, evalues)
      cmethod <- "evalues"
    }
    if (!is.null(Q)) {
      if (!is.numeric(Q) || length(Q) > 1L) 
        stop("Invalid Q argument")
      if (Q > length(tr)) 
        stop("Q larger than length of tr")
      Qres <- lagDistrImpacts(T, g, P, q = as.integer(Q))
      attr(res, "Qres") <- Qres
    }
    timings[["trace_impacts"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
    if (!is.null(R)) {
      if (is.null(samples)) {
        samples <- MASS::mvrnorm(n = R, mu = mu, Sigma = Sigma, 
                           tol = tol, empirical = empirical)
        if (mess) 
          samples[, irho] <- 1 - exp(samples[, irho])
      }
      if (!is.null(interval)) {
        check <- ((samples[, irho] > interval[1]) & 
                    (samples[, irho] < interval[2]))
        if (any(!check)) 
          samples <- samples[check, ]
      }
      timings[["impacts_samples"]] <- proc.time() - .ptime_start
      .ptime_start <- proc.time()
      sres <- apply(samples, 1, processSample, irho = irho, 
                    drop2beta = drop2beta, type = type, iicept = iicept, 
                    icept = icept, zero_fill = zero_fill, dvars = dvars, 
                    T = T, Q = Q, q = q, evalues = evalues)
      timings[["process_samples"]] <- proc.time() - .ptime_start
      .ptime_start <- proc.time()
      if (length(bnames) == 1L) {
        direct <- as.mcmc(t(matrix(sapply(sres, function(x) x$direct), 
                                   nrow = 1)))
        indirect <- as.mcmc(t(matrix(sapply(sres, function(x) x$indirect), 
                                     nrow = 1)))
        total <- as.mcmc(t(matrix(sapply(sres, function(x) x$total), 
                                  nrow = 1)))
      }
      else {
        direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
        indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
        total <- as.mcmc(t(sapply(sres, function(x) x$total)))
      }
      colnames(direct) <- bnames
      colnames(indirect) <- bnames
      colnames(total) <- bnames
      ssres <- list(direct = direct, indirect = indirect, 
                    total = total)
      if (!is.null(Q)) {
        Qdirect <- as.mcmc(t(sapply(sres, function(x) attr(x, 
                                                           "Qres")$direct)))
        Qindirect <- as.mcmc(t(sapply(sres, function(x) attr(x, 
                                                             "Qres")$indirect)))
        Qtotal <- as.mcmc(t(sapply(sres, function(x) attr(x, 
                                                          "Qres")$total)))
        Qnames <- c(sapply(bnames, function(x) paste(x, 
                                                     1:Q, sep = "__Q")))
        if (length(Qnames) == 1L) {
          Qdirect <- t(Qdirect)
          Qindirect <- t(Qindirect)
          Qtotal <- t(Qtotal)
        }
        colnames(Qdirect) <- Qnames
        colnames(Qindirect) <- Qnames
        colnames(Qtotal) <- Qnames
        Qmcmc <- list(direct = Qdirect, indirect = Qindirect, 
                      total = Qtotal)
        attr(ssres, "Qmcmc") <- Qmcmc
      }
      timings[["postprocess_samples"]] <- proc.time() - 
        .ptime_start
      res <- list(res = res, sres = ssres)
    }
    attr(res, "method") <- cmethod
  }
  else {
    stopifnot(length(listw$neighbours) == n)
    V <- listw2mat(listw)
    e <- eigen(V, only.values = TRUE)$values
    if (is.complex(e)) 
      interval <- 1/(range(Re(e)))
    else interval <- 1/(range(e))
    SW <- invIrW(listw, rho)
    if (type == "lag" || type == "sac") 
      res <- lagImpactsExact(SW, P, n)
    else if (type == "mixed" || type == "sacmixed") 
      res <- mixedImpactsExact(SW, P, n, listw)
    timings[["weights_impacts"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
    if (!is.null(R)) {
      samples <- mvrnorm(n = R, mu = mu, Sigma = Sigma, 
                         tol = tol, empirical = empirical)
      check <- ((samples[, irho] > interval[1]) & (samples[, 
                                                           irho] < interval[2]))
      if (any(!check)) 
        samples <- samples[check, ]
      timings[["impacts_samples"]] <- proc.time() - .ptime_start
      .ptime_start <- proc.time()
      sres <- apply(samples, 1, processXSample, drop2beta = drop2beta, 
                    type = type, iicept = iicept, icept = icept, 
                    n = n, listw = listw, irho = irho, zero_fill = zero_fill, 
                    dvars = dvars)
      timings[["process_samples"]] <- proc.time() - .ptime_start
      .ptime_start <- proc.time()
      if (length(bnames) == 1L) {
        direct <- as.mcmc(t(matrix(sapply(sres, function(x) x$direct), 
                                   nrow = 1)))
        indirect <- as.mcmc(t(matrix(sapply(sres, function(x) x$indirect), 
                                     nrow = 1)))
        total <- as.mcmc(t(matrix(sapply(sres, function(x) x$total), 
                                  nrow = 1)))
      }
      else {
        direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
        indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
        total <- as.mcmc(t(sapply(sres, function(x) x$total)))
      }
      colnames(direct) <- bnames
      colnames(indirect) <- bnames
      colnames(total) <- bnames
      timings[["postprocess_samples"]] <- proc.time() - 
        .ptime_start
      res <- list(res = res, sres = list(direct = direct, 
                                         indirect = indirect, total = total))
    }
    attr(res, "method") <- "exact"
  }
  if (!is.null(R)) 
    attr(res, "samples") <- list(samples = samples, irho = irho, 
                                 drop2beta = drop2beta)
  attr(res, "type") <- type
  attr(res, "bnames") <- bnames
  attr(res, "haveQ") <- !is.null(Q)
  attr(res, "timings") <- do.call("rbind", timings)[, c(1, 
                                                        3)]
  class(res) <- "lagImpact"
  res
}
