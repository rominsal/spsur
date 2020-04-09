# File containing auxiliary functions
get_array_E <- function(G){
  E <- array(0,dim=c(G,G,G,G))
  for(i in 1:G){
    for(j in 1:G){
      E[i,j,i,j] <- 1
      E[j,i,i,j] <- 1
    }
  }
  E
}

#######################################################
get_Sigma <- function(resids, N, G, Tm) {
  RR <- array(resids, dim = c(N, G, Tm))
  Sigma <- diag(rep(1,G))
  for (i in 1:G){
    for (j in 1:G){
      Sigma[i,j] <- cov(matrix(RR[,i,],ncol=1),
                        matrix(RR[,j,],ncol=1))
    }
  }
  Sigma_inv <- try(chol2inv(chol(Sigma)))
  if (inherits(Sigma_inv, "try-error"))
    Sigma_inv <- MASS::ginv(as.matrix(Sigma), tol = 1e-40)

  Sigma_corr <- diag(rep(1, G))
  for (i in 1:G) {
    for (j in 1:G) {
      Sigma_corr[i,j] <- cor(matrix(RR[,i,], ncol = 1),
                             matrix(RR[,j,], ncol = 1))
    }
  }
  Sigma_corr_inv <- try(solve(Sigma_corr))
  if (inherits(Sigma_corr_inv, "try-error"))
    Sigma_corr_inv <- MASS::ginv(as.matrix(Sigma_corr),
                                 tol = 1e-40)
  res <- list(Sigma = as.matrix(Sigma),
              Sigma_inv = as.matrix(Sigma_inv),
              Sigma_corr = as.matrix(Sigma_corr),
              Sigma_corr_inv = as.matrix(Sigma_corr_inv))
  res
}

#######################################################
get_ff_cf <- function(G){
  # Indexar los elementos diferentes de Sigma
  ff <- rep(0,G*(G+1)/2)
  cf <- rep(0,G*(G+1)/2)
  c1 <- 0;c2 <- 0;c3 <- 0;
  for (k1 in 1:G){
    c2 <- c2+1
    for (i in 1:(G-k1+1)){
      c1 <- c1+1
      c3 <- c3+1;
      ff[c1] <- c2
      cf[c1] <- c3
    }
    c3 <- c2
  }
  res <- list(ff=ff,cf=cf)
  res
}
