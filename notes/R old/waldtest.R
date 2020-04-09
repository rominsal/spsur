# Function to compute Wald tests about
# equals beta and delta in spatial SUR models
# with the same regressors in each equation.
# Falta incluir el caso sarar

wald_test <- function(betas,deltas,cov,p,G){
    H1 <- NULL
    H2 <- NULL
    for (i in 2:p){
      e1 <- matrix(0,nrow=1,ncol=p)
      e2 <- matrix(0,nrow=1,ncol=p)
      e3 <- matrix(0,nrow=1,ncol=p)
      e1[1,i] <- 1
      e2[1,i] <- -1
      block1 <- matrix(rep(e1,G-1),ncol=p,byrow=TRUE)
      block2 <- matrix(c(e2,rep(e3,G-2)),ncol=p,byrow=TRUE)
      H1 <- rbind(H1,block1)
      H2 <- rbind(H2,block2)
      H_betas <- cbind(H1,H2)
    }
    for(j in 1:(G-2)) {
      H2 <- rbind(0,H2[1:(nrow(H2)-1),])
      H_betas <- cbind(H_betas,H2)
    }
    df_betas <- nrow(H_betas)
    h_betas <- matrix(0,nrow=df_betas,ncol=1)
    holg_betas <- H_betas%*%betas - h_betas
    cov_betas <- cov[1:length(betas),1:length(betas)]
    wald_betas <- (t(holg_betas) %*%
              solve(H_betas %*% cov_betas %*% t(H_betas),
                  holg_betas)) / df_betas
    pval_betas <- pchisq(wald_betas,
                  df=df_betas,lower.tail=FALSE)

# Wald Test for equal deltas in G equations
    # Falta incluir el Sarar

    cov_deltas <- cov[(length(betas)+1):
                    (length(betas)+length(deltas)),
                    (length(betas)+1):
                      (length(betas)+length(deltas))]
    H_deltas <- NULL
    for (i in 1:(G-1)){
      e <- rep(0,G-1)
      e[i] <- -1
      H_deltas <- rbind(H_deltas,c(1,e))
    }
    # Cambia matriz para modelo SARAR
    if(length(deltas)==2*G){
      H_deltas <- as.matrix(Matrix::bdiag(H_deltas,H_deltas))
    }
    df_deltas <- nrow(H_deltas)
    h_deltas <- matrix(0,nrow=df_deltas,ncol=1)
    holg_deltas  <- H_deltas %*%deltas  - h_deltas
    wald_deltas  <- (t(holg_deltas ) %*%
                       solve(H_deltas  %*% cov_deltas  %*% t(H_deltas ),
                             holg_deltas )) / df_deltas
    pval_deltas  <- pchisq(wald_deltas,
                           df=df_deltas,lower.tail=FALSE)
    # Wald Test for equal betas and deltas in G equations
    cov_all <- cov[1:(length(betas)+length(deltas)),
                   1:(length(betas)+length(deltas))]
    H_all <- as.matrix(Matrix::bdiag(H_betas,H_deltas))
    h_all <- rbind(h_betas,h_deltas)
    par_all <- c(betas,deltas)
    df_all <- nrow(H_all)
    holg_all  <- H_all %*%par_all  - h_all
    wald_all  <- (t(holg_all ) %*%
            solve(H_all  %*% cov_all  %*% t(H_all),
                             holg_all )) / df_all
    pval_all  <- pchisq(wald_all,
                           df=df_all,lower.tail=FALSE)
  cat("H_0: Equal ",p-1," slopes in ",G," equations: \n")
  cat("Wald statistic (chi-squared distr.) ",
      round(wald_betas,3),"\n")
  cat("Degrees of freedom: ",df_betas,
      " p-value: (",round(pval_betas,3),") \n\n")
  cat("H_0: Equal spatial parameters in ",G," equations: \n")
  cat("Wald statistic (chi-squared distr.) ",
      round(wald_deltas,3),"\n")
  cat("Degrees of freedom: ",df_deltas,
      " p-value: (",round(pval_deltas,3),") \n\n")
  cat("H_0: Equal slopes and spatial parameters in ",
      G," equations: \n")
  cat("Wald statistic (chi-squared distr.) ",
      round(wald_all,3),"\n")
  cat("Degrees of freedom: ",df_all,
      " p-value: (",round(pval_all,3),") \n\n")
  res <- list(wald_betas = wald_betas,
              df_betas = df_betas,
              pval_betas = pval_betas,
              wald_deltas = wald_deltas,
              df_deltas = df_deltas,
              pval_deltas = pval_deltas
              )
}
