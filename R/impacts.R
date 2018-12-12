#' @name impacts
#' @rdname impacts
#'
#' @title Impacts estimation of spatial SUR model
#'
#' @description Returns the marginal effects of a Spatial SUR model
#'   (SUR-SAR; SUR-SDM; SUR-SARAR)
#'
#' @param spsurfit A fitted object of class \emph{spsur}.
#' @param nsim Number of simulations. Default = 1000.
#'
#' @details The marginal effects (impacts) of spatial autoregressive models
#'   (SAR; SDM; SARAR) are more complicated than usual measurements of impacts
#'   for non spatial models. Here we follow LeSage and Pace (2009) and propose
#'   the following summaries for impact measures for Spatial SUR models:
#'   \itemize{
#'     \item \strong{Average direct effects}: The average over all the
#'       observations of the effects of the change of an explanatory variable
#'       of a single observation on the choice probability of that same
#'       observation.
#'     \item \strong{Average indirect effects}: The average over all the
#'       observations of the effect of a change on a explanatory variable on
#'       the choice probability of the neighbouring observations.
#'     \item \strong{Average total effects}: The sum of direct and indirect
#'       impacts.
#'   }
#'
#' @return It returns the marginal effects of the estimated spatial SUR model
#'   and simulted values.
#'   \tabular{ll}{
#'   \code{table_dir_eff} \tab table of average direct effects. \cr
#'   \code{table_ind_eff} \tab table of average indirect effects. \cr
#'   \code{table_tot_eff} \tab table of average total effects. \cr
#'   \code{sim_dir_eff} \tab simulate direct effects. \cr
#'   \code{sim_ind_eff} \tab simulate indirect effects. \cr
#'   \code{sim_tot_eff} \tab simulate total effects. \cr
#'   }
#'
#' @references
#'   \itemize{
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'      \item LeSage, J., and Pace, R. K. (2009). \emph{Introduction to spatial
#'        econometrics}. Chapman and Hall/CRC.
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @seealso
#'
#' \code{\link{spsurml}},  \code{\link{spsur3sls}}
#'
#' @examples
#'
#' ####################################
#' ######## CROSS SECTION DATA ########
#' ####################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#'
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## A SUR-SAR model
#' spcSUR.sar <-spsurml(Form = Tformula, data = spc, type = "sar", W = Wspc)
#' summary(spcSUR.sar)
#' eff.spcSUR.sar <- impacts(spcSUR.sar)
#'
#' ## A SUR-SDM model
#' spcSUR.sdm <-spsurml(Form = Tformula, data = spc, type = "sdm", W = Wspc)
#' summary(spcSUR.sdm)
#' eff.spcSUR.sdm <- impacts(spcSUR.sdm, nsim = 300)
#'
#' ## A SUR-SARAR model
#' spcSUR.sarar <-spsurml(Form = Tformula, data = spc, type = "sarar", W = Wspc)
#' summary(spcSUR.sarar)
#' eff.spcSUR.sarar <- impacts(spcSUR.sarar, nsim = 300)
#'
#' ####################################
#' ######## PANEL DATA (nG>1; nT>1) ###
#' ####################################
#'
#' data(NAT)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SAR model
#' NATSUR.sar <-spsurml(Form = Tformula, data = NAT, type = "sar", W = W)
#' summary(NATSUR.sar)
#' eff.NAT.sar <- impacts(NATSUR.sar, nsim = 100)
#'
#' ####################################
#' ######## PANEL DATA (nG=1; nT>1) ###
#' ####################################
#'
#' data("unemp_it")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' unempitml_sar <- spsurtime(Form = form_un, data = unemp_it,
#'                            time = unemp_it$year, W = W_italy, type = "sar",
#'                            method = "ml")
#' summary(unempitml_sar)
#' eff.unempitml_sar <- impacts(unempitml_sar, nsim = 100)
#' @export
impacts <- function(spsurfit, nsim = 1000){
  z <- spsurfit
  type <- z$type
  if (!(type == "sar" | type == "sdm" | type == "sarar")){
    stop("type must be sar, sdm or sarar")
  }
  if (type == "sarar"){ # Treat sarar as a sar case...
    type <- "sar"
    z$deltas <- z$deltas[!grepl("rho",names(z$deltas))]
    z$cov <- z$cov[!grepl("rho",rownames(z$cov)),
                    !grepl("rho",colnames(z$cov))]
  }
  nR <- z$nR
  nG <- z$nG
  nT <- z$nT
  p <- z$p
  W <- z$W
  if(!is.null(W)) W <- Matrix::Matrix(W)
  param_fit <- c(z$betas,z$deltas)
  cov_betas_deltas <- z$cov[!grepl("sigma",rownames(z$cov)),
                            !grepl("sigma",colnames(z$cov))]
  demean <- z$demean
  # HAY QUE ELIMINAR INTERCEPTOS SI NO HAY DEMEAN
  if (!demean) {
    row_interc <- cumsum(p) - p[1] + 1
    param_fit <- c(z$betas[-c(row_interc)],z$deltas)
    cov_betas_deltas <- cov_betas_deltas[-c(row_interc),-c(row_interc)]
    p <- p - 1
  }
  fchol_cov <- t(chol(cov_betas_deltas))
  names_betas <- names(param_fit)[names(param_fit)!=names(z$deltas)]
  names_deltas <- names(z$deltas)
  betas_fit <- param_fit[names_betas]
  deltas_fit <- param_fit[names_deltas]

  # Direct and Indirect effects
  cump <- cumsum(p)
  lBetaD_fit <- lBetaI_fit <- lBetaT_fit <- lbeta_fit <- list()
  for (i in 1:nG) {
    if (i==1) {
      lbeta_fit[[i]] <- betas_fit[1:cump[1]]
    } else {
      lbeta_fit[[i]] <- betas_fit[(cump[i-1]+1):cump[i]]
    }
  }
  for (j in 1:nG){
    S0 <- Matrix::solve(Matrix::Diagonal(nR)-deltas_fit[j]*W)
    BetaD_fit_j <- NULL
    BetaI_fit_j <- NULL
    p_j <- p[j]
    if (type=="sdm") p_j <- p_j / 2
    for (i in 1:p_j){
      beta_ij <- lbeta_fit[[j]][i]
      name_beta_ij <- names(beta_ij)
      if (type=="sar") S_ij <- as.matrix(S0*beta_ij)
      if (type=="sdm"){
        name_theta_ij <- paste("W",name_beta_ij,sep="_")
        theta_ij <- lbeta_fit[[j]][name_theta_ij]
        S1_ij <- beta_ij*diag(nR) + (W*theta_ij)
        S_ij <- as.matrix(S0 %*% S1_ij)
      }
      BetaD_fit_j <- c(BetaD_fit_j,mean(diag(S_ij))) # average direct effect
      diag(S_ij) <- 0
      BetaI_fit_j <- c(BetaI_fit_j,mean(rowSums(S_ij))) # average indirect eff.
    }
    lBetaD_fit[[j]] <- BetaD_fit_j
    names(lBetaD_fit[[j]]) <- names(lbeta_fit[[j]])[1:p_j]
    lBetaI_fit[[j]] <- BetaI_fit_j
    names(lBetaI_fit[[j]]) <- names(lbeta_fit[[j]])[1:p_j]
    lBetaT_fit[[j]] <- lBetaD_fit[[j]] + lBetaI_fit[[j]]
    names(lBetaT_fit[[j]]) <- names(lbeta_fit[[j]])[1:p_j]
  }
  #################################3 Simulation for Direct & Indirect effects
  lBetaD_sim <- lBetaI_sim <- lBetaT_sim <- list()
  for (j in 1:nG) {
    lBetaD_sim[[j]] <- matrix(NA,nrow=nsim,ncol=length(lBetaD_fit[[j]]))
    colnames(lBetaD_sim[[j]]) <- names(lBetaD_fit[[j]])
    lBetaI_sim[[j]] <- matrix(NA,nrow=nsim,ncol=length(lBetaI_fit[[j]]))
    colnames(lBetaI_sim[[j]]) <- names(lBetaI_fit[[j]])
  }
  lbeta_sim <- list()

  for(k in 1:nsim){
    param_sim <- param_fit + fchol_cov %*% rnorm(n=length(param_fit))
    deltas_sim <- param_sim[names_deltas,]
    for (j in 1:nG)
    {
      lbeta_sim[[j]] <- param_sim[names(lbeta_fit[[j]]),]
       p_j <- p[j]
       if (type=="sdm") p_j <- p_j / 2
      S0 <- Matrix::solve(Matrix::Diagonal(nR)-deltas_sim[j]*W)
      BetaD_sim_j <- NULL
      BetaI_sim_j <- NULL
      for (i in 1:p_j)
      {
        beta_sim_ij <- lbeta_sim[[j]][i]
        name_beta_ij <- names(beta_sim_ij)
        if (type=="sar") S_ij <- as.matrix(S0*beta_sim_ij)
        if (type=="sdm"){
          name_theta_ij <- paste("W",name_beta_ij,sep="_")
          theta_sim_ij <- lbeta_sim[[j]][name_theta_ij]
          S1_ij <- beta_sim_ij*diag(nR) + (W*theta_sim_ij)
          S_ij <- as.matrix(S0 %*% S1_ij)
        }
        BetaD_sim_j <- c(BetaD_sim_j,mean(diag(S_ij))) # average direct effect
        diag(S_ij) <- 0
        BetaI_sim_j <- c(BetaI_sim_j,mean(rowSums(S_ij))) # average indirect eff
      }
      lBetaD_sim[[j]][k,] <- BetaD_sim_j
      lBetaI_sim[[j]][k,] <- BetaI_sim_j
    }
  }
  for(j in 1:nG){
    lBetaT_sim[[j]] <- lBetaD_sim[[j]] + lBetaI_sim[[j]]
    colnames(lBetaT_sim[[j]]) <- colnames(lBetaD_sim[[j]])
  }
  tableBetaD <- matrix(NA,nrow=length(betas_fit),ncol=4)
  rownames(tableBetaD) <- names(betas_fit)
  colnames(tableBetaD) <- c("mean","sd","t-stat","p-val")
  tableBetaI <- tableBetaT <- tableBetaD
  for(j in 1:nG) {
    meanD_j <- apply(lBetaD_sim[[j]],2,"mean")
    meanI_j <- apply(lBetaI_sim[[j]],2,"mean")
    meanT_j <- apply(lBetaT_sim[[j]],2,"mean")
    sdD_j <- apply(lBetaD_sim[[j]],2,"sd")
    sdI_j <- apply(lBetaI_sim[[j]],2,"sd")
    sdT_j <- apply(lBetaT_sim[[j]],2,"sd")
    names_j <- names(meanD_j)
    tableBetaD[names_j,"mean"] <- meanD_j
    tableBetaD[names_j,"sd"] <- sdD_j
    tableBetaI[names_j,"mean"] <- meanI_j
    tableBetaI[names_j,"sd"] <- sdI_j
    tableBetaT[names_j,"mean"] <- meanT_j
    tableBetaT[names_j,"sd"] <- sdT_j
  }
  tableBetaD[,"t-stat"] <- tableBetaD[,"mean"] / tableBetaD[,"sd"]
  tableBetaI[,"t-stat"] <- tableBetaI[,"mean"] / tableBetaI[,"sd"]
  tableBetaT[,"t-stat"] <- tableBetaT[,"mean"] / tableBetaT[,"sd"]
  tableBetaD[,"p-val"] <- 2*pnorm(abs(tableBetaD[,"t-stat"]),
                                          mean=0,sd=1,lower.tail=FALSE)
  tableBetaI[,"p-val"] <- 2*pnorm(abs(tableBetaI[,"t-stat"]),
                                           mean=0,sd=1,lower.tail=FALSE)
  tableBetaT[,"p-val"] <- 2*pnorm(abs(tableBetaT[,"t-stat"]),
                                           mean=0,sd=1,lower.tail=FALSE)
  if (type=="sdm") {
    tableBetaD <- tableBetaD[!grepl("W_",rownames(tableBetaD)),
                             !grepl("W_",colnames(tableBetaD))]
    tableBetaI <- tableBetaI[!grepl("W_",rownames(tableBetaI)),
                             !grepl("W_",colnames(tableBetaI))]
    tableBetaT <- tableBetaT[!grepl("W_",rownames(tableBetaT)),
                             !grepl("W_",colnames(tableBetaT))]
  }

  #########################   PRINT TABLES
  digits <- max(3L, getOption("digits") - 3L)
  if (spsurfit$type == "sarar") type <- "sarar" # Change sarar again...
  cat("\n","\n")
  cat("Spatial SUR model type: ",type,"\n")
  cat("\n\n Direct effects","\n\n")
  printCoefmat(tableBetaT, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n\n Indirect effects","\n\n")
  printCoefmat(tableBetaI, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n\n Total effects","\n\n")
  printCoefmat(tableBetaT, P.values = TRUE, has.Pvalue = TRUE)
  res <- list(table_dir_eff = tableBetaD,
              table_ind_eff = tableBetaI,
              table_tot_eff = tableBetaT,
              sim_dir_eff = lBetaD_sim,
              sim_ind_eff = lBetaI_sim,
              sim_tot_eff = lBetaI_sim)
  res
}

