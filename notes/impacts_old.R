#' @name impacts_old
#' @rdname impacts_old
#'
#' @title Direct, indirect and total effects estimated for a spatial SUR model
#'
#' @description
#' This function obtains the multiplier effects, on the explained variable, of a change in a regressor
#'  for the model that has been estimated. For reasons given below, this function
#'  only applies to models with an autoregressive structure (\strong{"slm"}, \strong{"sdm"} and \strong{"sarar"})
#'  or with spatial lags of the regressors (\strong{"slx"}, \strong{"sdem"}).\cr
#'  The measurement of the multiplier effects is a bit more complicated than in a pure time
#'  series context because, due to the spatial structure of the model, part of the impacts
#'  spills over non uniformly over the space. Using the notation introduced by LeSage and Pace (2009)
#'  we distinguish between:
#'
#'   \itemize{
#'     \item \strong{Average Direct effects}: The average over the \emph{N} spatial units and
#'      \emph{Tm} time periods of the effect of a unitary change in the value of a explanatory variable
#'      on the contemporaneous value of the corresponding explained variable, located in the same point of
#'      the intervened regressor. This calculus is solved for all the regressors that appear in the \emph{G}
#'      equations of the model.
#'     \item \strong{Average Indirect effects}: The average over the \emph{N} spatial units and
#'      \emph{Tm} time periods of the effects of a unitary change in the value of a explanatory variable
#'      on the contemporaneous value of the corresponding explained variable, located in a different spatial
#'      unit that that of the intervened regressor. This calculus is solved for all the regressors that appear in the \emph{G}
#'      equations of the model.
#'     \item \strong{Average total effects}: The sum of Direct and Indirect effects.
#'   }
#'
#'    The information on the three estimated effects is supplement with an indirect measure of statistical
#'    significance obtained from the randomization approach introduced in LeSage and Pace (2009).
#'
#' @param spsurfit A fitted object of class spsur.
#' @param nsim Number of simulations for the randomization procedure. Default = 1000.
#'
#' @details
#'  LeSage and Pace (2009) adapt the classical notion of \emph{'economic multiplier'} to the problem of
#'  measuring the impact that a unitary change in the value of a regressor, produced in a certain point in space,
#'  has on the explained variable. The question is interesting because, due to the spatial structure of the model,
#'  the impacts of such change spill non uniformly over the space. In fact, the reaction of the explained variable
#'  depends on its relative location in relation to the point of intervention.\cr
#'
#'  To simplify matters, LeSage and Pace (2009) propose to obtain aggregated multipliers for each regressor,
#'  just averaging the \eqn{N^{2}} impacts that results from intervening the value of each regressor on each of the
#'  N points in Space, on the explained variable, measured also in each of the \eqn{N} points in space.
#'  This aggregated average is the so-called \emph{Total effect}.\cr
#'
#'  Part of this impact will be absorved by the explained variable located in the same point of the regressor whose
#'  value has been changed (for example, the k-th regresor in the g-th equation, in the n-th spatial unit) or,
#'   in other words, we expect that \eqn{[d y_{tgn}]/[d x_{ktgn}] ne 0}. The aggregated average for the
#'  \emph{N} points in space (n=1,2,...,N) and \emph{Tm} time periods is the so-called \emph{Direct effect}.
#'  The difference between the \emph{Total effect} and the \emph{Direct effect} measures the portion of the impact
#'  on the explained variable that leakes to other points in space, \eqn{[d y_{tgn}]/[d x_{ktgm}] for n ne m};
#'  this is the \emph{Indirect effect}.
#'
#'  \code{\link{impacts}} obtains the three multipliers together with an indirect measure of statistical significance,
#'  according to the randomization approach described in Lesage and Pace (2009). Briefly, they suggest to obtain
#'  a sequence of \emph{nsim} random matrices of order \emph{(NTmxG)} from a multivariate normal distribution
#'  N(0; \strong{Sigma}), being \strong{Sigma} the estimated covariance matrix of the \emph{G} equations in the SUR
#'  model. These random matrices, combined with the observed values of the regressors and the estimated values of
#'  the parameters of the corresponding spatial SUR model, are used to obtain simulated values of the explained
#'  variables. Then, for each one of the \emph{nsim} experiments, the SUR model is estimated, and the effects
#'  are evaluated. The function \code{\link{impacts}} obtains the standard deviations of the \emph{nsim} estimated
#'  effects in the randomization procedure, which are used to test the significance of the estimated effects for the
#'  original data.
#'
#'  Finally, let us note that this is a SUR model where the \emph{G} equations are connected only through the error
#'  terms. This means that if we intervene a regressor in equation \emph{g}, in any point is space, only the explained
#'  variable of the same equation \emph{g} should react. The impacts do not spill over equations.
#'  Moreover, the impact of a regressor, intervened in the spatial unit \emph{n}, will
#'  cross the borders of this spatial unit only if in the right hand side of the equation there are spatial lags of the
#'  explained variables or of the regressors. In other words, the \emph{Indirect effect} is zero for the
#'  \strong{"sim"} and \strong{"sem"} models. \code{\link{impacts}} produces no output for these two models.
#'  Lastly, it is clear that all the impacts are contemporaneous because the equations in the SUR model
#'  have no time dynamics.
#'
#' @return Returns the Direct, Indirect and Total effects of the estimated spatial SUR model and
#'   simulted significance measures.
#'   \tabular{ll}{
#'   \code{table_dir_eff} \tab Table of average direct effects. \cr
#'   \code{table_ind_eff} \tab Table of average indirect effects. \cr
#'   \code{table_tot_eff} \tab Table of average total effects. \cr
#'   }
#'
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#'
#' @references
#'   \itemize{
#'      \item LeSage, J., and Pace, R. K. (2009). \emph{Introduction to spatial
#'        econometrics}. Chapman and Hall/CRC.
#'        \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'       \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'   }
#'
#'
#' @seealso
#'
#' \code{\link{spsurml}},  \code{\link{spsur3sls}}
#'
#' @examples
#'
#' ###############################################
#' ### PURE CROSS SECTIONAL DATA(G>1; Tm=1) ######
#' ###############################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## A SUR-SLM model.
#' spcsur.slm <-spsur3sls(Form = Tformula, data = spc, type = "slm", W = Wspc)
#' summary(spcsur.slm)
#' eff.spcsur.slm <- impacts(spcsur.slm, nsim = 30)
#' \donttest{
#' ## Each case usually takes 1-2 minutes maximum
#'
#' ## A SUR-SDM model
#' spcsur.sdm <-spsurml(Form = Tformula, data = spc, type = "sdm", W = Wspc)
#' summary(spcsur.sdm)
#' eff.spcsur.sdm <- impacts(spcsur.sdm, nsim = 300)
#'
#' ## A SUR-SLX model
#' spcsur.slx <-spsurml(Form = Tformula, data = spc, type = "slx", W = Wspc)
#' summary(spcsur.slx)
#' eff.spcsur.slx <- impacts(spcsur.slx, nsim = 300)
#'
#' ## A SUR-SDEM model
#' spcsur.sdem <-spsurml(Form = Tformula, data = spc, type = "sdem", W = Wspc)
#' summary(spcsur.sdem)
#' eff.spcsur.sdem <- impacts(spcsur.sdem, nsim = 300)
#'
#' ## A SUR-SARAR model
#' spcsur.sarar <-spsurml(Form = Tformula, data = spc, type = "sarar", W = Wspc)
#' summary(spcsur.sarar)
#' eff.spcsur.sarar <- impacts(spcsur.sarar, nsim = 300)
#'
#'
#' ####################################
#' ######## G=1; Tm>1               ###
#' ####################################
#'
#' ### Only execute if you have enough time...
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ## A SUR-SLM model
#' NCOVRSUR.slm <-spsurml(Form = Tformula, data = NCOVR, type = "slm", W = W)
#' summary(NCOVRSUR.slm)
#' eff.NCOVR.slm <- impacts(NCOVRSUR.slm, nsim = 100)
#'
#' }
#' @export
impacts_old <- function(spsurfit, nsim = 1000){
  z <- spsurfit
  type <- z$type
  if ((type == "sim" | type == "sem" )){
    stop("type must not be sim or sem")
  }
  N <- z$N
  G <- z$G
  Tm <- z$Tm
  p <- z$p
  W <- z$W
  if(!is.null(W)) W <- Matrix::Matrix(W)
  if (type == "sarar"){ # Treat sarar as a slm case...
    z$deltas <- z$deltas[!grepl("rho",names(z$deltas))]
    z$cov <- z$cov[!grepl("rho",rownames(z$cov)),
                    !grepl("rho",colnames(z$cov))]
    type <- "slm"
  }
  if (type == "slx" | type == "sdem"){
    if (type == "sdem"){
      z$deltas <- z$deltas[!grepl("rho",names(z$deltas))]
      z$cov <- z$cov[!grepl("rho",rownames(z$cov)),
                     !grepl("rho",colnames(z$cov))]
    }
    z$deltas <- rep(0,G) # set a sdm model with lambda's = 0
    for (i in 1:G){
      names(z$deltas)[i] <- paste("lambda",i,sep="_")
    }
    type <- "sdm"
  }
  param_fit <- c(z$betas,z$deltas)
  cov_betas_deltas <- z$cov[!grepl("sigma",rownames(z$cov)),
                            !grepl("sigma",colnames(z$cov))]
  demean <- z$demean
  if (!demean) {
    param_fit <- param_fit[!grepl("Intercept",names(param_fit))]
    cov_betas_deltas <- z$cov[!grepl("Intercept",rownames(cov_betas_deltas)),
                              !grepl("Intercept",colnames(cov_betas_deltas))]
    # row_interc <- 1
    # if (length(p) > 1){
    #   row_interc <- c(row_interc,(cumsum(p) + 1)[1:length(p)-1])
    # }
    # param_fit <- c(z$betas[-c(row_interc)],z$deltas)
    # cov_betas_deltas <- cov_betas_deltas[-c(row_interc),-c(row_interc)]
    p <- p - 1
  }
  fchol_cov <- t(chol(cov_betas_deltas))
  if (length(param_fit) > nrow(fchol_cov)){ # Fixed lambdas for type slx/sdem
    np <- nrow(fchol_cov)
    for (i in 1:G){
      fchol_cov <- cbind(rbind(fchol_cov,0),0)
      rownames(fchol_cov)[np+i] <- names(z$deltas)[i]
      colnames(fchol_cov)[np+i] <- names(z$deltas)[i]
    }
  }
  names_betas <- names(param_fit)[names(param_fit)!=names(z$deltas)]
  names_deltas <- names(z$deltas)
  betas_fit <- param_fit[names_betas]
  deltas_fit <- param_fit[names_deltas]

  # Direct and Indirect effects
  cump <- cumsum(p)
  lBetaD_fit <- lBetaI_fit <- lBetaT_fit <- lbeta_fit <- list()
  for (i in 1:G) {
    if (i==1) {
      lbeta_fit[[i]] <- betas_fit[1:cump[1]]
    } else {
      lbeta_fit[[i]] <- betas_fit[(cump[i-1]+1):cump[i]]
    }
  }
  for (j in 1:G){
    S0 <- Matrix::solve(Matrix::Diagonal(N)-deltas_fit[j]*W)
    BetaD_fit_j <- NULL
    BetaI_fit_j <- NULL
    p_j <- p[j]
    if (type=="sdm") p_j <- p_j / 2
    for (i in 1:p_j){
      beta_ij <- lbeta_fit[[j]][i]
      name_beta_ij <- names(beta_ij)
      if (type=="slm") S_ij <- as.matrix(S0*beta_ij)
      if (type=="sdm"){
        name_theta_ij <- paste("W",name_beta_ij,sep="_")
        theta_ij <- lbeta_fit[[j]][name_theta_ij]
        S1_ij <- beta_ij*diag(N) + (W*theta_ij)
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
  ################################# Simulation for Direct & Indirect effects
  lBetaD_sim <- lBetaI_sim <- lBetaT_sim <- list()
  for (j in 1:G) {
    lBetaD_sim[[j]] <- matrix(NA,nrow=nsim,ncol=length(lBetaD_fit[[j]]))
    colnames(lBetaD_sim[[j]]) <- names(lBetaD_fit[[j]])
    lBetaI_sim[[j]] <- matrix(NA,nrow=nsim,ncol=length(lBetaI_fit[[j]]))
    colnames(lBetaI_sim[[j]]) <- names(lBetaI_fit[[j]])
  }
  lbeta_sim <- list()

  for(k in 1:nsim){
    param_sim <- param_fit + fchol_cov %*% rnorm(n=length(param_fit))
    deltas_sim <- param_sim[names_deltas,]
    for (j in 1:G)
    {
      lbeta_sim[[j]] <- param_sim[names(lbeta_fit[[j]]),]
       p_j <- p[j]
       if (type == "sdm") p_j <- p_j / 2
      S0 <- Matrix::solve(Matrix::Diagonal(N) - deltas_sim[j]*W)
      BetaD_sim_j <- NULL
      BetaI_sim_j <- NULL
      for (i in 1:p_j)
      {
        beta_sim_ij <- lbeta_sim[[j]][i]
        name_beta_ij <- names(beta_sim_ij)
        if (type=="slm") S_ij <- as.matrix(S0*beta_sim_ij)
        if (type=="sdm"){
          name_theta_ij <- paste("W",name_beta_ij,sep="_")
          theta_sim_ij <- lbeta_sim[[j]][name_theta_ij]
          S1_ij <- beta_sim_ij*diag(N) + (W*theta_sim_ij)
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
  for(j in 1:G){
    lBetaT_sim[[j]] <- lBetaD_sim[[j]] + lBetaI_sim[[j]]
    colnames(lBetaT_sim[[j]]) <- colnames(lBetaD_sim[[j]])
  }
  tableBetaD <- matrix(NA,nrow=length(betas_fit),ncol=4)
  rownames(tableBetaD) <- names(betas_fit)
  colnames(tableBetaD) <- c("mean","sd","t-stat","p-val")
  tableBetaI <- tableBetaT <- tableBetaD
  for(j in 1:G) {
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
  if (spsurfit$type == "sarar") type <- "sarar" # Change types again...
  if (spsurfit$type == "slx") type <- "slx" # Change types again...
  if (spsurfit$type == "sdem") type <- "sdem" # Change types again...
  cat("\n")
  cat("Spatial SUR model type: ",type,"\n")
  cat("\n Direct effects","\n\n")
  printCoefmat(tableBetaD, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n Indirect effects","\n\n")
  printCoefmat(tableBetaI, P.values = TRUE, has.Pvalue = TRUE)
  cat("\n Total effects","\n\n")
  printCoefmat(tableBetaT, P.values = TRUE, has.Pvalue = TRUE)
  res <- list(table_dir_eff = tableBetaD,
              table_ind_eff = tableBetaI,
              table_tot_eff = tableBetaT,
              sim_dir_eff = lBetaD_sim,
              sim_ind_eff = lBetaI_sim,
              sim_tot_eff = lBetaI_sim)
  res
}

