#' @name impacts.spsur
#' @rdname impacts.spsur
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
#' @export
impacts.spsur <- function(obj, ..., tr = NULL, 
                           R = NULL, listw = NULL, 
                           evalues = NULL, 
                          #useHESS = NULL NO SE USA... 
                           tol = 1e-06, empirical = FALSE, 
                           Q = NULL) {
  if (obj$type == "sim" || obj$type == "sem") 
       stop("impact measures not for this model")
  
  if (is.null(listw) && !is.null(obj$listw_style) && 
      obj$listw_style != "W") 
    stop("Only row-standardised weights supported")
  
  deltas <- obj$deltas
  type <- obj$type
  G <- obj$G; N <- obj$N; Tm <- obj$Tm
  p <- obj$p; dvars <- obj$dvars
  Sigma <- obj$Sigma
  if (type == "slx" || type == "sdem") {
    ## CREAMOS UN SLX OBJECT PARA CADA ECUACIÓN
    res <- list()
    idxcoef <- 1
    idxres <- 1
    for (i in 1:G) {
      coeffi <- obj$coefficients[idxcoef:
                                   (idxcoef + p[i] - 1)]
      rest.sei <- obj$rest.se[idxcoef:
                                (idxcoef + p[i] - 1)]
      resi <- obj$residuals[idxres: 
                              (idxres + N - 1)]
      fiti <- obj$fitted.values[idxres: 
                                  (idxres + N - 1)]
      y <- obj$y[idxres:(idxres + N - 1),1]
      Xi <- obj$X[idxres: (idxres + N - 1), names(coeffi)]
      icept <- grep("(Intercept)", colnames(Xi))
      iicept <- length(icept) > 0L
      if (iicept) {
        lmi.model <- lm(formula(paste("y ~ ", 
                         paste(colnames(Xi)[c(-icept)], 
                         collapse = "+"))), 
                       data = as.data.frame(Xi))
      } else {
        lmi.model <- lm(formula(paste("y ~ 0 + ", 
                                paste(colnames(Xi), 
                                      collapse = "+"))), 
                        data = as.data.frame(Xi))
      }
      lmi.model$coefficients <- coeffi
      lmi.model$residuals <- resi
      lmi.model$fitted.values <- fiti
      lagcept <- grep("lag.", names(lmi.model$coefficients))
      nclti <- names(lmi.model$coefficients)[-lagcept]
      mixedImps <- NULL
      Ki <- ifelse(isTRUE(grep("Intercept", 
                  names(coefficients(lmi.model))[1]) == 
                         1L), 2, 1)
      mi <- length(coefficients(lmi.model))
      odd <- (mi%/%2) > 0
      if (odd) {
        m2i <- (mi - 1)/2
      }
      if (Ki == 1 && odd) {
        warning("model configuration issue: no total impacts")
      } else {
        cmi <- matrix(0, ncol = mi, nrow = m2i)
        if (Ki == 2) {
          if (odd) {
            rownames(cmi) <- nclti[2:(m2i + 1)]
          } else {
            rownames(cmi) <- nclti[1:m2i]
          }
          for (j in 1:m2i) cmi[j, c(j + 1, 
                                   j + (m2i + 1))] <- 1
          dirImpsi <- cbind(coeffi[2:(m2i + 1)],
                            rest.sei[2:(m2i + 1)])
          rownames(dirImpsi) <- rownames(cmi)
          indirImpsi <- cbind(coeffi[(m2i + 2):mi],
                            rest.sei[(m2i + 2):mi])
          rownames(indirImpsi) <- rownames(cmi)
        } else {
          rownames(cmi) <- nclti[1:m2i]
          for (j in 1:m2i) cmi[j, c(j, j + m2i)] <- 1
          dirImpsi <- cbind(coeffi[1:m2i],
                            rest.sei[1:m2i])          
          rownames(dirImpsi) <- rownames(cmi)
          indirImpsi <- cbind(coeffi[(m2i + 1):mi],
                            rest.sei[(m2i + 1):mi])          
          rownames(indirImpsi) <- rownames(cmi)
        }
        totImpsi <- as.matrix(gmodels::estimable(lmi.model, 
                                  cmi)[, 1:2, drop = FALSE])
      }
      mixedImps <- list(dirImps = dirImpsi, 
                        indirImps = indirImpsi, 
                        totImps = totImpsi)
      attr(lmi.model, "mixedImps") <- mixedImps
      attr(lmi.model, "dvars") <- dvars[i]
      class(lmi.model) <- c("SLX", class(lmi.model))
      impWXi <- spatialreg::impacts.SLX(lmi.model)
      res[[i]] <- impWXi
      idxcoef <- idxcoef + p[i] 
      idxres <- idxres + N
    }
    #return(res)
  } else {
    if (type == "slm") { type <- "lag" 
    } else if (type == "sdm") { type <- "mixed" 
    } else if (type == "sarar") { type <- "sac" } 
    if (any(type == c("lag","mixed"))) {
      rho <- deltas
    } else if (type == "sac") {
      rho <- deltas[1:G]
    }
    beta <- obj$coefficients
    Sigmaresids <- obj$Sigma
    resvar <- obj$resvar
    res <- list()
    idx <- 1
    for (i in 1:G) {
      #s2_i <- Sigmaresids[i,i]
      beta_i <- beta[idx:(idx + p[i] - 1)]
      rho_i <- rho[i]
      #lambda_i <- lambda[i]
      interval <- obj$interval
      if (is.null(interval)) interval <- c(-1, 0.999)
      irho <- 1  ## HIPÓTESIS: NO SE INCLUYEN SIGMAS EN RESVAR
      drop2beta <- 1
      #if (type == "sac") drop2beta <- c(drop2beta, 2)
      icept <- grep("(Intercept)", names(beta_i))
      iicept <- length(icept) > 0L
      #zero_fill <- NULL
      dvars_i <- dvars[i]
      if (type == "lag" || type == "sac") {
        if (iicept) {
          P_i <- matrix(beta_i[-icept], ncol = 1)
          bnames_i <- names(beta_i[-icept])
        } else {
          P_i <- matrix(beta_i, ncol = 1)
          bnames_i <- names(beta_i)
        }
        p_i <- length(beta_i)
      } else if (type == "mixed") {
    #if (!is.null(dvars_i)) zero_fill <- attr(dvars_i, "zero_fill")
        if (iicept) {
          b1_i <- beta_i[-icept]
        } else {
          b1_i <- beta_i
        }
        p_i <- length(b1_i)
        if (p_i %% 2 != 0)
          stop("non-matched coefficient pairs")
        P_i <- cbind(b1_i[1:(p_i/2)], b1_i[((p_i/2) + 1):p_i])
        bnames_i <- names(b1_i[1:(p_i/2)])
      }  
      n <- N
      mu_i <- NULL
      Sigma_i <- NULL
      if (!is.null(R)) {
        mu_i <- c(rho_i, beta_i)
        Sigma_i <- resvar[c(names(rho_i),names(beta_i)),
                          c(names(rho_i),names(beta_i))]
      }  
      res_i <- spatialreg::intImpacts(rho = rho_i, 
                                      beta = beta_i, 
                                      P = P_i, n = n, 
                                      mu = mu_i, 
                                      Sigma = Sigma_i, 
                                      irho = irho, 
                                      drop2beta = drop2beta, 
                                      bnames = bnames_i, 
                                      interval = interval, 
                                      type = type, 
                                      tr = tr, R = R, 
                                      listw = listw, 
                                      evalues = evalues, 
                                      tol = tol, 
                                      empirical = empirical, 
                                      Q = Q, 
                                      icept = icept, 
                                      iicept = iicept, 
                                      p = p_i, 
                                      zero_fill = NULL, 
                                      dvars = dvars_i)
      idx <- idx + p[i] ## Update Index
      attr(res_i, "iClass") <- class(obj)
      res[[i]] <- res_i
    }
  }
  return(res)
}
  
