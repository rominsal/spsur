#' @name impactspsur
#' @rdname impactspsur
#'
#' @title Direct, indirect and total effects estimated for a spatial SUR model
#'
#' @description
#'   This function is a wrapper for \code{\link[spatialreg]{impacts}} method
#'   used in \pkg{spatialreg} package. Nevertheless, in this case the same 
#'   method is used for both \code{\link[spatialreg]{lagsarlm}} and 
#'   \code{\link[spatialreg]{lmSLX}} objects. For details of implementation, 
#'   see the documentation of \code{\link[spatialreg]{impacts}} function in 
#'   \pkg{spatialreg} package. \cr
#'   The function obtains the multiplier effects, on the explained variable, 
#'   of a change in a regressor for the model that has been estimated. 
#'   For reasons given below, this function only applies to models with an 
#'   autoregressive structure ("slm", "sdm", "sarar" and "gnm") or with spatial lags 
#'   of the regressors ("slx", "sdem"). \cr
#'   The measurement of the multiplier effects is a bit more complicated than 
#'   in a pure time series context because, due to the spatial structure of 
#'   the model, part of the impacts spills over non uniformly over the space. 
#'   Using the notation introduced by LeSage and Pace (2009) we distinguish 
#'   between:
#'   \itemize{
#'     \item \strong{Average Direct effects}: The average over the \emph{N} 
#'     spatial units and \emph{Tm} time periods of the effect of a unitary 
#'     change in the value of a explanatory variable on the contemporaneous 
#'     value of the corresponding explained variable, located in the same 
#'     point of the intervened regressor. This calculus is solved for all the 
#'     regressors that appear in the \emph{G} equations of the model.
#'     \item \strong{Average Indirect effects}: The average over the \emph{N} 
#'     spatial units and \emph{Tm} time periods of the effects of a unitary 
#'     change in the value of a explanatory variable on the contemporaneous 
#'     value of the corresponding explained variable, located in a different 
#'     spatial unit that that of the intervened regressor. This calculus is 
#'     solved for all the regressors that appear in the \emph{G} equations of 
#'     the model.
#'     \item \strong{Average total effects}: The sum of Direct and 
#'     Indirect effects.
#'   }
#'
#'    The information on the three estimated effects is supplement with an 
#'    indirect measure of statistical significance obtained from the 
#'    randomization approach introduced in LeSage and Pace (2009).
#'    
#' @usage impactspsur (obj, ..., tr = NULL, R = NULL, listw = NULL, 
#'                       evalues = NULL,tol = 1e-06, 
#'                       empirical = FALSE, Q = NULL)     
#'
#' @param obj An \code{spsur} object created by \code{\link{spsurml}},
#'            \code{\link{spsur3sls}} or \code{\link{spsurtime}}.
#' @inheritParams spatialreg::impacts            
#'
#' @details
#'  LeSage and Pace (2009) adapt the classical notion of 
#'  \emph{'economic multiplier'} to the problem of measuring the impact that 
#'  a unitary change in the value of a regressor, produced in a certain point 
#'  in space, has on the explained variable. The question is interesting 
#'  because, due to the spatial structure of the model, the impacts of such 
#'  change spill non uniformly over the space. In fact, the reaction of the 
#'  explained variable depends on its relative location in relation to the 
#'  point of intervention. \cr
#'  To simplify matters, LeSage and Pace (2009) propose to obtain aggregated 
#'  multipliers for each regressor, just averaging the \eqn{N^{2}} impacts 
#'  that results from intervening the value of each regressor on each of the
#'  N points in Space, on the explained variable, measured also in each of 
#'  the \eqn{N} points in space. This aggregated average is the so-called 
#'  \emph{Total effect}. \cr
#'  Part of this impact will be absorved by the explained variable located in 
#'  the same point of the regressor whose value has been changed (for example, 
#'  the k-th regresor in the g-th equation, in the n-th spatial unit) or,
#'  in other words, we expect that \eqn{[d y_{tgn}]/[d x_{ktgn}] ne 0}. The 
#'  aggregated average for the \emph{N} points in space (n=1,2,...,N) and 
#'  \emph{Tm} time periods is the so-called \emph{Direct effect}.
#'  The difference between the \emph{Total effect} and the \emph{Direct effect} 
#'  measures the portion of the impact on the explained variable that leakes 
#'  to other points in space, \eqn{[d y_{tgn}]/[d x_{ktgm}] for n ne m};
#'  this is the \emph{Indirect effect}.
#'
#'  \code{\link{impacts}} obtains the three multipliers together with an 
#'  indirect measure of statistical significance, according to the 
#'  randomization approach described in Lesage and Pace (2009). Briefly, they 
#'  suggest to obtain a sequence of \emph{nsim} random matrices of order 
#'  \emph{(NTmxG)} from a multivariate normal distribution 
#'  N(0; \emph{Sigma}), being \emph{Sigma} the estimated covariance matrix 
#'  of the \emph{G} equations in the SUR model. These random matrices, 
#'  combined with the observed values of the regressors and the estimated 
#'  values of the parameters of the corresponding spatial SUR model, are used 
#'  to obtain simulated values of the explained variables. Then, for each one 
#'  of the \emph{nsim} experiments, the SUR model is estimated, and the effects
#'  are evaluated. The function \code{\link{impacts}} obtains the standard 
#'  deviations of the \emph{nsim} estimated effects in the randomization 
#'  procedure, which are used to test the significance of the estimated 
#'  effects for the original data.
#'
#'  Finally, let us note that this is a SUR model where the \emph{G} equations 
#'  are connected only through the error terms. This means that if we 
#'  intervene a regressor in equation \emph{g}, in any point is space, only 
#'  the explained variable of the same equation \emph{g} should react. The 
#'  impacts do not spill over equations. Moreover, the impact of a regressor, 
#'  intervened in the spatial unit \emph{n}, will cross the borders of this 
#'  spatial unit only if in the right hand side of the equation there are 
#'  spatial lags of the explained variables or of the regressors. In other 
#'  words, the \emph{Indirect effect} is zero for the "sim" and "sem" models. 
#'  \code{\link{impacts}} produces no output for these two models.
#'  Lastly, it is clear that all the impacts are contemporaneous because the 
#'  equations in the SUR model have no time dynamics.
#'
#' @return 
#' A list of \emph{G} objects either of class \code{lagImpact} 
#' or \code{WXImpact}.
#' 
#' For each of the G objects of the list, if no simulation is carried out 
#' the object returned is a list with:
#' \tabular{ll}{
#'   \code{direct} \tab numeric vector \cr
#'   \code{indirect} \tab numeric vector \cr
#'   \code{total} \tab numeric vector \cr
#' } 
#' and a matching \code{Qres} list attribute if {Q} was given.
#' 
#' On the other hand, for each of the G objects of the list, if 
#' simulation is carried out the object returned is a list with:
#' \tabular{ll}{  
#'   \code{res} \tab	a list with three components as for the non-simulation 
#'   case, with a matching \code{Qres} list attribute if \code{Q} was given. \cr
#'   \code{sres} \tab a list with three \code{mcmc} matrices, for the direct, 
#'   indirect and total impacts with a matching \code{Qmcmc}list attribute 
#'   if \code{Q} was given. \cr
#' }
#' 
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Bivand, R.S. and  Piras G. (2015). Comparing Implementations of 
#'        Estimation Methods for Spatial Econometrics. \emph{Journal of 
#'        Statistical Software}, 63(18), 1-36. 
#'        https://www.jstatsoft.org/v63/i18/.
#'      \item LeSage, J., and Pace, R. K. (2009). \emph{Introduction to spatial
#'        econometrics}. Chapman and Hall/CRC.
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'        
#'      \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'        effects in seemingly unrelated regressions.
#'        \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'        \url{https://doi.org/10.1080/17421772.2010.516443}
#'   }
#'
#' @seealso
#' \code{\link[spatialreg]{impacts}}, \code{\link{spsurml}}, 
#' \code{\link{spsur3sls}}
#'
#' @examples
#' 
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_impactspsur, package="spsur")
#' 
#' \donttest{
#' ###############################################
#' ### PURE CROSS SECTIONAL DATA(G>1; Tm=1) ######
#' ###############################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#'  rm(list = ls()) # Clean memory
#'  data(spc)
#'  lwspc <- spdep::mat2listw(Wspc, style = "W")
#'  Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## For SLM, SDM and SARAR models the output is a list of "lagImpact" objects
#' ## See spatialreg::impacts for details.
#'  spcsur_slm <-spsurml(formula = Tformula, data = spc, 
#'                       type = "slm", listw = lwspc)
#'  summary(spcsur_slm)
#'  impacts_slm <- impactspsur(spcsur_slm, listw = lwspc, R = 1000)
#' ## Impacts equation 1
#'  summary(impacts_slm[[1]], zstats = TRUE, short = TRUE)
#' ## Impacts equation 2
#'  summary(impacts_slm[[2]], zstats = TRUE, short = TRUE)
#' ## For SLX and SDEM models the output is a list of "WXImpact" objects
#' ## See spatialreg::impacts for details.
#' ## A SUR-SLX model
#'  spcsur_slx <-spsurml(formula = Tformula, data = spc, 
#'                       type = "slx", listw = lwspc)
#'  summary(spcsur_slx)
#'  impacts_slx <- impactspsur(spcsur_slx, listw = lwspc)
#'  summary(impacts_slx[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_slx[[2]], zstats = TRUE, short = TRUE)
#' 
#' ## A SUR-SDM model
#'  spcsur_sdm <-spsurml(formula = Tformula, data = spc, 
#'                       type = "sdm", listw = lwspc)
#'  impacts_sdm <- impactspsur(spcsur_sdm, listw = lwspc, R = 1000)
#'  summary(impacts_sdm[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_sdm[[2]], zstats = TRUE, short = TRUE)
#' ## A SUR-SDM model with different spatial lags in each equation
#'  TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
#'  spcsur_sdm2 <-spsurml(formula = Tformula, data = spc, type = "sdm", 
#'                        listw = lwspc, Durbin = TformulaD)
#'  summary(spcsur_sdm2)                       
#'  impacts_sdm2 <- impactspsur(spcsur_sdm2, listw = lwspc, R = 1000)
#'  summary(impacts_sdm2[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_sdm2[[2]], zstats = TRUE, short = TRUE)
#'  ## A SUR-SLX model with different spatial lags in each equation
#'  spcsur_slx2 <-spsurml(formula = Tformula, data = spc, 
#'                       type = "slx", listw = lwspc, Durbin = TformulaD)
#'  summary(spcsur_slx2)
#'  impacts_slx2 <- impactspsur(spcsur_slx2, listw = lwspc)
#'  summary(impacts_slx2[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_slx2[[2]], zstats = TRUE, short = TRUE)
#' ### A SUR-SDEM model
#'  spcsur_sdem <-spsurml(formula = Tformula, data = spc, 
#'                       type = "sdem", listw = lwspc)
#'  impacts_sdem <- impactspsur(spcsur_sdem, listw = lwspc)
#'  summary(impacts_sdem[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_sdem[[2]], zstats = TRUE, short = TRUE)
#' 
#' ### A SUR-SARAR model
#'  spcsur_sarar <-spsurml(formula = Tformula, data = spc, 
#'                       type = "sarar", listw = lwspc,
#'                       control = list(tol = 0.01))
#'  impacts_sarar <- impactspsur(spcsur_sarar, listw = lwspc, R = 1000)
#'  summary(impacts_sarar[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_sarar[[2]], zstats = TRUE, short = TRUE)
#'  
#' ## A SUR-GNM model
#'  spcsur_gnm <-spsurml(formula = Tformula, data = spc, 
#'                       type = "gnm", listw = lwspc,
#'                       control = list(tol = 0.1))
#'  impacts_gnm <- impactspsur(spcsur_gnm, listw = lwspc, R = 1000)
#'  summary(impacts_gnm[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_gnm[[2]], zstats = TRUE, short = TRUE)
#' ## A SUR-GNM model with different spatial lags in each equation
#'  TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
#'  spcsur_gnm2 <-spsurml(formula = Tformula, data = spc, type = "gnm", 
#'                        listw = lwspc, Durbin = TformulaD,
#'                        control = list(tol = 0.1))
#'  summary(spcsur_gnm2)                       
#'  impacts_gnm2 <- impactspsur(spcsur_gnm2, listw = lwspc, R = 1000)
#'  summary(impacts_gnm2[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_gnm2[[2]], zstats = TRUE, short = TRUE)
#'  
#' # ####################################
#' # ######## G=1; Tm>1               ###
#' # ####################################
#' #
#'  rm(list = ls()) # Clean memory
#'  data(NCOVR, package="spsur")
#'  nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ### Some regions with no links...
#'  lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#'  Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ### A SUR-SLM model
#'  NCOVRSUR_slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                         type = "slm", listw = lwncovr, 
#'                         method = "Matrix", zero.policy = TRUE, 
#'                         control = list(fdHess = TRUE))
#'  summary(NCOVRSUR_slm)
#' ### Use of trW to compute.
#'  Wncovr <- as(spdep::listw2mat(lwncovr), "CsparseMatrix")
#'  trwncovr <- spatialreg::trW(Wncovr, type = "MC")
#'  impacts_NCOVRSUR_slm <- impactspsur(NCOVRSUR_slm, tr = trwncovr,
#'                                  R = 1000)
#'  summary(impacts_NCOVRSUR_slm[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_NCOVRSUR_slm[[2]], zstats = TRUE, short = TRUE)
#' ### A SUR-SDM model
#'  NCOVRSUR_sdm <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                         type = "sdm", listw = lwncovr, 
#'                         method = "Matrix", zero.policy = TRUE, 
#'                         control = list(fdHess = TRUE))
#'  impacts_NCOVRSUR_sdm <- impactspsur(NCOVRSUR_sdm, tr = trwncovr,
#'                                  R = 1000)
#'  summary(impacts_NCOVRSUR_sdm[[1]], zstats = TRUE, short = TRUE)
#'  summary(impacts_NCOVRSUR_sdm[[2]], zstats = TRUE, short = TRUE)
#' }
#' 
#' @export
impactspsur <- function(obj, ..., tr = NULL, 
                           R = NULL, listw = NULL, 
                           evalues = NULL,
                           tol = 1e-06, empirical = FALSE, 
                           Q = NULL) {
  if (obj$type == "sim" || obj$type == "sem") 
       stop("impact measures not for this model")
  
  if (is.null(listw) && !is.null(obj$listw_style) && 
      obj$listw_style != "W") 
    stop("Only row-standardised weights supported")
  deltas <- obj$deltas
  type <- obj$type
  Durbin <- obj$Durbin
  G <- obj$G; N <- obj$N; Tm <- obj$Tm
  p <- obj$p
  dvars <- obj$dvars
  Sigma <- obj$Sigma
  if (type == "slx" || type == "sdem") {
    ## SLX object for each equation
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
      y <- obj$Y[idxres:(idxres + N - 1),1]
      Xi <- obj$X[idxres:(idxres + N - 1), names(coeffi)]
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
      mixedImps <- NULL
      Ki <- ifelse(isTRUE(grep("Intercept", 
                  names(coefficients(lmi.model))[1]) == 
                         1L), 2, 1)
      if (isTRUE(Durbin)) {
        lagcept <- grep("lag.", names(lmi.model$coefficients))
        nclti <- names(lmi.model$coefficients)[-lagcept]
        mi <- length(coefficients(lmi.model))
        odd <- (mi %/% 2) > 0
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
      } else if (inherits(Durbin, "formula")) {
        mi <- sum(dvars[[i]])
        m2i <- dvars[[i]][2]
        cmi <- matrix(0, ncol = mi, nrow = m2i)
        inds <- attr(dvars[[i]], "inds")
        for (j in 1:m2i) {
          cmi[j, c(inds[j], j + dvars[[i]][1])] <- 1
        }
        lagcept <- grep("lag.", names(lmi.model$coefficients))
        xni <- names(lmi.model$coefficients)[-lagcept]
        if (iicept) xni <- xni[-1]
        wxni <- names(lmi.model$coefficients)[lagcept]
        wxni <- substring(wxni, nchar("lag") + 2, nchar(wxni))
        zero_fill <- length(xni) + (which(!(xni %in% wxni)))
        rownames(cmi) <- wxni
        dirImpsi <- cbind(coeffi[2:dvars[[i]][1]],
                          rest.sei[2:dvars[[i]][1]])
        rownames(dirImpsi) <- xni
        indirImpsi <- cbind(coeffi[(dvars[[i]][1] + 1):mi],
                            rest.sei[(dvars[[i]][1] + 1):mi])
        if (!is.null(zero_fill)) {
          if (length(zero_fill) > 0L) {
            lres <- vector(mode = "list", length = 2L)
            for (j in 1:2) {
              jindirImpsi <- rep(as.numeric(NA), (dvars[[i]][1] - 
                                                   1))
              inds <- attr(dvars[[i]],"inds")
              for (k in seq(along = inds)) {
                jindirImpsi[(inds[k] - 1)] <- indirImpsi[k,j]
              }
              lres[[j]] <- jindirImpsi
            }
            indirImpsi <- do.call("cbind", lres)
          }
        }
        rownames(indirImpsi) <- xni

        totImpsi <- as.matrix(gmodels::estimable(lmi.model, 
                                                 cmi)[, 1:2, drop = FALSE])
        if (!is.null(zero_fill)) {
          if (length(zero_fill) > 0L) {
            lres <- vector(mode = "list", length = 2L)
            for (j in 1:2) {
              jtotImpsi <- dirImpsi[, j]
              for (k in seq(along = inds)) {
                jtotImpsi[(inds[k] - 1)] <- totImpsi[k, j]
              }
              lres[[j]] <- jtotImpsi
            }
            totImpsi <- do.call("cbind", lres)
          }
        }
        rownames(totImpsi) <- xni        
        
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
    } else if (any(type == c("sdm", "gnm"))) { 
      if (type == "gnm") deltas <- deltas[1:G] # supress lambda's in gnm
      type <- "mixed" 
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
      betai <- beta[idx:(idx + p[i] - 1)]
      rho_i <- rho[i]
      #lambda_i <- lambda[i]
      interval <- obj$interval
      if (is.null(interval)) interval <- c(-1, 0.999)
      irho <- 1  ## HIPÓTESIS: NO SE INCLUYEN SIGMAS EN RESVAR
      drop2beta <- 1
      #if (type == "sac") drop2beta <- c(drop2beta, 2)
      icept <- grep("(Intercept)", names(betai))
      iicept <- length(icept) > 0L
      zero_fill <- NULL
      dvarsi <- dvars[[i]]
      if (type == "lag" || type == "sac") {
        if (iicept) {
          Pi <- matrix(betai[-icept], ncol = 1)
          bnamesi <- names(betai[-icept])
        } else {
          Pi <- matrix(betai, ncol = 1)
          bnamesi <- names(betai)
        }
        pi <- length(betai)
      } else if (type == "mixed") {
        if (!is.null(dvarsi)) zero_fill <- attr(dvarsi, "zero_fill")
        if (iicept) {
          b1i <- betai[-icept]
        } else {
          b1i <- betai
        }
        if (!is.null(zero_fill)) {
          if (length(zero_fill) > 0L) {
            inds <- attr(dvarsi, "inds")
            b1i_long <- rep(0, 2 * (dvarsi[1] - 1))
            b1i_long[1:(dvarsi[1] - 1L)] <- b1i[1:
                                                (dvarsi[1] - 1)]
            names(b1i_long)[1:(dvarsi[1] - 1L)] <- names(b1i)[1:
                                                      (dvarsi[1] - 1)]
            for (j in seq(along = inds)) {
              b1i_long[(dvarsi[1] - 1L) + (inds[j] - 1L)] <- b1i[(dvarsi[1] - 
                                                                 1L) + j]
            }
            b1i <- b1i_long
          }
        }
        pi <- length(b1i)
        if (pi %% 2 != 0)
          stop("non-matched coefficient pairs")
        Pi <- cbind(b1i[1:(pi/2)], b1i[((pi/2) + 1):pi])
        bnamesi <- names(b1i[1:(pi/2)])
      }  
      #n <- N*Tm
      mu_i <- NULL
      Sigma_i <- NULL
      if (!is.null(R)) {
        mu_i <- c(rho_i, betai)
        Sigma_i <- resvar[c(names(rho_i),names(betai)),
                          c(names(rho_i),names(betai))]
      }  
      res_i <- spatialreg::intImpacts(rho = rho_i, 
                                      beta = betai, 
                                      P = Pi, n = N, 
                                      mu = mu_i, 
                                      Sigma = Sigma_i, 
                                      irho = irho, 
                                      drop2beta = drop2beta, 
                                      bnames = bnamesi, 
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
                                      p = pi, 
                                      zero_fill = zero_fill, 
                                      dvars = dvarsi)
      idx <- idx + p[i] ## Update Index
      attr(res_i, "iClass") <- class(obj)
      res[[i]] <- res_i
    }
  }
  return(res)
}
  
