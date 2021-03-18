#' @name lmtestspsur
#' @rdname lmtestspsur
#' @param ... further arguments passed to the methods.  
#' @export
#' 
lmtestspsur <- function(...) {
  UseMethod("lmtestspsur")
}

#' @name lmtestspsur
#' @rdname lmtestspsur
#' @title Testing for the presence of spatial effects in Seemingly 
#'  Unrelated Regressions
#' @description The function \code{\link{spsurml}}  reports a collection of 
#'  Lagrange Multipliers designed to test  for the presence of different 
#'  forms of spatial dependence in a \emph{SUR} model of the "sim" type. 
#'  That is, the approach of this function is from 
#'  \emph{'specific to general'}. As said, the model of the null hypothesis 
#'  is the "sim" model whereas the model of the alternative depends on 
#'  the effect whose omission we want to test.
#'
#' The collection of Lagrange Multipliers obtained by \code{lmtestspsur} 
#' are standard in the literature and take into account the multivariate 
#' nature of the \emph{SUR} model. As a limitation, note that each 
#' Multiplier tests for the omission of the same spatial effects in all 
#' the cross-sections of the \emph{G} equations.
#' 
#' @details \code{\link{lmtestspsur}} tests for the omission of spatial 
#'  effects in the "sim" version of the \emph{SUR} model: \cr
#'
#'     \deqn{y_{tg} = X_{tg} \beta_{g} + u_{tg}}
#'     \deqn{E[u_{tg}u_{th}']= \sigma_{gh}I_{N}  \quad E[u_{tg}u_{sh}']= 0 
#'           \mbox{ if } t ne s}
#'
#' where \eqn{y_{tg}} and \eqn{u_{tg}} are \emph{(Nx1)} vectors, 
#' corresponding to the g-th equation and time period t;
#' \eqn{X_{tg}} is the matrix of exogenous variables, of order 
#' \emph{\eqn{(Nxp_{g})}}. Moreover, \eqn{\beta_{g}} is an unknown
#' \emph{\eqn{(p_{g}x1)}} vector of coefficients and \eqn{\sigma_{gh}I_{N}} 
#' the covariance between equations \emph{g} and \emph{h},
#' being \eqn{\sigma_{gh}} and scalar and \eqn{I_{N}} the identity 
#' matrix of orden N.
#'
#'
#' The Lagrange Multipliers reported by this function are the followings:
#'
#'   \itemize{
#'     \item \strong{LM-SUR-LAG}: Tests for the omission of a spatial lag of 
#'      the explained variable in the right hand side of the "sim" equation. 
#'      The model of the alternative is: \cr
#'
#'     \eqn{y_{tg} = \rho_{g}Wy_{tg} + X_{tg} \beta_{g} + u_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'          \eqn{H_{0}: \rho_{g}=0 (forall g)} vs  
#'          \eqn{H_{A}: \rho_{g} ne 0 (exist g)}
#'
#'      \item \strong{LM-SUR-ERR}: Tests for the omission of spatial 
#'       dependence in the equation of the errors of the "sim" model. The 
#'       model of the alternative is:
#'
#'     \eqn{y_{tg} = X_{tg} \beta_{g} + u_{tg}}; 
#'     \eqn{u_{tg}= \lambda_{g}Wu_{tg}+\epsilon_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'          \eqn{H_{0}: \lambda_{g}=0 (forall g)} vs  
#'          \eqn{H_{A}: \lambda_{g}  ne 0 (exist g)}
#'
#'      \item \strong{LM-SUR-SARAR}: Tests for the simultaneous omission of 
#'      a spatial lag of the explained variable in the right hand side of 
#'      the "sim" equation and spatial dependence in the equation of the 
#'      errors. The model of the alternative is:
#'
#'
#'     \eqn{y_{tg} = \rho_{g}Wy_{tg}+X_{tg} \beta_{g} + u_{tg}}; 
#'     \eqn{u_{tg}= \lambda_{g}Wu_{tg}+\epsilon_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'      \eqn{H_{0}: \rho_{g}=\lambda_{g}=0 (forall g)} vs  
#'      \eqn{H_{A}: \rho_{g} ne 0 or \lambda_{g} ne 0 (exist g)}
#'
#'      \item
#'      \strong{LM*-SUR-SLM} and \strong{LM*-SUR-SEM}: These two test are 
#'       the robustifyed version of the original, raw Multipliers, 
#'       \strong{LM-SUR-SLM} and \strong{LM-SUR-SEM}, which can be severely 
#'       oversized if the respective alternative hypothesis is misspeficied 
#'       (this would be the case if, for example, we are testing for omitted 
#'       lags of the explained variable whereas the problem is that there is 
#'       spatial dependence in the errors, or viceversa). The null and 
#'       alternative hypotheses of both test are totally analogous to their
#'       twin non robust Multipliers.
#'     }
#' @param time time index for the spatial panel SUR data. 
#' @param ... further arguments passed to the method.  
#'  Default = \code{NULL}
#' @inheritParams spsurml 
#' @return A list of \code{htest} objects each one including the Wald
#'   statistic, the corresponding p-value and the degrees of
#'   freedom.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#
#'
#' @references
#'   \itemize{
#'     \item Mur, J., López, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions.
#'       \emph{Spatial Economic Analysis}, 5(4), 399-440.
#'       \url{https://doi.org/10.1080/17421772.2010.516443}
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'        \url{https://doi.org/10.1007/s00168-014-0624-2}
#'        \item Anselin, L. (1988) A test for spatial autocorrelation in seemingly unrelated 
#'        regressions \emph{Economics Letters} 28(4), 335-341.  
#'        \item Anselin, L. (1988) \emph{Spatial econometrics: methods and models} 
#'        Chap. 9 Dordrecht
#'        \item Anselin, L. (2016) Estimation and Testing in the Spatial Seemingly 
#'       Unrelated Regression (SUR). \emph{Geoda Center for Geospatial Analysis 
#'       and Computation, Arizona State University}. Working Paper 2016-01.
#'       \url{http://dx.doi.org/10.13140/RG.2.2.15925.40163}
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{anova}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1) # #######
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data("spc")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' lwspc <- spdep::mat2listw(Wspc, style = "W")
#' lmtestspsur(formula = Tformula, data = spc, listw = lwspc)
#'
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_lmtestspsur, package="spsur")
#' 
#' \donttest{
#' #################################################
#' ######## PANEL DATA (G>1; Tm>1)          ########
#' #################################################
#'
#' #### Example 2: Homicides & Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for
#' # continental U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' data(NCOVR, package="spsur")
#' nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ### Some regions with no links...
#' lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#' ### With different number of exogenous variables in each equation
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 | PS80 + UE80 +RD80 |
#'             PS90 + UE90 + RD90 + PO90
#' lmtestspsur(formula = Tformula, data = NCOVR.sf, 
#'             listw = lwncovr)
#'
#' #################################################################
#' ######### PANEL DATA: TEMPORAL CORRELATIONS (G=1; Tm>1) ########
#' #################################################################
#' 
#' ##### Example 3: NCOVR in panel data form
#' Year <- as.numeric(kronecker(c(1960,1970,1980,1990), 
#'                    matrix(1,nrow = dim(NCOVR.sf)[1])))
#' HR <- c(NCOVR.sf$HR60,NCOVR.sf$HR70,NCOVR.sf$HR80,NCOVR.sf$HR90)
#' PS <- c(NCOVR.sf$PS60,NCOVR.sf$PS70,NCOVR.sf$PS80,NCOVR.sf$PS90)
#' UE <- c(NCOVR.sf$UE60,NCOVR.sf$UE70,NCOVR.sf$UE80,NCOVR.sf$UE90)
#' NCOVRpanel <- as.data.frame(cbind(Year,HR,PS,UE))
#' Tformula <- HR ~ PS + UE
#' lmtestspsur(formula = Tformula, data = NCOVRpanel, time = Year, 
#' listw = lwncovr)
#' }
#' 
#' @export
lmtestspsur.formula <- function(formula, data, listw, na.action, 
                                time = NULL, Tm = 1, zero.policy = NULL, 
                                R = NULL, b = NULL, ...) {
  if (!inherits(formula, "Formula")) formula <- Formula::Formula(formula)
  cl <- match.call()
  if (is.null(time)) {
    mt <- terms(formula, data = data)
    mf <- lm(formula, data = data, na.action = na.action, 
             method = "model.frame")
    mf$drop.unused.levels <- TRUE
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
      subset <- !(1:length(listw$neighbours) %in% na.act)
      listw <- subset(listw, subset, zero.policy = zero.policy)
    }
    W <- Matrix::Matrix(spdep::listw2mat(listw))
    get_XY <- get_data_spsur(formula = formula, mf = mf, 
                             Durbin = FALSE,
                             listw = listw, 
                             zero.policy = zero.policy, 
                             N = N, Tm = Tm)
    Y <- get_XY$Y
    X <- get_XY$X
    G <- get_XY$G
    N <- get_XY$N
    Tm <- get_XY$Tm
    p <- get_XY$p
    dvars <- get_XY$dvars
    rm(get_XY)
    if (length(p) == 1) p <- rep(p,G)    
  } else { #G = 1 and Tm > 1 (Panel data) 
    if(!is.factor(time)) time <- as.factor(time)
    time <- droplevels(time)
    if (length(time) != nrow(data)) 
      stop("time must have same length than the number of rows in data")
    mt <- terms(formula)
    G <- length(levels(time))
    Ylist <- vector("list", G)
    Xlist <- vector("list", G)
    p <- NULL
    namesX <- NULL
    levels_time <- levels(time)
    for (i in 1:G) {
      data_i <- model.frame(mt, data = data[time == levels_time[i],])
      Ylist[[i]] <- data_i[, 1]
      Xlist[[i]] <- model.matrix(mt, data = data[time==levels_time[i],])
      p <- c(p,ncol(Xlist[[i]]))
      namesX <- c(namesX, paste(colnames(Xlist[[i]]), 
                                i, sep="_"))
    }
    Y <- matrix(unlist(Ylist), ncol=1)
    X <- as.matrix(Matrix::bdiag(Xlist))
    colnames(X) <- namesX
    N <- length(Ylist[[1]]); Tm <- 1
  }
  res <- lmtestspsur.default(Y = Y, X = X, G = G, 
                             N = N, Tm = Tm, listw = listw, 
                             p = p, R = R, b = b) 
  for (i in 1:length(res)) {
    res[[i]]$data.name <- cl[[3]]
  }
  return(res)
}

#' @name lmtestspsur
#' @rdname lmtestspsur
#' @param ... further arguments passed to the method.  
#' @inheritParams spsurml 
#' @export
lmtestspsur.default <- function(Y, X, G, N, Tm = 1, listw, 
                                p, R = NULL, b = NULL, ...) {
  if (!inherits(listw,c("listw","Matrix","matrix")))
    stop("listw format unknown")
  if (inherits(listw, "listw")) {
    W <- Matrix::Matrix(spdep::listw2mat(listw))
  }
  if (inherits(listw, "matrix")) {
    W <- Matrix::Matrix(listw)
    listw <- spdep::mat2listw(W)
  }  
  if (inherits(listw, "Matrix")) {
    W <- listw
    listw <- spdep::mat2listw(as.matrix(W))
  }
  if (Tm > 1 && G == 1){
    # Change dimensions in this case with matrix Data
    G <- Tm
    Tm <- 1
  }
  if (!is.null(R) && !is.null(b)) {
    restr <- X_restr(X = X, R = R, b = b, p = p)
    X <- restr$Xstar
    p <- restr$pstar
  }
  res <- sur3_spdiag(Tm = Tm, G = G, N = N, Y = Y, X = X, W = W)
  return(res)
}
