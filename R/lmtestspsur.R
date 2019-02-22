#' @name lmtestspsur
#' @rdname lmtestspsur
#'
#' @title Testing for the presence of spatial effects in Seemingly Unrelated Regressions
#'
#' @description The function \code{\link{spsurml}}  reports a collection of Lagrange Multipliers
#'  designed to test  for the presence of different forms of spatial dependence in a
#'  \emph{SUR} model of the \strong{"sim"} type. That is, the approach of this function is
#'  from \emph{'specific to general'}. As said, the model of the null hypothesis is the \strong{"sim"}
#'  model whereas the model of the alternative depends on the effect whose omission we want to test.
#'
#' The collection of Lagrange Multipliers obtained by \code{lmtestspsur} are standard in the
#' literature and take into account the multivariate nature of the \emph{SUR} model. As a limitation,
#' note that each Multiplier tests for the omission of the same spatial effects in all the cross-sections of
#' the \emph{G} equations.
#'
#' @inheritParams spsurml
#' @param print_table Logical value to print the output. Default = \code{TRUE}
#' @param time Time variable.
#'
#' @details \code{\link{lmtestspsur}} tests for the omission of spatial effects in the \strong{"sim"} version
#'  of the \emph{SUR} model: \cr
#'
#'     \deqn{y_{tg} = X_{tg} \beta_{g} + u_{tg}}
#'     \deqn{E[u_{tg}u_{th}']= \sigma_{gh}I_{N} & E[u_{tg}u_{sh}']= 0 if t ne s}
#'
#' where \eqn{y_{tg}} and \eqn{u_{tg}} are \emph{(Nx1)} vectors, corresponding to the g-th equation and time period t;
#' \eqn{X_{tg}} is the matrix of exogenous variables, of order\emph{\eqn{(Nxp_{g})}}. Moreover, \eqn{\beta_{g}} is an unknown
#' \emph{\eqn{(p_{g}x1)}} vector of coefficients and \eqn{\sigma_{gh}I_{N}} the covariance between equations \emph{g} and \emph{h},
#' being \eqn{\sigma_{gh}} and scalar and \eqn{I_{N}} the identity matrix of orden N.
#'
#'
#' The Lagrange Multipliers reported by this function are the followings:
#'
#'   \itemize{
#'     \item \strong{LM-SUR-LAG}: Tests for the omission of a spatial lag of the explained variable
#'      in the right hand side of the \strong{"sim"} equation. The model of the alternative is:\cr
#'
#'     \eqn{y_{tg} = \lambda_{g}Wy_{tg}+X_{tg} \beta_{g} + u_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'          \eqn{H_{0}: \lambda_{g}=0 (forall g)} vs  \eqn{H_{A}: \lambda_{g} ne 0 (exist g)}
#'
#'      \item \strong{LM-SUR-ERR}: Tests for the omission of spatial dependence in the equation of the errors
#'      of the \strong{"sim"} model. The model of the alternative is:
#'
#'     \eqn{y_{tg} = X_{tg} \beta_{g} + u_{tg}}; \eqn{u_{tg}= \rho_{g}Wu_{tg}+\epsilon_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'          \eqn{H_{0}: \rho_{g}=0 (forall g)} vs  \eqn{H_{A}: \rho_{g}  ne 0 (exist g)}
#'
#'      \item \strong{LM-SUR-SARAR}: Tests for the simultaneous omission of a spatial lag of the explained
#'       variable in the right hand side of the \strong{"sim"} equation and spatial dependence in the
#'       equation of the errors. The model of the alternative is:
#'
#'
#'     \eqn{y_{tg} = \lambda_{g}Wy_{tg}+X_{tg} \beta_{g} + u_{tg}}; \eqn{u_{tg}= \rho_{g}Wu_{tg}+\epsilon_{tg}}
#'
#'       The null and alternative hypotheses are:
#'
#'      \eqn{H_{0}: \lambda_{g}=\rho_{g}=0 (forall g)} vs  \eqn{H_{A}: \lambda_{g} ne 0 or \rho_{g} ne 0 (exist g)}
#'
#'      \item
#'      \strong{LM*-SUR-SLM} and \strong{LM*-SUR-SEM}: These two test are the robustifyed version of the original,
#'      raw Multipliers, \strong{LM-SUR-SLM} and \strong{LM-SUR-SEM}, which can be severely oversized if
#'      the respective alternative hypothesis is misspeficied (this would be the case if, for example, we are
#'      testing for omitted lags of the explained variable whereas the problem is that there is spatial dependence
#'      in the errors, or viceversa). The null and alternative hypotheses of both test are totally analogous to  their twin
#'      non robust Multipliers.
#'     }
#'
#' @return A list including:
#'   \tabular{ll}{
#'   \code{stat_names}   \tab Name of Lagrange Multiplier. \cr
#'   \code{stat}   \tab Value of the corresponding Lagrange Multiplier. \cr
#'   \code{df} \tab Degrees of freedom for each Multiplier.\cr
#'   }
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
#'      \item López, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1),
#'        197-220.
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{lrtestspsur}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' data("spc")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' LMs <- lmtestspsur(Form = Tformula, data = spc, W = Wspc)
#'
#' #################################################
#' ######## PANEL DATA (G>1; Tm>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides & Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for
#' # continental U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' data("NCOVR")
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 |PS80 + UE80 | PS90 + UE90
#' LMs <- lmtestspsur(Form = Tformula, data = NCOVR, W = W)
#'
#' # With different number of exogenous variables in each equation
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 | PS80 + UE80 +RD80 |
#'             PS90 + UE90 + RD90 + PO90
#' LMs <- lmtestspsur(Form = Tformula, data = NCOVR, W = W)
#'
#' ################################################################
#' ######## PANEL DATA: TEMPORAL CORRELATIONS (G=1; Tm>1) ########
#' ################################################################
#'
#' #### Example 3: NCOVR in panel data form
#' data("NCOVR")
#' Year <- as.numeric(kronecker(c(1960,1970,1980,1990),matrix(1,nrow = dim(NCOVR)[1])))
#' HR <- c(NCOVR$HR60,NCOVR$HR70,NCOVR$HR80,NCOVR$HR90)
#' PS <- c(NCOVR$PS60,NCOVR$PS70,NCOVR$PS80,NCOVR$PS90)
#' UE <- c(NCOVR$UE60,NCOVR$UE70,NCOVR$UE80,NCOVR$UE90)
#' NCOVRpanel <- as.data.frame(cbind(Year,HR,PS,UE))
#' Tformula <- HR ~ PS + UE
#' LM_time <- lmtestspsur(Form = Tformula, data = NCOVRpanel, time = Year, W = W)
#'
#' @export
lmtestspsur <- function(Form = NULL, data = NULL, W = NULL,
                        X = NULL, Y = NULL, time = NULL,
                        G = NULL, N = NULL, Tm = NULL,
                        print_table = TRUE) {
  if (is.null(W)) stop("W matrix is needed")
  #check for row-standardization of W
  if (!is.null(W)){
    if (class(W) != "matrix") W <- as.matrix(W)
    rsumW <- rowSums(W)
    rsumW[rsumW == 0] <- 1
    nW <- dim(W)[1]
    W <- W / matrix(rep(rsumW, each = nW),
                    nrow = nW, ncol = nW, byrow = TRUE)
    W <- Matrix::Matrix(W)
  }

  if (is.null(time)) {  # G > 1 (no temporal correlations are modelled)
    if(!is.null(Form) && !is.null(data)){
      # Lectura datos
      if (!any(class(Form) == "Formula")) Form <- Formula::Formula(Form)
      get_XY <- get_data_spsur(formula=Form,data=data,W=W)
      Y <- get_XY$Y
      X <- get_XY$X
      G <- get_XY$G
      N <- get_XY$N
      Tm <- get_XY$Tm
      p <- get_XY$p
      rm(get_XY)
    } else { ## Entry in matrix form.
      if(G == 1 && Tm > 1){  ## If G == 1 and Tm > 1 change values...
        G <- Tm
        Tm <- 1
      }
    }
  } else { #G = 1 and Tm > 1 (temporal correlations are modelled)
    if (class(time) != "factor") time <- as.factor(time)
    time <- droplevels(time)
    if (length(time) != nrow(data)) stop("time must have same length than the
                                         number of rows in data")
    mt <- terms(Form)
    G <- length(levels(time))
    Ylist <- vector("list",G)
    Xlist <- vector("list",G)
    p <- NULL
    namesX <- NULL
    levels_time <- levels(time)
    for (i in 1:G) {
      data_i <- model.frame(mt,data=data[time==levels_time[i],])
      Ylist[[i]] <- data_i[,1]
      Xlist[[i]] <- model.matrix(mt,data=data[time==levels_time[i],])
      p <- c(p,ncol(Xlist[[i]]))
      namesX <- c(namesX,paste(colnames(Xlist[[i]]),i,sep="_"))
    }
    Y <- matrix(unlist(Ylist),ncol=1)
    X <- as.matrix(Matrix::bdiag(Xlist))
    colnames(X) <- namesX
    N <- length(Ylist[[1]]); Tm <- 1
  }
  res <- sur3_spdiag(Tm=Tm,G=G,N=N,Y=Y,X=X,W=W)
  # cat(" \n\n")
  table_results <- cbind( res$stat, res$df,
                          pchisq(res$stat,res$df,lower.tail = FALSE) )
  rownames(table_results) <- res$stat_names
  colnames(table_results) <- c("LM-Stat.", "DF", "p-value")
  if(print_table)  printCoefmat(table_results, P.values = TRUE,
                                has.Pvalue = TRUE, digits = 4)
  return(res)
}
