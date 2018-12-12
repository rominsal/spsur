#' @name lmtestspsur
#' @rdname lmtestspsur
#'
#' @title Testing for spatial effects in Seemingly Unrelated Regression
#'
#' @description The function \code{\link{spsurml}} reports five statistis for
#'   testing spatial dependence in Semmingly Unrelated Regression (SUR) models.
#'   The five statistics are based on Lagrangian Multipliers principe.
#'   These tests can be used with one cross sectional data set or in a panel
#'   data framework.\cr
#'   The statistics are similar to the LM test of spatial autocorrelation for
#'   a single equation (see \code{\link[spdep]{lm.LMtests}} in package
#'   \pkg{spdep}). The LM test for error dependence (LM-SUR-err),
#'   the LM test for a missing spatially lagged dependent variable (LM-SUR-lag)
#'   and LM-SUR-SARAR for simultaneous spatial lagged dependent variable and
#'   spatial error structure. Variants of these tests robust to the presence
#'   of the other (RLM-SUR-err, RLM-SUR-lag) are also implemented.
#'
#' @inheritParams spsurml
#' @param print_table Logical value to print the output. Default = \code{TRUE}
#' @param time Variable including temporal periods.
#'
#' @details The model specification in the general case of \emph{nR}
#'   individuals, \emph{nT} cross sections and \emph{nG} equations, with
#'   spatial interaction mechanisms, is as follows: \cr
#'     \deqn{y_{gt} = \lambda_{g} W y_{gt} + X_{gt} \beta_{g} + u_{gt}}
#'     \deqn{u_{gt} = \rho_{g} W u_{gt} + \epsilon_{gt}}
#'   where \eqn{y_{gt}}, \eqn{u_{gt}} and \eqn{\epsilon_{gt}} are (\emph{nR}x1)
#'   vectors; \eqn{X_{gt}} is a matrix of exogenous variables of order
#'   (\emph{nR}x\eqn{p_g}); \eqn{\lambda_g} and \eqn{\rho_g} are parametres of
#'   spatial dependence; \eqn{W} is the \emph{nR}x\emph{nR} spatial weight
#'   matrix.
#'
#'   Five tests of diagnotic of spatial autocorrelation based on the principle
#'   of Lagrange Multiplier are developed in this code:
#'   \itemize{
#'     \item \strong{LM-SUR-lag}: \eqn{H_{0}: \lambda=0 (\forall g)} versus
#'       \eqn{H_{A}: No H_{0}}
#'      \item \strong{LM-SUR-err}: \eqn{H_{0}: \rho=0 (\forall g)} versus
#'        \eqn{H_{A}: No H_{0}}
#'      \item \strong{LM-SUR-sarar}: \eqn{H_{0}: \lambda=0} and
#'        \eqn{\rho=0 (\forall g)} versus \eqn{H_{A}: No H_{0}}
#'   }
#'   As the LM-SUR-lag and LM-SUR-err tests can be strongly oversized, two
#'   robust alternative tests are included under local misspecifications of
#'   the alternative: \strong{LM*-SUR-SAR} and \strong{LM*-SUR-SEM}.
#'
#' @return A list including:
#'   \tabular{ll}{
#'   \code{stat_names}   \tab the name of Lagrange Multiplier tests \cr
#'   \code{stat}   \tab the value of five Lagrange Multiplier tests \cr
#'   \code{df} \tab number of degrees of freedom for each test \cr
#'   }
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
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es}\cr
#'  }
#'
#' @seealso
#' \code{\link{spsurml}},  \code{\link{spsur3sls}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (nG>1; nT=1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' data("spc")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' LMs <- lmtestspsur(Form = Tformula, data = spc, W = Wspc)
#'
#' #################################################
#' ######## PANEL DATA (nG>1; nT>1)         ########
#' #################################################
#'
#' #### Example 2: Homicides & Socio-Economics (1960-90)
#' # Homicides and selected socio-economic characteristics for
#' # continental U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' data("NAT")
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 |PS80 + UE80 | PS90 + UE90
#' LMs <- lmtestspsur(Form = Tformula, data = NAT, W = W)
#'
#' # With different number of exogenous variables in each equation
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 | PS80 + UE80 +RD80 |
#'             PS90 + UE90 + RD90 + PO90
#' LMs <- lmtestspsur(Form = Tformula, data = NAT, W = W)
#'
#' ################################################################
#' ######## PANEL DATA: TEMPORAL CORRELATIONS (nG=1; nT>1) ########
#' ################################################################
#'
#' #### Example 3: Unemployment data in Italy
#' data("unemp_it")
#' form_un <- unrate  ~ empgrowth + partrate + agri + cons + serv
#' LM_time <- lmtestspsur(Form = form_un, data = unemp_it,
#'                        time = unemp_it$year, W = W_italy)
#'
#' @export
lmtestspsur <- function(Form = NULL, data = NULL, W = NULL,
                        X = NULL, Y = NULL, time = NULL,
                        nG = NULL, nR = NULL, nT = NULL,
                        print_table = TRUE) {
  if (is.null(W)) stop("W matrix is needed") else W <- Matrix::Matrix(W)
  if (is.null(time)) {  # nG > 1 (no temporal correlations are modelled)
    if(!is.null(Form) && !is.null(data)){
      # Lectura datos
      if (!class(Form) == "Formula") Form <- Formula::Formula(Form)
      get_XY <- get_data_spsur(formula=Form,data=data,W=W)
      Y <- get_XY$Y
      X <- get_XY$X
      nG <- get_XY$nG
      nR <- get_XY$nR
      nT <- get_XY$nT
      p <- get_XY$p
      rm(get_XY)
    }
  } else { #nG = 1 and nT > 1 (temporal correlations are modelled)
    if (class(time) != "factor") time <- as.factor(time)
    time <- droplevels(time)
    if (length(time) != nrow(data)) stop("time must have same length than the
                                         number of rows in data")
    mt <- terms(Form)
    nG <- length(levels(time))
    Ylist <- vector("list",nG)
    Xlist <- vector("list",nG)
    p <- NULL
    namesX <- NULL
    levels_time <- levels(time)
    for (i in 1:nG) {
      data_i <- model.frame(mt,data=data[time==levels_time[i],])
      Ylist[[i]] <- data_i[,1]
      Xlist[[i]] <- model.matrix(mt,data=data[time==levels_time[i],])
      p <- c(p,ncol(Xlist[[i]]))
      namesX <- c(namesX,paste(colnames(Xlist[[i]]),i,sep="_"))
    }
    Y <- matrix(unlist(Ylist),ncol=1)
    X <- as.matrix(Matrix::bdiag(Xlist))
    colnames(X) <- namesX
    nR <- length(Ylist[[1]]); nT <- 1
  }
  res <- sur2_spdiag(nT=nT,nG=nG,nR=nR,Y=Y,X=X,W=W)
  cat(" \n\n")
  table_results <- cbind( res$stat, res$df,
                          pchisq(res$stat,res$df,lower.tail = FALSE) )
  rownames(table_results) <- res$stat_names
  colnames(table_results) <- c("LM-Stat.", "DF", "p-value")
  if(print_table)  printCoefmat(table_results, P.values = TRUE,
                                has.Pvalue = TRUE, digits = 4)
  return(res)
}
