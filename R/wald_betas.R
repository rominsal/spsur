#' @name wald_betas
#' @rdname wald_betas
#'
#' @title Wald tests on the \emph{beta} coefficients of the equation of 
#'  the SUR model
#'
#'
#' @description
#' The function \code{\link{wald_betas}} can be seen as a complement 
#' to the restricted estimation procedures included in the functions 
#' \code{\link{spsurml}} and \code{\link{spsur3sls}}. 
#' \code{\link{wald_betas}} obtains Wald tests for sets of linear 
#' restrictions on the coefficients \eqn{\beta} of the SUR model.
#' The restrictions may involve coefficients of the same equation or 
#' coefficients from different equations. The function has great flexibility 
#' in this respect. Note that \code{\link{wald_betas}} is more general than 
#' \code{\link{lr_betas}} in the sense that the last function
#' only allows to test for restrictions of homogeneity of subsets of 
#' \eqn{\beta} coefficients among the different equations in the SUR model, 
#' and in a maximum-likelihood framework.
#'
#'  In order to work with \code{\link{wald_betas}}, the model on which the 
#'  linear restrictions are to be tested needs to exists as an \emph{spsur} 
#'  object.  Using the information contained in the object, 
#'  \code{\link{wald_betas}} obtains the corresponding Wald estatistic 
#'  for the null hypotheses specified by the user through the \emph{R} row 
#'  vector and \emph{b} column vector, used also in \code{\link{spsurml}} 
#'  and \code{\link{spsur3sls}}. The function shows the value of the Wald test
#'  statistics and its associated p-values.
#'
#' @usage wald_betas (obj , R , b)
#'
#' @param obj An \code{spsur} object created by \code{\link{spsurml}},
#'            \code{\link{spsur3sls}} or \code{\link{spsurtime}}.
#' @param R  A row vector of order \eqn{(1xPr)} showing  the set
#'  of \emph{r} linear constraints on the \eqn{\beta} parameters. 
#'  The \emph{first} restriction appears in the first \emph{K} terms 
#'  in \emph{R}, the \emph{second} restriction in the next \emph{K} terms 
#'  and so on. 
#' @param b A column vector of order \emph{(rx1)} with the values of the 
#'   linear restrictions on the \eqn{\beta} parameters.
#'
#' @return Object of \code{htest} class including the Wald
#'   statistic, the corresponding p-value, the degrees of
#'   freedom and the values of the sample estimates.
#'
#' @author
#'   \tabular{ll}{
#'   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#'   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#' @references
#'   \itemize{
#'     \item Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model
#'        selection strategies in a SUR framework. The case of regional
#'        productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
#'        <doi:10.1007/s00168-014-0624-2>
#'        
#'     \item Mur, J., Lopez, F., and Herrera, M. (2010). Testing for spatial
#'       effects in seemingly unrelated regressions. \emph{Spatial Economic 
#'       Analysis}, 5(4), 399-440. 
#'       <doi:10.1080/17421772.2010.516443>
#'       
#'      \item Anselin, L. (2016) Estimation and Testing in the Spatial Seemingly 
#'       Unrelated Regression (SUR). \emph{Geoda Center for Geospatial Analysis 
#'       and Computation, Arizona State University}. Working Paper 2016-01.
#'       <doi:10.13140/RG.2.2.15925.40163>
#'   }
#'
#'
#' @seealso
#'  \code{\link{spsurml}}, \code{\link{spsur3sls}}, \code{\link{lr_betas}}
#'
#' @examples
#' ## VIP: The output of the whole set of the examples can be examined 
#' ## by executing demo(demo_wald_betas, package="spsur")
#' 
#' #################################################
#' ######## CROSS SECTION DATA (G=1; Tm>1) ########
#' #################################################
#'
#' ##### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' rm(list = ls()) # Clean memory
#' data(spc)
#' lwspc <- spdep::mat2listw(Wspc, style = "W")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ### Estimate SUR-SLM model
#' spcsur.slm <- spsurml(formula = Tformula, data = spc, 
#'                       type = "slm", listw = lwspc)
#' summary(spcsur.slm)
#' ### H_0: equality between SMSA coefficients in both equations.
#' R1 <- matrix(c(0,0,0,1,0,0,0,-1), nrow=1)
#' b1 <- matrix(0, ncol=1)
#' 
#' wald_betas(spcsur.slm, R = R1, b = b1)
#' 
#' \donttest{
#' ### Estimate restricted SUR-SLM model
#' spcsur.slmr <- spsurml(formula = Tformula, data = spc, 
#'                       type = "slm", listw = lwspc,
#'                       R = R1, b = b1)
#' summary(spcsur.slmr)
#' 
#' ### H_0: equality between intercepts and SMSA coefficients in both equations.
#' R2 <- matrix(c(1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,-1),
#'              nrow = 2, ncol = 8, byrow = TRUE)
#' b2 <- matrix(c(0,0),ncol=1)
#' wald_betas(spcsur.slm, R = R2, b = b2)
#' ### Estimate restricted SUR-SLM model
#' spcsur.slmr2 <- spsurml(formula = Tformula, data = spc, 
#'                       type = "slm", listw = lwspc,
#'                       R = R2, b = b2)
#' 
#' #####################################
#' #########  G=1; Tm>1         ########
#' #####################################
#'
#' ##### Example 2: Homicides + Socio-Economics (1960-90)
#' #
#' rm(list = ls()) # Clean memory
#' ### Read NCOVR.sf object
#' data(NCOVR, package = "spsur")
#' nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
#' ### Some regions with no links...
#' lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' ##################################
#' ### A SUR-SLM model
#' NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                        type = "slm", listw = lwncovr,
#'                        method = "Matrix", zero.policy = TRUE, 
#'                        control = list(fdHess = TRUE))
#' summary(NCOVRSUR.slm)
#' R1 <- matrix(c(0,1,0,0,-1,0), nrow=1)
#' b1 <- matrix(0, ncol=1)
#' wald_betas(NCOVRSUR.slm, R = R1, b = b1)
#' NCOVRSUR.slmr <-spsurml(formula = Tformula, data = NCOVR.sf, 
#'                        type = "slm", listw = lwncovr,
#'                        method = "Matrix", zero.policy = TRUE, 
#'                        control = list(fdHess = TRUE),
#'                        R = R1, b = b1)
#' summary(NCOVRSUR.slmr)    
#' }
#' @export
 wald_betas <- function(obj , R , b){
  z <- obj 
  betas <- Matrix::Matrix(matrix(z$coefficients, ncol = 1))
  rownames(betas) <- names(z$coefficients)
  cov_betas <- Matrix::Matrix(z$resvar[rownames(betas),
                                    rownames(betas)])
  R <- Matrix::Matrix(R)
  colnames(R) <- rownames(betas)
  b <- Matrix::Matrix(matrix(b, ncol=1))
  holg <- (R %*% betas) - b
  parameter <- nrow(as.matrix(R))
  attr(parameter, "names") <- "df"
  statistic <- as.numeric( Matrix::t(holg) %*%
        Matrix::solve(R %*% cov_betas %*% Matrix::t(R),holg) )
  attr(statistic, "names") <- "Wald test"
  method <- paste("Wald test on beta parameters")
  p.value <- pchisq(statistic, df = parameter, 
                    lower.tail = FALSE)
  estimate <- as.numeric(betas)
  names(estimate) <- rownames(betas)
  data.name <- z$call[[3]]
  res <- list(statistic = statistic, parameter = parameter, 
              p.value = p.value, # estimate = estimate, 
              method = method, data.name = data.name)
  class(res) <- "htest"
  res  
 }
