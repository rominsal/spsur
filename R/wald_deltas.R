#' @name wald_deltas
#' @rdname wald_deltas
#'
#' @title Wald tests for spatial parameters coefficients.
#'
#' @description
#' Function \code{\link{wald_deltas}} obtains Wald tests for linear restrictions
#'  on the spatial coefficients of a SUR model that has been estimated previously through the function
#'  \code{\link{spsurml}}. The restrictions can affect to coefficients of the same equation
#'  (i.e., \eqn{\lambda_{g}=\rho_{g} forall g}) or can involve coefficients from different equations
#'  (i.e., \eqn{\lambda_{g}=\lambda_{h}}). The function has great flexibility in this respect.
#'  Note that \code{\link{wald_deltas}} only works in a maximum-likelihood framework.
#'
#'  In order to work with \code{\link{wald_betas}}, the model on which the linear restrictions are
#'  to be tested needs to exists as an \emph{spsur} object. Using the information contained in the object,
#'  \code{\link{wald_deltas}} obtains the corresponding Wald estatistic for the null hypotheses
#'  specified by the user through the \emph{R} row vector and \emph{b} column vector discussed, used also
#'  in \code{\link{spsurml}}. The function shows the  resulting Wald test statistics
#'  and their corresponding p-values.
#'
#'
#' @param    results : An object created with \code{\link{spsurml}} or \code{\link{spsur3sls}}. This argument
#'  serves the user to indicate the spatial SUR model, previously estimated by maximum-likelihood or 3sls,
#'  where the set of linear restrictions are to be tested.
#' @param    R       : A row vector of order \emph{(1xGr)} or \emph{(1x2Gr)} showing  the set of \emph{r} linear
#'  constraints on the  spatial parameters. The last case is reserved to \strong{"sarar"} models where there appear
#'  \emph{G} parameters \eqn{\lambda_{g}} and \emph{G} parameters \eqn{\rho_{g}}, \emph{2G} spatial parameters
#'  in total. The \emph{first} restriction appears in the first \emph{G} terms in \emph{R} (\emph{2G} for the
#'  \strong{"sarar"} case),   the second restriction in the next \emph{G} terms (\emph{2G} for the
#'  \strong{"sarar"} case) and so on. Default = NULL.
#' @param    b       : A column vector of order \emph{(rx1)} with the values of the linear restrictions on the
#'  \eqn{\beta} parameters. Default = NULL.
#'
#' @return
#' The output of the function is very simple and consists of two pieces of information, the value of the Wald
#' statistic and the corresponding p-value, plus the degrees of freedom of the test.
#'
#'   \tabular{ll}{
#'   \code{Wald stat} \tab The value of Wald test. \cr
#'   \code{p_val}    \tab The p-value of Wald test. \cr
#'   \code{q}    \tab Degrees of freedom of the corresponding \eqn{\chi^{2}} distribution. \cr
#'   }
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
#'   }
#'
#' @seealso
#'  \code{\link{spsurml}}, \code{\link{spsur3sls}}
#' @examples
#'
#' #################################################
#' ######## CROSS SECTION DATA (G>1; Tm=1) ########
#' #################################################
#' rm(list = ls()) # Clean memory
#' data(spc)
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' #################################
#' ## Estimate SUR-SLM model
#' spcsur.slm <-spsurml(Form = Tformula, data = spc, type = "slm", W= Wspc)
#' summary(spcsur.slm)
#' ## H_0: equality of the lambda parameters of both equations.
#' R1 <- matrix(c(1,-1), nrow=1)
#' b1 <- matrix(0, ncol=1)
#' wald_deltas(results = spcsur.slm, R = R1, b = b1)
#'
#' #################################
#' ## Estimate SUR-SEM model
#' spcsur.sem <-spsurml(Form = Tformula, data = spc, type = "sem", W = Wspc)
#' summary(spcsur.sem)
#' ## H_0: equality of the rho parameters of both equations.
#' R2 <- matrix(c(1,-1), nrow=1)
#' b2 <- matrix(0, ncol=1)
#' wald_deltas(results = spcsur.sem, R = R2, b = b2)
#'
#' #################################
#' ## Estimate SUR-SARAR model
#' spcsur.sarar <-spsurml(Form = Tformula, data = spc,
#'                        type = "sarar", W = Wspc)
#' summary(spcsur.sarar)
#' ## H_0: equality of the lambda and rho parameters of both equations.
#' R3 <- matrix(c(1,-1,0,0,0,0,1,-1),nrow=2,ncol=4,byrow=TRUE)
#' b3 <- matrix(c(0,0), ncol=1)
#' wald_deltas(results = spcsur.sarar, R = R3, b = b3)
#'
#' ####################################
#' ########  G=1; Tm>1         ########
#' ####################################
#'
#' #### Example 2: Homicides + Socio-Economics (1960-90)
#' \dontrun{
#' # It could make an error out-of-memory in some computers
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#' #################################
#' ## A SUR-SLM model
#' NCOVRSUR.slm <-spsurml(Form = Tformula, data = NCOVR, type = "slm", W = W)
#' summary(NCOVRSUR.slm)
#' ## H_0: equality of the lambda parameters of both equations.
#' R1 <- matrix(c(1,-1), nrow=1)
#' b1 <- matrix(0, ncol=1)
#' wald_deltas(results = NCOVRSUR.slm, R = R1, b = b1)
#' #################################
#' ## Estimate SUR-SEM model
#' NCOVRSUR.sem <-spsurml(Form = Tformula, data = NCOVR, type = "sem", W = W)
#' summary(NCOVRSUR.sem)
#' ## H_0: equality of the rho parameters of both equations.
#' R2 <- matrix(c(1,-1), nrow=1)
#' b2 <- matrix(0, ncol=1)
#' wald_deltas(results = NCOVRSUR.sem, R = R2, b = b2)
#' }
#' @export
wald_deltas <- function(results , R , b){
  z <- results # OBJETO QUE INCLUYE ESTIMACIÓN EN Rbetas <- z$betas
  deltas <- Matrix::Matrix(matrix(z$deltas,ncol=1))
  rownames(deltas) <- names(z$deltas)
  cov_deltas <- Matrix::Matrix(z$cov[rownames(deltas),rownames(deltas)])
  R <- Matrix::Matrix(R)
  colnames(R) <- rownames(deltas)
  b <- Matrix::Matrix(matrix(b,ncol=1))
  holg <- R %*% deltas - b
  q <- nrow(as.matrix(R))
  Wald <- as.numeric(Matrix::t(holg) %*%
                    Matrix::solve(R %*% cov_deltas %*% Matrix::t(R),holg)  )
  p_val <- pchisq(Wald,df=q,lower.tail=FALSE)
  #cat("\n statistical discrepancies: "); print(as.matrix(holg))
  cat("\n Wald stat.: ",round(Wald,3)," (",round(p_val,3),")",sep="")
  res <- list(stat = Wald,
              p_val = p_val,
              q = q,
              R = as.matrix(R),
              b = as.matrix(b),
              discr = as.matrix(holg) )
}
