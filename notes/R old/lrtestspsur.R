#' @name lrtestspsur
#' @rdname lrtestspsur
#'
#' @title Likelihood Ratio tests for the specification of spatial SUR models.
#'
#' @description The function computes a set of Likelihood Ratio tests, LR, that help
#'  the user to select the spatial structure of the SUR model. To achieve this goal, \code{\link{lrtestspsur}}
#'  needs to estimate the SUR models \strong{"sim"}, \strong{"slm"}, \strong{"sem"}, \strong{"sdm"},
#'  and \strong{"sarar"}, using the function \code{\link{spsurml}}.
#'
#'  The five models listed above are related by a nesting sequence, so they can be compared using
#'  the adequate LR tests. The function shows the log-likelihood corresponding to the maximum-likelihood
#'  estimates and the sequence of LR tests.
#'
#' @inheritParams spsurml
#' @param time Time variable.
#'
#' @details  A fundamental result in maximum-likelihood estimation shows that if \emph{model A} is nested
#' in \emph{model B}, by a set of \emph{n} restrictions on the parameters of \emph{model B}, then,
#' as the sample size increases, the test statistic: \emph{\eqn{-2log[l(H_{0}) / l(H_{A})]}}
#' is a \eqn{\chi^{2}(n)}, being l(H_{0} the estimated likelihood under the null hypothesis
#' (\emph{model A}) and  l(H_{A} the estimated likelihood under the alternative hypothesis (\emph{model B}).
#'
#'  The list of (spatial) models that can be estimated with the function \code{\link{spsurml}} includes
#'   the following (in addition to the \strong{"slx"} and \strong{"sdem"}):
#'
#'  \itemize{
#'     \item \strong{"sim"}: SUR model with no spatial effects
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"slm"}: SUR model with spatial lags of the explained variables
#'       \deqn{y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
#'     \item \strong{"sem"}: SUR model with spatial errors
#'       \deqn{ y_{tg} = X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \rho_{g} Wu_{tg} + \epsilon_{tg} }
#'     \item \strong{"sdm"}: SUR model of the Spatial Durbin type
#'       \deqn{ y_{tg} = \lambda_{g} Wy_{tg} + X_{tt} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
#'     \item \strong{"sarar"}: SUR model with spatial lags of the explained variables and spatial
#'       errors
#'       \deqn{ y_{tg} = \lambda_{g} Wy_{tg} + X_{tg} \beta_{g} + u_{tg} }
#'       \deqn{ u_{tg} = \rho_{g} W u_{tg} + \epsilon_{tg} }
#'   }
#'
#'
#'   This collection of models can be compared, on objective bases, using the LR principle  and the
#'    following  nesting relations:
#'
#'   \itemize{
#'     \item  \strong{"sim"} vs \strong{"sem"}, where the null hypotheses, in the \strong{"sem"} equation, are:
#'
#'   \deqn{ H_{0}: \rho_{g}=0 forall g vs  H_{A}: \rho_{g} ne 0 exist g}
#'
#'     \item  \strong{"sim"} vs \strong{"slm"}, where the null hypotheses, in the \strong{"slm"} equation, are:
#'
#'   \deqn{ H_{0}: \lambda_{g}=0 forall g vs  H_{A}: \lambda_{g} ne 0 exist g}
#'
#'     \item  \strong{"sim"} vs \strong{"sarar"}, where the null hypotheses, in the \strong{"sarar"} equation, are:
#'
#'   \deqn{ H_{0}: \rho_{g}=\lambda_{g}=0 forall g vs  H_{A}: \rho_{g} ne 0 or \lambda_{g} ne 0 exist g}
#'
#'     \item  \strong{"sem"} vs \strong{"sarar"}, where the null hypotheses, in the \strong{"sarar"} equation, are:
#'
#'   \deqn{ H_{0}: \lambda_{g}=0 forall g vs  H_{A}: \lambda_{g} ne 0 exist g}
#'
#'     \item  \strong{"slm"} vs \strong{"sarar"}, where the null hypotheses, in the \strong{"sarar"} equation, are:
#'
#'   \deqn{ H_{0}: \rho_{g}=0 forall g vs  H_{A}: \rho_{g} ne 0 exist g}
#'
#'     \item  \strong{"sem"} vs \strong{"sdm"}, also known as \emph{LR-COMFAC}, where the null hypotheses, in the \strong{"sdm"}
#'      equation, are:
#'
#'    \deqn{ H_{0}: -\lambda_{g}\beta_{g}=\theta_{g} forall g vs  H_{A}: -\lambda_{g}\beta_{g} ne \theta_{g} exist g}
#'
#'   }
#'
#'  The degrees of freedom of the corresponding \eqn{\chi^{2}} distribution is \emph{G} in the cases of \strong{"sim"}
#'  vs \strong{"sem"}, \strong{"sim"} vs \strong{"slm"}, \strong{"sem"} vs \strong{"sarar"}, \strong{"slm"} vs
#'  \strong{"sarar"}  and \strong{"sem"} vs  \strong{"sdm"} and \emph{2G} in the case of \strong{"sim"} vs \strong{"sarar"}.
#'   Moreover, function \code{\link{lrtestspsur}} also returns the p-values associated to the corresponding LR.
#'
#'
#' @return
#'    \code{\link{lrtestspsur}}, first, prints the value of the estimated log-likelihood for
#'    the major spatial specifications. Then, the function shows the values of the LR statistics corresponding to the nested
#'    and nesting models compared, together with their associated p-value.
#'
#'
#' @author
#'   \tabular{ll}{
#'   Fernando López  \tab \email{fernando.lopez@@upct.es} \cr
#'   Román Mínguez  \tab \email{roman.minguez@@uclm.es} \cr
#'   Jesús Mur  \tab \email{jmur@@unizar.es} \cr
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
#'   }
#'
#' @seealso
#' \code{\link{spsurml}}, \code{\link{lmtestspsur}}
#'
#' @examples
#' #################################################
#' ######## CROSS SECTION DATA (nG=1; nT>1) ########
#' #################################################
#'
#' #### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
#' \donttest{
#' ## It usually requires 1-2 minutes maximum
#' rm(list = ls()) # Clean memory
#' data("spc")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' LRs <- lrtestspsur(Form = Tformula, data = spc, W = Wspc)
#' }
#'
#' #################################################
#' ######## CROSS SECTION DATA (nG>1; nT=1) ########
#' #################################################
#'
#' #### Example 2: Homicides & Socio-Economics (1960-90)
#  # Different number of exogenous variables in each equation
#' # Homicides and selected socio-economic characteristics for
#' # continental U.S. counties.
#' # Data for four decennial census years: 1960, 1970, 1980 and 1990.
#' # https://geodacenter.github.io/data-and-lab/ncovr/
#' \donttest{
#' ## It could require some minutes
#' rm(list = ls()) # Clean memory
#' data("NCOVR")
#' Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 | PS80 + UE80 + RD80 |
#'             PS90 + UE90 + RD90 + PO90
#' LRs <- lrtestspsur(Form = Tformula, data = NCOVR, W = W)
#' }
#'
#' ################################################################
#' ######## PANEL DATA: TEMPORAL CORRELATIONS (nG=1; nT>1) ########
#' ################################################################
#'
#' #### Example 3: Classic panel data
#' \donttest{
#' ## It could require some minutes
#' rm(list = ls()) # Clean memory
#' data(NCOVR)
#' N <- nrow(NCOVR)
#' Tm <- 4
#' index_time <- rep(1:Tm, each = N)
#' index_indiv <- rep(1:N, Tm)
#' pHR <- c(NCOVR$HR60, NCOVR$HR70, NCOVR$HR80, NCOVR$HR90)
#' pPS <- c(NCOVR$PS60, NCOVR$PS70, NCOVR$PS80, NCOVR$PS90)
#' pUE <- c(NCOVR$UE60, NCOVR$UE70, NCOVR$UE80, NCOVR$UE90)
#' pNCOVR <- data.frame(indiv = index_indiv, time = index_time, HR = pHR, PS = pPS, UE = pUE)
#' rm(NCOVR,pHR,pPS,pUE,index_time,index_indiv)
#' form_pHR <- HR ~ PS + UE
#' LRs <- lrtestspsur(Form = form_pHR, data = pNCOVR, W = W, time = pNCOVR$time)
#' }
#' @export

lrtestspsur <- function(Form = NULL, data = NULL, W = NULL,
                        X = NULL, Y = NULL, time = NULL,
                        G = NULL, N = NULL, Tm = NULL) {
  ## PURPOSE:
  # Realiza todos los test LR.
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
  W <- as(W,"dgCMatrix")
  if (is.null(time)){  # G > 1 (no temporal correlations are modelled)
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

  ## Lik modelo SIM
  model_sim <- fit_spsursim(Tm=Tm,G=G,N=N,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sim <- model_sim$llsur
  cat("LogLik SUR-SIM:  ", round(llsur_sim,3),"\n")
  ## Lik modelo SLM
  model_slm <- fit_spsurslm(Tm=Tm,G=G,N=N,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_slm <- model_slm$llsur
  cat("LogLik SUR-SLM:  ", round(llsur_slm,3),"\n")
  ## Lik modelo SEM
  model_sem <- fit_spsursem(Tm=Tm,G=G,N=N,
                            Y=Y,X=X,W=W,trace=FALSE)
  llsur_sem <- model_sem$llsur
  cat("LogLik SUR-SEM:  ", round(llsur_sem,3),"\n")
  ## Lik modelo SARAR
  model_sarar <- fit_spsursarar(Tm=Tm,G=G,N=N,
                                Y=Y,X=X,W=W,trace=FALSE)
  llsur_sarar <- model_sarar$llsur
  cat("LogLik SUR-SARAR: ", round(llsur_sarar,3),"\n")

  if(!is.null(Form) && !is.null(data)){
    if (!is.null(time)){
      model_sdm <- spsurtime(Form = Form, data = data, time = time, type = "sdm", W = W,
                             cov = FALSE, trace = FALSE)
    } else {
      model_sdm <- spsurml(Form = Form, data = data, type = "sdm", W = W,
                           cov = FALSE, control = list(tol = 0.1, maxit = 200, trace = FALSE))
    }
    llsur_sdm <- model_sdm$llsur
    cat("LogLik SUR-SDM:  ", round(llsur_sdm,3),"\n" )
  }

  lr_sim_slm <- -2*(llsur_sim - llsur_slm)
  pval_lr_sim_slm <- pchisq(lr_sim_slm,df=G,lower.tail=FALSE)

  lr_sim_sem <- -2*(llsur_sim - llsur_sem)
  pval_lr_sim_sem <- pchisq(lr_sim_sem,df=G,lower.tail=FALSE)

  lr_sim_sarar <- -2*(llsur_sim - llsur_sarar)
  pval_lr_sim_sarar <- pchisq(lr_sim_sarar,df=2*G,lower.tail=FALSE)

  lr_slm_sarar <- -2*(llsur_slm - llsur_sarar)
  pval_lr_slm_sarar <- pchisq(lr_slm_sarar,df=G,lower.tail=FALSE)

  lr_sem_sarar <- -2*(llsur_sem - llsur_sarar)
  pval_lr_sem_sarar <- pchisq(lr_sem_sarar,df=G,lower.tail=FALSE)

  if(!is.null(Form) && !is.null(data)){
    lr_slm_sdm <- -2*(llsur_slm - llsur_sdm)
    pval_lr_slm_sdm <- pchisq(lr_slm_sdm,df=G,lower.tail=FALSE)

    lr_sem_sdm <- -2*(llsur_sem - llsur_sdm)
    pval_lr_sem_sdm <- pchisq(lr_sem_sdm,df=G,lower.tail=FALSE)
  }


  cat("\n LR test SUR-SIM versus SUR-SLM: \n")
  cat("statistic: ",round(lr_sim_slm,3),
      " p-value: (",round(pval_lr_sim_slm,3),") \n")
  cat(" LR test SUR-SIM versus SUR-SEM: \n")
  cat("statistic: ",round(lr_sim_sem,3),
      " p-value: (",round(pval_lr_sim_sem,3),") \n")
  cat(" LR test SUR-SIM versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_sim_sarar,3),
      " p-value: (",round(pval_lr_sim_sarar,3),") \n")
  cat(" LR test SUR-SLM versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_slm_sarar,3),
      " p-value: (",round(pval_lr_slm_sarar,3),") \n")
  cat(" LR test SUR-SEM versus SUR-SARAR: \n")
  cat("statistic: ",round(lr_sem_sarar,3),
      " p-value: (",round(pval_lr_sem_sarar,3),") \n")

  if(!is.null(Form) && !is.null(data)){
    cat("\n")
    cat(" LR testSUR-SLM versus SUR-SDM: \n")
    cat("statistic: ",round(lr_slm_sdm,3),
        " p-value: (",round(pval_lr_slm_sdm,3),") \n")
    cat(" COMFAC LR test SUR-SEM versus SUR-SDM: \n")
    cat("statistic: ",round(lr_sem_sdm,3),
        " p-value: (",round(pval_lr_sem_sdm,3),") \n")
  }


  res <- list(lr_sim_slm = lr_sim_slm,
              pval_lr_sim_slm=pval_lr_sim_slm,
              lr_sim_sem = lr_sim_sem,
              pval_lr_sim_sem = pval_lr_sim_sem,
              lr_sim_sarar  = lr_sim_sarar,
              pval_lr_sim_sarar = pval_lr_sim_sarar,
              lr_slm_sarar = lr_slm_sarar,
              pval_lr_slm_sarar = pval_lr_slm_sarar,
              lr_sem_sarar = lr_sem_sarar,
              pval_lr_sem_sarar = pval_lr_sem_sarar
              )
 res
}



