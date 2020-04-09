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
#' rm(list = ls()) # Clean memory
#' data("spc")
#' Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
#' ## It usually requires 1-2 minutes maximum
#' ## LRs <- lrtestspsur(Form = Tformula, data = spc, W = Wspc)
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

lrtestspsur <- function(objectr,  objectu = NULL, ...) {
                        #formula = NULL, data = NULL,
                        #listw = NULL, vtypes = NULL,
                        #method = "Matrix",
  # LR tests of model specification.
  class(objectr) <- "sarlm" ## ANOVA for sarlm class
  if(is.null(objectu)) {
    anova_table <- spatialreg::anova.sarlm(objectr)
    attr(anova_table, "row.names") <- paste(objectr$type,
                                    "model", sep = " ")
  } else {
    class(objectu) <- "sarlm"
    anova_table <- spatialreg::anova.sarlm(objectr, objectu)
    attr(anova_table, "row.names") <- c(paste(objectr$type,
                                        "model", sep = " "),
                                        paste(objectu$type,
                                        "model", sep = " "))
  }
  res <- anova_table
  res
}
  #### CODE TO FIX ########
  #### IT WOULD ALSO ALLOW TO COMPARE MODELS INCLUDING
  #### FORMULA, DATA, LISTW AND VTYPES (VECTOR OF TYPES)
  #### PROBLEMS WITH FORMULA AND CLASS SARLM
#   if(!is.null(objectr)) {
#     model <- objectr
#     class(model) <- "sarlm"
#     name_model <- paste("model", objectr$type, sep = "_")
#     assign(name_model, model)
#     vtypes <- c(vtypes, objectr$type)
#   }
#   if(!is.null(objectu)) {
#     model <- objectu
#     class(model) <- "sarlm"
#     name_model <- paste("model", objectu$type, sep = "_")
#     assign(name_model, model)
#     vtypes <- c(vtypes, objectu$type)
#   }
#   if (is.null(objectr) && is.null(objectu)) {
#     alltypes <- c("sim", "slx", "slm", "sdm", 
#                   "sem", "sdem", "sarar")
#     if (is.null(vtypes)) vtypes <- alltypes
#     if (is.null(formula) || is.null(data) || is.null(listw)) {
#       stop("Either restricted and unrestricted fitted models 
#             or formula, data and listw 
#             must be specified as arguments")
#     }
#     cl <- match.call()
#     name_form <- eval(quote(cl[[2]]))
#     form <- as.formula(eval(cl[[2]]))
#     assign(as.character(name_form), form, envir = .GlobalEnv)
#     for (i in seq_along(vtypes)) {
#       name_model <- paste("model", vtypes[i], sep = "_")
#       model <- spsurml(formula = formula, 
#                        data = data,
#                        listw = listw, type = vtypes[i], 
#                        method = method, 
#                        con = list(fdHess = TRUE))
#       class(model) <- "sarlm"
#       assign(name_model, model)
#     }
#   }
#   ### Print ANOVA models
#   if (length(vtypes) == 1) {
#     cat("\n ANOVA Table for ", vtypes[1], " model \n")
#     an_table <- spatialreg::anova.sarlm(get(name_model))
#     attr(an_table, "row.names") <- name_model
#     print(an_table)
#     res <- list()
#     model <- get(name_model)
#     class(model) <- "spsur"
#     name_an_table <- paste("an", vtypes[1], sep="_")
#     res[[name_an_table]] <- an_table
#     res[[name_model]] <- model
#     return(res)
#   }
#   if (exists("model_sim")) {
#     if (exists("model_slx")) {
#       cat("\n sim (restricted) versus slx (unrestricted) \n")
#       an_sim_slx <- spatialreg::anova.sarlm(model_sim, 
#                                             model_slx)
#       print(an_sim_slx)
#     }  
#     if (exists("model_slm")) {
#       cat("\n sim (restricted) versus slm (unrestricted) \n")
#       an_sim_slm <- spatialreg::anova.sarlm(model_sim, 
#                                             model_slm)
#       print(an_sim_slm)
#     }
#     if (exists("model_sdm")) {
#       cat("\n sim (restricted) versus sdm (unrestricted) \n")
#       an_sim_sdm <- spatialreg::anova.sarlm(model_sim, 
#                                             model_sdm)
#       print(an_sim_sdm)
#     }  
#     if (exists("model_sem")) {
#       cat("\n sim (restricted) versus sem (unrestricted) \n")
#       an_sim_sem <- spatialreg::anova.sarlm(model_sim, 
#                                             model_sem)
#       print(an_sim_sem)
#     }  
#     if (exists("model_sdem")) {
#       cat("\n sim (restricted) versus sdem (unrestricted) \n")
#       an_sim_sdem <- spatialreg::anova.sarlm(model_sim, 
#                                             model_sdem)
#       print(an_sim_sdem)
#     }  
#     if (exists("model_sarar")) {
#       cat("\n sim (restricted) versus sarar (unrestricted) \n")
#       an_sim_sarar <- spatialreg::anova.sarlm(model_sim, 
#                                              model_sarar)
#       print(an_sim_sarar)
#     }  
#   }  
#   if (exists("model_slx")) {
#     if (exists("model_sdm")) {
#       cat("\n slx (restricted) versus sdm (unrestricted) \n")
#       an_slx_sdm <- spatialreg::anova.sarlm(model_slx, 
#                                             model_sdm)
#       print(an_slx_sdm)
#     }
#     if (exists("model_sdem")) {
#       cat("\n slx (restricted) versus sdem (unrestricted) \n")
#       an_slx_sdem <- spatialreg::anova.sarlm(model_slx, 
#                                             model_sdem)
#       print(an_slx_sdem)
#     }
#   }  
#   
#   if (exists("model_slm")) {
#     if (exists("model_sdm")) {
#       cat("\n slm (restricted) versus sdm (unrestricted) \n")
#       an_slm_sdm <- spatialreg::anova.sarlm(model_slm, 
#                                             model_sdm)
#       print(an_slm_sdm)
#     }
#     if (exists("model_sarar")) {
#       cat("\n slm (restricted) versus sarar (unrestricted) \n")
#       an_slm_sarar <- spatialreg::anova.sarlm(model_slm, 
#                                             model_sarar)
#       print(an_slm_sarar)
#     }
#   }
#   
#   if (exists("model_sem")) {
#     if (exists("model_sdem")) {
#       cat("\n sem (restricted) versus sdem (unrestricted) \n")
#       an_sem_sdem <- spatialreg::anova.sarlm(model_sem, 
#                                             model_sdem)
#       print(an_sem_sdm)
#     }
#     if (exists("model_sarar")) {
#       cat("\n sem (restricted) versus sarar (unrestricted) \n")
#       an_sem_sarar <- spatialreg::anova.sarlm(model_sem, 
#                                               model_sarar)
#       print(an_sem_sarar)
#     }
#   }
#   lanova <- list(
#        if (exists("an_sim_slx")) anova_sim_slx = an_sim_slx,
#        if (exists("an_sim_slm")) anova_sim_slm = an_sim_slm,
#        if (exists("an_sim_sdm")) anova_sim_sdm = an_sim_sdm,
#        if (exists("an_sim_sem")) anova_sim_sem = an_sim_sem,
#        if (exists("an_sim_sdem")) anova_sim_sdem = an_sim_sdem,
#        if (exists("an_sim_sarar")) anova_sim_sarar = an_sim_sarar,
#        if (exists("an_slx_sdm")) anova_slx_sdm = an_slx_sdm,
#        if (exists("an_slx_sdem")) anova_slx_sdem = an_slx_sdem,
#        if (exists("an_slm_sdm")) anova_slm_sdm = an_slm_sdm,
#        if (exists("an_slm_sarar")) anova_slm_sarar = an_slm_sarar,
#        if (exists("an_sem_sdm")) anova_sem_sdm = an_sem_sdm,
#        if (exists("an_sem_sdem")) anova_sem_sdem = an_sem_sdem,
#        if (exists("an_sem_sarar")) anova_sem_sarar = an_sem_sarar)
#   lmodels <- list(
#      if (exists("model_sim")) model_sim = model_sim,
#      if (exists("model_slx")) model_slx = model_slx,
#      if (exists("model_slm")) model_slm = model_slm,
#      if (exists("model_sdm")) model_sdm = model_sdm,
#      if (exists("model_sem")) model_sem = model_sem,
#      if (exists("model_sdem")) model_sdem = model_sdem,
#      if (exists("model_sarar")) model_sarar = model_sarar)
#   for (i in seq_along(lmodels)) {
#     if (!is.null(lmodels[[i]])) class(lmodels[[i]]) <- "spsur"
#   }
#   res <- c(lanova, lmodels)
#   res
# }



