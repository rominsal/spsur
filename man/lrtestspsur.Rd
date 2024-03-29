% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrtestspsur.R
\name{lrtestspsur}
\alias{lrtestspsur}
\title{Likelihood Ratio tests for the specification of spatial SUR models.}
\usage{
lrtestspsur (objectr,  objectu = NULL)
}
\arguments{
\item{objectr}{An \code{spsur} object created by \code{\link{spsurml}},
\code{\link{spsur3sls}} or \code{\link{spsurtime}}.}

\item{objectu}{An \code{spsur} object created by \code{\link{spsurml}},
\code{\link{spsur3sls}} or \code{\link{spsurtime}} which
nests \code{objectr}. Default = \code{NULL}}
}
\value{
Object of \code{anova} class including, the list of models and, 
 for each model, the logLik, degrees of freedom and AIC. 
 If two nested models have been included as arguments, 
 the corresponding LR-test and its p-value associated.
}
\description{
The function computes a set of Likelihood Ratio tests, LR, 
 that help the user to select the spatial structure of the SUR model. 
 To achieve this goal, \code{\link{lrtestspsur}}
 needs to estimate the nested SUR models (options: "sim", "slx", "slm", 
 "sem", "sdm", "sdem", and "sarar"), using the function 
 \code{\link{spsurml}}.

 The nested models listed above are compared using ANOVA tables and 
 LR tests.
}
\details{
A fundamental result in maximum-likelihood estimation shows that 
if \emph{model A} is nested in \emph{model B}, by a set of \emph{n} 
restrictions on the parameters of \emph{model B}, then,
as the sample size increases, the test statistic: 
\emph{\eqn{-2log[l(H_{0}) / l(H_{A})]}}
is a \eqn{\chi^{2}(n)}, being l(H_{0} the estimated likelihood under 
the null hypothesis (\emph{model A}) and l(H_{A} the estimated likelihood 
under the alternative hypothesis (\emph{model B}).

 The list of (spatial) models that can be estimated with the function 
 \code{\link{spsurml}} includes the following (in addition to the 
 "slx" and "sdem"):

 \itemize{
    \item "sim": SUR model with no spatial effects
      \deqn{ y_{tg} = X_{tg} \beta_{g} + \epsilon_{tg} }
    \item "slm": SUR model with spatial lags of the explained variables
      \deqn{y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + \epsilon_{tg} }
    \item "sem": SUR model with spatial errors
      \deqn{ y_{tg} = X_{tg} \beta_{g} + u_{tg} }
      \deqn{ u_{tg} = \lambda_{g} Wu_{tg} + \epsilon_{tg} }
    \item "sdm": SUR model of the Spatial Durbin type
      \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tt} \beta_{g} + WX_{tg} \theta_{g} + \epsilon_{tg} }
    \item "sarar": SUR model with spatial lags of the explained variables and spatial
      errors
      \deqn{ y_{tg} = \rho_{g} Wy_{tg} + X_{tg} \beta_{g} + u_{tg} }
      \deqn{ u_{tg} = \lambda_{g} W u_{tg} + \epsilon_{tg} }
  }
  This collection of models can be compared, on objective bases, using the LR principle  and the
   following  nesting relations:

  \itemize{
    \item  "sim" vs either "slx",  slm", "sem", "sdm", "sarar"
    \item  "slm" vs either  "sdm", "sarar"
    \item  "sem" vs either  "sdm", "sdem", "sarar"
    \item  "slx" vs either  "sdm", "sdem"
  }
}
\examples{
#################################################
######## CROSS SECTION DATA (nG=1; nT>1) ########
#################################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data("spc", package = "spsur")
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
spcsur.slm <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc)
## ANOVA Table SLM model
lrtestspsur(spcsur.slm)    
## Test ANOVA SIM versus SLM
spcsur.sim <- spsurml(formula = Tformula, data = spc, 
                      type = "sim", listw = lwspc)
lrtestspsur(spcsur.sim, spcsur.slm)

## VIP: If you want to examine a particular example eliminate '#' and 
## execute the code of the example (they have been commented to 
## pass the checking time in CRAN)

## VIP: The output of the whole set of the examples can be examined 
## by executing demo(demo_lrtestspsur, package="spsur")

## Test ANOVA SLM vs SDM
#spcsur.sdm <- spsurml(formula = Tformula, data = spc, 
#                      type = "sdm", listw = lwspc)
#lrtestspsur(spcsur.slm, spcsur.sdm)
## Test ANOVA SEM vs SDM
#spcsur.sem <- spsurml(formula = Tformula, data = spc, 
#                      type = "sem", listw = lwspc)
#lrtestspsur(spcsur.sem, spcsur.sdm)
}
\references{
\itemize{
    \item Mur, J., Lopez, F., and Herrera, M. (2010). Testing for spatial
      effects in seemingly unrelated regressions.
      \emph{Spatial Economic Analysis}, 5(4), 399-440.
     \item Lopez, F.A., Mur, J., and Angulo, A. (2014). Spatial model
       selection strategies in a SUR framework. The case of regional
       productivity in EU. \emph{Annals of Regional Science}, 53(1),
       197-220.
  }
}
\seealso{
\code{\link{spsurml}}, \code{\link{lmtestspsur}}
}
\author{
\tabular{ll}{
  Fernando Lopez  \tab \email{fernando.lopez@upct.es} \cr
  Roman Minguez  \tab \email{roman.minguez@uclm.es} \cr
  Jesus Mur  \tab \email{jmur@unizar.es} \cr
  }
}
