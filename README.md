
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/f8l5h9/spqdep/workflows/R-CMD-check/badge.svg)](https://github.com/f8l5h9/spqdep/actions)
\[![CRAN status](https://www.r-pkg.org/badges/version-ago/spsur)
\[![CRAN
downloads-last-month](https://cranlogs.r-pkg.org/badges/last-month/spsur)
\[![CRAN
downloads-grand-total](https://cranlogs.r-pkg.org/badges/grand-total/spsur)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

# spsur

The **spsur** package allows for the estimation of the most popular
Spatial Seemingly Unrelated Regression models by maximum likelihood or
instrumental variable procedures (SUR-SLX; SUR-SLM; SUR-SEM;SUR-SDM;
SUR-SDEM and SUR-SARAR) and non spatial SUR model (SUR-SIM). Moreover,
**spsur** implements a collection of Lagrange Multipliers and Likelihood
Ratios to test for misspecifications in the SUR. Additional functions
allow for the estimation of the so-called spatial impacts (direct,
indirect and total effects) and also obtains random data sets, of a SUR
nature, with the features decided by the user. An important aspect of
**spsur** is that it operates both in a pure cross-sectional setting or
in panel data sets.

## Installation

You can install the released version of spsur from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spsur")
```

## Main functionalities of `spsur`

A few functions are necessary to test spatial autocorrelation in SUR and
estimate the spatial SUR models. Figure show the main functionalities of
the **spsur** package. The **spsur** package has one function to test
spatial autocorrelation on the residuals of a SUR model (lmtestspsur());
two functions to estimate SUR models with several spatial structures, by
maximum likelihood (spsurml()) and instrumental variables (spsur3sls());
three functions to help to the user to select the correct espeficication
(lrtestspsur(); wald\_betas() and wald\_deltas(); one function to get
the impacts (impacts()). Finally another function has been included in
this package to help to the user to develop Monte Carlo exercices
(dgp\_spsur()).

<img src="C:/Users/Fernando-pc/Dropbox/spSUR/notes/Journal of Statisitical Software/Figure1_JSS.png" title="\label{Fig1} Main functionatities of spsur package" alt="\label{Fig1} Main functionatities of spsur package" width="60%" style="display: block; margin: auto;" />

## Data sets in `spSUR`

The `spSUR` package include two data
sets:

#### The spc (Spatial Phillips-Curve). A classical data set from Anselin (1988, p.203)

A total of N=25 observations and Tm=2 time
periods

![](C:/Users/Fernando-pc/Dropbox/spSUR/notes/vignettes/spc.png)

| COUNTY   |   WAGE83 |     UN83 |      NMR83 | SMSA |   WAGE82 |   WAGE81 |     UN80 |      NMR80 |   WAGE80 |
| :------- | -------: | -------: | ---------: | ---: | -------: | -------: | -------: | ---------: | -------: |
| UNION    | 1.003127 | 0.080500 | \-0.002217 |    1 | 1.108662 | 1.146178 | 0.130375 | \-0.010875 | 1.084886 |
| DELAWARE | 1.039972 | 0.122174 |   0.018268 |    1 | 1.071271 | 1.104241 | 0.189603 |   0.041886 | 1.110426 |
| LICKING  | 1.050196 | 0.095821 | \-0.013681 |    1 | 1.058375 | 1.094732 | 0.124125 | \-0.004158 | 1.069776 |

#### Homicides + Socio-Economics characteristics for U.S. counties (1960-90)

from \[<https://geodacenter.github.io/data-and-lab/ncovr/>\]

Homicides and selected socio-economic characteristics for continental
U.S. counties. Data for four decennial census years: 1960, 1970, 1980
and 1990.  
A total of N=3,085 US
counties

![](C:/Users/Fernando-pc/Dropbox/spSUR/notes/vignettes/NAT.png)

0

| NAME              | STATE\_NAME | STATE\_FIPS | CNTY\_FIPS | FIPS  | STFIPS | COFIPS | FIPSNO | SOUTH |
| :---------------- | :---------- | :---------- | :--------- | :---- | -----: | -----: | -----: | ----: |
| Lake of the Woods | Minnesota   | 27          | 077        | 27077 |     27 |     77 |  27077 |     0 |
| Ferry             | Washington  | 53          | 019        | 53019 |     53 |     19 |  53019 |     0 |
| Stevens           | Washington  | 53          | 065        | 53065 |     53 |     65 |  53065 |     0 |

## How to specify multiple equations: The `Formula` package

By example: two equations with different number of
> regressors

> \(Y_{1} = \beta_{10} + \beta_{11} \ X_{11}+\beta_{12}X_{12}+\epsilon_{1}\)  
> \(Y_{2} = \beta_{20} + \beta_{21} \ X_{21}+\epsilon_{2}\)

> formula \<- \(Y_{1}\) | \(Y_{2}\) \~ \(X_{11}\) + \(X_{12}\)   |  
> \(X_{21}\)

Note that in the left side of the formula, two dependent variables has
been included separated by the symbol |. In right side, the independent
variables for each equation are included separated newly by a vertical
bar |, keeping the same order that in the left side.

-----

# The `spSUR` package step by step

### Step 1: Testing for spatial effects

### Step 2: Estimation of the Spatial SUR models

### Step 3: Looking for the correct especification

### Step 4: Impacts: Directs, Indirects and Total effects

### Step 5: `spSUR` in a panel data framework

### Step 6: Additional functionalities

### Step 7: Conclusion and work to do

-----

# Step 1: Testing for: `lmtestspsur`

The function `lmtestspsur()` obtain five LM statistis for testing
spatial dependence in Seemingly Unrelated Regression models  
(Mur J, López FA, Herrera M, 2010: Testing for spatial effect in
Seemingly Unrelated Regressions. *Spatial Economic Analysis* 5(4)
399-440).

> \(H_{0}:\) No spatial autocorrelation  
> \(H_{A}:\) SUR-SAR or  
> \(H_{A}:\) SUR-SEM or  
> \(H_{A}:\) SUR-SARAR

  - **LM-SUR-SAR**
  - **LM-SUR-SEM**
  - **LM-SUR-SARAR**

and two robust LM
tests

  - **LM\*-SUR-SAR**  
  - **LM\*-SUR-SEM**

#### Example 1: with Anselin’s data we can test spatial effects in the SUR model:

> \(WAGE_{83} = \beta_{10} + \beta_{11} \ UN_{83} + \beta_{12} \  NMR_{83} + \beta_{13} \ SMSA + \epsilon_{83}\)  
> \(WAGE_{81} = \beta_{20} + \beta_{21} \ UN_{80} + \beta_{22} \ NMR_{80} + \beta_{23} \ SMSA+ \epsilon_{81}\)  
> \(Corr(\epsilon_{83},\epsilon_{81}) \neq 0\)

``` r
library("spsur")
data("spc")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
LMs <- lmtestspsur(Form=Tformula,data=spc,W=Wspc)
#>              LM-Stat. DF p-value  
#> LM-SUR-SLM     5.2472  2  0.0725 .
#> LM-SUR-SEM     3.3050  2  0.1916  
#> LM*-SUR-SLM    2.1050  2  0.3491  
#> LM*-SUR-SEM    0.1628  2  0.9218  
#> LM-SUR-SARAR   5.7703  4  0.2170  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

In this example no spatial autocorrelation is identify\! (25
observations).

# Step 2: Estimation of a Spatial SUR

Two alternative estimation methods are implemented:

## 2.1. Maximum likelihood estimation: `spsurml`

Maximum likelihood estimation for differents spatial SUR models using
`spsurml`. The models are:

>   - **SUR-SIM:** with out spatial autocorrelation
>   - **SUR-SLX:** Spatial Lag of X SUR model
>   - **SUR-SLM:** Spatial Autorregresive SUR model
>   - **SUR-SEM:** Spatial Error SUR model
>   - **SUR-SDM:** Spatial Durbin SUR model
>   - **SUR-SDEM:** Spatial Durbin Error SUR model
>   - **SUR-SARAR:** Spatial Autorregresive with Spatial Error SUR model

### 2.1.1 Anselin data set

-----

**SUR-SAR:** Spatial autorregresive
model:  
(\(y_{t} = \lambda_{t} Wy_{t} + X_{t} \beta_{t} + \epsilon_{t}; \ t=1,...,T\)
)

> \(WAGE_{83} = \lambda_{83} WWAGE_{83} + \beta_{10} + \beta_{11} UN_{83} + \beta_{12} NMR_{83} + \beta_{13} SMSA + \epsilon_{83}\)  
> \(WAGE_{81} = \lambda_{81} WWAGE_{81} + \beta_{20} + \beta_{21} UN_{80} + \beta_{22} NMR_{80} + \beta_{23} SMSA+ \epsilon_{81}\)  
> \(Corr(\epsilon_{83},\epsilon_{81}) \neq 0\)

``` r
spcSUR.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc)
#> Initial point:   log_lik:  113.197  lambdas:  -0.472 -0.446 
#> Iteration:  1   log_lik:  114.085  lambdas:  -0.506 -0.482 
#> Iteration:  2   log_lik:  114.096  lambdas:  -0.506 -0.482 
#> Time to fit the model:  3.11  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.4  seconds
summary(spcSUR.slm)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "slm")
#> 
#>  
#> Spatial SUR model type:  slm 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1  1.4955217  0.2467240  6.0615 5.183e-07 ***
#> UN83_1         0.8070029  0.2557439  3.1555  0.003179 ** 
#> NMR83_1       -0.5194114  0.2590550 -2.0050  0.052318 .  
#> SMSA_1        -0.0073247  0.0118519 -0.6180  0.540347    
#> lambda_1      -0.5057334  0.2405734 -2.1022  0.042401 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6224 
#>   Equation  2 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2  1.7094414  0.2925620  5.8430 1.024e-06 ***
#> UN80_2        -0.6745562  0.3870737 -1.7427   0.08969 .  
#> NMR80_2        0.7502934  0.3842670  1.9525   0.05847 .  
#> SMSA_2         0.0014181  0.0241859  0.0586   0.95356    
#> lambda_2      -0.4821428  0.2557758 -1.8850   0.06730 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.4743 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0003085954 -0.0003561928
#>  -0.0003561928  0.0015864976
#> Correlation Matrix of inter-equation residuals:                    
#>   1.000000 -0.509062
#>  -0.509062  1.000000
#> 
#>  R-sq. pooled: 0.6603 
#>  Log-Likelihood:  114.096
#>  Breusch-Pagan: 6.516  p-value: (0.0107) 
#>  LMM: 0.50489  p-value: (0.477)
```

-----

Only change the **‘type’** argument in `spsurml` function it is possible
to estimate several spatial model

**SUR-SEM:** Spatial error
model:  
(\(y_{t} = X_{t}\ \beta_{t} + u_{t}\ ; u_{t}=\rho u_{t}+ \epsilon_{t}\ t=1,...,T\))

> \(WAGE_{83} = \beta_{10} + \beta_{11} UN_{83} + \beta_{12} NMR_{83} + \beta_{13} SMSA + u_{83}; \ u_{83}=\rho W \ u_{83} + \epsilon_{83}\)  
> \(WAGE_{81} =\beta_{20} + \beta_{21} UN_{80} + \beta_{22} NMR_{80} + \beta_{23} SMSA+ u_{81}; \ u_{81}=\rho W \ u_{81} + \epsilon_{81}\)  
> \(Corr(\epsilon_{83},\epsilon_{81}) \neq 0\)

``` r
spcSUR.sem <-spsurml(Form=Tformula,data=spc,type="sem",W=Wspc)
#> Initial point:   log_lik:  112.821  deltas:  -0.556 -0.477 
#> Iteration:  1  log_lik:  113.695  rhos:  -0.618 -0.537 
#> Iteration:  2  log_lik:  113.719  rhos:  -0.628 -0.548 
#> Time to fit the model:  3.99  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.36  seconds
summary(spcSUR.sem)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "sem")
#> 
#>  
#> Spatial SUR model type:  sem 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1  0.9805218  0.0151899 64.5511 < 2.2e-16 ***
#> UN83_1         0.7383349  0.2277247  3.2422  0.002513 ** 
#> NMR83_1       -0.4859228  0.2550377 -1.9053  0.064535 .  
#> SMSA_1        -0.0132403  0.0099122 -1.3358  0.189790    
#> rho_1         -0.6280610  0.2774391 -2.2638  0.029541 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6607 
#>   Equation  2 
#>                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)_2  1.1479484  0.0386961 29.6657  < 2e-16 ***
#> UN80_2        -0.4406330  0.3614882 -1.2189  0.23058    
#> NMR80_2        0.8223976  0.4062173  2.0245  0.05018 .  
#> SMSA_2         0.0041942  0.0204639  0.2050  0.83873    
#> rho_2         -0.5480668  0.2817155 -1.9455  0.05935 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.5092 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0002970481 -0.0003217158
#>  -0.0003217158  0.0015512097
#> Correlation Matrix of inter-equation residuals:                      
#>   1.0000000 -0.4739403
#>  -0.4739403  1.0000000
#> 
#>  R-sq. pooled: 0.673 
#>  Log-Likelihood:  113.719
#>  Breusch-Pagan: 5.512  p-value: (0.0189) 
#>  LMM: 1.4742  p-value: (0.225)
```

# Step 3: Testing for misspecification in spatial SUR

## 3.1 Testing for the diagonality of \(\Sigma\)

The Breush-Pagan test of diagonality of \(\Sigma\)

> \(H_{0}: \Sigma = \sigma^2 I_{R}\)  
> \(H_{A}: \Sigma \neq \sigma^2 I_{R}\)

## 3.2 Marginal tests: \(LM(\rho|\lambda)\) & \(LM(\lambda|\rho)\)

The Marginal Multiplier tests (LMM) are used to test for no correlation
in one part of the model allowing for spatial correlation in the other.

-----

>   - The \(LM(\rho|\lambda)\) is the test for spatial error correlation
>     in a model with subtantive spatial correlation (SUR-SAR; SUR-SDM).

> \(H_{0}: SUR-SAR\)  
> \(H_{A}: SUR-SARAR\)

-----

>   - The \(LM(\lambda|\rho)\) is the test for subtantive spatial
>     autocorrelation in a model with spatial autocorrelation in error
>     term (SUR-SEM; SUR-SDEM).

> \(H_{0}: SUR-SEM\)  
> \(H_{A}: SUR-SARAR\)

``` r
summary(spcSUR.sem)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "sem")
#> 
#>  
#> Spatial SUR model type:  sem 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1  0.9805218  0.0151899 64.5511 < 2.2e-16 ***
#> UN83_1         0.7383349  0.2277247  3.2422  0.002513 ** 
#> NMR83_1       -0.4859228  0.2550377 -1.9053  0.064535 .  
#> SMSA_1        -0.0132403  0.0099122 -1.3358  0.189790    
#> rho_1         -0.6280610  0.2774391 -2.2638  0.029541 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6607 
#>   Equation  2 
#>                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)_2  1.1479484  0.0386961 29.6657  < 2e-16 ***
#> UN80_2        -0.4406330  0.3614882 -1.2189  0.23058    
#> NMR80_2        0.8223976  0.4062173  2.0245  0.05018 .  
#> SMSA_2         0.0041942  0.0204639  0.2050  0.83873    
#> rho_2         -0.5480668  0.2817155 -1.9455  0.05935 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.5092 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0002970481 -0.0003217158
#>  -0.0003217158  0.0015512097
#> Correlation Matrix of inter-equation residuals:                      
#>   1.0000000 -0.4739403
#>  -0.4739403  1.0000000
#> 
#>  R-sq. pooled: 0.673 
#>  Log-Likelihood:  113.719
#>  Breusch-Pagan: 5.512  p-value: (0.0189) 
#>  LMM: 1.4742  p-value: (0.225)
```

-----

## 3.3 Coefficient stability/homogeneity

### 3.3.1 Wald tests for beta coefficients: `wald_betas`

In a **SUR-SAR** the
> model:

> \(WAGE_{83} = \lambda_{83} W \ WAGE_{83} + \beta_{10} + \beta_{11} UN_{83} + \beta_{12} NMR_{83} + \boldsymbol{\beta_{13}} SMSA + \epsilon_{83}\)  
> \(WAGE_{81} = \lambda_{81} W \ WAGE_{81} + \beta_{20} + \beta_{21} UN_{80} + \beta_{22} NMR_{80} + \boldsymbol{\beta_{23}} SMSA + \epsilon_{81}\)

It’s possible to test equality between SMSA coefficients in both
equations:

> \(H_{0}: \beta_{13} = \beta_{23}\)  
> \(H_{A}: \beta_{13} \neq \beta_{23}\)

``` r
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
spcSUR.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc)
#> Initial point:   log_lik:  113.197  lambdas:  -0.472 -0.446 
#> Iteration:  1   log_lik:  114.085  lambdas:  -0.506 -0.482 
#> Iteration:  2   log_lik:  114.096  lambdas:  -0.506 -0.482 
#> Time to fit the model:  2.89  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.2  seconds
R1 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
b1 <- matrix(0,ncol=1)
Wald_beta <- wald_betas(results=spcSUR.slm,R=R1,b=b1)
#> Wald stat.: 0.079 p-value: (0.779)
```

More complex hypothesis about \(\beta\) coefficients could be tested
using R1 vector

> \(H_{0}: \beta_{13} = \beta_{23}\) and \(\beta_{12} = \beta_{22}\)  
> \(H_{A}: \beta_{13} \neq \beta_{23}\) or
> \(\beta_{12} \neq \beta_{22}\)

``` r
# Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
# spcSUR.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc)
R1 <- t(matrix(c(0,0,0,1,0,0,0,-1,0,0,1,0,0,0,-1,0),ncol=2))
b1 <- matrix(0,ncol=2)
print(R1)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    0    0    0    1    0    0    0   -1
#> [2,]    0    0    1    0    0    0   -1    0
Wald_beta <- wald_betas(results=spcSUR.slm,R=R1,b=b1)
#> Wald stat.: 6.179 p-value: (0.046)
```

### Estimate the restricted model

> In case don’t reject the null, it’s possible to estimate the model
> with equal coefficient in both
> equations:

> \(WAGE_{83} = \lambda_{83} W \ WAGE_{83} + \beta_{10} + \beta_{11} UN_{83} + \beta_{12} NMR_{83} + \boldsymbol{\beta_{13}} SMSA + \epsilon_{83}\)  
> \(WAGE_{81} = \lambda_{81} W \ WAGE_{81} + \beta_{20} + \beta_{21} UN_{80} + \beta_{22} NMR_{80} + \boldsymbol{\beta_{13}} SMSA + \epsilon_{81}\)

``` r
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
R1 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
b1 <- matrix(0,ncol=1)
spcSUR.sar.restring <- spsurml(Form=Tformula, data=spc, type="slm", W=Wspc,R=R1,b=b1)
#> Initial point:   log_lik:  113.161  lambdas:  -0.428 -0.421 
#> Iteration:  1   log_lik:  114.01  lambdas:  -0.482 -0.465 
#> Iteration:  2   log_lik:  114.049  lambdas:  -0.495 -0.476 
#> Time to fit the model:  2.8  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.21  seconds
summary(spcSUR.sar.restring)
#> Call:
#> spsurml(Form = Tformula, data = spc, R = R1, b = b1, W = Wspc, 
#>     type = "slm")
#> 
#>  
#> Spatial SUR model type:  slm 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1  1.4876833  0.2444552  6.0857 4.343e-07 ***
#> UN83_1         0.7486660  0.2290706  3.2683  0.002301 ** 
#> NMR83_1       -0.4565917  0.2559136 -1.7842  0.082384 .  
#> SMSA_1        -0.0041128  0.0081857 -0.5024  0.618257    
#> lambda_1      -0.4950681  0.2392546 -2.0692  0.045377 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6115 
#>   Equation  2 
#>               Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2  1.69869    0.28885  5.8809 8.294e-07 ***
#> UN80_2        -0.61498    0.32522 -1.8910   0.06627 .  
#> NMR80_2        0.68793    0.37391  1.8398   0.07362 .  
#> NA                  NA         NA      NA        NA    
#> lambda_2      -0.47556    0.25016 -1.9010   0.06490 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.4691 
#>   Variance-Covariance Matrix of inter-equation residuals:                           
#>   0.0003220152 -0.000394289
#>  -0.0003942890  0.001616933
#> Correlation Matrix of inter-equation residuals:                      
#>   1.0000000 -0.5464249
#>  -0.5464249  1.0000000
#> 
#>  R-sq. pooled: 0.6531 
#>  Log-Likelihood:  114.049
#>  Breusch-Pagan: 7.318  p-value: (0.00683) 
#>  LMM: 0.36616  p-value: (0.545)
```

### 3.3.2 Wald test for ‘spatial’ coefficients homogeneity: `wald_deltas`

In same way a test for equal spatial autocorrelation coefficients can be
obtain with `wald_deltas` function:  
In the
> model:

> \(WAGE_{83} = \boldsymbol{\lambda_{83}} W \ WAGE_{83} + \beta_{10} + \beta_{11} UN_{83} + \beta_{12} NMR_{83} + \beta_{13} SMSA + \epsilon_{83}\)  
> \(WAGE_{81} = \boldsymbol{\lambda_{81}} W \ WAGE_{81} + \beta_{20} + \beta_{21} UN_{80} + \beta_{22} NMR_{80} + \beta_{23} SMSA + \epsilon_{81}\)

In this case the null
is:

> \(H_{0}: \lambda_{83} = \lambda_{81}\)  
> \(H_{A}: \lambda_{83} \neq \lambda_{81}\)

``` r
# Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
# spcSUR.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc,trace=F)
R1 <- matrix(c(1,-1),nrow=1)
b1 <- matrix(0,ncol=1)
res1 <- wald_deltas(results=spcSUR.slm,R=R1,b=b1)
#> 
#>  Wald stat.: 0.006 (0.939)
```

### 3.3.3 Likekihood ratio tests `lr_betas_spsur`

Alternatively to wald test, the Likelihoo Ration (LR) tests can be
obtain using the `lr_betas_spsur` function.

``` r
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
R1 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
b1 <- matrix(0,ncol=1)
LR_SMSA <-  lr_betas_spsur(Form=Tformula,data=spc,W=Wspc,type="slm",R=R1,b=b1,trace=F,printmodels=F)
#> 
#>  Fitting unrestricted model ... 
#> 
#>  Time to fit unrestricted model:  2.88  seconds
#> 
#>  Fitting restricted model ... 
#> Time to fit restricted model:  2.83  seconds
#> 
#>  LR-Test 
#> 
#>  Log-likelihood unrestricted model:  114.096
#>  Log-likelihood restricted model:  114.049
#>  LR statistic:  0.095  degrees of freedom:  1  p-value: ( 0.7576548 )
```

# 4\. Step 4: Marginal Effects: `impacts`

The marginal effects `impacts` of spatial autoregressive models
(SUR-SAR; SUR-SDM; SUR-SARAR) has been calculated following the propose
of LeSage and Pace (2009).

``` r
eff.spcSUR.sar <-impacts(spcSUR.slm,nsim=299)
#> 
#> Spatial SUR model type:  slm 
#> 
#>  Direct effects 
#> 
#>                mean          sd  t-stat     p-val    
#> UN83_1   8.6479e-01  2.5919e-01  3.3365 0.0008484 ***
#> NMR83_1 -5.7166e-01  2.9096e-01 -1.9647 0.0494447 *  
#> SMSA_1  -7.8388e-03  1.2291e-02 -0.6378 0.5236346    
#> UN80_2  -6.8648e-01  4.0118e-01 -1.7112 0.0870514 .  
#> NMR80_2  7.9222e-01  4.2757e-01  1.8528 0.0639067 .  
#> SMSA_2   1.7235e-05  2.4841e-02  0.0007 0.9994464    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  Indirect effects 
#> 
#>                mean          sd  t-stat   p-val  
#> UN83_1  -0.31238169  0.16813634 -1.8579 0.06318 .
#> NMR83_1  0.20488538  0.14047116  1.4586 0.14469  
#> SMSA_1   0.00224931  0.00480408  0.4682 0.63964  
#> UN80_2   0.22953581  0.19511769  1.1764 0.23944  
#> NMR80_2 -0.25803529  0.21139368 -1.2206 0.22222  
#> SMSA_2   0.00093014  0.00929599  0.1001 0.92030  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>  Total effects 
#> 
#>                mean          sd  t-stat    p-val   
#> UN83_1   0.55240913  0.19422962  2.8441 0.004454 **
#> NMR83_1 -0.36677246  0.20738882 -1.7685 0.076973 . 
#> SMSA_1  -0.00558949  0.00835723 -0.6688 0.503609   
#> UN80_2  -0.45694388  0.28367671 -1.6108 0.107225   
#> NMR80_2  0.53418244  0.30849200  1.7316 0.083346 . 
#> SMSA_2   0.00094737  0.01705941  0.0555 0.955713   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

# 5\. Step 5: The `spSUR` in a panel data framework

## 5.1 The `spSUR` with G equation and T periods

Case of T temporal cross-sections and G equations

![](/Users/roman/Dropbox/SpSUR/spsur2/vignettes/PanelGT.png)

By example with NAT data set:

>   - T = 4 (Four temporal periods)  
>   - G = 2 (Two equations with different numbers of independent
>     variables)  
>   - R = 3085 (Spatial observations)

A **SUR-SLM-PANEL**
model  
\(y_{gt} = \lambda_{g} Wy_{gt} + X_{gt} \beta_{g} + \epsilon_{gt};\)  
\(\ g=1,...,G; \ t=1,...,T\)  
\(Corr(\epsilon_{gt},\epsilon_{g't})=Corr(\epsilon_{g},\epsilon_{g'}) \neq 0\)
for \((\forall t)\)

``` r
Tformula <- HR80 | HC80 ~ UE80 + RD80 | UE80
spSUR.sar.panel <- spsurml(Form = Tformula,data=NCOVR,W=W,type="slm",N=3085,G=2,Tm=4)
#> Initial point:   log_lik:  -26267.63  lambdas:  0.376 0.145 
#> Iteration:  1   log_lik:  -26253.87  lambdas:  0.389 0.146 
#> Iteration:  2   log_lik:  -26253.87  lambdas:  0.389 0.146 
#> Time to fit the model:  6.3  seconds
#> Computing marginal test... 
#> Time to compute covariances:  37.11  seconds
summary(spSUR.sar.panel)
#> Call:
#> spsurml(Form = Tformula, data = NCOVR, W = W, G = 2, N = 3085, 
#>     Tm = 4, type = "slm")
#> 
#>  
#> Spatial SUR model type:  slm 
#> 
#> Equation  1 
#>               Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1 3.664135   0.267777 13.6835 < 2.2e-16 ***
#> UE80_1        0.081705   0.029855  2.7367  0.006224 ** 
#> RD80_1        2.537307   0.115395 21.9881 < 2.2e-16 ***
#> lambda_1      0.389281   0.021784 17.8697 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.3779 
#>   Equation  2 
#>               Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2 5.602462   2.232275  2.5098   0.01211 *  
#> UE80_2        0.110066   0.294648  0.3736   0.70875    
#> lambda_2      0.145692   0.028112  5.1825 2.259e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.01993 
#>   Variance-Covariance Matrix of inter-equation residuals:                    
#>  29.07190   60.48062
#>  60.48062 2940.02430
#> Correlation Matrix of inter-equation residuals:                    
#>  1.0000000 0.2068731
#>  0.2068731 1.0000000
#> 
#>  R-sq. pooled: 0.02445 
#>  Log-Likelihood:  -26253.9
#>  Breusch-Pagan:   132  p-value: (1.48e-30) 
#>  LMM: 371.03  p-value: (1.12e-82)
```

### 5.1.1 Unobserved effects: The `demaining` option

`spsur` package offers the possibility to transform the original data in
order to remove potential unobserved effects.  
The most popular transformation is demeaning the data: To subtract the
sampling averages of each individual, in every equation, from the
corresponding observation.

``` r
Tformula <- HR80 | HC80 ~ UE80 + RD80 | UE80
spSUR.sar.panel <- spsurml(Form = Tformula,data=NCOVR,W=W,type="slm",N=3085,G=2,Tm=4,demean = T)
#> Initial point:   log_lik:  -26267.63  lambdas:  0.376 0.145 
#> Iteration:  1   log_lik:  -26253.87  lambdas:  0.389 0.146 
#> Iteration:  2   log_lik:  -26253.87  lambdas:  0.389 0.146 
#> Time to fit the model:  6.51  seconds
#> Computing marginal test... 
#> Time to compute covariances:  35.82  seconds
summary(spSUR.sar.panel)
#> Call:
#> spsurml(Form = Tformula, data = NCOVR, W = W, G = 2, N = 3085, 
#>     Tm = 4, demean = T, type = "slm")
#> 
#>  
#> Spatial SUR model type:  slm 
#> 
#> Equation  1 
#>               Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1 3.664135   0.267777 13.6835 < 2.2e-16 ***
#> UE80_1        0.081705   0.029855  2.7367  0.006224 ** 
#> RD80_1        2.537307   0.115395 21.9881 < 2.2e-16 ***
#> lambda_1      0.389281   0.021784 17.8697 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.3779 
#>   Equation  2 
#>               Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2 5.602462   2.232275  2.5098   0.01211 *  
#> UE80_2        0.110066   0.294648  0.3736   0.70875    
#> lambda_2      0.145692   0.028112  5.1825 2.259e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.01993 
#>   Variance-Covariance Matrix of inter-equation residuals:                    
#>  29.07190   60.48062
#>  60.48062 2940.02430
#> Correlation Matrix of inter-equation residuals:                    
#>  1.0000000 0.2068731
#>  0.2068731 1.0000000
#> 
#>  R-sq. pooled: 0.02445 
#>  Log-Likelihood:  -26253.9
#>  Breusch-Pagan:   132  p-value: (1.48e-30) 
#>  LMM: 371.03  p-value: (1.12e-82)
```

# 6\. Aditional functionalities: `dgp_spSUR`

A Data Generating Process of a spatial SUR models is avalible using the
function `dgp_spSUR`

``` r
nT <- 1 # Number of time periods
nG <- 3 # Number of equations
nR <- 500 # Number of spatial elements
p <- 3 # Number of independent variables
Sigma <- matrix(0.3,ncol=nG,nrow=nG)
diag(Sigma)<-1
Coeff <- c(2,3)
rho <- 0.5 # level of spatial dependence
lambda <- 0.0 # spatial autocorrelation error term = 0
# ramdom coordinates
# co <- cbind(runif(nR,0,1),runif(nR,0,1))
# W <- spdep::nb2mat(spdep::knn2nb(spdep::knearneigh(co,k=5,longlat=F)))
# DGP <- dgp_spSUR(Sigma=Sigma,Coeff=Coeff,rho=rho,lambda=lambda,nT=nT,nG=nG,nR=nR,p=p,W=W)
```

# 7\. Conclusion & work to do

>   - `spSUR` is a powerful R-package to test, estimate and looking for
>     the correct specification  
>   - More functionalities and estimation algorithm coming soon
>       - GMM estimation  
>       - ML estimation with equal level of spatial dependence
>         (\(\lambda\)/\(\rho\)=constant)
>       - Orthogonal deamining for space-time SUR models
>       - …..
>   - The spSUR is avaliable in GitHub
>     \[<https://github.com/rominsal/spSUR/>\]

-----

# References

  - López, F.A., P. J. Martínez-Ortiz, and J.G. Cegarra-Navarro (2017).
    Spatial spillovers in public expenditure on a municipal level in
    spain. *The Annals of Regional Science* 58 (1), 39–65.  
  - López, F.A., J. Mur, and A. Angulo (2014). Spatial model selection
    strategies in a sur framework. the case of regional productivity in
    eu. *The Annals of Regional Science* 53 (1), 197–220.  
  - Mur, J., F. López, and M. Herrera (2010). Testing for spatial
    effects in seemingly unrelated regressions. *Spatial Economic
    Analysis* 5 (4), 399–440.

## Example

This is a basic example which shows you how test spatial structure in
the residual of a SUR model:

``` r
## Testing for spatial effects in SUR model
library("spsur")
data("spc")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
LMs <- lmtestspsur(Form = Tformula, data = spc, W = Wspc)
#>              LM-Stat. DF p-value  
#> LM-SUR-SLM     5.2472  2  0.0725 .
#> LM-SUR-SEM     3.3050  2  0.1916  
#> LM*-SUR-SLM    2.1050  2  0.3491  
#> LM*-SUR-SEM    0.1628  2  0.9218  
#> LM-SUR-SARAR   5.7703  4  0.2170  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Several Spatial SUR model can be estimated by maximun likelihhod:

``` r
## A SUR-SLM model
spcsur.slm <-spsurml(Form = Tformula, data = spc, type = "slm", W = Wspc)
#> Initial point:   log_lik:  113.197  lambdas:  -0.472 -0.446 
#> Iteration:  1   log_lik:  114.085  lambdas:  -0.506 -0.482 
#> Iteration:  2   log_lik:  114.096  lambdas:  -0.506 -0.482 
#> Time to fit the model:  3.04  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.21  seconds
summary(spcsur.slm)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "slm")
#> 
#>  
#> Spatial SUR model type:  slm 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_1  1.4955217  0.2467240  6.0615 5.183e-07 ***
#> UN83_1         0.8070029  0.2557439  3.1555  0.003179 ** 
#> NMR83_1       -0.5194114  0.2590550 -2.0050  0.052318 .  
#> SMSA_1        -0.0073247  0.0118519 -0.6180  0.540347    
#> lambda_1      -0.5057334  0.2405734 -2.1022  0.042401 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6224 
#>   Equation  2 
#>                 Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2  1.7094414  0.2925620  5.8430 1.024e-06 ***
#> UN80_2        -0.6745562  0.3870737 -1.7427   0.08969 .  
#> NMR80_2        0.7502934  0.3842670  1.9525   0.05847 .  
#> SMSA_2         0.0014181  0.0241859  0.0586   0.95356    
#> lambda_2      -0.4821428  0.2557758 -1.8850   0.06730 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.4743 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0003085954 -0.0003561928
#>  -0.0003561928  0.0015864976
#> Correlation Matrix of inter-equation residuals:                    
#>   1.000000 -0.509062
#>  -0.509062  1.000000
#> 
#>  R-sq. pooled: 0.6603 
#>  Log-Likelihood:  114.096
#>  Breusch-Pagan: 6.516  p-value: (0.0107) 
#>  LMM: 0.50489  p-value: (0.477)
```
