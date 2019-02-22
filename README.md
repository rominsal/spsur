# spsur
Spatial Seemingly Unrelated Regressions
The **spsur** package allows for the estimation of the most popular Spatial Seemingly Unrelated Regression models 
by maximum likelihood or instrumental variable procedures (SUR-SLX; SUR-SLM; SUR-SEM;SUR-SDM; SUR-SDEM and SUR-SARAR) 
and non spatial SUR model (SUR-SIM). Moreover, **spsur** implements a collection of Lagrange Multipliers and Likelihood 
Ratios to test for misspecifications in the SUR. Additional functions allow for the estimation of the so-called spatial 
impacts (direct, indirect and total effects) and also obtains random data sets, of a SUR nature, with the features decided 
by the user. An important aspect of **spsur** is that it operates both in a pure cross-sectional setting or in panel data sets.

Installation
------------

You can install the released version of spsur from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("spsur")
```

Main functionalities of **spsur**
-------------------------------

A few functions are necessary to test spatial autocorrelation in SUR and estimate the spatial SUR models. 
The **spsur** package has one function to test spatial autocorrelation on the residuals of a SUR model (`lmtestspsur()`); 
two functions to estimate SUR models with several spatial structures, by maximum likelihood (spsurml()) and instrumental 
variables (`spsur3sls()`); three functions to help to the user to select the correct espeficication (`lrtestspsur()`; 
`wald_betas()` and `wald_deltas()`; one function to get the impacts (`impacts()`). Finally another function has been 
included in this package to help to the user to develop Monte Carlo exercices (`dgp_spsur()`).


Data sets in **spsur**
--------------------

The **spsur** package include two data sets:

#### The spc (Spatial Phillips-Curve). A classical data set from Anselin (1988, p.203)

A total of N=25 observations and Tm=2 time periods


| COUNTY   |    WAGE83|      UN83|      NMR83|  SMSA|    WAGE82|    WAGE81|      UN80|      NMR80|    WAGE80|
|:---------|---------:|---------:|----------:|-----:|---------:|---------:|---------:|----------:|---------:|
| UNION    |  1.003127|  0.080500|  -0.002217|     1|  1.108662|  1.146178|  0.130375|  -0.010875|  1.084886|
| DELAWARE |  1.039972|  0.122174|   0.018268|     1|  1.071271|  1.104241|  0.189603|   0.041886|  1.110426|
| LICKING  |  1.050196|  0.095821|  -0.013681|     1|  1.058375|  1.094732|  0.124125|  -0.004158|  1.069776|

#### Homicides + Socio-Economics characteristics for U.S. counties (1960-90)

from \[<https://geodacenter.github.io/data-and-lab/ncovr/>\]

Homicides and selected socio-economic characteristics for continental U.S. counties. Data for four decennial census years: 
1960, 1970, 1980 and 1990. A total of N=3,085 US counties


| NAME              | STATE\_NAME | STATE\_FIPS | CNTY\_FIPS | FIPS  |  STFIPS|  COFIPS|  FIPSNO|  SOUTH|
|:------------------|:------------|:------------|:-----------|:------|-------:|-------:|-------:|------:|
| Lake of the Woods | Minnesota   | 27          | 077        | 27077 |      27|      77|   27077|      0|
| Ferry             | Washington  | 53          | 019        | 53019 |      53|      19|   53019|      0|
| Stevens           | Washington  | 53          | 065        | 53065 |      53|      65|   53065|      0|

How to specify multiple equations: The `Formula` package
--------------------------------------------------------

By example: two equations with different number of regressors

> *Y*<sub>1</sub> = *β*<sub>10</sub> + *β*<sub>11</sub> *X*<sub>11</sub> + *β*<sub>12</sub>*X*<sub>12</sub> + *ϵ*<sub>1</sub>
> *Y*<sub>2</sub> = *β*<sub>20</sub> + *β*<sub>21</sub> *X*<sub>21</sub> + *ϵ*<sub>2</sub>

> formula &lt;- *Y*<sub>1</sub> | *Y*<sub>2</sub> ~ *X*<sub>11</sub> + *X*<sub>12</sub>   |   *X*<sub>21</sub>

Note that in the left side of the formula, two dependent variables has been included separated by the symbol |. In right side, the independent variables for each equation are included separated newly by a vertical bar |, keeping the same order that in the left side.

------------------------------------------------------------------------

The **spsur** package step by step
================================

### Step 1: Testing for spatial effects

### Step 2: Estimation of the Spatial SUR models

### Step 3: Looking for the correct especification

### Step 4: Impacts: Directs, Indirects and Total effects

### Step 5: `spSUR` in a panel data framework

### Step 6: Additional functionalities

### Step 7: Conclusion and work to do

------------------------------------------------------------------------

Step 1: Testing for: `lmtestspsur`
==================================

The function `lmtestspsur()` obtain five LM statistis for testing spatial dependence in Seemingly Unrelated Regression models
(Mur J, López FA, Herrera M, 2010: Testing for spatial effect in Seemingly Unrelated Regressions. *Spatial Economic Analysis* 5(4) 399-440).

> *H*<sub>0</sub>: No spatial autocorrelation
> *H*<sub>*A*</sub>: SUR-SAR or
> *H*<sub>*A*</sub>: SUR-SEM or
> *H*<sub>*A*</sub>: SUR-SARAR

-   **LM-SUR-SAR**
-   **LM-SUR-SEM**
-   **LM-SUR-SARAR**

and two robust LM tests

-   **LM\*-SUR-SAR**
-   **LM\*-SUR-SEM**

#### Example 1: with Anselin's data we can test spatial effects in the SUR model:

> *W**A**G**E*<sub>83</sub> = *β*<sub>10</sub> + *β*<sub>11</sub> *U**N*<sub>83</sub> + *β*<sub>12</sub> *N**M**R*<sub>83</sub> + *β*<sub>13</sub> *S**M**S**A* + *ϵ*<sub>83</sub>
> *W**A**G**E*<sub>81</sub> = *β*<sub>20</sub> + *β*<sub>21</sub> *U**N*<sub>80</sub> + *β*<sub>22</sub> *N**M**R*<sub>80</sub> + *β*<sub>23</sub> *S**M**S**A* + *ϵ*<sub>81</sub>
> *C**o**r**r*(*ϵ*<sub>83</sub>, *ϵ*<sub>81</sub>)≠0

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

In this example no spatial autocorrelation is identify! (25 observations).

Step 2: Estimation of a Spatial SUR
===================================

Two alternative estimation methods are implemented:

2.1. Maximum likelihood estimation: `spsurml`
---------------------------------------------

Maximum likelihood estimation for differents spatial SUR models using `spsurml`. The models are:

> -   **SUR-SIM:** with out spatial autocorrelation
> -   **SUR-SLX:** Spatial Lag of X SUR model
> -   **SUR-SAR:** Spatial Autorregresive SUR model
> -   **SUR-SEM:** Spatial Error SUR model
> -   **SUR-SDM:** Spatial Durbin SUR model
> -   **SUR-SDEM:** Spatial Durbin Error SUR model
> -   **SUR-SARAR:** Spatial Autorregresive with Spatial Error SUR model / (SUR-SAC)

### 2.1.1 Anselin data set

------------------------------------------------------------------------

**SUR-SIM:** SUR model without spatial effects
(*y*<sub>*t*</sub> = *X*<sub>*t*</sub>*β*<sub>*t*</sub> + *ϵ*<sub>*t*</sub>; *t* = 1, ..., *T*)

> *W**A**G**E*<sub>83</sub> = *β*<sub>10</sub> + *β*<sub>11</sub>*U**N*<sub>83</sub> + *β*<sub>12</sub>*N**M**R*<sub>83</sub> + *β*<sub>13</sub>*S**M**S**A* + *ϵ*<sub>83</sub>
> *W**A**G**E*<sub>81</sub> = *β*<sub>20</sub> + *β*<sub>21</sub>*U**N*<sub>80</sub> + *β*<sub>22</sub>*N**M**R*<sub>80</sub> + *β*<sub>23</sub>*S**M**S**A* + *ϵ*<sub>81</sub>
> *C**o**r**r*(*ϵ*<sub>83</sub>, *ϵ*<sub>81</sub>)≠0

``` r
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
spcSUR.sim <-spsurml(Form=Tformula,data=spc,type="sim",W=Wspc)
#> Initial point:  
#> log_lik:  110.423 
#> Iteration:  1  log_lik:  111.201 
#> Iteration:  2  log_lik:  111.348 
#> Iteration:  3  log_lik:  111.378 
#> Time to fit the model:  0.09  seconds
#> Time to compute covariances:  0.02  seconds
summary(spcSUR.sim)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "sim")
#> 
#>  
#> Spatial SUR model type:  sim 
#> 
#> Equation  1 
#>                 Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)_1  0.9845592  0.0204381 48.1728  < 2e-16 ***
#> UN83_1         0.6595414  0.2878689  2.2911  0.02744 *  
#> NMR83_1       -0.3772648  0.2699660 -1.3975  0.17018    
#> SMSA_1        -0.0085762  0.0131588 -0.6517  0.51839    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.4224 
#>   Equation  2 
#>                Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)_2  1.159370   0.045728 25.3535   <2e-16 ***
#> UN80_2        -0.559787   0.414434 -1.3507   0.1846    
#> NMR80_2        0.583248   0.381367  1.5294   0.1342    
#> SMSA_2         0.010439   0.026166  0.3990   0.6921    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.2757 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0004539866 -0.0006769556
#>  -0.0006769556  0.0021133889
#> Correlation Matrix of inter-equation residuals:                      
#>   1.0000000 -0.6911129
#>  -0.6911129  1.0000000
#> 
#>  R-sq. pooled: 0.5328 
#>  Log-Likelihood:  111.378
#>  Breusch-Pagan: 12.07  p-value: (0.000513)
```

------------------------------------------------------------------------

Only change the **'type'** argument in `spsurml` function it is possible to estimate several spatial model

**SUR-SAR:** Spatial autorregresive model:
(*y*<sub>*t*</sub> = *λ*<sub>*t*</sub>*W**y*<sub>*t*</sub> + *X*<sub>*t*</sub>*β*<sub>*t*</sub> + *ϵ*<sub>*t*</sub>;  *t* = 1, ..., *T* )

> *W**A**G**E*<sub>83</sub> = *λ*<sub>83</sub>*W**W**A**G**E*<sub>83</sub> + *β*<sub>10</sub> + *β*<sub>11</sub>*U**N*<sub>83</sub> + *β*<sub>12</sub>*N**M**R*<sub>83</sub> + *β*<sub>13</sub>*S**M**S**A* + *ϵ*<sub>83</sub>
> *W**A**G**E*<sub>81</sub> = *λ*<sub>81</sub>*W**W**A**G**E*<sub>81</sub> + *β*<sub>20</sub> + *β*<sub>21</sub>*U**N*<sub>80</sub> + *β*<sub>22</sub>*N**M**R*<sub>80</sub> + *β*<sub>23</sub>*S**M**S**A* + *ϵ*<sub>81</sub>
> *C**o**r**r*(*ϵ*<sub>83</sub>, *ϵ*<sub>81</sub>)≠0

``` r
spcSUR.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc)
#> Initial point:   log_lik:  113.197  lambdas:  -0.472 -0.446 
#> Iteration:  1   log_lik:  114.085  lambdas:  -0.506 -0.482 
#> Iteration:  2   log_lik:  114.096  lambdas:  -0.506 -0.482 
#> Time to fit the model:  3.09  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.36  seconds
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

------------------------------------------------------------------------

**SUR-SEM:** Spatial error model:
(*y*<sub>*t*</sub> = *X*<sub>*t*</sub> *β*<sub>*t*</sub> + *u*<sub>*t*</sub> ; *u*<sub>*t*</sub> = *ρ**u*<sub>*t*</sub> + *ϵ*<sub>*t*</sub> *t* = 1, ..., *T*)

> *W**A**G**E*<sub>83</sub> = *β*<sub>10</sub> + *β*<sub>11</sub>*U**N*<sub>83</sub> + *β*<sub>12</sub>*N**M**R*<sub>83</sub> + *β*<sub>13</sub>*S**M**S**A* + *u*<sub>83</sub>;  *u*<sub>83</sub> = *ρ**W* *u*<sub>83</sub> + *ϵ*<sub>83</sub>
> *W**A**G**E*<sub>81</sub> = *β*<sub>20</sub> + *β*<sub>21</sub>*U**N*<sub>80</sub> + *β*<sub>22</sub>*N**M**R*<sub>80</sub> + *β*<sub>23</sub>*S**M**S**A* + *u*<sub>81</sub>;  *u*<sub>81</sub> = *ρ**W* *u*<sub>81</sub> + *ϵ*<sub>81</sub>
> *C**o**r**r*(*ϵ*<sub>83</sub>, *ϵ*<sub>81</sub>)≠0

``` r
spcSUR.sem <-spsurml(Form=Tformula,data=spc,type="sem",W=Wspc)
#> Initial point:   log_lik:  112.821  deltas:  -0.556 -0.477 
#> Iteration:  1  log_lik:  113.695  rhos:  -0.618 -0.537 
#> Iteration:  2  log_lik:  113.719  rhos:  -0.628 -0.548 
#> Time to fit the model:  2.74  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.33  seconds
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

------------------------------------------------------------------------

**SUR-SARAR:** Spatial autoregressive model with spatial autoregressive error term:
(*y*<sub>*t*</sub> = *λ*<sub>*t*</sub>*W**y*<sub>*t*</sub> + *X*<sub>*t*</sub> *β*<sub>*t*</sub> + *u*<sub>*t*</sub> ; *u*<sub>*t*</sub> = *ρ*<sub>*t*</sub>*u*<sub>*t*</sub> + *ϵ*<sub>*t*</sub> *t* = 1, ..., *T*)

``` r
spcSUR.sarar <-spsurml(Form=Tformula,data=spc,type="sarar",W=Wspc)
#> Initial point:   log_lik:  113.406  lambdas:  -0.343 -0.526  rhos:  -0.28 0.113 
#> Iteration:  1  log_lik:  114.574  lambdas:  -0.384 -0.783  rhos:  -0.307 0.411 
#> Iteration:  2  log_lik:  114.824  lambdas:  -0.394 -0.905  rhos:  -0.303 0.553 
#> Iteration:  3  log_lik:  115.019  lambdas:  -0.395 -0.995  rhos:  -0.296 0.657 
#> Iteration:  4  log_lik:  115.187  lambdas:  -0.388 -1  rhos:  -0.271 0.686 
#> Iteration:  5  log_lik:  115.237  lambdas:  -0.381 -1  rhos:  -0.254 0.699 
#> Iteration:  6  log_lik:  115.26  lambdas:  -0.376 -1  rhos:  -0.243 0.707 
#> Time to fit the model:  45.42  seconds
#> Time to compute covariances:  0.28  seconds
summary(spcSUR.sarar)
#> Call:
#> spsurml(Form = Tformula, data = spc, W = Wspc, type = "sarar")
#> 
#>  
#> Spatial SUR model type:  sarar 
#> 
#> Equation  1 
#>                Estimate Std. Error t value Pr(>|t|)   
#> (Intercept)_1  1.366431   0.506234  2.6992  0.01063 * 
#> UN83_1         0.776669   0.250638  3.0988  0.00382 **
#> NMR83_1       -0.477594   0.243241 -1.9635  0.05758 . 
#> SMSA_1        -0.011898   0.011087 -1.0731  0.29058   
#> lambda_1      -0.375539   0.498019 -0.7541  0.45585   
#> rho_1         -0.242517   0.568388 -0.4267  0.67223   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6488 
#>   Equation  2 
#>                Estimate Std. Error t value  Pr(>|t|)    
#> (Intercept)_2  2.261235   0.298670  7.5710  7.11e-09 ***
#> UN80_2        -0.600803   0.329069 -1.8258 0.0764273 .  
#> NMR80_2        0.422431   0.280809  1.5043 0.1414666    
#> SMSA_2         0.020267   0.026378  0.7683 0.4474484    
#> lambda_2      -1.000000   0.275726 -3.6268 0.0009050 ***
#> rho_2          0.707014   0.173853  4.0667 0.0002573 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> R-squared: 0.6785 
#>   Variance-Covariance Matrix of inter-equation residuals:                            
#>   0.0003134325 -0.0003735885
#>  -0.0003735885  0.0012071470
#> Correlation Matrix of inter-equation residuals:                      
#>   1.0000000 -0.6073533
#>  -0.6073533  1.0000000
#> 
#>  R-sq. pooled: 0.7421 
#>  Log-Likelihood:  115.26
#>  Breusch-Pagan: 9.432  p-value: (0.00213)
```

------------------------------------------------------------------------

Example
-------

This is a basic example which shows you how test spatial structure in the residual of a SUR model:

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
#> Time to fit the model:  2.98  seconds
#> Computing marginal test... 
#> Time to compute covariances:  0.24  seconds
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
