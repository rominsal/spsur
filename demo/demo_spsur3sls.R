### demo with the whole set of examples of spsur3sls() ########

#################################################
######## CROSS SECTION DATA (G=1; Tm>1) ########
#################################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
## A SUR model without spatial effects
rm(list = ls()) # Clean memory
data(spc)
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA

## A SUR-SLM model (3SLS Estimation)
spcsur_slm_3sls <-spsur3sls(formula = Tformula, data = spc,
                            type = "slm", listw = lwspc)
summary(spcsur_slm_3sls)

## A SUR-SDM model (3SLS Estimation)
spcsur_sdm_3sls <-spsur3sls(formula = Tformula, data = spc,
                            type = "sdm", listw = lwspc)
summary(spcsur_sdm_3sls)

# A SUR-SDM model with different spatial lags in each equation
TformulaD <-  ~ UN83 + NMR83 + SMSA | UN80 + NMR80  
spcsur_sdm2_3sls <-spsur3sls(formula = Tformula, data = spc,
                            type = "sdm", listw = lwspc,
                            Durbin = TformulaD)
summary(spcsur_sdm2_3sls)


#################################################
######## PANEL DATA (G>1; Tm>1)         #########
#################################################
#'
#### Example 2: Homicides + Socio-Economics (1960-90)
# Homicides and selected socio-economic characteristics for continental
# U.S. counties.
# Data for four decennial census years: 1960, 1970, 1980 and 1990.
# https://geodacenter.github.io/data-and-lab/ncovr/
rm(list = ls()) # Clean memory
## Read NCOVR.sf object
data(NCOVR, package = "spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
## A SUR-SLM model
NCOVRSUR_slm_3sls <-spsur3sls(formula = Tformula, data = NCOVR.sf, 
                              type = "slm", zero.policy = TRUE,
                              listw = lwncovr, maxlagW = 2)
summary(NCOVRSUR_slm_3sls)

