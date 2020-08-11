### demo with the whole set of examples of impactspsur() ########

###############################################
### PURE CROSS SECTIONAL DATA(G>1; Tm=1) ######
###############################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data(spc)
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
## For SLM, SDM and SARAR models the output is a list of "lagImpact" objects
## See spatialreg::impacts for details.
spcsur_slm <-spsurml(formula = Tformula, data = spc, 
                     type = "slm", listw = lwspc)
summary(spcsur_slm)
impacts_slm <- impactspsur(spcsur_slm, listw = lwspc, R = 1000)
## Impacts equation 1
summary(impacts_slm[[1]], zstats = TRUE, short = TRUE)
## Impacts equation 2
summary(impacts_slm[[2]], zstats = TRUE, short = TRUE)
## For SLX and SDEM models the output is a list of "WXImpact" objects
## See spatialreg::impacts for details.
## A SUR-SLX model
spcsur_slx <-spsurml(formula = Tformula, data = spc, 
                     type = "slx", listw = lwspc)
summary(spcsur_slx)
impacts_slx <- impactspsur(spcsur_slx, listw = lwspc)
summary(impacts_slx[[1]], zstats = TRUE, short = TRUE)
summary(impacts_slx[[2]], zstats = TRUE, short = TRUE)

## A SUR-SDM model
spcsur_sdm <-spsurml(formula = Tformula, data = spc, 
                     type = "sdm", listw = lwspc)
impacts_sdm <- impactspsur(spcsur_sdm, listw = lwspc, R = 1000)
summary(impacts_sdm[[1]], zstats = TRUE, short = TRUE)
summary(impacts_sdm[[2]], zstats = TRUE, short = TRUE)
# A SUR-SDM model with different spatial lags in each equation
TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
spcsur_sdm2 <-spsurml(formula = Tformula, data = spc, type = "sdm", 
                      listw = lwspc, Durbin = TformulaD)
summary(spcsur_sdm2)                       
impacts_sdm2 <- impactspsur(spcsur_sdm2, listw = lwspc, R = 1000)
summary(impacts_sdm2[[1]], zstats = TRUE, short = TRUE)
summary(impacts_sdm2[[2]], zstats = TRUE, short = TRUE)
# A SUR-SLX model with different spatial lags in each equation
spcsur_slx2 <-spsurml(formula = Tformula, data = spc, 
                     type = "slx", listw = lwspc, Durbin = TformulaD)
summary(spcsur_slx2)
impacts_slx2 <- impactspsur(spcsur_slx2, listw = lwspc)
summary(impacts_slx2[[1]], zstats = TRUE, short = TRUE)
summary(impacts_slx2[[2]], zstats = TRUE, short = TRUE)

## A SUR-SDEM model
spcsur_sdem <-spsurml(formula = Tformula, data = spc, 
                     type = "sdem", listw = lwspc)
impacts_sdem <- impactspsur(spcsur_sdem, listw = lwspc)
summary(impacts_sdem[[1]], zstats = TRUE, short = TRUE)
summary(impacts_sdem[[2]], zstats = TRUE, short = TRUE)

## A SUR-SARAR model
spcsur_sarar <-spsurml(formula = Tformula, data = spc, 
                     type = "sarar", listw = lwspc,
                     control = list(tol = 0.01))
impacts_sarar <- impactspsur(spcsur_sarar, listw = lwspc, R = 1000)
summary(impacts_sarar[[1]], zstats = TRUE, short = TRUE)
summary(impacts_sarar[[2]], zstats = TRUE, short = TRUE)

####################################
######## G=1; Tm>1               ###
####################################

rm(list = ls()) # Clean memory
data(NCOVR, package="spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
## A SUR-SLM model
NCOVRSUR_slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "slm", listw = lwncovr, 
                       method = "Matrix", zero.policy = TRUE, 
                       control = list(fdHess = TRUE))
summary(NCOVRSUR_slm)
## Use of trW to compute.
Wncovr <- as(spdep::listw2mat(lwncovr), "CsparseMatrix")
trwncovr <- spatialreg::trW(Wncovr, type = "MC")
impacts_NCOVRSUR_slm <- impactspsur(NCOVRSUR_slm, tr = trwncovr,
                                R = 1000)
summary(impacts_NCOVRSUR_slm[[1]], zstats = TRUE, short = TRUE)
summary(impacts_NCOVRSUR_slm[[2]], zstats = TRUE, short = TRUE)
## A SUR-SDM model
NCOVRSUR_sdm <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "sdm", listw = lwncovr, 
                       method = "Matrix", zero.policy = TRUE, 
                       control = list(fdHess = TRUE))
impacts_NCOVRSUR_sdm <- impactspsur(NCOVRSUR_sdm, tr = trwncovr,
                                R = 1000)
summary(impacts_NCOVRSUR_sdm[[1]], zstats = TRUE, short = TRUE)
summary(impacts_NCOVRSUR_sdm[[2]], zstats = TRUE, short = TRUE)
