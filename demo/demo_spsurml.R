### demo with the whole set of examples of spsurml() ########

#################################################
######## CROSS SECTION DATA (G>1; Tm=1) ########
#################################################
#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data(spc, package = "spsur")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
spcsur.sim <- spsurml(formula = Tformula, data = spc, type = "sim")
summary(spcsur.sim)
## A SUR-SLX model 
## (listw argument can be either a matrix or a listw object )
spcsur.slx <-spsurml(formula = Tformula, data = spc, type = "slx", 
  listw = Wspc)
summary(spcsur.slx)

## A SUR-SLM model
spcsur.slm <-spsurml(formula = Tformula, data = spc, type = "slm", 
  listw = Wspc)
summary(spcsur.slm)

## A SUR-SEM model
spcsur.sem <-spsurml(formula = Tformula, data = spc, type = "sem", 
  listw = Wspc)
summary(spcsur.sem)
rm(spcsur.sem) # remove

## A SUR-SDM model
spcsur.sdm <-spsurml(formula = Tformula, data = spc, type = "sdm", 
  listw = Wspc)
summary(spcsur.sdm)
rm(spcsur.sdm) # remove

# A SUR-SDM model with different spatial lags in each equation
TformulaD <- ~ UN83 + NMR83 + SMSA | UN80
spcsur.sdm2 <-spsurml(formula = Tformula, data = spc, type = "sdm", 
                      listw = Wspc, Durbin = TformulaD)
summary(spcsur.sdm2)
rm(spcsur.sdm2) # remove

## A SUR-SDEM model
spcsur.sdem <-spsurml(formula = Tformula, data = spc, type = "sdem", 
  listw = Wspc)
summary(spcsur.sdem)
rm(spcsur.sdem) # remove

## A SUR-SARAR model
spcsur.sarar <-spsurml(formula = Tformula, data = spc, type = "sarar", 
  listw = Wspc, control = list(tol = 0.1))
summary(spcsur.sarar)
rm(spcsur.sarar) # remove

#################################################
########  G=1; Tm>1         ########
#################################################

#### Example 2: Homicides + Socio-Economics (1960-90)
# Homicides and selected socio-economic characteristics for continental
# U.S. counties.
# Data for four decennial census years: 1960, 1970, 1980 and 1990.
# \url{https://geodacenter.github.io/data-and-lab/ncovr/}

## It usually requires 1-2 minutes maximum...
rm(list = ls()) # Clean memory
## Read NCOVR.sf object
data(NCOVR, package = "spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90

## A SUR-SIM model
NCOVRSUR.sim <-spsurml(formula = Tformula, data = NCOVR.sf, type = "sim")
summary(NCOVRSUR.sim)
rm(NCOVRSUR.sim)

## A SUR-SLX model
NCOVRSUR.slx <-spsurml(formula = Tformula, data = NCOVR.sf, type = "slx", 
  listw = lwncovr, zero.policy = TRUE)
summary(NCOVRSUR.slx)
rm(NCOVRSUR.slx)

## A SUR-SLM model
## method = "Matrix" (Cholesky) instead of "eigen"
## (fdHess = TRUE to compute numerical covariances )
NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
  type = "slm", listw = lwncovr, method = "Matrix", 
  zero.policy = TRUE, control = list(fdHess = TRUE))
summary(NCOVRSUR.slm)
rm(NCOVRSUR.slm)

## A SUR-SDM model with different spatial lags in each equation
## Analytical covariances (default)
TformulaD <- ~ PS80 + UE80 | PS90 
NCOVRSUR.sdm <-spsurml(formula = Tformula, data = NCOVR.sf, 
  type = "sdm", listw = lwncovr, method = "Matrix",
  Durbin = TformulaD, zero.policy = TRUE)
summary(NCOVRSUR.sdm)
rm(NCOVRSUR.sdm)

## A SUR-SEM model
NCOVRSUR.sem <- spsurml(formula = Tformula, data = NCOVR.sf, 
  type = "sem", listw = lwncovr, method = "Matrix",
  zero.policy = TRUE, control = list(fdHess = TRUE))
summary(NCOVRSUR.sem)

## A SUR-SDEM model
NCOVRSUR.sdem <-spsurml(formula = Tformula, data = NCOVR.sf, 
  type = "sdem", listw = lwncovr, method = "Matrix",
  zero.policy = TRUE, control = list(fdHess = TRUE))
summary(NCOVRSUR.sdem)
rm(NCOVRSUR.sdem)

## A SUR-SARAR model
NCOVRSUR.sarar <-spsurml(formula = Tformula, data = NCOVR.sf,
  type = "sarar", listw = lwncovr, method = "Matrix",
  zero.policy = TRUE, control = list(fdHess = TRUE))
summary(NCOVRSUR.sarar)
rm(NCOVRSUR.sarar)

##############################################
######## SUR with G>1; Tm>1  ###
##############################################
#### Reshape NCOVR in panel format
N <- nrow(NCOVR.sf)
Tm <- 4
index_time <- rep(1:Tm, each = N)
index_indiv <- rep(1:N, Tm)
pHR <- c(NCOVR.sf$HR60, NCOVR.sf$HR70, NCOVR.sf$HR80, NCOVR.sf$HR90)
pPS <- c(NCOVR.sf$PS60, NCOVR.sf$PS70, NCOVR.sf$PS80, NCOVR.sf$PS90)
pUE <- c(NCOVR.sf$UE60, NCOVR.sf$UE70, NCOVR.sf$UE80, NCOVR.sf$UE90)
pDV <- c(NCOVR.sf$DV60, NCOVR.sf$DV70, NCOVR.sf$DV80, NCOVR.sf$DV90)
pFP <- c(NCOVR.sf$FP59, NCOVR.sf$FP70, NCOVR.sf$FP80, NCOVR.sf$FP90)
pSOUTH <- rep(NCOVR.sf$SOUTH, Tm)
pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
                     HR = pHR, PS = pPS, UE = pUE, DV = pDV,
                     FP = pFP, SOUTH = pSOUTH)
pform <- HR | DV | FP ~ PS + UE | PS + UE + SOUTH | PS

## SIM 
## Remark: It is necessary to provide Tm value as argument 
## when G>1 && Tm>1
pNCOVRSUR.sim <- spsurml(formula = pform, data = pNCOVR, 
  type = "sim", Tm = Tm)
summary(pNCOVRSUR.sim)

pNCOVRSUR.slm <- spsurml(formula = pform, data = pNCOVR, 
  listw = lwncovr, type = "slm", method = "Matrix", Tm = Tm,
  zero.policy = TRUE, control= list(fdHess = TRUE))
summary(pNCOVRSUR.slm)   

pNCOVRSUR.sem <- spsurml(formula = pform, data = pNCOVR, 
  listw = lwncovr, type = "sem", method = "Matrix", Tm = Tm,
  zero.policy = TRUE, control= list(fdHess = TRUE))
summary(pNCOVRSUR.sem)  
