### demo with the whole set of examples of wald_betas() ########

#################################################
######## CROSS SECTION DATA (G=1; Tm>1) ########
#################################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data(spc)
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
## Estimate SUR-SLM model
spcsur.slm <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc)
summary(spcsur.slm)
## H_0: equality between SMSA coefficients in both equations.
R1 <- matrix(c(0,0,0,1,0,0,0,-1), nrow=1)
b1 <- matrix(0, ncol=1)

wald_betas(spcsur.slm, R = R1, b = b1)
## Estimate restricted SUR-SLM model
spcsur.slmr <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc,
                      R = R1, b = b1)
summary(spcsur.slmr)

## H_0: equality between intercepts and SMSA coefficients in both equations.
R2 <- matrix(c(1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,-1),
             nrow = 2, ncol = 8, byrow = TRUE)
b2 <- matrix(c(0,0),ncol=1)
wald_betas(spcsur.slm, R = R2, b = b2)
## Estimate restricted SUR-SLM model
spcsur.slmr2 <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc,
                      R = R2, b = b2)
summary(spcsur.slmr2)

####################################
########  G=1; Tm>1         ########
####################################

#### Example 2: Homicides + Socio-Economics (1960-90)

rm(list = ls()) # Clean memory
## Read NCOVR.sf object
data(NCOVR, package = "spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
#################################
## A SUR-SLM model
NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "slm", listw = lwncovr,
                       method = "Matrix", zero.policy = TRUE, 
                       control = list(fdHess = TRUE))
summary(NCOVRSUR.slm)
R1 <- matrix(c(0,1,0,0,-1,0), nrow=1)
b1 <- matrix(0, ncol=1)
wald_betas(NCOVRSUR.slm, R = R1, b = b1)