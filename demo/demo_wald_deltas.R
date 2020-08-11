### demo with the whole set of examples of wald_deltas() ########


#################################################
######## CROSS SECTION DATA (G>1; Tm=1) ########
#################################################
rm(list = ls()) # Clean memory
data(spc, package = "spsur")
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA

#################################
## Estimate SUR-SLM model
spcsur.slm <- spsurml(formula = Tformula, data = spc, 
                       type = "slm", listw = lwspc)
summary(spcsur.slm)
## H_0: equality of the lambda parameters of both equations.
R1 <- matrix(c(1,-1), nrow = 1)
b1 <- matrix(0, ncol = 1)
wald_deltas(spcsur.slm, R = R1, b = b1)

#################################
## Estimate SUR-SEM model
spcsur.sem <- spsurml(form = Tformula, data = spc, 
                     type = "sem", listw = lwspc)
summary(spcsur.sem)
## H_0: equality of the rho parameters of both equations.
R2 <- matrix(c(1,-1), nrow=1)
b2 <- matrix(0, ncol=1)
wald_deltas(spcsur.sem, R = R2, b = b2)

#################################
## Estimate SUR-SARAR model
## It usually requires 2-3 minutes maximum
spcsur.sarar <-spsurml(formula = Tformula, data = spc,
                       type = "sarar", listw = lwspc,
                       control = list(tol=0.1))
summary(spcsur.sarar)
## H_0: equality of the lambda and rho parameters of both equations.
R3 <- matrix(c(1,-1,0,0,0,0,1,-1), nrow=2, ncol=4, byrow=TRUE)
b3 <- matrix(c(0,0), ncol=1)
wald_deltas(spcsur.sarar, R = R3, b = b3)

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
## H_0: equality of the lambda parameters of both equations.
R1 <- matrix(c(1,-1), nrow=1)
b1 <- matrix(0, ncol=1)
wald_deltas( NCOVRSUR.slm, R = R1, b = b1)

#################################
## Estimate SUR-SEM model
NCOVRSUR.sem <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "sem", listw = lwncovr,
                       method = "Matrix", zero.policy = TRUE, 
                       control = list(fdHess = TRUE))
summary(NCOVRSUR.sem)
## H_0: equality of the rho parameters of both equations.
R2 <- matrix(c(1,-1), nrow=1)
b2 <- matrix(0, ncol=1)
wald_deltas(NCOVRSUR.sem, R = R2, b = b2)
