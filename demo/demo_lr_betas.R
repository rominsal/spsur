### demo with the whole set of examples of lr_betas() ########

#################################################
######## CROSS SECTION DATA (G>1; Tm=1)  ########
#################################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data(spc)
lwspc <- spdep::mat2listw(Wspc, style = "W")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
## H0: equal beta for SMSA in both equations.
R <- matrix(c(0,0,0,1,0,0,0,-1), nrow=1)
b <- matrix(0, ncol=1)
spcsur.slm <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc)
summary(spcsur.slm)
lr_betas(spcsur.slm, R = R, b = b)

## Estimate restricted SUR-SLM model
spcsur.slmr <- spsurml(formula = Tformula, data = spc, 
                      type = "slm", listw = lwspc,
                      R = R, b = b)
summary(spcsur.slmr)
################################################
####### PANEL DATA (G>1; Tm>1)          ########
################################################

#### Example 2: Homicides + Socio-Economics (1960-90)
# Homicides and selected socio-economic characteristics for continental
# U.S. counties.
# Data for four decennial census years: 1960, 1970, 1980 and 1990.
# \url{https://geodacenter.github.io/data-and-lab/ncovr/}
rm(list = ls()) # Clean memory
data(NCOVR, package="spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
## H0: equal beta for PS80 and PS90 in both equations.
NCOVRSUR.slm <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "slm", listw = lwncovr,
                       method = "LU", zero.policy = TRUE, 
                       control = list(fdHess = TRUE))
summary(NCOVRSUR.slm)
R <- matrix(c(0, 1, 0, 0, -1, 0), nrow = 1)
b <- matrix(0, ncol = 1)
lr_betas(NCOVRSUR.slm, R = R, b = b)
## Restricted model
NCOVRSUR.slmr <-spsurml(formula = Tformula, data = NCOVR.sf, 
                       type = "slm", listw = lwncovr,
                       method = "LU", zero.policy = TRUE, 
                       control = list(fdHess = TRUE),
                       R = R, b = b)
summary(NCOVRSUR.slmr)                        
################################################################
######## PANEL DATA: TEMPORAL CORRELATIONS (nG=1; nT>1) ########
################################################################
## Example 3: with classical panel data set. Database is
##            a spatio-temporal panel

N <- nrow(NCOVR.sf)
Tm <- 4
index_time <- rep(1:Tm, each = N)
index_indiv <- rep(1:N, Tm)
pHR <- c(NCOVR.sf$HR60, NCOVR.sf$HR70, NCOVR.sf$HR80, NCOVR.sf$HR90)
pPS <- c(NCOVR.sf$PS60, NCOVR.sf$PS70, NCOVR.sf$PS80, NCOVR.sf$PS90)
pUE <- c(NCOVR.sf$UE60, NCOVR.sf$UE70, NCOVR.sf$UE80, NCOVR.sf$UE90)
pNCOVR <- data.frame(indiv = index_indiv, time = index_time,
                     HR = pHR, PS = pPS, UE = pUE)
form_pHR <- HR ~ PS + UE
# H0: equal PS beta coefficient in equations 1, 3, and 4
# Fit with spsurtime and fit_method = "ml"...
pHR.slm <-spsurtime(formula = form_pHR, data = pNCOVR,
                    time = pNCOVR$time, 
                    type = "slm", listw = lwncovr,
                    zero.policy = TRUE,
                    fit_method = "ml", method = "LU", 
                    control = list(fdHess = TRUE))
summary(pHR.slm)
## H0: equal betas for PS in equations 1, 3 and 4.
R <- matrix(0, nrow = 2, ncol = 12) 
# nrow = number of restrictions 
# ncol = number of beta parameters
R[1, 2] <- 1; R[1, 8] <- -1 # PS beta coefficient in equations 1 equal to 3
R[2, 2] <- 1; R[2, 11] <- -1 # PS beta coefficient in equations 1 equal to 4
b <- matrix(0, nrow = 2, ncol = 1)
lr_betas(pHR.slm, R = R, b = b)


