### demo with the whole set of examples of lrtestspsur() ########

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
## Test ANOVA SLM vs SDM
spcsur.sdm <- spsurml(formula = Tformula, data = spc, 
                      type = "sdm", listw = lwspc)
lrtestspsur(spcsur.slm, spcsur.sdm)
## Test ANOVA SEM vs SDM
spcsur.sem <- spsurml(formula = Tformula, data = spc, 
                      type = "sem", listw = lwspc)
lrtestspsur(spcsur.sem, spcsur.sdm)

