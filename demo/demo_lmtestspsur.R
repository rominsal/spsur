### demo with the whole set of examples of lmtestspsur() ########


#################################################
######## CROSS SECTION DATA (G>1; Tm=1) # #######
#################################################

#### Example 1: Spatial Phillips-Curve. Anselin (1988, p. 203)
rm(list = ls()) # Clean memory
data("spc")
Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
lwspc <- spdep::mat2listw(Wspc, style = "W")
lmtestspsur(formula = Tformula, data = spc, 
             listw = lwspc)

#################################################
######## PANEL DATA (G>1; Tm>1)          ########
#################################################

#### Example 2: Homicides & Socio-Economics (1960-90)
# Homicides and selected socio-economic characteristics for
# continental U.S. counties.
# Data for four decennial census years: 1960, 1970, 1980 and 1990.
# https://geodacenter.github.io/data-and-lab/ncovr/
data(NCOVR, package="spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
## With different number of exogenous variables in each equation
Tformula <- HR70 | HR80  | HR90 ~ PS70 + UE70 | PS80 + UE80 +RD80 |
            PS90 + UE90 + RD90 + PO90
lmtestspsur(formula = Tformula, data = NCOVR.sf, 
            listw = lwncovr)

################################################################
######## PANEL DATA: TEMPORAL CORRELATIONS (G=1; Tm>1) ########
################################################################

#### Example 3: NCOVR in panel data form
Year <- as.numeric(kronecker(c(1960,1970,1980,1990), 
                   matrix(1,nrow = dim(NCOVR.sf)[1])))
HR <- c(NCOVR.sf$HR60,NCOVR.sf$HR70,NCOVR.sf$HR80,NCOVR.sf$HR90)
PS <- c(NCOVR.sf$PS60,NCOVR.sf$PS70,NCOVR.sf$PS80,NCOVR.sf$PS90)
UE <- c(NCOVR.sf$UE60,NCOVR.sf$UE70,NCOVR.sf$UE80,NCOVR.sf$UE90)
NCOVRpanel <- as.data.frame(cbind(Year,HR,PS,UE))
Tformula <- HR ~ PS + UE
lmtestspsur(formula = Tformula, data = NCOVRpanel, 
            time = Year, listw = lwncovr)

