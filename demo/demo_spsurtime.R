### demo with the whole set of examples of spsurtime() ########

####################################
######## PANEL DATA (G=1; Tm>1) ###
####################################

## Example 1:
rm(list = ls()) # Clean memory
data(spc)
lwspc <- spdep::mat2listw(Wspc, style = "W")
N <- nrow(spc)
Tm <- 2
index_time <- rep(1:Tm, each = N)
index_indiv <- rep(1:N, Tm)
WAGE <- c(spc$WAGE83, spc$WAGE81)
UN <- c(spc$UN83, spc$UN80)
NMR <- c(spc$NMR83, spc$NMR80)
SMSA <- c(spc$SMSA, spc$SMSA)
pspc <- data.frame(index_indiv, index_time, WAGE, UN,
                    NMR, SMSA)
form_pspc <- WAGE ~ UN + NMR + SMSA

# SLM 
pspc_slm <- spsurtime(formula = form_pspc, data = pspc, 
                      listw =lwspc,
                      time = pspc$index_time, 
                      type = "slm", fit_method = "ml")
summary(pspc_slm)
## Example 2:
rm(list = ls()) # Clean memory
## Read NCOVR.sf object
data(NCOVR, package="spsur")
nbncovr <- spdep::poly2nb(NCOVR.sf, queen = TRUE)
## Some regions with no links...
lwncovr <- spdep::nb2listw(nbncovr, style = "W", zero.policy = TRUE)
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

# SIM

pHR_sim <- spsurtime(formula = form_pHR, data = pNCOVR, 
                    time = pNCOVR$time, type = "sim", fit_method = "ml")
summary(pHR_sim)

# SLM by 3SLS. 

pHR_slm <- spsurtime(formula = form_pHR, data = pNCOVR, listw = lwncovr,
                     time = pNCOVR$time, type = "slm", 
                     fit_method = "3sls")
summary(pHR_slm)

############################ Wald tests about betas in spatio-temporal models
## H0: equal betas for PS in equations 1, 3 and 4.
R <- matrix(0, nrow = 2, ncol = 12) 
# nrow = number of restrictions 
# ncol = number of beta parameters
R[1, 2] <- 1; R[1, 8] <- -1 # PS beta coefficient in equations 1 equal to 3
R[2, 2] <- 1; R[2, 11] <- -1 # PS beta coefficient in equations 1 equal to 4
b <- matrix(0, nrow = 2, ncol = 1)
wald_betas(pHR_sim , R = R , b = b) # SIM model
wald_betas(pHR_slm , R = R , b = b) # SLM model

############################ Wald tests about spatial-parameters in
############################ spatio-temporal models
## H0: equal rhos in slm model for equations 1 and 2.
R2 <- matrix(0, nrow = 1, ncol = 4)
R2[1, 1] <- 1; R2[1, 2] <- -1
b2 <- matrix(0, nrow = 1, ncol = 1)
wald_deltas(pHR_slm, R = R2, b = b2)
