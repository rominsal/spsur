

# library(spsur)
# library(spdep); library(sf); library(spatialreg)
# #data(boston, package = "spData")
# #lw <- spdep::nb2listw(boston.soi)
data(house, package = "spData")
lw <- spdep::nb2listw(LO_nb, style = "W")
Tformula <- formula(log(price)  ~ TLA + beds)
house[c(1:3,5:7,20:22),c("price")] <- NA

sim_ho <- spsurml(formula = Tformula, data = house, type = "sim" )
summary(sim_ho)

# Check with lm
lm_ho <- lm(formula = Tformula, data = house)
summary(lm_ho)

slx_ho <- spsurml(formula = Tformula, data = house, listw = lw, 
                  zero.policy = TRUE,
                  type = "slx", method = "Matrix", covtype = "num")
summary(slx_ho)

## HACER UN BENCHMARK DEL TIEMPO. 
## TARDA MUCHO MÁS QUE LA FUNCIÓN SPATIALREG::LAGSARLM
slm_ho <- spsurml(formula = Tformula, data = house, listw = lw, 
                  zero.policy = TRUE,
                  type = "slm", method = "Matrix", covtype = "num")
summary(slm_ho)


lagsarlm_slm_h0 <- spatialreg::lagsarlm(formula = Tformula, data = house, 
                              zero.policy = TRUE,          
                              listw = lw, method = "Matrix", type = "lag")
summary(lagsarlm_slm_h0)

## COMPARACIÓN CON SPSE::SPSEML

spseml_slm_h0 <- spse::spseml(list(Tformula), data = house, panel = FALSE,
                              w = lw, 
                              method = "Matrix",
                              model = "lag", zero.policy = TRUE)
summary(spseml_slm_h0)


sdm_ho <- spsurml(formula = Tformula, data = house, listw = lw, 
                  zero.policy = TRUE, type = "sdm", 
                  method = "Matrix", covtype = "num")
summary(sdm_ho)

lagsarlm_sdm_h0 <- spatialreg::lagsarlm(formula = Tformula, data = house, 
                                        listw = lw, zero.policy = TRUE,
                                        method = "Matrix", 
                                        type = "Durbin")
summary(lagsarlm_sdm_h0)


sem_ho <- spsurml(formula = Tformula, data = house, listw = lw, 
                  type = "sem", method = "Matrix", covtype = "num")
summary(sem_ho)

lagsarlm_sem_h0 <- spatialreg::errorsarlm(formula = Tformula, data = house, 
                                           listw = lw, method = "Matrix")
summary(lagsarlm_sem_h0)

## TARDA MUCHO
#sarar_ho <- spsurml(formula = Tformula, data = house, listw = lw, 
#                  type = "sarar", method = "Matrix", covtype = "num")
#summary(sarar_ho)

 
lagsarlm_sarar_h0 <- spatialreg::sacsarlm(formula = Tformula, data = house, 
                                          listw = lw, method = "MC") 
### VIP: NO PUEDE CALCULARLO CON method = "Matrix"
summary(lagsarlm_sarar_h0)


## TARDA MUCHO
#sarar_ho_MC <- spsurml(formula = Tformula, data = house, listw = lw, 
#                    type = "sarar", method = "MC", covtype = "num")
#summary(sarar_ho_MC)


############################## EXAMPLE WITH NAT FILE ##################
ncovr <- sf::st_read("C:/Users/Roman.Minguez/OneDrive/spsurdev/notes/ncovr/NAT.shp") 


attr(ncovr, "sf_column")
(ncovr_geom <- sf::st_geometry(ncovr))
ncovr_geom[[1]]
par(mar = c(0,0,1,0))
#plot(ncovr[1,1], col = 'grey')
# PLOT OF USA
plot(ncovr[1], reset = FALSE) 

# some of the polygons in this dataset have multiple exterior rings; 
# they can be identified by
par(mar = c(0,0,1,0))
(w <- which(sapply(ncovr_geom, length) > 1))
plot(ncovr[w,1], col = 2:7)
attributes(ncovr_geom)

# neighbours
ncovr_nb <- spdep::poly2nb(ncovr, queen = TRUE)
#class(ncovr_nb)
summary(ncovr_nb)
spdep::is.symmetric.nb(ncovr_nb)
ncovr_lw <- spdep::nb2listw(ncovr_nb, style = "W", zero.policy = TRUE)

Tformula <- HR80  | HR90 ~ PS80 + UE80 | PS90 + UE90
Tformulalist <- list(eq1 = HR80 ~ PS80 + UE80, 
                     eq2 = HR90 ~ PS90 + UE90)
NCOVRSUR.sim <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "sim")
summary(NCOVRSUR.sim)

NCOVRSUR.slx <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "slx")
summary(NCOVRSUR.slx)
ncovrW <- as(ncovr_lw, "CsparseMatrix")
tr <- spatialreg::trW(ncovrW, type = "MC")
NCOVRSUR.slx.impacts <- impacts.spsur(NCOVRSUR.slx,
                                      tr = tr, R = 1000)
class(NCOVRSUR.slx.impacts[[1]])
summary(NCOVRSUR.slx.impacts[[1]])

NCOVRSUR.slm <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "slm", 
                        con = list(fdHess = FALSE))
summary(NCOVRSUR.slm)
ncovrW <- as(ncovr_lw, "CsparseMatrix")
tr <- spatialreg::trW(ncovrW, type = "MC")
NCOVRSUR.slm.impacts <- impacts(NCOVRSUR.slm, 
                                tr = tr, R = 1000)
summary(NCOVRSUR.slm.impacts[[1]], zstats = TRUE)
summary(NCOVRSUR.slm.impacts[[2]], zstats = TRUE)

## COMPARATIVA CON SPSEML
ncov_spseml.slm <- spse::spseml(Tformulalist, data = ncovr, panel = FALSE,
                                w = ncovr_lw, method = "Matrix",
                                model = "lag")

NCOVRSUR.sdm <- spsurml(formula = Tformula, data = ncovr, 
                        listw = ncovr_lw, zero.policy = TRUE,
                        method = "Matrix", type = "sdm", 
                        con = list(fdHess = TRUE))
summary(NCOVRSUR.sdm)
NCOVRSUR.sdm.impacts <- impacts(NCOVRSUR.sdm, tr = tr, R = 1000)
summary(NCOVRSUR.sdm.impacts[[1]], zstats = TRUE)
summary(NCOVRSUR.sdm.impacts[[2]], zstats = TRUE)

NCOVRSUR.sem <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "sem", covtype = "num")
summary(NCOVRSUR.sem)


NCOVRSUR.sdem <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "sdem", covtype = "num")
summary(NCOVRSUR.sdem)

NCOVRSUR.sarar <- spsurml(formula = Tformula, data = ncovr, listw = ncovr_lw, 
                        method = "Matrix", type = "sarar", 
                        con = list(fdHess = TRUE))
summary(NCOVRSUR.sarar)

NCOVRSUR.sarar.impacts <- impacts.spsur(NCOVRSUR.sarar, tr = tr, R = 1000)
summary(NCOVRSUR.sarar.impacts[[1]], zstats = TRUE)
summary(NCOVRSUR.sarar.impacts[[2]], zstats = TRUE)


## Check with W instead of listw argument...
ncovr_W <- Matrix::Matrix(spdep::listw2mat(ncovr_lw))
NCOVRSUR.slm2 <- spsurml(formula = Tformula, data = ncovr, W = ncovr_W, 
                        method = "LU", type = "slm", covtype = "num")
summary(NCOVRSUR.slm2)
## Not allowed with method == c("Matrix","Matrix_J","MC")
NCOVRSUR.slm3 <- spsurml(formula = Tformula, data = ncovr, W = ncovr_W, 
                         method = "Matrix", type = "slm", covtype = "num")


## Check without W instead of listw argument... ONLY SIM
NCOVRSUR.sim2 <- spsurml(formula = Tformula, data = ncovr, 
                         type = "sim", covtype = "num", N=nrow(ncovr))
summary(NCOVRSUR.sim2)
# Don't should work with type != ""sim"
NCOVRSUR.slm4 <- spsurml(formula = Tformula, data = ncovr, 
                         type = "slm", covtype = "num", N=nrow(ncovr))




# ## COMPROBACIÓN CÁLCULOS DETERMINANTE CON DISTINTOS MÉTODOS...
# env <- new.env()
# #assign("listw", listw, envir = env)
# assign("listw", lw, envir = env)
# assign("n", length(lw$neighbours), envir = env)
# assign("similar", FALSE, envir = env)
# assign("can.sim", spatialreg::can.be.simmed(lw), envir = env)
# quiet <- TRUE
# assign("verbose", !quiet, envir = env)
# assign("family", "SAR", envir = env) # CHEQUEAR OTRAS OPCIONES
# ## method <- c("eigen",  # NO USAR EIGEN PARA MUCHOS VECINOS...
# method <-c("Matrix", "Matrix_J", "spam",
#             "spam_update", "LU", "MC", "Chebyshev")
# det_jac <- vector("list",length(method))
# names(det_jac) <- method
# for (i in 1:length(method)) {
#   meth_i <- method[i]
#   print(meth_i)
#   #interval <-
#   ## OJO: CAMBIA EN ENVIRONMENT AL HACER jacobianSetup...
#   spatialreg::jacobianSetup(meth_i, env, con, pre_eig = con$pre_eig,
#                                         tr = trs)
#   #print(interval)
#   #lambda <- seq(interval[1]+0.001, interval[2]-0.001, length.out = 1000)
#   lambda <- seq(-0.999, 0.999, length.out = 1000)
#   print(system.time( for(j in 1:length(lambda)){
#     det_jac[[i]][j] <- spatialreg::do_ldet(lambda[j], env = env)
#   }
#   ))
# 
# }
# 
# # all.equal(det_jac$eigen, det_jac$Matrix)
# # all.equal(det_jac$eigen, det_jac$Matrix_J)
# # all.equal(det_jac$eigen, det_jac$spam, det_jac$spam_update)
# # all.equal(det_jac$eigen, det_jac$MC)
# # all.equal(det_jac$eigen, det_jac$Chebyshev)
# 
# all.equal(det_jac$Matrix, det_jac$Matrix_J)
# all.equal(det_jac$Matrix, det_jac$spam, det_jac$spam_update)
# all.equal(det_jac$Matrix, det_jac$MC)
# all.equal(det_jac$Matrix, det_jac$Chebyshev)
# 
# # LA MATRIZ DEL EJEMPLO SPC EN SPSURML NO PERMITE HACERSE SIMILAR...
