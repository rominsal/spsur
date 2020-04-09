# LA MATRIZ DEL EJEMPLO SPC EN SPSURML NO PERMITE HACERSE SIMILAR...
# # ARGUMENTOS "A MANO" PARA SPSURML FUNCTION
# data(spc)
# lWspc <- mat2listw(Wspc)
# Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA


# library(spsur)
# library(spdep); library(sf); library(spatialreg)
# #data(boston, package = "spData")
# #lw <- spdep::nb2listw(boston.soi)
# data(house, package = "spData")
# lw <- spdep::nb2listw(LO_nb, style = "W")
# Tformula <- Formula::Formula(log(price) | age ~ TLA + beds | TLA + beds)
# Tformula <- formula(log(price)  ~ TLA + beds)

# pr <- spsurml(Form = Tformula, data = house, listw = lw, 
#                           type = "sim", method = "Matrix",
#                           cov = TRUE, covtype = "num")
#   
# 
# # HAY QUE CARGAR LA FUNCIÓN SPSURML PREVIAMENTE EN EL ENVIRONMENT
# Form <- Tformula
# data <- house
# listw <- lw
# type <- "slm"
# method <- "Matrix"
# quiet <- TRUE
# zero.policy <- NULL; interval <- NULL; trs <- NULL
# 
# slm1eq <- spatialreg::lagsarlm(formula=Form, data = data, listw = lw,
#                            method = method, type = "lag")
# summary(slm1eq)
# 
# 
# sem1eq <- spatialreg::errorsarlm(formula=Form, data = data, listw = lw,
#                                method = method)
# summary(sem1eq)
# 
# sac1eq <- spatialreg::sacsarlm(formula = Form, data = data, listw = lw,
#                                  method = "MC")
# summary(sac1eq)
# 
# # Inicia function spsurml
# 
# 
# 
# 
# 
# 
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
