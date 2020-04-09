#library(spdep)
library(spatialreg)
tab_out <- list()
 dsets <- c("USC", "Walde4900", "wrld", "LO", "USZC", "UST")
 for (dset in dsets) {
     load(paste(dset, "_lw.RData", sep = ""))
     nb <- get(paste(dset, "_lw", sep = ""))$neighbours
     tab_out[[dset]] <- vector(mode = "list", length = 4)
     cnb <- card(nb)
     tab_out[[dset]][[1]] <- c(table(cnb))
     tab_out[[dset]][[2]] <- length(nb)
     tab_out[[dset]][[3]] <- n.comp.nb(nb)$nc
     tab_out[[dset]][[4]] <- sum(cnb)
     }
 save(tab_out, file = "tab_out.RData")

 library(spData)
 data(house)
 class(LO_nb)

 #library(spdep)
 #load("USC_lw.RData")
 #USC_nb <- USC_lw$neighbours
 set.ZeroPolicyOption(TRUE)
 eigs_out <- vector(mode = "list", length = 6)
 eigs_out[[1]] <- eigenw(nb2listw(LO_nb, style = "B"))
 eigs_out[[2]] <- eigenw(nb2listw(LO_nb, style = "C"))
 eigs_out[[3]] <- eigenw(nb2listw(LO_nb, style = "S"))
 eigs_out[[4]] <- eigenw(nb2listw(LO_nb, style = "W"))
 eigs_out[[5]] <- eigenw(similar.listw(nb2listw(LO_nb, style = "S")))
 eigs_out[[6]] <- eigenw(similar.listw(nb2listw(LO_nb, style = "W")))
 eig_res <- sapply(eigs_out, function(x) 1/range(Re(x)))
 save(eig_res, file = "eigs_out_res.RData")
 library(spdep)
 library(spam)
 set.ZeroPolicyOption(TRUE)
 lambda <- seq(-0.9, 0.99, 0.01)
 dsets <- c("USC", "Walde4900", "wrld", "LO", "USZC", "UST")
 output <- list()
 for (dset in dsets) {
     load(paste(dset, "_lw.RData", sep = ""))
     listw <- get(paste(dset, "_lw", sep = ""))
     can.sim <- FALSE
     if (listw$style %in% c("W", "S"))
         can.sim <- spdep:::can.be.simmed(listw)
         res <- list()
         if (length(listw$neighbours) < 5000) {
             env <- new.env(parent = globalenv())
             assign("listw", listw, envir = env)
             assign("verbose", FALSE, envir = env)
             assign("can.sim", can.sim, envir = env)
             assign("family", "SAR", envir = env)
             setTime <- system.time(eigen_setup(env))
             type <- get("method", envir = env)
             res[[type]] <- vector(mode = "list", length = 3)
             res[[type]][[1]] <- setTime
             res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                    function(x) do_ldet(x, env)))
             res[[type]][[3]] <- out
             rm(env)
             }
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("family", "SAR", envir = env)
         set.seed(length(listw$neighbours))
         setTime <- system.time(mcdet_setup(env))
         type <- get("method", envir = env)
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
        rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(cheb_setup(env, q = 2))
         type <- paste(get("method", envir = env), "2", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(cheb_setup(env, q = 5))
         type <- paste(get("method", envir = env), "5", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_setup(env, Imult = 2, super = FALSE))
         type <- paste(get("method", envir = env), "simp", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_setup(env, Imult = 2, super = TRUE))
         type <- paste(get("method", envir = env), "sup", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda, function(x) do_ldet(x,
                                                                                     env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_setup(env, Imult = 2,
                                             super = as.logical(NA)))
         type <- paste(get("method", envir = env), "NA", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_J_setup(env, super = FALSE))
         type <- paste(get("method", envir = env), "simp", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_J_setup(env, super = TRUE))
         type <- paste(get("method", envir = env), "sup", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(Matrix_J_setup(env, super = as.logical(NA)))
         type <- paste(get("method", envir = env), "NA", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(spam_setup(env))
         type <- get("method", envir = env)
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(spam_setup(env, pivot = "RCM"))
         type <- paste(get("method", envir = env), "RCM", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(spam_update_setup(env))
         type <- get("method", envir = env)
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("can.sim", can.sim, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(spam_update_setup(env, pivot = "RCM"))
         type <- paste(get("method", envir = env), "RCM", sep = "_")
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         env <- new.env(parent = globalenv())
         assign("listw", listw, envir = env)
         assign("n", length(listw$neighbours), envir = env)
         assign("family", "SAR", envir = env)
         setTime <- system.time(LU_setup(env))
         type <- get("method", envir = env)
         res[[type]] <- vector(mode = "list", length = 3)
         res[[type]][[1]] <- setTime
         res[[type]][[2]] <- system.time(out <- sapply(lambda,
                                                function(x) do_ldet(x, env)))
         res[[type]][[3]] <- out
         rm(env)
         output[[dset]] <- res
         }
 save(output, file = "output_Jacobian.RData")
