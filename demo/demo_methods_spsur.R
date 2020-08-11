### demo with the whole set of examples of methods_spsur ########
 rm(list = ls()) # Clean memory
 data(spc)
 Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
 spcsur.sim <-spsurml(formula = Tformula, data = spc, type = "sim")
 spcsur.slm <-spsurml(formula = Tformula, data = spc, type = "slm", 
                      listw = Wspc)
 # ANOVA table and LR test for nested models:
 anova(spcsur.sim, spcsur.slm)
 ## Print Table       
 print(spcsur.slm)
 ## Plot spatial and beta coefficients
 # Interactive plot
 plot(spcsur.slm)
 # Non-interactive plot
 if (require(gridExtra)) {
   pl <- plot(spcsur.slm, viewplot = FALSE) 
   grid.arrange(pl$lplbetas[[1]], pl$lplbetas[[2]], 
                pl$pldeltas, nrow = 3) }
