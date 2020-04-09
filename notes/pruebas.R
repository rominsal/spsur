library(spatialreg)
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1], quiet=TRUE)
#require("spdep", quietly=TRUE)
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])
listw <- spdep::nb2listw(col.gal.nb)
ev <- eigenw(listw)
lobjHess <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw,
                 control=list(pre_eig=ev, fdHess = TRUE)) 
lobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw,
                     control=list(pre_eig = ev)) 
## usa fdHESS para posterior impacts... 
summary(lobj)
mobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE,
                 control=list(pre_eig=ev))
summary(mobj)
mobj1 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, Durbin= ~ INC,
                  control=list(pre_eig=ev))
summary(mobj1)
sacobj <- sacsarlm(CRIME ~ INC + HOVAL, columbus, listw)
summary(sacobj)

sacmobj <- sacsarlm(CRIME ~ INC + HOVAL, columbus, listw, Durbin=TRUE)
summary(sacmobj)



lobj$dvars; mobj$dvars; mobj1$dvars
lobj$icept
