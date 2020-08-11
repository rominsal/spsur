### demo with the whole set of examples of dgp_spsur() ########

####################################
######## CROSS SECTION DATA ########
####################################

####################################
#### Example 1: DGP SLM model
####################################
rm(list = ls()) # Clean memory
Tm <- 1 # Number of time periods
G <- 3 # Number of equations
N <- 200 # Number of spatial elements
p <- 3 # Number of independent variables
Sigma <- matrix(0.3, ncol = G, nrow = G)
diag(Sigma) <- 1
Betas <- c(1, 2, 3, 1, -1, 0.5, 1, -0.5, 2)
rho <- 0.5 # level of spatial dependence
lambda <- 0.0 # spatial autocorrelation error term = 0
#  random coordinates
co <- cbind(runif(N,0,1),runif(N,0,1))
lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(co, k = 5,
                                                   longlat = FALSE)))
DGP <- dgp_spsur(Sigma = Sigma, Betas = Betas,
                 rho = rho, lambda = lambda, Tm = Tm,
                 G = G, N = N, p = p, listw = lw)
SLM <- spsurml(X = DGP$X, Y = DGP$Y, Tm = Tm, N = N, G = G, 
               p = c(3, 3, 3), listw = lw, type = "slm",
               method = "LU", control = list(fdHess = TRUE)) 
summary(SLM)
####################################
### Example 2: DGP SEM model with Tm>1; G=1 and
### different p for each equation
####################################
rm(list = ls()) # Clean memory
Tm <- 3 # Number of time periods
G <- 1 # Number of equations
N <- 500 # Number of spatial elements
p <- c(2,3,4) # Number of independent variables
Sigma <- matrix(0.8, ncol = Tm, nrow = Tm)
diag(Sigma) <- 1
Betas <- c(1,2,1,2,3,1,2,3,4)
rho <- 0 # level of spatial dependence = 0
lambda <- c(0.2,0.5,0.8) 
# spatial autocorrelation error terms for each equation
# random coordinates
co <- cbind(runif(N,0,1),runif(N,0,1))
lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(co, k = 5,
                                                   longlat = FALSE)))
DGP2 <- dgp_spsur(Sigma = Sigma, Betas = Betas, rho = rho, 
                  lambda = lambda, Tm = Tm, G = G, N = N, p = p, 
                  listw = lw)
SLM2 <- spsurml(X = DGP2$X, Y = DGP2$Y, Tm = Tm, N = N, G = G,
               p = c(2,3,4), listw = lw, type = "slm", 
               control = list(fdHess = TRUE))
summary(SLM2)
SEM2 <- spsurml(X = DGP2$X, Y = DGP2$Y, Tm = Tm, N = N, G = G,
               p = c(2,3,4), listw = lw, type = "sem",
               method = "LU", control = list(fdHess = TRUE))
summary(SEM2)

####################################
#### Example with G>1 and Tm>>1
####################################
## It usually requires 1-2 minutes maximum
rm(list = ls()) # Clean memory
Tm <- 10 # Number of time periods
G <- 3 # Number of equations
N <- 100 # Number of spatial elements
p <- 3 # Number of independent variables
Sigma <- matrix(0.5, ncol = G, nrow = G)
diag(Sigma) <- 1
Betas <- rep(1:3, G)
rho <- c(0.5, 0.1, 0.8)
lambda <- 0.0 # spatial autocorrelation error term = 0
# random coordinates
co <- cbind(runif(N,0,1),runif(N,0,1))
lw <- spdep::nb2listw(spdep::knn2nb(spdep::knearneigh(co, k = 5,
                                                   longlat = FALSE)))
DGP3 <- dgp_spsur(Sigma = Sigma, Betas = Betas, rho = rho, 
                  lambda = lambda, Tm = Tm, G = G, N = N, p = p, 
                  listw = lw)
SLM3  <-spsurml(Y = DGP3$Y, X = DGP3$X, G = G, N = N, Tm = Tm,
                p = p, listw = lw, type = "slm", 
                method = "LU", control = list(fdHess = TRUE))
summary(SLM3)

