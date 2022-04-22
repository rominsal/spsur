# @name betas_restr
# @rdname betas_restr
#
# @title betas restricted
#
# @description
# Restricted estimation of beta coefficients
#
# @param    results : An object create with \code{\link{spsurml}}.
# @param    R       : Coefficient matrix for betas.
# @param    r       : Vector of independent terms.
#
# @details
# ¿?¿?
#
# @return
# betas restricted and covariance matrix
#
# @references
#
# CAMBIAR
# J. LeSage and R.K. Pace. \emph{Introduction to Spatial Econometrics}, CRC Press, chapter 10.1.6, 2009
# Mur, J., Lopez, F., & Herrera, M. (2010). Testing for spatial effects in seemingly unrelated regressions. \emph{Spatial Economic Analysis}, 5(4), 399-440.
# \cr
# \cr
# Lopez, F.A., Mur, J., & Angulo, A. (2014). Spatial model selection strategies in a SUR framework. The case of regional productivity in EU. \emph{Annals of Regional Science}, 53(1), 197-220.
# \cr
# \cr
# Lopez, F.A., Martinez-Ortiz, P.J., & Cegarra-Navarro, J.G. (2017). Spatial spillovers in public expenditure on a municipal level in Spain. \emph{Annals of Regional Science}, 58(1), 39-65.
#
# @author
#   \tabular{ll}{
#   Fernando Lopez  \tab \email{fernando.lopez@@upct.es} \cr
#   Roman Minguez  \tab \email{roman.minguez@@uclm.es} \cr
#   Jesus Mur  \tab \email{jmur@@unizar.es} \cr
#   }
# @seealso
# \code{\link{spsur}}
# @examples
# data(spc)
# Tformula <- WAGE83 | WAGE81 ~ UN83 + NMR83 + SMSA | UN80 + NMR80 + SMSA
# ## Estimate SUR-SLM model
# spcsur.slm <-spsurml(Form=Tformula,data=spc,type="slm",W=Wspc)
# summary(spcsur.slm)
# ## H_0: equality between SMSA coefficients in both equations.
# R1 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
# r1 <- matrix(0,ncol=1)
# wald_betas(results=spcsur.slm,R=R1,r=r1)
# betas_rest1 <- betas_restr(spcsur.slm,R=R1,r=r1)
# ## Estimate SUR-SEM model
# spcsur.sem <-spsurml(Form=Tformula,data=spc,type="sem",W=Wspc)
# summary(spcsur.sem)
# ## H_0: equality between intercepts and SMSA coefficients in both equations.
# R2 <- matrix(c(1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,-1),nrow=2,ncol=8,byrow=TRUE)
# r2 <- matrix(c(0,0),ncol=1)
# res2 <- wald_betas(results=spcsur.sem,R=R2,r=r2)
# betas_rest2 <- betas_restr(spcsur.sem,R=R2,r=r2)
# ## Estimate SUR-SARAR model
# spcsur.sarar <-spsurml(Form=Tformula,data=spc,type="sarar",W=Wspc)
# summary(spcsur.sarar)
# ## H_0: equality between SMSA coefficients in both equations.
# R3 <- matrix(c(0,0,0,1,0,0,0,-1),nrow=1)
# r3 <- matrix(0,ncol=1)
# wald_betas(results=spcsur.sarar,R=R3,r=r3)
# betas_rest3 <- betas_restr(spcsur.sarar,R=R3,r=r3)


betas_restr <- function(results , R , r){
  z <- results # OBJETO QUE INCLUYE ESTIMACIÓN EN Rbetas <- z$betas
  betas <- Matrix::Matrix(matrix(z$betas,ncol=1))
  rownames(betas) <- names(z$betas)
  cov_betas <- Matrix::Matrix(z$cov[rownames(betas),rownames(betas)])
  R <- Matrix::Matrix(R)
  colnames(R) <- rownames(betas)
  r <- Matrix::Matrix(matrix(r,ncol=1))
  holg <- R %*% betas - r
  q <- nrow(R)

  X <- Matrix::Matrix(z$X)
  W <- Matrix::Matrix(z$W)
  Sigma <- Matrix::Matrix(z$Sigma)
  Sigma_inv <- Matrix::solve(Sigma)
  Tm <- z$Tm
  N <- z$N
  G <- z$G
  IT <- Matrix::Diagonal(Tm)
  IR <- Matrix::Diagonal(N)
  IGR <- Matrix::Diagonal(G*N)
  OME <- kronecker(IT,kronecker(Sigma,IR))
  OMEinv <- kronecker(IT,kronecker(Sigma_inv,IR))
  type <- z$type
  deltas <- Matrix::Matrix(matrix(z$deltas,ncol=1))
  rownames(deltas) <- names(z$deltas)
  cov_deltas <- Matrix::Matrix(z$cov[rownames(deltas),rownames(deltas)])
  if (type=="sem" || type=="sdem" || type=="sarar")
  {
    rhos <- deltas[grepl("rho",rownames(deltas))]
    rho_matrix <- Matrix::Matrix(diag(as.vector(rhos)))
    B <- kronecker(IT,(IGR - kronecker(rho_matrix,W)))
    X <- B %*% X
  }
  XtOMEinvX <- Matrix::t(X) %*% (OMEinv %*% X)
  rownames(XtOMEinvX) <- rownames(betas)
  colnames(XtOMEinvX) <- rownames(betas)
  XtOMEinvX_inv <- Matrix::solve(XtOMEinvX)
  rownames(XtOMEinvX_inv) <- rownames(betas)
  colnames(XtOMEinvX_inv) <- rownames(betas)

  betas_r <- betas + XtOMEinvX_inv %*% (Matrix::t(R) %*%
                           Matrix::solve(R %*% XtOMEinvX_inv %*%
                                           Matrix::t(R),-holg))

  cov_betas_r <- XtOMEinvX_inv - XtOMEinvX_inv %*%
                   (Matrix::t(R) %*%
                      Matrix::solve(R %*% XtOMEinvX_inv %*%
                                     Matrix::t(R),R %*% XtOMEinvX_inv))
  rownames(cov_betas_r) <- rownames(betas_r)
  colnames(cov_betas_r) <- rownames(betas_r)
  se_betas_r <- sqrt(Matrix::diag(cov_betas_r))
  t_betas_r <- betas_r / se_betas_r
  ##### Print Table of beta restricted
  p <- z$p
  rdf <- z$df.residual+q
  coef_table <- list(NULL)
  # Build coefficients table by Equation
  for (i in 1:G)
  {
    if(i==1){
      coef_table[[i]] <- cbind(betas_r[1:p[i]], se_betas_r[1:p[i]],
                                 t_betas_r[1:p[i]],
                                 2 * pt(abs(t_betas_r[1:p[i]]),rdf,
                                        lower.tail = FALSE))
      colnames(coef_table[[i]]) <- c("Estimate", "Std. Error",
                                        "t value", "Pr(>|t|)")
      rownames(coef_table[[i]]) <- rownames(betas_r)[1:p[i]]
    } else {
      coef_table[[i]] <-  cbind(betas_r[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                  se_betas_r[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                  t_betas_r[(cumsum(p)[i-1]+1):cumsum(p)[i]],
                                  2 * pt(abs(t_betas_r[(cumsum(p)[i-1]+1):cumsum(p)[i]]),rdf,
                                         lower.tail = FALSE))

      rownames(coef_table[[i]]) <- rownames(betas_r)[(cumsum(p)[i-1]+1):cumsum(p)[i]]
    }

    colnames(coef_table[[i]]) <- c("Estimate", "Std. Error",
                                     "t value", "Pr(>|t|)")
  }

  digits <- max(3L, getOption("digits") - 3L)

  cat("\n Betas Restricted: \n\n")
  for (i in 1:length(coef_table)){
    cat("Equation ",i,"\n")
    printCoefmat(coef_table[[i]], P.values = TRUE, has.Pvalue = TRUE)
  }
  res <- list(betas_restr = as.matrix(betas_r),
             cov_betas_restr = as.matrix(cov_betas_r),
             betas_unrestr = as.matrix(betas),
             cov_betas_unrestr = as.matrix(cov_betas))

}
