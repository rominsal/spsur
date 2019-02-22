X_restr <- function(X,R,b,p){
  # Function to fit beta restricted imposing the restrictions in matrix X
  Rt <- t(R)
  # Localizar coeficientes no nulos de Rt y b
  Rt_nonull <- Rt != 0
  rt_nonull <- b != 0
  # Primer coeficiente no nulo de cada columna en Rt es variable referencia
  # (mantiene signo). El resto de los coeficientes no nulos de la columna
  #cambiant signo
  Xstar <- X
  index_var_ref <- integer(ncol(Rt))
  index_var_del <- integer(ncol(Rt)) # All columns are deleted at the end...
  coef_var_ref <- rep(0,ncol(Rt))
  for (i in 1:ncol(Rt)) {
   Rt_col <- Rt[,i]
   index_Rt_col <- which(Rt_col != 0)
   nvar_Rt_col <- length(index_Rt_col)
   index_var_ref[i] <- index_Rt_col[1]
   coef_var_ref[i] <- Rt_col[index_var_ref[i]]
   for (j in 2:nvar_Rt_col) {
     index_var_del[i] <- index_Rt_col[j]
     if (!any(index_var_del[1:i] == index_Rt_col[j])) break
   }
   Rt_col_neg <- (-1)*Rt_col
   Rt_col_neg[index_var_ref[i]] <- coef_var_ref[i]
   newcolXstar <- Xstar %*% Rt_col_neg
   Xstar[,index_var_ref[i]] <- newcolXstar
  }
  # Si b es no nulo cambia el intercepto
  for (i in 1:nrow(b)) {
    if (rt_nonull[i]) X_star[,1] <- X_star[,1]+b[i]*X[,index_var_ref[i]]
  }
  # Elimina variables dependientes linealmente en Xstar y adapta p de cada ecuación
  Xstar <- Xstar[,-c(index_var_del)]
  colnames(Xstar) <- colnames(X)[-c(index_var_del)]
  pstar <- p
  for (i in 1:length(index_var_del))
  {
    index_pstar_i <- which(index_var_del[i] < cumsum(p))[1]
    pstar[index_pstar_i] <- pstar[index_pstar_i] - 1
  }
  res <- list(Xstar = Xstar, pstar = pstar)
  return (res)
}
  # CÓDIGO ANTIGUO
  #
  #   for (j in 1:nrow(Rt)) {
  #     if (Rt_nonull[j,i]) {
  #       var_ref[i] <- j
  #       coef_var_ref[i] <- Rt[j,i]
  #       # Chequea que no coincide con columna anterior CAMBIAR
  #       # if (i>1) {
  #       #   if (var_ref[i] != var_ref[i-1]) break
  #       # } else break
  #     }
  #   }
  # }
  #
  #
  #
  # # Último coeficiente de cada columna en Rt es variable a eliminar
  # # (hay un parámetro menos por cada restricción impuesta)
  # var_del <- integer(ncol(Rt))
  # for (i in 1:ncol(Rt)) {
  #   for (j in nrow(Rt):1) {
  #     if (Rt_nonull[j,i]) {
  #       var_del[i] <- j
  #       # Chequea que no coincide con columna anterior
  #       # if (i>1) {
  #       #   if (var_del[i]!=var_del[i-1]) break
  #       # } else break
  #     }
  #   }
  # }
  #
  # # Chequeo para comprobar la compatibilidad de las variables elegidas
  # for (i in 1:ncol(Rt)) {
  #   if (var_ref[i]==var_del[i]) stop("Restrictions can not be substituted in matrix X")
  # }
  #
  #
  #
  # # Cambia el signo de todos los coeficientes excepto las variables de referencia
  # Rt_star <- Rt*(-1)
  # for (i in 1:length(var_ref)) Rt_star[var_ref[i],i] <- coef_var_ref[i]
  #
  # # Creación variables restringidas
  # X_star <- X
  # for (i in 1:ncol(Rt_star)) {
  #   xnew_star <- X %*% Rt_star[,i]
  #   X_star[,var_ref[i]] <- xnew_star
  # }
  # X_star <- X_star[,-c(var_del)]
  #
  # # Si b es no nulo cambia el intercepto
  # for (i in 1:nrow(b)) {
  #   if (rt_nonull[i]) X_star[,1] <- X_star[,1]+b[i]*X[,var_ref[i]]
  # }
