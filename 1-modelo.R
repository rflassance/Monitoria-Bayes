################################################################################
################################################################################
##################____________________________________________##################
##################                                            ##################
##################        Drops 1: Modelagem bayesiana        ##################
##################____________________________________________##################
##################                                            ##################
################################################################################
################################################################################

# Parametros da priori
## Crença - Ignorância completa
alfa <- 1
beta <- 1
curve(dbeta(x, shape1 = alfa, shape2 = beta), ylab = "Densidade",
      ylim = c(0,8), main = "Prioris", lwd = 2)
## Consenso - Priori de Jeffreys
alfa <- 1/2
beta <- 1/2
curve(dbeta(x, shape1 = alfa, shape2 = beta), add = T, col = "blue", lwd = 2)
## Crítica - O que queremos testar?
alfa <- 20
beta <- 20
curve(dbeta(x, shape1 = alfa, shape2 = beta), add = T, col = "red", lwd = 2)
abline(v = 0.5, lty = 2)
legend("topright", c("Crença", "Consenso", "Crítica"), lty = 1, lwd = 2,
       col = c("black", "blue", "red"))

# Posterioris
# Tamanho da amostra
n <- 20
# Dados observados
X <- 13
plot_post <- function(n, X, alfa, beta, titulo = NULL, ylim = c(0,1)){
  curve(dbeta(x, shape1 = alfa, shape2 = beta), ylab = "Densidade",
        main = titulo, lwd = 2, ylim = ylim)
  curve(dbeta(x, shape1 = X + alfa, shape2 = n - X + beta), add = T,
        col = "blue", lwd = 2)
  abline(v = 0.5, lty = 2)
  legend("topright", c("Priori", "Posteriori"), lty = 1, lwd = 2,
         col = c("black", "blue"))
}
par(mfrow = c(1,3), mar = c(3,3,3,3))
## Crença
plot_post(n, X, alfa = 1, beta = 1, titulo = "Crença", ylim = c(0,8))
## Consenso
plot_post(n, X, alfa = 1/2, beta = 1/2, titulo = "Consenso", ylim = c(0,8))
## Crítica
plot_post(n, X, alfa = 20, beta = 20, titulo = "Crítica", ylim = c(0,8))
