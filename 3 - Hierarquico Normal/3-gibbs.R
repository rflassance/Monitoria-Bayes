################################################################################
################################################################################
##################____________________________________________##################
##################                                            ##################
##################        Drops 3: Amostrador de Gibbs        ##################
##################____________________________________________##################
##################                                            ##################
################################################################################
################################################################################

library(scales) #Melhora na visualização dos gráficos

# Contexto
## Xij|mu_i ~ N(mu_i, 1); i = 1, 2; j = 1, ..., ni
## mu_i|theta ~ N(theta, 1); i = 1, 2
## theta ~ N(0, 1)

# Condicionais completas
## mu_i|X, theta ~ N[(ni*Xi_barra + theta)/2, 1/2]
## theta|mu1, mu2 ~ N[(mu1 + mu2)/3, 1/3]

# Dados
n1 <- n2 <- 100
X1_barra <- 10.2
X2_barra <- 5.7

# Gibbs sampler
Gibbs <- function(X1, n1, X2, n2, theta0, B = 1000, seed = NULL){
  set.seed(seed)
  # Definindo os parametros
  mu1 <- mu2 <- theta <- numeric(B + 1)
  theta[1] <- theta0
  mu1[1] <- rnorm(1, (n1*X1 + theta[1])/(n1+1), sqrt(1/(n1+1)))
  mu2[1] <- rnorm(1, (n2*X2 + theta[1])/(n2+1), sqrt(1/(n2+1)))
  # Amostrador de Gibbs
  for(i in 2:(B+1)){ #B simulações
    theta[i] <- rnorm(1, (mu1[i-1] + mu2[i-1])/3, sqrt(1/3))
    mu1[i] <- rnorm(1, (n1*X1 + theta[i])/(n1+1), sqrt(1/(n1+1)))
    mu2[i] <- rnorm(1, (n2*X2 + theta[i])/(n2+2), sqrt(1/(n2+1)))
  }
  return(list(theta = theta, mu = cbind(mu1, mu2)))
}

# Amostra da posteriori
B <- 10000
parametros <- Gibbs(X1 = X1_barra, n1 = n1, X2 = X2_barra, n2 = n2,
                    theta0 = 0, B = B, seed = 42)
summary(parametros$theta)
summary(parametros$mu)

plot(0:B,0:B, type = "n", xlab = "Simulação", ylab = "Valor",
     ylim = c(min(parametros$theta, parametros$mu),
              max(parametros$theta, parametros$mu) + 4))
lines(0:B, parametros$theta, type = 'l')
lines(0:B, parametros$mu[,1], type = 'l', col = 'red')
lines(0:B, parametros$mu[,2], type = 'l', col = 'blue')
legend("topright", c(expression(theta), expression(mu[1]), expression(mu[2])),
       lty = 1, lwd = 2, col = c("black", "red", "blue"), bty = 'n')

minimo <- min(parametros$theta, parametros$mu)
maximo <- max(parametros$theta, parametros$mu)
opacidade <- 0.5
plot(density(parametros$theta, from = minimo, to = maximo), xlab = "Valores",
     ylab = "Densidade", main = NA, ylim = c(0,4), lwd = 2,
     col = alpha("black", .5))
lines(density(parametros$mu[,1], from = minimo, to = maximo), lwd = 2,
      col = alpha("red", opacidade))
lines(density(parametros$mu[,2], from = minimo, to = maximo), lwd = 2,
      col = alpha("blue", opacidade))
legend("topleft", c(expression(theta), expression(mu[1]), expression(mu[2])),
       lty = 1, lwd = 2, col = alpha(c("black", "red", "blue"), opacidade),
       bty = 'n')
