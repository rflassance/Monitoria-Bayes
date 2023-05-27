################################################################################
################################################################################
##################____________________________________________##################
##################                                            ##################
##################       Drops 2: Simulação condicional       ##################
##################____________________________________________##################
##################                                            ##################
################################################################################
################################################################################

# Para uso da T de Student com média e escala
library(ggdist)

# Parametros a priori
m <- 0
v <- 10
curve(dnorm(x,m,v), from = -40, to = 40, main = expression(Priori ~ para ~ mu),
      ylab = "Densidade", xlab = expression(mu))
a <- 2
b <- .5
curve(dgamma(x,a,b), from = 0, to = 15, main = expression(Priori ~ para ~ phi),
      ylab = "Densidade", xlab = expression(phi))

# Valores reais
mu <- 8
phi <- 0.5 # Variância  2

# Dados
n <- 50
set.seed(42)
Y <- rnorm(n, mean = mu, sd = sqrt(1/phi))

# Máxima verossimilhança
Ybarra <- mean(Y)
S2 <- sum((Y - Ybarra)^2)/n

# Parametros a posteriori
a_novo <- n/2 + a
b_novo <- 1/2*(n*S2 + (Ybarra - m)^2/(n+v) + b)
v_novo <- 1/(n + 1/v)
m_novo <- (n*Ybarra + m/v)*v_novo

# A marginal de mu é mesmo T de Student?
nu <- a_novo * 2

B <- 10000
phi_sim <- rgamma(B, nu/2, nu/2)
mu_sim <- rnorm(B, m_novo, sd = sqrt(v_novo * b_novo/nu * phi_sim))

curve(dstudent_t(x, df = nu, mu = m_novo, sigma = sqrt(v_novo * b_novo/nu)),
      from = 7, to = 9, ylab = "Densidade",
      xlab = expression(mu ~ symbol("|") ~ y))
hist(mu_sim, add = T, freq = F)

# Comparação dos resultados

## mu|Y
curve(dstudent_t(x, df = nu, mu = m_novo, sigma = sqrt(v_novo * b_novo/nu)),
      from = 7, to = 9, ylab = "Densidade",
      xlab = expression(mu ~ symbol("|") ~ y))
abline(v = Ybarra, col = 'red', lwd = 2)
abline(v = mu, col = 'blue', lwd = 2)
legend("topright", c("Verdadeiro", "EMV"), lty = 1, lwd = 2,
       col = c("blue", "red"))

## phi|Y
curve(dgamma(x, a_novo, b_novo), from = 0, to = 1, ylab = "Densidade",
      xlab = expression(phi ~ symbol("|") ~ y))
abline(v = 1/S2, col = 'red', lwd = 2)
abline(v = phi, col = 'blue', lwd = 2)
legend("topright", c("Verdadeiro", "EMV"), lty = 1, lwd = 2,
       col = c("blue", "red"))
