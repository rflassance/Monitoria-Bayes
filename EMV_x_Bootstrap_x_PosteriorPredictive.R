library(VGAM)

dados <- c(rep(0, 92), rep(1, 8))
alfa <- 9
beta <- 90

# Posteriori
alfa_n <- alfa + sum(dados)
beta_n <- beta + length(dados) - sum(dados)
dbbinom <- dbetabinom.ab(x = 0:50, size = 50, shape1 = alfa_n, shape2 = beta_n)

# Bootstrap
B <- 100000
Ttil <- numeric(B)
for(i in 1:B){
  theta <- sum(sample(dados, length(dados), replace = T))/length(dados)
  Ttil[i] <- rbinom(1, size = 50, theta)
}

hist(Ttil, col = NULL, border = 'blue', freq = F, ylab = "Probablilidade",
     xlab = expression(tilde(T)), lwd = 2, main = NULL, ylim = c(0,0.2))
lines(0:50, dbinom(x = 0:50, size = 50, sum(dados)/length(dados)), col = "green",
     lwd = 2, main = NULL, type = 's')
lines(0:50, dbbinom, col = 'red', type = 's', lwd = 2)

legend("topright", c("Posterior predictive", "Bootstrap", "EMV"), lwd = 1,
       col = c("red", "blue", "green"))
