################################################################################
################################################################################
##################____________________________________________##################
##################                                            ##################
##################        Drops 5: Introdução ao Stan         ##################
##################____________________________________________##################
##################                                            ##################
################################################################################
################################################################################

#Modelo: Y = X*beta + e, e ~ T_nu(0, sigma) => e|lambda ~ N(0, sigma2/lambda)
#beta|sigma2 ~ N(0, v2*sigma2)
#sigma2 ~ GamaInv(a, b)
#lambda|nu ~ Gama(nu/2, nu/2)
#nu ~ Fonseca et al (2008)

library(rstan) #Exige Rtools
library(parallel)

#Geração dos dados
set.seed(666)
n <- 200
beta <- c(2, .8, .2)
sigma2 <- .5
nu <- 5
X <- matrix(c(rep(1, n), runif(n), runif(n)), ncol = 3)
y <- rt(n, nu)*sqrt(sigma2) + X%*%beta

hist(y)

#Aplicação do modelo Stan
aux = stan("5 - Regressao Linear Stan/regressao.stan",
           data = list(N = n, M = dim(X)[2], y = y[,,drop=T], X = X,
                       v2 = 100, a = .01, b = .01),
           cores = 4, iter = 100, chains = 4, seed = 42)
aux = stan("5 - Regressao Linear Stan/regressao.stan",
           data = list(N = n, M = dim(X)[2], y = y[,,drop=T], X = X,
                       v2 = 100, a = .01, b = .01),
           cores = 4, iter = 1000, warmup = 100, chains = 4, seed = 42)

summary(aux)
plot(aux)

aux = stan("5 - Regressao Linear Stan/regressao.stan",
           data = list(N = n, M = dim(X)[2], y = y[,,drop=T], X = X,
                       v2 = 100, a = .01, b = .01),
           cores = 4, iter = 6000, warmup = 3000, chains = 4, seed = 42)
plot(aux)

posteriori <- extract(aux)

names(posteriori)

par(mfrow = c(2,3), mar = c(4,2,2,2))
for(i in 1:length(beta)){
  plot(posteriori$beta[,i], type = 'l', xlab = bquote(beta[.(i)]), ylab = '')
  abline(h = beta[i], col = 'red', lty = 2)
}

plot(posteriori$sigma2, type = 'l', xlab = expression(sigma^2), ylab = '')
abline(h = sigma2, col = 'red', lty = 2)

plot(posteriori$nu, type = 'l', xlab = expression(nu), ylab = '')
abline(h = nu, col = 'red', lty = 2)
