################################################################################
################################################################################
##################____________________________________________##################
##################                                            ##################
##################        Drops 4: Metropolis-Hastings        ##################
##################____________________________________________##################
##################                                            ##################
################################################################################
################################################################################

#Modelo: Y = X*beta + e, e ~ T_nu(0, sigma) => e|lambda ~ N(0, sigma2/lambda)
#beta|sigma2 ~ N(0, v2*sigma2)
#sigma2 ~ GamaInv(a, b)
#lambda|nu ~ Gama(nu/2, nu/2)
#nu ~ Fonseca et al (2008)

library(MASS) #ginv para a inversa

################################################################################
####################     Parte 1: nu fixado e conhecido     ####################
################################################################################
#B - simulações
#y - variável resposta
#X - covariáveis
#nu - graus de liberdade
#v2, a, b - parâmetros da priori
Gibbs <- function(B, y, X, nu, v2, a, b, seed = NULL){
  n <- length(y) #Tamanho da amostra
  n_cov <- dim(X)[2] #Número de covariáveis
  #Definindo os parâmetros de interesse (lambda não é de interesse!)
  beta <- matrix(NA, nrow = B+1, ncol = n_cov)
  sigma2 <- numeric(B+1)
  #Chutes iniciais (regressão linear frequentista)
  beta[1, ] <- ginv(t(X)%*%X)%*%t(X)%*%y
  sigma2[1] <- sum((y - X%*%beta[1,])^2)/(n - n_cov)
  #Amostrador de Gibbs
  set.seed(seed)
  for(i in 2:(B+1)){
    e2 <- (y - X%*%beta[i-1,])^2
    #lambda
    lambda <- rgamma(n, (nu+1)/2, 1/2*(nu + e2/sigma2[i-1]))
    #sigma2
    sigma2[i] <- 1/rgamma(1, (n+n_cov)/2 + a,
                         sum(lambda*e2)/2 + sum(beta[i-1,]^2)/(2*v2) + b)
    #beta
    V <- ginv(t(X)%*%X + diag(1/(v2*sigma2[i]), n_cov))
    m <- V%*%t(X)%*%y
    beta[i,] <- m + t(chol(V))%*%rnorm(n_cov)
  }
  return(list(beta = beta, sigma2 = sigma2, nu = nu))
}

set.seed(666)
n <- 200
beta <- c(2, .8, .2)
sigma2 <- .5
nu <- 5
X <- matrix(c(rep(1, n), runif(n), runif(n)), ncol = 3)
y <- rt(n, nu)*sqrt(sigma2) + X%*%beta

hist(y)

cadeias <- Gibbs(B = 1000, y = y, X = X, nu = nu,
                 v2 = 100, a = .01, b = .01, seed = 42)

par(mfrow = c(2, 2), mar = c(4,4,2,2))
for(i in 1:3){
  plot(cadeias$beta[,i], type = 'l', xlab = bquote(beta[.(i)]))
  abline(h = beta[i], col = 'red', lty = 2)
}
plot(cadeias$sigma2, type = 'l', xlab = expression(sigma^2))
abline(h = sigma2, col = 'red', lty = 2)

################################################################################
####################        Parte 2: nu desconhecido        ####################
################################################################################

priori.nu <- function(nu, log = T){
  den <- 1/2*log(nu/(nu+3)) +
    1/2*log(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/(nu*(nu+1)^2))
  if(!log) den <- exp(den)
  return(den)
}
MCMC <- function(B, y, X, nu0, v.nu, v2, a, b, seed = NULL){
  n <- length(y) #Tamanho da amostra
  n_cov <- dim(X)[2] #Número de covariáveis
  #Definindo os parâmetros de interesse (lambda não é de interesse!)
  beta <- matrix(NA, nrow = B+1, ncol = n_cov)
  sigma2 <- nu <- numeric(B+1)
  #Chutes iniciais (regressão linear frequentista)
  beta[1, ] <- ginv(t(X)%*%X)%*%t(X)%*%y
  sigma2[1] <- sum((y - X%*%beta[1,])^2)/(n - n_cov)
  nu[1] <- nu0
  #Taxa de aceitação
  aceita <- 0
  set.seed(seed)
  for(i in 2:(B+1)){
    #Resíduo ao quadrado
    e2 <- (y - X%*%beta[i-1,])^2
    #lambda (Gibbs)
    lambda <- rgamma(n, (nu[i-1]+1)/2, 1/2*(nu[i-1] + e2/sigma2[i-1]))
    
    #nu (MH)
    nu.prop <- exp(rnorm(1, log(nu[i-1]), sd = sqrt(v.nu)))
    log_a <-
      sum(dgamma(lambda, nu.prop/2, nu.prop/2, log = T) -
            dgamma(lambda, nu[i-1]/2, nu[i-1]/2, log = T)) + #log-vero
      priori.nu(nu.prop) - priori.nu(nu[i-1]) + #log-priori
      dnorm(log(nu[i-1]), mean = log(nu.prop), sd = sqrt(v.nu), log = T) -
      dnorm(log(nu.prop), mean = log(nu[i-1]), sd = sqrt(v.nu), log = T) #prop
    ifelse(log_a >= log(runif(1)),
           c(nu[i] <- nu.prop, aceita <- aceita + 1),
           nu[i] <- nu[i-1])
    #sigma2 (Gibbs)
    sigma2[i] <- 1/rgamma(1, (n+n_cov)/2 + a,
                          sum(lambda*e2)/2 + sum(beta[i-1,]^2)/(2*v2) + b)
    #beta (Gibbs)
    V <- ginv(t(X)%*%X + diag(1/(v2*sigma2[i]), n_cov))
    m <- V%*%t(X)%*%y
    beta[i,] <- m + t(chol(V))%*%rnorm(n_cov)
  }
  return(list(beta = beta, sigma2 = sigma2, nu = nu, aceita = aceita/B))
}

cadeias <- MCMC(B = 10000, y, X, nu0 = 3, v.nu = .055,
                v2 = 100, a = .01, b = .01, seed = 42)
cadeias$aceita

par(mfrow = c(2, 2), mar = c(4,4,2,2))
for(i in 1:3){
  plot(cadeias$beta[,i], type = 'l', xlab = bquote(beta[.(i)]))
  abline(h = beta[i], col = 'red', lty = 2)
}
plot(cadeias$sigma2, type = 'l', xlab = expression(sigma^2))
abline(h = sigma2, col = 'red', lty = 2)

par(mfrow = c(1, 1), mar = c(4,4,2,2))
plot(cadeias$nu, type = 'l', xlab = expression(nu))
abline(h = nu, col = 'red', lty = 2)


################################################################################
##################        Parte 3: outros diagnósticos        ##################
################################################################################

acf(cadeias$nu[-(1:7000)], lag = 150)

MCMC2 <- function(B, y, X, nu0, v.nu, v2, a, b, seed = NULL){
  n <- length(y) #Tamanho da amostra
  n_cov <- dim(X)[2] #Número de covariáveis
  #Definindo os parâmetros de interesse (lambda não é de interesse!)
  beta <- matrix(NA, nrow = B+1, ncol = n_cov)
  sigma2 <- nu <- numeric(B+1)
  #Chutes iniciais (regressão linear frequentista)
  beta[1, ] <- ginv(t(X)%*%X)%*%t(X)%*%y
  sigma2[1] <- sum((y - X%*%beta[1,])^2)/(n - n_cov)
  nu[1] <- nu0
  #Taxa de aceitação
  aceita <- 0
  set.seed(seed)
  for(i in 2:(B+1)){
    #Resíduo ao quadrado
    e <- (y - X%*%beta[i-1,])
    #lambda (Gibbs)
    lambda <- rgamma(n, (nu[i-1]+1)/2, 1/2*(nu[i-1] + e^2/sigma2[i-1]))
    
    #nu (MH)
    nu.prop <- exp(rnorm(1, log(nu[i-1]), sd = sqrt(v.nu)))
    
    log_a <-
      sum(dt(e/sqrt(sigma2[i-1]), nu.prop, log = T) -
            dt(e/sqrt(sigma2[i-1]), nu[i-1], log = T)) + #log-vero
      priori.nu(nu.prop) - priori.nu(nu[i-1]) + #log-priori
      dnorm(log(nu[i-1]), mean = log(nu.prop), sd = sqrt(v.nu), log = T) -
      dnorm(log(nu.prop), mean = log(nu[i-1]), sd = sqrt(v.nu), log = T) #prop
    
    ifelse(log_a >= log(runif(1)),
           c(nu[i] <- nu.prop, aceita <- aceita + 1),
           nu[i] <- nu[i-1])
    #sigma2 (Gibbs)
    sigma2[i] <- 1/rgamma(1, (n+n_cov)/2 + a,
                          sum(lambda*e^2)/2 + sum(beta[i-1,]^2)/(2*v2) + b)
    #beta (Gibbs)
    V <- ginv(t(X)%*%X + diag(1/(v2*sigma2[i]), n_cov))
    m <- V%*%t(X)%*%y
    beta[i,] <- m + t(chol(V))%*%rnorm(n_cov)
  }
  return(list(beta = beta, sigma2 = sigma2, nu = nu, aceita = aceita/B))
}

cadeias <- MCMC2(B = 10000, y, X, nu0 = 3, v.nu = .5,
                 v2 = 100, a = 0.01, b = 0.01, seed = 42)
cadeias$aceita

plot(cadeias$nu, type = 'l', xlab = expression(nu))
abline(h = nu, col = 'red', lty = 2)

acf(cadeias$nu, lag = 150)
acf(cadeias$nu, lag = 30)

acf(cadeias$beta[,1], lag = 30)
acf(cadeias$beta[,2], lag = 30)
acf(cadeias$beta[,3], lag = 30)
acf(cadeias$sigma2, lag = 30)

library(coda)
matMCMC <- cbind(cadeias$beta, cadeias$sigma2, cadeias$nu)
colnames(matMCMC) <- c("beta0", "beta1", "beta2", "sigma2", "nu")
matMCMC <- mcmc(data = matMCMC, start = 3000, thin = 30)
acfplot(matMCMC)
crosscorr.plot(matMCMC)
densplot(matMCMC)
effectiveSize(matMCMC)
geweke.plot(matMCMC)
HPDinterval(matMCMC, prob = 0.95)
rejectionRate(matMCMC)
