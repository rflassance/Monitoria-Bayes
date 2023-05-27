//
// Modelo:
//    Y = X*beta + e, e ~ T_nu(0, sqrt(sigma2))
//    beta|sigma2 ~ N(0, v2*sigma2)
//    sigma2 ~ GamaInv(a, b)
//    nu ~ Fonseca et al (2008)
//
// Mais sobre o Stan:
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// Data: inserção de informações externas ao modelo
data {
  // Tamanho da amostra
  int<lower=0> N;
  // Número de covariáveis
  int<lower=0> M;
  // Variável resposta
  vector[N] y;
  // Matriz das covariáveis
  matrix[N,M] X;
  // Informações das prioris
  real<lower=0> v2;
  real<lower=0> a;
  real<lower=0> b;
}

// Parameters: Detalhamento dos parâmetros existentes no modelo
parameters {
  // Vetor de betas
  vector[M] beta;
  // Escala da T
  real<lower=0> sigma2;
  // Grau de liberdade
  real<lower=0> nu;
}

// Model: Suposições acerca da distribuição dos dados e dos parâmetros
model {
  // Verossimilhança
  for(i in 1:N){
    y[i] ~ student_t(nu, dot_product(beta, X[i,:]), sqrt(sigma2));
  }
  // Prioris
  beta ~ normal(0, sqrt(v2*sigma2));
  sigma2 ~ inv_gamma(a, b);
  //// Priori para nu é de Jeffreys, exige uma construção própria
  target += 0.5*log(nu/(nu+3)) + 0.5*log(trigamma(nu/2)-trigamma((nu+1)/2)-2*(nu+3)/(nu*(nu+1)^2));
}

