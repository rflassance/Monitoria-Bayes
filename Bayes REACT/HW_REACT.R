################################################################################
#################______________________________________________#################
#################                                              #################
#################  Exemplo do REACT bayesiano para Dirichlet   #################
#################______________________________________________#################
#################                                              #################
################################################################################

# Carregando bibliotecas
library(dplyr)
library(ggplot2)
library(ggpubr)
library(DirichletReg)

################################################################################
######################## Carregando dados de alelos A,B ########################
################################################################################

# Contexto: Avaliação de se alelos estão passando por pressões evolutivas
# Combinações de alelos possíveis: {A,A}, {A,B}, {B,B}
# Sob a hipótese nula, se P(A) = p = 1 - P(B), então
##      X|theta ~  Multinomial(n, theta)
##      theta = [P(A,A), P(A,B), P(B,B)] = [p^2, 2*p*(1-p), (1-p)^2]
data = "Bayes REACT/brentani.csv" %>%
  read.csv() %>%
  as_tibble()
View(data)


################################################################################
######################## Avaliação visual da pragmática ########################
################################################################################
# theta3 = 1 - theta1 - theta2, logo só preciso avaliar dois parâmetros

# Grid para construção da pragmática
B = 500 # Termos de cada theta para a grade
thetas = seq(0.0000001, 0.99999999, length.out = B)
color_grid = tibble(theta1 = rep(thetas, each = B), #Valores para theta1 (livre)
                    theta2 = rep(thetas, B), #Valores para theta2 (livre)
                    color = rep(FALSE, B^2)) %>% #Cor (gráfico da pragmática)
  mutate(theta3 = 1 - theta1 - theta2) %>% #Aplicando a restrição
  filter(theta3 > 0, theta3 < 1) #Removendo casos que não fazem sentido

# Altera a grade para gerar cores quando a combinação estiver na pragmática
generate_grid = function(color_grid, dissimilarity, p0, eps){
  # Calcula a dissimilatidade para cada combinação
  diss <- with(color_grid,
               dissimilarity(theta1, theta2, theta3, p0)
               )
  color_grid %>%
    mutate(diss = diss, color = (diss < eps))
}


# Distância L1: |theta1 - p^2| + |theta2 - 2p(1-p)| + |theta3 - (1-p)^2|
ABS_diss = function(t1, t2, t3, p0) {
  abs(t1 - p0^2) + abs(t2 - 2*p0*(1-p0)) + abs(t3 - (1-p0)^2)
}

# Euclidiana: sqrt((theta1 - p^2)^2 + (theta2 - 2p(1-p))^2 + (theta3 - (1-p)^2)^2)
L2_diss = function(t1, t2, t3, p0) {
  sqrt((t1 - p0^2)^2 + (t2 - 2*p0*(1-p0))^2 + (t3 - (1-p0)^2)^2)
}

# Distância Linf: max(|theta1 - p^2|, |theta2 - 2p(1-p)|, |theta3 - (1-p)^2|)
MAX_diss = function(t1, t2, t3, p0) {
  pmax(abs(t1 - p0^2), abs(t2 - 2*p0*(1-p0)), abs(t3 - (1-p0)^2))
}

plot_grid = function(grid, col = "lightblue", alpha = 0.9){
  #Grade apenas para os casos na pragmática e agrupados
  grid %<>% filter(color) %>% group_by(theta1)
  #Menor e maior theta3 para cada theta1
  poly_min = grid %>% summarise(theta3 = min(theta3))
  poly_max = grid %>% 
    summarise(theta3 = max(theta3)) %>%
    arrange(desc(theta1))
  poly_grid = bind_rows(poly_min, poly_max)
  poly_grid %>%
    ggplot(aes(x = theta1, y = theta3)) +
    theme_minimal() +
    geom_polygon(color = NA, fill = col, alpha = alpha) +
    #Limites do espaço paramétrico
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0)) + 
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1)) +
    geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0)) +
    #H0
    stat_function(fun = function(x) (1-sqrt(x))^2, n = 101,
                  linetype="dashed", color = "black") +
    xlim(c(0,1)) + ylim(c(0,1)) +
    xlab(expression(theta[1])) +
    ylab(expression(theta[3]))
}

# Pragmática para cada dissimilaridade, p = 0.5
simple_grid_a = color_grid %>% generate_grid(ABS_diss, p0 = 0.5, eps = 0.1)
ga <- plot_grid(simple_grid_a)

simple_grid_l <- color_grid %>% generate_grid(L2_diss, p0 = 0.5, eps = 0.1)
gl <- plot_grid(simple_grid_l)

simple_grid_m <- color_grid %>% generate_grid(MAX_diss, p0 = 0.5, eps = 0.1)
gm <- plot_grid(simple_grid_m)

ggarrange(ga + labs(title = "L1"), gl + labs(title = "Euclidiana"),
          gm + labs(title = "Máximo"), ncol = 3)
  

#Grade para p livre
generate_join_grid = function(color_grid, dissimilarity, eps){
  #Pragmática para a primeira combinação de thetas
  join_grid = color_grid %>% generate_grid(dissimilarity, thetas[1], eps)
  for(theta in thetas){
    new_grid = color_grid %>% generate_grid(dissimilarity, theta, eps)
    #Atualiza a pragmática se a nova grade contém uma combinação na pragmática
    join_grid$color = (join_grid$color | new_grid$color)
  }
  return(join_grid)
}

# Pragmática para cada dissimilaridade, p livre
join_grid_a = color_grid %>% generate_join_grid(ABS_diss, eps = 0.1)
ga <- plot_grid(join_grid_a)

join_grid_l = color_grid %>% generate_join_grid(L2_diss, eps = 0.1)
gl <- plot_grid(join_grid_l)

join_grid_m = color_grid %>% generate_join_grid(MAX_diss, eps = 0.1)
gm <- plot_grid(join_grid_m)

ggarrange(ga + labs(title = "L1"), gl + labs(title = "Euclidiana"),
          gm + labs(title = "Máximo"), ncol = 3)

# Plot para várias pragmáticas simultaneamente
plot_grid = function(grid_list, col_list = "lightblue", label = NULL,
                     alpha = 0.9){
  poly_grid <- data.frame()
  
  for(i in 1:length(grid_list)){
    grid <- grid_list[[i]]
    grid %<>% filter(color) %>%
      group_by(theta1)
    poly_min = grid %>% summarise(theta3 = min(theta3)) %>%
      mutate(group_id = i)
    poly_max = grid %>% 
      summarise(theta3 = max(theta3)) %>%
      arrange(desc(theta1)) %>%
      mutate(group_id = i)
    poly_grid <- bind_rows(poly_grid, poly_min, poly_max)
  }
  g <- poly_grid %>%
    ggplot(aes(x = theta1, y = theta3)) +
    theme_minimal() +
    geom_polygon(color = NA, aes(fill = as.factor(group_id)), alpha = alpha) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 0)) + 
    geom_segment(aes(x = 0, xend = 0, y = 0, yend = 1)) +
    geom_segment(aes(x = 0, xend = 1, y = 1, yend = 0)) +
    stat_function(fun = function(x) (1-sqrt(x))^2, n = 101,
                  linetype="dashed", color = "black") +
    xlab(expression(theta[1])) +
    ylab(expression(theta[3]))
  if(length(col_list) == length(grid_list)){
    if(!is.null(label)){
      g <- g + scale_fill_manual(name = label[1], labels = label[-1],
                                   values = col_list)
    } else g <- g + scale_fill_manual(values = col_list)
  }
  g
}
tema <- theme(legend.position = c(.954, 0.645),
              legend.justification = c("right", "bottom"),
              legend.box.just = "right",
              legend.margin = margin(6, 6, 6, 6),
              legend.box.background = element_rect(colour = "black", size = .5),
              legend.text.align = 0,
              legend.title = element_text(size = 30, face = "bold"),
              legend.text = element_text(size=25),
              axis.text = element_text(size = 20),
              axis.title = element_text(size = 25))

# p = 0.5
gs <- list(simple_grid_m, simple_grid_l, simple_grid_a) %>%
  plot_grid(col_list = c("lightblue1", "lightblue2", "lightblue3"),
            label = c("Dissimilaridade", expression(Dist. ~ L[infinity]),
                      expression(Dist. ~ L[2]),
                      expression(Dist. ~ L[1])), alpha = .7) + tema
# p livre
gl <- list(join_grid_m, join_grid_l, join_grid_a) %>%
  plot_grid(col_list = c("lightblue1", "lightblue2", "lightblue3"),
            label = c("Dissimilaridade", expression(Dist. ~ L[infinity]),
                      expression(Dist. ~ L[2]),
                      expression(Dist. ~ L[1])), alpha = .7) + tema

ggarrange(gs + labs(title = expression(bold(H[0]: p == 0.5))),
          gl + labs(title = expression(bold(H[0]: p ~ livre))))


################################################################################
####################### Modelagem de theta via Dirichlet #######################
################################################################################

# Verossimilhança: X|theta ~ Multinomial(n, theta)
# Priori: theta ~ Dirichlet(alpha), alpha = (1, 1, 1)
# Posteriori: theta|X ~ Dirichlet(X + alpha)

# Parâmetros da posteriori da Dirichlet
posterior = function(data, prior_par = c(1, 1, 1)) data + prior_par

# Objetivo: obter o HPD de x% da posteriori
# HPD: menor região tal que P(theta no HPD) = x%
## 1) Identificar o valor da densidade tal que (1-x)% da amostra será maior 
hpd_post_cut = function(posterior_par, probability = 0.95, B = 10^5){
  quantile = posterior_par %>%
    rdirichlet(n = B, alpha = .) %>%
    ddirichlet(posterior_par, log = TRUE) %>%
    quantile(probs = 1 - probability) %>%
    as.numeric()
  return(quantile)
}
## 2) Obter o HPD para cada experimento
hpd_labels = tibble(idx = 1:nrow(data), 
                    x = rep(NA, nrow(data)), 
                    y = rep(NA, nrow(data)))
hpd_chull_grid = NULL
for(ind in 1:nrow(data)){
  counts = data[ind,]
  post = counts %>% posterior() %>% as.numeric() #Posteriori do experimento
  quant = hpd_post_cut(post, probability = 0.95) #Valor para o HPD
  this_hpd_grid = color_grid %>% #Valores no grid maiores que a referência
    mutate(color = (ddirichlet(x = cbind(theta1, theta2, theta3),
                               alpha = post, log = T) > quant)) %>%
    filter(color) %>% #Apenas os valores dentro do HPD
    select(theta1, theta3)
  this_chull = chull(this_hpd_grid) #Fecho convexo do HPD
  this_hpd_chull = this_hpd_grid[this_chull,] %>% #Pontos que formam o fecho
    mutate(idx = ind)
  hpd_chull_grid <- bind_rows(hpd_chull_grid,
                              this_hpd_chull) #Adiciona o HPD à lista
  hpd_labels[ind,] = list(idx = ind, x = (post/sum(post))[1],
                         y = (post/sum(post))[3]) #Média a posteriori de theta
}
hpd_chull_grid <- hpd_chull_grid %>%
  mutate(idx = as.factor(idx)) #Transformando índice em fator


################################################################################
######################## Teste hip. pragmática pelo HPD ########################
################################################################################

g1 <- gl
all_idx = unique(hpd_chull_grid$idx) #Índices dos experimentos
for(ind in all_idx){
  fecho = hpd_chull_grid %>% filter(idx == ind)
  g1 <- g1 + geom_polygon(aes(x = theta1, y = theta3),
                            fill = ind, col = ind,
                        data = fecho, alpha = .5)
}
g1 <- g1 + geom_text(aes(x = x, y = y, label = idx), data = hpd_labels)

g2 <- gs
all_idx = unique(hpd_chull_grid$idx) #Índices dos experimentos
for(ind in all_idx){
  fecho = hpd_chull_grid %>% filter(idx == ind)
  g2 <- g2 + geom_polygon(aes(x = theta1, y = theta3),
                        fill = ind, col = ind,
                        data = fecho, alpha = .5)
}
g2 <- g2 + geom_text(aes(x = x, y = y, label = idx), data = hpd_labels)

ggarrange(g1, g2)
