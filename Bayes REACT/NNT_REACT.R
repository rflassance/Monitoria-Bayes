################################################################################
#################______________________________________________#################
#################                                              #################
#################  Exemplo do REACT bayesiano para metanálise  #################
#################______________________________________________#################
#################                                              #################
################################################################################

# Carregando bibliotecas
library(dplyr)
library(ggplot2)
library(coda)

################################################################################
######################## Carregando dados de metanálise ########################
################################################################################

# Contexto: Estudo de aderência de pacientes a abordagens para combate ao fumo
## studlab: Referência do estudo realizado
## evente, eventc: Número de sucessos dos grupos de tratamento e controle
## ne, nc: Número de pacientes dos grupos de tratamento e controle
## mean: Estimativa pontual da diferença entre as médias dos grupos (e - c)
## study_class: Tipo do estudo
###             1 - Tecnologias de saúde com assistência computadorizada
###             2 - Adição do uso de fármacos ao tratamento
meta <- read.csv("Bayes REACT/meta.csv", header = T)
View(meta)


################################################################################
###################### Modelagem da diferença entre médias #####################
################################################################################

# Verossimilhança: Yi ~ Ber(theta)
# Priori: theta ~ Beta(a0, b0)
# Posteriori: theta | Y ~ Beta(sum(Yi) + a0, n - sum(Yi) + b0)
# Priori de Jeffreys: a0 = 1/2, b0 = 1/2

# Parâmetros da posteriori
alphae <- meta$evente + 1/2
betae <- meta$ne - meta$evente + 1/2
alphac <- meta$eventc + 1/2
betac <- meta$nc - meta$eventc + 1/2

# Simulação da posteriori da diff de médias via Monte Carlo
N <- 10000
diff_sim <- matrix(NA, ncol = dim(meta)[1], nrow = N)
colnames(diff_sim) <- meta$studlab

set.seed(42)
for(i in 1:dim(meta)[1]){
  diff_sim[,i] <- rbeta(N, alphae[i], betae[i]) - rbeta(N, alphac[i], betac[i])
}

diff_MC <- mcmc(diff_sim)

densplot(diff_MC)

HPD_int <- HPDinterval(diff_MC, prob = 0.95)
HPD_int #Que conclusão tomar? Tratamento difere do controle?


################################################################################
#################### REACT bayesiano a partir da medida NNT ####################
################################################################################

# NNT: Number Necessary to Treat, NNT = 1/[P(tratamento) - P(controle)]
## Número médio de pacientes para que tratamento se mostre melhor que controle
## Referência: https://thennt.com/nnt/varenicline-smoking-cessation/

NNT <- c(-10, 6)
# True positive: P(tratamento) - P(controle) > 1/6 
# True negative: P(tratamento) - P(controle) < -1/10


meta <- bind_cols(meta, HPD_int)

#Identificar decisão para cada caso
regiao <- colSums(
  apply(meta[c("lower", "upper")], 1, function(x) findInterval(1/NNT, x))
  )
cols <- ifelse(regiao %in% c(1,3), 1/2,
               ifelse(regiao %in% c(0,4), 1, 0))
cols <- factor(cols, levels = c(0, 1/2, 1),
               labels = c("Aceitar", "Agnóstico", "Rejeitar"))

meta_data <- data.frame(CIlower = meta$lower,
                        CIupper = meta$upper,
                        points = meta$mean,
                        color = cols,
                        study_class = factor(meta$study_class,
                                             levels = c(1,2),
                                             labels = c("Seguimento",
                                                        "Farmacoterapia")),
                        studlab = factor(meta$studlab,
                                         levels = rev(sort(unique(meta$studlab)))
                                         )
                        )

meta_data %>%
  ggplot(aes(y = studlab, col = color)) +
  geom_errorbarh(aes(xmin = CIlower, xmax = CIupper), height=.2) +
  geom_point(ggplot2::aes(x = points)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate('rect', xmin = 1/NNT[1], xmax = 1/NNT[2], ymin = 0,
           ymax = nrow(meta_data), alpha=.2, fill='dodgerblue3') +
  scale_color_manual(values=c("Aceitar" = "darkgreen",
                                       "Agnóstico" = "goldenrod",
                                       "Rejeitar" = "darkred")) +
  theme_bw() +
  labs(title = "(Forest plot) REACT baseado no NNT",
                x = "Diferença de médias", y = "Estudo", color = "Decisão") +
  facet_grid(~study_class, scales="free")

################################################################################