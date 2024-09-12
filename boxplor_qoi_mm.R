# install.packages("reshape2")
library(reshape2) 

# Caminho para ler o arquivo de entrada
data <- read.csv("C:/Users/Lucas Teixeira/Desktop/Sobol/2024-08-29_18-30-43.728964_1M_result_30_50_hcl_emualtor.csv", header=FALSE)

# Nomeando as colunas do dataframe
names(data) <- c("u_u", "theta_v", "tau_v_plus", "tau_w_plus", "tau_fi", "tau_so1", "tau_so2", "tau_si", 
                 "minAPD", "maxAPD", "MAXmaxdVdt", "MINmaxdVdt", "maxDerivada", "maxDerivada2")

# Nomes das quantidades de interesses que serão feitos os boxplots
yNames <- c("minAPD", "maxAPD", "MAXmaxdVdt", "MINmaxdVdt", "maxDerivada", "maxDerivada2")

# Retirando a primeira linha dos dados, que é uma string com os nomes
data <- data[-1,]

# Transformando os dados em dados numéricos
data <- apply(data, 2, as.numeric)



###########################################################
#                INÍCIO DO GRÁFICO INDIVIDUAL             #
###########################################################

# Definindo o grip da área de plotagem
par(mfrow = c(3, 2))

# Recursão entre que gerará os boxplots do vetor das QoIs (yNames)
lapply(1:length(yNames), function(i) {
  # Extrair os dados da coluna
  values <- data[,yNames[i]]
  
  # Calcular o boxplot (sem plotar) para obter os limites
  bp <- boxplot.stats(values)
  
  # Contar o número de outliers
  num_outliers <- length(bp$out)
  
  # Calcular o percentual de outliers
  total_points <- length(values)
  outlier_percent <- (num_outliers / total_points) * 100
  
  # Criar o boxplot com personalizações
  boxplot(values, # Valores para gerar o boxplot
          main = paste(yNames[i], "\nOutliers:", num_outliers, "(", round(outlier_percent, 2), "%)"), # Título
          col = "darkGrey", # Cor principal do boxplot
          border = "black", # Cor das linhas e bordas do boxplot
          cex.main = 1.5, # Tamanho da fonte do título
          cex.axis = 1.2, # Tamanho da fonte dos eixos
          outpch = 16, # Formato do dado marcado como outlier
          outcol = "red", # Cor principal do outlier
          ylim = range(values, na.rm = TRUE) # Limites dos gráficos
  )
})


###########################################################
#                 INÍCIO DO GRÁFICO CONJUNTO              #
###########################################################

# Padronizando os dados (todos seguem uma distribuição Normal Padrão, também chamada de distribuição Z)
data_zscore_normalized <- scale(data)

# Gerando um tipo de dataframe diferentes para permitir uma plotagem única das seis QoIs
data_mod <- melt(data_zscore_normalized[,9:14], measure.vars=yNames) 

# Criando vetor para armazenar a quantidade de outlier das QoI
outlier_percentages <- numeric(length(yNames))

# Calculando a porcentagem de outlier de cada QoI
outlier_percentages <- sapply(1:length(yNames), function(i) {
  values <- data[, yNames[i]]
  
  bp <- boxplot.stats(values)
  
  outlier_percent <- (length(bp$out) / length(values)) * 100
  
  outlier_percent
})

# Criando os nomes do eixo X de cada gráfico 
labels = sapply(1:length(outlier_percentages), function(i) {
  paste(yNames[i], "\nOutliers:" , round(outlier_percentages[i], 3), "%", sep = "")
})

# Definindo o grip da área de plotagem com um gráfico único
par(mfrow = c(1,1))

boxplot(value ~ Var2, data = data_mod,
  main = "Z-SCORED QOI", # Título
  col = "darkGrey", # Cor principal do boxplot
  border = "black", # Cor das linhas e bordas do boxplot
  cex.main = 1.5, # Tamanho da fonte do título
  cex.axis = 1, # Tamanho da fonte dos eixos
  outpch = 16, # Formato do dado marcado como outlier
  outcol = "red", # Cor principal do outlier
  names = labels, # Nomes de cada gráfico do eixo X
  xlab = "Variables", # Título do eixo X
  ylab = "Values" # Título do eixo Y
)
