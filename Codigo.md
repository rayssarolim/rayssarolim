```r
rm(list=ls(all=T))

library(readxl)
library(dplyr)
library(forecast)
library(ggplot2)
library(randtests)
library(lubridate)
library(tseries)
library(FinTS)
library(gridExtra)
library(rugarch)
library(skedastic)
library(rugarch)
library(purrr)

# --- DADOS ---
dados_diarios <- read_excel("10628227/dados diarios.xlsx", 
    col_types = c("numeric", "numeric", "text", "numeric", "text"))

dados_diarios$Data <- as.Date(dados_diarios$Data)

ibovespa_diario_ts <- ts(dados_diarios$Valor, start = c(2020, 1), 
    frequency = 252)
log_precos <- log(ibovespa_diario_ts)

# --- ANÁLISE EXPLORATÓRIA DA SÉRIE ---
autoplot(log_precos) + geom_line(linewidth=0.5) + 
  labs(y="Log-Ibovespa", x="Período") + theme_bw(base_size=14) +
  theme(axis.text=element_text(face="bold")+
  panel.border=element_rect(linewidth=1.5))

# ---  ACF e PACF da série  ---
par(mfrow = c(1, 2))
acf(log_precos, lag.max = 40, main="Função de Autocorrelação (ACF)")
pacf(log_precos, lag.max = 40, main="Autocorrelação Parcial (PACF)")

meu_tema <- theme_bw() +
  theme(axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 14))

# ---  Teste de Tendência ---
cox.stuart.test(log_precos)       #H0: Não existe tendência

# --- Testes de Raízes Unitárias (Estacionariedade) ---
# Teste Augmented Dickey-Fuller (ADF)
adf.test(log_precos)

# Teste Phillips-Perron (PP)
pp.test(log_precos)

# Teste KPSS (para estacionariedade em nível)
kpss.test(log_precos, null = "Level")

# Teste KPSS (para estacionariedade em tendência)
kpss.test(log_precos, null = "Trend")

# --- ANÁLISE DA SÉRIE DIFERENCIADA ---
serie_dif <- diff(log_precos) # série diferenciada

plot_dif <- autoplot(serie_dif),labs(title = "",
       y = "Log-Ibovespa", x = "Tempo") + theme_bw()

# ACF e PACF da série diferenciada
par(mfrow = c(1, 2))
acf(serie_dif, main = "ACF dos Log-Retornos", lag.max = 40)
pacf(serie_dif, main = "PACF dos Log-Retornos", lag.max = 40)
par(mfrow = c(1, 2))

#  ACF
acf <- ggAcf(serie_dif,lag.max = 40),labs(y ="Autocorrelação",x ="Defasagens")+
  ggtitle(""), theme_bw()

#  PACF
pacf <- ggPacf(serie_dif, lag.max = 40) + 
  labs(y = "Autocorrelação Parcial", x = "Defasagens") +
  ggtitle("")+ theme_bw()

# teste sazonalidade
log_retornos_vetor <- as.numeric(diff(log_precos))
dados_retornos <- data.frame(Data = tail(dados_diarios$Data, -1),
  LogRetorno = log_retornos_vetor) %>%
  mutate( DiaDaSemana = wday(Data, label = TRUE, week_start = 1),
    Mes = month(Data, label = TRUE))
    
# H0: Não há diferença de retorno entre os dias da semana.
kruskal.test(LogRetorno ~ DiaDaSemana, data = dados_retornos)

# AJUSTE DO MODELO ARIMA PARA A MÉDIA
modelo_media <- auto.arima(log_precos, trace = TRUE)
summary(modelo_media)

# DIAGNÓSTICO DE RESÍDUOS DO MODELO ARIMA
residuos_arima <- residuals(modelo_media)
residuos_df <- data.frame(Residuos = as.vector(residuos_arima))

# Teste de Autocorrelação (Ljung-Box)
checkresiduals(modelo_media) # H0: Resíduos são independentes.

# Teste de Heterocedasticidade para os resíduos do Modelo auto.arima
teste_engle_auto <- ArchTest(residuos_arima, lags = 12)

# teste de white
modelo_lm_auto <- lm(residuos_arima ~ fitted(modelo_media))
teste_white_bptest_auto <- bptest(modelo_lm_auto, ~ fitted(modelo_lm_auto) + 
    I(fitted(modelo_lm_auto)^2))

cf<- acf(residuos_arima, main = "ACF dos Log-Retornos", lag.max = 40)
pa <- pacf(residuos_arima, main = "PACF dos Log-Retornos", lag.max = 40)

# Teste de Heterocedasticidade (ARCH) / p-valor baixo, temos efeitos ARCH/GARCH.
FinTS::ArchTest(residuos_arima, lag = 12) # H0: Não há efeitos ARCH.
checkresiduals(residuos_arima) # H0: Os resíduos são "ruído branco"

teste_ljung_box <- Box.test(residuos_arima^2, lag = 10, type = "Ljung-Box")

#  ACF
ggAcf(residuos_arima,lag.max = 40),labs(y ="Autocorrelação", x ="Defasagens")+
  ggtitle(""),theme_bw()

#  PACF
ggPacf(residuos_arima, lag.max = 40) + 
labs(y = "Autocorrelação Parcial", x = "Defasagens") +
  ggtitle(""), theme_bw()

# Gráficos de Normalidade
plot_hist <- ggplot(residuos_df, aes(x = Residuos)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "white", 
  color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(residuos_df$Residuos), 
  sd = sd(residuos_df$Residuos)), color = "red", linewidth = 1) +
  labs(title = "  ", x = "Resíduos", y = "Densidade") + 
  theme_bw()

# QQ-Plot de Normalidade
plot_qq <- ggplot(residuos_df, aes(sample=Residuos)) + stat_qq() + 
           stat_qq_line(col="red", linewidth=1) + theme_bw() + 
           labs(x="Q. Teóricos", y="Q. Amostra")+
grid.arrange(plot_hist, plot_qq, ncol=2)

# Teste de Normalidade
shapiro.test(residuos_arima) #H0: Resíduos são normais

# Modelo ARIMA(2,1,2) 
# O PACF sugere p=2 e o ACF sugere q=2.
modelo_212 <- Arima(log_precos, order = c(2, 1, 2))
summary(modelo_212)

# RESÍDUOS Modelo ARMA(2,2)

residuos_212 <- residuals(modelo_212)
residuos_dff <- data.frame(Residuos = as.vector(residuos_212))

# Teste de Heterocedasticidade para os resíduos do Modelo ARIMA(2,1,2)
teste_engle_212 <- ArchTest(residuos_212, lags = 12)
# teste de white
modelo_lm_212 <- lm(residuos_212 ~ fitted(modelo_212))
teste_white_bptest_212 <- bptest(modelo_lm_212, ~ fitted(modelo_lm_212) +
    I(fitted(modelo_lm_212)^2))

# Teste de Autocorrelação (Ljung-Box)
checkresiduals(modelo_212) # H0: Resíduos são independentes.

acf_212<- acf(residuos_212, main = "ACF dos Log-Retornos", lag.max = 40)
pacf_212 <- pacf(residuos_212, main = "PACF dos Log-Retornos", lag.max = 40)

# Teste de Heterocedasticidade (ARCH)
FinTS::ArchTest(residuos_212, lag = 12) # H0: Não há efeitos ARCH.

checkresiduals(residuos_212)

teste_ljung_box <- Box.test(residuos_212^2, lag = 10, type = "Ljung-Box")

#  ACF
ggAcf(residuos_212, lag.max = 40), labs(y = "Autocorrelação", x = "Defasagens") +
  ggtitle(""),theme_bw()

#  PACF
ggPacf(residuos_212, lag.max = 40) + 
  labs(y = "Autocorrelação Parcial", x = "Defasagens"),ggtitle(""),theme_bw()

# Gráficos de Normalidade
hist <- ggplot(residuos_dff, aes(x = Residuos)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "white" + 
  color = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(residuos_df$Residuos),+
  sd = sd(residuos_df$Residuos)), color = "red", linewidth = 1) +
  labs(title = "  ", x = "Resíduos", y = "Densidade") + 
  theme_bw()

# QQ-Plot de Normaliresiduos_df# QQ-Plot de Normalidade
QQ <- ggplot(residuos_dff, aes(sample=Residuos)) + stat_qq() + 
      stat_qq_line(col="red", linewidth=1) + theme_bw() + 
      labs(x="Teóricos", y="Amostra")+
grid.arrange(hist, QQ, ncol=2)

# Teste de Normalidade
shapiro.test(residuos_arima) #H0: Resíduos são normais

# 3. AJUSTE DO MODELO ARIMA-GARCH
ordem_arima <- arimaorder(modelo_media)
ordem_arma <- c(ordem_arima["p"], ordem_arima["q"]) # Pegamos apenas p e q

# Especificar o modelo GARCH.
spec_garch <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = ordem_arma), distribution.model = "std" )
serie_dif <- diff(log_precos, differences = ordem_arima["d"])
modelo_final_garch <- ugarchfit(spec = spec_garch, data = serie_dif)

# 4. DIAGNÓSTICO DO MODELO GARCH FINAL
# Resíduos padronizados
residuos_padronizados_garch <- residuals(modelo_final_garch, standardize = TRUE)

# Teste de Efeitos ARCH nos resíduos do GARCH
FinTS::ArchTest(as.numeric(residuos_padronizados_garch), lag = 12)

# Gráficos de diagnóstico do próprio pacote rugarch
plot(modelo_final_garch, which = "all")
# RESULTADO: efeitos ARCH NÃO SUMIRAM

# --- Tentando o modelo GJR-GARCH ---
spec_garch_gjr <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)), 
  mean.model = list(armaOrder = ordem_arma), distribution.model = "std")
modelo_final_gjr <- ugarchfit(spec = spec_garch_gjr, data = serie_dif)

# diagnóstico para o novo modelo
print(modelo_final_gjr)
plot(modelo_final_gjr, which="all")
residuos_gjr <- residuals(modelo_final_gjr, standardize = TRUE)
FinTS::ArchTest(as.numeric(residuos_gjr), lag = 12)
# RESULTADO: efeitos ARCH SUMIRAM

# AJUSTE E DIAGNÓSTICO DO MODELO 2,1,2 ARIMA-GARCH
# A ordem da média agora é (2,2) e vamos manter o gjrGARCH(1,1) para a variância
spec_garch_212 <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(2, 2)), # <- A única mudança está aqui!
  distribution.model = "std")

# Ajuste do novo modelo GARCH
modelo_garch_212 <- ugarchfit(spec = spec_garch_212, data = serie_dif)
print(modelo_garch_212)

# 2. DIAGNOSTICAR O NOVO MODELO
residuos_212_padronizados <- residuals(modelo_garch_212, standardize = TRUE)

# Teste ARCH-LM. O objetivo é que o p-valor seja ALTO (maior que 0.05)
FinTS::ArchTest(as.numeric(residuos_212_padronizados), lag = 12)
plot(modelo_garch_212, which = "all")

# GRÁFICOS PARA O MODELO 1: ARIMA(1,1,3)-GJR(1,1)
residuos_modelo1 <- as.numeric(residuals(modelo_final_gjr, standardize = TRUE))
par(mfrow = c(3, 1), mar = c(4.5, 4.5, 1.5, 1.5)) 

# Gráfico 1: Série dos Resíduos Padronizados vs. Tempo
plot(residuos_modelo1, type = 'l', main = "", xlab = "Tempo", 
     ylab = "Resíduos Padronizados", col = "darkgray")
abline(h = 0, col = "black", lty = "dotted")

# Gráfico 2: ACF dos Resíduos Padronizados
acf_obj1 <- acf(residuos_modelo1, lag.max = 24, plot = FALSE)
plot(acf_obj1$lag[-1], acf_obj1$acf[-1], type = 'h', main = "", 
     xlab = "Defasagens", ylab = "ACF", lwd = 5, lend = "butt")
abline(h = c(-1.96/sqrt(length(residuos_modelo1)), 
    1.96/sqrt(length(residuos_modelo1))), col = "red", lty = "dashed")
abline(h = 0)

# Gráfico 3: P-valores do Teste Ljung-Box (nos resíduos ao quadrado)
lags_a_testar <- 1:24
p_valores1 <- purrr::map_dbl(lags_a_testar, function(h) {
  teste <- Box.test(residuos_modelo1^2, lag = h, type = "Ljung-Box")
  return(teste$p.value)})

plot(lags_a_testar, p_valores1, main = "",  xlab = "Defasagens", 
    ylab = "P-valor", ylim = c(0, 1), pch = 19, col = "blue")
abline(h = 0.05, col = "red", lty = "dashed", lwd = 2)

# GRÁFICOS PARA O MODELO 2: ARIMA(2,1,2)-GJR(1,1)
residuos_modelo2 <- as.numeric(residuals(modelo_garch_212, standardize = TRUE))
par(mfrow = c(3, 1), mar = c(4.5, 4.5, 1.5, 1.5))

# Gráfico 1: Série dos Resíduos Padronizados vs. Tempo
plot(residuos_modelo2, type = 'l', main = "", xlab = "Tempo", 
    ylab = "Resíduos Padronizados", col = "darkgray")
abline(h = 0, col = "black", lty = "dotted")

# Gráfico 2: ACF dos Resíduos Padronizados
acf_obj2 <- acf(residuos_modelo2, lag.max = 24, plot = FALSE)
plot(acf_obj2$lag[-1], acf_obj2$acf[-1], type = 'h', main = "", 
    xlab = "Defasagens", ylab = "ACF", lwd = 5, lend = "butt")
abline(h = c(-1.96/sqrt(length(residuos_modelo2)), 
    1.96/sqrt(length(residuos_modelo2))), col = "red", lty = "dashed")
abline(h = 0)

# Gráfico 3: P-valores do Teste Ljung-Box (nos resíduos ao quadrado)
p_valores2 <- purrr::map_dbl(lags_a_testar, function(h) {
  teste <- Box.test(residuos_modelo2^2, lag = h, type = "Ljung-Box")
  return(teste$p.value)})
  
plot(lags_a_testar, p_valores2, main = "", xlab = "Defasagens", ylab = "P-valor", 
    ylim = c(0, 1), pch = 19, col = "blue")
abline(h = 0.05, col = "red", lty = "dashed", lwd = 2)
par(mfrow = c(1, 1))

# --- Teste para o Modelo 1 GARCH ---
residuos_padronizados_1 <- as.numeric(residuals(modelo_final_gjr, 
    standardize = TRUE))
teste_shapiro_garch_1 <- shapiro.test(residuos_padronizados_1)

# --- Teste para o Modelo 2 GARCH ---
residuos_padronizados_2 <- as.numeric(residuals(modelo_garch_212, 
    standardize = TRUE))
teste_shapiro_garch_2 <- shapiro.test(residuos_padronizados_2)
