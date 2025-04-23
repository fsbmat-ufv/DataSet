rm(list = ls())
cat("\014")
library(tidyverse)
library(sampleSelection)

# 0) Ler o CSV
df <- read.csv("ess_dive/data/ForestGEO1.csv",
               na.strings = c("", "NA"))

# 1) Conversões essenciais (note hom)
df <- df %>%
        mutate(
                dbh.full.census      = log(as.numeric(dbh.full.census)),
                meanWD               = as.numeric(meanWD),
                hom                  = as.numeric(hom),            # coluna correta
                H_considering_damage = as.numeric(H_considering_damage),
                weights.ind          = as.numeric(weights.ind),
                site                 = factor(site),
                sel                  = status == "A",
                logH                 = log1p(H_considering_damage)
        )

# 2) Filtrar **apenas** as linhas que serão usadas nas fórmulas
df_mod <- df %>% 
        filter(
                !is.na(sel),
                !is.na(dbh.full.census),
                !is.na(meanWD),
                !is.na(hom),
                !is.na(site),
                !is.na(weights.ind)      # <-- garante que pesos e amostra tenham o mesmo n
        )

# 3) Verificação rápida (opcional)
stopifnot(nrow(df_mod) == sum(!is.na(df_mod$weights.ind)))

# 4) Ajuste do Heckman — agora pesos e dados têm o mesmo tamanho
model <- heckit(
        selection = sel ~ dbh.full.census + meanWD + hom,
        outcome   = logH ~ dbh.full.census + meanWD + site,
        data      = df_mod,
        weights   = df_mod$weights.ind         # passe o vetor já filtrado
)

summary(model)
str(df)
unique(df$hom)
quantile(df_mod$weights.ind, c(.01, .25, .5, .75, .99))
library(lubridate)   # facilita datas

df_std <- df_mod %>%                       # comece do objeto já limpo
        mutate(
                date.ams         = ymd(date.ams),
                date.full.census = ymd(date.full.census),
                time_since_census = as.numeric(date.ams - date.full.census, units = "days")/365.25,
                log_dbh          = scale(log(dbh.full.census))
        ) %>% 
        filter(!is.na(time_since_census))        # descarta linhas sem as duas datas
heck_ml <- selection(
        selection = sel ~ log_dbh + meanWD + hom + time_since_census,
        outcome   = logH   ~ log_dbh + meanWD + site,
        data      = df_std,
        weights   = df_std$weights.ind,  # <- esta linha garante que o vetor tenha o tamanho correto
        method    = "ml"
)


summary(heck_ml, part = "full")

library(ssmodels)

# 1. Filtrar apenas os dados com H > 0 (BS requer resposta positiva)
df_bs <- df_std %>%
        filter(!is.na(H_considering_damage),
               H_considering_damage > 0)

# 2. Fórmulas
selectEq <- sel ~ log_dbh + meanWD + hom + time_since_census
outcomeEq <- H_considering_damage ~ log_dbh + meanWD + site

# 3. Ajustar modelo
selectEq <- sel ~ log_dbh
outcomeEq <- H_considering_damage ~ log_dbh

modelo_teste <- HeckmanBS(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_bs
)

modelo_bs <- HeckmanBS(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_bs
)

# 4. Resultados
str(modelo_bs)      # estrutura geral
modelo_bs$Coefficients
modelo_bs$loglik
modelo_bs$aic


####
summary(df_bs$H_considering_damage)
range(df_bs$H_considering_damage, na.rm = TRUE)
any(!is.finite(df_bs$H_considering_damage))         # deve ser FALSE
any(is.na(df_bs$H_considering_damage))              # deve ser FALSE
dim(df_bs)
df_bs <- df_bs %>%
        filter(H_considering_damage > 0.01)
summary(df_bs$log_dbh)
any(!is.finite(df_bs$log_dbh))      # deve ser FALSE

modelo_bs <- HeckmanBS(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_bs,
        start     = rep(0.1, ncol(model.frame(selectEq, df_bs)) +
                                ncol(model.frame(outcomeEq, df_bs)) + 1)  # +1 para rho
)

library(sampleSelection)

# Modelo Heckman com normalidade (ml)
modelo_normal <- selection(
        selection = sel ~ log_dbh,
        outcome   = H_considering_damage ~ log_dbh,
        data      = df_bs,
        method    = "ml"
)

summary(modelo_normal, part = "full")
