rm(list = ls())
cat("\014")
library(tidyverse)
library(sampleSelection)
library(ssmodels)

# Lê o CSV com atenção para strings vazias ou lixo
df <- read.csv("ess_dive/data/ForestGEO1.csv",
               na.strings = c("", "NA", "?", "NULL", "null", "–", " ", "-", "N/A"))
# Quais valores únicos aparecem?
unique(df$dbh.full.census[!grepl("^[0-9\\.]+$", df$dbh.full.census)])
df <- df %>%
        mutate(
                dbh.full.census = as.character(dbh.full.census),
                dbh.full.census.num = suppressWarnings(as.numeric(dbh.full.census))
        )
sum(is.na(df$dbh.full.census.num))   # total de NAs após conversão
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


# 4) Ajuste do Heckman — agora pesos e dados têm o mesmo tamanho
selectEq <- sel ~ dbh.full.census + hom
outcomeEq <- logH ~ dbh.full.census + meanWD + site
model <- heckit(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_mod,
        method = "ml")
summary(model)
model <- HeckmanCL(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_mod
)


# 1. Filtrar apenas os dados com H > 0 (BS requer resposta positiva)
df_bs <- df_mod %>%
        mutate(
                log_dbh = log(dbh.full.census.num),
                H_considering_damage = ifelse(sel, H_considering_damage, 0)
        ) %>%
        filter(
                !is.na(H_considering_damage),
                H_considering_damage >= 0,
                is.finite(log_dbh),
                sel %in% c(TRUE, FALSE)
        )

df_bs <- df_bs %>%
        mutate(
                date.ams = as.Date(date.ams),
                date.full.census = as.Date(date.full.census),
                time_since_census = as.numeric(date.ams - date.full.census)
        )

df_bs <- df_bs %>%
        filter(
                !is.na(log_dbh),
                !is.na(meanWD),
                !is.na(hom),
                !is.na(site),
                !is.na(time_since_census),
                !is.na(H_considering_damage),
                is.finite(H_considering_damage),
                is.finite(log_dbh),
                H_considering_damage >= 0,
                sel %in% c(TRUE, FALSE)
        )
df_bs$sel <- 1*(df_bs$sel)
# 2. Fórmulas
selectEq <- sel ~ dbh.full.census + meanWD + hom
outcomeEq <- logH ~ dbh.full.census + meanWD + site



modelo_bs <- HeckmanCL(
        selection = selectEq,
        outcome   = outcomeEq,
        data      = df_bs
)

table(df_bs$sel, useNA = "ifany")
table(1*(df_bs$H_considering_damage<0))
table(1*is.na(df_bs$H_considering_damage))
lm_test <- lm(H_considering_damage ~ log_dbh + meanWD + site, 
              data = df_bs[df_bs$sel == TRUE, ])
summary(lm_test)
