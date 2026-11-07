#ANALYSIS 1 Comparing SR across Arousal State with one minute window
# Recode arousal labels throughout the analysis

library(lmerTest)
library(dplyr)
library(broom.mixed)
library(tidyverse)
library(reshape2)


recode_arousal <- function(df) {
  df %>%
    mutate(across(
      any_of(c("Arousal_State", "arousal_state", "arousal")),
      ~ dplyr::recode(
        as.character(.x),
        "drowsy" = "Low Arousal",
        "Low Arousal" = "Low Arousal",
        "alert" = "High Arousal",
        "High Arousal" = "High Arousal"
      )
    ))
}

#ANALYSIS 1 Networks
df_net<-read.csv("Net_Flex_MSTR_one_minute.csv")
df_net <- recode_arousal(df_net)

# Variables to test
vars <- c("a_sal", "p_sal", "d_dmn", "v_dmn", "r_cen", "l_cen")

test<-lmer(a_sal ~ Arousal_State + (1 | ID), data = df_net)

summary(test)

# Function to fit model and extract results
run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ Arousal_State + (1 |ID)")
  )
  
  model <- lmer(formula, data = df_net)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "Arousal_StateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#ANALYSIS 2 Parcels

df_parcel<- read.csv("Net_Parcel_MSTR_1min.csv")
all_cols <- colnames(df_parcel)

parcels <- all_cols[grep("dmn|sal|cen", all_cols)]
parcels <- parcels[2:51]

vars<- parcels

run_lmm_parcel <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ Arousal_State + (1 |ID)")
  )
  
  model <- lmer(formula, data = df_parcel)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "Arousal_Statedrowsy") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}


# Run all models
results <- bind_rows(lapply(vars, run_lmm_parcel))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"HCP-7T-SR-parcels-11-6-26-one-minute.csv")


#ANALAYSS 3 Community Assignment

community_assign_mat_fin<-read.csv("HCP_7T_community_allegience_one_min.csv")

all_cols <- colnames(community_assign_mat_fin)

values <- all_cols[grep("ASAL|PSAL", all_cols)]

vars<- values
run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal_state + (1 |ID)")
  )
  
  model <- lmer(formula, data = community_assign_mat_fin)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousal_statedrowsy") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"HCP-7T_community_aliegence_across_one_min.csv")



