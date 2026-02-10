#This code will see if cog-flex relates to network-switch

#Load df

library(tidyverse)
library(dplyr)
library(stats)

df_net<-read.csv("/Users/roggeokk/Desktop/Projects/HCP_7T_SwitchingRate/data/Net_Flex_MSTR.csv")

df_behavior<-read.csv("/Users/roggeokk/Desktop/Projects/HCP_7T_SwitchingRate/data/HCP-YA_allSubjects_behavorialData.csv")

#remove the subjects that have both alert and drowsy

subjects_to_remove <- df_net %>%
  group_by(ID) %>%
  filter(all(c("alert", "drowsy") %in% Arousal_State)) %>%
  distinct(ID)

# Remove only rows where these subjects have "alert" in the state column that have both alert and drowsy
df_filtered <- df_net %>%
  filter(!(ID %in% subjects_to_remove$ID & Arousal_State == "alert"))

df_filtered$ID<-as.factor(df_filtered$ID)

##See if our network flexibility relates to hours of sleep the night before

global_network_flexibility<-df_filtered %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

global_network_flexibility$global_flex<-rowMeans(global_network_flexibility[,-1,2],na.rm=TRUE)

df_behavior$ID<-as.factor(df_behavior$ID)

df_unique <- df_filtered %>%
  distinct(ID, .keep_all = TRUE)

gnf_behav<-left_join(global_network_flexibility,df_behavior,by="ID")

gnf_behav$arousal_state<-df_unique$Arousal_State


# Fit models (keep lm objects)
m1 <- lm(Relational_Task_Acc ~ global_flex * arousal_state, data = gnf_behav)
m2 <- lm(WM_Task_Acc ~ global_flex * arousal_state, data = gnf_behav)

get_interaction <- function(model) {
  ct <- coef(summary(model))
  term <- grep("^global_flex:arousal_state", rownames(ct), value = TRUE)
  ct[term, c("Estimate", "Std. Error", "Pr(>|t|)")]
}

lm_results <- data.frame(
  Task = c("Relational", "WorkingMemory"),
  Beta = c(
    get_interaction(m1)[1],
    get_interaction(m2)[1]
  ),
  SE = c(
    get_interaction(m1)[2],
    get_interaction(m2)[2]
  ),
  p_value = c(
    get_interaction(m1)[3],
    get_interaction(m2)[3]
  )
)

lm_results$q_value <- p.adjust(lm_results$p_value, method = "BH")

write.csv(lm_results,"/Users/roggeokk/Desktop/Projects/Net-flex-proj-VU-HCP7T/outputs/global_sr_task_moderation_2_9_26.csv")


##Test if this moderation works for any of the global signal network correlation or static correlation

hcp_df <- read.csv("/Users/roggeokk/Desktop/Projects/Net-flex-proj-VU-HCP7T/data/HCP_7T_measures_of_interest.csv")

df_behavior<-read.csv("/Users/roggeokk/Desktop/Projects/HCP_7T_SwitchingRate/data/HCP-YA_allSubjects_behavorialData.csv")

#remove the subjects that have both alert and drowsy

hcp_clean <- hcp_df %>%
  group_by(ID) %>%
  filter(!all(c("alert", "drowsy") %in% arousal)) %>%
  ungroup()

hcp_means <- hcp_clean %>%
  group_by(ID, arousal) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

hcp_mean_behav_df<-left_join(hcp_means,df_behavior, by="ID")

# Put all models into a named list
models <- list(
  Rel_asal  = summary(lm(Relational_Task_Acc ~ gs_asal + arousal + gs_asal:arousal, data=hcp_mean_behav_df)),
  Rel_ddmn  = summary(lm(Relational_Task_Acc ~ gs_ddmn + arousal + gs_ddmn:arousal, data=hcp_mean_behav_df)),
  Rel_vdmn  = summary(lm(Relational_Task_Acc ~ gs_vdmn + arousal + gs_vdmn:arousal, data=hcp_mean_behav_df)),
  Rel_mean  = summary(lm(Relational_Task_Acc ~ mean_stat_corr + arousal + mean_stat_corr:arousal, data=hcp_mean_behav_df)),
  
  WM_asal   = summary(lm(WM_Task_Acc ~ gs_asal + arousal + gs_asal:arousal, data=hcp_mean_behav_df)),
  WM_ddmn   = summary(lm(WM_Task_Acc ~ gs_ddmn + arousal + gs_ddmn:arousal, data=hcp_mean_behav_df)),
  WM_vdmn   = summary(lm(WM_Task_Acc ~ gs_vdmn + arousal + gs_vdmn:arousal, data=hcp_mean_behav_df)),
  WM_mean   = summary(lm(WM_Task_Acc ~ mean_stat_corr + arousal + mean_stat_corr:arousal, data=hcp_mean_behav_df))
)

# Initialize output
interaction_results <- data.frame()

for (model_name in names(models)) {
  
  mod <- models[[model_name]]
  coef_tbl <- as.data.frame(mod$coefficients)
  
  coef_tbl$Model <- model_name
  coef_tbl$Term <- rownames(coef_tbl)
  
  # Keep ONLY interaction terms
  interaction_only <- coef_tbl[grepl(":", coef_tbl$Term), ]
  
  interaction_results <- rbind(interaction_results, interaction_only)
}

# Rename columns
colnames(interaction_results)[1:4] <- c("Beta", "SE", "t_value", "p_value")

# FDR correction on interaction p-values only
interaction_results$FDR_p <- p.adjust(interaction_results$p_value, method = "fdr")

# Final ordering
interaction_results <- interaction_results[, c("Model", "Term", "Beta", "SE", "t_value", "p_value", "FDR_p")]

write.csv(interaction_results,"/Users/roggeokk/Desktop/Projects/Net-flex-proj-VU-HCP7T/outputs/static_gs_net_task_moderation_test_2_9_26.csv")
