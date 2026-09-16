#Extra code to meet reviewer comments

#Extra Moderation analysis to increase robustness of results 

#With task sr data #RELATIONAL TASK

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)
df_net<-read.csv("HCP_3T_REL_SR_ALL_SCANS_1min.csv")

df_net$global_sr <- rowMeans(df_net[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)

df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")

df_behavior$ID<-as.numeric(df_behavior$ID)

df_hr<-read.csv("relational_task_hr_measures.csv")

df_hr$hr_log<-log(df_hr$mean_hr)

df_mstr<-left_join(df_net,df_hr,by=c("ID","Scan"))

df_mean<-df_mstr %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

df_rel<-left_join(df_mean,df_behavior,by="ID")

df$PSQI_Score

r1 <- lm(Relational_Task_Acc ~ global_sr * PSQI_Score+ Gender + Age, data = df_rel)
summary(r1)

r2 <- lm(Relational_Task_Acc ~ global_sr * mean_hr+ Gender + Age + PSQI_Score, data = df_rel)
summary(r2)


#With task sr data #Working Memory TASK
df_net<-read.csv("HCP_3T_WM_SR_ALL_SCANS_1min.csv")

#df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")
df_hr<-read.csv("fixed_working_memory_task_hr_measures.csv")

df_hr$hr_log<-log(df_hr$mean_hr)

df_mstr<-left_join(df_net,df_hr,by=c("ID","Scan"))

df_mean<-df_mstr %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

df_wm<-left_join(df_mean,df_behavior,by="ID")

df_wm$PSQI_Score

wm1 <- lm(WM_Task_Acc ~ global_sr * PSQI_Score + Gender + Age , data = df_wm)
summary(wm1)

wm2 <- lm(WM_Task_Acc ~ global_sr * mean_hr+ Gender + Age + PSQI_Score, data = df_wm)
summary(wm2)
#remove the subjects that have both alert and drowsy

subjects_to_remove <- df_net %>%
  group_by(ID) %>%
  filter(all(c("alert", "drowsy") %in% Arousal_State)) %>%
  distinct(ID)

# Remove only rows where these subjects have "alert" in the state column that have both alert and drowsy
df_filtered <- df_net %>%
  filter(!(ID %in% subjects_to_remove$ID & Arousal_State == "alert"))

df_filtered$ID<-as.factor(df_filtered$ID)

global_network_flexibility<-df_filtered %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

global_network_flexibility$global_flex<-rowMeans(global_network_flexibility[,-1,2],na.rm=TRUE)

df_behavior$ID<-as.factor(df_behavior$ID)

df_unique <- df_filtered %>%
  distinct(ID, .keep_all = TRUE)

gnf_behav<-left_join(global_network_flexibility,df_behavior,by="ID")

gnf_behav$arousal_state<-df_unique$Arousal_State

gnf_behav$arousal_state <- factor(
  gnf_behav$arousal_state,
  levels = c("drowsy", "alert"),
  labels = c("Low arousal", "High arousal")
)

# Fit models (keep lm objects)
m1 <- lm(Relational_Task_Acc ~ global_flex * arousal_state+Gender+Age+PSQI_Score, data = gnf_behav)
m2 <- lm(WM_Task_Acc ~ global_flex * arousal_state+Gender+Age+PSQI_Score, data = gnf_behav)

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

lm_results$q_value <- p.adjust(lm_results$p_value, method = "fdr")

observed_beta<-lm_results$Beta[1] #Will use for next test

#write.csv(lm_results,"global_sr_task_moderation_2_9_26.csv")

#========================================================
#Compute analysis while iteration by net sr value used
#========================================================
df_net$ID<-as.factor(df_net$ID)
df_net$global_sr <- rowMeans(df_net[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)



gnf_behav2<-left_join(df_net,df_behavior,by="ID")

gnf_behav2<-gnf_behav2 %>%
  filter(
    !is.na(Relational_Task_Acc),
    !is.na(global_sr),
    !is.na(Arousal_State)
  )


set.seed(123)

n_iter <- 1000

interaction_results <- map_dfr(1:n_iter, function(i) {
  
  sampled_data <- gnf_behav2 %>%
    group_by(ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  model <- lm(
    Relational_Task_Acc ~ global_sr * Arousal_State,
    data = sampled_data
  )
  
  coef_table <- summary(model)$coefficients
  
  interaction_row <- grep(
    "global_sr:Arousal_State",
    rownames(coef_table)
  )
  
  tibble(
    iteration = i,
    beta = coef_table[interaction_row, "Estimate"],
    p_value = coef_table[interaction_row, "Pr(>|t|)"]
  )
})

interaction_results %>%
  summarise(
    mean_beta = mean(beta),
    median_beta = median(beta),
  )

#======================================
#Compute permutated p-value from betas
#=======================================
mean(abs(interaction_results$beta) >= abs(observed_beta))

#==================================================
#Making supplementary figure over statistic results
#==================================================
library(ggplot2)

ggplot(interaction_results, aes(x = beta)) +
  geom_histogram(
    bins = 50
  ) +
  geom_vline(
    xintercept = observed_beta,
    linewidth = 1,
    color="blue"
  ) +
  labs(
    x = "Interaction coefficient (β)",
    y = "Count",
    title = "Permutation Null Distribution of Moderation Effect",
    subtitle = paste0("Observed β = ", round(observed_beta, 3))
  ) +
  theme_classic(base_size = 16)
#==============================================
#Rerun analysis while having PSQ as a covariate
#===============================================

gnf_behav$PSQI_Score

# Fit models (keep lm objects)
m1 <- lm(Relational_Task_Acc ~ global_flex * arousal_state + Age +Gender +PSQI_Score, data = gnf_behav)
m2 <- lm(WM_Task_Acc ~ global_flex * arousal_state + Age + Gender+ PSQI_Score, data = gnf_behav)

get_interaction <- function(model) {
  ct <- coef(summary(model))
  term <- grep("^global_sr:arousal_state", rownames(ct), value = TRUE)
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

lm_results$q_value <- p.adjust(lm_results$p_value, method = "fdr")
#write.csv(lm_results,"global_sr_task_moderation_w_PSQI_Score_8_6_26.csv")

#=======================================
##Test if this moderation works for any 
#of the global signal network correlation 
#or static correlation
#=========================================

hcp_df <- read.csv("HCP_7T_measures_of_interest.csv")

#df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")

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

#write.csv(interaction_results,"static_gs_net_task_moderation_test_2_9_26.csv")


#==============================
# MAKE Manuscript PLOTS Figure 5
#==============================

library(interactions)

F1<-interact_plot(
  r1,
  pred = global_sr,
  modx = PSQI_Score,
  plot.points = TRUE,
  interval = FALSE,
  colors = c("purple","grey","hotpink"),
  x.label = "Global Switching",
  y.label = "Relational Task Accuracy",
  legend.main = "PSQI Score"
)+theme_bw()+theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

F1


p<-interact_plot(
  r2,
  pred = global_sr,
  modx = mean_hr,
  plot.points = TRUE,
  interval = FALSE,
  modx.values = c(40,60,80),
  colors = c("purple","grey","hotpink"),
  x.label = "Global Switching",
  y.label = "Relational Task Accuracy",
  legend.main = "Mean HR"
)+theme_bw()+theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

F2<-p +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "*",
    color = "green2",
    size = 10,
    hjust = 1.5,
    vjust = 1.5
  )

F2

f4<-interact_plot(
  wm1,
  pred = global_sr,
  modx = PSQI_Score,
  plot.points = TRUE,
  interval = FALSE,
  colors = c("purple","grey","hotpink"),
  x.label = "Global Switching",
  y.label = "Relational Task Accuracy",
  legend.main = "PSQI Score"
)+  scale_x_continuous(
  breaks = c(0.005, 0.010, 0.015, 0.020),
  limits = c(0.005, 0.020)) +
  theme_bw()+theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

F4<-f4 +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "*",
    color = "green2",
    size = 10,
    hjust = 1.5,
    vjust = 1.5
  )

F4

F5<-interact_plot(
  wm2,
  pred = global_sr,
  modx = mean_hr,
  plot.points = TRUE,
  interval = FALSE,
  modx.values = c(40,60,80),
  colors = c("purple","grey","hotpink"),
  x.label = "Global Switching",
  y.label = "Working Memory Task Accuracy",
  legend.main = "Mean HR"
) +  scale_x_continuous(
  breaks = c(0.005, 0.010, 0.015, 0.020),
  limits = c(0.005, 0.020))+
  theme_bw()+theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

F5


#===============================
# HCP-7T with Arousal State
#===============================


df_net<-read.csv("Net_Flex_MSTR.csv")

df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")

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

gnf_behav$arousal_state <- factor(
  gnf_behav$arousal_state,
  levels = c("drowsy", "alert"),
  labels = c("Low arousal", "High arousal")
)

library(ggplot2)

F3 <- ggplot(gnf_behav, aes(x = global_flex, y = Relational_Task_Acc, color = arousal_state)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se=FALSE) +
  scale_color_manual(
    values = c("High arousal" = "hotpink", "Low arousal" = "purple")
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) +
  labs(
    x = "Global Brain Switching",
    y = "Relational Task Accuracy",
    color = "Arousal State"
  ) +
  annotate(
    "text",
    x = 0.011,
    y = 97,
    label = "*",
    size = 10,
    color = "green2"
  )


F6 <- ggplot(gnf_behav, aes(x = global_flex, y = WM_Task_Acc, color = arousal_state)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",se=FALSE) +
  scale_color_manual(
    values = c("High arousal" = "hotpink", "Low arousal" = "purple")
  ) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) +
  labs(
    x = "Global Brain Switching",
    y = "Working Memory Task Accuracy",
    color = "Arousal State"
  )

#((F1 | F2 | F3)/( F4 | F5 | F6)) + plot_annotation(tag_levels = "a")

((F1 | F4) / (F2 | F5) / (F3 | F6)) +
  plot_annotation(tag_levels = "a") &
  theme(axis.title.y = element_blank())



