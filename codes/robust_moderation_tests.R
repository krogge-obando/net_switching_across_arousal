#Extra code to meet reviewer comments

#Extra Moderation analysis to increase robustness of results 

#With task sr data #RELATIONAL TASK

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)

df_net<-read.csv("HCP_3T_REL_SR_ALL_SCANS_1min.csv")

#df_net<-read.csv("HCP_3T_REL_SR_ALL_SCANS.csv")

df_net$global_sr <- rowMeans(df_net[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)

df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")

df_behavior$ID<-as.numeric(df_behavior$ID)

df_hr<-read.csv("relational_task_hr_measures.csv")

df_hr$hr_log<-log(df_hr$mean_hr)

df_mstr<-left_join(df_net,df_hr,by=c("ID","Scan"))

df_mean <- df_mstr %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x[.x != 0], na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_rel<-left_join(df_mean,df_behavior,by="ID")

df_rel$PSQI_Score

r1 <- lm(Relational_Task_Acc ~ global_sr * PSQI_Score+ Gender + Age, data = df_rel)
summary(r1)

r2 <- lm(Relational_Task_Acc ~ global_sr * mean_hr+ Gender + Age + PSQI_Score, data = df_rel)
summary(r2)


#With task sr data #Working Memory TASK
df_net<-read.csv("HCP_3T_WM_SR_ALL_SCANS_1min.csv")

df_net$global_sr <- rowMeans(df_net[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)

df_behavior<-read.csv("HCP-YA_allSubjects_behavorialData.csv")

df_behavior$ID<-as.numeric(df_behavior$ID)

df_hr<-read.csv("working_memory_task_hr_measures_fixed.csv")

df_hr$hr_log<-log(df_hr$mean_hr)

df_mstr<-left_join(df_net,df_hr,by=c("ID","Scan"))

df_mean <- df_mstr %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x[.x != 0], na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_wm<-left_join(df_mean,df_behavior,by="ID")

df_wm$PSQI_Score

wm1 <- lm(WM_Task_Acc ~ global_sr * PSQI_Score + Gender + Age , data = df_wm)
summary(wm1)

wm2 <- lm(WM_Task_Acc ~ global_sr * mean_hr+ Gender + Age + PSQI_Score, data = df_wm)
summary(wm2)

#=======================================
##Test if this moderation works for
#static correlation
#=========================================


rel_hcp_df<- read.csv("REL_static_corr.csv")
rel_hcp_df<- rel_hcp_df %>% rename(ID=Var1)
df_mean <- rel_hcp_df %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x[.x != 0], na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_rel_1<-left_join(df_mean,df_rel,by="ID")

wm_hcp_df<- read.csv("WM_static_corr.csv")

wm_hcp_df<- wm_hcp_df %>% rename(ID=Var1)

df_mean <- wm_hcp_df %>%
  group_by(ID) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x[.x != 0], na.rm = TRUE)
    ),
    .groups = "drop"
  )

df_wm_1<-left_join(df_mean,df_wm,by="ID")

# Put all models into a named list
models <- list(
   Rel_mean  = summary(lm(Relational_Task_Acc ~ stat_corr_list * mean_hr +Gender + Age, data=df_rel_1)),
   WM_mean   = summary(lm(WM_Task_Acc ~ stat_corr_list *mean_hr +Gender + Age, data=df_wm_1))
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

#write.csv(interaction_results,"/Users/roggeokk/Desktop/Projects/Net-flex-proj-VU-HCP7T/outputs/static_gs_net_task_moderation_test_2_9_26.csv")


#==============================
# MAKE Manuscript PLOTS Figure 5
#==============================

library(interactions)

SF1<-interact_plot(
  r1,
  pred = global_sr,
  modx = PSQI_Score,
  plot.points = TRUE,
  interval = FALSE,
  colors = c("blue","grey","green4"),
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

SF1


SF2<-interact_plot(
  wm1,
  pred = global_sr,
  modx = PSQI_Score,
  plot.points = TRUE,
  interval = FALSE,
  colors = c("blue","grey","green4"),
  x.label = "Global Switching",
  y.label = "Working Memory Task Accuracy",
  legend.main = "PSQI Score"
)+theme_bw()+theme(
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 20),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20)
)

SF2

SF2<-SF2 +
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

F1<-p +
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

F1


F2<-interact_plot(
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


((F1 | F2)) +
  plot_annotation(tag_levels = "a")

((SF1 | SF2)) +
  plot_annotation(tag_levels = "a")



