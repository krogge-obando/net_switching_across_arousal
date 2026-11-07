#This code will make the table that computes sd of eyeclosures for alert and drowsy scans

df_eye_alert<-read.csv("hcp_rs_eye_measures_alert.csv")
df_eye_alert$Arousal_State<-rep("High-Arousal",286)

df_eye_drowsy<-read.csv("hcp_rs_eye_measures_drowsy.csv")
df_eye_drowsy$Arousal_State<-rep("Low-Arousal",60)

df_eye<-rbind(df_eye_alert,df_eye_drowsy)

library(dplyr)

eyeclosure_summary <- df_eye %>%
  group_by(Arousal_State) %>%
  summarise(
    WB1_mean = mean(eyeclosure_wb_1, na.rm = TRUE),
    WB1_sd   = sd(eyeclosure_wb_1, na.rm = TRUE),
    
    WB2_mean = mean(eyeclosure_wb_2, na.rm = TRUE),
    WB2_sd   = sd(eyeclosure_wb_2, na.rm = TRUE),
    
    WB3_mean = mean(eyeclosure_wb_3, na.rm = TRUE),
    WB3_sd   = sd(eyeclosure_wb_3, na.rm = TRUE),
    
    WB4_mean = mean(eyeclosure_wb_4, na.rm = TRUE),
    WB4_sd   = sd(eyeclosure_wb_4, na.rm = TRUE)
  )

eyeclosure_summary

eyeclosure_summary_long <- eyeclosure_summary %>%
  tidyr::pivot_longer(
    cols = -Arousal_State,
    names_to = c("window", ".value"),
    names_pattern = "WB(\\d+)_(mean|sd)"
  ) %>%
  mutate(
    window = paste0("WB", window)
  )

eyeclosure_summary_long

write.csv(eyeclosure_summary_long,"eyeclosure_summary_across_windows.csv")

df<-read.csv("vigall_dynamic_sr_2.5min_4windows_a_d.csv")

df_a_list<-read.csv("FIND_atlas_net_sr_72TR_5_27_25 3.csv")

df_a<-select(df_a_list, c("ID","arousal_state"))

df_repeat <- df_a[rep(seq_len(nrow(df_a)), each = 4), ]

df<-cbind(df,df_repeat)

df$window<-rep(c("1","2","3","4"),22)

df <- df %>%
  rename(Arousal_State = arousal_state) %>%
  mutate(
    Arousal_State = recode(
      Arousal_State,
      "alert" = "High-Arousal",
      "drowsy" = "Low-Arousal"
    )
  )

summary_df <- df %>%
  group_by(Arousal_State, window) %>%
  summarise(
    mean = mean(vigall, na.rm = TRUE),
    sd = sd(vigall, na.rm = TRUE),
    n = sum(!is.na(vigall)),
    .groups = "drop"
  )

summary_df
write.csv(summary_df,"vigall_summary_across_windows.csv")

df_unique<-as.data.frame(df_unique)

df_unique <-df_unique %>% rename(ID=df_unique)

data<-left_join(df_unique,df_behavior,by="ID")



