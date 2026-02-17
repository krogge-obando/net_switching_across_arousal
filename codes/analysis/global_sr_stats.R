#This code will compute the statistics for global switching across arousal state

#Conduct for HCP-7T

df_net<-read.csv("Net_Flex_MSTR.csv")


df_net$global_flex <- apply(df_net[,4:17], 1, mean, na.rm = TRUE)

df_net$ranks <- rank(df_net$global_flex)

df_drowsy<-subset(df_net,df_net$Arousal_State=="drowsy")

df_alert<-subset(df_net, df_net$Arousal_State=="alert")

rank_mean_drowsy<-mean(df_drowsy$ranks)

rank_mean_alert<-mean(df_alert$ranks)

test<-wilcox.test(df_drowsy$global_flex,df_alert$global_flex)

test$p.value

n1 <- nrow(df_drowsy)
n2 <- nrow(df_alert)

R1 <- sum(df_drowsy$ranks, na.rm = TRUE)

# U statistic
U1 <- R1 - (n1 * (n1 + 1)) / 2

# Mean and SD of U
mu_U <- (n1 * n2) / 2
sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)

# z-score
z <- (U1 - mu_U) / sd_U

# Effect size r
r <- z / sqrt(n1 + n2)


HCP_df_wilcox_summary <- data.frame(
  U=U1,
  z = z,
  r = r,
  p_value = test$p.value,
  rank_mean_drowsy = rank_mean_drowsy,
  rank_mean_alert = rank_mean_alert
)

write.csv(HCP_df_wilcox_summary,"HCP_global_sr_results.csv")

df_net <- read.csv("FIND_atlas_net_sr_72TR_5_27_25.csv")

df_net$arousal_state<-as.factor(df_net$arousal_state)

df_net$global_flex <- apply(df_net[,6:19], 1, mean, na.rm = TRUE)

df_drowsy<-subset(df_net,df_net$arousal_state=="drowsy")

df_alert<-subset(df_net, df_net$arousal_state=="alert")

df_net$global_flex <- apply(df_net[,6:19], 1, mean, na.rm = TRUE)

test<-wilcox.test(df_drowsy$global_flex,df_alert$global_flex)

df_net$ranks <- rank(df_net$global_flex)

df_drowsy<-subset(df_net,df_net$arousal_state=="drowsy")

df_alert<-subset(df_net, df_net$arousal_state=="alert")

rank_mean_drowsy<-mean(df_drowsy$ranks)

rank_mean_alert<-mean(df_alert$ranks)


W <- as.numeric(test$statistic)

n1 <- nrow(df_drowsy)
n2 <- nrow(df_alert)

  R1 <- sum(df_drowsy$ranks, na.rm = TRUE)
  
  # U statistic
  U1 <- R1 - (n1 * (n1 + 1)) / 2
  
  # Mean and SD of U
  mu_U <- (n1 * n2) / 2
  sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
  
  # z-score
  z <- (U1 - mu_U) / sd_U
  
  # Effect size r
  r <- z / sqrt(n1 + n2)


# Summary dataframe
vu_df_wilcox_summary <- data.frame(
  U=U1,
  z = z,
  r = r,
  p_value = test$p.value,
  rank_mean_drowsy = rank_mean_drowsy,
  rank_mean_alert = rank_mean_alert
)

write.csv(vu_df_wilcox_summary,"VU_EEG_fMRI_global_sr_results.csv")







