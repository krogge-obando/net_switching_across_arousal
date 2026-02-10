#This code will compute a U-mann whitney test to identify if switching rates are altered across arousal state

#Computing the stats for HCP-7T data

df_net<-read.csv("/Net_Flex_MSTR.csv")

mw_results<-data.frame()

networks <- c("a_sal","p_sal","d_dmn","v_dmn","l_cen","r_cen")


for (network in networks) {
  
  # Extract values
  values <- df_net[[network]]
  group  <- df_net$Arousal_State
  
  # Rank across ALL subjects for this network
  ranks <- rank(values, na.last = "keep")
  
  # Split ranks by arousal state
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  W <- as.numeric(test$statistic)
  
  # Rank across all subjects
  
  # Rank sum for drowsy (group 1)
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # U1
  U1 <- R1 - (n1 * (n1 + 1)) / 2
  
  # Mean and SD of U
  mu_U <- (n1 * n2) / 2
  sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
  
  # z-score
  z <- (U1 - mu_U) / sd_U
  
  # Effect size r
  r <- z / sqrt(n1 + n2)
  
  # Store results
  mw_results <- rbind(
    mw_results,
    data.frame(
      Network = network,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value )
  )
}

mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

HCP_mw_results<-mw_results

write.csv(HCP_mw_results,"HCP_7T_Network_Wilcoxn_Results.csv")


##Now do the parcels

df_parcel<-read.csv("/Net_Parcel_MSTR.csv")
all_cols <- colnames(df_parcel)

# Select only parcels with dmn, sal, cen
parcels <- all_cols[grep("dmn|sal|cen", all_cols)]
parcels <- parcels[2:51]


mw_parcel_results <- data.frame()

for (parcel in parcels) {
  
  # Extract values
  values <- df_parcel[[parcel]]
  group  <- df_parcel$Arousal_State
  
  # Rank across ALL subjects for this network
  ranks <- rank(values, na.last = "keep")
  
  # Split ranks by arousal state
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  W <- as.numeric(test$statistic)
  
  # Rank across all subjects
  
  # Rank sum for drowsy (group 1)
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # U1
  U1 <- R1 - (n1 * (n1 + 1)) / 2
  
  # Mean and SD of U
  mu_U <- (n1 * n2) / 2
  sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
  
  # z-score
  z <- (U1 - mu_U) / sd_U
  
  # Effect size r
  r <- z / sqrt(n1 + n2)
  
  # Store results
  mw_parcel_results <- rbind(
    mw_parcel_results,
    data.frame(
      Parcel = parcel,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value)
  )
}

# FDR correction
mw_parcel_results$q_value <- p.adjust(mw_parcel_results$p_value, method = "BH")

HCP_parcel_results<-mw_parcel_results

write.csv(HCP_parcel_results,"/HCP_7T_Parcel_Wilcoxn_Results.csv")

##Repeat for VU-EEG-fMRI data

df_net <- read.csv("/FIND_atlas_net_sr_72TR_5_27_25.csv")

mw_results<-data.frame()

networks <- c("a_sal","p_sal","d_dmn","v_dmn","l_cen","r_cen")


for (network in networks) {
  
  # Extract values
  values <- df_net[[network]]
  group  <- df_net$arousal_state
  
  # Rank across ALL subjects for this network
  ranks <- rank(values, na.last = "keep")
  
  # Split ranks by arousal state
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  W <- as.numeric(test$statistic)
  
  # Rank across all subjects
  
  # Rank sum for drowsy (group 1)
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # U1
  U1 <- R1 - (n1 * (n1 + 1)) / 2
  
  # Mean and SD of U
  mu_U <- (n1 * n2) / 2
  sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
  
  # z-score
  z <- (U1 - mu_U) / sd_U
  
  # Effect size r
  r <- z / sqrt(n1 + n2)
  
  # Store results
  mw_results <- rbind(
    mw_results,
    data.frame(
      Network = network,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value )
  )
}

mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

VU_mw_results<-mw_results 

write.csv(VU_mw_results,"/VU_EEG_fMRI_Wilcoxn_Results.csv")

df_parcel<-read.csv("/FIND_atlas_parcel_sr_72TR_5_27_25.csv")

all_cols <- colnames(df_parcel)

# Select only parcels with dmn, sal, cen
parcels <- all_cols[grep("dmn|sal|cen", all_cols)]
parcels <- parcels[2:51]


mw_parcel_results <- data.frame()

for (parcel in parcels) {
  
  # Extract values
  values <- df_parcel[[parcel]]
  group  <- df_parcel$Arousal_State
  
  # Rank across ALL subjects for this network
  ranks <- rank(values, na.last = "keep")
  
  # Split ranks by arousal state
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  n1 <- length(group1)
  n2 <- length(group2)
  
  W <- as.numeric(test$statistic)
  
  # Rank across all subjects
  
  # Rank sum for drowsy (group 1)
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # U1
  U1 <- R1 - (n1 * (n1 + 1)) / 2
  
  # Mean and SD of U
  mu_U <- (n1 * n2) / 2
  sd_U <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
  
  # z-score
  z <- (U1 - mu_U) / sd_U
  
  # Effect size r
  r <- z / sqrt(n1 + n2)
  
  # Store results
  mw_parcel_results <- rbind(
    mw_parcel_results,
    data.frame(
      Parcel = parcel,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value)
  )
}

mw_parcel_results$q_value <- p.adjust(mw_parcel_results$p_value, method = "BH")

VU_parcel_results<-mw_parcel_results 

write.csv(VU_parcel_results,"/VU_EEG_fMRI_parcels_Wilcoxn_Results.csv")
