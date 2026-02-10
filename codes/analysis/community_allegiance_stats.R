#This code will compute the stats for comparing Community alleigience across arousal state.

#Load data
community_assign_mat_fin<-read.csv("/HCP_7T_community_allegience.csv")

# Get only the network-pair columns
network_columns <- setdiff(colnames(community_assign_mat_fin), c("Subject", "arousal_state"))

mw_results <- data.frame()

# Get only the network-pair columns
network_columns <- setdiff(
  colnames(community_assign_mat_fin),
  c("Subject", "arousal_state")
)

mw_results <- data.frame()

for (col in network_columns) {
  
  # Pull values and group
  values <- community_assign_mat_fin[[col]]
  group  <- community_assign_mat_fin$arousal_state
  
  # Rank across ALL subjects
  ranks <- rank(values, na.last = "keep")
  
  # Rank means
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for wilcox)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  # Sample sizes
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # Rank sum for drowsy
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  # U statistic
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
      NetworkPair = col,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value
    )
  )
}

# FDR correction
mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

hcp_mw_results<-mw_results

write.csv(hcp_mw_results,"/HCP_7T_community_allegiance_results.csv",row.names=FALSE)


#VU-EEG-fMRI


community_assign_mat_fin<-read.csv("eeg_fmri_community_allegience_2_9_26.csv")

# Get only the network-pair columns
network_columns <- setdiff(colnames(community_assign_mat_fin), c("Subject", "arousal_state"))

mw_results <- data.frame()

for (col in network_columns) {
  
  # Pull values and group
  values <- community_assign_mat_fin[[col]]
  group  <- community_assign_mat_fin$arousal_state
  
  # Rank across ALL subjects
  ranks <- rank(values, na.last = "keep")
  
  # Rank means
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for wilcox)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney test
  test <- wilcox.test(group1, group2, exact = FALSE, correct = FALSE)
  
  # Sample sizes
  n1 <- sum(group == "drowsy")
  n2 <- sum(group == "alert")
  
  # Rank sum for drowsy
  R1 <- sum(ranks[group == "drowsy"], na.rm = TRUE)
  
  # U statistic
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
      NetworkPair = col,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = test$p.value
    )
  )
}

# FDR correction
mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

write.csv(mw_results,"/EEG_fMRI_VU_community_allegiance.csv",row.names=FALSE)










