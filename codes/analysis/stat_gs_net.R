#This code will derive the statistics of static correlation and network to global mean correlation across arousal state

vu_df <- read.csv("/eeg_fmri_vu_measures_of_interest.csv")

# Map network names to column names in vu_df
network_cols <- c(
  DDMN     = "gs_ddmn",
  VDMN     = "gs_vdmn",
  ASAL     = "gs_asal",
  STAT_COR = "mean_stat_corr"
)

mw_results <- data.frame()

for (network in names(network_cols)) {
  
  # Pull raw values
  values <- vu_df[[ network_cols[network] ]]
  group  <- vu_df$arousal_state
  
  # Rank across ALL subjects
  ranks <- rank(values, na.last = "keep")
  
  # Rank means
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
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
  
  # p-value from your precomputed test
  

  p_val <- test$p.value
  
  # Store
  mw_results <- rbind(
    mw_results,
    data.frame(
      Network = network,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = p_val
    )
  )
}

mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

write.csv(mw_results,"/EEG_fMRI_VU_global_static_corr.csv")


hcp_df <- read.csv("/HCP_7T_measures_of_interest.csv")

# Map network names to column names in vu_df
network_cols <- c(
  DDMN     = "gs_ddmn",
  VDMN     = "gs_vdmn",
  ASAL     = "gs_asal",
  STAT_COR = "mean_stat_corr"
)

mw_results <- data.frame()

for (network in names(network_cols)) {
  
  # Pull raw values
  values <- hcp_df[[ network_cols[network] ]]
  group  <- hcp_df$arousal
  
  # Rank across ALL subjects
  ranks <- rank(values, na.last = "keep")
  
  # Rank means
  rank_mean_drowsy <- mean(ranks[group == "drowsy"], na.rm = TRUE)
  rank_mean_alert  <- mean(ranks[group == "alert"],  na.rm = TRUE)
  
  # Values by group (for MW test)
  group1 <- values[group == "drowsy"]
  group2 <- values[group == "alert"]
  
  # Mann–Whitney U test
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
  
  # p-value from your precomputed test
  
  
  p_val <- test$p.value
  
  # Store
  mw_results <- rbind(
    mw_results,
    data.frame(
      Network = network,
      U_value = U1,
      RankMean_Drowsy = rank_mean_drowsy,
      RankMean_Alert  = rank_mean_alert,
      z_value = z,
      r_value = r,
      p_value = p_val
    )
  )
}


mw_results$q_value <- p.adjust(mw_results$p_value, method = "BH")

hcp_mw_results<-mw_results

write.csv(hcp_mw_results,"/HCP_7T_global_static_corr.csv")

