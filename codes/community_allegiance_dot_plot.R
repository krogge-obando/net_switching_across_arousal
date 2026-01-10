#This code will compute t-test on the community assignment pairs, correct for bonferonii and lastly
# we will make a dot plot comparing the alert vs drowsy values for alert vs drowsy community assignment 
#and put a star on top of the ones that are significant

library(reshape2)
library(ggplot2)

compute_co_assignment_matrix <- function(df_com) {
  num_time_points <- ncol(df_com)
  n <- nrow(df_com)
  
  net_list<-c("asal","psal","ddmn","vdmn","lcen","rcen","aud","basg",
              "hvis","lang","prec","pvis","sens","viss")
  
  # Initialize the co-assignment matrix with zeros
  co_assignment_matrix <- matrix(0, nrow = n, ncol = n,dimnames = list(net_list,net_list))
  
  for (t in 1:num_time_points) {
    community_assignments <- df_com[, t]
    
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (community_assignments[i] == community_assignments[j]) {
          co_assignment_matrix[i, j] <- co_assignment_matrix[i, j] + 1
          co_assignment_matrix[j, i] <- co_assignment_matrix[j, i] + 1
        }
      }
    }
  }
  co_assignment_matrix <- co_assignment_matrix / num_time_points
  return(co_assignment_matrix)
}


#compute t-test investigations
alert_scans<-c("vcon02_scan01","vcon07_scan01","vcon15_scan01","vcon15_scan02","vcon17_scan02","vcon19_scan02","vcon20_scan01", "vcon20_scan02","vcon29_scan01","vcon29_scan02")

full_path <- "[enter full path]"

df_com_first<-read.csv("vcon02_scan01_community.csv",0)


df_paths<-as.data.frame(alert_scans)

colnames(df_paths) <-"file_name"

# Define the network pairs you want to extract
network_pairs <- list(
  c("asal", "ddmn"),
  c("asal", "vdmn"),
  c("asal", "lcen"),
  c("asal", "rcen"),
  c("psal", "ddmn"),
  c("psal", "vdmn"),
  c("psal", "lcen"),
  c("psal", "rcen")
)

# Define labels for each pair (used as column names)
pair_labels <- sapply(network_pairs, function(p) paste(toupper(p[1]), toupper(p[2]), sep = "–"))

# Initialize results data frame
results_wide <- data.frame()



# Loop through each subject
for (i in 1:nrow(df_paths)) {
  file_to_load <- paste0(full_path, df_paths$file_name[i], '_community.csv')
  
  # Load community assignment matrix
  df_com <- read.csv(file_to_load, 0)  # use row.names = 1 if there's a row header
  
  # Compute co-assignment matrix
  sub_matrix <- compute_co_assignment_matrix(df_com)
  
  # Extract values for specified pairs
  values <- sapply(network_pairs, function(p) sub_matrix[p[1], p[2]])
  
  # Create a one-row data frame for this subject
  subject_row <- data.frame(Subject = df_paths$file_name[i], t(values), stringsAsFactors = FALSE)
  
  # Bind to main results
  results_wide <- rbind(results_wide, subject_row)
}

# Set column names (first is subject, rest are network pairs)
colnames(results_wide)[-1] <- pair_labels

results_wide$arousal_state <- "alert"

community_assign_mat_a<-results_wide


#repeat for drowsy scans 
n_nodes<-14

drowsy_scans<-c( "vcon05_scan01", "vcon09_scan01", "vcon09_scan02", "vcon12_scan01", "vcon14_scan01",
                 "vcon22_scan01", "vcon23_scan01", "vcon24_scan01", "vcon25_scan02", "vcon32_scan01", "vcon36_scan01", "vcon37_scan01")

df_paths<-as.data.frame(drowsy_scans)
colnames(df_paths)<-"file_names"


results_wide <- data.frame()

# Loop through each subject
for (i in 1:nrow(df_paths)) {
  file_to_load <- paste0(full_path, df_paths$file_name[i], '_community.csv')
  
  # Load community assignment matrix
  df_com <- read.csv(file_to_load, 0)  # use row.names = 1 if there's a row header
  
  # Compute co-assignment matrix
  sub_matrix <- compute_co_assignment_matrix(df_com)
  
  # Extract values for specified pairs
  values <- sapply(network_pairs, function(p) sub_matrix[p[1], p[2]])
  
  # Create a one-row data frame for this subject
  subject_row <- data.frame(Subject = df_paths$file_name[i], t(values), stringsAsFactors = FALSE)
  
  # Bind to main results
  results_wide <- rbind(results_wide, subject_row)
}

# Set column names (first is subject, rest are network pairs)
colnames(results_wide)[-1] <- pair_labels

results_wide$arousal_state <- "drowsy"
community_assign_mat_d<-results_wide

community_assign_mat_fin<-rbind(community_assign_mat_a,community_assign_mat_d)

#Initialize a data frame to store results

# Get only the network-pair columns
network_columns <- setdiff(colnames(results_wide), c("Subject", "arousal_state"))

mw_results <- data.frame(
  Network_Pair = character(),
  p_value = numeric(),
  statistic = numeric(),
  stringsAsFactors = FALSE
)

# Loop over each network-pair column
for (col in network_columns) {
  # Get values for each arousal state
  group1 <- community_assign_mat_fin[community_assign_mat_fin$arousal_state == "drowsy", col]
  group2 <- community_assign_mat_fin[community_assign_mat_fin$arousal_state == "alert", col]
  
  # Run Wilcoxon rank-sum test (Mann–Whitney U)
  test <- wilcox.test(group1, group2, exact = FALSE, conf.int= TRUE)
  
  # Store results
  mw_results <- rbind(mw_results, data.frame(
    Network_Pair = col,
    p_value = test$p.value,
    statistic = test$statistic
  ))
}

mw_results$p_adj_FDR <- p.adjust(mw_results$p_value, method = "BH")

write.csv(mw_results,"community_assignment_mw_results_bfcorr.csv",row.names=FALSE)

##MAKE FIGURE 4 A)




df<-read.csv("coassignment_network_pairs_wide_VU.csv")

df_melt <- melt(df,id.vars=c("arousal_state","Subject"), variable.name="Community_Assignment_Pair",value.name = "Fraction")

# Wilcoxon tests for each community pair
library(dplyr)

results <- df_melt %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(
    p_value = wilcox.test(Fraction ~ arousal_state)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         significance = ifelse(p_adj < 0.05, "*", ""))  # You can add "**" for lower thresholds if you like


df_plot <- df_melt %>%
  left_join(results, by = "Community_Assignment_Pair")

# Compute Y positions for stars (just above max of each group)
label_positions <- df_plot %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(y_pos = max(Fraction) + 0.07)  # Adjust the offset if needed

# Combine with significance results
star_labels <- left_join(results, label_positions, by = "Community_Assignment_Pair")

ggplot(df_plot, aes(x = Community_Assignment_Pair, y = Fraction, color = arousal_state)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  stat_summary(
    fun = mean,
    geom = "point",
    aes(group = arousal_state),
    position = position_dodge(width = 0.5),
    color = "yellow3",
    size = 4,
    shape = 18
  ) +
  geom_text(data = star_labels,
            aes(x = Community_Assignment_Pair, y = y_pos, label = significance),
            inherit.aes = FALSE,
            size = 9, color = "black") +
  scale_color_manual(values = c("drowsy" = "purple", "alert" = "hotpink")) +
  theme_classic() +
  theme(
    legend.title= element_text(size=14),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  labs(x = "Community Assignment Pair", y = "Fraction of Community Allegiance")+ylim(0,1.09)



