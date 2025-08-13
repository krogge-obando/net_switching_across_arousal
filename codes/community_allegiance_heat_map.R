#This script will generate a count of how many times a network was in the same community


library(reshape2)
library(ggplot2)

# Make function to compute co-assignment-matrix
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


alert_scans<-c("vcon02_scan01","vcon07_scan01","vcon15_scan01","vcon15_scan02","vcon17_scan02","vcon19_scan02","vcon20_scan01", "vcon20_scan02","vcon29_scan01","vcon29_scan02")

full_path <- "/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/community_assignment/"

df_com_first<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/community_assignment/vcon02_scan01_community.csv",0)
n_nodes<-nrow(df_com_first)

sum_matrix<-matrix(0 ,nrow= n_nodes, ncol = n_nodes)

df_paths<-as.data.frame(alert_scans)

colnames(df_paths) <-"file_name"

for (i in 1:nrow(df_paths)) {
  
  file_to_load <- paste0(full_path, df_paths$file_name[i],'_community.csv')
  
  df_com <- read.csv(file_to_load,0)

  #df_com<-read.csv("/Users/roggeokk/Desktop/Projects/HCP_7T_SwitchingRate/data/100610-REST1_PA_community.csv")

  sub_matrix<-compute_co_assignment_matrix(df_com)
  
  sum_matrix<-sum_matrix + sub_matrix

}

alert_average_matrix <- sum_matrix/ length(df_paths$file_name)


##Repeat for drowsy nodes

# Initialize the co-assignment matrix with zeros
n_nodes<-14

drowsy_scans<-c( "vcon05_scan01", "vcon09_scan01", "vcon09_scan02", "vcon12_scan01", "vcon14_scan01",
                 "vcon22_scan01", "vcon23_scan01", "vcon24_scan01", "vcon25_scan02", "vcon32_scan01", "vcon36_scan01", "vcon37_scan01")

df_paths<-as.data.frame(drowsy_scans)
colnames(df_paths)<-"file_names"

sum_matrix_d<-matrix(0 ,nrow= n_nodes, ncol = n_nodes)

for (i in 1:nrow(df_paths)) {
  
  file_to_load <- paste0(full_path, df_paths$file_name[i],'_community.csv')
  
  df_com <- read.csv(file_to_load,0)
  
  #df_com<-read.csv("/Users/roggeokk/Desktop/Projects/HCP_7T_SwitchingRate/data/100610-REST1_PA_community.csv")
  
  sub_matrix_d<-compute_co_assignment_matrix(df_com)
  
  sum_matrix_d<-sum_matrix_d + sub_matrix_d
  
}

drowsy_average_matrix <- sum_matrix_d/ length(df_paths$file_name)


# Define the subset of networks you want to keep
subset_networks <- c("asal", "psal", "ddmn", "vdmn", "rcen", "lcen")

network_labels_sub <- toupper(subset_networks)



# Subset your matrices by rows and columns for these networks only
alert_subset <- alert_average_matrix[subset_networks, subset_networks]
drowsy_subset <- drowsy_average_matrix[subset_networks, subset_networks]

# Create network order and labels for subset
network_order_sub <- rownames(alert_subset)
network_labels_sub <- network_order_sub
network_indices_sub <- setNames(seq_along(network_order_sub), network_order_sub)

# Melt the subset matrices
alert_melted_sub <- melt(alert_subset)
colnames(alert_melted_sub) <- c("Network1", "Network2", "Fraction")
alert_melted_sub$Percentage <- round(alert_melted_sub$Fraction * 100, 1)

drowsy_melted_sub <- melt(drowsy_subset)
colnames(drowsy_melted_sub) <- c("Network1", "Network2", "Fraction")
drowsy_melted_sub$Percentage <- round(drowsy_melted_sub$Fraction * 100, 1)

# Map networks to numeric positions for plotting
alert_melted_sub$x <- network_indices_sub[alert_melted_sub$Network1]
alert_melted_sub$y <- network_indices_sub[alert_melted_sub$Network2]

drowsy_melted_sub$x <- network_indices_sub[drowsy_melted_sub$Network1]
drowsy_melted_sub$y <- network_indices_sub[drowsy_melted_sub$Network2]

# Keep only triangles for plotting
alert_tri_sub <- alert_melted_sub[alert_melted_sub$x < alert_melted_sub$y, ]
drowsy_tri_sub <- drowsy_melted_sub[drowsy_melted_sub$x > drowsy_melted_sub$y, ]

# Now plot
ggplot() +
  geom_tile(data = drowsy_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = drowsy_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white",
                      high = "purple", 
                      name = "Drowsy",
                      limits=c(0,0.9), 
                      breaks=c(0, 0.30,0.60,0.90), 
                      labels=c("0%","30%","60%","90%")) +
  
  new_scale_fill() +
  
  geom_tile(data = alert_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = alert_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white", 
                      high = "hotpink", 
                      name = "Alert", 
                      limits=c(0,0.9), 
                      breaks=c(0,  0.30,  0.60, 0.90),
                      labels=c("0%","30%","60%","90%")) +
  
  scale_x_continuous(breaks = seq_along(network_labels_sub), labels = network_labels_sub, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(network_labels_sub), labels = network_labels_sub, expand = c(0, 0)) +
  
  coord_fixed() +
  theme_minimal() +
  theme(
    panel.border=element_rect(color = "black", fill = NA, size = 1),
    legend.text=element_text(size=12),
    legend.title=element_text(size=20),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  labs(
    x = "Networks", y = "Networks"
  )

  
  
  
  
 