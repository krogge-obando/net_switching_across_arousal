#ANALYSIS 1 Comparing SR across Arousal State 
# Recode arousal labels throughout the analysis

library(lmerTest)
library(dplyr)
library(broom.mixed)
library(tidyverse)
library(reshape2)
library(ggplot2)

library(ggnewscale)
library(patchwork)

recode_arousal <- function(df) {
  df %>%
    mutate(across(
      any_of(c("Arousal_State", "arousal_state", "arousal")),
      ~ dplyr::recode(
        as.character(.x),
        "drowsy" = "Low Arousal",
        "Low Arousal" = "Low Arousal",
        "alert" = "High Arousal",
        "High Arousal" = "High Arousal"
      )
    ))
}

df_net<-read.csv("Net_Flex_MSTR.csv")
df_net <- recode_arousal(df_net)

# Variables to test
vars <- c("a_sal", "p_sal", "d_dmn", "v_dmn", "r_cen", "l_cen")

test<-lmer(a_sal ~ Arousal_State + (1 | ID), data = df_net)

summary(test)

# Function to fit model and extract results
run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ Arousal_State + (1 |ID)")
  )
  
  model <- lmer(formula, data = df_net)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "Arousal_StateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"HCP-7T_SR_across_arousal-9-11-26.csv")

#MAKE FIGURE 1

df_hcp_sub <- df_net %>%
  select(Arousal_State, a_sal, p_sal, d_dmn, v_dmn, l_cen, r_cen)

# Reshape to long format (modern replacement for melt)
df_hcp_melt <- df_hcp_sub %>%
  pivot_longer(
    cols = c(a_sal, p_sal, d_dmn, v_dmn, l_cen, r_cen),
    names_to = "network",
    values_to = "flexibility"
  )

# Clean network names
df_hcp_melt$network <- toupper(gsub("_", "", df_hcp_melt$network))

df_hcp_melt$network <- factor(
  df_hcp_melt$network,
  levels = c("ASAL", "PSAL", "DDMN", "VDMN", "LCEN", "RCEN")
)

# Star positions (IMPORTANT: use factor positions directly, not numeric x)
star_positions <- data.frame(
  network = c("ASAL", "DDMN", "VDMN"),
  y = 0.037,
  label = c("*", "**", "*")
)

y_limits <- c(-0.01, 0.045)
y_breaks <- seq(0, 0.04, by = 0.01)

# Plot
f1_a <- ggplot(df_hcp_melt,
               aes(x = network, y = flexibility, fill = Arousal_State)) +
  
  geom_violin(trim = FALSE, alpha = 0.5,
              position = position_dodge(0.9)) +
  
  geom_boxplot(width = 0.2,
               position = position_dodge(0.9),
               outlier.shape = NA) +
  
  scale_fill_manual(values = c("High Arousal" = "#FF1493",
                               "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  
  theme_classic(base_size = 16) +
  
  labs(title = "HCP-7T: Network-Level Switching",
       x = "Network",
       y = "Switching") +
  
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    expand = expansion(mult = c(0.25, 0.15))
  ) +
  
  # Stars (now uses factor positions correctly)
  geom_text(
    data = star_positions,
    aes(x = network, y = y, label = label),
    inherit.aes = FALSE,
    size = 10,
    color = "green3"
  ) +
  
  # Horizontal brace (FIXED: linewidth instead of size)
  geom_segment(
    data = star_positions,
    aes(x = network, xend = network,
        y = y - 0.001, yend = y - 0.001),
    inherit.aes = FALSE,
    linewidth = 1,
    color = "green3"
  ) +
  
  geom_segment(
    data = star_positions,
    aes(x = network, xend = network,
        y = y - 0.0019, yend = y - 0.0019),
    inherit.aes = FALSE,
    linewidth = 1,
    color = "black"
  )


df_parcel<-read.csv("Net_Parcel_MSTR.csv")
df_parcel <- recode_arousal(df_parcel)
all_cols <- colnames(df_parcel)

# Select only parcels with dmn, sal, cen
parcels <- all_cols[grep("dmn|sal|cen", all_cols)]
parcels <- parcels[2:51]

vars<- parcels

run_lmm_parcel <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ Arousal_State + (1 |ID)")
  )
  
  model <- lmer(formula, data = df_parcel)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "Arousal_StateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_parcel))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"HCP-7T_SR_across_arousal_parcels-9-11-26.csv")

##MAKE FIGURES

# Subset highlighted parcels (significant and non-significant)
df_hcp_sal <- df_parcel %>%
  select(Arousal_State, a_sal5, a_sal6, a_sal7, p_sal3, p_sal5, p_sal7,p_sal8, p_sal10, p_sal11, p_sal12)

# Subset highlighted parcels (significant and non-significant)
df_hcp_dmn <- df_parcel %>%
  select(Arousal_State,d_dmn1,d_dmn2,d_dmn4,d_dmn6,d_dmn7,d_dmn8, v_dmn1, v_dmn2, v_dmn5, v_dmn7)

# Subset highlighted parcels (significant and non-significant)
df_hcp_cen <- df_parcel %>%
  select(Arousal_State,lcen_1,lcen_6,rcen_5, rcen_6)


# Reshape to long format for sal, dmn, cen
df_hcp_melt <- melt(df_hcp_sal, id = "Arousal_State")
colnames(df_hcp_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_hcp_melt$parcel <- toupper(gsub("_", "", df_hcp_melt$parcel))


# Map parcels to brain regions
df_hcp_melt$region <- dplyr::recode(df_hcp_melt$parcel,
                                    "ASAL5" = "r-insula-IFO",
                                    "ASAL6" = "l-a-cerebellum",
                                    "ASAL7" = "r-a-cerebellum",
                                    "PSAL3" = "l-precuneus",
                                    "PSAL5" = "r-SPPP",
                                    "PSAL7" = "l-thal-hip",
                                    "PSAL8" = "l-p-cerebellum",
                                    "PSAL10" = "r-thalamus",
                                    "PSAL11" = "r-p-cerebellum",
                                    "PSAL12" = "r-insula"
                                    
                                    
                                    
)

# Ensure your main data's region column is a factor with proper order
df_hcp_melt$region <- factor(df_hcp_melt$region, levels = unique(df_hcp_melt$region))

# Check the exact factor levels
levels(df_hcp_melt$region)

# Define star positions (lifted above violins)

# Define star positions (lifted above violins)
star_positions <- data.frame(
  region = c("r-insula-IFO", "l-a-cerebellum", "r-a-cerebellum", "l-precuneus", 
             "r-SPPP", "l-thal-hip","l-p-cerebellum",
             "r-thalamus","r-p-cerebellum","r-insula"),
  y = 0.08,
  label = c("*", "**", "**", "*", "*","**","**","**","*","*"),
  bar = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE) ,
  bar2=c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
)

# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_hcp_melt$region))




# Plot
f2_a<-ggplot(df_hcp_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  theme_classic(base_size = 20) +
  labs( x = "Parcel", y = "SAL Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  # Stars (now using region directly)
  geom_text(data = star_positions, aes(x = region, y = y, label = label),
            inherit.aes = FALSE, size = 10, color="green3") +
  
  # Horizontal brace
  geom_segment(data = star_positions, 
               aes(x = as.numeric(factor(region)) - 0.2,
                   xend = as.numeric(factor(region)) + 0.2,
                   y = y - 0.004, yend = y - 0.004),
               inherit.aes = FALSE, size = 1, color="green3") +
  
  # Horizontal brace
  geom_segment(data = subset(star_positions, bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, linewidth = 1,color="black") +
  
  geom_segment(data = subset(star_positions, bar2), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, linewidth = 1,color="red") 

f2_a

##REPEAT FOR DMN
df_hcp_melt <- melt(df_hcp_dmn, id = "Arousal_State")
colnames(df_hcp_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_hcp_melt$parcel <- toupper(gsub("_", "", df_hcp_melt$parcel))


# Map parcels to brain regions
df_hcp_melt$region <- dplyr::recode(df_hcp_melt$parcel,
                                    "DDMN1" = "FMS-AC-R",
                                    "DDMN2" = "l-angular",
                                    "DDMN4" = "bi-Prec-Cing-Calc",
                                    "DDMN6" = "r-angular",
                                    "DDMN7" = "bi-thalamus",
                                    "DDMN8" = "l-hippo-fusiform",
                                    "VDMN1" = "l-Pre-Calc-Ling",
                                    "VDMN2" = "l-Frontal-S-M",
                                    "VDMN5" = "r-Pre-Calc-Ling",
                                    "VDMN7" = "r-Frontal-S-M"
                                    
                                    
                                    
)

# Ensure your main data's region column is a factor with proper order
df_hcp_melt$region <- factor(df_hcp_melt$region, levels = unique(df_hcp_melt$region))

# Check the exact factor levels
levels(df_hcp_melt$region)



# Define star positions (lifted above violins)
star_positions <- data.frame(
  region = c("FMS-AC-R", "l-angular", "bi-Prec-Cing-Calc", "r-angular", 
             "bi-thalamus", "l-hippo-fusiform","l-Pre-Calc-Ling","l-Frontal-S-M","r-Pre-Calc-Ling",
             "r-Frontal-S-M"),
  y = 0.085,
  label = c("*", "*", "**", "**", "**", "*","*","*","*","*"),
  bar = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
  bar2 = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
)
# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_hcp_melt$region))

# Plot

f2_c<-ggplot(df_hcp_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  theme_classic(base_size = 20) +
  labs(title = "HCP-7T: Parcel-Level Switching", x = "Parcel", y = " DMN Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  # Stars (now using region directly)
  geom_text(data = star_positions, aes(x = region, y = y, label = label),
            inherit.aes = FALSE, size = 10, color="green3") +
  
  # Horizontal brace
  geom_segment(data = star_positions, 
               aes(x = as.numeric(factor(region)) - 0.2,
                   xend = as.numeric(factor(region)) + 0.2,
                   y = y - 0.004, yend = y - 0.004),
               inherit.aes = FALSE, linewidth = 1,color="green3") +
  
  # Horizontal brace
  geom_segment(data = subset(star_positions, bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, linewidth = 1,color="black")+
  
  geom_segment(data = subset(star_positions, bar2), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, linewidth = 1,color="red") 


f2_c



df_hcp_melt <- melt(df_hcp_cen, id = "Arousal_State")
colnames(df_hcp_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_hcp_melt$parcel <- toupper(gsub("_", "", df_hcp_melt$parcel))

# Map parcels to brain regions
df_hcp_melt$region <- dplyr::recode(df_hcp_melt$parcel,
                                    "LCEN1" = "l-frontal-pre",
                                    "LCEN6" = "l-thalamus",
                                    "RCEN5" = "l-cerebellum",
                                    "RCEN6" = "r-thal-caudate"
                                    
                                    
                                    
)

# Ensure your main data's region column is a factor with proper order
df_hcp_melt$region <- factor(df_hcp_melt$region, levels = unique(df_hcp_melt$region))

# Check the exact factor levels
levels(df_hcp_melt$region)



# Define star positions (lifted above violins)
star_positions <- data.frame(
  region = c("l-frontal-pre", "l-thalamus", "l-cerebellum", 
             "r-thal-caudate"),
  y = 0.085,
  label = c("*", "*", "**", "**"),
  bar = c(FALSE, TRUE, TRUE, FALSE),
  bar2= c(FALSE,FALSE,FALSE,TRUE)
)
# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_hcp_melt$region))

# Plot

f2_e<-ggplot(df_hcp_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  theme_classic(base_size = 20) +
  labs( x = "Parcel", y = "CEN Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  # Stars (now using region directly)
  geom_text(data = star_positions, aes(x = region, y = y, label = label),
            inherit.aes = FALSE, size = 10, color="green3") +
  
  # Horizontal brace
  geom_segment(data = star_positions, 
               aes(x = as.numeric(factor(region)) - 0.2,
                   xend = as.numeric(factor(region)) + 0.2,
                   y = y - 0.004, yend = y - 0.004),
               inherit.aes = FALSE, size = 1,color="green3") +
  
  # Horizontal brace
  geom_segment(data = subset(star_positions, bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, size = 1,color="black") +
  
  
  geom_segment(data = subset(star_positions, bar2), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, size = 1, color="red3") 



f2_e

((f2_c ) / (f2_a)/ (f2_e)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )



##Repeat for VU-EEG-fMRI data

df_net <- read.csv("FIND_atlas_net_sr_72TR_5_27_25.csv")
df_net <- recode_arousal(df_net)

df_net<- df_net %>% rename(c(Arousal_State=arousal_state, scan_ID=ID,ID = subj))
df_net <- recode_arousal(df_net)

# Run all models
vars <- c("a_sal", "p_sal", "d_dmn", "v_dmn", "r_cen", "l_cen")
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"VU-EEG-fMRI-SR-networks-9-11-26.csv")

df_VU_sub <- df_net %>%
  select(Arousal_State, a_sal, p_sal, d_dmn, v_dmn, l_cen, r_cen)

# Reshape to long format
df_VU_melt <- melt(df_VU_sub, id = "Arousal_State")
colnames(df_VU_melt) <- c("Arousal_State", "network", "flexibility")

# Clean up network names
df_VU_melt$network <- toupper(gsub("_", "", df_VU_melt$network))
df_VU_melt$network <- factor(df_VU_melt$network, 
                             levels = c("ASAL", "PSAL", "DDMN", "VDMN", "LCEN", "RCEN"))

# Define significant networks for stars
star_positions <- data.frame(
  network = c("ASAL", "DDMN", "VDMN"),
  x = match(c("ASAL", "DDMN", "VDMN"), levels(df_VU_melt$network)),
  y = 0.041,  # lifted so stars/braces clear the violins
  bar= c(FALSE,TRUE,FALSE)
)

# Plot
f1_b<-ggplot(df_VU_melt, aes(x = network, y = flexibility, fill = Arousal_State)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  theme_classic(base_size = 16) +
  labs(title = "VU-EEG-fMRI: Network-Level Switching", 
       x = "Network", y = "Switching") +
  theme(text=element_text(size=20), axis.text.x = element_text( angle = 45, hjust = 1)) +
  scale_y_continuous(
    limits = y_limits,
    breaks = y_breaks,
    expand = expansion(mult = c(0.25, 0.15)))+
  geom_segment(data = subset(star_positions,bar), 
               aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y), 
               inherit.aes = FALSE, size = 1,color="black") 


library(patchwork)

((f1_a ) / (f1_b)) +
  plot_annotation(tag_levels = "a")
#Repeat for parcels


df_parcel<-read.csv("FIND_atlas_parcel_sr_72TR_5_27_25.csv")
df_parcel <- recode_arousal(df_parcel)

df_parcel<- df_parcel %>% rename(c( scan_ID=ID,ID = subj))
df_parcel <- recode_arousal(df_parcel)

all_cols <- colnames(df_parcel)
parcels <- all_cols[grep("dmn|sal|cen", all_cols)]
parcels <- parcels[2:51]

vars<- parcels

# Run all models
results <- bind_rows(lapply(vars, run_lmm_parcel))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"VU-EEG-fMRI-SR-parcels-9-11-26.csv")

#Make final parcel figure

# Subset highlighted parcels (significant and non-significant)
df_VU_sal <- df_parcel %>%
  select(Arousal_State, a_sal_5, a_sal_6, a_sal_7, p_sal_3, p_sal_5, p_sal_7,p_sal_8, p_sal_10, p_sal_11, p_sal_12)

# Subset highlighted parcels (significant and non-significant)
df_VU_dmn <- df_parcel %>%
  select(Arousal_State,d_dmn_1,d_dmn_2,d_dmn_4,d_dmn_6,d_dmn_7,d_dmn_8,v_dmn_1,v_dmn_2,v_dmn_5, v_dmn_7)

# Subset highlighted parcels (significant and non-significant)
df_VU_cen <- df_parcel %>%
  select(Arousal_State,lcen_1,lcen_6,rcen_5, rcen_6)


# Reshape to long format for sal, dmn, cen
df_VU_melt <- melt(df_VU_sal, id = "Arousal_State")
colnames(df_VU_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_VU_melt$parcel <- toupper(gsub("_", "", df_VU_melt$parcel))


# Map parcels to brain regions
df_VU_melt$region <- dplyr::recode(df_VU_melt$parcel,
                                   "ASAL5" = "r-insula-IFO",
                                   "ASAL6" = "l-a-cerebellum",
                                   "ASAL7" = "r-a-cerebellum",
                                   "PSAL3" = "l-precuneus",
                                   "PSAL5" = "r-SPPP",
                                   "PSAL7" = "l-thal-hip",
                                   "PSAL8" = "l-p-cerebellum",
                                   "PSAL10" = "r-thalamus",
                                   "PSAL11" = "r-p-cerebellum",
                                   "PSAL12" = "r-insula"
                                   
                                   
                                   
)

# Ensure your main data's region column is a factor with proper order
df_VU_melt$region <- factor(df_VU_melt$region, levels = unique(df_VU_melt$region))

# Check the exact factor levels
levels(df_VU_melt$region)

star_positions <- data.frame(
  region = c("r-insula-IFO", "l-a-cerebellum", "r-a-cerebellum", "l-precuneus", "r-SPPP", "l-thal-hip","l-p-cerebellum",
             "r-thalamus","r-p-cerebellum","r-insula"),
  y = 0.08,
  label = c("", "", "", "", "","","","","",""),
  bar = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
)

# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_VU_melt$region))


# Plot
f2_b<-ggplot(df_VU_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  theme_classic(base_size = 20) +
  labs(title = "VU-EEG-fMRI : Parcel-Level Switching", x = "Parcel", y = "SAL Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  
  # Horizontal brace
  geom_segment(data = subset(star_positions, bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, size = 1,color="black") 

##REPEAT FOR DMN
df_VU_melt <- melt(df_VU_dmn, id = "Arousal_State")
colnames(df_VU_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_VU_melt$parcel <- toupper(gsub("_", "", df_VU_melt$parcel))


# Map parcels to brain regions
df_VU_melt$region <- dplyr::recode(df_VU_melt$parcel,
                                   "DDMN1" = "FMS-AC-R",
                                   "DDMN2" = "l-angular",
                                   "DDMN4" = "bi-Prec-Cing-Calc",
                                   "DDMN6" = "r-angular",
                                   "DDMN7" = "bi-thalamus",
                                   "DDMN8" = "l-hippo-fusiform",
                                   "VDMN1" = "l-Pre-Calc-Ling",
                                   "VDMN2" = "l-Frontal-S-M",
                                   "VDMN5" = "r-Pre-Calc-Ling",
                                   "VDMN7" = "r-Frontal-S-M"
                                   
                                   
                                   
)

# Ensure your main data's region column is a factor with proper order
df_VU_melt$region <- factor(df_VU_melt$region, levels = unique(df_VU_melt$region))

# Check the exact factor levels
levels(df_VU_melt$region)



# Define star positions (lifted above violins)
star_positions <- data.frame(
  region = c("FMS-AC-R", "l-angular", "bi-Prec-Cing-Calc", "r-angular", 
             "bi-thalamus", "l-hippo-fusiform","l-Pre-Calc-Ling","l-Frontal-S-M","r-Pre-Calc-Ling",
             "r-Frontal-S-M"),
  y = 0.085,
  label = c("", "", "", "", "*", "","","","",""),
  bar = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE,FALSE, FALSE, FALSE) 
)
# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_VU_melt$region))

# Plot

f2_d<-ggplot(df_VU_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  
  theme_classic(base_size = 20) +
  labs(title = "VU-EEG-fMRI: Parcel-Level Switching", x = "Parcel", y = "DMN Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  # Stars (now using region directly)
  geom_text(data = star_positions, aes(x = region, y = y, label = label),
            inherit.aes = FALSE, size = 10, color="green3") +
  
  # Horizontal brace
  geom_segment(data = subset(star_positions,bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.004, yend = y - 0.004),
               inherit.aes = FALSE, size = 1,color="green3") +
  
  # Horizontal brace
  geom_segment(data = subset(star_positions, bar), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, size = 1,color="red") 






df_VU_melt <- melt(df_VU_cen, id = "Arousal_State")
colnames(df_VU_melt) <- c("arousal_state", "parcel", "flexibility")

# Clean parcel labels
df_VU_melt$parcel <- toupper(gsub("_", "", df_VU_melt$parcel))

# Map parcels to brain regions
df_VU_melt$region <- dplyr::recode(df_VU_melt$parcel,
                                   "LCEN1" = "l-frontal-pre",
                                   "LCEN6" = "l-thalamus",
                                   "RCEN5" = "l-cerebellum",
                                   "RCEN6" = "r-thal-caudate"
                                   
                                   
                                   
)

# Ensure your main data's region column is a factor with proper order
df_VU_melt$region <- factor(df_VU_melt$region, levels = unique(df_VU_melt$region))

# Check the exact factor levels
levels(df_VU_melt$region)


# Define star positions (lifted above violins)
star_positions <- data.frame(
  region = c("l-frontal-pre", "l-thalamus", "l-cerebellum", 
             "r-thal-caudate"),
  y = 0.085,
  label = c("", "", "", ""), 
  bar2= c(FALSE,TRUE,FALSE,TRUE)
)
# Add numeric x positions for plotting
star_positions$x <- match(star_positions$region, levels(df_VU_melt$region))

# Plot

f2_f<-ggplot(df_VU_melt, aes(x = region, y = flexibility, fill = arousal_state)) +
  geom_violin(trim = FALSE, alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA) +
  scale_fill_manual(values = c("High Arousal" = "#FF1493", "Low Arousal" = "purple2"),
                    name = "Arousal State") +
  scale_y_continuous(
    breaks = seq(-0.03, 0.09, by = 0.03),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  coord_cartesian(ylim = c(-0.03, 0.09))+
  theme_classic(base_size = 20) +
  labs(title = "VU-EEG-fMRI: Parcel-Level Switching", x = "Parcel", y = "CEN Switching") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  
  geom_segment(data = subset(star_positions, bar2), 
               aes(x = x - 0.2,
                   xend = x + 0.2,
                   y = y - 0.006, yend = y - 0.006),
               inherit.aes = FALSE, size = 1, color="red3") 


library(patchwork)

f2_a <- f2_a + theme(legend.position = "none")
f2_a <- f2_a + theme(plot.title = element_blank())

f2_b <- f2_b + theme(legend.position = "none")
f2_b <- f2_b + theme(plot.title = element_blank())

f2_c <- f2_c + theme(legend.position = "none")
#f2_c <- f2_c + theme(plot.title = element_blank())

f2_d <- f2_d + theme(legend.position = "none")
#f2_d <- f2_d + theme(plot.title = element_blank())

f2_e <- f2_e + theme(legend.position = "none")
f2_e <- f2_e + theme(plot.title = element_blank())

#f2_f <- f2_f + theme(legend.position = "none")
f2_f <- f2_f + theme(plot.title = element_blank())



#( f2_a | f2_b )+plot_annotation("a")

#( f2_c | f2_d )+plot_annotation("b")

#( f2_e | f2_f )+plot_annotation("c")


((f2_c | f2_d) / ( f2_a | f2_b)/ (f2_e | f2_f)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

#ANALYSIS 2 static correlation of network to global mean correlation across arousal state

hcp_df <- read.csv("HCP_7T_measures_of_interest.csv")
hcp_df <- recode_arousal(hcp_df)

vars<-c("gs_ddmn","gs_vdmn","gs_asal","mean_stat_corr")

run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal+ (1 |ID)")
  )
  
  model <- lmer(formula, data = hcp_df)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousalLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"HCP_7T_global_static_corr_9_11_26.csv")

#REDO FOR VU-EEG-FMRI

vu_df <- read.csv("eeg_fmri_vu_measures_of_interest.csv")
vu_df <- recode_arousal(vu_df)

vars<-c("gs_ddmn","gs_vdmn","gs_asal","mean_stat_corr")

run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal_state+ (1 |ID)")
  )
  
  model <- lmer(formula, data = vu_df)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousal_stateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"VU_eeg_fMRI_global_static_corr_9_11_26.csv")

#Now we will conduct the community assignment analysis

community_assign_mat_fin<-read.csv("HCP_7T_community_allegience.csv")
community_assign_mat_fin <- recode_arousal(community_assign_mat_fin)
all_cols <- colnames(community_assign_mat_fin)
values <- all_cols[grep("ASAL|PSAL", all_cols)]


vars<- values
run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal_state + (1 |ID)")
  )
  
  model <- lmer(formula, data = community_assign_mat_fin)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousal_stateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"HCP-7T_community_aliegence_across_arousal-9-11-26.csv")

#MAKE FIGURES

df_melt <- melt(community_assign_mat_fin,id.vars=c("arousal_state","ID","scan"), variable.name="Community_Assignment_Pair",value.name = "Fraction")

results <- df_melt %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(
    p_value = {
      model <- lmer(Fraction ~ arousal_state + (1 | ID), data = cur_data())
      summary(model)$coefficients["arousal_stateLow Arousal", "Pr(>|t|)"]
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = ifelse(p_adj < 0.05, "*", "")
  )# You can add "**" for lower thresholds if you like


df_plot <- df_melt %>%
  left_join(results, by = "Community_Assignment_Pair")

# Compute Y positions for stars (just above max of each group)
label_positions <- df_plot %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(y_pos = max(Fraction) + 0.07)  # Adjust the offset if needed

# Combine with significance results
star_labels <- left_join(results, label_positions, by = "Community_Assignment_Pair")

f4_a <-ggplot(df_plot, aes(x = Community_Assignment_Pair, y = Fraction, color = arousal_state)) +
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
            size = 12, color = "green3") +
  scale_color_manual(values = c("Low Arousal" = "purple", "High Arousal" = "hotpink")) +
  theme_classic() +
  theme(
    legend.title= element_text(size=14),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  labs(x = "Community Assignment Pair", y = "Fraction of Community Allegiance", title="HCP-7T")+ ylim(0,1.09)

#Repeat for EEG-fMRI

community_assign_mat_fin<-read.csv("eeg_fmri_community_allegience_2_9_26.csv")
community_assign_mat_fin <- recode_arousal(community_assign_mat_fin)
all_cols <- colnames(community_assign_mat_fin)
values <- all_cols[grep("ASAL|PSAL", all_cols)]


vars<- values
run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal_state + (1 |ID)")
  )
  
  model <- lmer(formula, data = community_assign_mat_fin)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousal_stateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"VU_EEG_fMRI_community_across_arousal-7-20-26.csv")

#MAKE FIGURES

df_melt <- melt(community_assign_mat_fin,id.vars=c("arousal_state","ID","scan"), variable.name="Community_Assignment_Pair",value.name = "Fraction")

results <- df_melt %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(
    p_value = {
      model <- lmer(Fraction ~ arousal_state + (1 | ID), data = cur_data())
      summary(model)$coefficients["arousal_stateLow Arousal", "Pr(>|t|)"]
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = ifelse(p_adj < 0.05, "*", "")
  )# You can add "**" for lower thresholds if you like


df_plot <- df_melt %>%
  left_join(results, by = "Community_Assignment_Pair")

# Compute Y positions for stars (just above max of each group)
label_positions <- df_plot %>%
  group_by(Community_Assignment_Pair) %>%
  summarise(y_pos = max(Fraction) + 0.07)  # Adjust the offset if needed

# Combine with significance results
star_labels <- left_join(results, label_positions, by = "Community_Assignment_Pair")

f4_b <-ggplot(df_plot, aes(x = Community_Assignment_Pair, y = Fraction, color = arousal_state)) +
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
            size = 12, color = "green3") +
  scale_color_manual(values = c("Low Arousal" = "purple", "High Arousal" = "hotpink")) +
  theme_classic() +
  theme(
    legend.title= element_text(size=14),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  labs(x = "Community Assignment Pair", y = "Fraction of Community Allegiance",title ="VU-EEG-fMRI")+ylim(0,1.09)

#Create Heatmap Figures

alert_average_matrix<-as.matrix(read.csv("alertcommunity_assignment_HCP_7T.csv",row.names = 1))
drowsy_average_matrix<-as.matrix(read.csv("drowsy_community_assignment_HCP_7T.csv",row.names= 1))


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
f4_c<-ggplot() +
  geom_tile(data = drowsy_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = drowsy_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white",
                      high = "purple", 
                      name = "Low Arousal",
                      limits=c(0,0.9), 
                      breaks=c(0, 0.30,0.60,0.90), 
                      labels=c("0%","30%","60%","90%")) +
  
  new_scale_fill() +
  
  geom_tile(data = alert_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = alert_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white", 
                      high = "hotpink", 
                      name = "High Arousal", 
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

#REPEAT FOR VU

alert_average_matrix<-as.matrix(read.csv("VU_Alert_Average_Matrix",row.names = 1))
drowsy_average_matrix<-as.matrix(read.csv("VU_Drowsy_Average_Matrix",row.names= 1))


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
f4_d<-ggplot() +
  geom_tile(data = drowsy_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = drowsy_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white",
                      high = "purple", 
                      name = "Low Arousal",
                      limits=c(0,0.9), 
                      breaks=c(0, 0.30,0.60,0.90), 
                      labels=c("0%","30%","60%","90%")) +
  
  new_scale_fill() +
  
  geom_tile(data = alert_tri_sub, aes(x = x, y = y, fill = Fraction), color = "black") +
  geom_text(data = alert_tri_sub, aes(x = x, y = y, label = Percentage), color = "black", size = 6) +
  scale_fill_gradient(low = "white", 
                      high = "hotpink", 
                      name = "High Arousal", 
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

((f4_a | f4_b) / ( f4_c | f4_d)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

#Global brain switching 
df_net_hcp<-read.csv("Net_Flex_MSTR.csv")
df_net_hcp <- recode_arousal(df_net_hcp)

df_net_hcp$global_sr <- apply(df_net_hcp[,4:17], 1, mean, na.rm = TRUE)

vars<- c("global_sr")

run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ Arousal_State + (1 |ID)")
  )
  
  model <- lmer(formula, data = df_net_hcp)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "Arousal_StateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

write.csv(results,"HCP_7T_GlobalSR_9-11-26.csv")

#Repeat for VU-EEG-fMRI

df_net_vu<-read.csv("FIND_atlas_net_sr_72TR_5_27_25.csv")
df_net_vu <- recode_arousal(df_net_vu)
df_net_vu$global_sr <- apply(df_net_vu[,6:19], 1, mean, na.rm = TRUE)

vars<- c("global_sr")

run_lmm_net <- function(var){
  
  formula <- as.formula(
    paste0(var, " ~ arousal_state + (1 |subj)")
  )
  
  model <- lmer(formula, data = df_net_vu)
  
  broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "arousal_stateLow Arousal") %>%   # replace with your actual factor level if different
    transmute(
      Variable = var,
      Beta = estimate,
      SE = std.error,
      df = df,
      t = statistic,
      p = p.value
    )
}

# Run all models
results <- bind_rows(lapply(vars, run_lmm_net))

# FDR correction
results$p_fdr <- p.adjust(results$p, method = "fdr")

results

#write.csv(results,"VU_EEG_fMRI_GlobalSR-9-11-26.csv")



