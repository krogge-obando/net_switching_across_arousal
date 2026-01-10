#This code will see if cog- relates to network-switch

#Load df

library(tidyverse)
library(dplyr)
library(stats)

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

global_network_flexibility<-df_filtered %>%
  group_by(ID) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

#derive global switching
global_network_flexibility$global_flex<-rowMeans(global_network_flexibility[,-1,2],na.rm=TRUE)

df_behavior$ID<-as.factor(df_behavior$ID)

df_unique <- df_filtered %>%
  distinct(ID, .keep_all = TRUE)

gnf_behav<-left_join(global_network_flexibility,df_behavior,by="ID")

gnf_behav$arousal_state<-df_unique$Arousal_State

# Make plots
library(ggplot2)

ggplot(gnf_behav, aes(x = global_flex, y = Relational_Task_Acc, color = arousal_state)) +
  geom_point(size=2) +  # Scatter plot
  geom_smooth(method = "lm") +  # Linear regression trend line
  scale_color_manual(values = c("alert" = "hotpink", "drowsy" = "purple")) +  # Custom colors
  theme_classic() +  # Clean theme
   theme(legend.text=element_text(size=14),
         legend.title=element_text(size=20),
    axis.text.x = element_text(size = 16),  # Bigger x-axis text
    axis.text.y = element_text(size = 16),  # Bigger y-axis text
    axis.title.x = element_text(size = 20), # Bigger x-axis title
    axis.title.y = element_text(size = 20)  # Bigger y-axis title
  ) + labs(x = "Global Brain Flexibility", y = "Relational Task Accuracy") + annotate(
    "text",
    x = mean(range(gnf_behav$global_flex)),    # center x
    y = 103,  # center y
    label = "*",
    size = 10,    # size of star
    color = "black"
  ) # Axis labels



model1<-summary(lm(Relational_Task_Acc ~ global_flex + arousal_state + global_flex:arousal_state, data=gnf_behav ))

ggplot(gnf_behav, aes(x = global_flex, y = WM_Task_Acc, color = arousal_state)) +
  geom_point(size=2) +  # Scatter plot
  geom_smooth(method = "lm") +  # Linear regression trend line
  scale_color_manual(values = c("alert" = "hotpink", "drowsy" = "purple")) +  # Custom colors
  theme_classic() +  # Clean theme
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=20),
    axis.text.x = element_text(size = 16),  # Bigger x-axis text
    axis.text.y = element_text(size = 16),  # Bigger y-axis text
    axis.title.x = element_text(size = 20), # Bigger x-axis title
    axis.title.y = element_text(size = 20)  # Bigger y-axis title
  ) + labs(x = "Global Brain Flexibility", y = "Working Memory Task Accuracy")  # Axis labels

model2<-summary(lm(WM_Task_Acc ~ global_flex + arousal_state + global_flex:arousal_state, data=gnf_behav ))


# Extract the p-values of the interaction term and conduct corrections
pvals <- c(
  model1$coefficients[4, 4],
  model2$coefficients[4, 4]
)

p.adjust(pvals, method = "BH")

summary(model1)$coefficients




