# ============================================================
# Supplementary Code to meet Reviewer 1 comments
# LM analysis looking at dynamic relationships between
# network switching rate and arousal levels
# ============================================================

# Packages
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(lme4)
library(lmerTest)


# ============================================================
# LINEAR MIXED-EFFECTS MODEL
# ============================================================

lm_test <- function(data, x, y, id = "ID") {
  
  # Keep only variables needed for this model
  d <- data %>%
    select(all_of(c(id, x, y))) %>%
    drop_na()
  
  # Make sure ID is a factor
  d[[id]] <- factor(d[[id]])
  
  # Formula:
  # outcome ~ predictor + random intercept for subject
  formula <- as.formula(
    paste0(y, " ~ ", x, " + (1 | ", id, ")")
  )
  
  # Fit model
  model <- lmer(
    formula,
    data = d,
    REML = FALSE
  )
  
  # Extract coefficient table
  coef_table <- summary(model)$coefficients
  
  tibble(
    outcome = y,
    predictor = x,
    N = nrow(d),
    n_subjects = n_distinct(d[[id]]),
    beta = coef_table[x, "Estimate"],
    SE = coef_table[x, "Std. Error"],
    t = coef_table[x, "t value"],
    p = coef_table[x, "Pr(>|t|)"]
  )
}


# ============================================================
# QUADRATIC MIXED-EFFECTS MODEL
# ============================================================

quad_test <- function(data, x, y, id = "ID") {
  
  # Keep only variables needed
  d <- data %>%
    select(all_of(c(id, x, y))) %>%
    drop_na()
  
  # Make ID a factor
  d[[id]] <- factor(d[[id]])
  
  # Center predictor
  d$x_centered <- d[[x]] - mean(d[[x]], na.rm = TRUE)
  
  # Quadratic model
  formula <- as.formula(
    paste0(
      y,
      " ~ x_centered + I(x_centered^2) + (1 | ",
      id,
      ")"
    )
  )
  
  # Fit model
  model <- lmer(
    formula,
    data = d,
    REML = FALSE
  )
  
  # Extract coefficients
  coef_table <- summary(model)$coefficients
  
  tibble(
    outcome = y,
    predictor = x,
    N = nrow(d),
    n_subjects = n_distinct(d[[id]]),
    
    # Linear component
    linear_beta = coef_table["x_centered", "Estimate"],
    linear_SE = coef_table["x_centered", "Std. Error"],
    linear_t = coef_table["x_centered", "t value"],
    linear_p = coef_table["x_centered", "Pr(>|t|)"],
    
    # Quadratic component
    quadratic_beta = coef_table["I(x_centered^2)", "Estimate"],
    quadratic_SE = coef_table["I(x_centered^2)", "Std. Error"],
    quadratic_t = coef_table["I(x_centered^2)", "t value"],
    quadratic_p = coef_table["I(x_centered^2)", "Pr(>|t|)"]
  )
  
  
}


# ============================================================
# LINEAR VS. QUADRATIC MODEL COMPARISON
# ============================================================

lm_vs_quad_test <- function(data, x, y, id = "ID") {
  
  # Use exactly the same observations for both models
  d <- data %>%
    select(all_of(c(id, x, y))) %>%
    drop_na()
  
  # Make ID a factor
  d[[id]] <- factor(d[[id]])
  
  # Center predictor
  d$x_centered <- d[[x]] - mean(d[[x]], na.rm = TRUE)
  
  # -------------------------
  # Linear model
  # -------------------------
  
  linear_formula <- as.formula(
    paste0(
      y,
      " ~ x_centered + (1 | ",
      id,
      ")"
    )
  )
  
  linear_model <- lmer(
    linear_formula,
    data = d,
    REML = FALSE
  )
  
  
  # -------------------------
  # Quadratic model
  # -------------------------
  
  quadratic_formula <- as.formula(
    paste0(
      y,
      " ~ x_centered + I(x_centered^2) + (1 | ",
      id,
      ")"
    )
  )
  
  quadratic_model <- lmer(
    quadratic_formula,
    data = d,
    REML = FALSE
  )
  
  
  # -------------------------
  # Likelihood-ratio test
  # -------------------------
  
  comparison <- anova(
    linear_model,
    quadratic_model
  )
  
  comparison_p <- comparison$`Pr(>Chisq)`[2]
  
  
  # -------------------------
  # AIC
  # -------------------------
  
  linear_AIC <- AIC(linear_model)
  quadratic_AIC <- AIC(quadratic_model)
  
  
  # -------------------------
  # Preferred model
  # -------------------------
  
  preferred_model <- case_when(
    comparison_p < 0.05 ~ "Quadratic",
    TRUE ~ "Linear"
  )
  
  
  # -------------------------
  # Return results
  # -------------------------
  
  tibble(
    outcome = y,
    predictor = x,
    N = nrow(d),
    n_subjects = n_distinct(d[[id]]),
    
    linear_AIC = linear_AIC,
    quadratic_AIC = quadratic_AIC,
    
    LR_chisq = comparison$Chisq[2],
    df = comparison$Df[2],
    model_comparison_p = comparison_p,
    
    preferred_model = preferred_model
  )
}


# ============================================================
# LOAD HCP-7T-DATA
# ============================================================

##Eyeclosures n/blinks
df_sr<-read.csv("HCP_7T_dynamic_sr_2.5min_4windows_all.csv")
df_sr_intermed<-read.csv("Inter_HCP_7T_dynamic_sr_2.5min_4windows_all.csv")
df_sr<-rbind(df_sr,df_sr_intermed)
df_sr$global_sr <- rowMeans(df_sr[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)

df_eye_intermed<-read.csv("hcp_rs_eye_measures_intermed.csv")
df_eye_alert<-read.csv("hcp_rs_eye_measures_alert.csv")
df_eye_drowsy<-read.csv("hcp_rs_eye_measures_drowsy.csv")

df_eye<-rbind(df_eye_alert,df_eye_drowsy,df_eye_intermed)

df_long <- df_eye %>%
  pivot_longer(
    cols = c(starts_with("eyeclosure_wb_"),starts_with("eyeclosure_nb_"),starts_with("pupil_")),
    names_to = c(".value", "window"),
    names_pattern = "(eyeclosure_wb|eyeclosure_nb|pupil)_(\\d+)"
  ) %>%
  mutate(
    window = as.numeric(window)
  )

df<-left_join(df_long,df_sr,by=c("ID","Scan","window"))



# ============================================================
# CALCULATE GLOBAL SWITCHING RATE
# ============================================================

df$global_sr <- rowMeans(
  df[, c(
    "a_sal",
    "p_sal",
    "d_dmn",
    "v_dmn",
    "l_cen",
    "r_cen",
    "aud",
    "bas_g",
    "h_vis",
    "lang",
    "prec",
    "p_vis",
    "sens",
    "vis_s"
  )],
  na.rm = TRUE
)



# ============================================================
# VARIABLES TO TEST
# ============================================================

# Predictor
x_vars <- "eyeclosure_wb"

# Outcomes
y_vars <- c(
  "global_sr",
  "a_sal",
  "d_dmn",
  "v_dmn"
)


# ============================================================
# 1. LINEAR MIXED-EFFECTS MODELS
# ============================================================

linear_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
linear_results

linear_results$q<-p.adjust(linear_results$p)

#write.csv(linear_results,"eye_wb_net_lm_test_9_3_26.csv")


# ============================================================
# 2. QUADRATIC MIXED-EFFECTS MODELS
# ============================================================

quadratic_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
quadratic_results <- quadratic_results %>%
  mutate(
    linear_p_fdr = p.adjust(linear_p, method = "fdr"),
    quadratic_p_fdr = p.adjust(quadratic_p, method = "fdr")
  )

#write.csv(quadratic_results,"eye_wb_net_quad_test_9_3_26.csv")

# ============================================================
# 3. LINEAR VS QUADRATIC MODEL COMPARISON
# ============================================================

model_comparison_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_vs_quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
model_comparison_results

model_comparison_results$q<-p.adjust(model_comparison_results$model_comparison_p)

write.csv(model_comparison_results,"eye_wb_net_lm_vs_quad_test_9_3_26.csv")

# ============================================================
# VARIABLES TO TEST EYECLOSURE NO BLINKS
# ============================================================

# Predictor
x_vars <- "eyeclosure_nb"

# Outcomes
y_vars <- c(
  "global_sr",
  "a_sal",
  "d_dmn",
  "v_dmn"
)


# ============================================================
# 1. LINEAR MIXED-EFFECTS MODELS
# ============================================================

linear_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
linear_results

linear_results$q<-p.adjust(linear_results$p)

#write.csv(linear_results,"eye_nb_net_lm_test_9_3_26.csv")


# ============================================================
# 2. QUADRATIC MIXED-EFFECTS MODELS
# ============================================================

quadratic_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
quadratic_results

quadratic_results <- quadratic_results %>%
  mutate(
    linear_p_fdr = p.adjust(linear_p, method = "fdr"),
    quadratic_p_fdr = p.adjust(quadratic_p, method = "fdr")
  )

#write.csv(quadratic_results,"eye_nb_net_quad_test_9_3_26.csv")

# ============================================================
# 3. LINEAR VS QUADRATIC MODEL COMPARISON
# ============================================================

model_comparison_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_vs_quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
model_comparison_results

model_comparison_results$q<-p.adjust(model_comparison_results$model_comparison_p)

#write.csv(model_comparison_results,"eye_nb_net_lm_vs_quad_test_9_3_26.csv")


# ============================================================
# VARIABLES TO TEST EYECLOSURE NO BLINKS
# ============================================================

# Predictor
x_vars <- "eyeclosure_nb"

# Outcomes
y_vars <- c(
  "global_sr",
  "a_sal",
  "d_dmn",
  "v_dmn"
)


# ============================================================
# 1. LINEAR MIXED-EFFECTS MODELS
# ============================================================

linear_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
linear_results

linear_results$q<-p.adjust(linear_results$p)

#write.csv(linear_results,"eye_nb_net_lm_test_9_3_26.csv")


# ============================================================
# 2. QUADRATIC MIXED-EFFECTS MODELS
# ============================================================

quadratic_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
quadratic_results

quadratic_results <- quadratic_results %>%
  mutate(
    linear_p_fdr = p.adjust(linear_p, method = "fdr"),
    quadratic_p_fdr = p.adjust(quadratic_p, method = "fdr")
  )

#write.csv(quadratic_results,"eye_nb_net_quad_test_9_3_26.csv")

# ============================================================
# 3. LINEAR VS QUADRATIC MODEL COMPARISON
# ============================================================

model_comparison_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_vs_quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
model_comparison_results

model_comparison_results$q<-p.adjust(model_comparison_results$model_comparison_p)

#write.csv(model_comparison_results,"eye_nb_net_lm_vs_quad_test_9_3_26.csv")



# ============================================================
# LOAD VU-EEG-fMRI DATA
# ============================================================

df <- read.csv(
  "vigall_dynamic_sr_2.5min_4windows_all.csv"
)


# ============================================================
# CALCULATE GLOBAL SWITCHING RATE
# ============================================================

df$global_sr <- rowMeans(
  df[, c(
    "a_sal",
    "p_sal",
    "d_dmn",
    "v_dmn",
    "l_cen",
    "r_cen",
    "aud",
    "bas_g",
    "h_vis",
    "lang",
    "prec",
    "p_vis",
    "sens",
    "vis_s"
  )],
  na.rm = TRUE
)


# ============================================================
# CREATE SUBJECT ID
# ============================================================
df_ID<-read.table("scan_list.txt")
df$ID<-rep(df_ID$V1,each=4)

# ============================================================
# VARIABLES TO TEST
# ============================================================

# Predictor
x_vars <- "vigall"

# Outcomes
y_vars <- c(
  "global_sr",
  "a_sal",
  "d_dmn",
  "v_dmn"
)


# ============================================================
# 1. LINEAR MIXED-EFFECTS MODELS
# ============================================================

linear_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
linear_results

linear_results$q<-p.adjust(linear_results$p)

#write.csv(linear_results,"vigall_cont_net_lm_test_9_3_26.csv")


# ============================================================
# 2. QUADRATIC MIXED-EFFECTS MODELS
# ============================================================

quadratic_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
quadratic_results

quadratic_results <- quadratic_results %>%
  mutate(
    linear_p_fdr = p.adjust(linear_p, method = "fdr"),
    quadratic_p_fdr = p.adjust(quadratic_p, method = "fdr")
  )

#write.csv(quadratic_results,"vigall_cont_net_quad_test_9_3_26.csv")


# ============================================================
# 3. LINEAR VS QUADRATIC MODEL COMPARISON
# ============================================================

model_comparison_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_vs_quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
model_comparison_results

model_comparison_results$q<-p.adjust(model_comparison_results$model_comparison_p)

#write.csv(model_comparison_results,"vigall_cont_net_lm_vs_quad_test_9_3_26.csv")

#===============================================
#MAKE PLOTS

#================================================
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# ============================================================
# FUNCTION TO CREATE MIXED-MODEL PREDICTIONS
# ============================================================

make_predictions <- function(data, x, y, id = "ID") {
  
  # Keep complete observations
  d <- data %>%
    select(all_of(c(id, x, y))) %>%
    drop_na()
  
  # Make ID a factor
  d[[id]] <- factor(d[[id]])
  
  # Center predictor
  d$x_centered <- d[[x]] - mean(d[[x]], na.rm = TRUE)
  
  # Linear mixed model
  linear_model <- lmer(
    as.formula(
      paste0(y, " ~ x_centered + (1 | ", id, ")")
    ),
    data = d,
    REML = FALSE
  )
  
  # Quadratic mixed model
  quadratic_model <- lmer(
    as.formula(
      paste0(y, " ~ x_centered + I(x_centered^2) + (1 | ", id, ")")
    ),
    data = d,
    REML = FALSE
  )
  
  # Prediction range
  prediction_data <- data.frame(
    x_centered = seq(
      min(d$x_centered, na.rm = TRUE),
      max(d$x_centered, na.rm = TRUE),
      length.out = 200
    )
  )
  
  # Linear predictions
  prediction_data$linear <- predict(
    linear_model,
    newdata = prediction_data,
    re.form = NA
  )
  
  # Quadratic predictions
  prediction_data$quadratic <- predict(
    quadratic_model,
    newdata = prediction_data,
    re.form = NA
  )
  
  # Convert centered x back to original scale
  x_mean <- mean(d[[x]], na.rm = TRUE)
  
  prediction_data[[x]] <- prediction_data$x_centered + x_mean
  
  prediction_data
}

#======================================
#VIGALL PLOTS
#======================================

library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


# ============================================================
# FUNCTION TO CREATE MIXED-MODEL PREDICTIONS
# ============================================================

make_predictions <- function(data, x, y, id = "ID") {
  
  # Keep complete observations
  d <- data %>%
    select(all_of(c(id, x, y))) %>%
    drop_na()
  
  # Make ID a factor
  d[[id]] <- factor(d[[id]])
  
  # Center predictor
  d$x_centered <- d[[x]] - mean(d[[x]], na.rm = TRUE)
  
  # Linear mixed model
  linear_model <- lmer(
    as.formula(
      paste0(y, " ~ x_centered + (1 | ", id, ")")
    ),
    data = d,
    REML = FALSE
  )
  
  # Quadratic mixed model
  quadratic_model <- lmer(
    as.formula(
      paste0(y, " ~ x_centered + I(x_centered^2) + (1 | ", id, ")")
    ),
    data = d,
    REML = FALSE
  )
  
  # Prediction range
  prediction_data <- data.frame(
    x_centered = seq(
      min(d$x_centered, na.rm = TRUE),
      max(d$x_centered, na.rm = TRUE),
      length.out = 200
    )
  )
  
  # Linear predictions
  prediction_data$linear <- predict(
    linear_model,
    newdata = prediction_data,
    re.form = NA
  )
  
  # Quadratic predictions
  prediction_data$quadratic <- predict(
    quadratic_model,
    newdata = prediction_data,
    re.form = NA
  )
  
  # Convert centered x back to original scale
  x_mean <- mean(d[[x]], na.rm = TRUE)
  
  prediction_data[[x]] <- prediction_data$x_centered + x_mean
  
  prediction_data
}

#=====================================================
df_vigall <- read.csv("vigall_dynamic_sr_2.5min_4windows_all.csv")

df_vigall$global_sr <- rowMeans(
  df_vigall[, c(
    "a_sal", "p_sal", "d_dmn", "v_dmn",
    "l_cen", "r_cen", "aud", "bas_g",
    "h_vis", "lang", "prec", "p_vis",
    "sens", "vis_s"
  )],
  na.rm = TRUE
)

# Subject IDs: 4 observations per subject
df_ID <- read.table(
  "scan_list.txt")

df_vigall$ID <- rep(df_ID$V1, each = 4)

pred_asal <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "a_sal",
  id = "ID"
)

pred_ddmn <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "d_dmn",
  id = "ID"
)

pred_vdmn <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "v_dmn",
  id = "ID"
)

pred_global <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "global_sr",
  id = "ID"
)

pred_asal <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "a_sal",
  id = "ID"
)

pred_ddmn <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "d_dmn",
  id = "ID"
)

pred_vdmn <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "v_dmn",
  id = "ID"
)

pred_global <- make_predictions(
  df_vigall,
  x = "vigall",
  y = "global_sr",
  id = "ID"
)


F5_B <- ggplot(df_vigall, aes(x = vigall, y = a_sal)) +
  
  geom_jitter(
    width = 0.05,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_asal,
    aes(x = vigall, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_asal,
    aes(x = vigall, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.4)
  ) +
  
  annotate(
    "text",
    x = 4.0,
    y = 0.040,
    label = "*",
    size = 10,
    color = "blue"
  ) +
  
  labs(
    x = "EEG Arousal-VIGALL",
    y = "ASAL"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_E <- ggplot(df_vigall, aes(x = vigall, y = d_dmn)) +
  
  geom_jitter(
    width = 0.05,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_ddmn,
    aes(x = vigall, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_ddmn,
    aes(x = vigall, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.4)
  ) +
  
  labs(
    x = "EEG Arousal-VIGALL",
    y = "DDMN"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_H <- ggplot(df_vigall, aes(x = vigall, y = v_dmn)) +
  
  geom_jitter(
    width = 0.05,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_vdmn,
    aes(x = vigall, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_vdmn,
    aes(x = vigall, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.3)
  ) +
  
  labs(
    x = "EEG Arousal-VIGALL",
    y = "VDMN"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_K <- ggplot(df_vigall, aes(x = vigall, y = global_sr)) +
  
  geom_jitter(
    width = 0.05,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_global,
    aes(x = vigall, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_global,
    aes(x = vigall, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.5)
  ) +
  
  labs(
    x = "EEG Arousal-VIGALL",
    y = "Global-SR"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))



df_sr <- read.csv("HCP_7T_dynamic_sr_2.5min_4windows_all.csv")

df_sr_intermed <- read.csv("Inter_HCP_7T_dynamic_sr_2.5min_4windows_all.csv")

df_sr <- rbind(df_sr,df_sr_intermed)

df_sr$global_sr <- rowMeans(
  df_sr[, c(
    "a_sal", "p_sal", "d_dmn", "v_dmn",
    "l_cen", "r_cen", "aud", "bas_g",
    "h_vis", "lang", "prec", "p_vis",
    "sens", "vis_s"
  )],
  na.rm = TRUE
)


df_eye_intermed <- read.csv(
  "hcp_rs_eye_measures_intermed.csv"
)

df_eye_alert <- read.csv(
  "hcp_rs_eye_measures_alert.csv"
)

df_eye_drowsy <- read.csv(
  "hcp_rs_eye_measures_drowsy.csv"
)

df_eye <- rbind(
  df_eye_alert,
  df_eye_drowsy,
  df_eye_intermed
)


df_long <- df_eye %>%
  pivot_longer(
    cols = c(
      starts_with("eyeclosure_wb_"),
      starts_with("eyeclosure_nb_"),
      starts_with("pupil_")
    ),
    names_to = c(".value", "window"),
    names_pattern = "(eyeclosure_wb|eyeclosure_nb|pupil)_(\\d+)"
  ) %>%
  mutate(
    window = as.numeric(window)
  )


df_eye_final <- left_join(
  df_long,
  df_sr,
  by = c("ID", "Scan", "window"))
  
pred_eye_asal <- make_predictions(
  df_eye_final,
  x = "eyeclosure_nb",
  y = "a_sal",
  id = "ID"
)

pred_eye_ddmn <- make_predictions(
  df_eye_final,
  x = "eyeclosure_nb",
  y = "d_dmn",
  id = "ID"
)

pred_eye_vdmn <- make_predictions(
  df_eye_final,
  x = "eyeclosure_nb",
  y = "v_dmn",
  id = "ID"
)

pred_eye_global <- make_predictions(
  df_eye_final,
  x = "eyeclosure_nb",
  y = "global_sr",
  id = "ID"
)

F5_A <- ggplot(
  df_eye_final,
  aes(x = eyeclosure_nb, y = a_sal)
) +
  
  geom_jitter(
    width = 0.01,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_eye_asal,
    aes(x = eyeclosure_nb, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_eye_asal,
    aes(x = eyeclosure_nb, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.4)
  ) +
  
  labs(
    x = "Percent of Eye Closures",
    y = "ASAL"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_D <- ggplot(
  df_eye_final,
  aes(x = eyeclosure_nb, y = d_dmn)
) +
  
  geom_jitter(
    width = 0.01,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_eye_ddmn,
    aes(x = eyeclosure_nb, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_eye_ddmn,
    aes(x = eyeclosure_nb, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.4)
  ) +
  
  annotate(
    "text",
    x = 0.55,
    y = 0.04,
    label = "*",
    size = 10,
    color = "green"
  ) +
  
  labs(
    x = "Percent of Eye Closures",
    y = "DDMN"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_G <- ggplot(
  df_eye_final,
  aes(x = eyeclosure_nb, y = v_dmn)
) +
  
  geom_jitter(
    width = 0.01,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_eye_vdmn,
    aes(x = eyeclosure_nb, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_eye_vdmn,
    aes(x = eyeclosure_nb, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.5)
  ) +
  
  annotate(
    "text",
    x = 0.50,
    y = 0.045,
    label = "*",
    size = 10,
    color = "blue"
  ) +
  
  labs(
    x = "Percent of Eye Closure",
    y = "VDMN"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))


F5_J <- ggplot(
  df_eye_final,
  aes(x = eyeclosure_nb, y = global_sr)
) +
  
  geom_jitter(
    width = 0.01,
    height = 0,
    size = 3,
    alpha = 0.7
  ) +
  
  geom_line(
    data = pred_eye_global,
    aes(x = eyeclosure_nb, y = linear),
    linewidth = 2,
    linetype = "dashed",
    color = scales::alpha("green3", 0.6)
  ) +
  
  geom_line(
    data = pred_eye_global,
    aes(x = eyeclosure_nb, y = quadratic),
    linewidth = 2,
    color = scales::alpha("blue", 0.5)
  ) +
  
  labs(
    x = "Percent of Eye Closures",
    y = "Global-SR"
  ) +
  
  theme_classic() +
  theme(text = element_text(size = 20))

library(patchwork)

((F5_A | F5_B) /
    (F5_D | F5_E) /
    (F5_G | F5_H) /
    (F5_J | F5_K)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

#Supplementary test with pupil


# ============================================================
# Check Pupil (Only Eyes Open)
# ============================================================
df_eye_alert<-read.csv("hcp_rs_eye_measures_alert.csv")
df_long <- df_eye_alert %>%
  pivot_longer(
    cols = c(starts_with("eyeclosure_wb_"),starts_with("eyeclosure_nb_"),starts_with("pupil_")),
    names_to = c(".value", "window"),
    names_pattern = "(eyeclosure_wb|eyeclosure_nb|pupil)_(\\d+)"
  ) %>%
  mutate(
    window = as.numeric(window)
    
  )

df_sr<-read.csv("HCP_7T_dynamic_sr_2.5min_4windows_all.csv")

df_sr$global_sr <- rowMeans(df_sr[, c("a_sal", "p_sal", "d_dmn","v_dmn","l_cen","r_cen","aud","bas_g","h_vis","lang","prec","p_vis","sens","vis_s")], na.rm = TRUE)

df<-left_join(df_long,df_sr,by=c("ID","Scan","window"))

# Predictor
x_vars <- "pupil"

# Outcomes
y_vars <- c(
  "global_sr",
  "a_sal",
  "d_dmn",
  "v_dmn"
)


# ============================================================
# 1. LINEAR MIXED-EFFECTS MODELS
# ============================================================

linear_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
linear_results

linear_results$q<-p.adjust(linear_results$p)

#write.csv(linear_results,"pupil_net_lm_test_9_3_26.csv")


# ============================================================
# 2. QUADRATIC MIXED-EFFECTS MODELS
# ============================================================

quadratic_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
quadratic_results
quadratic_results <- quadratic_results %>%
  mutate(
    linear_p_fdr = p.adjust(linear_p, method = "fdr"),
    quadratic_p_fdr = p.adjust(quadratic_p, method = "fdr")
  )

#write.csv(quadratic_results,"pupil_net_quad_test_9_3_26.csv")

# ============================================================
# 3. LINEAR VS QUADRATIC MODEL COMPARISON
# ============================================================

model_comparison_results <- expand.grid(
  x = x_vars,
  y = y_vars,
  stringsAsFactors = FALSE
) %>%
  pmap_dfr(
    ~ lm_vs_quad_test(
      data = df,
      x = ..1,
      y = ..2,
      id = "ID"
    )
  )


# View results
model_comparison_results

model_comparison_results$q<-p.adjust(model_comparison_results$model_comparison_p)

#write.csv(model_comparison_results,"pupil_net_lm_vs_quad_test_9_3_26.csv")

