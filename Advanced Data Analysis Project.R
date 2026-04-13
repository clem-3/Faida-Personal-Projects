# Install 
required_packages <- c(
  "tidyverse", "lme4", "lmerTest", "performance", "patchwork",
  "gtsummary", "janitor", "skimr", "DHARMa", "broom.mixed",
  "emmeans", "ggpubr", "tableone", "naniar", "GGally",
  "ggcorrplot", "effectsize", "car", "reshape2", "viridis"
)

installed <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

# Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(performance)
library(patchwork)
library(gtsummary)
library(janitor)
library(skimr)
library(DHARMa)
library(broom.mixed)
library(emmeans)
library(ggpubr)
library(tableone)
library(naniar)
library(GGally)
library(ggcorrplot)
library(effectsize)
library(car)
library(reshape2)
library(viridis)

# Set global ggplot theme for publication quality
theme_set(
  theme_bw(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
      strip.background = element_rect(fill = "#2E75B6", colour = NA),
      strip.text    = element_text(colour = "white", face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
)

# Colour palette (treatment groups)
pal_trt <- c("0" = "#E05252", "1" = "#2E75B6")   # red = control, blue = eSMC

cat("\n DATA IMPORT & STRUCTURE \n")

# Import 
data <- read_csv("ADAM PROJECT/malaria_anemia_rct.csv", show_col_types = FALSE)

# Inspection
str(data)
cat("\n glimpse() \n"); glimpse(data)
cat("\n--- skimr::skim() ---\n"); print(skim(data))

# Missing value summary 
cat("\n Missing values \n")
print(data |> summarise(across(everything(), ~sum(is.na(.)))))

# SECTION 2 — DATA CLEANING

cat("\n DATA CLEANING \n")

df <- data   # working copy

# Factor conversion 
df <- df |>
  mutate(
    treatment = factor(treatment, levels = c(0, 1),
                       labels = c("Control", "eSMC")),
    itn_use   = factor(itn_use,   levels = c(0, 1),
                       labels = c("No", "Yes")),
    month_f   = factor(month, levels = c(0, 2, 4, 6),
                       labels = c("M0", "M2", "M4", "M6")),
    child_id  = factor(child_id)
  )

cat("\nData structure after factor conversion:\n"); glimpse(df)
str(df)
summary(df)

# Age distribution
age_hist <- ggplot(df |> distinct(child_id, .keep_all = TRUE),
                   aes(x = age_baseline)) +
  geom_histogram(binwidth = 6, fill = "#2E75B6", colour = "white") +
  labs(title = "Distribution of Baseline Age (months)",
       x = "Age (months)", y = "Count")

print(age_hist)

# Logical consistency checks
# Duplicate child × month combinations
dupl <- df |> group_by(child_id, month) |> filter(n() > 1)
cat("\nDuplicate child-month rows:", nrow(dupl), "\n")

# Check each child has ≤ 4 visits
visit_counts <- df |> count(child_id, name = "n_visits")
cat("Visit count distribution:\n"); print(table(visit_counts$n_visits))

# BASELINE CHARACTERISTICS TABLE

cat(" BASELINE CHARACTERISTICS \n")

df_baseline <- df |> filter(month == 0)

# Publication-ready table using gtsummary
tbl_baseline <- df_baseline |>
  select(treatment, age_baseline, itn_use, hemoglobin) |>
  tbl_summary(
    by = treatment,
    label = list(
      age_baseline ~ "Age at baseline (months)",
      itn_use      ~ "ITN use",
      hemoglobin   ~ "Baseline hemoglobin (g/dL)"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(all_continuous() ~ 1)
  ) |>
  add_p(test = list(all_continuous()  ~ "t.test",
                    all_categorical() ~ "chisq.test")) |>
  add_overall() |>
  modify_header(label ~ "**Characteristic**") |>
  modify_caption("**Table 1. Baseline Characteristics by Treatment Group**") |>
  bold_labels()

print(tbl_baseline)

# Hemoglobin Trajectories Over Time by Treatment Group

ggplot(df, aes(x = month, y = hemoglobin, group = child_id, color = treatment)) +
  geom_line(alpha = 0.2) +
  geom_smooth(aes(group = treatment), se = FALSE, size = 1.5) +
  labs(
    title = "Hemoglobin Trajectories Over Time by Treatment Group",
    x = "Month",
    y = "Hemoglobin (g/dL)",
    color = "Treatment"
  ) +
  theme_minimal()

# Histograms of hemoglobin by group and time

hist_overall <- ggplot(df, aes(x = hemoglobin, fill = treatment)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5, position = "identity", alpha = 0.6) +
  geom_density(aes(colour = treatment), linewidth = 1, fill = NA) +
  scale_fill_manual(values  = pal_trt, name = "Treatment") +
  scale_colour_manual(values = pal_trt, name = "Treatment") +
  labs(title = "Hemoglobin Distribution by Treatment Group (All Timepoints)",
       x = "Hemoglobin (g/dL)", y = "Density")

hist_time <- ggplot(df, aes(x = hemoglobin, fill = treatment)) +
  geom_histogram(binwidth = 0.5, position = "identity", alpha = 0.6) +
  facet_wrap(~ month_f, ncol = 2) +
  scale_fill_manual(values = pal_trt, name = "Treatment") +
  labs(title = "Hemoglobin Distribution by Visit and Treatment",
       x = "Hemoglobin (g/dL)", y = "Count")

print(hist_overall)
print(hist_time)

# Bar charts for categorical variables

bar_trt <- df_baseline |>
  count(treatment) |>
  ggplot(aes(x = treatment, y = n, fill = treatment)) +
  geom_col(width = 0.5, show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.4, fontface = "bold") +
  scale_fill_manual(values = pal_trt) +
  labs(title = "Sample Size by Treatment Group",
       x = "Treatment", y = "N")

bar_itn <- df_baseline |>
  count(treatment, itn_use) |>
  ggplot(aes(x = treatment, y = n, fill = itn_use)) +
  geom_col(position = "fill", width = 0.5) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2", name = "ITN Use") +
  labs(title = "ITN Use by Treatment Group",
       x = "Treatment", y = "Proportion")

print(bar_trt + bar_itn)

cat("\n DISTRIBUTION ANALYSIS \n")

p_dist_hist <- ggplot(df, aes(x = hemoglobin)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5, fill = "#2E75B6", colour = "white", alpha = 0.8) +
  geom_density(linewidth = 1.2, colour = "navy") +
  stat_function(fun = dnorm,
                args = list(mean = mean(df$hemoglobin, na.rm = TRUE),
                            sd   = sd(df$hemoglobin,   na.rm = TRUE)),
                linewidth = 1, colour = "red", linetype = "dashed") +
  labs(title = "Hemoglobin Distribution (All Observations)",
       subtitle = "Blue = empirical density | Red dashed = theoretical normal",
       x = "Hemoglobin (g/dL)", y = "Density")

p_qq <- ggplot(df, aes(sample = hemoglobin)) +
  stat_qq(colour = "#2E75B6", alpha = 0.5) +
  stat_qq_line(colour = "red", linewidth = 1) +
  labs(title = "Normal Q-Q Plot of Hemoglobin",
       subtitle = "Points close to line = good normality",
       x = "Theoretical Quantiles", y = "Sample Quantiles")

print(p_dist_hist + p_qq)

# QQ plot by treatment and time
p_qq_facet <- ggplot(df, aes(sample = hemoglobin, colour = treatment)) +
  stat_qq(alpha = 0.6) +
  stat_qq_line(aes(colour = treatment)) +
  facet_wrap(~ month_f) +
  scale_colour_manual(values = pal_trt, name = "Treatment") +
  labs(title = "Q-Q Plots by Visit and Treatment",
       x = "Theoretical Quantiles", y = "Sample Quantiles")
print(p_qq_facet)

### Model
# -----------------------------
# 11. EMPTY MODEL / ICC
# -----------------------------
cat("\n EMPTY MODEL / ICC \n")

model0 <- lmer(
  hemoglobin ~ 1 + (1 | child_id),
  data = df,
  REML = TRUE
)

print(summary(model0))

var_comp0 <- as.data.frame(VarCorr(model0))
sigma_child <- var_comp0$vcov[var_comp0$grp == "child_id"]
sigma_resid <- attr(VarCorr(model0), "sc")^2
ICC <- sigma_child / (sigma_child + sigma_resid)

cat("\nIntraclass Correlation Coefficient (ICC):\n")
print(ICC)


### Model with no interaction

model1 <- lmer(
  hemoglobin ~ treatment + month + age_baseline + itn_use + (1 | child_id),
  data = df,
  REML = TRUE    # REML preferred for variance component estimation
)

cat("\n LMM Summary \n")
print(summary(model1))

### Model with interaction 

model_main <- lmer(
  hemoglobin ~ treatment * month + age_baseline + itn_use + (1 | child_id),
  data = df,
  REML = TRUE    # REML preferred for variance component estimation
)

cat("\n  LMM Summary \n")
print(summary(model_main))


anova_comp0 <- anova(model1, model_main)
print(anova_comp0)

# Tidy output
fixed_effects <- tidy(model_main, conf.int = TRUE, effects = "fixed") |>
  mutate(across(where(is.numeric), ~round(., 4)))

cat("\n Fixed Effects (broom.mixed) \n")
print(fixed_effects)

random_effects <- tidy(model_main, effects = "ran_pars")
cat("\n Random Effects \n")
print(random_effects)

# Random slope model 
model_slope <- lmer(
  hemoglobin ~ treatment * month + age_baseline + itn_use +
    (month | child_id),
  data = df, REML = TRUE,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 2e5))
)

cat("\n Random Slope Model Summary \n")
print(summary(model_slope))

# Model comparison (AIC / BIC / LRT)
cat("\n Model Comparison (random intercept vs. random slope) \n")
anova_comp <- anova(model_main, model_slope)
print(anova_comp)

# MODEL DIAGNOSTICS

cat("\n MODEL DIAGNOSTICS \n")

#   residuals and random effects
resid_vals  <- residuals(m_final)
fitted_vals <- fitted(m_final)
re          <- ranef(m_final)$child_id

#  Residuals vs Fitted
p1 <- ggplot(data.frame(fitted = fitted_vals, resid = resid_vals),
             aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue",
              linewidth = 0.8, formula = y ~ x) +
  labs(x = "Fitted Values", y = "Residuals",
       title = "(A) Residuals vs. Fitted") +
  theme_minimal(base_size = 11)

# Q-Q — observation-level residuals
p2 <- ggplot(data.frame(resid = resid_vals), aes(sample = resid)) +
  stat_qq(alpha = 0.3, size = 0.7) +
  stat_qq_line(color = "red") +
  labs(title = "(B) Q-Q Plot: Level-1 Residuals") +
  theme_minimal(base_size = 11)

# Q-Q — random intercepts
p3 <- ggplot(data.frame(intercept = re[, 1]), aes(sample = intercept)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(color = "red") +
  labs(title = "(C) Q-Q Plot: Random Intercepts") +
  theme_minimal(base_size = 11)

#  Panel D: Q-Q — random slopes
p4 <- ggplot(data.frame(slope = re[, 2]), aes(sample = slope)) +
  stat_qq(alpha = 0.5) +
  stat_qq_line(color = "red") +
  labs(title = "(D) Q-Q Plot: Random Slopes") +
  theme_minimal(base_size = 11)

(p1 + p2) / (p3 + p4)

# --- 7.6 performance::check_model() 
cat("\n--- performance::check_model() ---\n")
check_model(model_slope)  # Opens interactive diagnostics panel

# Individual checks
cat("\n--- Normality of residuals:\n"); print(check_normality(model_slope))
cat("\n--- Homogeneity of variance:\n"); print(check_heteroscedasticity(model_slope))
cat("\n--- Random effects normality:\n"); print(check_normality(model_slope, effects = "random"))
