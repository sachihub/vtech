# =============================================================
# R translation of Stata script for Figure 1
# =============================================================

# ----------------- Setup
library(haven)
library(fixest)
library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(tidyr)
library(data.table)
library(broom)
library(patchwork)
library(stringr)

# -------------------------------------------------------------
# Load main data
df <- read_dta("output/full_data_for_regression.dta")
df$modate <- as.factor(df$modate)

# -------------------------------------------------------------
# Run fixed-effects model
m_linear <- feols(tdths_ttpop ~ . - modate | modate, data = df)
saveRDS(m_linear, "output/beta/m_linear.rds")

# -------------------------------------------------------------
# Predict values
df$y_imt_TC3 <- predict(m_linear)
df$modate_lbl <- make_date(df$year, df$month, 1)

# -------------------------------------------------------------
# Plot for NJ and FL
plot_pred_vs_actual <- function(state_code) {
  state_data <- df %>%
    filter(state_abbrev == state_code, year >= 1950)

  ggplot(state_data, aes(x = modate_lbl)) +
    geom_line(aes(y = tdths_ttpop), color = "orange", alpha = 0.5) +
    geom_point(aes(y = y_imt_TC3), color = "darkred", size = 0.5, alpha = 0.5) +
    labs(
      x = NULL,
      y = "Monthly All Cause Mortality Rate (per 100,000)",
      title = state_code
    ) +
    theme_minimal()
}

p_nj <- plot_pred_vs_actual("NJ")
p_fl <- plot_pred_vs_actual("FL")

ggsave("figures/state_NJ_pred.png", p_nj, width = 6, height = 4)
ggsave("figures/state_FL_pred.png", p_fl, width = 6, height = 4)

# -------------------------------------------------------------
# Collapse to yearly predictions
yearly <- df %>%
  group_by(year) %>%
  summarise(
    tdths_ttpop_year = weighted.mean(tdths_ttpop, totalpop, na.rm = TRUE),
    y_imt_TC3_year = weighted.mean(y_imt_TC3, totalpop, na.rm = TRUE)
  )

df <- left_join(df, yearly, by = "year")

# One obs per year
df <- df %>%
  group_by(year) %>%
  mutate(
    y_imt_TC3_year = if_else(row_number() == 1, y_imt_TC3_year, NA_real_),
    tdths_ttpop_year = if_else(row_number() == 1, tdths_ttpop_year, NA_real_)
  )

# R² value
r2 <- summary(m_linear)$r.squared

# -------------------------------------------------------------
# True vs Predicted Scatterplot
p_tp <- ggplot(df, aes(x = y_imt_TC3, y = tdths_ttpop)) +
  geom_point(alpha = 0.1, color = "darkred") +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    x = "Predicted Mortality Rate (per 100,000)",
    y = "True Mortality Rate (per 100,000)",
    caption = paste0("Model R² = ", round(r2, 3))
  ) +
  theme_minimal()

ggsave("figures/prediction_true_plot.png", p_tp, width = 6, height = 4)

# -------------------------------------------------------------
# Load and reshape simulation results
load_sim_data <- function(path, model_beta_path) {
  sim_df <- read_dta(path)
  setDT(sim_df)
  sim_df[, simulation := .I]
  sim_df <- melt(sim_df, id.vars = "simulation", variable.name = "var", value.name = "sim_coeff")
  sim_df[, lag := as.integer(gsub("[^0-9]", "", var))]

  # Load model estimates
  m_linear_beta <- read_dta(model_beta_path)
  sim_df <- merge(sim_df, m_linear_beta, by = "lag", all.x = TRUE)

  sim_df
}

# Example use (you can repeat this per simulation set)
# sim_df <- load_sim_data("output/randomization/total_randomization_output/total_1000_simulation_output.dta", 
#                         "output/m_linear_beta.dta")

# -------------------------------------------------------------
# p-value curve plot (simplified)
# Use sim_df from the previous function
plot_pval_curve <- function(sim_df, out_path) {
  sim_df[, p_indicator := ifelse(B_s < abs(sim_coeff), 1, 0)]
  pval_df <- sim_df[, .(p_value = mean(p_indicator, na.rm = TRUE)), by = lag]
  pval_df[lag == 0, p_value := 0]

  ggplot(pval_df, aes(x = lag, y = p_value)) +
    geom_line() +
    geom_hline(yintercept = c(0.01, 0.05, 0.1), linetype = "dashed", color = "grey") +
    labs(
      x = "Months Since Tropical Cyclone",
      y = "p-value"
    ) +
    theme_minimal()

  ggsave(out_path, width = 6, height = 4)
}

# plot_pval_curve(sim_df, "figures/appendix/total_rand_pvalue.png")

# -------------------------------------------------------------
# Rainfall data
# Clean and merge rainfall datasets

rainfall <- read_csv("data/rainfall_idw_state_storm.csv")
storm_meta <- read_delim("data/storm_list.txt", delim = ",")

# Clean strings
storm_meta <- storm_meta %>%
  rename(storm_name = v1, year = v9, ep_name = v21) %>%
  mutate(across(c(ep_name, storm_name), ~ str_replace_all(., " ", "")))

rainfall <- rainfall %>%
  mutate(
    storm_name1 = str_sub(storm_name, 1, nchar(storm_name) - 5),
    storm_name1 = str_replace_all(storm_name1, " ", "_"),
    ep_name = str_replace_all(toupper(storm_name1), "TROPICAL_DEPRESSION_", ""),
    filename = storm_name,
    year = as.integer(str_sub(filename, -4))
  )

rainfall_merged <- left_join(rainfall, storm_meta, by = c("ep_name", "year"))

# -------------------------------------------------------------
# Regress rainfall ~ windspeed
damage_data <- read_csv("data/panel_by_storm__NA_USA_density_8_yr_1930_2018.csv")

# Clean and merge with rainfall
damage_data <- damage_data %>%
  mutate(across(c(storm_name, adm_name), ~ str_replace_all(., " ", "")))

rain_df <- inner_join(damage_data, rainfall_merged, by = c("storm_name", "adm_name", "year")) %>%
  mutate(average_rain = replace_na(average_rain, 0))

# Regress
rain_model <- lm(average_rain ~ maxs, data = rain_df)
summary(rain_model)

# Plot
ggplot(rain_df, aes(x = maxs, y = average_rain)) +
  geom_point(alpha = 0.2, color = "purple") +
  geom_smooth(method = "lm") +
  labs(
    x = "Max Wind Speed (m/s)",
    y = "Average Rainfall (inches)",
    caption = paste0("R² = ", round(summary(rain_model)$r.squared, 3))
  ) +
  theme_minimal()

ggsave("figures/scatter_rain_wind_full.png", width = 6, height = 4)

# -------------------------------------------------------------
# Final combine (Patchwork)
final_fig <- (p_fl | p_nj) / p_tp
ggsave("figures/figure1.png", final_fig, width = 10, height = 8)

