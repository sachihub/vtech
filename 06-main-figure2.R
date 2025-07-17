# 04-main-figure2.R
# Conversion of Stata script for Figure 2 to R
# Requires: tidyverse, fixest, haven, patchwork, broom, sandwich

library(haven)
library(dplyr)
library(tidyr)
library(stringr)
library(fixest)
library(ggplot2)
library(patchwork)
library(broom)
library(purrr)
library(modelr)

message("---------- STARTING FIGURE 2 ----------")

#---------------------------------------------
# Helper function to reshape simulation output
process_simulation <- function(file, varname) {
  read_dta(file) %>%
    rename(sim_coeff_L240D_maxs = sim_coeff_L240_maxs) %>%
    select(matches("sim_coeff_L.*D_maxs")) %>%
    filter(!is.na(sim_coeff_L1D_maxs)) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(everything(), names_to = "term", values_to = varname) %>%
    mutate(lag = as.numeric(str_extract(term, "\\d+"))) %>%
    select(lag, all_of(varname))
}

# Load and process each simulation dataset
sim_total  <- process_simulation("output/randomization/total_randomization_output/total_1000_simulation_output.dta", "total_rand")
sim_within <- process_simulation("output/randomization/whithin_panel_randomization_output/whithin_panel_1000_simulation_output.dta", "within_panel")
sim_month  <- process_simulation("output/randomization/whithin_monthyear_randomization_output/whithin_monthyear_1000_simulation_output.dta", "within_month")
sim_cross  <- process_simulation("output/randomization/cross_panel_randomization_output/cross_panel_1000_simulation_output.dta", "cross_panel")

# Merge all
randomization_means <- list(sim_total, sim_within, sim_month, sim_cross) %>%
  reduce(full_join, by = "lag")

write_dta(randomization_means, "data/randomization_mean.dta")

#-----------------------------------------------------
# Load main regression data
df <- read_dta("output/full_data_for_regression.dta") %>%
  mutate(lag = row_number() - n() + 240,
         lag = ifelse(lag < -72, NA, lag)) %>%
  filter(!is.na(lag))

# Fit FE model (modate absorbed)
fe_model <- feols(tdths_ttpop ~ . - modate | modate, data = df)

# Get lag terms
lags_ff <- paste0("_FF", 1:72, "_maxs")
lags_ll <- paste0("_LL", 0:240, "_maxs")

# Extract cumulative effect of lags
coef_all <- coef(fe_model)
vcov_all <- vcov(fe_model, cluster = "modate")

get_cumulative <- function(lags, coef_vec, vcov_mat) {
  idx <- names(coef_vec) %in% lags
  beta <- sum(coef_vec[idx])
  se <- sqrt(sum(vcov_mat[idx, idx]))
  ci <- 1.96 * se
  return(c(beta, beta - ci, beta + ci))
}

# Calculate B_leads and CI for each lag
calc_betas <- function(lag_seq, prefix) {
  map_dfr(lag_seq, function(i) {
    lag_name <- paste0(prefix, i, "_maxs")
    if (lag_name %in% names(coef_all)) {
      val <- coef_all[lag_name]
      se <- sqrt(vcov_all[lag_name, lag_name])
      ci <- 1.96 * se
      tibble(lag = ifelse(prefix == "_FF", -i, i),
             b = val,
             ci1 = val - ci,
             ci2 = val + ci)
    } else {
      tibble(lag = ifelse(prefix == "_FF", -i, i), b = NA, ci1 = NA, ci2 = NA)
    }
  })
}

b_leads <- bind_rows(
  calc_betas(1:72, "_FF"),
  calc_betas(0:240, "_LL")
)

# Merge with randomization
b_leads <- left_join(b_leads, randomization_means, by = "lag")

#-----------------------------------------------------
# Plot main figure (Figure 2A)
plot_main <- ggplot(b_leads, aes(x = lag)) +
  geom_ribbon(aes(ymin = ci1, ymax = ci2), fill = "forestgreen", alpha = 0.3) +
  geom_line(aes(y = b), color = "forestgreen") +
  geom_line(aes(y = cross_panel), color = "orange", linetype = "dashed") +
  geom_line(aes(y = total_rand), color = "orange", linetype = "dotted") +
  geom_line(aes(y = within_panel), color = "orange", linetype = "dotdash") +
  geom_line(aes(y = within_month), color = "orange", linetype = "longdash") +
  geom_vline(xintercept = 172, linetype = "dashed") +
  geom_vline(xintercept = 0) +
  labs(x = "Months Since Tropical Cyclone",
       y = "Cumulative All Cause Mortality\n(per 100,000 per m/s)") +
  theme_minimal()

ggsave("figures/figure2_panelA.pdf", plot_main, width = 9, height = 6)

#-----------------------------------------------------
# Regress b on lag and lag^2
b_leads <- b_leads %>% mutate(lag2 = lag^2)
reg_quad <- lm(b ~ lag + lag2, data = filter(b_leads, lag > 0 & lag < 172))

#-----------------------------------------------------
# Loop over outcomes and extract cumulative effects
outcomes <- c("tdths_ttpop", "dths_cvd_pop", "dths_neo_pop", "dths_rpd_pop", 
              "dths_ifd_pop", "dths_mva_pop", "dths_other_pop", "dths_0000_pop", 
              "dths_0144_pop", "dths_4564_pop", "dths_6599_pop", "dths_0000_apop", 
              "dths_0144_apop", "dths_4564_apop", "dths_6599_apop", 
              "dths_black_pop", "dths_white_pop", "dths_black_bpop", "dths_white_wpop")

for (out in outcomes) {
  fml <- as.formula(paste(out, "~ . - modate | modate"))
  mod <- feols(fml, data = df)
  coefs <- coef(mod)
  vcovs <- vcov(mod, cluster = "modate")

  # Extract cumulative _LL0_maxs to _LL240_maxs
  terms <- paste0("_LL", 0:240, "_maxs")
  idx <- names(coefs) %in% terms
  est <- sum(coefs[idx], na.rm = TRUE)
  se <- sqrt(sum(vcovs[idx, idx], na.rm = TRUE))
  ci <- 1.96 * se
  message(paste(out, ": Estimate =", round(est, 4), "CI [", round(est - ci, 4), ",", round(est + ci, 4), "]"))
}

#-----------------------------------------------------
# Adaptation models
df <- df %>%
  group_by(stfips) %>%
  mutate(mean_exposure = mean(maxs, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    pct_4 = ntile(mean_exposure, 4),
    pct_2 = if_else(ntile(mean_exposure, 4) >= 2, 2, 1)
  )

# Interaction models
adapt4 <- feols(tdths_ttpop ~ c(_LL*_maxs) * factor(pct_4) | modate, data = df)
adapt2 <- feols(tdths_ttpop ~ c(_LL*_maxs) * factor(pct_2) | modate, data = df)

# Export models
saveRDS(adapt4, file = "output/beta/m_adapt_linear_groups4.rds")
saveRDS(adapt2, file = "output/beta/m_adapt_linear_groups2.rds")

#-----------------------------------------------------
# Stacked plots (simplified example)
# Replace below with actual CI variables per outcome group
gg_age <- ggplot(df, aes(x = lag)) +
  geom_line(aes(y = B_dths_0000_apop), color = "blue") +
  geom_line(aes(y = B_dths_0144_apop), color = "red") +
  geom_line(aes(y = B_dths_4564_apop), color = "purple") +
  geom_line(aes(y = B_dths_6599_apop), color = "green") +
  labs(x = "Months Since Tropical Cyclone",
       y = "Cumulative Mortality by Age",
       title = "Age Specific Effects") +
  theme_minimal()

ggsave("figures/appendix/est_by_age.pdf", gg_age, width = 7, height = 5)

#-----------------------------------------------------
# Final Combined Plot (Figure 2)
final_figure <- plot_main / gg_age + plot_layout(ncol = 1)
ggsave("figures/figure2.pdf", final_figure, width = 10, height = 10)
