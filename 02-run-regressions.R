# ----------------------- 1. Setup -----------------------
library(fixest)
library(dplyr)
library(data.table)
library(haven)
library(purrr)
library(broom)
library(stringr)
library(tidyr)

# Load dataset
df <- read_dta("output/full_data_for_regression.dta")
setDT(df)

# Create state-level average of maxs and quartiles
df[, maxs_mean_st := mean(maxs, na.rm = TRUE), by = stfips]
df[, pct_4 := cut(maxs_mean_st, 
                  breaks = quantile(maxs_mean_st, probs = seq(0, 1, 0.25), na.rm = TRUE), 
                  labels = FALSE, include.lowest = TRUE)]
df[, pct_2 := pct_4]
df[pct_2 >= 2, pct_2 := 2]

# --------------------- 2. Define Outcomes and Models ---------------------
outcome_list <- c(
  "tdths_ttpop", "dths_0144_pop", "dths_4564_pop", "dths_6599_pop",
  "dths_0000_apop", "dths_0144_apop", "dths_4564_apop", "dths_6599_apop",
  "dths_black_bpop", "dths_white_wpop", "dths_0000_pop", "dths_black_pop", "dths_white_pop"
)

model_list <- c(
  "m_linear", "m_nonlinear", "m_cubic", "m_linear_leads",
  "m_linear_adapt_4", "m_linear_adapt_pool24",
  "m_nonlinear_adapt_4", "m_nonlinear_adapt_pool24", "m_cubic_adapt_pool24"
)

# --------------------- 3. Regression Runner ---------------------
run_model <- function(df, outcome, model) {
  fe_vars <- c("modate")  # Fixed effects
  
  # Base predictors
  X_base <- c(
    grep("^_iy_", names(df), value = TRUE),
    grep("^_imt", names(df), value = TRUE),
    grep("^_mi_", names(df), value = TRUE),
    grep("^_iT_", names(df), value = TRUE),
    grep("^_LL\\d+_maxs$", names(df), value = TRUE)
  )
  
  # Additional predictors based on model
  X_extra <- character(0)
  if (model == "m_nonlinear") {
    X_extra <- grep("^_LL\\d+_maxs_2$", names(df), value = TRUE)
  } else if (model == "m_cubic") {
    X_extra <- c(
      grep("^_LL\\d+_maxs_2$", names(df), value = TRUE),
      grep("^_LL\\d+_maxs_3$", names(df), value = TRUE)
    )
  } else if (model == "m_linear_leads") {
    X_extra <- grep("^_FF\\d+_maxs$", names(df), value = TRUE)
  }
  
  # Combine and escape variable names
  all_vars <- c(X_base, X_extra)
  all_vars_escaped <- paste0("`", all_vars, "`")
  formula_rhs <- paste(all_vars_escaped, collapse = " + ")
  
  # Adaptation interaction models
  if (model %in% c("m_linear_adapt_4", "m_nonlinear_adapt_4")) {
    formula_rhs <- paste0("(", formula_rhs, ") * factor(pct_4)")
  } else if (model %in% c("m_linear_adapt_pool24", "m_nonlinear_adapt_pool24", "m_cubic_adapt_pool24")) {
    formula_rhs <- paste0("(", formula_rhs, ") * factor(pct_2)")
  }
  
  # Build and run the regression formula
  fml <- as.formula(paste0(outcome, " ~ ", formula_rhs, " | ", paste(fe_vars, collapse = "+")))
  
  model_fit <- feols(fml, data = df, weights = ~totalpop)
  return(model_fit)
}

# --------------------- 4. Loop Over All Models and Outcomes ---------------------
dir.create("output/beta", showWarnings = FALSE, recursive = TRUE)

for (Y in outcome_list) {
  for (m in model_list) {
    cat("Running model:", m, "on outcome:", Y, "\n")
    
    # Run regression
    model_fit <- run_model(df, Y, m)
    
    # Save full model
    saveRDS(model_fit, file = paste0("output/beta/", m, "_", Y, ".rds"))
    
    # Save coefficients as CSV
    coefs <- broom::tidy(model_fit)
    write.csv(coefs, file = paste0("output/beta/", m, "_", Y, "_coefs.csv"), row.names = FALSE)
  }
}

# --------------------- 5. Extract Beta Coefficients (Optional) ---------------------
extract_beta_coefficients <- function(model_fit, pattern = "_LL\\d+_maxs") {
  coef_df <- broom::tidy(model_fit)
  lag_rows <- coef_df %>% filter(str_detect(term, pattern))
  
  beta_matrix <- lag_rows %>%
    mutate(lag = str_extract(term, "\\d+")) %>%
    select(lag, estimate) %>%
    pivot_wider(names_from = lag, values_from = estimate, names_prefix = "beta_L") %>%
    as.data.frame()
  
  return(beta_matrix)
}


