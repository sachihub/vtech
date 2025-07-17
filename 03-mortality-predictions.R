library(dplyr)
library(readr)
library(haven)

# Define outcomes and models
outcome_list <- c("_0000_apop", "_0144_apop", "_4564_apop", "_6599_apop", "_black_bpop", "_white_wpop")
model_list_basic <- c("m_linear", "m_cubic_adapt_pool24")
model_list_full <- c("m_linear", "m_nonlinear", "m_linear_leads", "m_linear_adapt_4", 
                     "m_linear_adapt_pool24", "m_nonlinear_adapt_4", "m_nonlinear_adapt_pool24",
                     "m_cubic_adapt_pool24", "m_cubic")

#------------------ Run attribute() for all outcome-model combinations
for (n in model_list_basic) {
  for (z in outcome_list) {
    attribute(Y = z, m = n, fixedpop = "4")
  }
}

#------------------ Run attribute() for tdths_ttpop for all models
for (n in model_list_full) {
  attribute(Y = "tdths_ttpop", m = n, fixedpop = "4")
}

#------------------ Run attribute() for quartiles of m_cubic_adapt_pool24
for (j in 1:4) {
  attribute(Y = "tdths_ttpop", m = "m_cubic_adapt_pool24", fixedpop = as.character(j))
}

#------------------ Combine storm-state-month data
files <- list.files("output", pattern = "^mortality_predict_m_cubic_adapt_pool24_tdths_ttpop_\\d+\\.dta$", full.names = TRUE)
mortality_data <- lapply(files, read_dta) %>% bind_rows()
write_dta(mortality_data, "output/mortality_predict_m_cubic_adapt_pool24_tdths_ttpop.dta")

#------------------ Collapse and merge prediction data (by date and by state)
outcome_list <- c("tdths_ttpop", "_0000_apop", "_0144_apop", "_4564_apop", "_6599_apop", "_black_bpop", "_white_wpop")

for (Y in outcome_list) {

  model_list <- if (Y == "tdths_ttpop") model_list_full else model_list_basic
  mortmodel <- character(0)

  # Collapse each model's date-level prediction
  for (m in model_list) {
    df <- read_dta(paste0("output/mort_TC_date_total_", m, "_", Y, ".dta"))
    collapsed <- df %>%
      group_by(modate) %>%
      summarise(mort = sum(mort, na.rm = TRUE),
                mort_sd = sd(mort, na.rm = TRUE),
                mort_n = n())
    write_dta(collapsed, paste0("output/mort_TC_date_total_", m, "_", Y, ".dta"))
  }

  # Merge all models into one date-level dataset
  df_main <- read_dta(paste0("output/mort_TC_date_total_m_linear_", Y, ".dta"))
  for (m in model_list) {
    df_temp <- read_dta(paste0("output/mort_TC_date_total_", m, "_", Y, ".dta")) %>%
      rename(!!paste0("mort_", m) := mort)
    df_main <- left_join(df_main, df_temp, by = "modate")
    mortmodel <- c(mortmodel, paste0("mort_", m))
  }

  write_dta(df_main, paste0("output/mort_TC_date_allmodels_", Y, ".dta"))

  # Merge with observed death data
  full_data <- read_dta("output/full_data_for_regression.dta") %>%
    filter(year >= 1930) %>%
    group_by(year, month) %>%
    summarise(dths_tot = sum(dths_tot, na.rm = TRUE)) %>%
    mutate(modate = (year - 1900) * 12 + month)

  merged <- left_join(full_data, df_main, by = "modate") %>%
    filter(year >= 1950)

  # Compute mortality proportions
  for (m in model_list) {
    mort_col <- paste0("mort_", m)
    merged[[mort_col]] <- na_if(merged[[mort_col]], 0)
    merged[[paste0("p_total_", m)]] <- merged[[mort_col]] / merged$dths_tot
  }

  # Summary statistics
  mort_summary <- merged %>%
    summarise(across(all_of(mortmodel), list(sum = sum, mean = mean), na.rm = TRUE))

  # Save merged file
  write_dta(merged, paste0("output/mort_TC_date_allmodels_", Y, ".dta"))

  #------------------ Collapse and merge prediction data by state
  # Total observed deaths by state
  full_data_state <- read_dta("output/full_data_for_regression.dta") %>%
    group_by(adm_name) %>%
    summarise(dths_tot = sum(dths_tot, na.rm = TRUE))

  for (m in model_list) {
    state_model <- read_dta(paste0("output/mort_TC_state_total_", m, "_", Y, ".dta")) %>%
      rename(!!paste0("mort_", m) := mort)
    full_data_state <- left_join(full_data_state, state_model, by = "adm_name")
  }

  write_dta(full_data_state, paste0("output/mort_TC_state_total_allmodels_", Y, ".dta"))

  # Summarize by state
  state_summary <- full_data_state %>%
    summarise(across(starts_with("mort_"), sum, na.rm = TRUE))

  write_dta(full_data_state, paste0("output/mort_TC_state_allmodels_", Y, ".dta"))
}

