library(dplyr)
library(tidyr)
library(haven)
library(data.table)
library(zoo)
library(stringr)
library(tidyverse)

#-------------------------------------------------------------------------------
# Load main dataset (1930–2015)
#-------------------------------------------------------------------------------
df <- read_dta("data/DATA_hurricane_mortality_temp_month_state_19302015.dta")

#-------------------------------------------------------------------------------
# Create lag and lead variables
#-------------------------------------------------------------------------------

# Convert to data.table for efficiency
setDT(df)
setkey(df, adm_name, modate)
show_col_types = FALSE

# Create leads: F72.maxs to F1.maxs
for (i in 72:1) {
  df[, paste0("_FF", i, "_maxs") := shift(maxs, n = -i, type = "lead"), by = adm_name]
}

# Create lags: L0.maxs to L240.maxs
for (i in 0:240) {
  df[, paste0("_LL", i, "_maxs") := shift(maxs, n = i, type = "lag"), by = adm_name]
}

# Nonlinear lags (squared, cubed)
df[, maxs_2 := maxs^2]
df[, maxs_3 := maxs^3]

for (i in 0:240) {
  df[, paste0("_LL", i, "_maxs_2") := shift(maxs_2, n = i, type = "lag"), by = adm_name]
  df[, paste0("_LL", i, "_maxs_3") := shift(maxs_3, n = i, type = "lag"), by = adm_name]
}

#-------------------------------------------------------------------------------
# Fixed effects and polynomial trends
#-------------------------------------------------------------------------------

# Center modate (mean since 1950)
modate_mean <- df[year > 1949, mean(modate, na.rm = TRUE)]
df[, modate := modate - round(modate_mean)]

# Create polynomial terms
for (i in 1:8) {
  df[, paste0("modate_", i) := modate^i]
}

# Temperature nonlinearity
df[, tavg_2 := tavg^2]
df[, tavg_3 := tavg^3]

# Generate group variables for fixed effects
df[, `_im_` := .GRP, by = .(stfips, month)]

# One-hot encoding (equivalent of xi: i.var)
df <- df %>%
  fastDummies::dummy_cols(select_columns = c("stfips", "month", "year"))

# Optional: interaction terms for fixed effects
for (i in 1:8) {
  df <- df %>%
    mutate(!!paste0("stfips_modate_", i) := as.numeric(stfips) * get(paste0("modate_", i)))
}

df <- df %>%
  mutate(imt_modate_1 = `_im_` * modate_1,
         stfips_tavg = as.numeric(stfips) * tavg,
         stfips_tavg_2 = as.numeric(stfips) * tavg_2)

# Drop redundant encodings
df <- df %>%
  select(-matches("imt_im__|iy_stfips_|iT_stfips_"))

# Save main regression dataset
write_dta(df, "output/full_data_for_regression.dta")

#-------------------------------------------------------------------------------
# Save population data for plotting
#-------------------------------------------------------------------------------

df_pop <- read_dta("data/DATA_hurricane_mortality_temp_month_state_19302015.dta")

# Annual total population
pop_annual <- df_pop %>%
  group_by(year) %>%
  summarise(totalpop = mean(totalpop, na.rm = TRUE)) %>%
  filter(year >= 1950 & year <= 2015) %>%
  arrange(year)

# Compute population change over 65 years
pop_annual <- pop_annual %>%
  mutate(t_p_change_1950_2015 = (totalpop - lag(totalpop, 65)) / lag(totalpop, 65))

# Carryforward
pop_annual <- pop_annual %>%
  arrange(desc(totalpop)) %>%
  mutate(p_change_1950_2015 = na.locf(t_p_change_1950_2015, na.rm = FALSE)) %>%
  select(year, p_change_1950_2015)

# Merge with regression data
df_plot <- read_dta("output/full_data_for_regression.dta") %>%
  select(month, year, adm_name, totalpop, pop_black, pop_white,
         pop_0000, pop_0144, pop_4564, pop_6599) %>%
  filter(year >= 1930) %>%
  mutate(modate = (year - 1900) * 12 + month,
         adm_name = str_to_upper(adm_name)) %>%
  arrange(adm_name, modate)

# Rename population variables
df_plot <- df_plot %>%
  rename(black_bpop = pop_black,
         white_wpop = pop_white,
         `_0000_apop` = pop_0000,
         `_0144_apop` = pop_0144,
         `_4564_apop` = pop_4564,
         `_6599_apop` = pop_6599) %>%
  filter(adm_name != "")

# Merge population growth data
df_plot <- left_join(df_plot, pop_annual, by = "year")

# Save output
write_dta(df_plot, "output/population_data_for_plotting.dta")

#-------------------------------------------------------------------------------
# Compute state-level maxs quartiles for attribution
#-------------------------------------------------------------------------------

df_q <- read_dta("output/full_data_for_regression.dta")

# Average maxs by state
df_q <- df_q %>%
  group_by(stfips) %>%
  mutate(maxs_mean_st = mean(maxs, na.rm = TRUE)) %>%
  ungroup()

# Quartiles
df_q <- df_q %>%
  mutate(pct = ntile(maxs_mean_st, 4)) %>%
  group_by(adm_name) %>%
  summarise(pct = mean(pct, na.rm = TRUE))

# Save
write_dta(df_q, "output/pct_maxs.dta")

#-------------------------------------------------------------------------------
# Final plotting dataset prep
#-------------------------------------------------------------------------------

df_storm <- read_csv("data/panel_by_storm__NA_USA_density_8_yr_1930_2018.csv") %>%
  select(-pddi) %>%
  mutate(modate = (year - 1900) * 12 + month,
         adm_name = str_to_upper(adm_name)) %>%
  filter(year >= 1930 & year <= 2015)

# Merge with quartiles
df_storm <- left_join(df_storm, df_q, by = "adm_name")

# Flag first modate per quartile
df_storm <- df_storm %>%
  group_by(pct) %>%
  mutate(seq = row_number(),
         firstdate_flag = if_else(modate == min(modate), 1, NA_real_)) %>%
  ungroup()

# Drop records
df_storm <- df_storm %>%
  filter(!(maxs == 0 & is.na(firstdate_flag))) %>%
  select(-firstdate_flag, -seq)

# Save final plotting dataset
write_dta(df_storm, "output/plotting_deaths.dta")

print("01-Prep job completed")

