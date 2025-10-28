# # ADR Risk Simulation Script with Sensitivity Analysis

# Author: Alfie Lloyd
# Date: 23/10/2025

# Purpose:
# Simulates the number of individual dogs that may have been treated with 
# bedinvetmab to estimate realistic adverse drug event risks for the drug with
# uncertainty quantification through sensitivity analysis.

# 1. Setup ---------------------------------------------------------------------

# Load required packages
required_packages <- c(
  "tidyverse", "readxl", "here", "survival", "truncnorm", "forcats", "viridis",
  "tictoc", "janitor", "scales", "progressr", "gt", "jsonlite", "lubridate",
  "gtExtras", "htmltools", "future", "furrr","lhs"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[!packages %in% rownames(installed.packages())]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
}

install_if_missing(required_packages)
invisible(lapply(required_packages, library, character.only = TRUE))


# HYPERPARAMETER CONFIGURATION 

# Sensitivity analysis settings
N_SIMS <- 1000                  
N_WORKERS <- 8                   
CUTOFF_DATE <- ymd("2025-06-30")  

# Latin Hypercube
lhs_matrix <- randomLHS(N_SIMS, 5)

# Hyperparameter ranges (uniform distributions unless noted)
CORRELATION_MIN <- -0.3          # between age and weight. -.2 across full dataset in Montoya but ~0 if restricted to dogs >8yo
CORRELATION_MAX <- -0.0          

DROPOUT_RATE_MIN <- 0.005        # Min monthly discontinuation rate
DROPOUT_RATE_MAX <- 0.03         # Max monthly discontinuation rate

DOSE_BUFFER_MIN <- 0.05          # Min dose wastage/stock (5%)
DOSE_BUFFER_MAX <- 0.30          # Max dose wastage/stock (30%)

TOTAL_DOSES_MIN <- 3e7           # Min total doses (30M)
TOTAL_DOSES_MAX <- 3.1e7         # Max total doses (31M)

# Underreporting correction (Hazell & Shakir 2006)
# Use beta distribution + min/max to get heavy left tail (less underreporting) 
# with weight of the distribution ~ 80% corresponding to Hazell
UNDERREPORTING_MIN <- 1          # No under-reporting
UNDERREPORTING_MAX <- 20         # 95% under-reporting

# Reporting rate temporal drift (since Monteiro 2025)
REPORTING_ADJ_MIN <- 0.7         
REPORTING_ADJ_MAX <- 1.3         

# Configure parallel processing
plan(multisession, workers = N_WORKERS)
cat(sprintf("✓ Parallel processing enabled with %d workers\n", N_WORKERS))

# Load adverse events data
ae_df <- read_rds(here("data", "processed", "matched_data.rds")) %>%
  filter(drug == "bedinvetmab",
         date <= CUTOFF_DATE)

# Comparison to Monteiro
monteiro <- ae_df %>%
  filter(between(date, ymd("2021-02-01"),ymd("2024-06-30")))

n_distinct(monteiro$dog_id)

# 55966 individual PTs in 16,175 ADR reports compared to   17,162

# 2. Country-specific reporting rates ------------------------------------------

reporting_rates <- tribble(
  ~country_group, ~ae_per_10k_doses,
  "United States", 18.18,
  "United Kingdom", 14.61,
  "Germany", 6.06,
  "Spain", 5.41,
  "France", 3.22,
  "Italy", 3.22,
  "Canada", 18.72,
  "Australia", 13.86,
  "EU_combined", mean(c(3.22, 5.41, 3.22)),
  "Rest_of_World", mean(c(3.22, 5.41, 3.22))
)

eu_countries <- c(
  "Austria", "Belgium", "Bulgaria", "Croatia", "Czech Republic",
  "Denmark", "Estonia", "Finland", "France", "Germany", "Greece",
  "Hungary", "Iceland", "Ireland", "Italy", "Lithuania", "Luxembourg",
  "Netherlands", "Norway", "Poland", "Portugal", "Slovakia", "Slovenia",
  "Spain", "Sweden", "Switzerland", "Andorra"
)

# 3. Life expectancy tables ----------------------------------------------------

# Load from Montoya Supplementary material Table 8

life_table_by_size <- read_rds(here("data", "life_table_by_size"))

# 4. Helper Functions ----------------------------------------------------------

assign_size_group <- function(weights) {
  case_when(
    weights < 5.5 ~ "Toy",
    weights < 11 ~ "Small",
    weights < 26 ~ "Medium",
    weights < 45 ~ "Large",
    TRUE ~ "Giant"
  )
}

survival_months_from_size_and_age_vectorized <- function(ages, size_groups) {
  n <- length(ages)
  survival_months <- numeric(n)
  
  dogs_df <- tibble(
    idx = 1:n,
    age = ages,
    size_group = size_groups
  )
  
  for (size in unique(size_groups)) {
    size_rows <- life_table_by_size %>% filter(size_group == size)
    
    for (i in 1:nrow(size_rows)) {
      mask <- (dogs_df$size_group == size) & 
        (dogs_df$age >= size_rows$age_lower[i]) & 
        (dogs_df$age < size_rows$age_upper[i])
      
      n_in_group <- sum(mask)
      
      if (n_in_group > 0) {
        survival_months[mask] <- rtruncnorm(
          n_in_group,
          a = size_rows$ci_lower[i] * 12,
          b = size_rows$ci_upper[i] * 12,
          mean = size_rows$e_x_months[i],
          sd = size_rows$sd_months[i]
        )
      }
    }
  }
  survival_months <- pmin(survival_months, 216)
  return(pmax(survival_months, 1))
}

calculate_K <- function(current_cumulative, months_available) {
  t_saturation <- 72
  r_assumed <- 0.1
  t0_assumed <- t_saturation / 2
  expected_prop <- 1 / (1 + exp(-r_assumed * (months_available - t0_assumed)))
  K_estimated <- current_cumulative / expected_prop
  return(K_estimated)
}

fit_logistic_fixed_K <- function(data, K_fixed) {
  max_time <- max(data$months_since_launch)
  tryCatch({
    nls(
      cumulative_doses ~ K_fixed / (1 + exp(-r * (months_since_launch - t0))),
      data = data,
      start = list(r = 0.1, t0 = max_time * 0.5),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    )
  }, error = function(e) NULL)
}

# 5. Create Output Directories -------------------------------------------------

dir.create(here("output"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "sensitivity"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)

# 6. Define Simulation Function ------------------------------------------------


run_single_simulation <- function(sim_id, ae_df, reporting_rates, eu_countries, 
                                  life_table_by_size, cutoff_date) {
  
  # Sample hyperparameters from configured ranges 
  correlation <- qunif(lhs_matrix[sim_id, 1], CORRELATION_MIN, CORRELATION_MAX)
  dropout_rate <- qunif(lhs_matrix[sim_id, 2], DROPOUT_RATE_MIN, DROPOUT_RATE_MAX)
  dose_buffer_pct <- qunif(lhs_matrix[sim_id, 3], DOSE_BUFFER_MIN, DOSE_BUFFER_MAX)
  beta_val <- qbeta(lhs_matrix[sim_id, 4], shape1 = 2.5, shape2 = 7) 
  underreporting_mult <- UNDERREPORTING_MIN + beta_val * UNDERREPORTING_MAX
  known_global_doses_sold <- round(qunif(lhs_matrix[sim_id, 5], TOTAL_DOSES_MIN, TOTAL_DOSES_MAX))
  target_doses_consumed <- known_global_doses_sold * (1 - dose_buffer_pct) 
  
  # Country-specific reporting rate adjustments
  reporting_rates_adjusted <- reporting_rates %>%
    mutate(
      reporting_rate_adj = runif(n(), min = REPORTING_ADJ_MIN, max = REPORTING_ADJ_MAX),
      ae_per_10k_doses = ae_per_10k_doses * reporting_rate_adj
      ) 
  
# Market share estimation 
  analysis_data <- ae_df %>%
    group_by(country) %>%
    summarise(n_aes = n(), .groups = "drop") %>%
    mutate(
      country_group = case_when(
        country == "United States" ~ "United States",
        country == "United Kingdom" ~ "United Kingdom",
        country == "Canada" ~ "Canada",
        country == "France" ~ "France",
        country == "Germany" ~ "Germany", 
        country == "Italy" ~ "Italy", 
        country == "Spain" ~ "Spain",
        country == "Australia" ~ "Australia",
        country %in% eu_countries ~ "EU_combined",
        TRUE ~ "Rest_of_World" 
      ) 
    ) %>% 
    group_by(country_group) %>%
    summarise(total_aes = sum(n_aes), .groups = "drop") %>%
    left_join(reporting_rates_adjusted, by = "country_group") %>% 
    mutate( 
      license_date = case_when(
        country_group == "United States" ~ as_date("2023-05-05"), 
        country_group == "United Kingdom" ~ as_date("2020-11-10"),
        country_group == "Canada" ~ as_date("2023-03-08"), 
        country_group == "Australia" ~ as_date("2022-09-01"), 
        country_group %in% c("France", "Germany", "Italy", "Spain", "EU_combined") ~ as_date("2020-11-10"), 
        country_group == "Rest_of_World" ~ as_date("2022-01-01") 
      ), 
      estimated_doses = total_aes / (ae_per_10k_doses / 10000),
      market_share = estimated_doses / sum(estimated_doses) ) 
  
  # Logistic growth modeling # 
empirical_uptake <- ae_df %>% 
  mutate( 
    year_month = floor_date(date, "month"),
    country_group = case_when( 
      country == "United States" ~ "United States", 
      country == "United Kingdom" ~ "United Kingdom",
      country == "Canada" ~ "Canada", 
      country == "France" ~ "France",
      country == "Germany" ~ "Germany",
      country == "Italy" ~ "Italy",
      country == "Spain" ~ "Spain", 
      country == "Australia" ~ "Australia",
      country %in% eu_countries ~ "EU_combined", 
      TRUE ~ "Rest_of_World" 
    ) 
  ) %>%
  count(country_group, year_month, name = "n_aes") %>%
  left_join(analysis_data %>%
              dplyr::select(country_group, license_date), by = "country_group") %>% 
  filter(year_month >= license_date) %>%
  left_join(reporting_rates_adjusted %>% 
              dplyr::select(country_group, ae_per_10k_doses), by = "country_group") %>% 
  mutate( 
    months_since_launch = as.numeric(interval(license_date, year_month) / months(1)),
    estimated_doses = n_aes / (ae_per_10k_doses / 10000)
  ) %>% 
  group_by(country_group) %>%
  arrange(year_month) %>%
  mutate(cumulative_doses = cumsum(estimated_doses)) %>% 
  ungroup() 

logistic_params <- empirical_uptake %>%
  group_by(country_group) %>%
  summarise(
    months_available = max(months_since_launch), 
    current = max(cumulative_doses), 
    K = calculate_K(current, months_available),
    model = list(fit_logistic_fixed_K(pick(everything()), K)),  # ← CHANGE THIS 
    r = ifelse(!is.null(model[[1]]), coef(model[[1]])["r"], 0.1), 
    t0 = ifelse(!is.null(model[[1]]), coef(model[[1]])["t0"], months_available * 0.5), 
    pct_capacity = 100 * current / K, .groups = "drop" ) 

analysis_data <- analysis_data %>% left_join(logistic_params, by = "country_group") 

# Define parameterized functions 

sample_dog_age_and_weight_param <- function(n, corr, a_sd, w_sd) {
  
  mu <- c(0, 0) 
  sigma <- matrix(c(1, corr, corr, 1), nrow = 2)
  z <- MASS::mvrnorm(n, mu = mu, Sigma = sigma)
  
  u1 <- pnorm(z[, 1]) 
  u2 <- pnorm(z[, 2])
  
  age_at_start <- qbeta(u1, shape1 = 4.5, shape2 = 2.7) * 17 + 1 
  weight_kg <- pmin(qlnorm(u2, meanlog = log(26.06), sdlog = 0.57), 100)

  tibble(age_at_start, weight_kg) 
} 

treatment_persistence_param <- function(survival_months, lambda) { 
  
  n <- length(survival_months)
  treatment_months <- rexp(n, rate = lambda)
  pmin(ceiling(treatment_months), ceiling(survival_months))
} 

# Run simulation 

dogs_list <- list() 
cumulative_doses <- 0 
batch_num <- 0 

while (cumulative_doses < target_doses_consumed) { 
  batch_num <- batch_num + 1
  
  remaining_pct <- cumulative_doses / target_doses_consumed 
  batch_size <- case_when( 
    remaining_pct < 0.90 ~ 50000,
    remaining_pct < 0.96 ~ 10000,
    TRUE ~ 2000 
  ) 
  
  dog_demographics <- sample_dog_age_and_weight_param(
    batch_size, correlation, age_sd, weight_sd 
  )
  
  batch_dogs <- tibble(
    dog_id = (batch_num - 1) * 50000 + 1:batch_size,
    country_group = sample( 
      analysis_data$country_group,
      size = batch_size, 
      replace = TRUE,
      prob = analysis_data$market_share
    ) 
  ) %>% 
    left_join(logistic_params %>%
                dplyr::select(country_group, K, r, t0), 
              by = "country_group") %>% 
    left_join(analysis_data %>% 
                dplyr::select(country_group, license_date),
              by = "country_group") %>%
    bind_cols(dog_demographics) %>%
    mutate(
      total_months = as.numeric(interval(license_date, cutoff_date) / months(1)),
      N_total = K / (1 + exp(-r * (total_months - t0))),
      u = runif(n()) * N_total, 
      months_to_entry = t0 - (1 / r) * log(pmax(K / u - 1, 0.001)), 
      months_to_entry = pmin(pmax(months_to_entry, 0), total_months), 
      entry_date = license_date + months(round(months_to_entry))
    ) 
  
  batch_dogs$size_group <- assign_size_group(batch_dogs$weight_kg) 
  
  batch_dogs$survival_months <- survival_months_from_size_and_age_vectorized(
    batch_dogs$age_at_start, 
    batch_dogs$size_group
  ) 
  
  batch_dogs$treatment_months <- treatment_persistence_param( 
    batch_dogs$survival_months,
    dropout_rate
  )
  
  batch_dogs <- batch_dogs %>%
    mutate( 
      months_to_cutoff = as.numeric(interval(entry_date, cutoff_date) / months(1)),
      doses_by_cutoff = pmin(treatment_months, pmax(0, months_to_cutoff)) %>% ceiling()
    ) %>% 
    filter(entry_date <= cutoff_date, doses_by_cutoff > 0)
  
  dogs_list[[batch_num]] <- batch_dogs 
  cumulative_doses <- cumulative_doses + sum(batch_dogs$doses_by_cutoff)
  
  if (batch_num > 300) { 
    warning("Exceeded maximum iterations")
    break 
  } 
}

dogs_all <- bind_rows(dogs_list)

# Trim excess doses 

if (cumulative_doses > target_doses_consumed) {
  dogs_sorted <- dogs_all %>%
    arrange(desc(doses_by_cutoff)) 
  dogs_to_remove <- 0 
  running_total <- cumulative_doses 
  
  for (i in 1:nrow(dogs_sorted)) {
    if (running_total - dogs_sorted$doses_by_cutoff[i] >= target_doses_consumed) { 
      running_total <- running_total - dogs_sorted$doses_by_cutoff[i] 
      dogs_to_remove <- dogs_to_remove + 1
    
      } else {
        
      break 
      }
  }
  
  dogs_final <- dogs_all %>%
    arrange(desc(doses_by_cutoff)) %>%
    slice(-(1:dogs_to_remove)) 
} else { 
    dogs_final <- dogs_all 
} 

# Calculate ADR risks (as percentages) 

exposure_summary <- dogs_final %>%
  count(country_group, name = "n_exposed") 

adr_summary <- ae_df %>%
  mutate(country_group = case_when(
    country == "United States" ~ "United States", 
    country == "United Kingdom" ~ "United Kingdom",
    country == "Canada" ~ "Canada",
    country == "France" ~ "France",
    country == "Germany" ~ "Germany", 
    country == "Italy" ~ "Italy",
    country == "Spain" ~ "Spain",
    country == "Australia" ~ "Australia",
    country %in% eu_countries ~ "EU_combined", TRUE ~ "Rest_of_World"
  )) %>% 
  count(country_group, pt, name = "n_observed") %>%
  mutate(n_corrected = n_observed * underreporting_mult)

adr_risks <- adr_summary %>%
  left_join(exposure_summary, by = "country_group") %>% 
  mutate( 
    risk_pct_observed = 100 * n_observed / n_exposed,
    risk_pct_corrected = 100 * n_corrected / n_exposed 
  ) 

# Return results 

result <- adr_risks %>% 
  mutate(
    sim_id = sim_id, 
    correlation = correlation, 
    dropout_rate = dropout_rate,
    dose_buffer_pct = dose_buffer_pct,
    total_doses = known_global_doses_sold,
    underreporting_mult = underreporting_mult,
    total_dogs = nrow(dogs_final),
    total_doses_consumed = sum(dogs_final$doses_by_cutoff) 
  ) 

return(result) 
}


# 7. Run Simulations -----------------------------------------------------------

cat("STARTING PARALLEL SENSITIVITY ANALYSIS\n")


start_time <- Sys.time()
tic()
# Run simulations in parallel
results_list <- with_progress({
  p <- progressor(steps = N_SIMS)
  
  future_map(1:N_SIMS, function(sim_id) {
    start_time <- Sys.time()
    res <- run_single_simulation(
      sim_id = sim_id,
      ae_df = ae_df,
      reporting_rates = reporting_rates,
      eu_countries = eu_countries,
      life_table_by_size = life_table_by_size,
      cutoff_date = CUTOFF_DATE
    )
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    p(message = sprintf("Sim %d (%.1fs)", sim_id, elapsed))
    res
  }, .options = furrr_options(seed = TRUE))
})
toc()
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))

print(elapsed)
cat("SENSITIVITY ANALYSIS COMPLETE\n")
cat(sprintf("Total: %.1f min | Average: %.2f sec per sim\n", 
            elapsed, elapsed*60 / N_SIMS))
# 150 mins, average 9s per simulation

# 8. Aggregate and Save Results ------------------------------------------------

all_results <- bind_rows(results_list)

saveRDS(all_results, here("output", "simulation", "full_sim_dataset.rds"))
# write_csv(all_results, here("output", "simulation", "full_sim_dataset.csv"))

all_results <- readRDS(here("output","simulation","full_sim_dataset.rds"))

# 9. Summary Statistics --------------------------------------------------------

# ADR risk summary with 95% CIs

# Weighted quantiles to reflect different market share
wq <- function(x, w, p) {
  Hmisc::wtd.quantile(x, weights = w, probs = p, na.rm = TRUE)
}

final_summary <- all_results %>%
  group_by(pt, country_group) %>%
  summarise(
    n_sims = n(),
    median_risk_observed   = median(risk_pct_observed,  na.rm = TRUE),
    lower_95_observed      = quantile(risk_pct_observed,  0.025, na.rm = TRUE),
    upper_95_observed      = quantile(risk_pct_observed,  0.975, na.rm = TRUE),
    median_risk_corrected  = median(risk_pct_corrected, na.rm = TRUE),
    lower_95_corrected     = quantile(risk_pct_corrected, 0.025, na.rm = TRUE),
    upper_95_corrected     = quantile(risk_pct_corrected, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(country_group, desc(median_risk_corrected)) %>%
  bind_rows(
    all_results %>%
      group_by(pt) %>%
      summarise(
        n_sims = n(),
        country_group = "Total",
        median_risk_corrected = wq(risk_pct_corrected, w = n_exposed, p = 0.5),
        lower_95_corrected    = wq(risk_pct_corrected, w = n_exposed, p = 0.025),
        upper_95_corrected    = wq(risk_pct_corrected, w = n_exposed, p = 0.975),
        .groups = "drop"
      )
  ) 

# 9. Display Top ADRs ---------------------------------------------------------

cat("TOP 20 ADRs (Global) - CORRECTED FOR UNDERREPORTING\n")
final_summary %>%
  filter(pt != "Lack of efficacy", country_group == "Total") %>%
  slice_max(median_risk_corrected, n = 20) %>%
  mutate(
    risk_display = sprintf("%.3f%% [%.3f%% - %.3f%%]", 
                           median_risk_corrected,
                           lower_95_corrected,
                           upper_95_corrected  
  )) %>%
  select(pt, risk_display) %>%
  print(n = 20)

# Parameter summary
param_summary <- all_results %>%
  distinct(sim_id, .keep_all = TRUE) %>%
  select(correlation:total_doses_consumed) %>%
  summarise(across(everything(), list(
    mean = mean, median = median, sd = sd, min = min, max = max
  ), .names = "{.col}_{.fn}"))

print(param_summary, n = Inf)
write_csv(param_summary, here("output", "simulation", "parameter_summary.csv"))

# 11. Plots --------------------------------------------------------------------

# Plot 1: Parameter distributions
param_dist_plot <- all_results %>%
  distinct(sim_id, .keep_all = TRUE) %>%
  select(sim_id, correlation, dropout_rate, dose_buffer_pct,
         underreporting_mult, total_dogs) %>%
  pivot_longer(-sim_id, names_to = "parameter", values_to = "value") %>%
  mutate(
    parameter = recode(parameter,
                       "correlation" = "Age-Weight Correlation",
                       "dropout_rate" = "Monthly Dropout Rate",
                       "dose_buffer_pct" = "Dose Buffer %",
                       "underreporting_mult" = "Underreporting Multiplier",
                       "total_dogs" = "Total Dogs Simulated")
  ) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, fill = "#9b4b16", alpha = 0.7, color = "#2d2d2d") +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  labs(
    title = "Distribution of Sampled Parameters Across Simulations",
    subtitle = sprintf("n = %d simulations", N_SIMS),
    x = "Parameter Value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    strip.background = element_rect(fill = "#f2e6dd", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "#555555")
  )

param_dist_plot

ggsave(here("output", "figures", "parameter_distributions.png"), 
       param_dist_plot, width = 12, height = 8, dpi = 300, bg = "#fff1e5")


# Plot 2: Top ADRs Globally

top_adrs_plot <- final_summary %>%
  filter(pt != "Lack of efficacy", country_group == "Total") %>%
  slice_max(median_risk_corrected, n = 20) %>%
  mutate(
    pt = str_to_title(pt) %>% str_replace_all("_", " "),
    pt = recode(pt,
                "Increased Blood Urea Nitrogen (Bun) Or Creatinine" = "↑ BUN/Creatinine"),
    pt = fct_reorder(pt, median_risk_corrected)
  ) %>%
  ggplot(aes(y = pt, x = median_risk_corrected)) +
  geom_errorbar(aes(xmin = lower_95_corrected, xmax = upper_95_corrected),
                 height = 0.3, color = "#555555", linewidth = 0.8) +
  geom_vline(xintercept = c(0.1, 1),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  annotate("text", x = 0.55, y = 20.5, label = "Uncommon",
           vjust = -0.5, hjust = 0.5, size = 5, fontface = "italic") +
  annotate("text", x = 2.5, y = 20.5, label = "Common",
           vjust = -0.5, hjust = 1, size = 5, fontface = "italic") +
  
  scale_x_continuous(
    limits = c(0,6.5),
    breaks = c(0, 1, 2, 3, 4, 5, 6)
  ) +
  geom_point(size = 3, color = "#9b4b16") +
  labs(
    title = "Top 20 ADRs Rates",
    subtitle = "Median with 95% confidence interval | Dotted lines show CIOMS frequency thresholds",
    x = "ADR Risk (%)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  scale_y_discrete(expand = expansion(mult = c(0.02, 0.15))) +
  theme_minimal(base_size = 16, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 16, color = "#555555"),
    axis.text.y = element_text(size = 16)
  )

top_adrs_plot

ggsave(here("output", "figures", "top_adrs_with_uncertainty.png"), 
       top_adrs_plot, width = 12, height = 8, dpi = 600, bg = "#fff1e5")

write_csv(final_summary, here("output", "simulation", "summary_sim_results.csv"))


# Plot 3. Country-specific Heatmap plot

top_pts <- final_summary %>%
  filter(pt != "Lack of efficacy", country_group == "Total") %>%
  slice_max(median_risk_corrected, n = 20) %>%
  pull(pt)

adr_top20_clean <- final_summary %>%
  filter(pt %in% top_pts) %>%
  mutate(
    pt = str_to_title(pt) |> str_replace_all("_", " "),
    country_group = case_match(
      country_group,
      "EU_combined"   ~ "EU (other)",
      "United Kingdom" ~ "UK",
      "United States" ~ "US",
      "Rest_of_World" ~ "Rest of World",
      .default = country_group
    ),
    country_group = fct_reorder(country_group, median_risk_corrected, .desc = TRUE))

adr_rate_plot <- adr_top20_clean %>%
  arrange(desc(median_risk_corrected)) %>%
  mutate(
    pt = recode(pt,
                "Increased Blood Urea Nitrogen (Bun) Or Creatinine" = "↑ BUN/Creatinine")
  ) %>%
  ggplot(aes(x = country_group, y = pt, fill = median_risk_corrected)) +
  geom_tile(color = "#d1bfa7", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f%%", median_risk_corrected)),
    size = 3,
    family = "Merriweather",
    fontface = "bold"
    )+
  scale_color_identity() +
  scale_fill_gradient(
    low = "#f2e6dd", high = "#c04a0b",
    name = "ADR rate (%)",
    labels = label_percent(accuracy = 0.01)
  ) +
  scale_x_discrete(position = "top", name = "Country or Region") +
  scale_y_discrete(
    name = "Preferred Term",
    labels = \(x) str_wrap(x, width = 25)
  ) +
  labs(
    title = "Top 20 ADR Rates by Country or Region",
    subtitle = "Median Rates calculated using simulated exposure denominators (%)"
  ) +
  theme_minimal(base_size = 13, base_family = "Merriweather") +
  theme(
    plot.background  = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid       = element_blank(),
    axis.text.x      = element_text(face = "bold", size = 11, vjust = 0.5, color = "#2d2d2d"),
    axis.text.y      = element_text(face = "bold", size = 10, color = "#2d2d2d"),
    axis.title.x     = element_text(face = "bold", size = 12, color = "#333333", margin = margin(b = 10)),
    axis.title.y     = element_text(face = "bold", size = 12, color = "#333333", margin = margin(r = 10)),
    legend.position  = "none",
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5, color = "#2d2d2d"),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "#555555"),
    plot.margin      = margin(15, 15, 15, 15)
  )

adr_rate_plot

ggsave(here("output", "figures", "adr_country_rate_heatmap.png"), 
       adr_rate_plot, width = 14, height = 8, dpi = 600, bg = "#fff1e5")

# Plot 4. Country contribution boxplot 
country_contrib <- all_results %>%
  group_by(sim_id, country_group) %>%
  slice(1) %>%
  summarise(country_n = sum(n_exposed),
            sim_doses = total_doses,
            .groups = "drop_last") %>%
  mutate(proportion = country_n / sum(country_n)) %>%
  ungroup() %>%
  mutate(
    country_group = case_match(
      country_group,
      "EU_combined"   ~ "EU (other)",
      "United Kingdom" ~ "UK",
      "United States" ~ "US",
      "Rest_of_World" ~ "Rest of World",
      .default = country_group
    )
  )

total_doses <- country_contrib %>%
  group_by(sim_id) %>%
  summarise(total_n = sum(as.numeric(country_n)), 
            total_dose = first(sim_doses),
            doses_per_dog = total_dose/total_n,
            .groups = "drop_last") %>%
  mutate(median_n = median(total_n),
         lower_n = quantile(total_n, 0.025),
         upper_n = quantile(total_n, 0.975),
         median_doses = median(doses_per_dog),
         lower_doses = quantile(doses_per_dog, 0.025),
         upper_doses = quantile(doses_per_dog, 0.975))

country_contrib_plot <- ggplot(country_contrib, 
                               aes(x = fct_reorder(country_group, proportion, median), 
                                   y = proportion)) +
  geom_boxplot(fill = "#c04a0b", alpha = 0.7, color = "#2d2d2d") +
  scale_y_continuous(labels = scales::percent, name = "Contribution to Total") +
  scale_x_discrete(name = "Country or Region") +
  coord_flip() +
  labs(
    title = "Country/Region Contributions to Total Number of Dogs Exposed"
  ) +
  theme_minimal(base_size = 13, base_family = "Merriweather") +
  theme(
    plot.background  = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid.major.x = element_line(color = "#d1bfa7", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x      = element_text(face = "bold", size = 11, color = "#2d2d2d"),
    axis.text.y      = element_text(face = "bold", size = 11, color = "#2d2d2d"),
    axis.title.x     = element_text(face = "bold", size = 12, color = "#333333", margin = margin(t = 10)),
    axis.title.y     = element_text(face = "bold", size = 12, color = "#333333", margin = margin(r = 10)),
    plot.title       = element_text(face = "bold", size = 15, hjust = 0.5, color = "#2d2d2d"),
    plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "#555555"),
    plot.margin      = margin(15, 15, 15, 15)
  )

country_contrib_plot

ggsave(here("output", "figures", "country_contribution_boxplot.png"), 
       country_contrib_plot, width = 10, height = 6, dpi = 600, bg = "#fff1e5")

sum <- final_summary %>%
  filter(country_group == "Total") %>%
  arrange(desc(median_risk_corrected))

dose_summary <- all_results %>%
  group_by(sim_id) %>%
  summarise(
    doses_per_dog = first(total_doses) / first(total_dogs),
    .groups = "drop"
  )

# Core statistics
med <- median(dose_summary$doses_per_dog)
ci <- quantile(dose_summary$doses_per_dog, c(0.025, 0.975))

library(ggtext)

plot_dose_summary <- dose_summary %>%
  ggplot(aes(x = doses_per_dog)) +
  geom_area(
    stat = "density",
    fill = "#9b4b16",
    alpha = 0.8,
    colour = NA
  ) +
  scale_x_continuous(limits = c(10,25)) +
  geom_vline(xintercept = med, colour = "#222", linewidth = 0.8) +
  geom_vline(xintercept = ci, colour = "#222", linetype = "dashed", linewidth = 0.6) +
  annotate(
    "richtext",
    x = med, y = 0.9 * ggplot_build(ggplot(dose_summary, aes(x = doses_per_dog)) +
                                      geom_density())$data[[1]]$y %>% max(),
    label = glue::glue("<b>Median:</b> {round(med, 2)}<br><span style='font-size:10pt'>95% CI: {round(ci[1], 2)} – {round(ci[2], 2)}</span>"),
    hjust = -0.05, vjust = 1, fill = NA, label.color = NA, family = "Merriweather"
  ) +
  labs(
    title = "Estimated Doses per Dog",
    subtitle = "Distribution across all simulations (median and 95% credible interval)",
    x = "Doses per dog",
    y = NULL
  ) +
  theme_minimal(base_family = "Merriweather") +
  theme(
    text = element_text(colour = "#222"),
    plot.title = element_text(size = 15, face = "bold", colour = "#4d3319"),
    plot.subtitle = element_text(size = 11, colour = "#4d3319"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#d9c7b8", linewidth = 0.3),
    plot.background = element_rect(fill = "#fff1e5", colour = NA),
    panel.background = element_rect(fill = "#fff1e5", colour = NA),
    axis.text = element_text(colour = "#4d3319"),
    axis.title.x = element_text(margin = margin(t = 8))
  )

ggsave(here("output", "figures", "simulation_dose_distribution.png"), 
       plot_dose_summary, width = 10, height = 6, dpi = 600, bg = "#fff1e5")
