# ===============================================================================
# COMPREHENSIVE PHARMACOVIGILANCE SIGNAL DETECTION ANALYSIS
# ===============================================================================
# 
# PURPOSE:
# This script performs comprehensive statistical signal detection for adverse
# drug reactions using multiple validated disproportionality methods.
# Starts with clean 'matched' data from the cleaning pipeline.
#
# CLINICAL CONTEXT:
# - Signal detection ≠ proof of causation
# - Statistical associations may reflect confounding, reporting bias, or chance
# - All signals require clinical evaluation and additional data sources
# - Multiple methods provide validation and confidence in findings
#
# METHODS IMPLEMENTED:
# 1. Proportional Reporting Ratio (PRR) - frequentist approach
# 2. Multi-item Gamma Poisson Shrinker (MGPS) - Bayesian shrinkage
# 3. Bayesian Confidence Propagation Neural Network (BCPNN) - Bayesian IC
# 4. Organ-system level analysis
# 5. Terms of interest analysis
#
# AUTHOR: [Your Name]
# DATE: [Current Date]
# VERSION: 2.0
# ===============================================================================

# ===============================================================================
# 1. ENVIRONMENT SETUP AND DATA LOADING
# ===============================================================================

# Clear workspace and set reproducible seed
rm(list = ls())
set.seed(42)

# Load required packages
required_packages <- c(
  "tidyverse",    # Data manipulation and visualization
  "readxl",       # Excel reading (if needed)
  "gt",           # Publication-quality tables
  "scales",       # Number formatting
  "openEBGM",     # Empirical Bayes signal detection
  "here",         # Robust file paths
  "plotly",       # Interactive plots
  "corrplot",     # Correlation visualization
  "ggplot2"       # Enhanced plotting
)

# Load packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set options
options(
  scipen = 999,
  digits = 3,
  tibble.print_max = 20
)

# Create output directories
dir.create("output", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)

# Set up analysis logging
log_file <- paste0("logs/analysis_log_", Sys.Date(), ".txt")
cat("Signal detection analysis started at:", as.character(Sys.time()), "\n", file = log_file)

message("✓ Environment setup complete")

# ===============================================================================
# 2. LOAD CLEANED DATA AND VALIDATION
# ===============================================================================

message("📁 Loading cleaned matched dataset...")

# Load the matched dataset from cleaning pipeline
matched_path <- here("data", "processed", "matched_data.rds")

if (!file.exists(matched_path)) {
  stop("❌ Matched data file not found at: ", matched_path,
       "\nPlease run the data cleaning script first.")
}

matched <- read_rds(matched_path)
message("✓ Loaded ", nrow(matched), " reaction records for ", n_distinct(matched$dog_id), " dogs")

# Validate essential columns exist
essential_columns <- c("dog_id", "drug", "pt", "llt", "organ", "hlt", "sex", "year")
missing_columns <- setdiff(essential_columns, names(matched))

if (length(missing_columns) > 0) {
  stop("❌ Essential columns missing from matched data: ", paste(missing_columns, collapse = ", "))
}

# Display data overview
cat("\n=== DATASET OVERVIEW ===\n")
cat("Total reaction records:", format(nrow(matched), big.mark = ","), "\n")
cat("Unique dogs:", format(n_distinct(matched$dog_id), big.mark = ","), "\n")
cat("Mapped to PT:", format(sum(!is.na(matched$pt)), big.mark = ","), 
    "(", round(mean(!is.na(matched$pt)) * 100, 1), "%)\n")
cat("Drugs included:", paste(sort(unique(matched$drug)), collapse = ", "), "\n")

# ===============================================================================
# 3. ANALYSIS SETUP AND CONFIGURATION
# ===============================================================================

message("⚙️ Setting up analysis parameters...")

# Define target drug for analysis (modify as needed)
target_drug <- "bedinvetmab"  # Primary drug of interest

# Alternative drug names mapping (if needed)
drug_aliases <- c(
  "bedinvetmab" = "librela",
  "librela" = "bedinvetmab"
)

# Check if target drug exists in data
available_drugs <- unique(matched$drug)
if (!target_drug %in% available_drugs) {
  # Check for aliases
  alias_match <- drug_aliases[target_drug]
  if (!is.na(alias_match) && alias_match %in% available_drugs) {
    target_drug <- alias_match
    message("ℹ️  Using alias: ", target_drug)
  } else {
    message("Available drugs: ", paste(available_drugs, collapse = ", "))
    stop("❌ Target drug '", target_drug, "' not found in dataset")
  }
}

# Set analysis parameters
min_reports <- 3  # Minimum number of reports for signal detection

# Create dog-level dataset (one row per dog per PT per drug)
dog_pt <- matched %>%
  filter(!is.na(pt)) %>%
  arrange(dog_id, drug, pt) %>%
  distinct(dog_id, drug, pt, .keep_all = TRUE)

# Calculate denominators
denominators <- dog_pt %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target) %>%
  summarise(
    n_target = sum(is_target),
    n_comparator = sum(!is_target),
    .groups = "drop"
  )

cat("\n=== ANALYSIS POPULATION ===\n")
cat("Target drug (", target_drug, "):", format(denominators$n_target, big.mark = ","), "dogs\n")
cat("Comparator drugs:", format(denominators$n_comparator, big.mark = ","), "dogs\n")
cat("Total analysis population:", format(denominators$n_target + denominators$n_comparator, big.mark = ","), "dogs\n")

message("✓ Analysis setup complete")

# ===============================================================================
# 4. PROPORTIONAL REPORTING RATIO (PRR) ANALYSIS
# ===============================================================================

message("📊 Performing PRR analysis...")

# Enhanced PRR calculation function
calculate_prr_comprehensive <- function(a, n_target, b, n_comparator) {
  # Convert to numeric to prevent overflow
  a <- as.numeric(a)
  b <- as.numeric(b)
  n_target <- as.numeric(n_target)
  n_comparator <- as.numeric(n_comparator)
  
  # Calculate 2x2 table
  c <- n_target - a
  d <- n_comparator - b
  
  # Calculate proportions
  prop_target <- a / n_target
  prop_comparator <- b / n_comparator
  
  # PRR
  prr <- prop_target / prop_comparator

  # Chi-square test (using expected frequencies to prevent overflow)
  n_total <- n_target + n_comparator
  expected_a <- (a + b) * n_target / n_total
  expected_b <- (a + b) * n_comparator / n_total
  expected_c <- (c + d) * n_target / n_total
  expected_d <- (c + d) * n_comparator / n_total
  
  chi_square <- (a - expected_a)^2 / expected_a + 
    (b - expected_b)^2 / expected_b +
    (c - expected_c)^2 / expected_c +
    (d - expected_d)^2 / expected_d
  
  
  list(
    cell_a = a, cell_b = b, cell_c = c, cell_d = d,
    n_target_total = n_target, n_comparator_total = n_comparator,
    prop_target = prop_target,
    prop_comparator = prop_comparator, 
    prr = prr,
    chi_square = chi_square
  )
}

# Perform PRR analysis
prr_results <- dog_pt %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, pt) %>%
  group_by(pt) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, denominators$n_target, 
                                             reports_comparator, denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  # Apply FDR correction
  mutate(
    disproportionality_signal = prr >= 2 & chi_square >= 4
  ) %>%
  arrange(desc(prr))

cat("\n=== PRR ANALYSIS RESULTS ===\n")
cat("Total PTs analyzed:", nrow(prr_results), "\n")
cat("disproportionality signals (PRR ≥2, χ² ≥4):", sum(prr_results$disproportionality_signal), "\n")

message("✓ PRR analysis complete")

# ===============================================================================
# 5. MULTI-ITEM GAMMA POISSON SHRINKER (MGPS) ANALYSIS
# ===============================================================================

message("🎯 Performing unstratified MGPS analysis...")

# Prepare MGPS data with stratification
mgps_data <- dog_pt %>%
  transmute(
    id = dog_id,
    var1 = drug,
    var2 = pt)

# Process data for MGPS
mgps_processed <- processRaw(
  data = mgps_data,
  stratify = FALSE,
  zeroes = FALSE
)

# Hyperparameter estimation
theta_initial <- data.frame(
  alpha1 = c(0.5, 1.0),
  beta1 = c(0.5, 1.0), 
  alpha2 = c(2, 3),
  beta2 = c(2, 3),
  p = c(0.1, 0.2)
)

mgps_hyperparameters <- autoHyper(
  data = mgps_processed,
  theta_init = theta_initial,
  squashed = FALSE
)

# Calculate EBGM scores
mgps_results <- ebScores(
  processed = mgps_processed,
  hyper_estimate = mgps_hyperparameters,
  quantiles = c(5, 95)
)

# Extract target drug results
mgps_signals <- mgps_results$data %>%
  filter(var1 == target_drug) %>%
  transmute(
    pt = var2,
    n_observed = N,
    n_expected = E,
    relative_risk = RR,
    prr_mgps = PRR,
    ebgm = EBGM,
    eb05 = QUANT_05,
    eb95 = QUANT_95,
    mgps_signal = n_observed >= min_reports & eb05 >= 2
  ) %>%
  arrange(desc(eb05))

cat("\n=== MGPS ANALYSIS RESULTS ===\n")
cat("MGPS signals (N≥", min_reports, ", EB05≥2):", sum(mgps_signals$mgps_signal), "\n")

message("✓ MGPS analysis complete")

# ===============================================================================
# 6. BAYESIAN CONFIDENCE PROPAGATION NEURAL NETWORK (BCPNN)
# ===============================================================================

message("🧠 Performing BCPNN analysis...")

# Enhanced BCPNN calculation
calculate_bcpnn_robust <- function(n11, n1p, np1, n_total, n_simulations = 50000) {
  
  # Calculate 2x2 table
  n10 <- n1p - n11
  n01 <- np1 - n11
  n00 <- n_total - n11 - n10 - n01
  
  # Validate inputs
  if (any(c(n11, n10, n01, n00) < 0)) {
    return(list(ic = NA, ic_lower = NA, ic_upper = NA))
  }
  
  tryCatch({
    # Dirichlet sampling with Jeffreys prior
    alpha_vec <- c(n11, n10, n01, n00) + 0.5
    
    # Generate Dirichlet samples using gamma property
    gamma_samples <- matrix(
      rgamma(n_simulations * 4, shape = rep(alpha_vec, each = n_simulations)),
      nrow = n_simulations
    )
    
    # Normalize to get probabilities
    p_samples <- gamma_samples / rowSums(gamma_samples)
    
    # Calculate Information Component
    p11 <- p_samples[, 1]
    p_drug <- p_samples[, 1] + p_samples[, 2]
    p_event <- p_samples[, 1] + p_samples[, 3]
    
    ic_samples <- log2(p11 / (p_drug * p_event))
    ic_samples <- ic_samples[is.finite(ic_samples)]
    
    if (length(ic_samples) == 0) {
      return(list(ic = NA, ic_lower = NA, ic_upper = NA))
    }
    
    list(
      ic = mean(ic_samples),
      ic_lower = quantile(ic_samples, 0.025),
      ic_upper = quantile(ic_samples, 0.975)
    )
    
  }, error = function(e) {
    list(ic = NA, ic_lower = NA, ic_upper = NA)
  })
}

# Prepare BCPNN analysis
n_total_dogs <- n_distinct(dog_pt$dog_id)
n_target_dogs <- denominators$n_target

pt_totals <- dog_pt %>%
  group_by(pt) %>%
  summarise(np1 = n_distinct(dog_id), .groups = "drop")

pt_target <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n11 = n_distinct(dog_id), .groups = "drop")

# Calculate BCPNN
bcpnn_results <- pt_totals %>%
  left_join(pt_target, by = "pt") %>%
  mutate(n11 = replace_na(n11, 0)) %>%
  filter(n11 >= min_reports) %>%
  rowwise() %>%
  mutate(
    bcpnn_stats = list(calculate_bcpnn_robust(n11, n_target_dogs, np1, n_total_dogs))
  ) %>%
  unnest_wider(bcpnn_stats) %>%
  ungroup() %>%
  mutate(
    bcpnn_signal = n11 >= min_reports & ic_lower > 0
  ) %>%
  arrange(desc(ic_lower))

cat("\n=== BCPNN ANALYSIS RESULTS ===\n")
cat("BCPNN signals (N≥", min_reports, ", IC025>0):", sum(bcpnn_results$bcpnn_signal, na.rm = TRUE), "\n")

message("✓ BCPNN analysis complete")

# ===============================================================================
# 7. SIGNAL INTEGRATION AND RANKING
# ===============================================================================

message("🔄 Integrating signals across all methods...")

# Combine all PT-level results
integrated_results <- prr_results %>%
  select(pt, prr_signal = disproportionality_signal, prr_value = prr, chi_square) %>%
  full_join(
    mgps_signals %>% select(pt, mgps_signal, ebgm, eb05, eb95),
    by = "pt"
  ) %>%
  full_join(
    bcpnn_results %>% select(pt, bcpnn_signal, ic, ic_lower, ic_upper),
    by = "pt"
  ) %>%
  mutate(
    # Count methods detecting signal
    methods_detecting = (prr_signal %||% FALSE) + 
      (mgps_signal %||% FALSE) + 
      (bcpnn_signal %||% FALSE),
    
    # Consensus signal (≥2 methods)
    consensus_signal = methods_detecting >= 2,
    
    # Any method signal  
    any_method_signal = methods_detecting >= 1,
    
    # Calculate composite signal strength score for ranking
    signal_strength = case_when(
      methods_detecting >= 2 ~ 100 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      methods_detecting == 1 ~ 50 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(signal_strength))

# Signal overlap summary
signal_summary <- tibble(
  method = c("PRR", "MGPS", "BCPNN", "Consensus (≥2)", "Any method"),
  n_signals = c(
    sum(integrated_results$prr_signal, na.rm = TRUE),
    sum(integrated_results$mgps_signal, na.rm = TRUE),
    sum(integrated_results$bcpnn_signal, na.rm = TRUE),
    sum(integrated_results$consensus_signal, na.rm = TRUE),
    sum(integrated_results$any_method_signal, na.rm = TRUE)
  )
)

cat("\n=== SIGNAL INTEGRATION SUMMARY ===\n")
print(signal_summary)

# Priority signals for clinical review
priority_signals <- integrated_results %>%
  filter(consensus_signal) %>%
  head(10)

cat("\nTop consensus signals for clinical review:", nrow(priority_signals), "\n")

message("✓ Signal integration complete")

# ===============================================================================
# 8. COMPREHENSIVE RESULTS OUTPUT
# ===============================================================================


# Create main results table with numerical values
results_table <- integrated_results %>%
  filter(any_method_signal) %>%
  mutate(
    prr_display = ifelse(is.na(prr_value), "—", as.character(round(prr_value, 2))),
    chi_display = ifelse(is.na(chi_square), "—", as.character(round(chi_square, 1))),
    ebgm_display = ifelse(is.na(ebgm), "—", as.character(round(ebgm, 2))),
    eb05_display = ifelse(is.na(eb05), "—", as.character(round(eb05, 2))),
    ic_display = ifelse(is.na(ic), "—", as.character(round(ic, 2))),
    ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2)))
  ) %>%
  select(
    PT = pt,
    PRR = prr_display,
    `Chi-square` = chi_display,
    EB05 = eb05_display,
    `IC Floor` = ic_lower_display,
  ) %>%
  slice_head(n = 20)  # Top 50 signals

# Create publication table
top_20_pts_table <- results_table %>%
  gt() %>%
  tab_header(
    title = paste("Pharmacovigilance Signal Detection Results:", str_to_title(target_drug)),
    subtitle = "Multi-Method Disproportionality Analysis - Ranked by Signal Strength"
  ) %>%
  tab_source_note(
    source_note = paste("Minimum reports:", min_reports, 
                        "| Signal thresholds: PRR≥2 & χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  cols_label(
    PT = "Preferred Term",
    PRR = "PRR",
    `Chi-square` = "χ²", 
    EB05 = "EB05",
    `IC Floor` = "IC Floor",
  ) %>%
  tab_style(style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), locations = cells_body()) %>%
  tab_style(style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), locations = cells_column_labels()) %>%
  # Bold PTs with signals
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = PT, rows = which(results_table_data$has_signal))) %>%
  # Uniform red highlighting for all significant values
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = PRR, rows = which(results_table_data$prr_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `Chi-square`, rows = which(results_table_data$chi_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = EB05, rows = which(results_table_data$eb05_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `IC Floor`, rows = which(results_table_data$ic_significant))) %>%
  cols_align(align = "center", columns = c(Reports, PRR, `Chi-square`, EB05, `IC Floor`)) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12) %>%
  tab_footnote(footnote = "Red highlighting indicates significant values", 
               locations = cells_column_labels(columns = c(PRR, `Chi-square`, EB05, `IC Floor`)))
  


tab_style(
    style = cell_fill(color = "#FFE6E6"),
    locations = cells_body(rows = Consensus == TRUE)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = PT)
  ) %>%
  cols_align(
    align = "center",
    columns = c(PRR, `Chi-square`, EB05, `IC Floor`,)
  )

top_20_pts_table

# Save results
write_csv(integrated_results, here("output", "integrated_signal_results.csv"))
write_csv(prr_results, here("output", "prr_detailed_results.csv"))
write_csv(mgps_signals, here("output", "mgps_detailed_results.csv"))
write_csv(bcpnn_results, here("output", "bcpnn_detailed_results.csv"))

if (nrow(organ_dpa) > 0) {
  write_csv(organ_dpa, here("output", "organ_system_analysis.csv"))
}

if (nrow(terms_dpa) > 0) {
  write_csv(terms_dpa, here("output", "terms_of_interest_analysis.csv"))
}

gtsave(publication_table, here("output", "tables", "signal_detection_summary.html"))

message("✓ Results output complete")

# ===============================================================================
# 9. TOP REPORTED PREFERRED TERMS TABLE
# ===============================================================================

message("📋 Creating top reported PTs table...")

# # Create table of top 20 most reported PTs for target drug
# top_pts_data <- dog_pt %>%
#   filter(drug == target_drug) %>%
#   group_by(pt) %>%
#   summarise(
#     n_dogs = n_distinct(dog_id),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_dogs)) %>%
#   slice_head(n = 20) %>%
#   mutate(
#     rank = row_number(),
#     percentage = round(n_dogs / denominators$n_target * 100, 1),
#     pt_clean = str_to_title(pt)  # Clean up PT formatting
#   )
# 
# # Create publication-quality table
# top_pts_table <- top_pts_data %>%
#   select(
#     Rank = rank,
#     `Preferred Term` = pt_clean,
#     `Dogs Affected` = n_dogs,
#     `Percentage` = percentage
#   ) %>%
#   gt() %>%
#   tab_header(
#     title = paste("Top 20 Reported Adverse Reactions:", str_to_title(target_drug)),
#     subtitle = paste("Most frequently reported Preferred Terms (N =", 
#                      format(denominators$n_target, big.mark = ","), "dogs)")
#   ) %>%
#   tab_source_note(
#     source_note = paste("Analysis date:", Sys.Date(), "| Percentage of total dogs receiving", str_to_title(target_drug))
#   ) %>%
#   cols_label(
#     Rank = "Rank",
#     `Preferred Term` = "Preferred Term",
#     `Dogs Affected` = "Dogs Affected",
#     Percentage = "% of Dogs"
#   ) %>%
#   fmt_number(
#     columns = `Dogs Affected`,
#     decimals = 0,
#     use_seps = TRUE
#   ) %>%
#   fmt_number(
#     columns = Percentage,
#     decimals = 1,
#     pattern = "{x}%"
#   ) %>%
#   # Financial Times inspired styling
#   tab_style(
#     style = list(
#       cell_fill(color = "#fff1e5"),
#       cell_text(color = "#2d2d2d")
#     ),
#     locations = cells_body()
#   ) %>%
#   tab_style(
#     style = list(
#       cell_fill(color = "#f2e6dd"),
#       cell_text(weight = "bold", color = "#2d2d2d")
#     ),
#     locations = cells_column_labels()
#   ) %>%
#   tab_style(
#     style = list(
#       cell_fill(color = "#e8ddd4"),
#       cell_text(weight = "bold", color = "#1a1a1a")
#     ),
#     locations = cells_title()
#   ) %>%
#   tab_style(
#     style = cell_text(color = "#555555"),
#     locations = cells_source_notes()
#   ) %>%
#   # Highlight top 5 rows
#   tab_style(
#     style = list(
#       cell_fill(color = "#ffe6e6"),
#       cell_text(weight = "bold")
#     ),
#     locations = cells_body(rows = 1:5)
#   ) %>%
#   cols_align(
#     align = "center",
#     columns = c(Rank, `Dogs Affected`, Percentage)
#   ) %>%
#   cols_align(
#     align = "left", 
#     columns = `Preferred Term`
#   ) %>%
#   tab_options(
#     table.font.names = "Arial",
#     table.font.size = 12,
#     heading.title.font.size = 16,
#     heading.subtitle.font.size = 12,
#     column_labels.font.weight = "bold",
#     table.border.top.style = "none",
#     table.border.bottom.style = "none",
#     column_labels.border.bottom.width = 2,
#     column_labels.border.bottom.color = "#d62728"
#   )
# 
# # Display the table
# top_pts_table
# 
# # Save the table
# gtsave(top_pts_table, here("output", "tables", "top_20_dpa.html"))
# 
# # Also create a simple CSV version for further analysis

message("✓ Top PTs table created and saved")

message("Creating main results table...")

# Top 20 PTs by report count
top_pts_by_reports <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n_reports = n_distinct(dog_id), .groups = "drop") %>%
  arrange(desc(n_reports)) %>%
  slice_head(n = 20)

# Join with signal results
results_table_data <- top_pts_by_reports %>%
  left_join(integrated_results, by = "pt") %>%
  mutate(
    prr_display = ifelse(is.na(prr_value), "—", as.character(round(prr_value, 2))),
    chi_display = ifelse(is.na(chi_square), "—", as.character(round(chi_square, 1))),
    eb05_display = ifelse(is.na(eb05), "—", as.character(round(eb05, 2))),
    ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2))),
    has_signal = any_method_signal %||% FALSE,
    prr_significant = !is.na(prr_value) & prr_value >= 2,
    chi_significant = !is.na(chi_square) & chi_square >= 4,
    eb05_significant = !is.na(eb05) & eb05 >= 2,
    ic_significant = !is.na(ic_lower) & ic_lower > 0,
    pt_clean = str_to_title(pt)
  )

# Create publication table with uniform red highlighting
publication_table <- results_table_data %>%
  select(`Preferred Term` = pt_clean, Reports = n_reports, PRR = prr_display, `Chi-square` = chi_display, 
         EB05 = eb05_display, `IC Floor` = ic_lower_display) %>%
  gt() %>%
  tab_header(
    title = paste("Top 20 Reported Adverse Reactions:", str_to_title(target_drug)),
    subtitle = "Ranked by frequency with signal detection results"
  ) %>%
  tab_source_note(
    source_note = paste("N =", format(denominators$n_target, big.mark = ","), 
                        "dogs | Significant values: PRR≥2, χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  fmt_number(columns = Reports, decimals = 0, use_seps = TRUE) %>%
  # Base styling
  tab_style(style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), locations = cells_body()) %>%
  tab_style(style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), locations = cells_column_labels()) %>%
  # Bold PTs with signals
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = `Preferred Term`, rows = which(results_table_data$has_signal))) %>%
  # Uniform red highlighting for all significant values
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = PRR, rows = which(results_table_data$prr_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `Chi-square`, rows = which(results_table_data$chi_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = EB05, rows = which(results_table_data$eb05_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `IC Floor`, rows = which(results_table_data$ic_significant))) %>%
  cols_align(align = "center", columns = c(Reports, PRR, `Chi-square`, EB05, `IC Floor`)) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12) %>%
  tab_footnote(footnote = "Red highlighting indicates significant values", 
               locations = cells_column_labels(columns = c(PRR, `Chi-square`, EB05, `IC Floor`)))

publication_table

gtsave(publication_table, here("output", "tables", "signal_detection_summary.html"))

message("✓ Results output complete")

# Top 20 PTs by significance
top_pts_by_dpa <- dog_pt %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n_reports = n_distinct(dog_id), .groups = "drop")

# Join with signal results
results_table_data <- top_pts_by_dpa %>%
  left_join(integrated_results, by = "pt") %>%
  arrange(desc(prr_value)) %>%
  slice_head(n = 20) %>%
  mutate(
    prr_display = ifelse(is.na(prr_value), "—", as.character(round(prr_value, 2))),
    chi_display = ifelse(is.na(chi_square), "—", as.character(round(chi_square, 1))),
    eb05_display = ifelse(is.na(eb05), "—", as.character(round(eb05, 2))),
    ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2))),
    has_signal = any_method_signal %||% FALSE,
    prr_significant = !is.na(prr_value) & prr_value >= 2,
    chi_significant = !is.na(chi_square) & chi_square >= 4,
    eb05_significant = !is.na(eb05) & eb05 >= 2,
    ic_significant = !is.na(ic_lower) & ic_lower > 0,
    pt_clean = str_to_title(pt)
  )

# Create publication table with uniform red highlighting
dpa_publication_table <- results_table_data %>%
  select(`Preferred Term` = pt_clean, Reports = n_reports, PRR = prr_display, `Chi-square` = chi_display, 
         EB05 = eb05_display, `IC Floor` = ic_lower_display) %>%
  gt() %>%
  tab_header(
    title = paste("Most Disproportionality:", str_to_title(target_drug)),
    subtitle = "Ranked by Proportional Reporting Ratio"
  ) %>%
  tab_source_note(
    source_note = paste("N =", format(denominators$n_target, big.mark = ","), 
                        "dogs | Significant values: PRR≥2, χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  fmt_number(columns = Reports, decimals = 0, use_seps = TRUE) %>%
  # Base styling
  tab_style(style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), locations = cells_body()) %>%
  tab_style(style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), locations = cells_column_labels()) %>%
  # Bold PTs with signals
  tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = `Preferred Term`, rows = which(results_table_data$has_signal))) %>%
  # Uniform red highlighting for all significant values
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = PRR, rows = which(results_table_data$prr_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `Chi-square`, rows = which(results_table_data$chi_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = EB05, rows = which(results_table_data$eb05_significant))) %>%
  tab_style(style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
            locations = cells_body(columns = `IC Floor`, rows = which(results_table_data$ic_significant))) %>%
  cols_align(align = "center", columns = c(Reports, PRR, `Chi-square`, EB05, `IC Floor`)) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 12) %>%
  tab_footnote(footnote = "Red highlighting indicates significant values", 
               locations = cells_column_labels(columns = c(PRR, `Chi-square`, EB05, `IC Floor`)))

dpa_publication_table

gtsave(dpa_publication_table, here("output", "tables", "top_20_dpa_summary.html"))

message("✓ Results output complete")
# ===============================================================================
# 10. VISUALIZATION AND PLOTS
# ===============================================================================

message("📈 Creating visualization plots...")

# Prepare data for volcano plot combining PRR and BCPNN results
volcano_data <- prr_results %>%
  select(pt, prr, reports_target, disproportionality_signal, chi_square) %>%
  left_join(
    bcpnn_results %>% select(pt, ic_lower),
    by = "pt"
  ) %>%
  mutate(
    log2_prr = log2(prr),
    sqrt_chi = sqrt(chi_square),
    # Color based on IC lower bound - SET FACTOR LEVELS FOR PROPER ORDER
    signal_strength = case_when(
      !is.na(ic_lower) & ic_lower > 0 & disproportionality_signal ~ "Strong Signal (IC Floor > 0 & PRR >= 2 & Chi-square >= 4)",
      disproportionality_signal ~ "Moderate Signal (PRR >= 2 & Chi-square >= 4)",
      !is.na(ic_lower) & ic_lower > 0  ~ "IC Floor Signal Only",
      TRUE ~ "No Signal"
    ),
    # Convert to factor with proper level ordering (Strong first)
    signal_strength = factor(signal_strength, levels = c(
      "Strong Signal (IC Floor > 0 & PRR >= 2 & Chi-square >= 4)",
      "Moderate Signal (PRR >= 2 & Chi-square >= 4)",
      "IC Floor Signal Only", 
      "No Signal"
    )),
    # Size based on number of reports (n11)
    n_reports = reports_target
  ) %>%
  # Remove infinite values that could cause plotting issues
  filter(is.finite(log2_prr), is.finite(sqrt_chi))

# Financial Times inspired volcano plot
volcano_plot <- volcano_data %>%
  ggplot(aes(x = log2_prr, y = sqrt_chi)) +
  geom_point(aes(color = signal_strength, size = n_reports), alpha = 0.7) +
  
  # Reference lines
  geom_hline(yintercept = sqrt(4), linetype = "dashed", 
             color = "#767676", size = 0.5) +
  geom_vline(xintercept = log2(2), linetype = "dashed", 
             color = "#767676", size = 0.5) +
  
  # Financial Times color palette
  scale_color_manual(
    values = c(
      "Strong Signal (IC Floor > 0 & PRR >= 2 & Chi-square >= 4)" = "#d62728",      # Strong red
      "Moderate Signal (PRR >= 2 & Chi-square >= 4)",    # Orange  
      "IC Floor Signal Only" = "#1f77b4",                # Blue
      "No Signal" = "#2ca02c"        # Green
    )
  ) +
  
  # Size scale for number of reports
  scale_size_continuous(
    name = "PT Reports (Bedinvetmab)",
    range = c(1, 15),
    breaks = c(3, 10, 25, 50, 100),
    labels = c("3", "10", "25", "50", "100+")
  ) +
  
  # Clean, professional theme inspired by Financial Times
  theme_minimal() +
  theme(
    # Background
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    
    # Grid
    panel.grid.major = element_line(color = "#f0f0f0", size = 0.5),
    panel.grid.minor = element_blank(),
    
    # Text
    text = element_text(family = "Arial", color = "#333333"),
    plot.title = element_text(size = 16, face = "bold", color = "#2d2d2d"),
    plot.subtitle = element_text(size = 12, color = "#555555", margin = margin(b = 20)),
    axis.title = element_text(size = 12, color = "#333333"),
    axis.text = element_text(size = 10, color = "#555555"),
    
    # Legend
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    legend.box = "horizontal",
    
    # Margins
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  
  # Labels - USE EXPRESSION() FOR PROPER MATHEMATICAL NOTATION
  labs(
    title = paste("Pharmacovigilance Signal Detection:", str_to_title(target_drug)),
    subtitle = "Proportional Reporting Ratio vs Statistical Significance (Square root of Chi-square value)\nSized by report count",
    x = expression(log[2]("Proportional Reporting Ratio")),
    y = expression(sqrt(chi^2)),
    color = "Signal Strength",
    caption = "Dashed lines: PRR ≥ 2, Chi-square ≥ 4"
  ) +
  
  # Guide specifications - CONTROL LEGEND LAYOUT
  guides(
    color = guide_legend(
      override.aes = list(size = 4, alpha = 1),
      title.position = "top",
      ncol = 1,  # 2 columns so Strong appears above Moderate
      byrow = TRUE  # Fill by row, not column
    ),
    size = guide_legend(
      title.position = "top",
      nrow = 1
    )
  )

# Display the plot
volcano_plot

ggsave(here("output", "figures", "volcano_plot.png"), volcano_plot)

# Signal comparison plot with Financial Times styling
method_comparison_data <- integrated_results %>%
  select(pt, prr_signal, mgps_signal, bcpnn_signal) %>%
  pivot_longer(cols = c(prr_signal, mgps_signal, bcpnn_signal),
               names_to = "method", values_to = "signal") %>%
  mutate(
    method = case_when(
      method == "prr_signal" ~ "Disproportionality Analysis",
      method == "mgps_signal" ~ "Multi-item Gamma Poisson Shrinker", 
      method == "bcpnn_signal" ~ "Bayesian Confidence Propogation Neural Network"
    ),
    signal_status = ifelse(signal %in% TRUE, "Signal Detected", "No Signal")
  ) %>%
  count(method, signal_status)

comparison_plot <- method_comparison_data %>%
  ggplot(aes(x = method, y = n, fill = signal_status)) +
  geom_col(position = "stack", width = 0.7, alpha = 0.8) +
  
  # Financial Times color palette
  scale_fill_manual(
    values = c(
      "Signal Detected" = "#d62728",  # FT red
      "No Signal" = "#c7c7c7"        # Light grey
    )
  ) +
  
  # Clean theme matching volcano plot
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#fff1e5", color = NA),
    panel.background = element_rect(fill = "#fff1e5", color = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#f0f0f0", size = 0.5),
    panel.grid.minor = element_blank(),
    
    text = element_text(family = "Arial", color = "#333333"),
    plot.title = element_text(size = 16, face = "bold", color = "#2d2d2d"),
    plot.subtitle = element_text(size = 12, color = "#555555", margin = margin(b = 20)),
    axis.title = element_text(size = 12, color = "#333333"),
    axis.text = element_text(size = 10, color = "#555555"),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "#fff1e5", color = NA),
    
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  
  labs(
    title = "Signal Detection Method Comparison",
    subtitle = paste("Number of Preferred Terms flagged by each method -", str_to_title(target_drug)),
    x = "Detection Method",
    y = "Number of Preferred Terms",
    caption = paste("Analysis date:", Sys.Date())
  )

# Display the comparison plot  
comparison_plot


message("✓ Plots created and displayed in R environment")
message("ℹ️  Use ggsave() to save plots if desired:")
message("   ggsave('volcano_plot.png', volcano_plot, width = 12, height = 8, dpi = 300)")
message("   ggsave('comparison_plot.png', comparison_plot, width = 10, height = 6, dpi = 300)")



# ===============================================================================
# SENSITIVITY 1: ORGAN SYSTEM ANALYSIS
# ===============================================================================

message("Performing organ system analysis...")

organ_level <- matched %>%
  filter(!is.na(organ)) %>%
  distinct(dog_id, drug, organ) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = organ, values_from = value, values_fill = 0)

# PRR analysis function for organs
calc_organ_prr <- function(organ_col, df, test_drug = target_drug) {
  df_analysis <- df %>% mutate(has_reaction = .data[[organ_col]] == 1, is_test = drug == test_drug)
  
  a <- sum(df_analysis$is_test & df_analysis$has_reaction, na.rm = TRUE)
  b <- sum(!df_analysis$is_test & df_analysis$has_reaction, na.rm = TRUE)
  n_test <- sum(df_analysis$is_test, na.rm = TRUE)
  n_other <- sum(!df_analysis$is_test, na.rm = TRUE)
  
  if (a < min_reports || b < min_reports) return(NULL)
  
  prr <- (a / n_test) / (b / n_other)
  c <- n_test - a; d <- n_other - b; total <- a + b + c + d
  expected <- c((a + c) * (a + b) / total, (a + c) * (c + d) / total, 
                (b + d) * (a + b) / total, (b + d) * (c + d) / total)
  observed <- c(a, c, b, d)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  
  tibble(organ = organ_col, n_target = a, n_comparator = b, prr = prr, chi_square = chisq,
         prr_signal = a >= min_reports & prr >= 2 & chisq >= 4)
}

# BCPNN function for organs
calc_organ_bcpnn <- function(organ_col, df, test_drug = target_drug) {
  df_analysis <- df %>% mutate(has_reaction = .data[[organ_col]] == 1, is_test = drug == test_drug)
  
  n11 <- sum(df_analysis$is_test & df_analysis$has_reaction, na.rm = TRUE)
  n1p <- sum(df_analysis$is_test, na.rm = TRUE)
  np1 <- sum(df_analysis$has_reaction, na.rm = TRUE)
  n_total <- nrow(df_analysis)
  
  if (n11 < min_reports) return(NULL)
  
  bcpnn_stats <- calculate_bcpnn_robust(n11, n1p, np1, n_total)
  tibble(organ = organ_col, n11 = n11, ic = bcpnn_stats$ic, ic_lower = bcpnn_stats$ic_lower,
         bcpnn_signal = n11 >= min_reports & bcpnn_stats$ic_lower > 0)
}

# Analyze organs
organ_cols <- setdiff(names(organ_level), c("dog_id", "drug"))
organ_prr <- map_dfr(organ_cols, ~ calc_organ_prr(.x, organ_level))
organ_bcpnn <- map_dfr(organ_cols, ~ calc_organ_bcpnn(.x, organ_level))

# Integrate organ results
organ_integrated <- organ_prr %>%
  full_join(organ_bcpnn, by = "organ") %>%
  mutate(methods_detecting = (prr_signal %||% FALSE) + (bcpnn_signal %||% FALSE),
         consensus_signal = methods_detecting >= 2, any_signal = methods_detecting >= 1) %>%
  arrange(desc(methods_detecting), desc(prr))

cat("Organ systems with signals:", sum(organ_integrated$any_signal, na.rm = TRUE), "\n")

# Create organ table
if (nrow(organ_integrated) > 0) {
  organ_table <- organ_integrated %>%
    mutate(organ_clean = str_replace_all(organ, "_", " ") %>% str_to_title(),
           prr_display = ifelse(is.na(prr), "—", as.character(round(prr, 2))),
           ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2)))) %>%
    select(`Organ System` = organ_clean, `Target Reports` = n_target, PRR = prr_display, 
           `IC Lower` = ic_lower_display, Consensus = consensus_signal) %>%
    slice_head(n = 15) %>%
    gt() %>%
    tab_header(title = paste("Organ System Analysis:", str_to_title(target_drug)), 
               subtitle = "Disproportionality analysis by system organ class") %>%
    tab_style(style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), locations = cells_body()) %>%
    tab_style(style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), locations = cells_column_labels())
  
  organ_table
  write_csv(organ_integrated, here("output", "organ_system_analysis.csv"))
  gtsave(organ_table, here("output", "tables", "organ_system_analysis.html"))
}

# ===============================================================================
# SENSITIVITY 2: HIGH-LEVEL TERM ANALYSIS
# ===============================================================================

message("Performing High-Level Term analysis...")

# Create HLT dataset
dog_hlt <- matched %>%
  filter(!is.na(hlt)) %>%
  distinct(dog_id, drug, hlt, .keep_all = TRUE)

# HLT PRR analysis
hlt_prr_results <- dog_hlt %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, hlt) %>%
  group_by(hlt) %>%
  summarise(reports_target = sum(is_target), reports_comparator = sum(!is_target), .groups = "drop") %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(stats = list(calculate_prr_comprehensive(reports_target, denominators$n_target, 
                                                  reports_comparator, denominators$n_comparator))) %>%
  unnest_wider(stats) %>%
  mutate(prr_signal = prr >= 2 & chi_square >= 4) %>%
  arrange(desc(prr))

# HLT BCPNN analysis
hlt_totals <- dog_hlt %>% group_by(hlt) %>% summarise(np1 = n_distinct(dog_id), .groups = "drop")
hlt_target <- dog_hlt %>% filter(drug == target_drug) %>% group_by(hlt) %>% summarise(n11 = n_distinct(dog_id), .groups = "drop")

hlt_bcpnn_results <- hlt_totals %>%
  left_join(hlt_target, by = "hlt") %>%
  mutate(n11 = replace_na(n11, 0)) %>%
  filter(n11 >= min_reports) %>%
  rowwise() %>%
  mutate(bcpnn_stats = list(calculate_bcpnn_robust(n11, denominators$n_target, np1, n_distinct(dog_hlt$dog_id)))) %>%
  unnest_wider(bcpnn_stats) %>%
  mutate(bcpnn_signal = n11 >= min_reports & ic_lower > 0) %>%
  arrange(desc(ic_lower))

# Integrate HLT results
hlt_integrated <- hlt_prr_results %>%
  select(hlt, prr_signal, prr, chi_square) %>%
  full_join(hlt_bcpnn_results %>% select(hlt, bcpnn_signal, ic, ic_lower), by = "hlt") %>%
  mutate(methods_detecting = (prr_signal %||% FALSE) + (bcpnn_signal %||% FALSE),
         consensus_signal = methods_detecting >= 2, any_signal = methods_detecting >= 1) %>%
  arrange(desc(methods_detecting), desc(prr))

cat("HLTs with signals:", sum(hlt_integrated$any_signal, na.rm = TRUE), "\n")

# Create HLT table
if (nrow(hlt_integrated) > 0 && sum(hlt_integrated$any_signal) > 0) {
  hlt_table <- hlt_integrated %>%
    filter(any_signal) %>%
    mutate(hlt_clean = str_to_title(hlt),
           prr_display = ifelse(is.na(prr), "—", as.character(round(prr, 2))),
           ic_lower_display = ifelse(is.na(ic_lower), "—", as.character(round(ic_lower, 2)))) %>%
    select(`High-Level Term` = hlt_clean, PRR = prr_display, `IC Lower` = ic_lower_display, 
           Methods = methods_detecting, Consensus = consensus_signal) %>%
    slice_head(n = 15) %>%
    gt() %>%
    tab_header(title = paste("High-Level Term Analysis:", str_to_title(target_drug)), 
               subtitle = "Signal detection at High-Level Term classification") %>%
    tab_style(style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), locations = cells_body()) %>%
    tab_style(style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), locations = cells_column_labels())
  
  hlt_table
  write_csv(hlt_integrated, here("output", "hlt_analysis.csv"))
  gtsave(hlt_table, here("output", "tables", "hlt_analysis.html"))
}

# ===============================================================================
# SENSITIVITY 3: MONOTHERAPY ANALYSIS
# ===============================================================================

message("Performing monotherapy analysis...")

# Filter to monotherapy dogs
dog_pt_mono <- dog_pt %>% filter(polypharmacy == "Monotherapy")

mono_denominators <- dog_pt_mono %>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target) %>%
  summarise(
    n_target = sum(is_target),
    n_comparator = sum(!is_target),
    .groups = "drop"
  )

# Perform PRR analysis
prr_mono_results <- dog_pt_mono%>%
  mutate(is_target = drug == target_drug) %>%
  distinct(dog_id, is_target, pt) %>%
  group_by(pt) %>%
  summarise(
    reports_target = sum(is_target),
    reports_comparator = sum(!is_target),
    .groups = "drop"
  ) %>%
  filter(reports_target >= min_reports, reports_comparator >= min_reports) %>%
  rowwise() %>%
  mutate(
    stats = list(calculate_prr_comprehensive(reports_target, mono_denominators$n_target, 
                                             reports_comparator, mono_denominators$n_comparator))
  ) %>%
  unnest_wider(stats) %>%
  ungroup() %>%
  # Apply FDR correction
  mutate(
    disproportionality_signal = prr >= 2 & chi_square >= 4
  ) %>%
  arrange(desc(prr))

cat("\n=== PRR ANALYSIS RESULTS ===\n")
cat("Total PTs analyzed:", nrow(prr_mono_results), "\n")
cat("disproportionality signals (PRR ≥2, χ² ≥4):", sum(prr_mono_results$disproportionality_signal), "\n")

message("✓ PRR analysis complete")

# ===============================================================================
# 5. MULTI-ITEM GAMMA POISSON SHRINKER (MGPS) ANALYSIS
# ===============================================================================

message("🎯 Performing unstratified MGPS analysis...")

# Prepare MGPS data with stratification
mgps_mono_data <- dog_pt_mono %>%
  transmute(
    id = dog_id,
    var1 = drug,
    var2 = pt)

# Process data for MGPS
mgps_mono_processed <- processRaw(
  data = mgps_mono_data,
  stratify = FALSE,
  zeroes = FALSE
)

# Hyperparameter estimation
theta_initial <- data.frame(
  alpha1 = c(0.5, 1.0, 1.5),
  beta1 = c(0.5, 1.0, 1.5), 
  alpha2 = c(2, 3, 4),
  beta2 = c(2, 3, 4),
  p = c(0.1, 0.2, 0.3)
)

mgps_hyperparameters <- autoHyper(
  data = mgps_mono_processed,
  theta_init = theta_initial,
  squashed = FALSE
)

# Calculate EBGM scores
mgps_mono_results <- ebScores(
  processed = mgps_mono_processed,
  hyper_estimate = mgps_hyperparameters,
  quantiles = c(5, 95)
)

# Extract target drug results
mgps_mono_signals <- mgps_mono_results$data %>%
  filter(var1 == target_drug) %>%
  transmute(
    pt = var2,
    n_observed = N,
    n_expected = E,
    relative_risk = RR,
    prr_mgps = PRR,
    ebgm = EBGM,
    eb05 = QUANT_05,
    eb95 = QUANT_95,
    mgps_signal = n_observed >= min_reports & eb05 >= 2
  ) %>%
  arrange(desc(eb05))

cat("\n=== MGPS ANALYSIS RESULTS ===\n")
cat("MGPS signals (N ≥", min_reports, ", EB05≥2):", sum(mgps_mono_signals$mgps_signal), "\n")

message("✓ MGPS analysis complete")

# ===============================================================================
# 6. BAYESIAN CONFIDENCE PROPAGATION NEURAL NETWORK (BCPNN)
# ===============================================================================

message("🧠 Performing BCPNN analysis...")

# Enhanced BCPNN calculation
calculate_bcpnn_robust <- function(n11, n1p, np1, n_total, n_simulations = 50000) {
  
  # Calculate 2x2 table
  n10 <- n1p - n11
  n01 <- np1 - n11
  n00 <- n_total - n11 - n10 - n01
  
  # Validate inputs
  if (any(c(n11, n10, n01, n00) < 0)) {
    return(list(ic = NA, ic_lower = NA, ic_upper = NA))
  }
  
  tryCatch({
    # Dirichlet sampling with Jeffreys prior
    alpha_vec <- c(n11, n10, n01, n00) + 0.5
    
    # Generate Dirichlet samples using gamma property
    gamma_samples <- matrix(
      rgamma(n_simulations * 4, shape = rep(alpha_vec, each = n_simulations)),
      nrow = n_simulations
    )
    
    # Normalize to get probabilities
    p_samples <- gamma_samples / rowSums(gamma_samples)
    
    # Calculate Information Component
    p11 <- p_samples[, 1]
    p_drug <- p_samples[, 1] + p_samples[, 2]
    p_event <- p_samples[, 1] + p_samples[, 3]
    
    ic_samples <- log2(p11 / (p_drug * p_event))
    ic_samples <- ic_samples[is.finite(ic_samples)]
    
    if (length(ic_samples) == 0) {
      return(list(ic = NA, ic_lower = NA, ic_upper = NA))
    }
    
    list(
      ic = mean(ic_samples),
      ic_lower = quantile(ic_samples, 0.025),
      ic_upper = quantile(ic_samples, 0.975)
    )
    
  }, error = function(e) {
    list(ic = NA, ic_lower = NA, ic_upper = NA)
  })
}

# Prepare BCPNN analysis
n_total_dogs <- n_distinct(dog_pt_mono$dog_id)
n_target_dogs <- denominators$n_target

pt_totals <- dog_pt_mono %>%
  group_by(pt) %>%
  summarise(np1 = n_distinct(dog_id), .groups = "drop")

pt_target <- dog_pt_mono %>%
  filter(drug == target_drug) %>%
  group_by(pt) %>%
  summarise(n11 = n_distinct(dog_id), .groups = "drop")

# Calculate BCPNN
bcpnn_mono_results <- pt_totals %>%
  left_join(pt_target, by = "pt") %>%
  mutate(n11 = replace_na(n11, 0)) %>%
  filter(n11 >= min_reports) %>%
  rowwise() %>%
  mutate(
    bcpnn_stats = list(calculate_bcpnn_robust(n11, n_target_dogs, np1, n_total_dogs))
  ) %>%
  unnest_wider(bcpnn_stats) %>%
  ungroup() %>%
  mutate(
    bcpnn_signal = n11 >= min_reports & ic_lower > 0
  ) %>%
  arrange(desc(ic_lower))

cat("\n=== BCPNN ANALYSIS RESULTS ===\n")
cat("BCPNN signals (N≥", min_reports, ", IC025>0):", sum(bcpnn_mono_results$bcpnn_signal, na.rm = TRUE), "\n")

message("✓ BCPNN analysis complete")

# ===============================================================================
# 7. SIGNAL INTEGRATION AND RANKING
# ===============================================================================

message("🔄 Integrating signals across all methods...")

# Combine all PT-level results
integrated_mono_results <- prr_mono_results %>%
  select(pt, prr_signal = disproportionality_signal, prr_value = prr, chi_square) %>%
  full_join(
    mgps_mono_signals %>% select(pt, mgps_signal, ebgm, eb05, eb95),
    by = "pt"
  ) %>%
  full_join(
    bcpnn_mono_results %>% select(pt, bcpnn_signal, ic, ic_lower, ic_upper),
    by = "pt"
  ) %>%
  mutate(
    # Count methods detecting signal
    methods_detecting = (prr_signal %||% FALSE) + 
      (mgps_signal %||% FALSE) + 
      (bcpnn_signal %||% FALSE),
    
    # Consensus signal (≥2 methods)
    consensus_signal = methods_detecting >= 2,
    
    # Any method signal  
    any_method_signal = methods_detecting >= 1,
    
    # Calculate composite signal strength score for ranking
    signal_strength = case_when(
      methods_detecting >= 2 ~ 100 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      methods_detecting == 1 ~ 50 + (prr_value %||% 0) + (eb05 %||% 0) + 
        (pmax(0, ic_lower, na.rm = TRUE) * 10),
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(signal_strength))

# Signal overlap summary
mono_signal_summary <- tibble(
  method = c("PRR", "MGPS", "BCPNN", "Consensus (≥2)", "Any method"),
  n_signals = c(
    sum(integrated_mono_results$prr_signal, na.rm = TRUE),
    sum(integrated_mono_results$mgps_signal, na.rm = TRUE),
    sum(integrated_mono_results$bcpnn_signal, na.rm = TRUE),
    sum(integrated_mono_results$consensus_signal, na.rm = TRUE),
    sum(integrated_mono_results$any_method_signal, na.rm = TRUE)
  )
)

cat("\n=== SIGNAL INTEGRATION SUMMARY ===\n")
print(mono_signal_summary)

# Join with both full and monotherapy results
joint_results_table_data <- top_pts_by_dpa %>%
  left_join(integrated_results, by = "pt", suffix = c("", "_full")) %>%
  left_join(integrated_mono_results, by = "pt", suffix = c("_full", "_mono")) %>%
  arrange(desc(prr_value_mono)) %>%
  slice_head(n = 20) %>%
  mutate(
    # Full analysis displays
    prr_full_display = ifelse(is.na(prr_value_full), "—", as.character(round(prr_value_full, 2))),
    chi_full_display = ifelse(is.na(chi_square_full), "—", as.character(round(chi_square_full, 1))),
    eb05_full_display = ifelse(is.na(eb05_full), "—", as.character(round(eb05_full, 2))),
    ic_full_display = ifelse(is.na(ic_lower_full), "—", as.character(round(ic_lower_full, 2))),
    
    # Monotherapy analysis displays
    prr_mono_display = ifelse(is.na(prr_value_mono), "—", as.character(round(prr_value_mono, 2))),
    chi_mono_display = ifelse(is.na(chi_square_mono), "—", as.character(round(chi_square_mono, 1))),
    eb05_mono_display = ifelse(is.na(eb05_mono), "—", as.character(round(eb05_mono, 2))),
    ic_mono_display = ifelse(is.na(ic_lower_mono), "—", as.character(round(ic_lower_mono, 2))),
    
    # Signal flags
    has_signal_full = any_method_signal_full %||% FALSE,
    has_signal_mono = any_method_signal_mono %||% FALSE,
    
    # Full analysis significance
    prr_full_sig = !is.na(prr_value_full) & prr_value_full >= 2,
    chi_full_sig = !is.na(chi_square_full) & chi_square_full >= 4,
    eb05_full_sig = !is.na(eb05_full) & eb05_full >= 2,
    ic_full_sig = !is.na(ic_lower_full) & ic_lower_full > 0,
    
    # Monotherapy significance
    prr_mono_sig = !is.na(prr_value_mono) & prr_value_mono >= 2,
    chi_mono_sig = !is.na(chi_square_mono) & chi_square_mono >= 4,
    eb05_mono_sig = !is.na(eb05_mono) & eb05_mono >= 2,
    ic_mono_sig = !is.na(ic_lower_mono) & ic_lower_mono > 0,
    
    pt_clean = str_to_title(pt)
  )


# Create publication table with dual analysis columns
joint_dpa_publication_table <- joint_results_table_data %>%
  select(
    `Preferred Term` = pt_clean, 
    Reports = n_reports,
    # Full analysis
    `PRR Full` = prr_full_display, 
    `Chi-square Full` = chi_full_display, 
    `EB05 Full` = eb05_full_display, 
    `IC Floor Full` = ic_full_display,
    # Monotherapy analysis
    `PRR Mono` = prr_mono_display,
    `Chi-square Mono` = chi_mono_display,
    `EB05 Mono` = eb05_mono_display,
    `IC Floor Mono` = ic_mono_display
  ) %>%
  gt() %>%
  tab_header(
    title = paste("Most Disproportionality:", str_to_title(target_drug)),
    subtitle = "Ranked by Proportional Reporting Ratio (Full Analysis)"
  ) %>%
  tab_source_note(
    source_note = paste("N =", format(mono_denominators$n_target, big.mark = ","), 
                        "dogs (full) | Significant values: PRR≥2, χ²≥4, EB05≥2, IC₀₂₅>0")
  ) %>%
  # Add column spanners
  tab_spanner(
    label = "Full Analysis",
    columns = c(`PRR Full`, `Chi-square Full`, `EB05 Full`, `IC Floor Full`)
  ) %>%
  tab_spanner(
    label = "Monotherapy Analysis",
    columns = c(`PRR Mono`, `Chi-square Mono`, `EB05 Mono`, `IC Floor Mono`)
  ) %>%
  # Rename column labels to remove suffixes
  cols_label(
    `PRR Full` = "PRR",
    `Chi-square Full` = "Chi-square",
    `EB05 Full` = "EB05",
    `IC Floor Full` = "IC Floor",
    `PRR Mono` = "PRR",
    `Chi-square Mono` = "Chi-square",
    `EB05 Mono` = "EB05",
    `IC Floor Mono` = "IC Floor"
  ) %>%
  fmt_number(columns = Reports, decimals = 0, use_seps = TRUE) %>%
  # Base styling
  tab_style(
    style = list(cell_fill(color = "#fff1e5"), cell_text(color = "#2d2d2d")), 
    locations = cells_body()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#f2e6dd"), cell_text(weight = "bold", color = "#2d2d2d")), 
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#e8ddd4"), cell_text(weight = "bold", color = "#1a1a1a")),
    locations = cells_column_spanners()
  ) %>%
  # Bold PTs with signals in either analysis
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_body(
      columns = `Preferred Term`, 
      rows = which(joint_results_table_data$has_signal_full | joint_results_table_data$has_signal_mono)
    )
  ) %>%
  # Full analysis highlighting
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `PRR Full`, rows = which(joint_results_table_data$prr_full_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `Chi-square Full`, rows = which(joint_results_table_data$chi_full_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `EB05 Full`, rows = which(joint_results_table_data$eb05_full_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `IC Floor Full`, rows = which(joint_results_table_data$ic_full_sig))
  ) %>%
  # Monotherapy analysis highlighting
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `PRR Mono`, rows = which(joint_results_table_data$prr_mono_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `Chi-square Mono`, rows = which(joint_results_table_data$chi_mono_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `EB05 Mono`, rows = which(joint_results_table_data$eb05_mono_sig))
  ) %>%
  tab_style(
    style = list(cell_fill(color = "#ffcccc"), cell_text(weight = "bold", color = "#d62728")), 
    locations = cells_body(columns = `IC Floor Mono`, rows = which(joint_results_table_data$ic_mono_sig))
  ) %>%
  cols_align(
    align = "center", 
    columns = c(Reports, `PRR Full`, `Chi-square Full`, `EB05 Full`, `IC Floor Full`,
                `PRR Mono`, `Chi-square Mono`, `EB05 Mono`, `IC Floor Mono`)
  ) %>%
  cols_align(align = "left", columns = `Preferred Term`) %>%
  tab_options(table.font.names = "Arial", table.font.size = 11) %>%
  tab_footnote(
    footnote = "Red highlighting indicates significant values", 
    locations = cells_column_spanners(spanners = c("Full Analysis", "Monotherapy Analysis"))
  )

joint_dpa_publication_table

# ===============================================================================
# 12. EXECUTIVE SUMMARY AND RECOMMENDATIONS
# ===============================================================================

message("📋 Generating executive summary...")

# Create executive summary
executive_summary <- list(
  analysis_date = as.character(Sys.Date()),
  target_drug = target_drug,
  total_dogs_analyzed = denominators$n_target + denominators$n_comparator,
  target_drug_dogs = denominators$n_target,
  comparator_drug_dogs = denominators$n_comparator,
  total_reactions_analyzed = nrow(matched),
  pts_analyzed = nrow(integrated_results),
  
  # Signal counts by method
  prr_signals = sum(integrated_results$prr_signal, na.rm = TRUE),
  mgps_signals = sum(integrated_results$mgps_signal, na.rm = TRUE),
  bcpnn_signals = sum(integrated_results$bcpnn_signal, na.rm = TRUE),
  consensus_signals = sum(integrated_results$consensus_signal, na.rm = TRUE),
  
  # Top signals
  top_consensus_signals = if(nrow(priority_signals) > 0) {
    priority_signals$pt[1:min(5, nrow(priority_signals))]
  } else character(0),
  
  # Additional analyses
  organ_signals = if(exists("organ_dpa")) sum(organ_dpa$signal, na.rm = TRUE) else 0,
  terms_of_interest_signals = if(exists("terms_dpa")) sum(terms_dpa$signal, na.rm = TRUE) else 0,
  
  # Clinical recommendations
  immediate_review_recommended = sum(integrated_results$consensus_signal, na.rm = TRUE) > 0,
  
  # Analysis quality metrics
  veddra_mapping_rate = round(mean(!is.na(matched$pt)) * 100, 1),
  analysis_completion_status = "SUCCESS"
)

# Save executive summary
write_json(executive_summary, here("output", "executive_summary.json"), pretty = TRUE)

# Display key findings
cat("\n=== EXECUTIVE SUMMARY ===\n")
cat("Analysis Date:", executive_summary$analysis_date, "\n")
cat("Target Drug:", str_to_title(executive_summary$target_drug), "\n")
cat("Dogs Analyzed:", format(executive_summary$total_dogs_analyzed, big.mark = ","), "\n")
cat("Consensus Signals:", executive_summary$consensus_signals, "\n")
cat("Immediate Clinical Review:", ifelse(executive_summary$immediate_review_recommended, "YES", "NO"), "\n")

if (length(executive_summary$top_consensus_signals) > 0) {
  cat("Top Priority Terms:\n")
  cat(paste("-", executive_summary$top_consensus_signals), sep = "\n")
}

# Final completion log
final_log <- paste0(
  "\n=== ANALYSIS COMPLETION ===\n",
  "Completion time: ", Sys.time(), "\n",
  "Status: SUCCESS\n",
  "Consensus signals: ", executive_summary$consensus_signals, "\n",
  "Files generated in: output/\n"
)

cat(final_log, file = log_file, append = TRUE)
cat(final_log)

message("✅ Comprehensive pharmacovigilance analysis complete!")
message("📁 Results available in output/ directory")
message("📊 Review executive_summary.json and signal_detection_summary.html")
message("🏥 ", executive_summary$consensus_signals, " consensus signals require clinical review")

# ===============================================================================
# END OF ANALYSIS
# ===============================================================================