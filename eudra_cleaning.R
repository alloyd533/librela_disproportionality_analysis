# vet_data_cleaning.R
library(tidyverse)
library(readxl)
library(lubridate)

# Define drug names
drugs <- c("librela", "metacam", "previcox", "galliprant", "rimadyl", "onsior", "daxocox")

# List all yearly CSV files
file_list <- list.files("data", pattern = "^\\d{4}_eudra_.*\\.csv$", full.names = TRUE)

# Read and clean individual CSV file
read_and_clean <- function(path) {
  df <- read_csv(path, col_types = cols(.default = "c"), show_col_types = FALSE)
  
  # Standardise column names
  names(df) <- str_replace_all(names(df), " ", "_")
  
  # Parse dates
  parsed <- parse_date_time(df$Received_date, orders = c("ymd", "dmy", "ymd HMS", "dmy HMS"), quiet = TRUE)
  if (all(is.na(parsed))) {
    parsed <- suppressWarnings(as.Date(df$Received_date, format = "%Y-%m-%d"))
    if (all(is.na(parsed))) {
      parsed <- suppressWarnings(as.Date(df$Received_date, format = "%d/%m/%Y"))
    }
  }
  
  df <- df |>
    mutate(
      date = parsed,
      source_file = basename(path),
      year = str_extract(source_file, "^\\d{4}")
    ) |>
    select(-Received_date)
  
  return(df)
}

# Read and save cleaned yearly data per drug
for (drug in drugs) {
  drug_files <- keep(file_list, ~ str_detect(.x, paste0("eudra_", drug, "\\.csv$")))
  if (length(drug_files) > 0) {
    combined <- map_dfr(drug_files, read_and_clean)
    write_csv(combined, file = file.path("data", paste0(drug, "_combined.csv")))
  }
}

# Re-load cleaned drug data and apply filters
fx_filter_combine <- function(df, drug_name) {
  df |>
    filter(Species == "Dog") |>
    select(-c(Case_number, AER_form, source_file, Animals_treated)) |>
    mutate(drug = drug_name) |>
    select(drug, year, everything())
}

# Load each cleaned drug dataset, filter and assign
for (drug in drugs) {
  df <- read_csv(paste0("data/", drug, "_combined.csv"), show_col_types = FALSE)
  assign(drug, fx_filter_combine(df, drug))
}

# Combine all drug data
complete <- bind_rows(mget(drugs))

# Save combined dog-level, reaction-wide dataset
write_csv(complete, "data/complete_combined_clean.csv")
for (drug in drugs) {
  df <- get(drug)
  df_clean <- fx_filter_combine(df, drug)
  assign(drug, df_clean)
}

# Bind them all together

complete <- bind_rows(mget(drugs))

### Pivot longer to separate reactions

max_reactions <- complete$Reaction |>
  str_count(",") |>
  max(na.rm = TRUE) + 1

reactions <- paste0("reaction", seq_len(max_reactions))

complete_long <- complete |>
  mutate(dog_id = row_number()) |>
  separate_longer_delim(Reaction, delim = ",") |>
  mutate(llt = str_trim(Reaction),
         llt = str_to_lower(llt)) |>
  select(-Reaction)

# Now join complete with veddra to get the terms

matched <- complete_long |>
  left_join(veddra, by = "llt")

# Ensure each HLT maps to one organ system (if this is valid)
hlt_map <- matched |>
  select(hlt, organ) |>
  distinct()

# Now count occurrences per HLT × drug
hlt_counts <- matched |>
  count(hlt, drug, name = "n")

# Join organ system back on using HLT
hlt_props <- hlt_counts |>
  left_join(hlt_map, by = "hlt") |>
  
  # Get total reactions per drug (for proportions)
  group_by(drug) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

hlt_props_filtered <- hlt_props |>
  filter(organ%in% c("Musculoskeletal disorders", "Neurological disorders")) |>
  select(organ, hlt, drug, prop)



# Pivot wider: HLTs as rows, drugs as columns
hlt_prop_wide <- hlt_props_filtered |>
  mutate(drug = str_to_sentence(drug)) |>
  pivot_wider(
    names_from = drug,
    values_from = prop,
    values_fill = 0
  ) 
  
require(gt)
# Create nicely formatted gt table
hlt_prop_wide |>
  arrange(organ, hlt) |>
  mutate(across(where(is.numeric), ~ scales::percent(.x, accuracy = 0.1))) |>
  gt() |>
  tab_header(
    title = "Proportional Reporting Table by High-Level Term (HLT)",
    subtitle = "Each cell = % of total reactions for that drug with that HLT"
  ) |>
  tab_row_group(
    label = "MSK",
    rows = organ == "Musculoskeletal disorders"
  ) |>
  tab_row_group(
    label = "Neurological",
    rows = organ == "Neurological disorders"
  ) |>
  cols_label(
    hlt = "High-Level Term"
  ) |>
  cols_hide(columns = organ)|>
  
  # ✅ Add styling to the row group labels
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  )


# Create a binary for any of reactions are MSK or Neuro
organ_level <- matched |>
  group_by(dog_id, drug) |>
  summarise(
    MSK = as.integer(any(str_detect(organ, regex("musculoskeletal", ignore_case = TRUE)))),
    Neuro = as.integer(any(str_detect(organ, regex("neurological", ignore_case = TRUE)))),
    .groups = "drop"
  ) 
  
organ_table <- organ_level |>
  group_by(drug) |>
  summarise(
    MSK_prop   = mean(MSK, na.rm = TRUE),
    Neuro_prop = mean(Neuro, na.rm = TRUE),
    .groups = "drop"
  )

organ_table |>
  mutate(drug = str_to_sentence(drug)) |>
  mutate(across(where(is.numeric), ~ scales::percent(.x, accuracy = 0.1))) |>
  gt::gt() |>
  gt::tab_header(
    title = "Proportion of Dogs with MSK or Neuro Reactions",
    subtitle = "Binary indicator based on whether dog had any reaction in that organ system"
  ) |>
  gt::cols_label(
    MSK_prop = "Musculoskeletal",
    Neuro_prop = "Neurological",
    drug = "Drug"
  )

## Disproportionality analysis

calc_dpa_librela <- function(hlt_name, df) {
  df_bin <- df |>
    mutate(
      has_reaction = hlt == hlt_name,
      is_librela = drug == "librela"
    )
  
  # 2x2 counts
  a <- sum(df_bin$is_librela & df_bin$has_reaction, na.rm = TRUE)
  b <- sum(df_bin$is_librela & !df_bin$has_reaction, na.rm = TRUE)
  c <- sum(!df_bin$is_librela & df_bin$has_reaction, na.rm = TRUE)
  d <- sum(!df_bin$is_librela & !df_bin$has_reaction, na.rm = TRUE)
  
  # Drop if unstable (low counts)
  if (a < 5 || c < 5) return(NULL)
  
  # PRR
  prr <- (a / (a + b)) / (c / (c + d))
  
  a <- as.numeric(a)
  b <- as.numeric(b)
  c <- as.numeric(c)
  d <- as.numeric(d)
  
  # Chi-squared (without continuity correction)
  total <- a + b + c + d
  expected <- c(
    (a + b) * (a + c) / total,  # E[a]
    (a + b) * (b + d) / total,  # E[b]
    (c + d) * (a + c) / total,  # E[c]
    (c + d) * (b + d) / total   # E[d]
  )
  observed <- c(a, b, c, d)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  
  # Signal flag per FDA-style rule
  signals <- a >= 3 & prr >= 2 & chisq >= 4
  
  tibble(
    hlt = hlt_name,
    a, b, c, d,
    prr,
    chisq,
    signals
  )
}



dog_hlt_level <- matched |>
  select(dog_id, drug, hlt) |>
  distinct()

librela_hlts <- dog_hlt_level |>
  filter(drug == "librela") |>
  distinct(hlt) |>
  drop_na()

librela_dpa <- map_dfr(librela_hlts$hlt, ~ calc_dpa_librela(.x, dog_hlt_level))

library(ggplot2)

librela_dpa |>
  filter(prr > 1.5) |>
  mutate(
    hlt = fct_reorder(hlt, prr),
    sig_label = if_else(signals, "Signal", "No signal")
  ) |>
  ggplot(aes(x = prr, y = hlt, fill = sig_label)) +
  geom_col() +
  scale_fill_manual(values = c("Signal" = "red", "No signal" = "grey80")) +
  labs(
    title = "HLTs with Disproportional Reporting for Librela",
    x = "PRR",
    y = "High-Level Term (HLT)",
    fill = "Signal status"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


## DPA at dog level

organ_matrix <- matched |>
  filter(!is.na(organ)) |>
  distinct(dog_id, drug, organ) |>
  mutate(value = 1) |>
  pivot_wider(
    names_from = organ,
    values_from = value,
    values_fill = 0
  )

calc_prop_dpa <- function(organ_col, df, test_drug = "librela") {
  df_bin <- df |>
    mutate(
      has_reaction = .data[[organ_col]] == 1,
      is_test = drug == test_drug
    )
  
  a <- sum(df_bin$is_test & df_bin$has_reaction, na.rm = TRUE)
  b <- sum(!df_bin$is_test & df_bin$has_reaction, na.rm = TRUE)
  n_test <- sum(df_bin$is_test, na.rm = TRUE)
  n_other <- sum(!df_bin$is_test, na.rm = TRUE)
  
  # Drop if too few cases
  if (a < 5 || b < 5) return(NULL)
  
  prr <- (a / n_test) / (b / n_other)
  
  # Chi-squared (2x2)
  c <- n_test - a
  d <- n_other - b
  total <- a + b + c + d
  expected <- c(
    (a + c) * (a + b) / total,
    (a + c) * (c + d) / total,
    (b + d) * (a + b) / total,
    (b + d) * (c + d) / total
  )
  observed <- c(a, c, b, d)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  
  tibble(
    organ = organ_col,
    librela_prop = a / n_test,
    other_prop   = b / n_other,
    prr,
    chisq,
    signals = a >= 3 & prr >= 2 & chisq >= 4
  )
}


organ_cols <- setdiff(names(organ_matrix), c("dog_id", "drug"))

organ_dpa <- map_dfr(organ_cols, ~ calc_prop_dpa(.x, organ_matrix))

library(scales)

organ_dpa |>
  mutate(
    librela_prop = percent(librela_prop, accuracy = 0.1),
    other_prop   = percent(other_prop, accuracy = 0.1),
    prr = round(prr, 2),
    chisq = round(chisq, 1),
    signals = if_else(signals, "⚠️", "")
  ) |>
  arrange(desc(prr)) |>
  gt::gt() |>
  gt::tab_header(
    title = "Organ-Level Disproportionality: Librela vs All Others",
    subtitle = "Proportion of dogs with any reaction in each system"
  ) |>
  gt::cols_label(
    organ = "Organ System",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR",
    chisq = "Chi-sq",
    signals = "Signal"
  )


organ_dpa |>
  mutate(
    log_prr = log2(prr),
    sqrt_chisq = sqrt(chisq)
  ) |>
  ggplot(aes(x = log_prr, y = sqrt_chisq, label = organ)) +
  geom_point(aes(color = signals)) +
  geom_hline(yintercept = sqrt(4), linetype = "dashed") +
  geom_vline(xintercept = log2(2), linetype = "dashed") +
  geom_text(aes(label = ifelse(signals, organ, "")), hjust = 1.1, vjust = 0.5, size = 3) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  labs(
    x = "log2(PRR)",
    y = "sqrt(Chi-squared)",
    title = "Volcano Plot of Organ-Level Disproportionality",
    subtitle = "Highlighting organs with disproportional signals for Librela"
  ) +
  theme_minimal()


# Terms of interest from FDA study
pts_of_interest <- c(
  "Ataxia", "Convulsion", "Paresis", "Proprioception abnormality",
  "Paralysis", "Recumbency", "Muscle weakness", "Muscle tremors",
  "Lameness", "Collapse NOS",
  "Pancreatitis", "Death",
  "Immune mediated haemolytic anaemia", 
  "Immune mediated thrombocytopenia", 
  "Immune mediated polyarthritis"
) |> str_to_lower()

dog_pt <- matched |>
  select(dog_id, drug, pt) |>
  filter(!is.na(pt)) |>
  mutate(pt = str_to_lower(pt)) |>
  distinct()

llt_subset <- dog_llt |>
  filter(llt %in% llts_of_interest)

calc_dpa_pt <- function(term, df, test_drug = "librela") {
  df_bin <- df |>
    mutate(
      has_pt = pt == term,
      is_test = drug == test_drug
    ) |>
    group_by(dog_id, is_test) |>
    summarise(
      has_pt = any(has_pt),
      .groups = "drop"
    )
  
  a <- sum(df_bin$is_test & df_bin$has_pt, na.rm = TRUE)
  b <- sum(!df_bin$is_test & df_bin$has_pt, na.rm = TRUE)
  n_test <- sum(df_bin$is_test, na.rm = TRUE)
  n_other <- sum(!df_bin$is_test, na.rm = TRUE)
  
  if (a < 5 || b < 5) return(NULL)
  
  prr <- (a / n_test) / (b / n_other)
  
  c <- n_test - a
  d <- n_other - b
  total <- a + b + c + d
  expected <- c(
    (a + c) * (a + b) / total,
    (a + c) * (c + d) / total,
    (b + d) * (a + b) / total,
    (b + d) * (c + d) / total
  )
  observed <- c(a, c, b, d)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  
  tibble(
    pt = str_to_title(term),
    librela_prop = a / n_test,
    other_prop   = b / n_other,
    prr,
    chisq,
    signal = a >= 3 & prr >= 2 & chisq >= 4
  )
}


pt_dpa <- map_dfr(pts_of_interest, ~ calc_dpa_pt(.x, dog_pt))

pt_dpa |>
  mutate(
    librela_prop = scales::percent(librela_prop, accuracy = 0.1),
    other_prop   = scales::percent(other_prop, accuracy = 0.1),
    prr = round(prr, 2),
    chisq = round(chisq, 1),
    signal = if_else(signal, "⚠️", "")
  ) |>
  arrange(desc(prr)) |>
  gt::gt() |>
  gt::tab_header(
    title = "Disproportionality Analysis of Selected PTs",
    subtitle = "Librela vs Other Drugs (Dog-level, PT-based)"
  ) |>
  gt::cols_label(
    pt = "Preferred Term (PT)",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR",
    chisq = "Chi-sq",
    signal = "Signal"
  )
