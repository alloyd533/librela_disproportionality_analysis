# Load necessary libraries
library(tidyverse)  # Core packages
library(readxl)     # To read Excel files
library(gt)         # To build summary tables
library(scales)     # For percentage formatting

# --------------------- LOAD AND CLEAN DATA ---------------------

# Load the combined cleaned CSV file with all reported adverse events
complete <- read_csv("data/complete_combined_clean.csv", show_col_types = FALSE)

# Create a new column 'polypharmacy' indicating if the dog was on >1 drug at the point of having the ADR
complete <- complete |>
  mutate(
    polypharmacy = Drug |>
      # Drugs separated by commas
      str_split(",") |>
      map_chr(~ {
        # For each comma-separated entry, extract the first word (this is the drug name). Accept different doses of same drug
        drug_names <- str_trim(.x) |> str_extract("^[^\\s]+")
        # If more than one unique drug, flag as polypharmacy
        as.character(if (length(unique(drug_names)) > 1) 1 else 0)
      })
  )

# --------------------- CLEAN VEDDRA TERMINOLOGY ---------------------

# Load VEDDRA terminology to map LLTs to PT, HLT, and Organ System
veddra <- read_excel("data/combined-veddra-list-clinical-terms-reporting-suspected-adverse-events-animals-humans-veterinary-medicinal-products_en.xlsx") |>
  rename_with(~ str_replace_all(.x, " ", "_")) |>       # Replace spaces with underscores in column names
  filter(Current__Term_Type != "H") |>                  # Drop human terms
  mutate(across(where(is.character), as.factor)) |>     # Convert all character fields to factors
  rename(
    llt = `Current_Low_Level_Term_(LLT)`,
    pt = `Current_Preferred_Term_(PT)`,
    hlt = `Current_High_Level_Term_(HLT)`,
    organ = `Current_System_Organ_Class_(SOC)_Term`
  ) |>
  select(organ, hlt, pt, llt) |>                        # Keep only relevant columns
  mutate(llt = str_to_lower(llt))                       # Standardise LLTs to lowercase for joining

# --------------------- EXPAND MULTIPLE REACTIONS PER DOG ---------------------

# Expand rows so each LLT reaction gets its own row per dog
complete_long <- complete |>
  mutate(dog_id = row_number()) |>                      # Assign unique dog ID
  separate_longer_delim(Reaction, delim = ",") |>       # Split multiple reactions into separate rows
  mutate(llt = str_trim(Reaction),                      # Clean and lowercase reaction text
         llt = str_to_lower(llt)) |>
  select(-Reaction)                                     # Drop original Reaction column

# Join with VEDDRA to map each LLT to PT, HLT, and organ system
matched <- complete_long |>
  left_join(veddra, by = "llt")

# --------------------- 1) HLT-LEVEL DISPROPORTIONALITY ---------------------

# Prepare a dog-level dataset of each dog × HLT
dog_hlt <- matched |>
  select(dog_id, drug, hlt) |>
  filter(!is.na(hlt)) |>
  distinct()

# Function to calculate PRR for one HLT
calc_dpa_hlt <- function(hlt_name, df, test_drug = "librela") {
  df_bin <- df |>
    mutate(
      # Create logical for whether HLT/drug matches test drug
      has_hlt = hlt == hlt_name,
      is_test = drug == test_drug
    ) |>
    group_by(dog_id, is_test) |>
    summarise(has_hlt = any(has_hlt), .groups = "drop")
  
  # 2x2 table components for disporportionality
  a <- sum(df_bin$is_test & df_bin$has_hlt)
  b <- sum(!df_bin$is_test & df_bin$has_hlt)
  n_test <- sum(df_bin$is_test)
  n_other <- sum(!df_bin$is_test)
  c <- n_test  - a
  d <- n_other - b
  
  
  if (a < 5 || b < 5) return(NULL) # Drop sparse terms
  
  # Propoortional Reporting Ratio
  prr <- (a / n_test) / (b / n_other)
  chiSq = ((a*d - b*c)^2) * (n_test+n_other) / ((a+b)*(c+d)*(a+c)*(b+d))
  
  tibble(
    hlt = hlt_name,
    librela_prop = a / n_test,
    other_prop = b / n_other,
    prr = prr,
    chi = chiSq
  )
}

# Run disproportionality analysis across all HLTs seen in librela reports
hlt_list <- dog_hlt |>
  filter(drug == "librela") |>
  distinct(hlt) |>
  pull(hlt)

hlt_dpa <- map_dfr(hlt_list, ~ calc_dpa_hlt(.x, dog_hlt))

# Display HLT table
hlt_dpa |>
  arrange(desc(prr)) |>
  mutate(
    prr = round(prr, 2),
    librela_prop = percent(librela_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "High-Level Term (HLT) Disproportionality: Librela vs Others",
    subtitle = "Dog-level proportions and PRRs"
  ) |>
  cols_label(
    hlt = "HLT",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR"
  ) |>
  tab_style(
    style = cell_fill(color = "#FFDFDF"),
    locations = cells_body(rows = Highlight)
  ) |>
  cols_hide(columns = Highlight)

# --------------------- 2) ORGAN SYSTEM DISPROPORTIONALITY ---------------------

# Build wide-format matrix of each dog × organ system (binary 0/1)
organ_matrix <- matched |>
  filter(!is.na(organ)) |>
  distinct(dog_id, drug, organ, polypharmacy) |>
  mutate(value = 1) |>
  pivot_wider(names_from = organ, values_from = value, values_fill = 0)

# Function to calculate PRR for one organ system (similar as HLT one above)
calc_prop_dpa <- function(organ_col, df, test_drug = "librela") {
  df_bin <- df |>
    mutate(
      has_reaction = .data[[organ_col]] == 1,
      is_test = drug == test_drug
    )
  
  a <- sum(df_bin$is_test & df_bin$has_reaction)
  b <- sum(!df_bin$is_test & df_bin$has_reaction)
  n_test <- sum(df_bin$is_test)
  n_other <- sum(!df_bin$is_test)
  
  if (a < 10 || b < 10) return(NULL)
  
  prr <- (a / n_test) / (b / n_other)
  
  tibble(
    organ = organ_col,
    librela_prop = a / n_test,
    other_prop = b / n_other,
    prr = prr
  )
}

# Run disproportionality analysis for each organ system
organ_cols <- setdiff(names(organ_matrix), c("dog_id", "drug", "polypharmacy"))
organ_dpa <- map_dfr(organ_cols, ~ calc_prop_dpa(.x, organ_matrix))

# Display organ system results
organ_dpa |>
  arrange(desc(prr)) |>
  mutate(
    prr = round(prr, 2),
    librela_prop = percent(librela_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "Organ System Disproportionality: Librela vs Others",
    subtitle = "Proportion of dogs with any reaction in each system"
  ) |>
  cols_label(
    organ = "Organ System",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR"
  ) |>
  tab_style(
    style = cell_fill(color = "#FFDFDF"),
    locations = cells_body(rows = Highlight)
  ) |>
  cols_hide(columns = Highlight)

# --------------------- PT TERMS OF INTEREST ---------------------

# Define a list of key PTs based on FDA concerns
pts_of_interest <- c(
  "Ataxia", "Convulsion", "Paresis", "Proprioception abnormality",
  "Paralysis", "Recumbency", "Muscle weakness", "Muscle tremors",
  "Lameness", "Collapse NOS",
  "Pancreatitis", "Death",
  "Immune mediated haemolytic anaemia", 
  "Immune mediated thrombocytopenia", 
  "Immune mediated polyarthritis"
) |> str_to_lower()

# Build PT-level dataset per dog
dog_pt <- matched |>
  select(dog_id, drug, pt) |>
  filter(!is.na(pt)) |>
  mutate(pt = str_to_lower(pt)) |>
  distinct()

# Function to calculate PRR for individual PT
calc_dpa_pt <- function(term, df, test_drug = "librela") {
  df_bin <- df |>
    mutate(
      has_pt = pt == term,
      is_test = drug == test_drug
    ) |>
    group_by(dog_id, is_test) |>
    summarise(has_pt = any(has_pt), .groups = "drop")
  
  a <- sum(df_bin$is_test & df_bin$has_pt)
  b <- sum(!df_bin$is_test & df_bin$has_pt)
  n_test <- sum(df_bin$is_test)
  n_other <- sum(!df_bin$is_test)
  
  if (a < 10 || b < 10) return(NULL)
  
  prr <- (a / n_test) / (b / n_other)
  
  tibble(
    pt = str_to_title(term),
    librela_prop = a / n_test,
    other_prop = b / n_other,
    prr = prr
  )
}

# Run PT DPA
pt_dpa <- map_dfr(pts_of_interest, ~ calc_dpa_pt(.x, dog_pt))

pt_dpa |>
  arrange(desc(prr)) |>
  mutate(
    prr = round(prr, 2),
    librela_prop = percent(librela_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "Disproportionality Analysis of Selected Preferred Terms",
    subtitle = "Librela vs Other Drugs (Dog-level)"
  ) |>
  cols_label(
    pt = "Preferred Term (PT)",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR"
  ) |>
  tab_style(
    style = cell_fill(color = "#FFDFDF"),
    locations = cells_body(rows = Highlight)
  ) |>
  cols_hide(columns = Highlight)

# --------------------- SENSITIVITY: NO POLYPHARMACY ---------------------

# Subset organ matrix to dogs without polypharmacy
organ_matrix_nopoly <- organ_matrix |>
  filter(polypharmacy == "0")

# Re-run organ DPA on this subset
organ_dpa_nopoly <- map_dfr(organ_cols, ~ calc_prop_dpa(.x, organ_matrix_nopoly))

organ_dpa_nopoly |>
  arrange(desc(prr)) |>
  mutate(
    prr = round(prr, 2),
    librela_prop = percent(librela_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "Organ System Disproportionality (Single agent)",
    subtitle = "Subset analysis restricted to monotherapy reports"
  ) |>
  cols_label(
    organ = "Organ System",
    librela_prop = "Librela %",
    other_prop = "Other Drugs %",
    prr = "PRR"
  ) |>
  tab_style(
    style = cell_fill(color = "#FFDFDF"),
    locations = cells_body(rows = Highlight)
  ) |>
  cols_hide(columns = Highlight)

## assumes: tidyverse loaded; objects `matched` (your join result) available
## matched has columns: dog_id (report proxy), drug, pt, hlt, organ
## 1) FDA-style exclusions -----------------------------------------------

excl_pts   <- c("Emesis")
excl_hlt   <- c("Lack of efficacy")
excl_socs  <- c("Product defects","Med error","Med dev","IGA Animal","Other events")

pt_df <- matched |>
  filter(!is.na(pt)) |>
  mutate(
    pt   = str_to_title(pt),
    hlt  = str_to_title(hlt),
    organ= str_to_title(organ)
  ) |>
  filter(!(pt %in% excl_pts),
         !(hlt %in% excl_hlt),
         !(organ %in% str_to_title(excl_socs))) |>
  distinct(dog_id, drug, pt)        # report × drug × PT presence

## 2) Helper to compute 2×2, PRR, chi-square, EBGM/EB05 ------------------

dpa_one_pt <- function(df, term, test_drug = "librela") {
  x <- df |>
    mutate(is_test = drug == test_drug,
           has_pt  = pt == term) |>
    group_by(dog_id, is_test) |>
    summarise(has_pt = any(has_pt), .groups = "drop")
  a <- sum(x$is_test &  x$has_pt)
  c <- sum(x$is_test & !x$has_pt)
  b <- sum(!x$is_test &  x$has_pt)
  d <- sum(!x$is_test & !x$has_pt)
  N <- a + b + c + d
  if (N == 0) return(NULL)
  
  # PRR and chi-square (without Yates, as per many PV implementations)
  prr <- (a / (a + c)) / (b / (b + d))
  chi <- ( (a * d - b * c)^2 ) * N / ((a + b) * (c + d) * (a + c) * (b + d))
  
  # MGPS-style EBGM with simple Gamma-Poisson shrinker
  # Expected under independence:
  E   <- (a + c) * (a + b) / N
  # Prior hyperparams (weakly informative). You can tune (alpha, beta); alpha=1,beta=1 is common.
  alpha <- 1; beta <- 1
  post_alpha <- alpha + a
  post_beta  <- beta  + E
  ebgm <- post_alpha / post_beta
  # EB05 is lower 90% bound of (theta/E); theta ~ Gamma(post_alpha, post_beta)
  eb05 <- qgamma(0.05, shape = post_alpha, rate = post_beta) / E
  
  tibble(
    pt = term, A = a, B = b, C = c, D = d, N = N,
    PRR = prr, ChiSq = chi, EBGM = ebgm, EB05 = eb05
  )
}

## 3) Run across all PTs present for test drug ----------------------------

pts <- pt_df |> filter(str_to_lower(drug) == "librela") |> distinct(pt) |> pull(pt)

pt_dpa <- map_dfr(pts, ~ dpa_one_pt(pt_df, .x, test_drug = "librela")) |>
  mutate(
    signal_PRR = (A >= 3) & (PRR >= 2) & (ChiSq >= 4),
    signal_EB  = EB05 >= 2,
    SIGNAL     = signal_PRR | signal_EB
  ) |>
  arrange(desc(SIGNAL), desc(PRR))

## 4) View signals (replicating FDA “any of the algorithms” rule) --------

pt_dpa |> filter(SIGNAL)
