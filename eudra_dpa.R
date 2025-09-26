# Load necessary libraries
library(tidyverse)  # Core packages
library(readxl)     # To read Excel files
library(gt)         # To build summary tables
library(scales)     # For percentage formatting
require(openEBGM)   # For MGPS

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

# Join with VEDDRA to map each LLT to PT, HLT, and organ system then clean
matched <- complete_long |>
  left_join(veddra, by = "llt") |>
  rename(drug_hx = Drug,
         sex = Gender,
         death = Animals_died) |>
  mutate(across(where(is.character) & !any_of("drug_hx"), as.factor),
         year = factor(year),
         death = factor(death),
         sex = case_when(
           sex == "FEMALE" ~ "f",
           sex == "MALE" ~ "m",
           TRUE ~ "unknown"
         )) |>
  select(-Animals_affected) |>
  rename_with(tolower)
          
write_rds(matched, "data/matched")

# --------------------- 1) PT-LEVEL DISPROPORTIONALITY ---------------------

# Prepare a dog-level dataset of each dog × PT
dog_pt <- matched |>
  filter(!is.na(pt)) |>
  arrange(dog_id, drug, pt) %>%  
  distinct(dog_id, drug, pt, .keep_all = TRUE)

test_drug <- "bedinvetmab"

denoms <- dog_pt %>%
  mutate(is_test = drug == test_drug) %>%
  distinct(dog_id, is_test) %>%
  summarise(n_test = sum(is_test), n_other = sum(!is_test), .groups = "drop")

dpa_pt <- dog_pt %>%
  mutate(is_test = drug == test_drug) %>%
  distinct(dog_id, is_test, pt) %>%                 # presence of PT within test/other per dog
  group_by(pt) %>%
  summarise(a = sum(is_test),                       # test drug + PT
            b = sum(!is_test),                      # other drugs + PT
            .groups = "drop") %>%
  mutate(across(c(a,b), as.numeric),
         n_test = denoms$n_test,
         n_other = denoms$n_other,
         c = n_test - a,
         d = n_other - b,
         n_tot = n_test + n_other,
         bedinvetmab_prop = a / n_test, 
         other_prop = b / n_other,
         prr = (a / n_test) / (b / n_other),
         chi = ((a*d - b*c)^2) * n_tot / ((a+b)*(c+d)*(a+c)*(b+d))) %>%
  filter(a >= 3, b >= 3, is.finite(prr), is.finite(chi)) %>%
  arrange(desc(prr))

# Display pt table
dpa_pt |>
  arrange(desc(prr)) |>
  mutate(
    prr = round(prr, 2),
    chi = round(chi, 2),
    bedinvetmab_prop = percent(bedinvetmab_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2 & chi >=4,
  ) |>
  select(pt, a, bedinvetmab_prop, b, other_prop, prr, chi, Highlight) |>
  gt() |>
  tab_header(
    title = "Preferred Term (pt) Disproportionality: Bedinvetmab vs Comparator drugs",
    subtitle = "Dog-level proportions and PRRs"
  ) |>
  cols_label(
    pt = "Preferred Term",
    bedinvetmab_prop = "Bedinvetmab %",
    a = "Bedinvetmab Reports",
    b = "Comparator Reports",
    other_prop = "Comparator Drugs %",
    prr = "PRR",
    chi = "Chi-Square Value"
  ) |>
  tab_style(
    style = cell_fill(color = "#FFDFDF"),
    locations = cells_body(rows = Highlight)
  ) |>
  cols_hide(columns = Highlight)

# --------------------- 2) Shrinkage methods -----------------------------------
## ------------------------------------------------------------
## 1) MGPS (openEBGM) — EBGM / EB05 (stratified)
## ------------------------------------------------------------

dat_mgps <- dog_pt |>
  mutate(year = case_when(
         year %in% 2004:2014 ~ "2004–2014",
         year %in% 2015:2021 ~ "2015–2021",
         year == 2022        ~ "2022",
         year == 2023        ~ "2023",
         year == 2024        ~ "2024",
         year == 2025        ~ "2025",
         TRUE                ~ NA_character_
  )) |>
  mutate(year = factor(
    year,
    levels = c("2004–2014","2015–2021","2022","2023","2024","2025")
  )) |>
  transmute(
    id     = dog_id,
    var1   = drug,
    var2   = pt)

,
    # strat_year     = as.character(year),
    # strat_region  = countrycode::countrycode(country, "country.name", "continent"),
    # strat_sex    = toupper(sex),
    # strat_reporter = tolower(reporter)
  )

# Process raw (stratify is logical; strata detected by substring 'strat')
proc <- processRaw(
  data     = dat_mgps,
  # stratify = TRUE,   # uses all columns containing 'strat'
  zeroes   = FALSE
)

set.seed(1)

# startpoints for the initial hyperparameter guesses 
theta_init <- data.frame(
  alpha1 = c(0.5, 1), 
  beta1 = c(0.5, 1),
  alpha2 = c(2, 3),  
  beta2 = c(2, 3),
  p      = c(0.1, 0.2)
)

# estimate hyperparameters
hyp  <- autoHyper(data = proc, 
                  theta_init = theta_init, 
                  squashed = FALSE)    # estimates the 2-gamma mixture prior

# EBGM + quantiles (returns an openEBGM object; results in $data)
obj  <- ebScores(processed = proc, hyper_estimate = hyp, quantiles = c(5, 95))
res  <- obj$data

# Pull Librela PT signals (common gate: A≥3 & EB05≥2)
out <- res |>
  filter(var1 == "bedinvetmab") |>
  transmute(
    pt = var2, N, E, RR, PRR, EBGM, EB05 = QUANT_05, EB95 = QUANT_95
  ) |>
  arrange(desc(EB05))

sig_mgps <- filter(out, N >= 3, EB05 >= 2)

## ------------------------------------------------------------
## 2) BCPNN — EBGM / EB05 (stratified)
## ------------------------------------------------------------
 

# pooled 2×2 vs all other drugs, per PT
n_tot <- nrow(dog_pt)
n_lib <- sum(dog_pt$drug == "bedinvetmab")

pt_tot <- dog_pt |> count(pt, name = "np1")
pt_lib <- dog_pt |> 
  filter(drug == "bedinvetmab") |>
  count(pt, name = "n11")

bc <- pt_tot |>
  left_join(pt_lib, by = "pt") |>
  mutate(
    n11 = replace_na(n11, 0L),
    n1p = n_lib,
    n10 = n1p - n11,
    n01 = np1 - n11,
    n00 = n_tot - n11 - n10 - n01
  )

# helper: Dirichlet via Gamma (no extra pkg)
rdir <- function(n, alpha) {
  X <- matrix(rgamma(n*length(alpha), shape=alpha), nrow=n)
  X / rowSums(X)
}

# IC & 95% floor with Jeffreys prior Dirichlet(½,½,½,½)

bc_ic <- function(a,b,c,d, S=20000L){
  P  <- rdir(S, c(a,b,c,d) + 0.5)      # Jeffreys prior
  p11 <- P[,1]; p10 <- P[,2]; p01 <- P[,3]
  pd  <- p11 + p10; pe <- p11 + p01
  IC  <- log2(p11 / (pd * pe))
  tibble(
    IC       = mean(IC),
    IC_floor = quantile(IC, 0.025, na.rm = TRUE)
  )
}

bcpnn <- bc |>
  rowwise() |>
  mutate(stats = list(bc_ic(n11, n10, n01, n00))) |>
  unnest_wider(stats) |>
  ungroup() 

sig_bcpnn <- filter(bcpnn, n11 >= 3, IC_floor >= 0)

# --------------------- 2) ORGAN SYSTEM DISPROPORTIONALITY ---------------------

# Build wide-format matrix of each dog × organ system (binary 0/1)
organ_matrix <- matched |>
  filter(!is.na(organ)) |>
  distinct(dog_id, drug, organ, polypharmacy) |>
  mutate(value = 1) |>
  pivot_wider(names_from = organ, values_from = value, values_fill = 0)

# Function to calculate PRR for one organ system (similar as pt one above)
calc_prop_dpa <- function(organ_col, df, test_drug = "bedinvetmab") {
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
    bedinvetmab_prop = a / n_test,
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
    bedinvetmab_prop = percent(bedinvetmab_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "Organ System Disproportionality: bedinvetmab vs Others",
    subtitle = "Proportion of dogs with any reaction in each system"
  ) |>
  cols_label(
    organ = "Organ System",
    bedinvetmab_prop = "bedinvetmab %",
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
calc_dpa_pt <- function(term, df, test_drug = "bedinvetmab") {
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
    bedinvetmab_prop = a / n_test,
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
    bedinvetmab_prop = percent(bedinvetmab_prop, 0.1),
    other_prop = percent(other_prop, 0.1),
    Highlight = prr >= 2
  ) |>
  gt() |>
  tab_header(
    title = "Disproportionality Analysis of Selected Preferred Terms",
    subtitle = "bedinvetmab vs Other Drugs (Dog-level)"
  ) |>
  cols_label(
    pt = "Preferred Term (PT)",
    bedinvetmab_prop = "bedinvetmab %",
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
    bedinvetmab_prop = percent(bedinvetmab_prop, 0.1),
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
    bedinvetmab_prop = "bedinvetmab %",
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

dpa_one_pt <- function(df, term, test_drug = "bedinvetmab") {
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

pts <- pt_df |> filter(str_to_lower(drug) == "bedinvetmab") |> distinct(pt) |> pull(pt)

pt_dpa <- map_dfr(pts, ~ dpa_one_pt(pt_df, .x, test_drug = "bedinvetmab")) |>
  mutate(
    signal_PRR = (A >= 3) & (PRR >= 2) & (ChiSq >= 4),
    signal_EB  = EB05 >= 2,
    SIGNAL     = signal_PRR | signal_EB
  ) |>
  arrange(desc(SIGNAL), desc(PRR))

## 4) View signals (replicating FDA “any of the algorithms” rule) --------

pt_dpa |> filter(SIGNAL)
