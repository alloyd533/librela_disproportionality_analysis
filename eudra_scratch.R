## --- Packages ---
stopifnot(
  require(tidyverse),
  require(readxl),
  require(gt),
  require(scales),
  require(openEBGM)
)

## --- LOAD & POLYPHARMACY (unchanged but tightened) ---
complete <- read_csv("data/complete_combined_clean.csv", show_col_types = FALSE)

first_token <- function(x) str_extract(str_trim(x), "^[^\\s]+")
complete <- complete |>
  mutate(polypharmacy = Drug |>
           str_split(",") |>
           map_chr(~ as.integer(length(unique(first_token(.x))) > 1) |> as.character()))

## --- VEDDRA (tidy select/rename) ---
veddra <- read_excel(
  "data/combined-veddra-list-clinical-terms-reporting-suspected-adverse-events-animals-humans-veterinary-medicinal-products_en.xlsx"
) |>
  filter(`Current  Term Type` != "H") |>
  transmute(
    organ = `Current System Organ Class (SOC) Term`,
    hlt   = `Current High Level Term (HLT)`,
    pt    = `Current Preferred Term (PT)`,
    llt   = str_to_lower(`Current Low Level Term (LLT)`)
  )

## --- EXPAND REACTIONS & JOIN ---
matched <- complete |>
  mutate(dog_id = row_number()) |>
  tidyr::separate_longer_delim(Reaction, ",") |>
  mutate(llt = Reaction |> str_trim() |> str_to_lower()) |>
  select(-Reaction) |>
  left_join(veddra, by = "llt")

## --- Presence builders ---
presence_by <- function(df, level = c("pt","hlt","organ")) {
  level <- match.arg(level)
  df |>
    select(dog_id, drug, !!sym(level)) |>
    filter(!is.na(.data[[level]])) |>
    mutate(!!level := str_to_title(.data[[level]])) |>
    distinct()
}
dog_pt  <- presence_by(matched, "pt")
dog_hlt <- presence_by(matched, "hlt")
dog_org <- presence_by(matched, "organ")

## --- Core 2x2 + PRR + ChiSq ---
dpa_one_term <- function(df, term_col, term_val, test_drug = "librela") {
  x <- df |>
    mutate(is_test = str_to_lower(drug) == str_to_lower(test_drug),
           has_trm = .data[[term_col]] == term_val) |>
    group_by(dog_id, is_test) |>
    summarise(has_trm = any(has_trm), .groups = "drop")
  
  a <- sum( x$is_test &  x$has_trm)
  c <- sum( x$is_test & !x$has_trm)
  b <- sum(!x$is_test &  x$has_trm)
  d <- sum(!x$is_test & !x$has_trm)
  N <- a + b + c + d
  if (N == 0) return(NULL)
  
  prr <- (a/(a+c)) / (b/(b+d))
  chi <- ((a*d - b*c)^2) * N / ((a+b)*(c+d)*(a+c)*(b+d))
  
  tibble(term = term_val, A = a, B = b, C = c, D = d, N = N, PRR = prr, ChiSq = chi)
}

run_dpa_base <- function(df_level, term_col = "pt", test_drug = "librela") {
  terms <- df_level |>
    filter(str_to_lower(drug) == str_to_lower(test_drug)) |>
    distinct(.data[[term_col]]) |>
    pull(1)
  map_dfr(terms, ~ dpa_one_term(df_level, term_col, .x, test_drug))
}

## --- EBGM/EB05 with openEBGM (no squashing for small data) ---
ebgm_for_level <- function(df_level, term_col = "pt") {
  pairs <- df_level |>
    transmute(id = dog_id, var1 = str_to_title(drug), var2 = str_to_title(.data[[term_col]]))
  
  proc <- processRaw(pairs)
  
  set.seed(1)
  theta_init <- data.frame(
    alpha1 = c(0.5, 1), beta1 = c(0.5, 1),
    alpha2 = c(2, 3),   beta2 = c(2, 3),
    p      = c(0.1, 0.2)
  )
  
  hyp   <- autoHyper(data = proc, theta_init = theta_init, squashed = FALSE)
  scored <- ebScores(processed = proc, hyper_estimate = hyp, quantiles = c(5, 95))
  
  scored$data |>
    transmute(drug = var1,
              term = var2, 
              EBGM = EBGM,
              EB05 = QUANT_05, 
              EB95 = QUANT_95)
}

## --- PT analysis (repeat for HLT/Organ if needed) ---
base_pt <- run_dpa_base(dog_pt, "pt", "librela")
eb_pt   <- ebgm_for_level(dog_pt, "pt") %>%
  filter(drug == "Librela")


pt_dpa <- base_pt |>
  left_join(eb_pt, by = "term") |>
  mutate(
    signal_PRR = (A >= 3) | (PRR >= 2) | (ChiSq >= 4),
    signal_EB  = EB05 >= 2,
    SIGNAL     = coalesce(signal_PRR, FALSE) | coalesce(signal_EB, FALSE)
  ) |>
  arrange(desc(SIGNAL), desc(PRR))

## --- HLT / Organ (optional) ---
# base_hlt <- run_dpa_base(dog_hlt, "hlt", "librela")
# eb_hlt   <- ebgm_for_level(dog_hlt, "hlt")
# hlt_dpa  <- base_hlt |> left_join(eb_hlt, by = "term") |>
#   mutate(signal_PRR = (A>=3)&(PRR>=2)&(ChiSq>=4),
#          signal_EB = EB05>=2,
#          SIGNAL = coalesce(signal_PRR,FALSE)|coalesce(signal_EB,FALSE)) |>
#   arrange(desc(SIGNAL), desc(PRR))

## --- Pretty table (PT) ---
pt_dpa |>
  mutate(across(c(PRR, EBGM, EB05, EB95), ~ round(.x, 2))) |>
  filter(PRR >1, is.finite(PRR), SIGNAL == "TRUE") |>
  arrange(desc(PRR)) |>
  gt() |>
  tab_header(title = "PT Disproportionality: Librela vs All Other Drugs") |>
  cols_label(term = "PT", ChiSq = "Chi²") |>
  data_color(columns = SIGNAL, colors = col_factor(c("#FFDFDF","#FFFFFF"), domain = c(TRUE, FALSE)))

## --- Monotherapy sensitivity (polypharmacy == 0) ---
presence_by_nopoly <- function(level) matched |> filter(polypharmacy == "0") |> presence_by(level)
pt_dpa_nopoly <- {
  base_np <- run_dpa_base(presence_by_nopoly("pt"), "pt", "librela")
  eb_np   <- ebgm_for_level(presence_by_nopoly("pt"), "pt")
  base_np |>
    left_join(eb_np, by = "term") |>
    mutate(
      signal_PRR = (A >= 3) & (PRR >= 2) & (ChiSq >= 4),
      signal_EB  = EB05 >= 2,
      SIGNAL     = signal_PRR | signal_EB
    )
}

### checker
## assumes packages already loaded and `matched` exists

## 1) Build report-level presence at PT level (no functions)
dog_hlt <- matched |>
  select(dog_id, drug, hlt) |>
  filter(!is.na(hlt)) |>
  mutate(
    drug = str_to_title(drug),
    hlt   = str_to_title(hlt)
  ) |>
  distinct()

test_drug <- "Librela"

## 2) 2x2 counts for Librela vs pooled others
n_test  <- dog_hlt |> filter(drug == test_drug) |> distinct(dog_id) |> nrow()
n_other <- dog_hlt |> filter(drug != test_drug) |> distinct(dog_id) |> nrow()

pt_counts <- dog_hlt |>
  group_by(hlt) |>
  summarise(
    A = n_distinct(dog_id[drug == test_drug]),
    B = n_distinct(dog_id[drug != test_drug]),
    .groups = "drop"
  ) |>
  mutate(
    C = n_test  - A,
    D = n_other - B,
    N = A + B + C + D,
    PRR = (A / (A + C)) / (B / (B + D)),
    ChiSq = ((A*D - B*C)^2) * N / ((A+B)*(C+D)*(A+C)*(B+D))
  )

## 3) EBGM / EB05 with openEBGM (fit on ALL drugs, keep Librela rows only)
pt_pairs <- dog_hlt |>
  transmute(id = dog_id, var1 = drug, var2 = hlt)

proc_pt <- processRaw(pt_pairs)

set.seed(1)
theta_init <- data.frame(
  alpha1 = c(0.5, 1), beta1 = c(0.5, 1),
  alpha2 = c(2, 3),   beta2 = c(2, 3),
  p      = c(0.05, 0.15)
)

hyp_pt   <- autoHyper(data = proc_pt, theta_init = theta_init, squashed = FALSE)
scored_pt <- ebScores(processed = proc_pt, hyper_estimate = hyp_pt, quantiles = c(5, 95))

# column names for quantiles vary by version; pick what exists
q5_name  <- intersect(c("EB05","Q_5","QUANT_05"),  names(scored_pt$data))[1]
q95_name <- intersect(c("EB95","Q_95","QUANT_95"), names(scored_pt$data))[1]

eb_pt <- scored_pt$data |>
  filter(var1 == test_drug) |>
  transmute(
    pt   = var2,
    EBGM = EBGM,
    EB05 = .data[[q5_name]],
    EB95 = .data[[q95_name]]
  )

## 4) Merge and flag signals (Librela-only)
pt_dpa <- pt_counts |>
  left_join(eb_pt, by = "pt") |>
  mutate(
    signal_PRR = (A >= 3) & (PRR >= 2) & (ChiSq >= 4),
    signal_EB  = EB05 >= 2,
    SIGNAL     = coalesce(signal_PRR, FALSE) | coalesce(signal_EB, FALSE)
  ) |>
  arrange(desc(SIGNAL), desc(PRR))

## 5) Display
pt_dpa |>
  filter(B > 5, A > 5) |>
  arrange(desc(PRR)) |>
  mutate(across(c(PRR, EBGM, EB05, EB95), ~ round(.x, 2))) |>
  gt() |>
  tab_header(title = "PT Disproportionality: Librela vs All Other Drugs") |>
  cols_label(pt = "PT", ChiSq = "Chi²")

# Build the report-level input pvda expects (one row per report–drug–event)
pvda_pt <- matched |>
  select(report_id = dog_id, drug, event = pt) |>
  filter(!is.na(event)) |>
  mutate(drug = str_to_title(drug),
         event = str_to_title(event)) |>
  distinct()

# Run BCPNN/IC (da() computes obs & expected, and IC with shrinkage = 0.5)
da_pt <- da(pvda_pt)                     # default df_colnames already match

# Extract Librela vs all others, and flag on IC025 > 0
ic_pt_librela <- da_pt$da_df |>
  filter(drug == "Librela") |>
  transmute(drug, event, obs, exp_rrr, ic2.5, ic, ic97.5,
            SIGNAL = ic2.5 > 0)
