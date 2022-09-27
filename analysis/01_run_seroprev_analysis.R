# This script estimates seroprevalence
# 
# Preamble ----------------------------------------------------------------

library(tidyverse)
library(rstan)
library(uuid)
library(optparse)

source(here::here("analysis/functions_Stan_multinomial.R"))

set.seed(123)

option_list <- list(
  make_option(c("-i", "--input"), default = here::here("data", "processed", "synthetic_sero_sample.csv"),
              type = "character", help = "Input serology dataset"),
  make_option(c("-h", "--hh_effect"), default = TRUE,
              type = "logical", help = "Whether to use a household random effect or not"),
  make_option(c("-n", "--ID_note"), default = "_HH_multinomial_test",
              type = "character", help = "Note on file names")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Session ID
if (!opt$hh_effect) {
  opt$ID_note <- paste0("_NO", opt$ID_note)
}

# use a session ID for similar filenames from same code run
session_ID <- paste0(substr(uuid::UUIDgenerate(), 1, 8), opt$ID_note)

## Stan control settings that can be changed for testing
# (for production use 4 chains and at least 1500 iter, 250 warmup)
n_chains <- 4
n_iter <- 20
n_warmup <- 10
random_seed <- 1546
options(mc.cores = 4)
p_delta <- 0.99
n_treedepth <- 20
rstan_options(auto_write = TRUE)

# Load data ---------------------------------------------------------------

input_dat <- read_csv(opt$input, 
                      col_types = 
                        cols(uid = col_character(),
                             hh_id = col_character(),
                             sex = col_factor(),
                             age = col_number(),
                             age_cat_10y = col_factor(),
                             age_cat_UEP = col_factor(),
                             S_interp = col_character(),
                             N_interp = col_character(),
                             vaccinated = col_logical()
                        ))

## Define desired age cuts
age_cuts <- c(0, 6, 12, 18, 25, 35, 50, 65, 75, 105)

# Factor references
age_ref <- "[25,35)"
sex_ref <- "Female"

input_dat <- input_dat %>%
  mutate(
    age_cat = factor(cut(age, age_cuts, right = FALSE, include.lowest = TRUE)),
    Sex = fct_recode(sex, "Female" = "f", "Male" = "m"),
    UEP_S_result = fct_recode(
      S_interp,
      "0" = "neg", "1" = "pos"
    ) %>%
      as.character() %>%
      as.integer(),
    UEP_N_result = fct_recode(
      N_interp,
      "0" = "neg", "1" = "pos"
    ) %>%
      as.character() %>%
      as.integer(),
    vaccinated = as.integer(vaccinated)
  ) %>%
  droplevels() %>% 
  mutate(
    Sex = fct_relevel(Sex, ref = "Female"),
    age_cat = fct_relevel(age_cat, ref = age_ref)
  )

# Control data ------------------------------------------------------------
# Roche S data from lab validation and Roche website
# https://diagnostics.roche.com/global/en/products/params/elecsys-anti-sars-cov-2-s.html
# Clinical sensitivity and Clinical specificity
# 1,423 of the tested samples had a sampling date of 14 days or later after diagnosis with PCR. 1,406 of these 1,423 samples were determined with ≥0.8 U/mL in the Elecsys® Anti‑SARS‑CoV‑2 S assay and hence considered positive, resulting in a sensitivity of 98.8 % (95 % CI: 98.1 – 99.3 %) in this sample cohort
# A total of 5,991 samples from diagnostic routine and blood donors drawn before October 2019 were tested with the Elecsys® Anti-SARSCoV-2 S assay. Overall specificity in this cohort of pre-pandemic samples was 99.98 % (95 % CI: 99.91 – 100 %).
## number of positive controls from validation data
pos_control_S = c(172, 1423)
## number of negative controls
neg_control_S = c(185, 5991)
## number of true positives for cases
control_tp_S = c(159, 1406)
## number of false positives for controls (1-specificity)
control_fp_S = c(1, 1)

# Roche N data from lab validation and
# https://www.thelancet.com/action/showPdf?pii=S1473-3099%2820%2930634-4
# Table at bottom of pg 4 "Samples with complete data available taken ≥14 days post symptom onset"
## number of positive controls from validation data
pos_control_N = c(172, 561)
## number of negative controls
neg_control_N = c(185, 976)
## number of true positives for cases
control_tp_N = c(151, 543)
## number of false positives for controls (1-specificity)
control_fp_N = c(0, 2)

# Geneva population ages data ---------------------------------------------------------------

# !!NB: this would need to be changed if use with different population data
population_data <- getGE_age_cats(age_cuts) %>%
  mutate(Sex = fct_relevel(Sex, ref = "Female"),
         age_cat = fct_relevel(age_cat, ref = age_ref)) %>%
  filter(age_cat %in% levels(input_dat$age_cat),
         Sex %in% levels(input_dat$Sex)) %>%
  droplevels()

# Geneva population vaccination data -------------------------------------------
# !!NB: These fractions are based on the Geneva state's statistics. Vaccination 
# fractions would need to be adapted for other settings.
# Repeat for each age class.
GE_vacc <- rep(c(0.0217, 0.0486, 0.5113, 0.6845, 0.7559, 0.8052, 0.8587, 0.8915, 0.9383), 2)

# Save data for later use -------------------------------------------------

seroprev_run_data <- list(
  input_dat = input_dat,
  population_data = population_data,
  GE_vacc = GE_vacc,
  sero_controls = list(
    pos_control_S  = pos_control_S,
    neg_control_S = neg_control_S,
    control_tp_S = control_tp_S,
    control_fp_S = control_fp_S,
    pos_control_N  = pos_control_N,
    neg_control_N = neg_control_N,
    control_tp_N = control_tp_N,
    control_fp_N = control_fp_N
  )
)

saveRDS(seroprev_run_data, 
        file = here::here("data", "processed", paste0(session_ID, "_seroprev_run_data.rds")))

# Where to save output --------------------------------------------------------------
if (!dir.exists("output")) 
  dir.create("output")

if (!dir.exists("output/results")) 
  dir.create("output/results")

output_result_list_filename  <- here::here(
  "data", "processed", paste0(
    session_ID,
    "_results-list_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  ))

# Define regression model -------------------------------------------------

model <- "Sex + age_cat"

# Run models --------------------------------------------------------------

if (!opt$hh_effect) {
  stan_script <- here::here("analysis/stan/seroprev-hh-multinomial_pred.stan")
} else {
  stan_script <- here::here("analysis/stan/seroprev-NOhh-multinomial_pred.stan")
}

if (!opt$hh_effect) {
  ## model the overall seroprevalence with household effect
  fun <- run_multinomial_with_hh_pred
} else {
  ## model the overall seroprevalence (only neutralisation (index) participants - no household effect)
  fun <- run_multinomial_no_hh_pred
}

calc_seropos = fun(
  model_script = stan_script,
  dat = input_dat,
  coef_eqn = model,
  vaccinated = input_dat %>% pull(vaccinated),
  GE_vacc = GE_vacc,
  pos_control_S = pos_control_S,
  neg_control_S = neg_control_S,
  control_tp_S = control_tp_S,
  control_fp_S = control_fp_S,
  pos_control_N = pos_control_N,
  neg_control_N = neg_control_N,
  control_tp_N = control_tp_N,
  control_fp_N = control_fp_N,
  n_cores = getOption("mc.cores"),
  chains = n_chains,
  iter = n_iter,
  warmup = n_warmup,
  Stan_control = list(
    adapt_delta = p_delta,
    max_treedepth = n_treedepth
  ),
  seed = random_seed,
  pop_age_cats = population_data,
  session_ID = session_ID,
  sex_ref = sex_ref,
  age_ref = age_ref
)

# Calculate seroprevalence probabilities --------------------------------------

seroprev_probs <- calc_integrated_seroprev_probs_pred(calc_seropos$stan_posterior,
                                                      coef_eqn = model,
                                                      pop_age_cats = population_data,
                                                      session_ID = session_ID
)

pop_cats_resp_probs <- seroprev_probs$pop_cats_resp_probs
pop_cats_infect_probs <- seroprev_probs$pop_cats_infect_probs
pop_cats_vacc_probs <- seroprev_probs$pop_cats_vacc_probs
pop_cats_any_probs <- seroprev_probs$pop_cats_any_probs

subset_estimates_pospos <- compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S+N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_negpos <- compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S-N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_posneg <- compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S+N-`, age_cat, Sex, strata_pop, sim))
subset_estimates_negneg <- compute_weighted_estimates(pop_cats_resp_probs %>% select(seropos=`S-N-`, age_cat, Sex, strata_pop, sim))

subset_estimates_infect <- compute_weighted_estimates(pop_cats_infect_probs)
subset_estimates_vacc <- compute_weighted_estimates(pop_cats_vacc_probs)
subset_estimates_any <- compute_weighted_estimates(pop_cats_any_probs)

all_results_list <- c(
  model_matrix = list(calc_seropos$model_matrix),
  seroprev_probs,
  subset_estimates_pospos = list(subset_estimates_pospos),
  subset_estimates_negpos = list(subset_estimates_negpos),
  subset_estimates_posneg = list(subset_estimates_posneg),
  subset_estimates_negneg = list(subset_estimates_negneg),
  subset_estimates_infect = list(subset_estimates_infect),
  subset_estimates_vacc = list(subset_estimates_vacc),
  subset_estimates_any = list(subset_estimates_any)
)

saveRDS(all_results_list, output_result_list_filename)

cat("\n-------\n", output_result_list_filename, " saved at ",
    format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n\n", sep = "")

# Using GE vaccination data
seroprev_probs_GE_vacc <- calc_integrated_seroprev_probs_using_GE_vacc(calc_seropos$stan_posterior,
                                                                       coef_eqn = model,
                                                                       pop_age_cats = population_data,
                                                                       session_ID = session_ID
)

pop_cats_resp_probs_GE_vacc <- seroprev_probs_GE_vacc$pop_cats_resp_probs_GE_vacc
pop_cats_infect_probs_GE_vacc <- seroprev_probs_GE_vacc$pop_cats_infect_probs_GE_vacc
pop_cats_any_probs_GE_vacc <- seroprev_probs_GE_vacc$pop_cats_any_probs_GE_vacc

subset_estimates_pospos_GE_vacc <- compute_weighted_estimates(pop_cats_resp_probs_GE_vacc %>% select(seropos=`S+N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_negpos_GE_vacc <- compute_weighted_estimates(pop_cats_resp_probs_GE_vacc %>% select(seropos=`S-N+`, age_cat, Sex, strata_pop, sim))
subset_estimates_posneg_GE_vacc <- compute_weighted_estimates(pop_cats_resp_probs_GE_vacc %>% select(seropos=`S+N-`, age_cat, Sex, strata_pop, sim))
subset_estimates_negneg_GE_vacc <- compute_weighted_estimates(pop_cats_resp_probs_GE_vacc %>% select(seropos=`S-N-`, age_cat, Sex, strata_pop, sim))

subset_estimates_infect_GE_vacc <- compute_weighted_estimates(pop_cats_infect_probs_GE_vacc)
subset_estimates_any_GE_vacc <- compute_weighted_estimates(pop_cats_any_probs_GE_vacc)

all_results_list_GE_vacc <- c(
  model_matrix = list(calc_seropos$model_matrix),
  seroprev_probs_GE_vacc,
  subset_estimates_pospos_GE_vacc = list(subset_estimates_pospos_GE_vacc),
  subset_estimates_negpos_GE_vacc = list(subset_estimates_negpos_GE_vacc),
  subset_estimates_posneg_GE_vacc = list(subset_estimates_posneg_GE_vacc),
  subset_estimates_negneg_GE_vacc = list(subset_estimates_negneg_GE_vacc),
  subset_estimates_infect_GE_vacc = list(subset_estimates_infect_GE_vacc),
  subset_estimates_any_GE_vacc = list(subset_estimates_any_GE_vacc)
)

output_result_list_filename_GE_vacc <- here::here(
  "data", "processed", paste0(
    session_ID,
    "_results-list_GE_vacc_",
    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
    ".rds"
  ))

saveRDS(all_results_list_GE_vacc, output_result_list_filename_GE_vacc)

cat("\n-------\n", output_result_list_filename_GE_vacc, " saved at ",
    format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n\n", sep = "")

# Save final output table ----------------------------------------------------

join_cols <- c("var", "Category", "Obs", "Vaccinated", 
               "S Test positive", "S Test negative", 
               "N Test positive", "N Test negative")

# Table with seroprevalence estimates
final_table <- 
  full_join(
    make_final_table(pop_cats_any_probs, subset_estimates_any, 
                     input_dat, age_ref, sex_ref) %>%
      rename(Any_antibody_response=`Seroprevalence (95% CI)`) %>%
      select(-p, -`Median Seroprevalence (95% HDCI)`),
    make_final_table(pop_cats_infect_probs, subset_estimates_infect,
                     input_dat, age_ref, sex_ref) %>%
      select(-`Median Seroprevalence (95% HDCI)`) %>%
      rename(natural_infection=`Seroprevalence (95% CI)`) %>%
      relocate(p, .after = last_col()),
    by = join_cols) %>%
  full_join(
    make_final_table(pop_cats_vacc_probs, subset_estimates_vacc, 
                     input_dat, age_ref, sex_ref) %>% 
      select(-p, -`Median Seroprevalence (95% HDCI)`),
    by = join_cols) %>%
  rename(estimated_vaccination=`Seroprevalence (95% CI)`) %>%
  relocate(estimated_vaccination, .after = "Vaccinated") %>%
  full_join(
    make_final_table(pop_cats_any_probs_GE_vacc, subset_estimates_any_GE_vacc, 
                     input_dat, age_ref, sex_ref) %>% 
      select(-p, -`Median Seroprevalence (95% HDCI)`) %>%
      rename(Any_antibody_response_GE_vacc=`Seroprevalence (95% CI)`),
    by = join_cols) %>%
  relocate(Any_antibody_response_GE_vacc, .after = Any_antibody_response)

res_table_file <- paste0(
  "output/results/", 
  session_ID, "_results_table_", 
  format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")

write_excel_csv(final_table, res_table_file)

cat("\n---- Done ----\n", res_table_file, " written at ", format(Sys.time(), "%Y-%m-%d-%H:%M:%S"), "\n-------\n", sep = "")

# Tables with relative risks
rrisk_any <- calc_rrisk(all_results_list$subset_estimates_any, 
                        input_dat, age_ref, sex_ref)
rrisk_any %>%
  select(Category, Obs, `Relative risk (95% CI)`, p) %>%
  write_excel_csv(
    here::here("output", "results",
               paste0(session_ID, "_rrisk_any_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")))

# For Geneva vaccination
rrisk_any_GE_vacc <- calc_rrisk(all_results_list_GE_vacc$subset_estimates_any_GE_vacc,
                                input_dat, age_ref, sex_ref)
rrisk_any_GE_vacc %>%
  select(Category, Obs, `Relative risk (95% CI)`, p) %>%
  write_excel_csv(
    here::here("output", "results", 
               paste0(session_ID, "_rrisk_any_GE_vacc_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")))

# For infection
rrisk_infect <- calc_rrisk(all_results_list$subset_estimates_infect, 
                           input_dat, age_ref, sex_ref)
rrisk_infect %>%
  select(Category, Obs, `Relative risk (95% CI)`, p) %>%
  write_excel_csv(
    here::here("output", "results", 
               paste0(session_ID, "_rrisk_infect_", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".csv")))
