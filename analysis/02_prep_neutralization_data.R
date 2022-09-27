# This script prepares the neutralization data for modeling

# Preamble ----------------------------------------------------------------

library(tidyverse)
library(optparse)
library(cmdstanr)

source(here::here("analysis", "utils.R"))

option_list <- list(
  make_option(c("-v", "--variant"), default = "var-a", action ="store",
              type = "character", help = "Variant to model"),
  make_option(c("-a", "--age_cat"), default = "age_cat_uep", action ="store",
              type = "character", help = "age categories to use"),
  make_option(c("-m", "--model"), default = "additive", action ="store",
              type = "character", help = "Model type"),
  make_option(c("-s", "--suffix"), default = NULL, action ="store",
              type = "character", help = "Suffix in names"),
  make_option(c("-i", "--input"), default = here::here("data", "processed", "synthetic_neut_sample.csv"),
              type = "character", help = "Input raw neutralization dataset"),
  make_option(c("-x", "--full_data"), default = here::here("data", "processed", "synthetic_true_sample.rds"),
              type = "character", help = "Full dataset with vaccination/infection information")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load data ---------------------------------------------------------------

# Neutralization data
neut_data <- read_csv(opt$input,
                      col_types = 
                        cols(uid = col_character(),
                             sex = col_factor(),
                             age = col_number(),
                             age_cat_10y = col_factor(),
                             age_cat_UEP = col_factor(),
                             variant = col_character(),
                             ec50 = col_number(),
                             inf_latest = col_character(),
                             vacc_status_3way = col_character()
                        )) %>%
  janitor::clean_names() %>%
  mutate(
    uid = as.character(uid),
    # EC50 threshold
    ec50_thresh = ec50 > 50,
    inf_latest = factor(inf_latest) %>% forcats::fct_relevel(ref = "inf_bef_2022"),
    variant = flattenStr(variant) %>% factor(),
    # Relevel factors
    sex = factor(sex) %>% forcats::fct_relevel(ref = "f"),
    age_cat_uep = forcats::fct_relevel(age_cat_uep, ref = "[35,50)"),
    age_cat_10y = forcats::fct_relevel(age_cat_10y, ref = "[20,30)"),
    vacc_status_3way = factor(vacc_status_3way) %>% forcats::fct_relevel(ref = "Dose 1 or 2 only"),
  ) %>%
  # Select age cat for analysis
  rename(age_cat = !!rlang::sym(opt$age_cat))

if (!(opt$variant %in% unique(neut_data$variant))) {
  stop("Uknown variant name '", opt$variant, "', provide valid variant among: ",
       str_c(unique(neut_data$variant), collapse = ", "))
}

# Save all neutralization data
saveRDS(neut_data, file = here::here("data", "processed", "neutralization_data_all.rds"))

# Keep data only for variant to model
neut_data <- filter(neut_data, variant == opt$variant) %>%
  mutate(obs_id = row_number())

# Save for further use
saveRDS(neut_data, file = makeProcessedNeutDataName(opt))

# Regression matrix -------------------------------------------------------

# Model formula for regression
if (opt$model == "additive") {
  model_frml <- as.formula("~ age_cat + sex + inf_latest + vacc_status_3way")
} else if (opt$model == "interaction") {
  model_frml <- as.formula("~ age_cat + sex + inf_latest * vacc_status_3way")
} else {
  stop("Unknown model type.")
}

# Build regression matrix for participants for which we are confident about
# the self-reported infections (excluding N+ & no reported infection)
design_mat_full <- model.matrix(
  model_frml,
  data = neut_data %>% filter(!is.na(inf_latest))
)

# Regression matrix for participants for which we do not have infection dates
# and have N+ serology. Use infection by other variant than omicron as default.
# Uncertainty in the model will be accounted for through marginalization.
design_mat_missing <- model.matrix(
  model_frml,
  data = neut_data %>%
    filter(is.na(inf_latest)) %>%
    # Set levels to match other
    mutate(inf_latest = "inf_bef_2022",
           inf_latest = factor(inf_latest, levels = levels(neut_data$inf_latest)) %>%
             forcats::fct_relevel(ref = "uninfected"))
)

if (ncol(design_mat_full) != ncol(design_mat_missing))
  stop("Dimensions of design matrices not equal.")


# Prior probabilities of non-omicron infections ---------------------------

# !! This is arbitrary. In our analysis we used differences in cumulative
# seroincidence between our seroprevalence surveys. Values used in the analysis
# were:
# c(0.417, 0.478, 0.549, 0.598, 0.523, 0.546, 0.553, 0.528, 0.482)
prior_other_infection <- rep(.45, length(levels(neut_data$age_cat)))

# Columns for omicron/other infection classes
col_omicron <- grep("inf_latestinf_2022", colnames(design_mat_missing)) %>% as.array()
col_other <- grep("inf_latestinf_bef_2022", colnames(design_mat_missing)) %>% as.array()

# Priors ------------------------------------------------------------------

priors <- list(
  mu_alpha = 0,    # Mean of the prior of the intercept
  sd_alpha = 1,    # sd of the prior of the intercept
  mu_alpha_reg = 5,    # Mean of the prior of the intercept for regression
  sd_alpha_reg = 1,    # sd of the prior of the intercept for regression
  sd_beta = 1    # sd of the other regression coefficients
)

# Stan input data ---------------------------------------------------------

stan_data <- list(
  N_survey = nrow(neut_data),      # total number of participants
  N_missing = nrow(design_mat_missing),  # number of participants with N+ but missing infection date
  N_age_cats = length(levels(neut_data$age_cat)),
  N_col_omicron = length(col_omicron),
  ind_full = neut_data %>% filter(!is.na(inf_latest)) %>% pull(obs_id),  # indices of full observations
  ind_missing = neut_data %>% filter(is.na(inf_latest)) %>% pull(obs_id),  # indices of missing infection dates
  y = neut_data$ec50_thresh,       # observations of neutralization above threshold
  p_vars = ncol(design_mat_full),  # number of regression parameter
  X = design_mat_full,             # covariate model matrix for participants with full info
  X_missing = design_mat_missing,   # model matrix for participants with missing infection dates
  col_omicron = col_omicron,
  col_other =  col_other,
  map_missing_age_cat = neut_data %>% filter(is.na(inf_latest)) %>% pull(age_cat) %>% as.numeric(),
  prior_other = prior_other_infection,    # age-specific prior that infection was not omicron
  N_pos_control = 122,  # control data from Fenwick et al., 2021
  control_tp = 118,
  N_neg_control = 84,
  control_fp = 0,
  sd_sigma = 2,
  sd_od = 2,
  child_age_cat = neut_data$age < 12
) %>%
  append(priors)

# Write to json
cmdstanr::write_stan_json(stan_data, file = makeNeutStanDataName(opt))

# Save all to file
saveRDS(
  list(
    stan_data = stan_data,
    design_mat_full = design_mat_full,
    design_mat_missing = design_mat_missing
  ),
  file = makeNeutStanInputName(opt))


# Prepare post-stratification weights -------------------------------------

# Calculated GE vacc status data by age  
# !!NB: This would need to be modified for other settings
vacc_status_by_age = tibble::tribble(
  ~age_cat_UEP,  ~vacc_status_3way,  ~pop_total, ~strata_total, ~percentage,
  "[0,6)",          "Boosted",     31849.8,             0,           0,
  "[6,12)",          "Boosted", 32144.57143,             0,           0,
  "[12,18)",          "Boosted",       31902,        1955.4, 6.129396276,
  "[18,25)",          "Boosted",       43021,        9283.9, 21.57992608,
  "[25,35)",          "Boosted",     70489.5,       20486.5, 29.06319381,
  "[35,50)",          "Boosted",    114727.5,         42982, 37.46442658,
  "[50,65)",          "Boosted",       96005,         51167, 53.29618249,
  "[65,75)",          "Boosted",     42196.5,       28753.5,  68.1419075,
  "[75,105]",          "Boosted",     44096.5,       33756.5, 76.55142698,
  "all",          "Boosted",      506343,        190620, 37.64641755,
  "[0,6)", "Dose 1 or 2 only",     31849.8,         682.8,  2.14381252,
  "[6,12)", "Dose 1 or 2 only", 32144.57143,   1563.428571, 4.863740601,
  "[12,18)", "Dose 1 or 2 only",       31902,       14356.8, 45.00282114,
  "[18,25)", "Dose 1 or 2 only",       43021,       20166.3, 46.87547942,
  "[25,35)", "Dose 1 or 2 only",     70489.5,         32796, 46.52607835,
  "[35,50)", "Dose 1 or 2 only",    114727.5,       49392.5, 43.05201456,
  "[50,65)", "Dose 1 or 2 only",       96005,       31268.5, 32.56965783,
  "[65,75)", "Dose 1 or 2 only",     42196.5,          8864, 21.00648158,
  "[75,105]", "Dose 1 or 2 only",     44096.5,        7618.5, 17.27688139,
  "all", "Dose 1 or 2 only",      506343,        169307, 33.43721548,
  "[0,6)",          "No vacc",     31849.8,         31167, 97.85618748,
  "[6,12)",          "No vacc", 32144.57143,   30581.14286,  95.1362594,
  "[12,18)",          "No vacc",       31902,       15589.8, 48.86778258,
  "[18,25)",          "No vacc",       43021,       13570.8,  31.5445945,
  "[25,35)",          "No vacc",     70489.5,         17207, 24.41072784,
  "[35,50)",          "No vacc",    114727.5,         22353, 19.48355887,
  "[50,65)",          "No vacc",       96005,       13569.5, 14.13415968,
  "[65,75)",          "No vacc",     42196.5,          4579, 10.85161092,
  "[75,105]",          "No vacc",     44096.5,        2721.5, 6.171691631,
  "all",          "No vacc",      506343,        146416, 28.91636697
)

# Full sample inf by vacc by age (for post-stratification if required)
# !!NB: This is based on the full seroprevalence sample.
inf_latest_by_vacc_by_age_full_sample <- readRDS(opt$full_data) %>%
  drop_na(vacc_status_3way, inf_latest, age_cat_UEP) %>%
  count(vacc_status_3way, inf_latest, age_cat_UEP) %>%
  add_count(vacc_status_3way, age_cat_UEP, wt = n, name = "Total_per_vacc") %>%
  mutate(pc_inf_per_vacc = 100*n/Total_per_vacc)

post_strat_finegrained_weights_full_sample <- full_join(
  inf_latest_by_vacc_by_age_full_sample, 
  vacc_status_by_age) %>% 
  mutate(strata_size_N = pc_inf_per_vacc/100*pop_total*percentage/100) %>%
  select(age_cat_UEP, vacc_status_3way, inf_latest, strata_size_N) %>%
  drop_na() %>% 
  complete(age_cat_UEP, vacc_status_3way, inf_latest, fill = list(strata_size_N = 0))

saveRDS(post_strat_finegrained_weights_full_sample, 
        here::here("data", "processed", paste0(format(Sys.time(), "%Y-%m-%d-%H-%M"), "_neutralising_post_strat_finegrained_weights_full_sample.rds")))
