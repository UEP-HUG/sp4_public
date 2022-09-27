# This script runs the poststratification of neutralization

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
  make_option(c("-o", "--outcome"), default = "thresh", action ="store",
              type = "character", help = "Model outcome"),
  make_option(c("-s", "--suffix"), default = NULL, action ="store",
              type = "character", help = "Suffix in names"),
  make_option(c("-w", "--strata_weights"), 
              default = here::here("data", "processed", "2022-09-27-15-40_neutralising_post_strat_finegrained_weights_full_sample.rds"),
              action ="store", type = "character", 
              help = "Strata weight of age/sex/infection/vaccination based on analysis sample. 
              This file is produced by 02_prep_neutralization_data.R."),
  make_option(c("-x", "--seroprev_output"), 
              default = here::here("data", "processed", "55512dcf_HH_multinomial_test_results-list_GE_vacc_2022-09-27-15-16-53.rds"),
              action ="store", type = "character", 
              help = "Post-stratified seroprevalence estimates. 
              This file is produced by 01_run_seroprev_analysis.R.")
)

opt <- parse_args(OptionParser(option_list = option_list)) 

# Load data ---------------------------------------------------------------

# Model fit
fit <- readRDS(makeNeutStanOutputName(opt))

# Input data
stan_input <- readRDS(makeNeutStanInputName(opt))

# Population weights of infection/vaccination status.
# In the analysis these were based on our seroprevalence estimates and on
# sample.
pop_weights <- readRDS(opt$strata_weights) %>%
  janitor::clean_names() %>%
  ungroup() %>%
  rename(age_cat = !!rlang::sym(opt$age_cat)) %>%
  group_by(age_cat, vacc_status_3way) %>%
  mutate(n_inf_tot_sample = sum(strata_size_n[inf_latest != "uninfected"]),
         n_tot_sample = sum(strata_size_n))

# Load estimates of infections from SP4
sp4_inf <- readRDS(opt$seroprev_output)$pop_cats_infect_probs_GE_vacc %>%
  group_by(age_cat) %>%
  summarise(pc_inf_sp4 = weighted.mean(seropos, strata_pop) * 100)

neut_data <- readRDS(here::here("data", "processed", "neutralization_data_all.rds"))

inf_comp <- pop_weights %>%
  group_by(age_cat) %>%
  summarise(n_notinf_age = sum(strata_size_n[inf_latest == "uninfected"]),
            n_inf_before_age = sum(strata_size_n[inf_latest == "inf_bef_2022"]),
            n_inf_2022_age = sum(strata_size_n[inf_latest == "inf_2022"]),
            n_tot_age = sum(strata_size_n),
            n_inf_tot_age = n_inf_2022_age + n_inf_before_age,
            pc_inf_file_age = (1-n_notinf_age/n_tot_age) * 100,
            pc_inf_2022_age = n_inf_2022_age/(n_inf_tot_age)) %>%
  inner_join(sp4_inf, by = c("age_cat" = "age_cat"))


# Note: infection estimates seem to require adjustment to our SP4 inference
# Adjust infection estimates:
pop_weights2 <- inf_comp %>%
  mutate(pc_ratio = pc_inf_sp4/pc_inf_file_age) %>%
  inner_join(pop_weights, .) %>%
  group_by(age_cat, vacc_status_3way) %>%
  group_modify(function(x, y) {
    n_tot <- sum(x$strata_size_n)
    n_inf_old <- sum(x$strata_size_n[x$inf_latest != "uninfected"], na.rm = T)
    n_inf_new <- n_inf_old * x$pc_ratio[1]
    n_inf_2022_old <- x$strata_size_n[x$inf_latest == "inf_2022"]
    n_inf_bef_2022_old <- x$strata_size_n[x$inf_latest == "inf_bef_2022"]
    ratio_inf <- n_inf_2022_old / (n_inf_2022_old + n_inf_bef_2022_old)
    if (is.nan(ratio_inf)) {
      ratio_inf <- 0
    }
    n_inf_2022_new <- min(n_inf_new * ratio_inf, n_tot)
    n_inf_bef_2022_new <-  min(n_inf_new * (1-ratio_inf), n_tot - n_inf_2022_new)
    
    x$strata_size_n2[x$inf_latest == "inf_2022"] <- n_inf_2022_new
    x$strata_size_n2[x$inf_latest == "inf_bef_2022"] <- n_inf_bef_2022_new
    x$strata_size_n2[x$inf_latest == "uninfected"] <- n_tot - n_inf_2022_new - n_inf_bef_2022_new
    x
  }) %>%
  ungroup()

pop_weights_final <- bind_rows(
  pop_weights2 %>% mutate(sex = "m"),
  pop_weights2 %>% mutate(sex = "f")) %>%
  mutate(strata_size_final = strata_size_n2/2,
         inf_latest = factor(inf_latest, levels = levels(neut_data$inf_latest)),
         sex = factor(sex, levels = levels(neut_data$sex)),
         age_cat = factor(age_cat, levels = levels(neut_data$age_cat)),
         vacc_status_3way = factor(vacc_status_3way, levels(neut_data$vacc_status_3way)),
         full_cat = str_c("sex", sex, ":age_cat", age_cat, ":inf_latest", inf_latest, ":vacc", vacc_status_3way),
         row = row_number())

# Setup data for post-stratification --------------------------------------

# Model formula for regression
if (opt$model == "additive") {
  model_frml <- as.formula("~ age_cat + sex + inf_latest + vacc_status_3way")
} else if (opt$model == "interaction") {
  model_frml <- as.formula("~ age_cat + sex + inf_latest * vacc_status_3way")
} else {
  stop("Unknown model type.")
}

# Make unique combinations of covariates
X_post <- model.matrix(data = pop_weights_final, model_frml)

# Make postratitifaction numbers for each covariate subset of interest.
# Here we look at each level of each covariate
postrat_covars <- select_if(pop_weights_final, is.factor) %>%
  colnames() %>%
  map(., function(x) as.list(levels(pop_weights_final[[x]]) %>% set_names(x))) %>%
  unlist(recursive = F)

postrat_covars <- map(1:length(postrat_covars), function(x){
  list(postrat_covars[[x]]) %>% set_names(names(postrat_covars)[x])
})

# Inf/vacc combinations
inf_vacc_strata <- expand.grid(inf = levels(pop_weights_final$inf_latest),
                               vacc = levels(pop_weights_final$vacc_status_3way))

# Add combinations of infection and vaccination status
postrat_covars <- map(1:nrow(inf_vacc_strata),
                      function(x) {
                        list(inf_latest = inf_vacc_strata$inf[x],
                             vacc_status_3way = inf_vacc_strata$vacc[x])
                      }) %>%
  append(postrat_covars, .)

all_poststrat <- map_df(1:length(postrat_covars),
                        function(x) {
                          makePostStratMap(df = pop_weights_final,
                                           cov_list = postrat_covars[[x]]) %>%
                            mutate(postrat_name = str_replace_all(expr, "==", ":") %>%
                                     str_remove_all(" |'"))
                        })


saveRDS(all_poststrat, file = here::here("data", "processed", "all_postratification_covariates.rds"))

# Make stan data for generation
poststrat_data <- list(
  N_comb = nrow(X_post),
  X_comb = X_post,
  pop_strata = pop_weights_final$strata_size_final,
  N_poststrat = nrow(all_poststrat),
  M_poststrat = sum(all_poststrat$n),
  map_poststrat = all_poststrat$row_ind %>% unlist(),
  map_poststrat_start = c(1, cumsum(all_poststrat$n)[-nrow(all_poststrat)] + 1),
  map_poststrat_end = cumsum(all_poststrat$n)
) %>%
  append(rjson::fromJSON(file = makeNeutStanDataName(opt)))

# Gen quant ---------------------------------------------------------------

if (opt$outcome == "thresh") {
  gen_model <- cmdstanr::cmdstan_model(
    here::here("analysis", "stan", "neutralization_thresh_model_gen.stan"),
    force_recompile = F)
} else if (opt$outcome == "thresh_od") {
  gen_model <- cmdstanr::cmdstan_model(
    here::here("analysis", "stan", "neutralization_thresh_model_od_gen.stan"),
    force_recompile = F)
}

# Make draws
gen_quant <- gen_model$generate_quantities(fitted_params = fit,
                                           data = poststrat_data,
                                           parallel_chains = 4)

gen_quant$save_object(file = makeNeutStanGenName(opt))

p_postsrat <- gen_quant$summary(c("p", "p_poststrat"),
                                mean = mean, ~ quantile(. , c(0.025, .975))) %>%
  rename(q025 = `2.5%`,
         q975 = `97.5%`) %>%
  mutate(var = c("all", all_poststrat$postrat_name),
         var_cat = case_when(var == "all" ~ "all",
                             str_detect(var, "age") ~ "age_cat",
                             str_detect(var, "sex") ~ "sex",
                             str_detect(var, "age") ~ "age",
                             str_detect(var, "\\&") ~ "inf_latest:vacc_status_3way",
                             T ~ str_extract(var, ".*(?=:)")) %>%
           factor(levels = names(getVarDict())),
         var = case_when(str_detect(var_cat, "\\:", negate = T) ~ str_remove(var, str_c(var_cat, ":")),
                         T ~ str_remove_all(var, "inf_latest\\:|vacc_status_3way\\:") %>%
                           str_replace_all("\\&", ":")) %>%
           factor(levels = getVarLevels()),
         variant = opt$variant)

saveRDS(p_postsrat, file = makeNeutPostStratName(opt))
