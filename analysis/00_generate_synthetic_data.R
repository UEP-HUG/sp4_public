# This script generates the synthetic data on which to run the analysis
# 
# Data formats correspond to the origin datasets used in the analysis

# Preamble ----------------------------------------------------------------

library(tidyverse)

# Parameters --------------------------------------------------------------

N <- 1e3    # number of synthetic participants
N_hh_max <- round(N/1.5)   # number of households

# Event probabilities
p_vacc <- .7    # probability of being vaccinated
p_boost <- .3   # probability of having booster dose
p_inf <-  .6    # probability of being infected
p_test <- .5    # probability of being tested
p_inf_2022 <- .6    # probability of having a positive test in 2022

# Test specifications
sens_s <- .97
spec_s <- .98
sens_n <- .9
spec_n <- .98

set.seed(123)

if (!dir.exists(here::here("data", "processed")))
  dir.create(here::here("data", "processed"))

# Generate true states ----------------------------------------------------

true_pop <- map_df(1:N, function(x) {
  
  tibble(
    uid = x,
    hh_id = sample(1:N_hh_max, 1),
    age = round(runif(1, 0, 99)),
    sex = round(runif(1)),
    infected = runif(1) < p_inf,
    vaccinated = runif(1) < p_vacc
  ) %>% 
    mutate(
      # When infection occurred
      infected_when = case_when(infected & runif(1) < p_inf_2022 ~ "inf_2022",
                                infected ~ "inf_bef_2022",
                                T ~ "uninfected"),
      # Whether a test was available
      pos_test = case_when(infected ~ runif(1) < p_test,
                           T ~ FALSE),
      # Whether boosted or not
      boosted = case_when(vaccinated ~ runif(1) < p_boost,
                          T ~ FALSE),
      inf_latest = case_when(infected & runif(1) < p_test ~ infected_when,
                             infected ~ NA_character_,
                             T ~ infected_when),
      vacc_status_3way = case_when(!vaccinated  ~ "No vacc",
                                   runif(1) > p_boost ~ "Dose 1 or 2 only",
                                   T ~ "Boosted")
    )
}) %>% 
  mutate(
    uid = factor(uid),
    hh_id = factor(hh_id)
  )

# Set age categories and factor levels
true_pop <- true_pop %>% 
  mutate(
    sex = case_when(sex == 1 ~ "f", T ~ "m") %>% factor(),
    age_cat_10y = cut(age, c(seq(0, 80, 10), 105), 
                      right = FALSE, include.lowest = TRUE) %>% factor(),
    age_cat_UEP = cut(age, c(0, 6, 12, 18, 25, 35, 50, 65, 75, 105), 
                      right = FALSE, include.lowest = TRUE) %>% factor()
  )

saveRDS(true_pop, here::here("data", "processed", "synthetic_true_sample.rds"))

# Generate serological data -----------------------------------------------

sero_data <- true_pop %>% 
  rowwise() %>% 
  mutate(
    # anti-S and anti-N serology results
    S_interp = ifelse(runif(1) < 
                        (sens_s * (infected | vaccinated) +
                           (1-spec_s) * !(infected | vaccinated)),
                      "pos", "neg"),
    N_interp = ifelse(runif(1) < 
                        (sens_n * infected + (1-spec_n) * !infected),
                      "pos", "neg")
  ) %>% 
  ungroup() %>% 
  mutate(vaccinated = vacc_status_3way != "No vacc") %>%
  select(uid, hh_id, sex, contains("age"), S_interp, N_interp, vaccinated)

write_csv(sero_data, here::here("data", "processed", "synthetic_sero_sample.csv"))

# Generate neutralization data --------------------------------------------

# Make fake variants A-C
u_variants <- paste0("var-", c("a", "b", "c"))

# Mean of ec50 value for each variant
mean_variants <- log(c(3e2, 3e2, 1e2))
names(mean_variants) <- u_variants

neut_data <- true_pop %>% 
  group_by(uid) %>%
  group_modify(function(x, y) {
    
    map_df(u_variants, 
           ~ mutate(x,
                    ec50 = case_when(
                      infected ~ exp(rnorm(1, mean = mean_variants[.], sd = 1.5)),
                      T ~ rnorm(1, mean = .1, sd = .1)),
                    variant = .))
  }) %>% 
  ungroup() %>% 
  select(uid, sex, contains("age"), variant, ec50, inf_latest, vacc_status_3way)

write_csv(neut_data, here::here("data", "processed", "synthetic_neut_sample.csv"))
