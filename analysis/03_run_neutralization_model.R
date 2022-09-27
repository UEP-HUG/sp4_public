# This script runs the neutralization model

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
              type = "character", help = "Suffix in names")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Run stan ----------------------------------------------------------------

if (opt$outcome == "thresh") {
  model <- cmdstanr::cmdstan_model(
    here::here("analysis", "stan", "neutralization_thresh_model.stan"),
    force_recompile = F)
} else if (opt$outcome == "thresh_od") {
  model <- cmdstanr::cmdstan_model(
    here::here("analysis", "stan", "neutralization_thresh_model_od.stan"),
    force_recompile = F)
} else if (opt$outcome == "thresh_multiacc") {
  model <- cmdstanr::cmdstan_model(
    here::here("analysis", "stan", "neutralization_thresh_model_multiacc.stan"),
    force_recompile = F)
} else {
  stop("Unkown model")
}

# Make draws
fit <- model$sample(data = makeNeutStanDataName(opt),
                    seed = 1234,
                    init = 1,
                    chain_ids = seq_len(4),
                    chains = 4,
                    parallel_chains = 4,
                    threads_per_chain = NULL,
                    iter_warmup = 250,
                    iter_sampling = 1000,
                    max_treedepth = 12L,
                    metric = "unit_e",
                    adapt_delta = .95,
                    save_warmup = F,
                    refresh = 100)

fit$save_object(file = makeNeutStanOutputName(opt))
