# These are utility functions for the analysis
makeSuffix <- function(opt) {
  if (!is.null(opt$suffix)) {
    str_c("_", opt$suffix)
  } else {
    ""
  }
}

makeProcessedNeutDataName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  here::here("data", "processed", str_glue("neutralization_data_{str}_{opt$age_cat}{makeSuffix(opt)}.rds"))
}

makeNeutStanDataName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  here::here("data", "processed", str_glue("neut_stan_data_{str}_{opt$model}_{opt$age_cat}{makeSuffix(opt)}.json"))
}

makeNeutStanInputName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  here::here("data", "processed", str_glue("neut_stan_input_{str}_{opt$model}_{opt$age_cat}{makeSuffix(opt)}.rds"))
}

makeNeutStanOutputName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  here::here("data", "processed", str_glue("neut_stan_output_{opt$outcome}_{opt$model}_{str}_{opt$age_cat}{makeSuffix(opt)}.rds"))
}

makeNeutStanGenName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  here::here("data", "processed", str_glue("neut_stan_genquant_{opt$outcome}_{opt$model}_{str}_{opt$age_cat}{makeSuffix(opt)}.rds"))
}

makeNeutPostStratName <- function(opt) {
  str <- janitor::make_clean_names(opt$variant)
  if (opt$outcome == "thresh") {
    outcome <-""
  } else {
    outcome <- "_od"
  }
  here::here("data", "processed", str_glue("neut_poststrat_{opt$model}{outcome}_{str}_{opt$age_cat}{makeSuffix(opt)}.rds"))
}

flattenStr <- function(str){
  str %>%
    str_to_lower() %>%
    str_replace_all("\n", "_") %>%
    str_replace_all("/", "-") %>%
    str_replace_all("\\.", "_")
}

parseResFileName <- function(str) {
  model_type <- str_extract(str, "additive|interaction")
  tibble(
    model_type = model_type,
    outcome_type  = str_extract(str, "thresh_od|thresh|reg"),
    variant  = str_extract(str, str_glue("(?<={model_type}_).*(?=_age)")),
    age_cat = str_extract(str, "age_cat_uep|age_cat_10yr")
  )
}

extractOutput <- function(res_file) {
  res <- readRDS(res_file)
  res_info <- parseResFileName(res_file)

  stan_input <- readRDS(makeNeutStanInputName(list(variant = res_info$variant,
                                                   model = res_info$model_type,
                                                   age_cat = res_info$age_cat)))

  if (res_info$outcome_type == "reg") {
    variables <- c("beta", "beta_hurdle")
  } else {
    variables <- c("beta", "odd_ratios")
  }

  params <- try(res$summary(variables = variables, .cores = 4,
                            mean = mean, ~ quantile(. , c(0.025, .975))) %>%
                  rename(q025 = `2.5%`,
                         q975 = `97.5%`) %>%
                  mutate(var = rep(colnames(stan_input$design_mat_full), length(variables))))

  if (!inherits(params, "try-error")) {
    loo <- loo::loo(res$draws("log_lik"))

    params %>%
      mutate(loo = list(loo)) %>%
      bind_cols(res_info)

  } else {
    res_info
  }
}


getGEVaccUEPCat <- function() {
  tribble(
    ~Category, ~GE_vac_pc,
    "[0,6)",    2.17,
    "[6,12)",   4.86,
    "[12,18)",  51.13,
    "[18,25)",  68.45,
    "[25,35)",  75.59,
    "[35,50)",  80.52,
    "[50,65)",  85.87,
    "[65,75)",  89.15,
    "[75,105]", 93.83,
    "Male",     69.6,
    "Female",   72.48,
    "All",      71.08)
}


# Define post-stratifications
makePostStratMap <- function(df, cov_list) {

  # Make filter
  f <- map_chr(1:length(cov_list),
               ~ str_c(names(cov_list)[.], " =='", cov_list[[.]],
                       "'")) %>%
    str_c(collapse = " & ")

  # Subset to retain
  sub <- filter(df, eval(parse(text = f)))

  # Make numbers
  tibble(
    expr = f,
    row_ind = list(sub$row),
    n = nrow(sub)
  )
}

fracToPrettyPc <- function(x) {
  formatC(x * 100, digits = 1, format = "f")
}

getVariantLevels <- function() {
  # Variants in chronological order
  c('any',
    'd614g',
    'alpha',
    'beta',
    'gamma',
    'delta',
    'iota',
    'kappa',
    'lambda',
    'omicron_ba_1',
    'omicron_ba_2',
    'omicron_ba_2_12_1',
    'omicron_ba_4-ba_5')
}


getVariantDict <- function() {
  # Variants in chronological order
  c('any'   = 'any',
    'd614g' = 'D614G',
    'alpha' = 'Alpha',
    'beta' = 'Beta',
    'gamma' = 'Gamma',
    'delta' = 'Delta',
    'iota' = 'Iota',
    'kappa' = 'Kappa',
    'lambda' = 'Lambda',
    'omicron_ba_1' = 'Omicron BA.1',
    'omicron_ba_2' = 'Omicron BA.2',
    'omicron_ba_2_12_1' = 'Omicron BA.2.12.1',
    'omicron_ba_4-ba_5' = 'Omicron BA.4/BA.5')
}

getVariantColors <- function() {
  # Colors excluding kappa and iota
  c("gray15", "#D93434", "#DB7734", "#CDDB34", "#34DB3D", "#34BADB",
    "#343ADB", "#8534DB", "#BA34DB", "#D534DB", "#DB34A3")
}

getVOIColors <- function() {
  # Colors excluding kappa and iota
  getVariantColors()[getVariantLevels()[-c(7, 8)] %in% getVOI()]
}

getVarLevels <- function() {
  c("all", "f", "m", "[0,6)", "[6,12)", "[12,18)",
    "[18,25)", "[25,35)", "[35,50)", "[50,65)",
    "[65,75)", "[75,105]", "uninfected", "inf_bef_2022",
    "inf_2022", "Novacc", "Dose1or2only", "Boosted",
    "uninfected:Novacc",
    "inf_bef_2022:Novacc",
    "inf_2022:Novacc",
    "uninfected:Dose1or2only",
    "inf_bef_2022:Dose1or2only",
    "inf_2022:Dose1or2only",
    "uninfected:Boosted",
    "inf_bef_2022:Boosted",
    "inf_2022:Boosted")
}

getVOI <- function(all = T) {
  voi <- c('any',
           'd614g',
           'alpha',
           'delta',
           'omicron_ba_1',
           'omicron_ba_2',
           'omicron_ba_2_12_1',
           'omicron_ba_4-ba_5')

  if (all) {
    return(voi)
  } else {
    return(voi[-1])
  }
}

getVarDict <- function() {
  c(
    "all" = "All",
    "sex" = "Sex",
    "age_cat" = "Age",
    "inf_latest" = "Infection",
    "vacc_status_3way" = "Vaccination",
    "inf_latest:vacc_status_3way" = "Infection/Vaccination"
  )
}

getVarLabelDict <- function() {
  c(
    "all" = "all",
    "f" = "women",
    "m" = "men",
    "uninfected" = "uninfected",
    "inf_bef_2022" = "latest infection\nbefore 2022",
    "inf_2022" = "latest infection\nin 2022",
    "[0,6)" = "0-5",
    "[6,12)" = "6-11",
    "[12,18)" = "12-17",
    "[18,25)" = "18-24",
    "[25,35)" = "25-34",
    "[35,50)" = "35-49",
    "[50,65)" = "50-64",
    "[65,75)" = "65-74",
    "[75,105]" = "75+",
    "Novacc" = "unvaccinated",
    "Dose1or2only" = "no booster dose",
    "Boosted" = "booster dose",
    "uninfected:Novacc" = "unvaccinated\nuninfected",
    "inf_bef_2022:Novacc" = "unvaccinated\nlatest infection before 2022",
    "inf_2022:Novacc" = "unvaccinated\nlatest infection in 2022",
    "uninfected:Dose1or2only" = "vaccinated no booster\nuninfected",
    "inf_bef_2022:Dose1or2only" = "vaccinated no booster\nlatest infection before 2022",
    "inf_2022:Dose1or2only" = "vaccinated no booster\nlatest infection in 2022",
    "uninfected:Boosted" = "vaccinated and booster\nuninfected",
    "inf_bef_2022:Boosted" = "vaccinated and booster\nlatest infection before 2022",
    "inf_2022:Boosted" = "vaccinated and booster\nlatest infection in 2022"
    # "uninfected:Novacc" = "uninfected and\nunvaccinated",
    # "uninfected:Dose1or2only" = "uninfected and\nno booster dose",
    # "uninfected:Boosted" = "uninfected and\nbooster dose",
    # "inf_bef_2022:Novacc" = "latest infection\nbefore 2022 and\nunvaccinated",
    # "inf_bef_2022:Dose1or2only" = "latest infection\nbefore 2022 and\nno booster dose",
    # "inf_bef_2022:Boosted" = "latest infection\nbefore 2022 and\nbooster dose",
    # "inf_2022:Novacc" = "latest infection\nin 2022 and\nunvaccinated",
    # "inf_2022:Dose1or2only" = "latest infection\nin 2022 and\nno booster dose",
    # "inf_2022:Boosted" = "latest infection\nin 2022 and\nbooster dose"
  )
}


makeNeutPlot <- function(bars = F,
                         points = T,
                         all_variants = F,
                         pdw = .7,
                         df,
                         keep_var_cat = unique(df$var_cat)) {
  subdat <-  df %>%
    {
      if (!all_variants) {
        filter(., (variant %in% getVariantDict()[getVOI()] |
                     outcome != "neutralizing"))
      } else {
        .
      }
    } %>%
    filter(var_cat %in% keep_var_cat)

  if (points) {
    p <- subdat %>%
      ggplot(aes(x = var, ymin = q025, ymax = q975, y = mean))
  } else {
    p <- subdat %>%
      filter(outcome == "neutralizing") %>%
      ggplot(aes(x = var, ymin = q025, ymax = q975, y = mean))
  }

  if (bars) {

    p <- p + geom_bar(data = subdat %>%
                        filter(outcome != "neutralizing"),
                      aes(y = mean, fill = outcome_name),
                      stat = "identity", position = "identity", alpha = .5) +
      scale_fill_manual(values = c("gray70", "gray35")) +
      guides(fill = guide_legend("Antibody origin"))
  }

  if (all_variants) {
    variant_cols <- getVariantColors()
    keep_variants <- setdiff(getVariantDict()[-1], c("Kappa" ,"Iota"))
  } else {
    variant_cols <- getVOIColors()
    keep_variants <-  getVariantDict()[getVOI()][-1]
  }


  pd <- position_dodge(width = pdw)

  p <- p + geom_point(aes(color = variant, pch = outcome_name), position = pd, alpha = .7) +
    geom_errorbar(aes(color = variant, lty = outcome_name), width = 0, position = pd, alpha = .3)

  p <- p +
    facet_grid(~ var_cat, scales = "free_x", space = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(x = "Covariate", y = "Seroprevalence") +
    scale_linetype_manual(values = rep(1, 3)) +
    scale_shape_manual(values = c(17, 16, 15)) +
    scale_color_manual(values = variant_cols[-1], na.value = "gray15",
                       limits = keep_variants) +
    guides(shape = guide_legend("Antibody measure", order = 1),
           lty = guide_legend("Antibody measure", order = 1),
           color = guide_legend("Variant", override.aes = list(shape = 15)), order = 2) +
    theme(legend.title = element_text(size = 10),
          panel.grid.major.y = element_line(color = "gray80", size = .3, linetype = 1),
          panel.grid.major.x = element_blank())  +
    theme(legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.key.height = unit(.15, units = "in"),
          strip.text = element_text(size = 7))

  p
}


makeFinalTable <- function(p_for_table,
                           sample_numbers,
                           all_variants = F) {
  p_for_table %>%
    mutate(txt_obs_pos = str_c(n_pos, " (", pc_pos, ")"),
           var = as.character(var),
           var_cat = as.character(var_cat)) %>%
    select(variant, var, var_cat, n, txt_obs_pos, estimate) %>%
    {
      if (all_variants) {
        .
      } else {
        filter(., variant %in% getVariantDict()[getVOI()])
      }
    } %>%
    pivot_wider(values_from = c("txt_obs_pos", "estimate"),
                id_cols = c("var_cat", "var"),
                names_from = "variant") %>%
    inner_join(sample_numbers %>%
                 mutate(variable = getVarDict()[variable],
                        label = getVarLabelDict()[label]), .,
               by = c("variable" = "var_cat", "label" = "var")) %>%
    select(variable, label, N,
           one_of(map(getVariantDict(),
                      ~ str_c(c("txt_obs_pos", "estimate"), ., sep = "_")) %>%
                    unlist())) %>%
    mutate(label = factor(label, levels = getVarLabelDict()[getVarLevels()])) %>%
    filter(!is.na(label)) %>%
    arrange(label) %>%
    mutate(variable = factor(variable, levels = getVarDict())) %>%
    as_grouped_data(groups = "variable")
}

makeFinalFlextable <- function(final_table) {

  variants <- map(getVariantDict(), function(x) {
    if(any(str_detect(colnames(final_table), x))) {
      res <- as.character(x)
    } else {
      res <-  NA_character_
    }
    res
  }) %>% .[!is.na(.)]

  new_colnames <-  set_names(c("", "", "n", rep(c("Obs.", "Est."), length(variants)) %>%
                                 as.list()),
                             colnames(final_table))

  final_table_flex <- flextable(final_table) %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    set_header_labels(values = new_colnames) %>%
    add_header_row(values = c("", variants), colwidths = c(3, rep(2, length(variants)))) %>%
    add_header_row(top = FALSE,
                   values = c("", rep(c("n (%)", "% (95% CrI)"), length(variants))), colwidths = c(3, rep(1, 2 * length(variants)))) %>%
    align(align = "center", part = "all") %>%
    fontsize(i = 1, size = 8, part = "header") %>%
    fontsize(i = 2, size = 8, part = "header") %>%
    fontsize(i = 3, size = 6, part = "header") %>%
    fontsize(j = 1:2, size = 7) %>%
    fontsize(i = 2:nrow(final_table),
             j = 3:ncol(final_table), size = 7)

  final_table_flex

}

makeFinalFlextaleExport <- function(final_table_flex, final_table) {
  final_table_flex %>%
    padding(i = 1:nrow(final_table),
            j = 1:ncol(final_table),
            padding.top = 1,
            padding.bottom = 1,
            padding.right = 0,
            padding.left = 0)
}

selectBestModel <- function(df,
                            se_diff_thresh = 2,
                            set_diff_tol = .1) {

  loo_compare_df <- loo::loo_compare(df$loo) %>%
    { x <- .
    as.data.frame(x) %>%
      mutate(id = rownames(x) %>% str_extract("[0-9]") %>% as.numeric())
    } %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(df %>%
                        select(model_type, outcome_type) %>%
                        mutate(id = row_number()), by = "id") %>%
    dplyr::select(-id) %>%
    mutate(model_id = str_c(outcome_type, model_type, sep = "-")) %>%
    arrange(desc(elpd_diff)) %>%
    mutate(rank = row_number())

  if (str_detect(df$outcome_type[1], "reg")) {
    model_complexity <- c(1,2) %>% set_names(str_c("reg-", c("additive", "interaction")))
  } else {
    model_complexity <- c(1:4) %>% set_names(str_c(rep(c("thresh-", "thresh_od-"), each = 2),
                                                   c("additive", "interaction")))
  }


  # Get the models that are not significanlty different than the best
  # best_model <- model_ics$variant[model_ics$loo_rank == 1]
  best_model_set <- dplyr::filter(loo_compare_df, elpd_diff == 0 | abs(elpd_diff) < se_diff * (se_diff_thresh - set_diff_tol))

  # If multiple models fall within the SE thresh limit choose the simplest
  if (nrow(best_model_set) == 1){
    best_model <- best_model_set
  } else {
    # Get the simplest model
    best_model <- best_model_set %>%
      dplyr::mutate(complexity = model_complexity[model_id]) %>%
      dplyr::arrange(complexity) %>%
      dplyr::slice(1)
  }

  loo_compare_df %>%
    left_join(best_model %>% mutate(best = T),
              by = c("elpd_diff", "se_diff", "elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic",
                     "model_type", "outcome_type", "model_id", "rank")) %>%
    replace_na(list(best = FALSE))
}
