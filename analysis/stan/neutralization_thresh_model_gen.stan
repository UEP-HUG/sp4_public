//
// This Stan program defines a model for producing population-level estimates
// of the presence of neutralizing antibodies in the poppulation.
//

data {
  // Numbers
  int<lower=1> N_survey;    // numbber of participants in the survey
  int<lower=1> N_missing;   // numbber of participants with missing infection dates in the survey
  int<lower=0> N_age_cats;  // number of age categories (for priors on omicron vs. other infections)
  int<lower=1> N_col_omicron;  // number of columns where omicron infections are involved
  int<lower=1> N_comb; // number of unique covariate combinationss

  // Observations
  int<lower=0, upper=1> y[N_survey];  // observations for neutralization results (EC50 > 50)

  // Regression model
  int<lower=1> p_vars; //number of variables to adjust for

  matrix[N_comb, p_vars] X_comb;  // covariate model matrix (age and sex in these analyses)

  // Post-stratification
  row_vector[N_comb] pop_strata;  // number of people in each strata
  int<lower=0> N_poststrat;    // number of unique covariates across which to post-stratify
  int<lower=0> M_poststrat;    // number of mappings for post-stratitication
  int<lower=1, upper=M_poststrat> map_poststrat_start[N_poststrat];
  int<lower=1, upper=M_poststrat> map_poststrat_end[N_poststrat];
  int<lower=1, upper=N_comb> map_poststrat[M_poststrat];
}
transformed data {
  real tot_pop = sum(pop_strata);
  row_vector[N_comb] pop_weights = pop_strata/tot_pop;    /// weights overall
  real tot_postrats[N_poststrat];  // total in each poststratification starta

  if (N_poststrat>0) {
    for (i in 1:N_poststrat) {
      tot_postrats[i] = 0;
      for (m in map_poststrat_start[i]:map_poststrat_end[i]) {
        tot_postrats[i] += pop_strata[map_poststrat[m]];
      }
    }
  }
}
parameters {
  vector[p_vars] beta; // fixed regression coefficients
}
generated quantities{
  vector[N_comb] p_ucomb;  // probabilities for each stratat
  real p;                  // post-stratified estimate for all the population
  real p_poststrat[N_poststrat];    // post-stratified estimate for all the population

  p_ucomb = inv_logit(X_comb * beta);

  // Post-stratification across all covariates for total estimates
  p = pop_weights * p_ucomb;

  // Post-stratification for user-specific
  if (N_poststrat>0) {
    for (i in 1:N_poststrat) {
      p_poststrat[i]  = 0;
      for (m in map_poststrat_start[i]:map_poststrat_end[i]) {
        p_poststrat[i] += p_ucomb[map_poststrat[m]] * pop_strata[map_poststrat[m]]/tot_postrats[i];
      }
    }
  }
}
