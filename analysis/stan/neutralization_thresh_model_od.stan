//
// This Stan program defines a model for producing population-level estimates
// of the presence of neutralizing antibodies in the poppulation.
//

data {
  // Numbers
  int<lower=1> N_survey;    // numbber of participants in the survey
  int<lower=1> N_missing;   // numbber of participants with missing infection dates in the survey
  int<lower=0> N_age_cats;  // number of age categories (for priors on omicron vs. other infections)
  int<lower=1> N_col_omicron;

  // Observations
  int<lower=0, upper=N_survey> ind_full[N_survey-N_missing];  // indeces of full observations
  int<lower=0, upper=N_survey> ind_missing[N_missing];        // indeces of missing infection dates
  int<lower=0, upper=1> y[N_survey];  // observations for neutralization results (EC50 > 50)

  // Regression model
  int<lower=1> p_vars; //number of variables to adjust for
  matrix[N_survey-N_missing, p_vars] X;  // covariate model matrix (age and sex in these analyses)
  matrix[N_missing, p_vars] X_missing;   // model matrix for missing infection (default is other infection)
  int<lower=0, upper=p_vars> col_omicron[N_col_omicron];  // column index for omicron infection in model matrix
  int<lower=0, upper=p_vars> col_other[N_col_omicron];    // column index for other infection in model matrix
  int<lower=1, upper=N_age_cats> map_missing_age_cat[N_missing]; // map from observations with missing infections to age category

  // Priors
  real mu_alpha;
  real<lower=0> sd_alpha;
  real<lower=0> sd_beta;
  vector<lower=0, upper=1>[N_age_cats] prior_other;    // Prior that unknown infection dates are < Dec 25th 2021

  // Controls for sens/spec from Fenwick et al., 2021 (DOI: 10.1126/scitranslmed.abi8452)
  int<lower=0> N_pos_control;
  int<lower=0> control_tp;
  int<lower=0> N_neg_control;
  int<lower=0> control_fp;

  real<lower=0> sd_od;
}

transformed data {
  int<lower=0, upper=N_survey> N_full = N_survey - N_missing;    // number of full observations
  matrix[N_missing, p_vars] X_missing_omicron;
  vector[N_age_cats] log_prior_other = log(prior_other);
  vector[N_age_cats] log_prior_omicron = log(1-prior_other);

  // Initialize model matrix for missing infection dates assuming omicron infection
  for (i in 1:N_missing) {
    for (p in 1:p_vars) {
      X_missing_omicron[i,p] = X_missing[i,p];
    }
  }

  // Set omicron infection
  for (i in 1:N_missing) {
    for (j in 1:N_col_omicron) {
      if (X_missing[i, col_other[j]] == 1) {
        X_missing_omicron[i, col_other[j]] = 0;
        X_missing_omicron[i, col_omicron[j]] = 1;
      }
    }
  }
}
parameters {
  vector[p_vars] beta; // fixed regression coefficients
  real<lower=0.5, upper=1> sens;
  real<lower=0.5, upper=1> spec;
  real<lower=0> sigma;
  vector[N_survey] eta;
}
model {
  vector[N_full] p_full;   // probability of infection
  matrix[N_missing, 2] ll;

  // Logit-scale probabilities for participants with full observations
  p_full = inv_logit(X * beta + eta[ind_full]*sigma);

  {
    vector[N_missing] p_missing_other = inv_logit(X_missing * beta + eta[ind_missing]*sigma);
    vector[N_missing] p_missing_omicron = inv_logit(X_missing_omicron * beta + eta[ind_missing]*sigma);

    // Loglik of observations for missing infection dates: age-specif prior of omicron infection * bernoulli
    for (i in 1:N_missing) {
      ll[i, 1] = log_prior_other[map_missing_age_cat[i]] + bernoulli_lpmf(y[ind_missing[i]]| sens * p_missing_other[i] + (1-spec) * (1-p_missing_other[i]));
      ll[i, 2] = log_prior_omicron[map_missing_age_cat[i]] + bernoulli_lpmf(y[ind_missing[i]]| sens * p_missing_omicron[i] + (1-spec) * (1-p_missing_omicron[i]));
    }
  }

  y[ind_full] ~ bernoulli(sens * p_full + (1-spec) * (1-p_full));

  // Loglik for participant with missing infection dates
  for (i in 1:N_missing) {
    target += log_sum_exp(ll[i, ]);
  }

  // Priors for regression coefficients
  beta[1] ~ normal(mu_alpha, sd_alpha);
  beta[2:p_vars] ~ normal(0, sd_beta);

  // Priors for validation data
  control_tp ~ binomial(N_pos_control, sens);
  control_fp ~ binomial(N_neg_control, 1-spec);

  // Over-dispersion
  sigma ~ normal(0, sd_od);
  eta ~ std_normal();
}
generated quantities{
  real log_lik[N_survey];
  real odd_ratios[p_vars];

  for (i in 1:p_vars) {
    odd_ratios[i] = exp(beta[i]);
  }

  {
    vector[N_full] p_full;   // probability of infection

    // Logit-scale probabilities for participants with full observations
    p_full = inv_logit(X * beta + eta[ind_full]*sigma);

    // Log-likelihood for full obs
    for (i in 1:N_full) {
      log_lik[ind_full[i]] = bernoulli_lpmf(y[ind_full[i]]| sens*p_full[i] + (1-spec)*(1-p_full[i]));
    }
  }
  {
    vector[N_missing] p_missing_other = inv_logit(X_missing * beta + eta[ind_missing]*sigma);
    vector[N_missing] p_missing_omicron = inv_logit(X_missing_omicron * beta + eta[ind_missing]*sigma);

    // Log-lik for uncertain obs
    for (i in 1:N_missing) {
      vector[2] ll_tmp;
      ll_tmp[1] = log_prior_other[map_missing_age_cat[i]] + bernoulli_lpmf(y[ind_missing[i]]| sens*p_missing_other[i] + (1-spec)*(1-p_missing_other[i]));
      ll_tmp[2] = log_prior_omicron[map_missing_age_cat[i]] + bernoulli_lpmf(y[ind_missing[i]]| sens*p_missing_omicron[i] + (1-spec)*(1-p_missing_omicron[i]));
      log_lik[ind_missing[i]] = log_sum_exp(ll_tmp);
    }
  }
}

