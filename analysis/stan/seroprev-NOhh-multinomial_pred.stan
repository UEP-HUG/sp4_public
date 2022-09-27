//
// This Stan program defines a model for adjusting a predicted
// seroincidence by the sensitivity and specificity of the diagnostic using data from lab validation study.
// used for analyses of the samples sent to EPFL for neutralisation analysis

// There are no household random effects

// Hierarchical model for specificity and sensitivity following Gelman & Carpenter 2020
// For code see https://github.com/bob-carpenter/diagnostic-testing/blob/master/src/specificity-santa-clara/stan/santa-clara-hierarchical.stan
// and https://github.com/bob-carpenter/diagnostic-testing/blob/master/src/specificity-santa-clara/R/rstan-santa-clara.R
// For initial discussion see https://statmodeling.stat.columbia.edu/2020/05/01/simple-bayesian-analysis-inference-of-coronavirus-infection-rate-from-the-stanford-study-in-santa-clara-county/
// and for final paper see http://www.stat.columbia.edu/~gelman/research/published/specificity.pdf

// We have input data from a validation set (which determines the diagnostic performance) and the survey (from which we'd like to estimate seropos).
data {
  int<lower=1> N_survey; //number of participants in the survey
  matrix<lower=0, upper=1>[N_survey, 2] survey_pos; //observations for S and N results

  int<lower=1> p_vars; //number of variables to adjust for
  matrix[N_survey, p_vars] X; //covariate model matrix (age and sex in these analyses)
  int<lower=1> N_comb; // number of unique combinations
  matrix[N_comb, p_vars] X_comb; //model matrix with unique combinations
  vector<lower=0, upper=1>[N_comb] GE_vacc; // Geneva vaccination proportions per sex/age transformed from cantonal data

  int<lower=0, upper=1> vaccinated[N_survey]; // vaccination status
  vector<lower=0>[3] p_infect_prior; // Drichilet prior for the three types of natural infection immune response (++, -+, +-)

  int<lower=0> J_studies_S; // number of Roche S studies with validation data
  int<lower=0> J_studies_N; // number of Roche N studies with validation data

  int<lower=0> N_pos_control_S[J_studies_S]; //number of positive controls in the validation data of Roche S
  int<lower=0> control_tp_S[J_studies_S]; // number of true positive tests in the validation data of Roche S
  int<lower=0> N_neg_control_S[J_studies_S]; // number of negative controls in the validation data of Roche S
  int<lower=0> control_fp_S[J_studies_S];// number of false positives by the diagnostic test in the validation study of Roche S

  int<lower=0> N_pos_control_N[J_studies_N]; //number of positive controls in the validation data of Roche N
  int<lower=0> control_tp_N[J_studies_N]; // number of true positive tests in the validation data of Roche N
  int<lower=0> N_neg_control_N[J_studies_N]; // number of negative controls in the validation data of Roche N
  int<lower=0> control_fp_N[J_studies_N];// number of false positives by the diagnostic test in the validation study of Roche N
}
transformed data {
  int<lower=0, upper=1> y[N_survey, 4];// one-hot encoded matrix of observations for each participant

  // Build observations: ++, -+, +-, --
  // Initialize
  for (i in 1:N_survey) {
    for (j in 1:4) {
      y[i, j] = 0;
    }
  }

  // Fill - this must be the same order as in theta below
  for (i in 1:N_survey) {
    if (survey_pos[i, 1] == 1 && survey_pos[i, 2] == 1) {// {S+,N+}
    y[i, 1] = 1;
    } else if (survey_pos[i, 1] == 0 && survey_pos[i, 2] == 1) {// {S-,N+}
    y[i, 2] = 1;
    } else if (survey_pos[i, 1] == 1 && survey_pos[i, 2] == 0) {// {S+,N-}
    y[i, 3] = 1;
    } else if (survey_pos[i, 1] == 0 && survey_pos[i, 2] == 0) {// {S-,N-}
    y[i, 4] = 1;
    }
  }
}
parameters {
  vector[p_vars] beta; // fixed regression coefficients
  vector[p_vars] beta_vacc;

  // Roche S
  real mu_logit_spec_S; // mean for hierarchical specificity
  real mu_logit_sens_S; // mean for hierarchical sensitivity
  real<lower = 0> sigma_logit_spec_S; // variance for hierarchical specificity
  real<lower = 0> sigma_logit_sens_S; // variance for hierarchical sensitivity
  vector<offset = mu_logit_spec_S, multiplier = sigma_logit_spec_S>[J_studies_S] logit_spec_S; // normal dist on log odds scale for hierarchical specificity
  vector<offset = mu_logit_sens_S, multiplier = sigma_logit_sens_S>[J_studies_S] logit_sens_S; // normal dist on log odds scale for hierarchical sensitivity

  // Roche N
  real mu_logit_spec_N; // mean for hierarchical specificity
  real mu_logit_sens_N; // mean for hierarchical sensitivity
  real<lower = 0> sigma_logit_spec_N; // variance for hierarchical specificity
  real<lower = 0> sigma_logit_sens_N; // variance for hierarchical sensitivity
  vector<offset = mu_logit_spec_N, multiplier = sigma_logit_spec_N>[J_studies_N] logit_spec_N; // normal dist on log odds scale for hierarchical specificity
  vector<offset = mu_logit_sens_N, multiplier = sigma_logit_sens_N>[J_studies_N] logit_sens_N; // normal dist on log odds scale for hierarchical sensitivity

  // Additional infection parameters
  simplex[3] p_infect_response;    // the probability of having ++/-+/+- responses to natural infection
  real<lower=0, upper = 1> p_vacc_response;    // prob of immune response upon vaccination

}

transformed parameters {
  vector<lower=0, upper=1>[N_survey] p;   // probability of infection
  vector<lower=0, upper=1>[N_survey] p_vacc;   // probability of vaccination

  vector<lower=0, upper=1>[N_survey] P1;  // {S+,N+} natural infection (I)
  vector<lower=0, upper=1>[N_survey] P2;  // {S-,N+} I
  vector<lower=0, upper=1>[N_survey] P3;  // {S+,N-} I
  vector<lower=0, upper=1>[N_survey] P_pospos;  // {S+,N+} test result: {S+,N+} I OR ({S-,N+} I AND V)
  vector<lower=0, upper=1>[N_survey] P_negpos;  // {S-,N+} test result: {S-,N+} I AND ~V
  vector<lower=0, upper=1>[N_survey] P_posneg;  // {S+,N-} test result: {S+,N-} I OR (~I AND V)
  vector<lower=0, upper=1>[N_survey] P_negneg;  // {S-,N-} test result: ~I AND ~V

  matrix<lower=0, upper=1>[N_survey, 4] theta;

  vector[J_studies_S] spec_S = inv_logit(logit_spec_S); // unwind the logistic/log odds transform (I think so its normal value is used in the model??)
  vector[J_studies_S] sens_S = inv_logit(logit_sens_S);

  vector[J_studies_N] spec_N = inv_logit(logit_spec_N); // unwind the logistic/log odds transform (I think so its normal value is used in the model??)
  vector[J_studies_N] sens_N = inv_logit(logit_sens_N);
  vector<lower=0, upper=1>[N_survey] Pv; // Probability of S+,N- response to vaccination


  for (i in 1:N_survey) {
    Pv[i] = vaccinated[i] * p_vacc_response;
  }

  // Probability of natural infection
  p = inv_logit(X * beta);
  p_vacc = inv_logit(X * beta_vacc);

  // print("p_infect_response ", p_infect_response);
  // Probabilities of responses to natural infection
  P1 = p * p_infect_response[1]; // {S+,N+} I
  P2 = p * p_infect_response[2]; // {S-,N+} I
  P3 = p * p_infect_response[3]; // {S+,N-} I

  // print("p[1]= ", p[1], " P1[1]=", P1[1], " P2[1]=", P2[1], " P3[1]=", P3[1]);

  // Probabilities of "true" S/N status

  P_posneg = Pv .* (1 - p) + P3;
  P_negpos = P2 .* (1 - Pv);
  P_negneg = 1 - Pv .* (1 - p) - p;
  // True formula
  // P_pospos = P1 + P2 .* Pv;
  P_pospos = 1  - (P_posneg + P_negpos + P_negneg);

  // Checks
  // print("P_pospos ", P_pospos[1], " P_negpos ", P_negpos[1], " P_posneg ", P_posneg[1], " P_negneg ", P_negneg[1]);
  // print("Sum of P's ", P_pospos[1] + P_negpos[1] + P_posneg[1] + P_negneg[1]);
  // print("sens_S[1] ", sens_S[1], " sens_N[1] ", sens_N[1], " spec_S[1] ", spec_S[1], " spec_N[1] ", spec_N[1]);

  // Adjusted probabilities accounting for S/N test performance
  theta[, 1] = sens_S[1]*sens_N[1]*P_pospos + (1-spec_S[1])*sens_N[1]*P_negpos + sens_S[1]*(1-spec_N[1])*P_posneg + (1-spec_S[1])*(1-spec_N[1])*P_negneg;
  theta[, 2] = (1-sens_S[1])*sens_N[1]*P_pospos + spec_S[1]*sens_N[1]*P_negpos + (1-sens_S[1])*(1-spec_N[1])*P_posneg + spec_S[1]*(1-spec_N[1])*P_negneg;
  theta[, 3] = sens_S[1]*(1-sens_N[1])*P_pospos + (1-spec_S[1])*(1-sens_N[1])*P_negpos + sens_S[1]*spec_N[1]*P_posneg + (1-spec_S[1])*spec_N[1]*P_negneg;
  // This is the true formula, to avoid rounding errors we use the complement
  // theta[, 4] = (1-sens_S[1])*(1-sens_N[1])*P_pospos + spec_S[1]*(1-sens_N[1])*P_negpos + (1-sens_S[1])*spec_N[1]*P_posneg + spec_S[1]*spec_N[1]*P_negneg;
  theta[, 4] =  1 - (theta[, 1] +  theta[, 2] +  theta[, 3]);

  // Checks
  // print("theta1 ", theta[1,]);
  // print("Sum theta1 ", sum(theta[1,]));
}

//  We observe 'y' cases as a multinomial distribution based on the
//  N and S tests with observations coming as the 'theta' vector with the four
//  probs based on test accuracy above
model {

  // Individual status probabilities
  for (i in 1:N_survey) {
    target += multinomial_lpmf(y[i, ] | to_vector(theta[i, ]));
  }

  // probability of vaccination
  target += bernoulli_lpmf(vaccinated | p_vacc);

  // prior on different types of infection responses
  target += dirichlet_lpdf(p_infect_response | p_infect_prior);

  // Prior on probability of vaccine response
  target +=  beta_lpdf(p_vacc_response| 10, .1);

  target+= binomial_lpmf(control_tp_S | N_pos_control_S, sens_S);
  target+= binomial_lpmf(control_fp_S | N_neg_control_S, 1-spec_S);
  target+= binomial_lpmf(control_tp_N | N_pos_control_N, sens_N);
  target+= binomial_lpmf(control_fp_N | N_neg_control_N, 1-spec_N);

  target+= normal_lpdf(beta | 0, 1); // priors for coefficients
  target+= normal_lpdf(beta_vacc | 0, 5); // priors for coefficients

  target+= normal_lpdf(logit_spec_S | mu_logit_spec_S, sigma_logit_spec_S); // hierarchical prior for spec
  target+= normal_lpdf(logit_sens_S | mu_logit_sens_S, sigma_logit_sens_S); // hierarchical prior for sens
  target+= normal_lpdf(sigma_logit_spec_S | 0, 1);// weak truncated half-normal prior for spec variance
  target+= normal_lpdf(sigma_logit_sens_S | 0, 1);// weak truncated half-normal prior for sens variance
  target+= normal_lpdf(mu_logit_spec_S | 4, 2);  // weak prior on mean of distribution of spec
  target+= normal_lpdf(mu_logit_sens_S | 4, 2);  // weak prior on mean of distribution of sens

  target+= normal_lpdf(logit_spec_N | mu_logit_spec_N, sigma_logit_spec_N); // hierarchical prior for spec
  target+= normal_lpdf(logit_sens_N | mu_logit_sens_N, sigma_logit_sens_N); // hierarchical prior for sens
  target+= normal_lpdf(sigma_logit_spec_N | 0, 1);// weak truncated half-normal prior for spec variance
  target+= normal_lpdf(sigma_logit_sens_N | 0, 1);// weak truncated half-normal prior for sens variance
  target+= normal_lpdf(mu_logit_spec_N | 4, 2);  // weak prior on mean of distribution of spec
  target+= normal_lpdf(mu_logit_sens_N | 4, 2);  // weak prior on mean of distribution of sens
}

generated quantities {
  // vector[N_survey] log_lik;
  vector[N_comb] p_est;
  vector[N_comb] p_est_vacc;
  matrix[N_comb, 4] post_probs;
  vector[N_comb] p_any; // probability of any antibody response

  // Post-stratification using Geneva vaccination data
  // (at least as best as I could transform them to work with our age cats)
  matrix[N_comb, 4] post_probs_with_GE_vacc;// using GE vacc proportions
  vector[N_comb] p_any_with_GE_vacc;// using GE vacc proportions

  for (i in 1:N_comb) {
    real p_vacc_resp;
    real p_vacc_resp_GE;

    // infection
    p_est[i] = inv_logit(X_comb[i, ] * beta);

    // vaccination
    p_est_vacc[i] = inv_logit(X_comb[i, ] * beta_vacc);
    p_vacc_resp = p_est_vacc[i] * p_vacc_response;

    post_probs[i, 2] =  p_vacc_resp * (1 - p_est[i]) + p_est[i] * p_infect_response[3];
    post_probs[i, 3] = p_est[i] * p_infect_response[2] * (1 - p_vacc_resp);
    post_probs[i, 4] = 1 - p_vacc_resp * (1 -  p_est[i]) -  p_est[i] ;

    // True formula
    // P_pospos = P1 + P2 .* Pv;
    post_probs[i, 1] = 1  - sum(post_probs[i, 2:4]);

    // Using Geneva vaccination levels (cantonal estimates) passed as data
    p_vacc_resp_GE = GE_vacc[i] * p_vacc_response;

    post_probs_with_GE_vacc[i, 2] =  p_vacc_resp_GE * (1 - p_est[i]) + p_est[i] * p_infect_response[3];
    post_probs_with_GE_vacc[i, 3] = p_est[i] * p_infect_response[2] * (1 - p_vacc_resp_GE);
    post_probs_with_GE_vacc[i, 4] = 1 - p_vacc_resp_GE * (1 -  p_est[i]) -  p_est[i] ;

    // True formula
    // P_pospos = P1 + P2 .* Pv;
    post_probs_with_GE_vacc[i, 1] = 1  - sum(post_probs_with_GE_vacc[i, 2:4]);
  }

  p_any = 1 - post_probs[, 4];
  p_any_with_GE_vacc = 1 - post_probs_with_GE_vacc[, 4];

  // print("Done gq");

  // for(i in 1:N_survey){
  //   log_lik[i] = multinomial_lpmf(y[i, ] | to_vector(theta[i, ]));
  // }

}

