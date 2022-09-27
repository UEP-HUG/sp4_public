# Seroprevalence of anti-SARS-CoV-2 antibodies and cross-variant neutralization capacity after the Omicron BA. 2 wave in Geneva, Switzerland 

This is a repository containing code used in the analysis used in the manuscript:

Zaballa, M.E., Perez-Saez, J., de Mestral, C., Pullen, N., Lamour, J., Turelli, P., Raclot, C., Baysson, H., Pennacchio, F., Villers, J. Duc, J., Richard V., Dumont R., Semaani C., Jutta A. Loizeau, Graindorge C., Lorthe E., Balavoine J-F, Pittet D., Schibler M., Vuilleumier N., Chappuis F., Kherad O.,ORCID Azman A.S., Posfay-Barbe K.M., Kaiser L., Trono D., Stringhini S., Guessous S. for the Specchio-COVID19 study group 2022. *Seroprevalence of anti-SARS-CoV-2 antibodies and cross-variant neutralization capacity after the Omicron BA. 2 wave in Geneva, Switzerland*. medRxiv.  https://doi.org/10.1101/2022.07.27.22278126 

## Repository structure

- All analysis scripts are in `analysis/`, with stan model codes in `analysis/stan`.
- `00_generate_synthetic_data.R` generates synthetic data to test the analysis pipeline
- `01_run_seroprev_analysis.R` runs seroprevalence estimation, post-stratification and result table creation.
- `02_prep_neutralization_data.R` prepares neutralization data for analysis and post-stratification
- `03_run_neutralization_model.R` runs the regression model on neutralization capacity
- `04_run_neutralization_poststratification.R` produces post-stratified estimates of neutralization capacity in the population.
- The data folder contains demographic data of the state of Geneva used in the analysis

## Running with synthetic data

The analysis pipeline can be tested using a synthetic dataset. To do so:

1. Run `analysis/00_generate_synthetic_data.R` to generate the data.
2. Run analysis steps 1 to 3. 
3. Run analysis step 4 passing the population strata and seroprevalence estimate file names as produced in steps 1 and 3. 

## Running with your own data

To run the analysis pipeline with your own data:

1. Prepare serological and neutralization datasets in the same format as outputted by `00_generate_synthetic_data.R`
2. Modify inputs for seroprevalence estimation:
- Modify function `getGE_age_cats` in `analysis/functions_Stan_multinomial.R` to read in population counts by age and sex in your target population.
- Modify variable `GE_vacc` in `analysis/01_run_seroprev_analysis.R` with vaccination coverage by age category in your target population.
3. Modify inputs for neutralization analysis in `analysis/02_prep_neutralization_data.R`:
- Modify probabilities of non-Omicron infection (before 2022) by age category in variable `prior_other_infection` to your target population
- Modify dataframe `vacc_status_by_age` to population counts by vaccination status and age in your target population.
4. Run steps 2-3 as when using a synthetic dataset.
