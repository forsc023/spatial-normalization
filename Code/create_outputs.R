
### This script reads in the prepared inputs for one cluster / section and generates the
### normalized counts, either by running the MCMC if we are doing spatial
### normalization or running a negative binomial regression if we are doing
### non spatial normalization

### The Nimble package is used to run the MCMC model. For Mac, this requires having
### an installation of Xcode for the C++ compilation to work.

library(here)
library(readr)

source(here("Code", "procedure_functions.R"))

filenum <- 1
input_files <- list.files(path=here("Processed Data", "MCMC Inputs"),
                          pattern="prepared_inputs_cluster[0-9]+_section[0-9]+",
                          full.names=T)
prepared_inputs <- read_rds(input_files[[filenum]])

cluster_id <- prepared_inputs$cluster
section_id <- prepared_inputs$section

num_chains <- 3

## Run negative binomial MCMC model with spatial random effect in Nimble
if (prepared_inputs$normalization_type=="spatial") {
  output <- run_nimble_code(code=nimble_code_covariates_nb, prepared_inputs=prepared_inputs,
                            num_iterations=50000, num_burnin=40000, n_thin=1, 
                            num_chains=num_chains, seed=prepared_inputs$seed[1:num_chains])
} else {
  output <- list(prepared_inputs=prepared_inputs)
}


est_params <- construct_est_params(output, spatial_model="NB", nonspatial_model="NB")

## Output estimated parameters
write_rds(est_params, file=here("Processed Data", "Estimated Parameters", 
                                paste0("est_params_cluster", cluster_id, "_section", section_id, ".rds")))
