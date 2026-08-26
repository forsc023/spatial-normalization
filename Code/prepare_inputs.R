
### This script runs a modified SpatialPCA procedure and prepares inputs for the MCMC
### procedure for one of the gene clusters.
### As an example, we are going to run this on Cluster 4.

library(here)
library(readr)

library(parallel)
library(mcprogress)
ncores <- 6

## Gene cluster for which we will prepare the inputs
clustnum <- 4

source(here("Code", "procedure_functions.R"))

rawcounts <- read_rds(here("Processed Data", "rawcounts.rds"))
partition_assignments <- read_rds(here("Processed Data", "partition_assignments.rds"))
all_gs_matrices <- read_rds(here("Processed Data", "all_gs_matrices.rds"))
clust_gs_matrices <- all_gs_matrices[[clustnum]] # GS matrices for cluster 4

## Get covariate values: sequencing depths, for each section
## These will be the same regardless of gene cluster.
all_seq_depths <- colSums(as.matrix(rawcounts))
seq_depths_list <- lapply(1:length(partition_assignments$overlapping_sections),
                          function(i) {
                            all_seq_depths[partition_assignments$overlapping_sections[[i]]]
                          })

## Prepare inputs
## We are using parallelization to run this more quickly on all the sections.
## Using 6 cores, took about 5 minutes on my computer
inputs_clust <- pmclapply(1:length(clust_gs_matrices), function(i) {
  curr_matrix <- as.matrix(clust_gs_matrices[[i]])
  curr_coords <- as.matrix(partition_assignments$overlapping_coordinates[[i]])
  curr_seq_depths <- seq_depths_list[[i]]
  curr_X <- cbind(rep(1, length(curr_seq_depths)), curr_seq_depths)
  
  
  curr_inputs <- prepare_inputs(gs_matrix=curr_matrix, covariates=curr_X,
                                coords=curr_coords, max_pc_number=5, seed=c(i, i+10000, i+20000))
  return(curr_inputs)
},
mc.cores=ncores)

for (i in 1:length(inputs_clust)) {
  inputs_clust[[i]]$cluster <- clustnum
  inputs_clust[[i]]$section <- unname(inputs_clust[[i]]$original_coords[,3][1])
}

## Output all sections if you ran all the sections
for (i in 1:length(inputs_clust)) {
  section <- inputs_clust[[i]]$section
  write_rds(inputs_clust[[i]], file=here("Processed Data", "MCMC Inputs", 
                                         paste0("prepared_inputs_cluster", clustnum, 
                                                "_section", section, ".rds")))
}

## Could also just run on a single section
#  section <- 1
#  curr_matrix <- as.matrix(clust_gs_matrices[[section]])
#  curr_coords <- as.matrix(partition_assignments$overlapping_coordinates[[section]])
#  curr_seq_depths <- seq_depths_list[[section]]
#  curr_X <- cbind(rep(1, length(curr_seq_depths)), curr_seq_depths)
#  curr_inputs <- prepare_inputs(gs_matrix=curr_matrix, covariates=curr_X,
#                                coords=curr_coords, max_pc_number=5, seed=c(section, section+10000, section+20000))
#  curr_inputs$cluster <- clustnum
#  curr_inputs$section <- section
# 
# ## And output single section
# write_rds(curr_inputs, file=here("Processed Data", "MCMC Inputs", 
#                                   paste0("prepared_inputs_cluster", clustnum, "_section", section, ".rds")))
