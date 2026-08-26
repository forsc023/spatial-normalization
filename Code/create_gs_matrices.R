### This script (1) splits the tissue sample into unique and overlapping sections, and
### (2) creates smaller gene-spot matrices from the full raw gene-spot matrix, with
### each smaller matrix containing the genes from one cluster and the spots from one
### of the overlapping sections.

library(here)
library(readr)
library(dplyr)

source(here("Code", "modified_spatialpca_functions.R"))
source(here("Code", "procedure_functions.R"))

## Read data
rawcounts <- read_rds(here("Processed Data", "rawcounts.rds"))


## Read in clusters
ndclust <- data.frame(read_csv(here("Processed Data", "nd_constrained_clusters.csv"))) %>%
  rename(cluster=X0)

real_coords <- read_rds(here("Processed Data", "coords.rds"))
spots <- as.matrix(real_coords)

## Create the overlapping sections
partition_assignments <- moving_window(spots=spots, x_subsections=9, 
                                       y_subsections=12, proportion_overlap=0.9)

partition_coords <- do.call(rbind, partition_assignments$unique_coordinates)
overlap_coords <- do.call(rbind, partition_assignments$overlapping_coordinates)

## Create g-s matrices
n_clusters <- length(unique(ndclust$cluster))
n_sections <- length(partition_assignments$overlapping_sections)
all_gs_matrices <- vector(mode="list", length=n_clusters)
for (i in 0:(n_clusters-1)) {
  clust_indices <- which(ndclust==i)
  rawcounts_clust <- rawcounts[clust_indices,]
  clust_gs_matrices <- lapply(1:length(partition_assignments$overlapping_sections), function(j) {
    rawcounts_clust[,partition_assignments$overlapping_sections[[j]]]
  })
  all_gs_matrices[[i+1]] <- clust_gs_matrices
}


## Output results
write_rds(partition_assignments, file=here("Processed Data", "partition_assignments.rds"))
write_rds(all_gs_matrices, file=here("Processed Data", "all_gs_matrices.rds"))