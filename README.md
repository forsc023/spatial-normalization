# spatial-normalization
This repository contains sample code for our spatial normalization procedure for spatial transcriptomics data. We use a public dataset from 10x Genomics, obtained the data from the 10x Genomics website: https://www.10xgenomics.com/datasets/visium-cytassist-gene-expression-libraries-of-post-xenium-human-colon-cancer-ffpe-using-the-human-whole-transcriptome-probe-set-2-standard 

The raw data obtained from 10x Genomics is located in the "Raw Data" folder. The scripts in the "Code" folder will read in the raw data, perform the procedure, and output to the "Processed Data" folder. Due to the computational complexity of the MCMC step of the procedure, the scripts are set up to run this step on one gene cluster and spatial section.

The scripts are described below and should be run in the following order:

1. preprocessing.R: this script reads in the raw gene-spot matrix and features, implements the Seurat normalization functions NormalizeData and SCTransform, and does a global filtering on the genes. The three output gene-spot matrices (raw, NormalizeData, and SCTransform counts) are organized to include the same genes & spots in the same order. Outputs are: (a)	rawcounts.rds (raw gene-spot matrix); (b) ndcounts.rds (NormalizeData gene-spot matrix in compressed format); (c) sctcounts.rds (SCTransform gene-spot matrix in compressed format); (d) coords_all_info.rds (coordinates data frame that also includes pixel and array coordinates); (e) coords.rds (coordinates matrix); (f) features.rds (dataframe of gene names, filtered to the genes we kept and in the same order as the outputted gene-spot matrices. The field       “sct_index” refers to the row on the outputted gene-spot matrices and “feature_index” refers to the row on the original unfiltered raw gene-spot matrix); (g) ndcounts_as_matrix.rds (NormalizeData gene-spot matrix in matrix format); (h) sct_kept_genes.rds (list of gene indices maintained in the SCTransform filtering)

2. clustering.Rmd: this is an R markdown file that calls Python via the reticulate package and runs a clustering algorithm from the       k_means_constrained Python package. It uses the NormalizeData (i.e. log-transformed and multiplied by a scaling factor) counts to create the clusters. The file paths will need to be changed to reflect the user’s computer. Outputs are:
   * nd_constrained_clusters.csv (list of cluster assignments for each gene.)

3. create_gs_matrices.R: partitions the tissue into unique & overlapping sections and then splits the full gene-spot matrix into smaller overlapping matrices using the overlapping sections and the gene clusters created in the previous script. Outputs are:
   * partition_assignments.rds (contains coordinates and spot indices for the unique and overlapping sections)
   * all_gs_matrices.rds (list of all the smaller gene-spot matrices, each corresponding to a single gene cluster + spatial section combination)

4.	prepare_inputs.R: runs modified SpatialPCA procedure and creates MCMC inputs for one of the gene clusters. This also determines whether spatial or non-spatial normalization will be run for each section with the given cluster. Outputs are:
   * Files with names of the form prepared_inputs_cluster[]_section[] (e.g., prepared_inputs_cluster4_section1.rds). Each contains the MCMC inputs for one of the clusters + sections.

5.	create_outputs.R: runs the MCMC procedure given a set of inputs for one cluster + section. Outputs are:
  a.	Files with names of the form est_params_cluster[]_section[] (e.g., est_params_cluster4_section1.rds). Each contains the estimated parameters and normalized counts for the indicated cluster and section.

6.	construct_normalized_gs.R: reads in the normalized gene-spot matrices for each cluster & section and combines them together into the final normalized gene-spot matrix. Outputs are:
  a.	spcounts.rds (final normalized gene-spot matrix)
