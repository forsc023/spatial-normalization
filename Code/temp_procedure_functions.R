library(readr)
library(spdep)
library(nimble)
library(coda)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(janitor)
library(berryFunctions) # for if.error function
library(Rfast)
library(pscl)
library(MASS)

st_directory <- "/projects/standard/zhangl4/forsc023/spatial_transcriptomics/"
source(paste0(st_directory, "functions/Modified SpatialPCA Functions for MSI.R"))


## Function to assign labels to the gene spot matrix rows/columns and the
## rows of the spatial coordinates matrix so that they match and we can
## enter them into SpatialPCA
label_matrix <- function(gs_matrix, covariates=NA, coords) {
  colnames(gs_matrix) <- paste0("S", c(1:ncol(gs_matrix)))
  rownames(gs_matrix) <- paste0("G", c(1:nrow(gs_matrix)))
  rownames(coords) <- paste0("S", c(1:nrow(coords)))
  
  if (!identical(covariates, NA)) {
    rownames(covariates) <- paste0("S", c(1:ncol(gs_matrix)))
  }
  
  result <- list(gs=gs_matrix, covariates=covariates, coords=coords)
  return(result)
}

## Function to filter out genes expressed in less than a given number
## (threshold_num) of spots:
filter_lowly_expressed_genes <- function(gs_matrix, threshold_num=5) {
  # For each gene, compute number of spots with nonzero counts
  num_expressed_locations <- rowSums(gs_matrix != 0)
  
  exclude_indices <- which(num_expressed_locations < threshold_num)
  include_indices <- which(num_expressed_locations >= threshold_num)
  
  new_gs <- gs_matrix[include_indices,]
  
  result <- list(new_gs=new_gs,
                 include_indices=include_indices, exclude_indices=exclude_indices)
  return(result)
}


## Function to define neighborhood structure for the purpose of computing Moran's I.
## This uses k-nearest neighbors, with the number of neighbors equal to the floor
## of sqrt(S), where S is the number of spots.
define_neighbors <- function(coords) {
  num_neighbors <- floor(sqrt(nrow(coords)))
  nbk_temp <- knearneigh(x = coords, k = num_neighbors)
  nbk <- knn2nb(nbk_temp)
  nbw <- nb2listw(nbk)
  return(nbw)
}

## Function to identify the most variable gene (MVG) in a gene-spot matrix.
## Takes in a gene-spot matrix with genes as the rows, spots as the columns.
## Returns (1) a vector with the indices of the genes ordered from most to least
## variable, and (2) a vector containing the counts of the most variable gene.
rank_gene_variances <- function(gs_matrix) {
  ## Turn gs matrix into a dataframe with each "variable" (column) being a gene
  gs_transpose <- data.frame(t(gs_matrix))
  # Compute variance of each column (gene)
  gene_variances <- sapply(gs_transpose, var) 
  n_genes <- nrow(gs_matrix)
  ## Returns indices of genes from most to least variable
  index_order <- order(gene_variances, decreasing=T)
  mvg_quantity <- min(10, n_genes)
  mvg_indices <- index_order[1:mvg_quantity]
  mvg_counts <- vector(mode="list", length=mvg_quantity)
  for (i in 1:length(mvg_counts)) {
    mvg_counts[[i]] <- gs_matrix[index_order[i],]
  }
  result <- list(index_order = index_order, 
                 mvg_indices = mvg_indices,
                 mvg_count_list = mvg_counts)
  return(result)
}

## Computes Moran's I for the most variable genes and takes the maximum
get_moran_stat <- function(mvg_count_list, nbw) {
  moran_tests <- lapply(mvg_count_list, moran.test, listw=nbw)
  max_moran_stat <- max(sapply(moran_tests, "[[", "estimate")["Moran I statistic",])
  return(max_moran_stat)
}

get_median_moran_stat <- function(mvg_count_list, nbw) {
  moran_tests <- lapply(mvg_count_list, moran.test, listw=nbw)
  median_moran_stat <- median(sapply(moran_tests, "[[", "estimate")["Moran I statistic",])
  return(median_moran_stat)
}

## Function to successively run SpatialPCA and compute Moran's I to determine
## the number of spatial PCs to use in the MCMC procedure.
determine_pc_number <- function(gs_matrix, covariates=NA, coords, mvg_indices, 
                                nbw, max_pc_number=5) {
  
  ## I want to set a ceiling on the number of spatial PCs we use even if we don't
  ## manage to get Moran's I down to the threshold, because the MCMC procedure
  ## will simply become unwieldy if there are too many PCs.
  ## (Plus, the SPCA procedure itself will start to take too long if we do it too
  ## many times.) This could be a 
  ## (somewhat arbitrary) number, like 10, or maybe a certain fraction of the 
  ## number of genes in the cluster. For now, let's try just setting a number.
  
  ## Making sure our max PC number is not larger than the number of genes in the
  ## cluster:
  max_pc <- min(nrow(gs_matrix), max_pc_number)
  
  moran_stats <- rep(NA, max_pc)
  #median_moran_stats <- rep(NA, max_pc)
  tau_values <- rep(NA, max_pc)
  sigma2_values <- rep(NA, max_pc)
  gamma_values <- rep(NA, max_pc)
  W_values <- vector(mode="list", length=max_pc)
  Z_values <- vector(mode="list", length=max_pc)
  B_values <- vector(mode="list", length=max_pc)
  X <- NA # covariates - we will fill in if covariates are entered to the function
  
  ## Build the sPCA kernel, which won't change when we change the number of PCs:
  if (identical(covariates, NA)) { # if there are no covariates
    spca_kernel <- CreateSpatialPCAObject_MODIFIED(counts=gs_matrix,
                                                   location=coords)
  } else {
    spca_kernel <- CreateSpatialPCAObject_MODIFIED(counts=gs_matrix, 
                                                   covariate=covariates, location=coords)
    X <- covariates
  }
  
  ## I keep running into a situation where the buildKernel function does not work
  ## because the matrix is too sparse to find the bandwidth (even after filtering
  ## the original GS matrix). So to prevent the whole thing from failing when this
  ## happens, I will use the if.error function, which allows us to return something
  ## (e.g. NA) if the buildKernel function throws an error and return the SPCA
  ## object otherwise. If we get an NA return from this, we'll proceed with
  ## non-spatial normalization.
  spca_kernel <- if.error(SpatialPCA_buildKernel(spca_kernel, kerneltype="gaussian",
                                                 bandwidthtype="SJ"), 
                          error_true=NA, 
                          error_false=SpatialPCA_buildKernel(spca_kernel, kerneltype="gaussian",
                                                             bandwidthtype="SJ"))
  #spca_kernel <- SpatialPCA_buildKernel(spca_kernel, kerneltype="gaussian", 
  #                                      bandwidthtype="SJ")
  if (identical(spca_kernel, NA)) {
    return(NA)
  }
  
  for (i in 1:max_pc) {
    pc_number_to_use <- i
    
    ## Compute SpatialPCA on the original g-s matrix
    spca <- SpatialPCA_EstimateLoading(spca_kernel, fast=F, SpatialPCnum=pc_number_to_use)
    spca <- SpatialPCA_SpatialPCs(spca, fast=F)
    
    if (!identical(covariates, NA)) { # if there are covariates
      ## B hat, using formula from supplement of SpatialPCA paper
      B_hat <- solve(spca@params$XtX)%*%t(covariates)%*%(t(gs_matrix) - t(spca@SpatialPCs)%*%t(spca@W))
      B_values[[i]] <- B_hat
      ## Subtract (XB)^T + WZ from original matrix
      new_gs_matrix <- gs_matrix - t(covariates%*%B_hat) - spca@W%*%spca@SpatialPCs
    } else {
      ## Subtract WZ from original matrix
      new_gs_matrix <- gs_matrix - spca@W%*%spca@SpatialPCs
    }
    
    ## New counts with which to compute Moran's I:
    ## We will pull out the rows of the gene-spot matrix corresponding to the
    ## most variable genes that we determined earlier, and then make a list
    ## where each element of the list is a row of this reduced matrix
    ## (i.e., each element of the list is the counts of one of the most
    ## variable genes.)
    new_gs_matrix_mvg <- new_gs_matrix[mvg_indices,]
    new_mvg_count_list <- lapply(seq_len(nrow(new_gs_matrix_mvg)), 
                                 function(i) new_gs_matrix_mvg[i,])
    
    
    ## Compute Moran's I
    curr_moran_stat <- get_moran_stat(mvg_count_list = new_mvg_count_list, 
                                      nbw = nbw)
    #curr_median_moran_stat <- get_median_moran_stat(mvg_count_list = new_mvg_count_list,
    #                                                nbw = nbw)
    moran_stats[i] <- curr_moran_stat
    #median_moran_stats[i] <- curr_median_moran_stat
    tau_values[i] <- spca@tau
    sigma2_values[i] <- spca@sigma2_0
    gamma_values[i] <- spca@bandwidth
    W_values[[i]] <- spca@W
    Z_values[[i]] <- spca@SpatialPCs
    
    ## Break out of the for loop if Moran's I is low enough
    ## Or if the difference between the previous statistic and the current 
    ## statistic is < 0.2
    if (curr_moran_stat < 0.3) {
      break
    }
    
    if (i > 1) {
      if (abs(moran_stats[i-1] - curr_moran_stat) < 0.02) {
        break
      }
    }
  }
  
  ## Get rid of the NAs at the end of our output lists, which will occur if
  ## we achieved a low enough Moran's I and broke out of the loop early.
  tau_values <- tau_values[!is.na(tau_values)]
  sigma2_values <- sigma2_values[!is.na(sigma2_values)]
  gamma_values <- gamma_values[!is.na(gamma_values)]
  moran_stats <- moran_stats[!is.na(moran_stats)]
  #median_moran_stats <- median_moran_stats[!is.na(median_moran_stats)]
  W_values <- Filter(Negate(is.null), W_values)
  Z_values <- Filter(Negate(is.null), Z_values)
  B_values <- Filter(Negate(is.null), B_values)
  
  result <- list(pc_number_to_use=pc_number_to_use, moran_stats=moran_stats,
                 tau_values=tau_values, sigma2_values=sigma2_values,
                 gamma_values=gamma_values, W_values=W_values, Z_values=Z_values,
                 B_values=B_values, X=X)
  return(result)
}

## Gaussian covariance function
gaussian_cov <- function(dist, tau, gamma) {
  num_spots <- dim(dist)[1]
  result <- matrix(data=NA, nrow=num_spots, ncol=num_spots)
  for (i in 1:num_spots) {
    for (j in 1:num_spots) {
      result[i,j] <- tau*exp(-1*(dist[i,j]^2)/gamma)
    }
  }
  return(result)
}


#==============================================================================#
#==============================================================================#
#==============================================================================#


### Function that calls the other functions in order to compile the inputs we'll
### need for the MCMC procedure. Takes in a gene-spot matrix and spatial coordinates.
### This prepares all the inputs for one cluster and one section of the map.
### The seed will be used to generate rnorm() random initial values for B in the
### case that non-spatial normalization is selected. Otherwise, we will determine
### initial values for B from the SpatialPCA output.
prepare_inputs <- function(gs_matrix, coords, covariates=NA, filter_threshold=5,
                           max_pc_number=5, seed) {
  original_coords <- coords
  coords_xy <- coords[,c("x", "y")]
  coords_xy_scaled <- scale(coords_xy)
  coords[,"x"] <- coords_xy_scaled[,"x"]
  coords[,"y"] <- coords_xy_scaled[,"y"]
  
  labeled_data <- label_matrix(gs_matrix=gs_matrix, covariates=covariates, coords=coords) # label
  gs <- labeled_data$gs
  coords <- labeled_data$coords
  covariates <- labeled_data$covariates
  nbw <- define_neighbors(coords=coords) # define neighbors
  filtered_gs <- filter_lowly_expressed_genes(gs_matrix=gs, threshold_num=filter_threshold)$new_gs # filter
  
  ## If the filtered matrix has only 1 row, nrow() will return NULL, so we need
  ## to account for this in our if-else statement:
  if (is.null(nrow(filtered_gs))) {
    result <- list(filtered_gs=filtered_gs, original_gs=gs_matrix,
                   scaled_coords=coords, original_coords=original_coords,
                   nbw=nbw, normalization_type="none")
    return(result)
  }
  if (nrow(filtered_gs)==0) {
    result <- list(filtered_gs=filtered_gs, original_gs=gs_matrix,
                   scaled_coords=coords, original_coords=original_coords,
                   nbw=nbw, normalization_type="none")
    return(result)
  }
  gene_vars <- apply(filtered_gs, 1, var)
  ## Take the maximum of 1 and the largest empirical gene variance to use as a
  ## constant (the max of the prior on sigma) in the mcmc. This is to protect
  ## against artificially restricting the variance parameter *too* much in the
  ## case of very small empirical gene variances.
  max_gene_var <- max(1, max(gene_vars))
  ## Taking the log before identifying MVG and computing Moran's I, to be
  ## consistent with the next step where we will run SpatialPCA and compute
  ## Moran's I several times using the log counts.
  loggs <- log(filtered_gs + 0.1)
  
  g <- nrow(filtered_gs) # number of genes
  S <- ncol(filtered_gs) # number of spots
  y <- filtered_gs # Gene-spot matrix
  X <- covariates # covariates (sequencing depth and intercept)
  X_t <- t(X)
  
  variance_ranking <- rank_gene_variances(loggs)
  mvg_indices <- variance_ranking$mvg_indices
  mvg_count_list <- variance_ranking$mvg_count_list
  moran_stat <- get_moran_stat(mvg_count_list = mvg_count_list, nbw = nbw)
  
  ## If initial Moran's I is high enough, prepare inputs for spatial
  ## normalization
  if (moran_stat >= 0.3) {
    pc_determination <- determine_pc_number(gs_matrix=loggs, covariates=covariates,
                                            coords=coords,
                                            mvg_indices=mvg_indices, nbw=nbw,
                                            max_pc_number=max_pc_number)
    ## If the matrix wasn't too sparse to compute the bandwidth and we got a
    ## return from pc_determination, proceed:
    #if (!is.na(pc_determination)) {
    if (!identical(pc_determination, NA)) {
      pc_num <- pc_determination$pc_number_to_use
      tau_to_use <- pc_determination$tau_values[length(pc_determination$tau_values)] * pc_determination$sigma2_values[length(pc_determination$sigma2_values)]
      gamma_to_use <- pc_determination$gamma_values[length(pc_determination$gamma_values)]
      W_to_use <- pc_determination$W_values[[length(pc_determination$W_values)]]
      Z_to_use <- pc_determination$Z_values[[length(pc_determination$Z_values)]]
      if (!identical(covariates, NA)) { # if there are covariates
        B <- pc_determination$B_values[[length(pc_determination$B_values)]]
        p <- ncol(covariates) # dimension of beta (2 if we have intercept & seq depth)
      } else {
        B <- NA
        p <- NA
      }
      
      ## Assign data, constants, and initial values to put into the MCMC
      dist <- as.matrix(dist(coords))
      Sigma <- gaussian_cov(dist=dist, tau=tau_to_use, gamma=gamma_to_use)
      d <- pc_num # number of PCS, based on our procedure above
      
      Z_init <- Z_to_use
      W_init <- W_to_use
      B_t_init <- t(B)
      r_init <- rep(1, g)
      
      max_Z_var <- max(1, max(apply(Z_init, 1, var)))
      constants <- list(X_t=X_t, S=S, d=d, g=g, p=p, zero=rep(0, S), Sigma=Sigma,
                        max_gene_var=max_gene_var, max_Z_var=max_Z_var)
      data <- list(y=y)
      inits <- list(W=W_init, Z=Z_init, B_t=B_t_init, r=r_init)
      
      ## Note in our returns:
      ## filtered_gs is the gs matrix after filtering out lowly expressed genes
      ## original_gs is the original count matrix entered into the function
      ## log_gs is the log transformed count matrix, which is what is used in all
      ## the subsequent calculations (Moran's I, spatial PC estimation)
      result <- list(filtered_gs=filtered_gs, original_gs=gs_matrix, log_gs=loggs,
                     max_pc_number = max_pc_number,
                     scaled_coords=coords, original_coords=original_coords, nbw=nbw,
                     mvg_indices=mvg_indices,
                     mvg_count_list=mvg_count_list, initial_moran=moran_stat,
                     pc_determination=pc_determination, pc_num=pc_num,
                     tau_to_use=tau_to_use, gamma_to_use=gamma_to_use,
                     W_to_use=W_to_use, Z_to_use=Z_to_use, B_to_use=B, p=p,
                     data=data, constants=constants, inits=inits, seed=seed,
                     normalization_type="spatial")
      return(result)
      
    }
  }
  ### If initial Moran's I not high enough, return just non-spatial inputs to
  ### put into the non-spatial MCMC
  ## Instead of getting the initial values of B from SpatialPCA, we will
  ## randomly generate them. Dimension of B should be pxg
  if (!identical(covariates, NA)) { # if there are covariates
    p <- ncol(covariates) # dimension of beta (2 if we have intercept & seq depth)
    set.seed(seed)
    #B <- matrix(nrow=p, ncol=g, data=rnorm(n=p*g, mean=0, sd=1), byrow=T)
    # Initialize the intercept to standard normal r.v.'s and
    # initialize the coefficients for sequencing depth to be 1/sequencing depth.
    B <- matrix(nrow=p, ncol=g, data=c(rnorm(g, mean=0, sd=1),
                                       rep(1/max(covariates[,2]), g),
                                       rnorm((p-2)*g, 0, 1)),
                byrow=T)
    
  } else {
    B <- NA
  }
  
  B_t_init <- t(B)
  constants <- list(X_t=X_t, S=S, g=g, p=p, max_gene_var=max_gene_var)
  data <- list(y=y)
  inits <- list(B_t=B_t_init)
  
  result <- list(filtered_gs=filtered_gs, original_gs=gs_matrix, log_gs=loggs,
                 scaled_coords=coords, original_coords=original_coords,
                 nbw=nbw, mvg_indices=mvg_indices,
                 mvg_count_list=mvg_count_list, initial_moran=moran_stat,
                 B_to_use=B, p=p, data=data, constants=constants, inits=inits,
                 seed=seed, normalization_type="nonspatial")
  
  return(result)
  
}


#==============================================================================#
#==============================================================================#
#==============================================================================#

## Function to run the MCMC in Nimble
run_nimble_code <- function(code, prepared_inputs, num_iterations, num_burnin,
                            n_thin=100,
                            num_chains=3, seed) {
  library(nimble)
  data <- prepared_inputs$data
  constants <- prepared_inputs$constants
  inits <- prepared_inputs$inits
  
  poisModel <- nimbleModel(code, constants = constants, data = data, inits = inits)
  cpoisModel <- compileNimble(poisModel)
  
  # monitors
  # If we are not doing spatial normalization:
  if (prepared_inputs$normalization_type=="nonspatial") {
    poisconfMC <- configureMCMC(poisModel, monitors = c("B_t", "sigma"))
    # If spatial normalization without covariates
  } else if (identical(inits$B_t, NA)) {
    poisconfMC <- configureMCMC(poisModel, monitors = c("Z", "W"))
    # Otherwise, spatial normalization with covariates
  } else {
    poisconfMC <- configureMCMC(poisModel, monitors = c("Z", "W", "B_t", "r"))
    
    #poisconfMC$removeSamplers(c("Z", "W"))
    #poisconfMC$addSampler(
    #  target = c("Z"), type="RW_block", control=list(scale=0.01, adaptive=T, adaptInterval=20)
    #)
    #poisconfMC$addSampler(
    #  target = c("W"), type="RW_block", control=list(adaptInterval=20)
    #)
    
    #poisconfMC$addSampler(
    #  target = c("Z", "W"), type="AF_slice", control=list(adaptive=T)
    #)
    #poisconfMC$addSampler(
    #  target = "B_t", type="AF_slice", control=list(adaptive=T)
    #)
    #poisconfMC$printSamplers()
  }
  
  poisMCMC <- buildMCMC(poisconfMC)
  cpoisMCMC <- compileNimble(poisMCMC, project = cpoisModel)
  
  ## Run MCMC
  MCMC_out <- runMCMC(cpoisMCMC, niter=num_iterations, nburnin=num_burnin,
                      thin=n_thin,
                      nchains=num_chains, samplesAsCodaMCMC = TRUE, setSeed=seed)
  return(list(MCMC_out=MCMC_out,
              prepared_inputs=prepared_inputs,
              num_iterations=num_iterations,
              num_burnin=num_burnin,
              n_thin=n_thin,
              num_chains=num_chains,
              seed_used=seed))
}

## Nimble code without covariates -- written while developing the method and not
## intended to be used, since our method involves regressing out the impact of
## sequencing depth.  
nimble_code <- nimbleCode({
  # For each gene i
  for (i in 1:g) {
    # For each spot j
    for (j in 1:S) {
      y[i,j] ~ dpois(lambda=mu[i,j])
      # The WZ term gives the dot product of the ith row of W and the jth column of Z
      # Doing sum(W[i,1:d]*Z[1:d,j]) gives the same result as doing matrix multiplication
      # Idea is that this gives us the i,jth entry of Y, the gxS count matrix
      
      ## Nimble gets upset when trying to compute the inner product when there is 
      ## only one spatial PC because then W[i,1:d] and Z[1:d,j] are scalars. So
      ## we will include an if/else statement to deal with this.
      if (d==1) {
        log(mu[i,j]) <- W[i,1]*Z[1,j]
      } else {
        log(mu[i,j]) <- sum(W[i,1:d]*Z[1:d,j])
      }
      
    }
  }
  
  # Assign priors for Z, using the spatial covariance matrix that we'll enter
  # as a constant.
  for (k in 1:d) { # for each spatial PC
    Z[k,1:S] ~ dmnorm(mean=zero[1:S], cov=Sigma[1:S,1:S])
  }
  
  # Assign priors for W 
  for (i in 1:g) {
    for (k in 1:d) {
      W[i,k] ~ dnorm(mean=0, sd=1)
    }
  }
})


## Nimble code including the covariates
nimble_code_covariates <- nimbleCode({
  # For each gene i
  for (i in 1:g) {
    # For each spot j
    for (j in 1:S) {
      y[i,j] ~ dpois(lambda=mu[i,j])
      # The WZ term gives the dot product of the ith row of W and the jth column of Z
      # Doing sum(W[i,1:d]*Z[1:d,j]) gives the same result as doing matrix multiplication
      # Idea is that this gives us the i,jth entry of Y, the gxS count matrix
      
      ## Nimble gets upset when trying to compute the inner product when there is 
      ## only one spatial PC because then W[i,1:d] and Z[1:d,j] are scalars. So
      ## we will include an if/else statement to deal with this.
      if (d==1) {
        ## p = number of covariates
        log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j]) + W[i,1]*Z[1,j]
      } else {
        log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j]) + sum(W[i,1:d]*Z[1:d,j])
      }
      
    }
  }
  
  # Assign priors for Z, using the spatial covariance matrix that we'll enter
  # as a constant.
  for (k in 1:d) { # for each spatial PC
    Z[k,1:S] ~ dmnorm(mean=zero[1:S], cov=Sigma[1:S,1:S])
  }
  
  # Assign priors for W 
  for (i in 1:g) {
    for (k in 1:d) {
      W[i,k] ~ dnorm(mean=0, sd=1)
    }
  }
  
  # Priors for B 
  for (i in 1:g) {
    for (l in 1:p) {
      B_t[i,l] ~ dnorm(mean=0, sd=5)
      #B_t[i,l] ~ dnorm(mean=0, sd=100)
    }
  }
})



nimble_code_covariates_nb <- nimbleCode({
  for (i in 1:g) {
    for (j in 1:S) {
      ## Negative binomial likelihood with overdispersion parameter r
      ## E[y] = mu = r(1-p)/p
      ## p = r/(mu + r)
      y[i,j] ~ dnegbin(prob = pi[i,j], size = r[i])
      
      if (d==1) {
        ## p = number of covariates
        log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j]) + W[i,1]*Z[1,j]
      } else {
        log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j]) + sum(W[i,1:d]*Z[1:d,j])
      }
      pi[i,j] <- r[i]/(mu[i,j] + r[i])
    }
  }
  
  # Assign priors for Z, using the spatial covariance matrix that we'll enter
  # as a constant.
  for (k in 1:d) { # for each spatial PC
    Z[k,1:S] ~ dmnorm(mean=zero[1:S], cov=Sigma[1:S,1:S])
  }
  
  # Assign priors for W 
  for (i in 1:g) {
    for (k in 1:d) {
      W[i,k] ~ dnorm(mean=0, sd=1)
    }
  }
  
  ## Priors for beta
  for (i in 1:g) {
    for (l in 1:p) {
      B_t[i,l] ~ dnorm(mean=0, sd=5)
    }
  }
  ## Prior for overdispersion parameter
  for (i in 1:g) {
    #r[i] ~ dunif(0,1000)
    r[i] ~ dinvgamma(shape=2, scale=1)
  }
  
})


## Nimble code without the spatial PCs, so we are just estimating the coefficients
## associated with sequencing depth without using the spatial information.
## This then becomes a simple poisson model.
nimble_code_covariates_nonspatial <- nimbleCode({
  # For each gene i
  for (i in 1:g) {
    # For each spot j
    for (j in 1:S) {
      y[i,j] ~ dpois(lambda=mu[i,j])
      log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j]) + theta[i,j]
      theta[i,j] ~ dnorm(mean=0, sd=sigma[i])
    }
    ## each gene has its own error variance?
    sigma[i] ~ dunif(0, 3*max_gene_var)
    #sigma[i] ~ dinvgamma(shape=3, scale=1)
  }
  # Priors for B 
  for (i in 1:g) {
    for (l in 1:p) {
      B_t[i,l] ~ dnorm(mean=0, sd=5)
    }
  }
})

nimble_code_covariates_nonspatial_nb <- nimbleCode({
  for (i in 1:g) {
    for (j in 1:S) {
      ## Negative binomial likelihood with overdispersion parameter r
      ## E[y] = mu = r(1-p)/p
      ## p = r/(mu + r)
      y[i,j] ~ dnegbin(prob = pi[i,j], size = r[i])
      log(mu[i,j]) <- sum(B_t[i,1:p]*X_t[1:p,j])
      pi[i,j] <- r[i]/(mu[i,j] + r[i])
    }
  }
  ## Priors for beta
  for (i in 1:g) {
    for (l in 1:p) {
      B_t[i,l] ~ dnorm(mean=0, sd=5)
    }
  }
  ## Prior for overdispersion parameter
  for (i in 1:g) {
    #r[i] ~ dunif(0,1000)
    r[i] ~ dinvgamma(shape=2, scale=1)
  }
  
})



#==============================================================================#
#==============================================================================#
#==============================================================================#

### ANALYZING RESULTS

## Function to plot spatial PCs. i.e., given a matrix (z) and coordinates (coords)
## with the number of rows in "coords" matching the number of columns in the matrix,
## the function will spatially plot each row of the matrix.
plot_spatialpcs <- function(z, coords, labs=T) {
  z_plots <- vector(mode="list", length=nrow(z))
  for (i in 1:length(z_plots)) {
    curr_pc <- z[i,]
    curr_plotdata <- cbind(coords, curr_pc)
    num_spots <- nrow(coords)
    dot_size <- floor(100/sqrt(num_spots))
    
    curr_plot <- ggplot(curr_plotdata, aes(x=xcoords, y=ycoords, color=curr_pc)) +
      geom_point(size=dot_size) + theme_bw() + 
      #scale_color_continuous(type="viridis", name=paste0("Spatial PC ", i)) +
      #labs(title=paste0("Spatial PC ", i))
      scale_color_continuous(type="viridis") 
    
    if (labs) {
      curr_plot <- curr_plot + 
        labs(title=paste0("Spatial PC ", i), color=paste0("Spatial PC ", i))
    } 
    
    z_plots[[i]] <- curr_plot
  }
  return(z_plots)
}


### Negative binomial regression, to be used for non spatially normalized sections.
### Runs a separate regression for each gene and calculates the fitted mean and
### variance (using the overdispersion parameter)
run_nb_regression <- function(prepared_inputs, method="NB") {
  B_t_est <- matrix(nrow=prepared_inputs$constants$g, ncol=2, data=NA)
  mu_hat <- matrix(nrow=prepared_inputs$constants$g, ncol=prepared_inputs$constants$S,
                   data=NA)
  marginal_variance <- matrix(nrow=prepared_inputs$constants$g, ncol=prepared_inputs$constants$S,
                              data=NA)
  
  for (g in 1:prepared_inputs$constants$g) {
    curr_y <- prepared_inputs$filtered_gs[g,]
    curr_seq_depths <- prepared_inputs$constants$X_t[2,]
    
    if (method=="NB") {
      ## Check if we are able to run negative binomial regression for this gene. If
      ## we are not, e.g. due to too much sparsity, run quasipoisson regression instead.
      model_use <- if.error(glm.nb(curr_y ~ curr_seq_depths), 
                            error_true=list(model_use=qpois.reg(x=curr_seq_depths, y=curr_y),
                                            method="QP"),
                            error_false = list(model_use=glm.nb(curr_y ~ curr_seq_depths),
                                               method="NB"))
      if (model_use$method=="NB") { # use negative binominal regression
        nb_reg <- model_use$model_use
        B_t_est[g,] <- nb_reg$coefficients
        log_mu_hat <- nb_reg$coefficients["(Intercept)"] + nb_reg$coefficients["curr_seq_depths"]*curr_seq_depths
        curr_mu_hat <- exp(log_mu_hat)
        mu_hat[g,] <- curr_mu_hat
        marginal_variance[g,] <- curr_mu_hat + (curr_mu_hat^2)/nb_reg$theta
      } else { # use the quasipoisson regression
        ## Old code for quasi poisson regression to find coefficients and overdispersion
        ## parameter. We have Y ~ (mu, phi*mu),
        ## log(mu) = beta0 + beta1*seq_depth
        qpois_reg <- model_use$model_use
        B_t_est[g,] <- qpois_reg$be
        log_mu_hat <- qpois_reg$be["(Intercept)",] + qpois_reg$be["x",]*curr_seq_depths
        curr_mu_hat <- exp(log_mu_hat)
        mu_hat[g,] <- curr_mu_hat
        marginal_variance[g,] <- qpois_reg$phi*curr_mu_hat
      }

    }
    if (method=="ZINB") {
      ## Zero inflated negative binomial regression
      ## Y = (1-Z)X where X ~ NB(theta, mu), Z ~ Ber(pi), E(X) = mu, E(Z) = pi,
      ## Var(X) = mu + mu^2/theta, Var(Z) = pi(1-pi)
      ## E(Y) = mu(1-pi) where pi is the probability of being in the "zero" group (first process)
      ## Var(Y) = mu(1-pi)(1 + mu(pi + 1/theta))
      zinb_reg <- zeroinfl(curr_y ~ curr_seq_depths | curr_seq_depths, dist="negbin", 
                           link="logit")
      B_t_est[g,] <- zinb_reg$coefficients$count
      curr_mu_hat <- exp(zinb_reg$coefficients$count["(Intercept)"] + zinb_reg$coefficients$count["curr_seq_depths"]*curr_seq_depths)
      pi_hat <- expit(zinb_reg$coefficients$zero["(Intercept)"] + 
                        zinb_reg$coefficients$zero["curr_seq_depths"]*curr_seq_depths)
      fitted_vals <- curr_mu_hat*(1-pi_hat)
      mu_hat[g,] <- fitted_vals
      alpha_hat <- 1/zinb_reg$theta
      marginal_variance[g,] <- curr_mu_hat*(1-pi_hat)*(1 + curr_mu_hat*(pi_hat + alpha_hat))
    }
    
  }
  result <- list(mu_hat=mu_hat,
                 marginal_variance=marginal_variance,
                 B_t_est=B_t_est)
}

## Alt version that smooths the overdispersion parameter estimates as a function
## of the geometric mean of the gene expression
run_nb_regression_smoothed <- function(prepared_inputs, method="NB") {
  B_t_est <- matrix(nrow=prepared_inputs$constants$g, ncol=2, data=NA)
  mu_hat <- matrix(nrow=prepared_inputs$constants$g, ncol=prepared_inputs$constants$S,
                   data=NA)
  marginal_variance <- matrix(nrow=prepared_inputs$constants$g, ncol=prepared_inputs$constants$S,
                              data=NA)
  
  for (g in 1:prepared_inputs$constants$g) {
    curr_y <- prepared_inputs$filtered_gs[g,]
    curr_seq_depths <- prepared_inputs$constants$X_t[2,]
    
    if (method=="NB") {
      ## Gene geometric means
      geometric_means <- sapply(1:nrow(prepared_inputs$filtered_gs), function(gene) {
        gene_counts <- prepared_inputs$filtered_gs[gene,]
        gene_gmean <- exp(mean(log(gene_counts + 1))) - 1
        return(gene_gmean)
      }) 
      
      ## Get estimates from negative binomial regression. 
      ## If the negative binomial regression can't be run (i.e. due to sparsity),
      ## use the parameter estimates from a quasi poisson regression for the mean
      ## estimation, and we will use a Loess smoothing to estimate the theta parameters
      ## as a function of the gene geometric mean. Any missing thetas from the NB regression
      ## not being able to be run can then be filled in using the predictions from the
      ## smoothing.
      raw_reg_estimates <- lapply(1:nrow(prepared_inputs$filtered_gs), function(gene) {
        curr_y <- prepared_inputs$filtered_gs[gene,]
        curr_seq_depths <- prepared_inputs$constants$X_t[2,]
        
        model_use <- if.error(glm.nb(curr_y ~ curr_seq_depths), 
                              error_true=list(model_use=qpois.reg(x=curr_seq_depths, y=curr_y),
                                              reg_method="QP"),
                              error_false = list(model_use=glm.nb(curr_y ~ curr_seq_depths),
                                                 reg_method="NB"))
        
        if (model_use$reg_method=="NB") { # use negative binominal regression
          nb_reg <- model_use$model_use
          reg_coefficients <- nb_reg$coefficients
          log_mu_hat <- nb_reg$coefficients["(Intercept)"] + nb_reg$coefficients["curr_seq_depths"]*curr_seq_depths
          curr_mu_hat <- exp(log_mu_hat)
          theta <- nb_reg$theta
        } else {
          qpois_reg <- model_use$model_use
          reg_coefficients <- qpois_reg$be
          log_mu_hat <- qpois_reg$be["(Intercept)",] + qpois_reg$be["x",]*curr_seq_depths
          curr_mu_hat <- exp(log_mu_hat)
          theta <- NA
        }
        result <- list(
          mu_hat = curr_mu_hat,
          theta = theta,
          reg_coefficients = reg_coefficients
        )
        return(result)
      })
      
      theta_estimates <- sapply(raw_reg_estimates, "[[", "theta")
      gmean_data <- data.frame(
        gmean = geometric_means,
        theta = theta_estimates
      ) %>%
        filter(!is.na(theta)) %>%
        arrange(gmean)
      
      theta_loess <- loess(theta ~ gmean, data=gmean_data)
      theta_smoothed <- predict(theta_loess, newdata=data.frame(gmean=geometric_means))
      
      for (gene in 1:nrow(prepared_inputs$filtered_gs)) {
        curr_mu_hat <- raw_reg_estimates[[gene]]$mu_hat
        mu_hat[gene,] <- curr_mu_hat
        B_t_est[gene,] <- raw_reg_estimates[[gene]]$reg_coefficients
        if (is.na(theta_smoothed[gene])) {
          theta <- 1
        } else {
          theta <- theta_smoothed[gene]
        }
        marginal_variance[g,] <- curr_mu_hat + (curr_mu_hat^2)/theta
      }

    } # end of "if method==NB"
    if (method=="ZINB") {
      ## Zero inflated negative binomial regression
      ## Y = (1-Z)X where X ~ NB(theta, mu), Z ~ Ber(pi), E(X) = mu, E(Z) = pi,
      ## Var(X) = mu + mu^2/theta, Var(Z) = pi(1-pi)
      ## E(Y) = mu(1-pi) where pi is the probability of being in the "zero" group (first process)
      ## Var(Y) = mu(1-pi)(1 + mu(pi + 1/theta))
      zinb_reg <- zeroinfl(curr_y ~ curr_seq_depths | curr_seq_depths, dist="negbin", 
                           link="logit")
      B_t_est[g,] <- zinb_reg$coefficients$count
      curr_mu_hat <- exp(zinb_reg$coefficients$count["(Intercept)"] + zinb_reg$coefficients$count["curr_seq_depths"]*curr_seq_depths)
      pi_hat <- expit(zinb_reg$coefficients$zero["(Intercept)"] + 
                        zinb_reg$coefficients$zero["curr_seq_depths"]*curr_seq_depths)
      fitted_vals <- curr_mu_hat*(1-pi_hat)
      mu_hat[g,] <- fitted_vals
      alpha_hat <- 1/zinb_reg$theta
      marginal_variance[g,] <- curr_mu_hat*(1-pi_hat)*(1 + curr_mu_hat*(pi_hat + alpha_hat))
    }
    
  }
  result <- list(mu_hat=mu_hat,
                 marginal_variance=marginal_variance,
                 B_t_est=B_t_est)
}

## Function that constructs estimated matrices from the output that comes out of
## run_nimble_code: For genes that were filtered out of the current section before
## normalization, we can (1) input zero counts (filter_method=="zero") or
## (2) input the raw counts (filter_method=="raw").
## spatial_model refers to the MCMC model for spatial normalization: poisson or negative binomial.
## nonspatial_model refers to the regression model for non-spatial normalization: negative binomial or zero inflated negative binomial.
construct_est_params <- function(output, filter_method="zero", spatial_model="NB", nonspatial_model="NB") {
  coords <- output$prepared_inputs$original_coords %>%
    data.frame() %>%
    rename(xcoords=x, ycoords=y)
  normalization_type <- output$prepared_inputs$normalization_type
  
  #coords <- output$prepared_inputs$coords %>%
  #  data.frame() %>%
  #  rename(xcoords=x, ycoords=y)
  
  ## First, if we didn't perform any normalization due to sparsity,
  ## just input everything as 0's
  if (output$prepared_inputs$normalization_type=="none") {
    Y_original <- output$prepared_inputs$original_gs
    Y_original_filtered <- output$prepared_inputs$filtered_gs
    Y_normalized <- matrix(nrow=nrow(Y_original), ncol=ncol(Y_original), data=0)
    return(list(Y_original=Y_original,
                Y_original_filtered=Y_original_filtered,
                Y_normalized=Y_normalized,
                coords=coords,
                normalization_type=normalization_type))
  }
  
  X_t <- output$prepared_inputs$constants$X_t
  Y_true <- output$prepared_inputs$filtered_gs
  #Y_log <- output$prepared_inputs$log_gs
  
  n_genes <- nrow(Y_true)
  n_spots <- ncol(Y_true)
  
  if (output$prepared_inputs$normalization_type=="nonspatial") {
    regression_output <- run_nb_regression(output$prepared_inputs, method=nonspatial_model)
    mu_hat <- regression_output$mu_hat
    marginal_variance <- regression_output$marginal_variance
    Y_normalized_filtered <- (Y_true - mu_hat)/sqrt(marginal_variance)
    B_t_est <- regression_output$B_t_est
    Z_est <- NA
    W_est <- NA
    r_est <- NA
  }
  
  if (output$prepared_inputs$normalization_type=="spatial") {
    mcmc_output <- output$MCMC_out
    post_estimates <- data.frame(as.matrix(mcmc_output,
                                           chains=T)) %>%
      clean_names() %>%
      dplyr::select(-any_of(c("chain"))) %>%
      pivot_longer(cols = everything(), names_to="parameter", values_to="value") %>%
      mutate(parameter = gsub("b_t", "bt", parameter)) %>%
      group_by(parameter) %>%
      ## Get posterior means and 95% credible intervals
      summarize(mean = mean(value),
                low = quantile(value, 0.025),
                high = quantile(value, 0.975)) %>%
      data.frame()
    
    ## Reconstruct the B, W, and Z matrices
    estimates_matrices <- post_estimates %>%
      #filter(!str_detect(parameter, "^r")) %>%
      dplyr::select(parameter, mean) %>%
      separate(col=parameter, into=c("parameter", "row", "col"), sep="_") %>%
      arrange(parameter, as.numeric(col), as.numeric(row)) %>%
      pivot_wider(names_from=col, values_from=mean, names_prefix="col")
    
    B_t_est <- estimates_matrices %>%
      filter(parameter=="bt") %>%
      ## Select column if not all NAs
      select_if(function(column) sum(is.na(column)) == 0) %>%
      dplyr::select(-parameter, -row) %>%
      as.matrix()
    
    B_est <- t(B_t_est)
    
    sigma_est <- estimates_matrices %>%
      filter(parameter=="sigma") %>%
      select_if(function(column) sum(is.na(column)) == 0) %>%
      dplyr::select(-parameter, -row) %>%
      as.matrix() %>% c()
    
    # over dispersion parameter, if we used the negative binomial model
    r_est <- estimates_matrices %>%
      filter(parameter=="r") %>%
      select_if(function(column) sum(is.na(column))==0) %>%
      dplyr::select(-parameter, -row) %>%
      as.matrix()
    
    Z_est <- estimates_matrices %>%
      filter(parameter=="z") %>%
      select_if(function(column) sum(is.na(column)) == 0) %>%
      dplyr::select(-parameter, -row) %>%
      as.matrix()
    
    W_est <- estimates_matrices %>%
      filter(parameter=="w") %>%
      select_if(function(column) sum(is.na(column)) == 0) %>%
      dplyr::select(-parameter, -row) %>%
      as.matrix()
    
    marginal_variance <- matrix(nrow=n_genes, ncol=n_spots, data=NA) # Pois mixed model
    nb_marginal_variance <- matrix(nrow=n_genes, ncol=n_spots, data=NA) # NegBin mixed model
    tau <- output$prepared_inputs$tau_to_use
    
    for (g in 1:n_genes) {
      beta0_g <- B_est[1,g]
      beta1_g <- B_est[2,g]
      w_sq_sum <- sum(W_est[g,]^2)
      for (s in 1:n_spots) {
        x_s <- X_t[2,s]
        # Poisson marginal variance
        marginal_variance[g,s] <- exp(2*(beta0_g + beta1_g*x_s) + tau*w_sq_sum)*(exp(tau*w_sq_sum) - 1) + exp(beta0_g + beta1_g*x_s + 0.5*tau*w_sq_sum)
        
        if (spatial_model=="NB") {
          ## Negative binomial mixed model marginal variance
          nb_marginal_variance[g,s] <- exp(2*(beta0_g + beta1_g*x_s) + tau*w_sq_sum)*(exp(tau*w_sq_sum) - 1)*(r_est[g]+1)/r_est[g] + 
            exp(beta0_g + beta1_g*x_s + 0.5*tau*w_sq_sum) +
            exp(2*(beta0_g + beta1_g*x_s + 0.5*tau*w_sq_sum))/r_est[g]
        }
        
      }
    }
    if (spatial_model=="Pois") {
      Y_normalized_filtered <- (Y_true - exp(B_t_est%*%X_t))/sqrt(marginal_variance)
    }
    if (spatial_model=="NB") {
      Y_normalized_filtered <- (Y_true - exp(B_t_est%*%X_t))/sqrt(nb_marginal_variance)
    }
    
  }
  
  ## Old code for calculating marginal variance for non spatially normalized counts
  ## when we ran the Nimble poisson model with a random effect for extra variance
  #log_mu_hat <- B_t_est%*%X_t
  #mu_hat <- exp(log_mu_hat)
  #for (g in 1:n_genes) {
  #  beta0_g <- B_est[1,g]
  #  beta1_g <- B_est[2,g]
  #  for (s in 1:n_spots) {
  #    x_s <- X_t[2,s]
  #    summand1 <- exp(2*(beta0_g + beta1_g*x_s + sigma_est[g]^2))
  #    summand2 <- -1*exp(2*(beta0_g + beta1_g*x_s) + sigma_est[g]^2)
  #    summand3 <- exp(0.5*sigma_est[g]^2)
  #    marginal_variance[g,s] <- summand1 + summand2 + summand3
  #  }
  #}
  #Y_normalized_filtered <- (Y_true - mu_hat)/sqrt(marginal_variance)
  
  ## Original unfiltered gene spot matrix
  Y_original <- output$prepared_inputs$original_gs
  original_gene_number <- nrow(Y_original)
  normalized_gene_names <- rownames(Y_normalized_filtered)
  all_gene_names <- paste0("G", 1:original_gene_number)
  missing_genes <- which(!(all_gene_names %in% normalized_gene_names))
  missing_gene_names <- all_gene_names[missing_genes]
  
  
  ## Find the original counts for the genes that were filtered out prior to
  ## normalization
  if (length(missing_genes)==1) {
    missing_gene_rows <- t(as.matrix(Y_original[missing_genes,]))
  } else {
    missing_gene_rows <- Y_original[missing_genes,]
  }
  ## If we have selected "zero" as our filter method, replace all these
  ## original counts with 0's (note they will have been mostly 0's to begin with
  ## anyway)
  if (filter_method=="zero") {
    missing_gene_rows[] <- 0
  }
  
  ## If we had genes filtered out, insert those missing counts into our
  ## normalized Y matrix (in the appropriate positions).
  if (length(missing_genes) > 0) {
    rownames(missing_gene_rows) <- missing_gene_names
    ## Normalized Y, adding in the original counts for the genes we filtered and
    ## didn't perform the procedure on because they were too lowly expressed
    Y_normalized <- rbind(Y_normalized_filtered, missing_gene_rows)
    Y_normalized <- Y_normalized[str_sort(rownames(Y_normalized), numeric=T),]
  } else {
    Y_normalized <- Y_normalized_filtered
  }
  
  results <- list(Y_original=Y_original,
                  Y_original_filtered=Y_true,
                  Y_normalized_filtered=Y_normalized_filtered,
                  Y_normalized=Y_normalized,
                  X_t=X_t,
                  B_t_est=B_t_est,
                  Z_est=Z_est,
                  W_est=W_est,
                  r_est=r_est,
                  coords=coords,
                  missing_genes=missing_genes,
                  missing_gene_names=missing_gene_names,
                  normalization_type=normalization_type)
  
  return(results)
}


### Create various plots of the true, normalized, and fitted counts of the
### most variable gene, using the output from the construct_est_params function.
get_mvg_count_plots <- function(est_params) {
  coords_relabeled <- est_params$coords %>% data.frame() %>%
    rename(xcoords=x, ycoords=y) %>% as.matrix()
  
  ## Get most variable gene, based on the true log counts
  v_ranking <- rank_gene_variances(est_params$Y_log)
  mvg_index <- v_ranking$mvg_indices[1]
  true_log_mvg_counts <- t(as.matrix(v_ranking$mvg_count_list[[1]]))
  
  true_mvg_counts <- t(as.matrix(est_params$Y_true[mvg_index,]))
  normalized_mvg_counts <- t(as.matrix(est_params$Y_normalized[mvg_index,]))
  
  fitted_mvg_counts <- t(as.matrix(est_params$Y_fitted[mvg_index,])) 
  residuals_mvg_counts <- t(as.matrix((est_params$Y_true - est_params$Y_fitted)[mvg_index,])) 
  
  
  ## True log counts
  true_log_plots <- plot_spatialpcs(z=true_log_mvg_counts, 
                                    coords=coords_relabeled, labs=F)
  
  ## True counts
  true_plots <- plot_spatialpcs(z=true_mvg_counts, 
                                coords=coords_relabeled, labs=F) 
  
  ## Normalized counts
  normalized_plots <- plot_spatialpcs(z=normalized_mvg_counts, 
                                      coords=coords_relabeled, labs=F) 
  
  ## Fitted counts
  fitted_plots <- plot_spatialpcs(z=fitted_mvg_counts, 
                                  coords=coords_relabeled, labs=F) 
  
  ## Residuals: Y - fitted Y
  residual_plots <- plot_spatialpcs(z=residuals_mvg_counts, 
                                    coords=coords_relabeled, labs=F) 
  
  ## Edit labels and dot size
  for (i in 1:length(true_log_plots)) {
    true_log_plots[[i]] <- true_log_plots[[i]] +
      scale_color_continuous(type="viridis", option="inferno") +
      labs(title="True Log Counts of the Most Variable Gene", color="Counts")
    
    true_plots[[i]] <- true_plots[[i]] +
      scale_color_continuous(type="viridis", option="inferno") +
      labs(title="True Counts of the Most Variable Gene", color="Counts")
    
    normalized_plots[[i]] <- normalized_plots[[i]] +
      scale_color_continuous(type="viridis", option="inferno") +
      labs(title="Normalized Counts of the Most Variable Gene", color="Counts")
    
    fitted_plots[[i]] <- fitted_plots[[i]] +
      scale_color_continuous(type="viridis", option="inferno") +
      labs(title="Fitted Counts of the Most Variable Gene", color="Counts")
    
    residual_plots[[i]] <- residual_plots[[i]] +
      scale_color_continuous(type="viridis", option="inferno") +
      labs(title="Residuals of the Most Variable Gene", color="Counts")
  }
  
  results <- list(
    true_log_plots=true_log_plots,
    true_plots=true_plots,
    normalized_plots=normalized_plots,
    fitted_plots=fitted_plots,
    residual_plots=residual_plots
  )
  
  return(results)
}



#==============================================================================#
#==============================================================================#
#==============================================================================#

### MOVING WINDOW FUNCTION
moving_window <- function(spots, x_subsections, y_subsections, 
                          proportion_overlap=0.25) {
  
  # Range of x and y coordinates
  x_range <- c(min(spots[,"x"]), max(spots[,"x"])) #x
  y_range <- c(min(spots[,"y"]), max(spots[,"y"])) #y
  
  x_section_width <- (x_range[2] - x_range[1])/x_subsections
  y_section_width <- (y_range[2] - y_range[1])/y_subsections
  
  x_buffer <- x_section_width*proportion_overlap
  y_buffer <- y_section_width*proportion_overlap
  
  x_markers <- c(x_range[1])
  y_markers <- c(y_range[1])
  for (i in 1:x_subsections) {
    #if (length(x_markers)==i) {
    ## If the max x coordinate is more than section_length away from the
    ## previous marker, add a new marker
    if (x_range[2] - x_markers[i] >= x_section_width) {
      new_marker <- x_markers[i] + x_section_width
      x_markers <- c(x_markers, new_marker)
      #}
    }
  }
  if (length(x_markers) < (x_subsections+1)) {
    x_markers <- c(x_markers, x_range[2])
  }
  
  for (i in 1:y_subsections) {
    #if (length(y_markers)==i) {
    #remaining_y_range <- round(y_range[2] - y_markers[i], digits=3)
    if (y_range[2] - y_markers[i] >= y_section_width) {
      new_marker <- y_markers[i] + y_section_width
      y_markers <- c(y_markers, new_marker)
    }
    #}
  }
  if (length(y_markers) < (y_subsections+1)) {
    y_markers <- c(y_markers, y_range[2])
  }
  
  # for each section, define x coords as x_markers[i], x_markers[i]+section_width if 
  # x_markers[i]+section_width <= x_range[2] (max x coord)
  sections <- list()
  overlapping_sections <- list()
  for (i in 1:(length(y_markers)-1)) {
    #overlapping_y_lower <- y_markers[i] - y_buffer
    #overlapping_y_upper <- y_markers[i] + y_section_width + y_buffer
    
    # At a fixed y interval, loop over all x intervals
    if (y_markers[i] + y_section_width <= y_range[2]) {
      if (i==1) { # if on the first section
        curr_y_indices <- which( spots[,"y"] >= y_markers[i] & spots[,"y"] <= y_markers[i]+y_section_width )
      } else { # after that, don't include the last point from the previous section 
        ## (hence the greater than instead of greater than or equal to)
        curr_y_indices <- which( spots[,"y"] > y_markers[i] & spots[,"y"] <= y_markers[i]+y_section_width )
      }
    } else { # last section; go up to the max y range
      if (i==1) {
        curr_y_indices <- which(spots[,"y"] >= y_markers[i] & spots[,"y"] <= y_range[2])
      } else {
        curr_y_indices <- which(spots[,"y"] > y_markers[i] & spots[,"y"] <= y_range[2])
      }
      
    }
    
    for (j in 1:( length(x_markers)-1 )) {
      #overlapping_x_lower <- x_markers[j] - x_buffer
      #overlapping_x_upper <- x_markers[j] + x_section_width + x_buffer
      
      if (x_markers[j] + x_section_width <= x_range[2]) {
        if (j==1) {
          curr_x_indices <- which( spots[,"x"] >= x_markers[j] & spots[,"x"] <= x_markers[j]+x_section_width )
        } else {
          curr_x_indices <- which( spots[,"x"] > x_markers[j] & spots[,"x"] <= x_markers[j]+x_section_width )
        }
        
      } else {
        if (j==1) {
          curr_x_indices <- which(spots[,"x"] >= x_markers[j] & spots[,"x"] <= x_range[2])
        } else {
          curr_x_indices <- which(spots[,"x"] > x_markers[j] & spots[,"x"] <= x_range[2])
        }
      }
      
      curr_indices <- intersect(curr_x_indices, curr_y_indices)
      
      # if ( length(sections) > 0 ) {
      #   if (length(sections[[length(sections)]]) < 30 ) { # previous section has 1-99 spots
      #     sections[[length(sections)]] <- c(sections[[length(sections)]], curr_indices)
      #     sections[[length(sections)]] <- unique(sections[[length(sections)]]) # get rid of duplicate indices
      #   } else if (length(curr_indices) >= 30) { # previous section has >= 100 spots & so does current section
      #     sections[[length(sections)+1]] <- curr_indices
      #   } else { # previous section has >= 100 spots & current section has less
      #     sections[[length(sections)]] <- c(sections[[length(sections)]], curr_indices)
      #     sections[[length(sections)]] <- unique(sections[[length(sections)]]) # get rid of duplicate indices
      #   }
      # } else if (length(curr_indices) > 0) { # no sections have been defined yet & current section has > 0 spots
      #   sections[[length(sections)+1]] <- curr_indices
      # }
      if (length(curr_indices) > 0) {
        sections[[length(sections)+1]] <- curr_indices
      }
    }
  }
  
  ## Assign the coordinates corresponding to each section
  unique_coordinates <- vector(mode="list", length=length(sections))
  overlapping_sections <- vector(mode="list", length=length(sections))
  overlapping_coordinates <- vector(mode="list", length=length(sections))
  
  for (i in 1:length(sections)) {
    sections[[i]] <- sort(sections[[i]])
    if (length(sections[[i]])==1) {
      curr_spots <- data.frame(x=spots[sections[[i]],]["x"],
                               y=spots[sections[[i]],]["y"])
    } else {
      curr_spots <- spots[sections[[i]],] %>% data.frame()
    }
    curr_spots$section <- i
    unique_coordinates[[i]] <- curr_spots
    
    curr_x_range <- c(min(curr_spots$x), max(curr_spots$x))
    curr_y_range <- c(min(curr_spots$y), max(curr_spots$y))
    
    ## Define a wider set of coordinates: those within the range of the current 
    ## section plus a buffer defined earlier
    curr_x_overlap <- which(spots[,"x"] >= curr_x_range[1] - x_buffer & 
                              spots[,"x"] <= curr_x_range[2] + x_buffer)
    curr_y_overlap <- which(spots[,"y"] >= curr_y_range[1] - y_buffer & 
                              spots[,"y"] <= curr_y_range[2] + y_buffer)
    curr_indices_overlap <- intersect(curr_x_overlap, curr_y_overlap)
    
    overlapping_sections[[i]] <- sort(curr_indices_overlap)
    curr_spots_overlap <- spots[overlapping_sections[[i]],] %>% data.frame()
    curr_spots_overlap$section <- i
    overlapping_coordinates[[i]] <- curr_spots_overlap
    
  }
  
  midpoints <- matrix(nrow=length(sections), ncol=2, data=NA)
  for (i in 1:length(sections)) {
    
    if (length(sections[[i]])==1) {
      curr_coords <- matrix(data=spots[sections[[i]],], nrow=1, ncol=2)
      colnames(curr_coords) <- colnames(spots)
    } else {
      curr_coords <- spots[sections[[i]],] 
    }
    #curr_coords <- spots[sections[[i]],]
    x_mid <- (max(curr_coords[,"x"]) + min(curr_coords[,"x"]))/2
    y_mid <- (max(curr_coords[,"y"]) + min(curr_coords[,"y"]))/2
    midpoints[i,] <- c(x_mid, y_mid)
  }
  colnames(midpoints) <- c("x","y")
  midpoints_df <- data.frame(midpoints)
  
  
  mw <- list(unique_sections=sections, unique_coordinates=unique_coordinates,
             overlapping_sections=overlapping_sections, 
             overlapping_coordinates=overlapping_coordinates,
             x_markers=x_markers, y_markers=y_markers, midpoints=midpoints_df)
  return(mw)
}










### No longer used - old version of the moving window function
### Had a couple of issues - was sometimes not assigning all spots to a section.
# Function for creating moving window sub-section squares
# Parameters: x, y coordinates of all the spots (as a matrix), # of x subsections, # of y subsections, proportion overlap
# moving_window <- function(spots, x_subsections, y_subsections, proportion_overlap) {
#   
#   # Range of x and y coordinates
#   x_range <- c(min(spots[,1]), max(spots[,1])) #x
#   y_range <- c(min(spots[,2]), max(spots[,2])) #y
#   
#   # We will make our windows 1/(subsections^2) the size of the full map, so 1/subsections the size of the x range and 
#   # 1/subsections the size of the y range. We will overlap proportion_overlap of the previous window with each 
#   # subsequent window
#   overlap_multiplier <- (1/(1 - proportion_overlap))
#   x_increment <- (x_range[2] - x_range[1])/(x_subsections*overlap_multiplier)
#   y_increment <- (y_range[2] - y_range[1])/(y_subsections*overlap_multiplier)
#   
#   x_section_width <- (x_range[2] - x_range[1])/x_subsections
#   y_section_width <- (y_range[2] - y_range[1])/y_subsections
#   
#   x_markers <- c(x_range[1])
#   y_markers <- c(y_range[1])
#   for (i in 1:ceiling(overlap_multiplier*x_subsections)) {
#     if (length(x_markers)==i) {
#       if (x_range[2] - x_markers[i] > x_section_width) {
#         new_marker <- x_markers[i] + x_increment
#         x_markers <- c(x_markers, new_marker)
#       }
#     }
#   }
#   
#   for (i in 1:ceiling(overlap_multiplier*y_subsections)) {
#     if (length(y_markers)==i) {
#       remaining_y_range <- round(y_range[2] - y_markers[i], digits=3)
#       if (remaining_y_range >= y_section_width) {
#         new_marker <- y_markers[i] + y_increment
#         y_markers <- c(y_markers, new_marker)
#       }
#     }
#   }
#   
#   # for each section, define x coords as x_markers[i], x_markers[i]+section_width if 
#   # x_markers[i]+section_width <= x_range[2] (max x coord)
#   sections <- list()
#   for (i in 1:(length(y_markers)-1)) {
#     
#     # At a fixed y interval, loop over all x intervals
#     if (y_markers[i] + y_section_width <= y_range[2]) {
#       curr_y_indices <- which( spots[,2] >= y_markers[i] & spots[,2] < y_markers[i]+y_section_width )
#     } else {
#       curr_y_indices <- which(spots[,2] >= y_markers[i] & spots[,2] < y_range[2])
#     }
#     
#     for (j in 1:( length(x_markers)-1 )) {
#       if (x_markers[j] + x_section_width <= x_range[2]) {
#         curr_x_indices <- which( spots[,1] >= x_markers[j] & spots[,1] < x_markers[j]+x_section_width )
#       } else {
#         curr_x_indices <- which(spots[,1] >= x_markers[j] & spots[,1] < x_range[2])
#       }
#       
#       curr_indices <- intersect(curr_x_indices, curr_y_indices)
#       
#       if ( length(sections) > 0 ) {
#         if (length(sections[[length(sections)]]) < 100 ) { # previous section has 1-99 spots
#           sections[[length(sections)]] <- c(sections[[length(sections)]], curr_indices)
#           sections[[length(sections)]] <- unique(sections[[length(sections)]]) # get rid of duplicate indices
#         } else if (length(curr_indices) >= 100) { # previous section has >= 100 spots & so does current section
#           sections[[length(sections)+1]] <- curr_indices
#         } else { # previous section has >= 100 spots & current section has less
#           sections[[length(sections)]] <- c(sections[[length(sections)]], curr_indices)
#           sections[[length(sections)]] <- unique(sections[[length(sections)]]) # get rid of duplicate indices
#         }
#       } else if (length(curr_indices) > 0) { # no sections have been defined yet & current section has > 0 spots
#         sections[[length(sections)+1]] <- curr_indices
#       }
#       
#     }
#   }
#   
#   ## We also want to make a list of the indices that are unique to each section:
#   repeated_spot_indices <- unlist(sections)[duplicated(unlist(sections))]
#   unique_sections <- lapply(1:length(sections), function(i) {
#     curr_section <- sections[[i]]
#     #keep_indices <- curr_section[which(!(curr_section %in% repeated_spot_indices))]
#     keep_indices <- which(!(curr_section %in% repeated_spot_indices))
#     return(keep_indices)
#   })
#   
#   midpoints <- matrix(nrow=length(sections), ncol=2, data=NA)
#   for (i in 1:length(sections)) {
#     curr_coords <- spots[sections[[i]],]
#     x_mid <- (max(curr_coords[,1]) + min(curr_coords[,1]))/2
#     y_mid <- (max(curr_coords[,2]) + min(curr_coords[,2]))/2
#     midpoints[i,] <- c(x_mid, y_mid)
#   }
#   colnames(midpoints) <- c("x","y")
#   midpoints_df <- data.frame(midpoints)
#   
#   
#   mw <- list(sections=sections, unique_sections=unique_sections,
#              x_markers=x_markers, y_markers=y_markers, midpoints=midpoints_df)
#   return(mw)
# }