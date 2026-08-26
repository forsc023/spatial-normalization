
### PURPOSE:
### In this program, I modify the source code of the CreateSpatialPCAObject function
### from the SpatialPCA package to create a modified object that can still be
### entered into the subsequent SpatialPCA functions.

### The original SpatialPCA function takes in raw count data and automatically
### normalizes it using Seruat's SCTransform procedure. It also calls SPARK to
### detect spatially variable genes, and goes on to estimate the spatial PCs
### with just the spatially variable genes. 

### I want to be able to use the SpatialPCA procedure on non-count data (i.e.,
### gene spot matrices that have non whole-number entries).
### So in this script, I modify the function to skip the SCTransform and SPARK steps.
### The new function takes in a gene-spot matrix, which can be pre-normalized
### or not and does not need to contain whole number values, as well as the
### spatial coordinates corresponding to the columns (spots) of the matrix. 

### This also contains the other functions to run SpatialPCA, copied directly
### from the package documentation, which allows them to be used without installing
### the SpatialPCA package. We used this in our paper to be able to run our procedure
### in parallel in the Minnesota Supercomputing Institute (MSI), because there was
### some trouble with installing the SpatialPCA package on MSI.

## Setting the class:
## In the original SpatialPCA code, this defines a class "SpatialPCA".
## Here, the only thing I changed is to call the class "SpatialPCA_MODIFIED".
setClass("SpatialPCA_MODIFIED", slots=list(
  counts = "ANY",
  normalized_expr = "ANY",
  project = "character",
  covariate = "ANY",
  location = "matrix", 
  kernelmat = "ANY",
  kerneltype = "character",
  bandwidthtype = "character",
  bandwidth = "numeric",
  sparseKernel="logical",
  sparseKernel_tol = "numeric",
  sparseKernel_ncore = "numeric",
  fast = "logical",
  eigenvecnum = "numeric",
  SpatialPCnum = "numeric",
  tau = "numeric",
  sigma2_0 = "numeric",
  W = "ANY",
  SpatialPCs = "ANY",
  highPCs = "ANY",
  highPos = "ANY",
  expr_pred="ANY",
  params = "ANY"
) )

## Modified function:
## This code is taken directly from the CreateSpatialPCAObject() source code.
## The changes I made:
## 1) Remove the SCTransform step
## 2) Remove the SPARK step, which filters to only include spatially variable genes
## 3) Create an object of class SpatialPCA_MODIFIED instead of class SpatialPCA
CreateSpatialPCAObject_MODIFIED <- function(counts, location, covariate=NULL,
                                            project = "SpatialPCA", 
                                            min.loctions = 20,  min.features=20){
  
  ## check dimension
  if(ncol(counts)!=nrow(location)){
    stop("The number of cells in counts and location should be consistent (counts -- m genes x n locations; location -- n locations x d dimension).")
  }# end fi
  
  ## check data order should consistent
  if(!identical(colnames(counts), rownames(location))){
    stop("The column names of counts and row names of location should be should be matched (counts -- m genes x n locations; location -- n locations x d dimension).")
  }# end fi
  
  ## inheriting
  object <- new(
    Class = "SpatialPCA_MODIFIED",
    counts = counts,
    location = location,
    project = project
  )
  
  if(!is.null(covariate)){
    ## check data order should consistent
    if(!identical(rownames(covariate), rownames(location))){
      stop("The row names of covariate and row names of location should be should be matched (covariate -- n locations x q covariates; location -- n locations x d dimension).")
    }# end fi
    
    q=dim(covariate)[2]
    n_covariate=dim(covariate)[1]
    # remove the intercept if added by user, later intercept will add automatically
    if(length(unique(covariate[,1])) == 1){
      covariate = covariate[, -1]
      q=q-1
    }# end fi
    
    object@covariate = as.matrix(covariate,n_covariate,q)
    
  }# end fi
  
  
  object@counts <- counts # store count matrix in sparse matrix
  object@location <- location
  object@project <- project
  
  #### SET normalized_expr 
  object@normalized_expr <- counts
  
  rm(counts) # to save memory
  rm(location)
  
  #  covariates, i.e., confounding or batch effects
  if(!is.null(covariate)){
    object@covariate = object@covariate[match(colnames(object@normalized_expr), rownames(object@location)),1:q]
    object@covariate = as.matrix(object@covariate,dim(object@normalized_expr)[2],q )
  }
  
  ## store count matrix as a sparse matrix
  if(class(object@counts)[1] != "dgCMatrix" ){
    object@counts <- as(object@counts, "dgCMatrix")
  }# end fi
  
  object@params = list()
  
  return(object)
}# end function





#===============================================================================#
#===============================================================================#
#===============================================================================#
### The functions below are unchanged from the original SpatialPCA package. I am
### copying them here so I can use them in MSI without needing to install the
### SpatialPCA package, because the package installation was a difficult process
### on my local device that is likely to be extremely difficult to replicate in 
### the MSI environment.



SpatialPCA_buildKernel = function(object, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL,sparseKernel=FALSE,sparseKernel_tol=1e-20,sparseKernel_ncore=1) {
  
  ## extract the data from the slot of object, createSpatialPCAobject() function goes first
  if(length(object@counts) == 0) {
    stop("object@counts has not been set. Run CreateSpatialPCAObject() first and then retry.")
  }# end fi
  
  object@kerneltype = kerneltype
  object@bandwidthtype = bandwidthtype
  
  cat(paste("## Selected kernel type is: ", object@kerneltype," \n"))
  
  #************************************************************#
  #    Calculate the bandwidth for kernel matrix               #
  #************************************************************#
  
  #cat(paste("## Scale the expression of each gene. \n"))
  expr=object@normalized_expr
  for(i in 1:dim(object@normalized_expr)[1]){
    expr[i,] = scale(object@normalized_expr[i,])
  }
  
  if(is.null(bandwidth.set.by.user)){
    object@bandwidth = bandwidth_select(expr, method=object@bandwidthtype)
  }else{
    object@bandwidth = bandwidth.set.by.user
  }
  object@params$expr=expr
  
  rm(expr)
  
  cat(paste("## The bandwidth is: ", object@bandwidth," \n"))
  
  #************************************************************#
  #    Calculate the kernel matrix with above bandwidth        #
  #************************************************************#
  
  location_normalized = scale(object@location)
  
  if(sparseKernel==FALSE){
    cat(paste("## Calculating kernel matrix\n"))
    object@kernelmat = kernel_build(kerneltype=object@kerneltype,location=location_normalized, bandwidth=object@bandwidth)
    object@sparseKernel=sparseKernel
  }else if(sparseKernel==TRUE){
    cat(paste("## Calculating sparse kernel matrix\n"))
    object@sparseKernel=sparseKernel
    object@sparseKernel_tol = sparseKernel_tol
    object@sparseKernel_ncore = sparseKernel_ncore
    object@kernelmat = kernel_build_sparse(kerneltype=object@kerneltype,location=location_normalized, bandwidth=object@bandwidth,tol = object@sparseKernel_tol, ncores=object@sparseKernel_ncore)
    
  }
  
  cat(paste("## Finished calculating kernel matrix.\n"))
  
  # return results
  return(object)
}# end function



#' @title Select bandwidth in Gaussian kernel.
#' @description This function selects bandwidth in Gaussian kernel.
#' @param expr A m gene by n location matrix of normalized gene expression matrix.
#' @param method The method used in bandwidth selection, "SJ" usually for small sample size data, "Silverman" usually for large sample size data.
#' @return A numeric value of calculated bandwidth.
#' @export
bandwidth_select=function (expr, method)
{
  N = dim(expr)[2]
  if (method == "SJ") {
    
    bw_SJ = c()
    for (i in 1:dim(expr)[1]) {
      tryCatch({
        #print(i)
        bw_SJ[i] = bw.SJ(expr[i, ], method = "dpi")
      }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
    }
    
    beta = median(na.omit(bw_SJ))
  }
  else if (method == "Silverman") {
    bw_Silverman = c()
    for (i in 1:dim(expr)[1]) {
      tryCatch({
        bw_Silverman[i] = bw.nrd0(expr[i, ])
      }, error=function(e){cat("Gene",i," :",conditionMessage(e), "\n")})
    }
    beta = median(na.omit(bw_Silverman))
  }
}



#' @title Build kernel matrix.
#' @description This function calculates kernel matrix from spatial locations.
#' @param kerneltype The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel, and "delaunday" for gaussian kernel built with non-linear Delaunay triangulation based distance.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @return The kernel matrix for spatial relationship between locations.
#' @export
kernel_build = function (kerneltype = "gaussian", location, bandwidth)
{
  if (kerneltype == "gaussian") {
    K = exp(-1*as.matrix(dist(location)^2)/bandwidth)
  }
  else if (kerneltype == "cauchy") {
    K = 1/(1 + 1*as.matrix(dist(location)^2)/as.numeric(bandwidth))
  }
  else if (kerneltype == "quadratic") {
    ED2=1*as.matrix(dist(location)^2)
    K = 1 - ED2/(ED2 + as.numeric(bandwidth))
  }else if (kerneltype == "delaunday") {
    require(spatstat.geom)
    tmp <- ppp(location[,1], location[,2],window=owin(c(-3,3),c(-3,3)))
    Delaunay_dist=delaunayDistance(tmp)
    K = exp(-Delaunay_dist^2/30) 
  }
  return(K)
}




#' @title Build sparse kernel matrix.
#' @description This function calculates kernel matrix.
#' @param kerneltype The type of kernel to be used, either "gaussian", or "cauchy" for cauchy kernel, or "quadratic" for rational quadratic kernel.
#' @param location A n by d matrix of cell/spot location coordinates.
#' @param bandwidth A numeric value of bandwidth.
#' @param tol A numeric value of cut-off value when building sparse kernel matrix.
#' @param ncores A integer value of number of CPU cores to use when building sparse kernel matrix.
#' @return The sparse kernel matrix for spatial relationship between locations.
#'
#' @import parallel
#' @import MASS
#' @import pdist
#' @import tidyr
#'
#' @export
kernel_build_sparse = function(kerneltype,location, bandwidth,tol, ncores)
{
  
  # suppressMessages(require(tidyr))
  # suppressMessages(require(parallel))
  # suppressMessages(require(MASS))
  # suppressMessages(require(pdist))
  # suppressMessages(require(Matrix))
  
  if (kerneltype == "gaussian") {
    fx_gaussian <- function(i){
      line_i = rep(0,dim(location)[1])
      line_i[i] = 1
      line_i[-i] = exp(-(pdist(location[i,],location[-i,])@dist^2)/bandwidth)
      ind_i=which(line_i>=tol)
      return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
    }
    
    results = mclapply(1:dim(location)[1], fx_gaussian, mc.cores = ncores)
    tib = tibble(results)  %>%  unnest_wider(results)
    K_sparse = Matrix::sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
    #K = exp(-1*as.matrix(dist(location)^2)/bandwidth)
  } else if (kerneltype == "cauchy") {
    fx_cauchy <- function(i){
      line_i = rep(0,dim(location)[1])
      line_i[i] = 1
      line_i[-i] = 1/(1 + (pdist(location[i,],location[-i,])@dist^2)/as.numeric(bandwidth))
      ind_i=which(line_i>=tol)
      return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
    }
    
    results = mclapply(1:dim(location)[1], fx_cauchy, mc.cores = ncores)
    tib = tibble(results)  %>%  unnest_wider(results)
    K_sparse = Matrix::sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
    # K = 1/(1 + 1*as.matrix(dist(location)^2)/as.numeric(bandwidth))
  }else if (kerneltype == "quadratic"){
    
    fx_quadratic <- function(i){
      line_i = rep(0,dim(location)[1])
      line_i[i] = 1
      ED2=pdist(location[i,],location[-i,])@dist^2
      line_i[-i] = 1 - ED2/(ED2 + as.numeric(bandwidth))
      ind_i=which(line_i>=tol)
      return(list("ind_i"=ind_i,"ind_j"=rep(i,length(ind_i)),"val_i"=line_i[ind_i] ))
    }
    
    results = mclapply(1:dim(location)[1], fx_quadratic, mc.cores = ncores)
    tib = tibble(results)  %>%  unnest_wider(results)
    K_sparse = sparseMatrix(i =unlist(tib[[1]]), j= unlist(tib[[2]]), x= unlist(tib[[3]]),  dims = c(dim(location)[1],dim(location)[1] ))
    # ED2=1*as.matrix(dist(location)^2)
    # K = 1 - ED2/(ED2 + as.numeric(bandwidth))
    
  }else if (kerneltype == "delaunday"){
    
    require(spatstat.geom)
    tmp <- ppp(location[,1], location[,2],window=owin(c(-3,3),c(-3,3))) # scaled distance often ranges from -3 to 3
    Delaunay_dist=delaunayDistance(tmp)
    K = exp(-Delaunay_dist^2/30) 
    K[K<tol]=0
    K_sparse = as(K, "sparseMatrix")
  }
  
  return(K_sparse)
  
}



SpatialPCA_EstimateLoading = function(object, maxiter=300,initial_tau=1,fast=FALSE,eigenvecnum=NULL,SpatialPCnum=20){
  
  suppressMessages(require(RSpectra))
  set.seed(1234)
  param_ini=log(initial_tau)
  object@SpatialPCnum = SpatialPCnum
  object@fast = fast
  object@params$X = scale(object@location)
  object@params$n = dim(object@params$X)[1]
  object@params$p=dim(object@params$X)[2]
  
  if(is.null(object@covariate)){
    object@params$H = matrix(1, dim(object@params$X)[1],1)
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH = object@params$H%*%HH_inv%*%t(object@params$H)
    object@params$M=diag(object@params$n)-HH
    # Y=expr
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
    object@params$q=1
  }else{
    object@params$q = dim(object@covariate)[2]+1
    object@params$H = matrix(0, object@params$n,object@params$q)
    object@params$H[,1]=1
    object@params$H[,2:object@params$q] = object@covariate
    HH_inv=solve(t(object@params$H)%*%object@params$H,tol = 1e-40)
    HH=object@params$H%*%HH_inv%*%t(object@params$H)
    object@params$M=diag(object@params$n)-HH
    #Y=expr
    object@params$tr_YMY=sum(diag(object@params$expr%*%object@params$M%*%t(object@params$expr)))
    object@params$YM = object@params$expr%*%object@params$M
  }
  
  
  if(fast==FALSE){
    object@fast=fast
    print("Eigen decomposition on kernel matrix!")
    eigen_res = eigen(object@kernelmat)
    object@params$delta = eigen_res$values
    object@params$U = eigen_res$vectors
    print("Using all eigenvectors and eigenvalues in the Kernel matrix!")
  }else{
    object@fast=fast
    if(!is.null(eigenvecnum)){
      print("Eigen decomposition on kernel matrix!")
      object@eigenvecnum=eigenvecnum
      if(object@sparseKernel==TRUE){
        eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum)
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }else{
        eigen_res = eigs_sym(object@kernelmat, k=object@eigenvecnum, which = "LM")
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }
      
      print("Low rank approximation!")
      print(paste0("Using user selected top ",object@eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
    }else if(object@params$n>5000){
      print("Eigen decomposition on kernel matrix!")
      if(object@sparseKernel==TRUE){
        eigen_res = eigs_sym(object@kernelmat, k=20)
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }else{
        eigen_res = eigs_sym(object@kernelmat, k=20, which = "LM")
        object@params$delta = eigen_res$values
        object@params$U = eigen_res$vectors
      }
      print("Low rank approximation!")
      print("Large sample, using top 20 eigenvectors and eigenvalues in the Kernel matrix!")
    }else{
      eigen_res = eigen(object@kernelmat)
      delta_all = eigen_res$values
      U_all = eigen_res$vectors
      ind = which(cumsum(delta_all/length(delta_all))>0.9)[1]
      print("Low rank approximation!")
      print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
      object@params$delta = delta_all[1:ind]
      object@params$U = U_all[,1:ind]
      rm(U_all)
    }
  }
  
  
  object@params$MYt = object@params$M %*% t(object@params$expr)
  object@params$YMMYt = object@params$YM %*% object@params$MYt
  object@params$YMU = object@params$YM %*% object@params$U
  object@params$Xt = t(object@params$H)
  object@params$XtU = object@params$Xt %*% object@params$U
  object@params$Ut = t(object@params$U)
  object@params$UtX = object@params$Ut %*% object@params$H
  object@params$YMX = object@params$YM %*% object@params$H
  object@params$UtU = object@params$Ut %*% object@params$U
  object@params$XtX = object@params$Xt %*% object@params$H
  object@params$SpatialPCnum = SpatialPCnum
  
  
  optim_result =try(optim(param_ini, SpatialPCA_estimate_parameter,params=object@params,control = list(maxit = maxiter), lower = -10, upper = 10,method="Brent"),silent=T)
  
  object@tau = exp(optim_result$par)
  k = dim(object@params$expr)[1]
  n = dim(object@params$expr)[2]
  q=object@params$q
  tauD_UtU_inv = solve(object@tau*diag(object@params$delta) + object@params$UtU, tol = 1e-40)
  YMU_tauD_UtU_inv_Ut = object@params$YMU %*% tauD_UtU_inv %*% object@params$Ut
  YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% object@params$H
  XtU_inv_UtX = object@params$XtU %*% tauD_UtU_inv %*% object@params$UtX
  left = object@params$YMX - YMU_tauD_UtU_inv_UtX
  right = t(left)
  middle = solve(-XtU_inv_UtX, tol = 1e-40)
  G_each = object@params$YMMYt - YMU_tauD_UtU_inv_Ut %*% object@params$MYt - left %*% middle %*% right
  object@W = eigs_sym(G_each, k=SpatialPCnum, which = "LM")$vectors
  object@sigma2_0 = as.numeric((object@params$tr_YMY+F_funct_sameG(object@W,G_each))/(k*(n-q)))
  
  rm(eigen_res)
  rm(tauD_UtU_inv)
  rm(YMU_tauD_UtU_inv_Ut)
  rm(YMU_tauD_UtU_inv_UtX)
  rm(XtU_inv_UtX)
  rm(left)
  rm(right)
  rm(middle)
  rm(G_each)
  gc()
  
  return(object)
}


#' @import RSpectra
SpatialPCA_estimate_parameter = function(param_ini, params){
  # suppressMessages(require(RSpectra))
  set.seed(1234)
  tau=exp(param_ini[1])
  k = dim(params$expr)[1]
  n = dim(params$expr)[2]
  q=params$q
  PCnum=params$SpatialPCnum
  tauD_UtU_inv = solve(tau*diag(params$delta) + params$UtU, tol = 1e-40)
  YMU_tauD_UtU_inv_Ut = params$YMU %*% tauD_UtU_inv %*% params$Ut
  YMU_tauD_UtU_inv_UtX = YMU_tauD_UtU_inv_Ut %*% params$H
  XtU_inv_UtX = params$XtU %*% tauD_UtU_inv %*% params$UtX
  left = params$YMX - YMU_tauD_UtU_inv_UtX
  right = t(left)
  middle = solve(-XtU_inv_UtX, tol = 1e-40)
  G_each = params$YMMYt - YMU_tauD_UtU_inv_Ut %*% params$MYt - left %*% middle %*% right
  log_det_tauK_I = determinant(1/tau*diag(1/params$delta)+ params$UtU, logarithm=TRUE)$modulus[1] + determinant(tau*diag(params$delta), logarithm=TRUE)$modulus[1]
  Xt_invmiddle_X = params$XtX - params$XtU %*% solve(params$UtU + 1/tau *diag( 1/params$delta) , tol = 1e-40) %*% params$UtX
  log_det_Xt_inv_X = determinant(Xt_invmiddle_X, logarithm=TRUE)$modulus[1]
  sum_det=0
  sum_det=sum_det+(0.5*log_det_tauK_I+0.5*log_det_Xt_inv_X  )*PCnum
  
  rm(tauD_UtU_inv)
  rm(YMU_tauD_UtU_inv_Ut)
  rm(YMU_tauD_UtU_inv_UtX)
  rm(XtU_inv_UtX)
  rm(left)
  rm(middle)
  rm(right)
  rm(Xt_invmiddle_X)
  gc()
  
  W_est_here = eigs_sym(G_each, k=PCnum, which = "LM")$vectors
  -(-sum_det -(k*(n-q))/2*log(params$tr_YMY+F_funct_sameG(W_est_here,G_each)))
}


F_funct_sameG = function(X,G){ # G is a matrix
  return_val=0
  for(i in 1: dim(X)[2]){
    return_val=return_val+t(X[,i])%*%G%*%X[,i]
  }
  -return_val
}





SpatialPCA_SpatialPCs= function(object,fast=FALSE,eigenvecnum=NULL){
  
  # suppressMessages(require(RSpectra))
  
  n = object@params$n
  PCnum = object@SpatialPCnum
  Z_hat = matrix(0, PCnum, n)
  tau = object@tau
  W_hat = object@W
  
  if(fast==FALSE){
    U=object@params$U
    delta=object@params$delta
  }else if(fast==TRUE){
    
    if(!is.null(eigenvecnum)){
      print(paste0("Low rank approximation!"))
      print(paste0("Using user selected top ",eigenvecnum," eigenvectors and eigenvalues in the Kernel matrix!"))
      EIGEN = eigs_sym(object@kernelmat, k=eigenvecnum, which = "LM")
      U=EIGEN$vectors
      delta=EIGEN$values
      
    }else if(n>5000){
      fast_eigen_num = ceiling(n*0.1)
      print(paste0("Low rank approximation!"))
      print("Large sample, using top 10% sample size of eigenvectors and eigenvalues in the Kernel matrix!")
      EIGEN = eigs_sym(object@kernelmat, k=fast_eigen_num, which = "LM")
      U=EIGEN$vectors
      delta=EIGEN$values
    }else{
      U=object@params$U
      delta=object@params$delta
      ind=length(delta)
      #print(paste0("Low rank approximation!"))
      print(paste0("Small sample, using top ",ind," eigenvectors and eigenvalues in the Kernel matrix!"))
    }
  }
  object@params$U=U
  object@params$delta=delta
  
  W_hat_t = t(W_hat)
  WtYM = W_hat_t%*% object@params$YM
  WtYMK = WtYM %*% object@kernelmat
  WtYMU = WtYM %*% object@params$U
  Ut=t(object@params$U)
  UtM = Ut %*% object@params$M
  UtMK = UtM %*% object@kernelmat
  UtMU = UtM %*% object@params$U
  middle_inv = solve(1/tau * diag(1/delta) + UtMU, tol = 1e-40)
  
  object@SpatialPCs = tau*WtYMK - tau*WtYMU %*% middle_inv %*% UtMK
  object@SpatialPCs = as.matrix(object@SpatialPCs)
  rm(W_hat_t)
  rm(WtYM)
  rm(WtYMK)
  rm(WtYMU)
  rm(Ut)
  rm(UtM)
  rm(UtMK)
  rm(UtMU)
  rm(middle_inv)
  gc()
  
  
  return(object)
}


