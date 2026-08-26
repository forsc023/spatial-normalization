### Dataset: 10X Genomics - Visium CytAssist Gene Expression Libraries of Post-Xenium Human Colon Cancer (FFPE)
### Post-Xenium Replicate 1
### Analyzed with Space Ranger 2.1.0 (spaceranger count pipeline)
### Visium
### 6518 spots
### 18085 genes
### url: https://www.10xgenomics.com/datasets/visium-cytassist-gene-expression-libraries-of-post-xenium-human-colon-cancer-ffpe-using-the-human-whole-transcriptome-probe-set-2-standard

library(Seurat)
library(here)
library(tidyverse)

## Read in data as Seurat object and extract coordinates and the raw gene-spot matrix
data <- Load10X_Spatial(
  data.dir=here("Raw Data"),
  filename="CytAssist_FFPE_Human_Colon_Post_Xenium_Rep1_filtered_feature_bc_matrix.h5",
  assay="Spatial",
  slice="tissue_lowres_image", #this was 'tissue_lowres_image.png' from the Spatial Imaging Data folder
  filter.matrix=TRUE,
)
coords <- data@images$tissue_lowres_image@boundaries$centroids@coords
rawcounts <- data@assays$Spatial@layers$counts
#spot_counts <- colSums(rawcounts)
#which(spot_counts==0)

## Get the spot names
metadata <- data@meta.data
colnames(rawcounts) <- rownames(metadata)

## Read in feature names
features <- read_tsv(here("Raw Data", "filtered_feature_bc_matrix/features.tsv.gz"), col_names=F) %>%
  rename(gene_long=X1, gene=X2) %>%
  mutate(feature_index=row_number())
rownames(rawcounts) <- features$gene_long


### Run SCTransform, for future comparison
### This step is not necessary if you just want to run the spatial normalization procedure.
sct_data <- SCTransform(object=rawcounts, cell.attr=metadata, variable.features.n=nrow(rawcounts))

### Run NormalizeData, which scales the counts and does a log transformation
nd_data <- NormalizeData(object=rawcounts)

### Moving forward, only use the genes SCTransform saved so that we can directly
### compare everything:
sctcounts <- sct_data$y
sct_corrected_counts <- sct_data$umi_corrected
rawcounts_use <- rawcounts[sct_data$variable_features,]
ndcounts_use <- nd_data[sct_data$variable_features,]

## Reorder feature names to match the order outputted from sctransform
sct_features <- data.frame(sct_index=c(1:length(sct_data$variable_features)), 
                           gene_long=sct_data$variable_features)
features_use <- features %>%
  left_join(sct_features, by="gene_long") %>% arrange(sct_index) %>% 
  filter(!is.na(sct_index))
identical(features_use$gene_long, rownames(sctcounts))

## Check that we've put the genes & spots in the same order in all three datasets
identical(rownames(rawcounts_use), rownames(ndcounts_use))
identical(rownames(rawcounts_use), rownames(sctcounts))
identical(colnames(rawcounts_use), colnames(ndcounts_use))
identical(colnames(rawcounts_use), colnames(sctcounts))

spot_names <- data@images$tissue_lowres_image@boundaries$centroids@cells
identical(spot_names, colnames(rawcounts_use))

tissue_positions <- read.csv(here("Raw Data", "spatial/tissue_positions.csv"), header=T) %>%
  filter(in_tissue==1) %>%
  rename(spot=barcode, y_array=array_row, x_array=array_col, y_pixel=pxl_row_in_fullres, 
         x_pixel=pxl_col_in_fullres)

coord_data <- cbind(coords, spot=spot_names) %>% data.frame() %>%
  left_join(tissue_positions, by="spot") %>%
  mutate(x=as.numeric(x), y=as.numeric(y), x_pixel=as.numeric(x_pixel), y_pixel=as.numeric(y_pixel),
         x_array=as.numeric(x_array), y_array=-1*as.numeric(y_array))

identical(coord_data$spot, colnames(rawcounts_use))
coords_matrix <- coord_data %>% dplyr::select(x, y) %>% as.matrix()
rownames(coords_matrix) <- coord_data$spot

write_rds(rawcounts_use, file=here("Processed Data", "rawcounts.rds"))
write_rds(ndcounts_use, file=here("Processed Data", "ndcounts.rds"))
write_rds(sctcounts, file=here("Processed Data", "sctcounts.rds"))
write_rds(coord_data, file=here("Processed Data", "coords_all_info.rds"))
write_rds(coords_matrix, file=here("Processed Data", "coords.rds"))
write_rds(features_use, file=here("Processed Data", "features.rds"))

## Save ND counts as matrix instead of the compressed matrix format. We'll use this
## for the constrained k means clustering because Python can't recognize the compressed
## format.
write_rds(as.matrix(ndcounts_use), file=here("Processed Data", "ndcounts_as_matrix.rds"))

## Save which genes (ie gene indices) SCTransform kept, so we know which ones we
## filtered out of the raw and NormalizeData counts
write_rds(sct_data$umi_corrected, file=here("Processed Data", "sct_kept_genes.rds"))

