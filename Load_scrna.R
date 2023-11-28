#Script for reading different single cell matrices in different formats 
#From there we will be converting them to seurat objects 

setwd("/Users/michael/Documents/Bioinformatics/Personal_Projects/scRNA/Data/")

library(Seurat)
#library(SeuratDisk)

#For the .RDS format 

rds_obj <-readRDS("name_of_file.rds")
#str(rds_obj) to view the object 

#Loading the 10X CellRanger .HDF5 format
hdf5_obj <- Read10X_h5(filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)

seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj) #converting to a seurat object
str(seurat_hdf5) #Visualizing the seurat objects

#Loading the .mtx file
mtx_obj <- ReadMtx(mtx = "Path_to_file_matrix.mtx",
                   features = "Path_features.tsv",
                   cells = "Path_barcodes.tsv")


#loading the .loom file
loom_obj <- Connect(filename = "path_to_file", mode = "r")
seurat_loom <- as.Seurat(loom_obj)


#.h5ad format
Convert("adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE) #This converts the AnnData object to an h5seurat file
seurat_anndata <- LoadH5Seurat("adata_SS2_for_download.h5seurat")
