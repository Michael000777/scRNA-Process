#Analysis of Single cell data from 10x genomics 

#Loading needed libraries 
library(Seurat)
library(tidyverse)

#Loading the NSCLC dataset 
#ensure you're in the right working directory 

nsclc.sparse.m <- Read10X_h5(filename = "20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5")
str(nsclc.sparse.m)#to see the modalities present 
counts <- nsclc.sparse.m$`Gene Expression`

counts[1:10, 1:10]#to visualize the first 10 rows and 10 coulmns

nsclc.seurat.obj <- CreateSeuratObject(counts = counts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj


#QC
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm")

#Filtering out low quality scores 
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                             percent.mt < 5)
#Normalization step 
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)

#Identify highly variable features
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

#To visualize the top 20 variable features 
top20 <- head(VariableFeatures(nsclc.seurat.obj), 20)

#Plot and visualize the variable features
plot_a <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot_a, points = top20, repel = TRUE)

#Scaling 
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)


#PCA- dimension reduction 
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

#Visualization of PCA results 
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

#Determining dimensionality of the data: number of PC to pick out 
ElbowPlot(nsclc.seurat.obj)

#Clustering
nsclc.seurat.obj <-FindNeighbors(nsclc.seurat.obj, dims = 1:15) #15 comes from the Elbow plot as that represents the number of pca components that explain maximum variability in the dataset

#Understanding res
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1)) #You can further breakdown the resolutions
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.3", label = TRUE) #This resolution works 
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)
DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.7", label = TRUE) #Looks like the higher the resolution the more the clusters are broken down into smaller groups. So the higher you go smilar cell groups start to be broken into sub groups and that might not be biologically relevant to what we want.


#setting identity of clusters
Idents(nsclc.seurat.obj)#This gives you the default identity: which is number of clusters
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.3"


#Non-linear dimensionality reduction UMAP
reticulate::py_install(packages = 'umap-learn') #if needed
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

DimPlot(nsclc.seurat.obj, reduction = "umap")

