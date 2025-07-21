--------------

library(Seurat)
library(ggplot2)
library(Seurat)
library(ggplot2)
library(dplyr)

# Run commands to reduce diferentiation
seuratObj <- RunPCA(seuratObj)
seuratObj <- NormalizeData(seuratObj)
seuratObj <- FindVariableFeatures(seuratObj)
seuratObj <- ScaleData(seuratObj)
seuratObj <- RunUMAP(seuratObj, dims = 1:10)

# Both UMAP & PCA simplify complex drm(ata by reducing the number of variables while retaining important information

# Create UMAP plot side-by-side with annotations for stromalClass and celltype
DimPlot(seuratObj, reduction = "umap", group.by = "stromalClass") + 
  ggtitle("UMAP by Stromal Class") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

DimPlot(seuratObj, reduction = "umap", group.by = "celltype") + 
  ggtitle("UMAP by Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

# Now, we use FeaturePlot in order to create UMAP projections showing the expressions of 3 genes:

FeaturePlot(seuratObj, features = "Ptn", reduction = "umap")
FeaturePlot(seuratObj,features = "Pecam1", reduction = "umap")
FeaturePlot(seuratObj, features = "Acta2", reduction = "umap") 

# Create a Vln & DotPlot to show distributions of gene expression levels across cell types (Ptn,Acta2 & Pecam1)
VlnPlot(seuratObj, features = c("Ptn", "Pecam1", "Acta2"), group.by = "celltype", pt.size =0.3)
DotPlot(seuratObj, features = c("Ptn", "Pecam1", "Acta2"), group.by = "celltype")

# Identify marker genes

markers <- FindAllMarkers(seuratObj)

# Filter table to top 5 gene markers for each cell type
top_markers <- markers %>%
  group_by(cluster) %>%              
  top_n(5, wt = avg_log2FC) %>%                    
  arrange(cluster, desc(avg_log2FC))                

# Rearrange table 
top_markers <- top_markers %>%
  select(c(6, 7, everything()))  # Move columns 6 and 7 to the front

# Convert cell type to a factor & drop unused levels
seuratObj$celltype <- factor(seuratObj$celltype)
seuratObj$celltype <- droplevels(seuratObj$celltype)

# Create DotPlot showing expression of top cell type marker genes
DotPlot(seuratObj, features = top_markers$gene, group.by = "celltype") + RotatedAxis()
?ggsave
ggsave(filename = "expression.png",
       plot = last_plot(), 
       width = 20,
       height = 20)
  