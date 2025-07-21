library(limma)
library(tidyverse)

########################################### Put path to input file here
inputFile <- "GSE200431_all_dep_ctrl_noNA.csv"
###########################################

# Import log2CPM data into R
log2CPM_matrix <- read.csv(inputFile)  

# Create data.matrix from log2CPM_matrix, with gene symbols as rownames
log2CPM_matrix <-                     
  log2CPM_matrix %>%
  distinct(symbol, .keep_all = TRUE) %>%
  column_to_rownames("symbol") %>%
  dplyr::select(-X) %>%
  as.matrix()

# Create a vector of sample names
whichSamples <- colnames(log2CPM_matrix)  

########################################### Use brackets to subset samples to analyze
whichSamples <- whichSamples[-c(1:12, 19:24)]  

###########################################

# Extract treatment names from group names
group <- str_remove((str_extract(whichSamples, "[_][A-Za-z]+")), "_")
group <- as.factor(group)                                          

########################################### # Create design matrix for linear model 
design <- model.matrix(~ 0 + group)
########################################### 
colnames(design) <- levels(group)
design

# Fit linear model to log2CPM data
fit <- lmFit(log2CPM_matrix[, whichSamples], design)

########################################### Define contrast(s), e.g., treated vs control
contrast.matrix <- makeContrasts(CntrlvsSH= plko - shPTN,levels=design)
###########################################

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
topTable(fit2, coef="CntrlvsSH", lfc = 0, number = 200)

# Create new file
cooltable <- topTable(fit2, coef="CntrlvsSH", lfc = 0, number = 200)

# Filter tale to select genes with adjusted p-value < 0.5 & log2 fold change > 1 or < -1
library(dplyr)
coolfilteredgenes <- cooltable %>%
  filter(adj.P.Val < 0.5 & (logFC > 1 | logFC < -1))

# Write table to CSV file
write.csv(coolfilteredgenes, "GSEFiltered_PTN_Signature.csv", row.names = FALSE)

## Fltered PTN Genes Table
DT::datatable(coolfilteredgenes, options = list(pageLength = 10))

# Dowload pakages 
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(dplyr)

# Remove columns 1-12 and 19-25
log2CPM_matrix_filtered <- log2CPM_matrix[, -c(1:12, 19:25)]

# Filter top 50 genes (based on absolute logFC and adj.P.Val < 0.05)
top50_genes_filtered <- coolfilteredgenes %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1:50)

# Subset the log2CPM matrix for the top 50 genes
log2CPM_top50 <- log2CPM_matrix_filtered[rownames(log2CPM_matrix_filtered) %in% rownames(top50_genes_filtered), ]

# Mean centering the rows of the expression matrix(Log2CPM values)
log2CPM_top50_centered <- t(scale(t(log2CPM_top50), center = TRUE, scale = FALSE))

# Create Heatmap & Legend
Heatmap(log2CPM_top50_centered,
        col = colorRampPalette(c("blue", "white", "red"))(256),  # Custom color palette (blue-white-red)
        show_column_names = TRUE,  # Show column names (samples)
        show_row_names = TRUE,
        width = unit(14, "cm"),  
        height = unit(14, "cm"),
        name = "Top 50 Differentially\nExpressed Genes")

# Separate the posotive‐FC and negative‐FC gene lists
pos_genes <- rownames(top50_genes_filtered)[ top50_genes_filtered$logFC > 0 ]
neg_genes <- rownames(top50_genes_filtered)[ top50_genes_filtered$logFC < 0 ]
      
mat_pos <- log2CPM_top50_centered[ pos_genes, ]
mat_neg <- log2CPM_top50_centered[ neg_genes, ]

library(circlize)

# Create Heatmap for - and + values of the log2FC 
Heatmap(mat_neg,
        col = colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
        name = "Top 50\n(log2FC < 0)",
        width = unit(14, "cm"),  
        height = unit(14, "cm"))
      
Heatmap(mat_pos,
        col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
        width = unit(14, "cm"),  
        height = unit(14, "cm"),
        name = "Top 50\n(log2FC > 0)")

# Insert Libraries

library(Seurat)
library(ggplot2)

# Create DotPlot for GSE postive genes celltypes
DotPlot(seuratObj, features   = posvalid_genes, group.by   = "celltype") +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x = element_blank()) + ggtitle("Expression of Positive‐FC Genes by Cell Type")

# Create DotPlot for GSE negative genes celltypes
DotPlot(seuratObj, features   = neg_genes, group.by   = "celltype") +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x = element_blank()) + ggtitle("Expression of Negative‐FC Genes by Cell Type")
  
# Create DotPlot for GSE positive genes stromal content
DotPlot(seuratObj,features   = pos_genes, group.by   = "stromalClass") +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
axis.title.x = element_blank()) + ggtitle("Expression of Positive‐FC Genes by Stromal Content")



# Create UMAP for +FC stromal content, & label each UMAP plot with the gene name
FeaturePlot(seuratObj, features  = pos_genes, reduction = "umap", 
ncol      = 3, pt.size   = 0.3) &
ggtitle("Expression of +FC genes")
plot_list <- lapply(posvalid_genes, function(gene) 
{FeaturePlot(seuratObj, features = gene, reduction = "umap", pt.size = 0.3) + ggtitle(gene)})

# Arrange all the plots in a grid
gridExtra::grid.arrange(grobs = plot_list, ncol = 3)  # Adjust ncol as needed



# Create 5 groups for -FC stromal content
gp1 <- neg_genes[1:9]
gp2 <- neg_genes[10:18]
gp3 <- neg_genes[19:27]
gp4 <- neg_genes[28:36]
gp5 <- neg_genes[38:41]

#gp1-5 umap
FeaturePlot(seuratObj, features  = gp5, reduction = "umap", ncol = 3, pt.size   = 0.3) &
ggtitle("Expression of +FC genes")
plot_list <- lapply(gp5, function(gene) 
{FeaturePlot(seuratObj, features = gene, reduction = "umap", pt.size = 0.3) + ggtitle(gene)})

# Arrange all the plots in a grid
gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
