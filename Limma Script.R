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

# Filter list of DE genes to 50 of largest log2FC w adj. P-value < 0.05
top50_genes <- coolfilteredgenes %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC))) %>%
  slice(1:50)

# Dowload pakages 
install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(dplyr)

# Scale by row, so expression is centered per gene
log2CPM_top50_scaled <- t(scale(t(top50_genes)))

# Create heatmap
heatmap(log2CPM_top50_scaled,
        scale = "none",  
        col = heat.colors(256),   # change to your preferred color palette
        main = "Top 50 Differentially Expressed Genes",
        Rowv = TRUE,              # cluster genes (rows)
        Colv = TRUE)              # cluster samples (columns)

