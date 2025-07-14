library(limma)
library(tidyverse)

########################################### Put path to input file here
inputFile <- "Custom/GSE200431_all_dep_ctrl_noNA.csv"
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
whichSamples <- whichSamples[]
###########################################

# Extract treatment names from group names
group <- str_remove((str_extract(whichSamples, "[_][A-Za-z]+")), "_") # Ask ChatGPT what this line does
group <- as.factor(group)                                             # Ask ChatGPT what this line does

########################################### # Create design matrix for linear model 
design <- model.matrix()
########################################### 
colnames(design) <- levels(group)
design

# Fit linear model to log2CPM data
fit <- lmFit(log2CPM_matrix[, whichSamples], design)

########################################### Define contrast(s), e.g., treated vs control
contrast.matrix <- makeContrasts()
###########################################

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
topTable(fit2, coef="TreatedVsControl", lfc = 0, number = 20)

