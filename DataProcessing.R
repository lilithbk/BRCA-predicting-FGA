# Data source: https://www.cbioportal.org/study/clinicalData?id=brca_tcga_pan_can_atlas_2018 
# (both the clinical and genomic data)

# Load libraries
library(tidyverse)
library(glmnet)
library(faraway)

# First, need to get list of sample ID's that are present in the clinical data
# (Due to inconsistencies within cBioPortal, all samples in clin data are in genomic data but opposite not true)
clin <- read.delim("../Data/brca_tcga_pan_can_atlas_2018_clinical_data.tsv")
clin <- clin %>%
  mutate(sampleId = str_replace_all(sampleId, "-", ".")) # Format to match RNAseq data
samp_ids <- pull(clin, sampleId)

# Filter RNAseq data to contain only the samples in clin data and keep only one gene name column
rna <- read.delim("../Data/brca_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
rna <- rna %>%
  select(-Entrez_Gene_Id) %>% # Hugo symbol much easier to work with
  drop_na(Hugo_Symbol) %>%
  select(c(Hugo_Symbol, any_of(samp_ids))) %>% # Subset for relevant samples
  group_by(Hugo_Symbol) %>%
  filter(n() == 1) %>% # Removing instances of duplicated Hugo symbols
  ungroup()

# Transpose RNA so that Hugo symbol is the columns and samples are rows
rna_mtx <- t(as.matrix(rna[, -1]))
colnames(rna_mtx) <- rna$Hugo_Symbol
rna <- as.data.frame(rna_mtx) # Now the rownames are the sample ID's
rna <- rna %>%
  mutate(sampleId = rownames(rna), .before = UBE2Q2P2) %>% # Make first column sample ID's
  remove_rownames()

# Merge RNAseq data with FGA from clin data
fga <- clin %>%
  select(c(sampleId, FRACTION_GENOME_ALTERED)) %>%
  drop_na()
rna <- merge(fga, rna, by = "sampleId")

# Lasso regression of RNAseq data using FGA as the response
## Prepare data
x <- rna %>% 
  select(-c(FRACTION_GENOME_ALTERED, sampleId)) %>%
  as.matrix()
y <- rna$FRACTION_GENOME_ALTERED
## Do 10-fold CV lasso regression
set.seed(27) # For reproducibility
rna.lasso <- cv.glmnet(x, y, alpha = 1)
## Get Hugo symbols that are non-zero
hugo.lasso <- names(coef(rna.lasso)[, 1][coef(rna.lasso)[, 1] != 0][-1])
## Remove mRNAs from RNAseq data with coefficient of 0
rna <- rna %>%
  select(c(sampleId, FRACTION_GENOME_ALTERED, hugo.lasso))
## Save shrunk RNAseq data to file
write_tsv(rna, "ProcessedData/shrunk_lasso_reg_subset_rna_expression.tsv")
dim(rna) # Down to only 51 genes from 20,531

# Check RNAseq data for multicollinearity issues
rna.cor <- cor(rna[, -(1:2)])
rna.cor[abs(rna.cor) > 0.8 & rna.cor < 1] # No correlations exceed a magnitude of 0.8
rna.vif <- vif(rna[, -(1:2)])
rna.vif[rna.vif > 6] # No VIFs exceed 6, so no major collinearity problems
rna.vif[rna.vif > 2] # Several exceed 2, though, so there are some collinearity problems, which is to be expected (and is unavoidable with this type of data)
mean(rna.vif) # The mean of VIFs does exceed 1, but is not concerningly large

# Creating a version of the rna data with predictors from clin data that we may use
clin.2 <- clin %>%
  select(c(sampleId, AGE, GENETIC_ANCESTRY_LABEL, SUBTYPE))
clin.rna <- merge(clin.2, rna, by = "sampleId")
clin.rna <- clin.rna %>%
  relocate(sampleId, FRACTION_GENOME_ALTERED)
clin.rna <- clin.rna %>%
  filter(GENETIC_ANCESTRY_LABEL != "") %>%
  filter(SUBTYPE != "")
write_tsv(clin.rna, "ProcessedData/rna_expression_with_clinical_vars_of_interest.tsv")
