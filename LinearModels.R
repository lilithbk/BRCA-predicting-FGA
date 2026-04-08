library(tidyverse)
library(rms)

# Load data
df <- read.delim("ProcessedData/rna_expression_with_clinical_vars_of_interest.tsv")

# Convert categorical variables to factors
df <- df %>%
  mutate(SUBTYPE = as.factor(SUBTYPE)) %>%
  mutate(GENETIC_ANCESTRY_LABEL = as.factor(GENETIC_ANCESTRY_LABEL)) %>%
  select(-sampleId) # Unnecessary column

# Full linear model
full.lm <- lm(FRACTION_GENOME_ALTERED ~ ., df)
summary(full.lm)
write_rds(full.lm, "LinearRegressions/full_model.rds")
write_rds(summary(full.lm), "LinearRegressions/full_model_summary.rds")

# Backward stepwise selection using BIC as evaluation metric
reduced.lm <- step(full.lm, direction = "backward", k = log(nrow(df)))
summary(reduced.lm)
write_rds(reduced.lm, "LinearRegressions/reduced_model.rds")
write_rds(summary(reduced.lm), "LinearRegressions/reduced_model_summary.rds")
