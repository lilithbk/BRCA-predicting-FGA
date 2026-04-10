library(tidyverse)

# Load data
df <- read.delim("ProcessedData/rna_expression_with_clinical_vars_of_interest.tsv")

# Make different dataframe without sample ID's, FGA, and the two categorical variables
df.2 <- df %>%
  select(-c(sampleId, FRACTION_GENOME_ALTERED, SUBTYPE, GENETIC_ANCESTRY_LABEL))

# Initialize list in which to define a new observation
new.obs <- list()

# Iterate through each quantitative predictor
quant.vars <- names(df.2)
set.seed(27) # For reproducibility
for (i in 1:length(quant.vars)){
  df.temp <- df.2[, i]
  # Get new observation by sampling a uniform distribution where the min and max are the observed min and max values for that predictor
  new.obs[quant.vars[i]] <- runif(1, min(df.temp), max(df.temp))
}

# Add new observation for categorical variables by sampling from unique observed values
## BRCA subtype
subtypes <- unique(df$SUBTYPE)
new.obs["SUBTYPE"] <- sample(subtypes, 1)
## Genetic ancestry label
ancestry <- unique(df$GENETIC_ANCESTRY_LABEL)
new.obs["GENETIC_ANCESTRY_LABEL"] <- sample(ancestry, 1)

# Load linear models
full.lm <- readRDS("LinearRegressions/full_model.rds")
reduced.lm <- readRDS("LinearRegressions/reduced_model.rds")

# Use both models to get a prediction interval for FGA using the new observation
predict(full.lm, new.obs, interval = "prediction")
  # fit         lwr       upr
  # 0.3602199 -0.06896844 0.7894083
predict(reduced.lm, new.obs, interval = "prediction")
  # fit         lwr       upr
  # 0.2434185 -0.03925945 0.5260964

# Convert new.obs to data frame to make saving easier
new.obs <- data.frame(predictor = names(new.obs), new_observation = unlist(new.obs))
write.csv(new.obs, "LinearRegressions/new_observation.csv", row.names = FALSE)
