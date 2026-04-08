# Read full linear model from file
full.lm <- read_rds("LinearRegressions/full_model.rds")
# Read data from file
df <- read.delim("ProcessedData/rna_expression_with_clinical_vars_of_interest.tsv")

# Residual plot
full.resid <- resid(full.lm)
full.fitted <- fitted(full.lm)

df.2 <- df %>%
  mutate(resid = full.resid) %>%
  mutate(fitted = full.fitted)

ggplot(df.2, aes(fitted, resid, color = SUBTYPE)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 0) +
  xlab("Fitted Values") +
  ylab("Residual Values") +
  labs(color = "Subtype") +
  theme_bw()
ggsave("Plots/residual_vs_fitted_by_subtype.png", width = 6, height = 4)

ggplot(df.2, aes(fitted, resid, color = GENETIC_ANCESTRY_LABEL)) +
  geom_point(alpha = 0.75) +
  geom_hline(yintercept = 0) +
  xlab("Fitted Values") +
  ylab("Residual Values") +
  labs(color = "Genetic Ancestry Label") +
  theme_bw()
ggsave("Plots/residual_vs_fitted_by_ancestry.png", width = 6, height = 4)

# QQ Plot
png("Plots/QQ_plot.png", width = 6, height = 4, units = "in", res = 300)
qqnorm(full.resid); qqline(full.resid, col = "red")
dev.off()
