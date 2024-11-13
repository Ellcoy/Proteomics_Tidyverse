# Set seed for reproducibility
set.seed(42)

# Original sample data
proteins <- paste0("Protein_", 1:100)  # 100 proteins
conditions <- c("Control", "Treatment")
replicates <- 1:5  # 5 replicates per condition

df_proteomics <- data.frame(
  ProteinID = rep(proteins, each = length(conditions) * length(replicates)),
  Condition = rep(rep(conditions, each = length(replicates)), times = length(proteins)),
  Replicate = rep(replicates, times = length(proteins) * length(conditions)),
  Intensity = rlnorm(length(proteins) * length(conditions) * length(replicates), meanlog = 4, sdlog = 0.5)
)

# Add couple of "standard" proteins with NA in the Condition column
standard_proteins <- paste0("Protein_X_", 101:102)  # Simulate 2 proteins as standards

# Create a dataframe for these proteins with NA in Condition
df_standard <- data.frame(
  ProteinID = rep(standard_proteins, each = length(replicates) * length(conditions)),
  Condition = rep(NA, length(standard_proteins) * length(replicates) * length(conditions)),
  Replicate = rep(replicates, times = length(standard_proteins) * length(conditions)),
  Intensity = rlnorm(length(standard_proteins) * length(replicates) * length(conditions), meanlog = 4, sdlog = 0.5)
)

# Combine the standard proteins with the original data
df_proteomics <- bind_rows(df_proteomics, df_standard)

# Display the first few rows of the updated dataset
head(df_proteomics)


write_csv(df_proteomics,"proteomics_sample_data.csv")
