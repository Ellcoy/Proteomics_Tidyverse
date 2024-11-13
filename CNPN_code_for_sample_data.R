#Install Tidyverse Package
# install.packages("tidyverse")
# Load necessary libraries
library(tidyverse)

df_proteomics <- read_csv("proteomics_sample_data.csv")

# Step 1: Filter out proteins with "X" in the name
df_proteomics_filtered <- df_proteomics %>%
  filter(!str_detect(ProteinID, "X"))

# Step 2: Calculate mean intensity by Protein and Condition for simplicity
mean_intensity <- df_proteomics_filtered %>%
  group_by(ProteinID, Condition) %>%
  summarize(mean_intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Step 3: Calculate p-values for each protein (between conditions)
p_values <- df_proteomics_filtered %>%
  group_by(ProteinID) %>%
  summarize(
    p.value = t.test(Intensity[Condition == "Control"], Intensity[Condition == "Treatment"], var.equal = FALSE)$p.value,
    .groups = "drop"
  )

# Step 4: Merge summarized mean intensity and p-values, calculate log2 Fold Change and significance
mean_intensity_wide <- mean_intensity %>%
  pivot_wider(names_from = Condition, values_from = mean_intensity) %>%
  inner_join(p_values, by = "ProteinID") %>%
  mutate(
    log2FC = log2(Treatment / Control),
    significance = if_else(p.value < 0.05 & log2FC > 1, "Upregulated", 
                           if_else(p.value < 0.05 & log2FC < -1, "Downregulated", "Not Significant"))
  )

# Step 5: Create the volcano plot with cutoff lines and color coding
volcano_plot <- ggplot(mean_intensity_wide, aes(x = log2FC, y = -log10(p.value))) +
  geom_point(aes(color = significance)) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot of Simulated Data", x = "Log2 Fold Change", y = "-Log10(p-value)") +
  # Add labels for the most significant upregulated and downregulated proteins
  geom_text(data = mean_intensity_wide %>% filter(significance != "Not Significant"), aes(label = ProteinID), vjust = -1, size = 3, color = "black")

volcano_plot


# Save the volcano plot
ggsave("volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# Step 6: Create a detailed table with each replicate and the summary results
detailed_table <- df_proteomics_filtered %>%
  left_join(mean_intensity_wide %>% select(ProteinID, log2FC, p.value, significance), by = "ProteinID") %>%
  rename(p_value = p.value) %>%
  arrange(desc(significance), p_value) %>%
  select(ProteinID, Condition, Replicate, Intensity, log2FC, p_value, significance)

# View the detailed table
print(detailed_table)

# Optionally, export the detailed table
write_csv(detailed_table, "detailed_proteomics_table.csv")
