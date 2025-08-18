# Clear the workspace
rm(list = ls())

# Load required libraries
library(dplyr)      # For data manipulation
library(tidyr)      # For drop_na() function
library(mgcv)       # For GAM models
library(MuMIn)      # For all-subset regression (dredge)
library(visreg)     # For partial effect plots
library(ggplot2)    # For plotting
library(MASS)       # For Negative Binomial model
library(car)        # For variance inflation factor (VIF) and Anova tests
library(lmtest)     # For likelihood ratio tests

# -----------------------------
# PART 1: Data Preparation
# -----------------------------

# Read in the dataset (adjust the file path as needed)
data <- read.csv("Desktop/Virus_fasta/complete_Virus_model_unique_sci_name_aggregated.csv", 
                 stringsAsFactors = FALSE)

# Convert selected columns to factors and create virus_family
data <- data %>%
  mutate(
    infected_host_sample_size = factor(infected_host_sample_size),
    actual_host_sample_size   = factor(actual_host_sample_size),
    hosts                     = factor(hosts),
    location                  = factor(location),
    sample_type               = factor(sample_type),
    virus_family              = factor(family),
    host_order                = factor(host_order),
    host_habitat              = factor(host_habitat)
  )

# -----------------------------
# PART 2: Host-Level Data Aggregation and GLMs
# -----------------------------

# Aggregate data at the host level
host_data <- data %>%
  group_by(hosts, host_order, actual_host_sample_size, host_habitat, location, sample_type) %>%
  summarise(
    virus_species_richness = n_distinct(sci_name),
    virus_abundance        = sum(RPM),
    host_sample_size       = n(),
    .groups = "drop"
  )

# Compute a proxy for host geographic distribution: number of unique locations per host
geo_dist <- data %>%
  group_by(hosts) %>%
  summarise(
    host_geographic_distribution = n_distinct(location),
    .groups = "drop"
  )

# Merge geographic distribution into host_data
host_data <- left_join(host_data, geo_dist, by = "hosts")

# -----------------------------
# Boxplot: Virus Richness by Host
# -----------------------------

# Define color mapping for hosts
host_colors <- c('Chicken' = '#FF33A6',
                 'Bat'     = '#7B68EE',
                 'Dog'     = '#4D7D73',
                 'Goat'    = '#884936',
                 'Egret'   = '#A52A2A',
                 'Cat'     = '#FFD433',
                 'Sheep'   = '#338CFF',
                 'Rodent'  = 'grey',
                 'Lizard'  = '#33FFD4',
                 'Squirrel'='#FFC0CB',
                 'Pig'     = '#B833FF')

# Save the plot to a PNG file with increased dimensions
png("./Desktop/Virus_fasta/1Virus_Species_Richness_by_Host.png", width = 1200, height = 800, res = 150)

# Boxplot of virus species richness by host
boxplot(virus_species_richness ~ hosts, data = host_data,
        col = host_colors[levels(host_data$hosts)],
        main = "Virus Species Richness by Host",
        xlab = "Host",
        ylab = "Virus Species Richness")

dev.off()

# -----------------------------
# All-Subset Regression for Host-Level Responses
# -----------------------------

# Remove rows with missing values using tidyr's drop_na()
host_data_clean <- drop_na(host_data)

# Fit Negative Binomial model for total virus richness
nb_model_richness <- glm.nb(virus_species_richness ~  host_sample_size + host_order + host_habitat +
                              location + sample_type + 1,
                            data = host_data_clean,
                            na.action = na.fail)

# Summary of the Negative Binomial model
summary(nb_model_richness)

# Check overdispersion
overdisp_rich <- sum(residuals(nb_model_richness, type = "pearson")^2) / 
  nb_model_richness$df.residual
overdisp_rich

# All-subset regression (dredge) and best model selection
dredge_richness <- dredge(nb_model_richness)
best_model_richness <- get.models(dredge_richness, subset = 1)[[1]]
summary(best_model_richness)

# Likelihood Ratio Test
lrtest(best_model_richness)

# Wald Test for individual predictors
Anova(best_model_richness, type = "III")


# Assuming the 'host_data_clean' dataframe contains the relevant columns: 'virus_species_richness', 'host_sample_size', and 'host_order'


# -----------------------------
# Diagnostic Plots for the Negative Binomial Model
# -----------------------------

# Diagnostic plots for the best model
par(mfrow = c(2, 2))
plot(best_model_richness)


# Perform the Kruskal-Wallis test for virus richness across host_order
kruskal_test <- kruskal.test(virus_species_richness ~ host_order, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Load necessary libraries
library(ggplot2)
library(ggsignif)  # For adding significance bars

# Check the distribution of host orders
table(host_data_clean$host_order)

# Perform the Kruskal-Wallis test for virus richness across host_order
kruskal_test <- kruskal.test(virus_species_richness ~ host_order, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Create the boxplot to visualize virus richness across host orders
plot <- ggplot(host_data_clean, aes(x = host_order, y = virus_species_richness, fill = host_order)) +
  geom_boxplot() +  # Create the boxplot
  labs(
    title = "Virus Species Richness by Host Order", 
    x = "Host Order", 
    y = "Virus Species Richness"
  ) +
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(angle = 0, size = 20, face = "bold", hjust = 1),  # Increase size and adjust angle for x-axis labels
    axis.text.y = element_text(size = 20, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis titles size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",  # Remove legend for this case
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border to the panel (axis frame)
  ) +
  scale_fill_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936'))  # Customize colors for different host orders

# Add significance bar (adjust the comparisons based on your specific results)
# You may want to ensure that the pairwise comparisons have sufficient observations
plot <- plot + 
  geom_signif(
    comparisons = list(c("Herbivores", "Insectivores")),  # Specify groups to compare
    map_signif_level = TRUE,  # Add significance levels
    textsize = 8,  # Adjust text size of the p-value
    size = 1,  # Size of the significance bar
    y_position = max(host_data_clean$virus_species_richness) + 2  # Position the significance bar above the boxes
  )

# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_Species_Richness_by_Host_Order.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, you can save it in other formats, like PDF, by changing the file extension to .pdf
ggsave("./Desktop/Virus_fasta/Virus_Species_Richness_by_Host_Order.pdf", plot = plot, width = 12, height = 8)



# Load necessary libraries
library(ggplot2)
library(dplyr)
library(broom)  # For tidy() and glance() functions to extract model statistics

# Fit linear models for each host_order and extract p-values and R-squared values
model_stats <- host_data_clean %>%
  group_by(host_order) %>%
  do({
    model <- lm(log10(virus_species_richness) ~ host_sample_size, data = .)  # Fit model
    tidy_model <- tidy(model)  # Extract coefficients and p-values
    glance_model <- glance(model)  # Extract R-squared
    data.frame(
      host_order = unique(.$host_order),
      p_value = tidy_model$p.value[2],  # p-value for host_sample_size
      r_squared = glance_model$r.squared  # R-squared value
    )
  }) %>%
  ungroup()

# Print the model statistics
print(model_stats)

# Adjust color palette to match the number of host orders
plot <- ggplot(host_data_clean, aes(x = host_sample_size, y = virus_species_richness, color = host_order)) +
  geom_point(size = 3, alpha = 0.6) +  # Scatter plot for individual points with transparency
  geom_smooth(method = "lm", aes(group = host_order), se = TRUE, linetype = "solid", size = 1) +  # Add regression lines with confidence bands
  scale_y_log10() +  # Log-transform the y-axis
  labs(
    title = "Total Virus Richness by Host Order and Sample Size (Log Scale)", 
    x = "Host Sample Size", 
    y = "Virus Species Richness (Log Scale)",
    color = "Host Order"  # Add a legend title for host_order
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.title = element_text(size = 14, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 12)  # Customize legend text
  ) + 
  scale_color_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936', '#FFB6C1')) +  # Customize colors for 5 host orders
  geom_text(
    data = model_stats, 
    aes(
      x = Inf, y = Inf,  # Position annotations at the top-right corner
      label = paste0("R² = ", round(r_squared, 2), "\np = ", round(p_value, 3)),  # Add R² and p-value
      hjust = 1.1, vjust = 1.5,  # Adjust text position
      size = 5,  # Text size
      color = "black"  # Text color
    )
  )


# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_Richness_by_Host_Order_and_Sample_Size_Log_Scale_with_Stats.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/Virus_Richness_by_Host_Order_and_Sample_Size_Log_Scale_with_Stats.pdf", plot = plot, width = 12, height = 8)





# -----------------------------
# Partial Effect Plots
# -----------------------------

# Fit Poisson model for visualization purposes
model_for_vis <- glm(virus_species_richness ~ host_sample_size + host_order + 1, 
                     family = poisson(link = "log"), data = host_data_clean, 
                     na.action = na.fail)

# Save the partial effect plot for host_sample_size
png("./Desktop/Virus_fasta/Host_Sample_Size_vs_virus_richness.png", width = 1200, height = 800, res = 150)
visreg(model_for_vis, "host_sample_size", 
       main = "Partial Effect: Host Sample Size (Exploratory)",
       line = list(col = "blue", lwd = 1.5),         
       points = list(col = "darkblue", pch = 16),    
       fill = list(col = rgb(0.1, 0.1, 0.9, alpha = 0.2)))  
dev.off()

# Save the partial effect plot for host_order

png("./Desktop/Virus_fasta/Feedinghabitat_Effect_Host_Order.png", width = 1200, height = 800, res = 150)
visreg(model_for_vis, "host_order", 
       main = "Partial Effect: Host Order (Exploratory)",
       line = list(col = "red", lwd = 2),
       points = list(col = "darkred", pch = 16),
       fill = list(col = rgb(0.9, 0.1, 0.1, alpha = 0.2)))
dev.off()

# -----------------------------
# Boxplot: Virus Richness by Sample Type
# -----------------------------

# Define the sample type color map
sample_type_color_map <- c(
  "Liver"   = "blue",
  "Lung"    = "red",
  "Swab"    = "gold",
  "Spleen"  = "pink",
  "Oral"    = "purple",
  "Rectal"  = "green"
)

# Ensure that sample_type is a factor with levels in the desired order
host_data_clean$sample_type <- factor(host_data_clean$sample_type, 
                                      levels = names(sample_type_color_map))

# Save a boxplot of virus_species_richness by sample_type as a PNG file
png("./Desktop/Virus_fasta/Partial_Effect_sample_type_Boxplot.png", 
    width = 1200, height = 800, res = 150)
par(mar = c(7, 4, 4, 2) + 0.1)  # Adjust margins if necessary

boxplot(virus_species_richness ~ sample_type, data = host_data_clean,
        col = sample_type_color_map[levels(host_data_clean$sample_type)],
        main = "Virus Species Richness by Sample Type (Boxplot)",
        xlab = "Sample Type",
        ylab = "Virus Species Richness")
dev.off()

# -----------------------------
# Boxplot: Virus Richness by Host Order
# -----------------------------

# Save a boxplot of virus_species_richness by host_order as a PNG file
png("./Desktop/Virus_fasta/Partial_Effect_Host_Order_Boxplot.png", width = 1200, height = 800, res = 150)
par(mar = c(7, 4, 4, 2) + 0.1)  # Adjust margins if necessary

boxplot(virus_species_richness ~ host_order, data = host_data_clean,
        col = c("darkred", "#0E6251", "steelblue"),  # Customize colors for each category
        main = "Virus Species Richness by Host Order (Boxplot)",
        xlab = "Host Order",
        ylab = "Virus Species Richness")
dev.off()

# -----------------------------
# Partial Effect Plots for Total Virus Richness
# -----------------------------

# Visualizing partial effects for host factors
visreg(best_model_richness, "hosts", 
       main = "Partial Effect: Host Sample Size on Virus Richness")
visreg(best_model_richness, "location", 
       main = "Partial Effect: Host Order on Virus Richness")



library(MASS)
library(MuMIn)  # for dredge()

# Fit a negative binomial model for virus_abundance
full_model_abundance_nb <- glm.nb(
  virus_abundance ~ hosts + host_order + host_habitat +
    host_sample_size + host_geographic_distribution +
    location + sample_type,
  data = host_data_clean,
  na.action = na.fail
)

summary(full_model_abundance_nb)

# Check overdispersion
overdisp_abund_nb <- sum(residuals(full_model_abundance_nb, type="pearson")^2) /
  full_model_abundance_nb$df.residual
overdisp_abund_nb
dredge_abundance_nb <- dredge(full_model_abundance_nb, trace = FALSE)
head(dredge_abundance_nb)  # top models by AIC

best_model_abundance_nb <- get.models(dredge_abundance_nb, subset = 1)[[1]]
summary(best_model_abundance_nb)


# Fit Poisson model for visualization purposes
model_for_vis <- glm(virus_abundance ~  host_habitat + host_sample_size + 1, 
                     family = poisson(link = "log"), data = host_data_clean, 
                     na.action = na.fail)
# ----------------------------------------------
# BEAUTIFUL PARTIAL EFFECT PLOTS & SAVE TO PNG
# ----------------------------------------------

# 1. Partial Effect: Host Sample Size on Virus Abundance
png("./Desktop/Virus_fasta/Host_Sample_Size_on_Virus_Abundance.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_abundance_nb,
  "host_sample_size",
  main = "Partial Effect: Host Sample Size on Virus Abundance",
  line = list(col = "blue", lwd = 2),            # Customize line color & thickness
  points = list(col = "darkblue", pch = 16),     # Customize point color & symbol
  fill = list(col = rgb(0.1, 0.1, 0.9, alpha = 0.2))  # Transparent fill
)

dev.off()  # Closes the PNG device and saves the file

# 2. Partial Effect: Host Habitat on Virus Abundance
png("./Desktop/Virus_fasta/Host_Habitat_on_Virus_Abundance.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_abundance_nb,
  "host_habitat",
  main = "Partial Effect: Host Habitat on Virus Abundance",
  line = list(col = "red", lwd = 2),
  points = list(col = "darkred", pch = 16),
  fill = list(col = rgb(0.9, 0.1, 0.1, alpha = 0.2))
)

dev.off()

# Load necessary libraries
library(ggplot2)
library(ggsignif)  # For adding significance bars

# Check the distribution of host habitats
table(host_data_clean$host_habitat)

# Perform the Kruskal-Wallis test for log-transformed virus abundance across host_habitat
kruskal_test <- kruskal.test(log10(virus_abundance) ~ host_habitat, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Create the boxplot to visualize log-transformed virus abundance across host habitats
plot <- ggplot(host_data_clean, aes(x = host_habitat, y = log10(virus_abundance), fill = host_habitat)) +
  geom_boxplot() +  # Create the boxplot
  labs(
    title = "Virus Abundance by Host Habitat (Log Scale)", 
    x = "Host Habitat", 
    y = "Virus Abundance (Log Scale)"
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(angle = 0, size = 20, face = "bold", hjust = 1),  # Increase size and adjust angle for x-axis labels
    axis.text.y = element_text(size = 20, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",  # Remove legend for this case
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border to the panel (axis frame)
  ) + 
  scale_fill_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936'))  # Customize colors for different host habitats

# Add significance bar (adjust the comparisons based on your specific results)
# Ensure that the pairwise comparisons have sufficient observations
plot <- plot + 
  geom_signif(
    comparisons = list(c("peri-domesticated", "Wild")),  # Specify groups to compare
    map_signif_level = TRUE,  # Add significance levels
    textsize = 8,  # Adjust text size of the p-value
    size = 1,  # Size of the significance bar
    y_position = max(log10(host_data_clean$virus_abundance)) + 2  # Position the significance bar above the boxes
  )

# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_abundance_by_Host_habitat.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in other formats, like PDF, by changing the file extension to .pdf
ggsave("./Desktop/Virus_fasta/Virus_abundance_by_Host_habitat.pdf", plot = plot, width = 12, height = 8)





# Load necessary libraries
library(ggplot2)
library(dplyr)
library(broom)  # For tidy() and glance() functions to extract model statistics

# Fit linear models for each host_order and extract p-values and R-squared values
model_stats <- host_data_clean %>%
  group_by(host_habitat) %>%
  do({
    model <- lm(log10(virus_abundance) ~ host_sample_size, data = .)  # Fit model
    tidy_model <- tidy(model)  # Extract coefficients and p-values
    glance_model <- glance(model)  # Extract R-squared
    data.frame(
      host_habitat = unique(.$host_habitat),
      p_value = tidy_model$p.value[2],  # p-value for host_sample_size
      r_squared = glance_model$r.squared  # R-squared value
    )
  }) %>%
  ungroup()

# Print the model statistics
print(model_stats)

# Adjust color palette to match the number of host orders
plot <- ggplot(host_data_clean, aes(x = host_sample_size, y = virus_abundance, color = host_habitat)) +
  geom_point(size = 3, alpha = 0.6) +  # Scatter plot for individual points with transparency
  geom_smooth(method = "lm", aes(group = host_habitat), se = TRUE, linetype = "solid", size = 1) +  # Add regression lines with confidence bands
  scale_y_log10() +  # Log-transform the y-axis
  labs(
    title = "Total Virus Richness by Host Order and Sample Size (Log Scale)", 
    x = "Host Sample Size", 
    y = "Virus Species Richness (Log Scale)",
    color = "Host Order"  # Add a legend title for host_order
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.title = element_text(size = 14, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 12)  # Customize legend text
  ) + 
  scale_color_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936', '#FFB6C1')) +  # Customize colors for 5 host orders
  geom_text(
    data = model_stats, 
    aes(
      x = Inf, y = Inf,  # Position annotations at the top-right corner
      label = paste0("R² = ", round(r_squared, 2), "\np = ", round(p_value, 3)),  # Add R² and p-value
      hjust = 1.1, vjust = 1.5,  # Adjust text position
      size = 5,  # Text size
      color = "black"  # Text color
    )
  )


# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/virus_abundance_Host_habitat_and_Sample_Size_Log_Scale_with_Stats.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/virus_abundance_Host_habitat_and_Sample_Size_Log_Scale_with_Stats.pdf", plot = plot, width = 12, height = 8)





# Convert virus_host_count from factor to numeric
data$virus_host_count <- as.numeric(as.character(data$virus_host_count))


library(dplyr)

host_data <- data %>%
  group_by(hosts, host_order, host_habitat, location, sample_type) %>%
  summarise(
    virus_species_richness = n_distinct(sci_name),
    virus_abundance        = sum(RPM, na.rm = TRUE),
    host_sample_size       = n(),
    cross_species_count    = sum(virus_host_count > 1, na.rm = TRUE),  # newly created
    .groups = "drop"
  )


geo_dist <- data %>%
  group_by(hosts) %>%
  summarise(
    host_geographic_distribution = n_distinct(location),
    .groups = "drop"
  )

# Merge into host_data
host_data <- left_join(host_data, geo_dist, by = "hosts")

library(tidyr)

host_data_clean <- drop_na(host_data)
nrow(host_data_clean)

summary(host_data)
colSums(is.na(host_data))

names(host_data_clean)


library(MASS)

full_model_cross_nb <- glm.nb(
  cross_species_count ~ virus_species_richness + virus_abundance + 
    hosts + host_order + host_sample_size + host_geographic_distribution +
    location + sample_type,
  data = host_data_clean,
  na.action = na.fail
)

summary(full_model_cross_nb)



library(MuMIn)

dredge_cross_nb <- dredge(full_model_cross_nb, trace = FALSE)
head(dredge_cross_nb)
best_model_cross_nb <- get.models(dredge_cross_nb, subset = 1)[[1]]
summary(best_model_cross_nb)


# Fit the best-fit model (from dredge)
best_model_cross_nb <- glm.nb(cross_species_count ~ host_feeding + host_habitat + 
                                host_sample_size + virus_abundance + 1, 
                              data = host_data_clean, na.action = na.fail)

# Extract null deviance from the best-fit model
null_deviance_full_model <- best_model_cross_nb$null.deviance
print(null_deviance_full_model)  # Should print 234.138

# Function to calculate % deviance explained for a given predictor
calculate_deviance_explained <- function(null_formula, full_formula, predictor_name) {
  # Fit the null model (without the predictor)
  null_model <- glm.nb(null_formula, data = host_data_clean, na.action = na.fail)
  
  # Fit the full model (with the predictor)
  full_model <- glm.nb(full_formula, data = host_data_clean, na.action = na.fail)
  
  # Perform likelihood ratio test (LRT)
  lrt <- anova(null_model, full_model, test = "LRT")
  print(lrt)  # Print LRT results
  
  # Calculate deviance difference
  deviance_diff <- lrt$`Chisq`[2]
  
  # Calculate percentage of deviance explained
  percent_deviance_explained <- (deviance_diff / null_deviance_full_model) * 100
  print(paste("Percentage of deviance explained by", predictor_name, ":", percent_deviance_explained))
}

# Calculate % deviance explained for each predictor

# 1. For host_feeding
calculate_deviance_explained(
  null_formula = cross_species_count ~ host_habitat + host_sample_size + virus_abundance + 1,
  full_formula = cross_species_count ~ host_feeding + host_habitat + host_sample_size + virus_abundance + 1,
  predictor_name = "host_feeding"
)

# 2. For host_habitat
calculate_deviance_explained(
  null_formula = cross_species_count ~ host_feeding + host_sample_size + virus_abundance + 1,
  full_formula = cross_species_count ~ host_feeding + host_habitat + host_sample_size + virus_abundance + 1,
  predictor_name = "host_habitat"
)

# 3. For host_sample_size
calculate_deviance_explained(
  null_formula = cross_species_count ~ host_feeding + host_habitat + virus_abundance + 1,
  full_formula = cross_species_count ~ host_feeding + host_habitat + host_sample_size + virus_abundance + 1,
  predictor_name = "host_sample_size"
)

# 4. For virus_abundance
calculate_deviance_explained(
  null_formula = cross_species_count ~ host_feeding + host_habitat + host_sample_size + 1,
  full_formula = cross_species_count ~ host_feeding + host_habitat + host_sample_size + virus_abundance + 1,
  predictor_name = "virus_abundance"
)





# (F) Partial Effect Plots for Number of Cross-Species Transmitted Viruses ---------
# 1. Partial Effect: Host Sample Size on Virus Abundance
png("./Desktop/Virus_fasta/cross_species_count_effect_distribution_cross.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_cross_nb,
  "host_order",
  main = "Partial Effect: Host Feeding habit on Virus cross transmission",
  line = list(col = "red", lwd = 2),            # Customize line color & thickness
  points = list(col = "darkred", pch = 16),     # Customize point color & symbol
  fill = list(col = rgb(0.1, 0.1, 0.9, alpha = 0.2))  # Transparent fill
)

dev.off()  # Closes the PNG device and saves the file

# 2. Partial Effect: Host Habitat on Virus Abundance
png("./Desktop/Virus_fasta/Partial_Cross_Effect_on_cross_trans.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_cross_nb,
  "sample_type",
  main = "Partial Effect: Sample Type on Virus cross transmission",
  line = list(col = "blue", lwd = 2),
  points = list(col = "darkblue", pch = 16),
  fill = list(col = rgb(0.9, 0.1, 0.1, alpha = 0.2))
)

dev.off()

# 2. Partial Effect: Host Habitat on Virus Abundance
png("./Desktop/Virus_fasta/Effect_virus_species_richness_on_cross_trans.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_cross_nb,
  "virus_species_richness",
  main = "Partial Effect: virus abundance on Virus cross transmission",
  line = list(col = "red", lwd = 2),
  points = list(col = "darkred", pch = 16),
  fill = list(col = rgb(0.9, 0.1, 0.1, alpha = 0.2))
)

dev.off()


# Adjust color palette to match the number of host orders
plot <- ggplot(host_data_clean, aes(x = host_sample_size, y = virus_species_richness, color = host_order)) +
  geom_point(size = 3, alpha = 0.6) +  # Scatter plot for individual points with transparency
  geom_smooth(method = "lm", aes(group = host_order), se = TRUE, linetype = "solid", size = 1) +  # Add regression lines with confidence bands
  scale_y_log10() +  # Log-transform the y-axis
  labs(
    title = "Total Virus Richness by Host Order and Sample Size (Log Scale)", 
    x = "Host Sample Size", 
    y = "Virus Species Richness (Log Scale)",
    color = "Host Order"  # Add a legend title for host_order
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.title = element_text(size = 14, face = "bold"),  # Customize legend title
    legend.text = element_text(size = 12)  # Customize legend text
  ) + 
  scale_color_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936', '#FFB6C1')) +  # Customize colors for 5 host orders
  geom_text(
    data = model_stats, 
    aes(
      x = Inf, y = Inf,  # Position annotations at the top-right corner
      label = paste0("R² = ", round(r_squared, 2), "\np = ", round(p_value, 3)),  # Add R² and p-value
      hjust = 1.1, vjust = 1.5,  # Adjust text position
      size = 5,  # Text size
      color = "black"  # Text color
    )
  )


# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_Richness_by_Host_Order_and_Sample_Size_Log_Scale_with_Stats.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/Virus_Richness_by_Host_Order_and_Sample_Size_Log_Scale_with_Stats.pdf", plot = plot, width = 12, height = 8)







# 2. Partial Effect: Host Habitat on Virus Abundance
png("./Desktop/Virus_fasta/Partial_Effect_sample_type_on_cross_trans.png",
    width = 1200, height = 800, res = 150)

visreg(
  best_model_cross_nb,
  "sample_type",
  main = "Partial Effect: sample type on Virus cross transmission",
  line = list(col = "red", lwd = 2),
  points = list(col = "darkred", pch = 16),
  fill = list(col = rgb(0.9, 0.1, 0.1, alpha = 0.2))
)

dev.off()

overdisp_cross_nb <- sum(residuals(best_model_cross_nb, type="pearson")^2) /
  best_model_cross_nb$df.residual
overdisp_cross_nb


# Load necessary libraries
library(ggplot2)
library(ggsignif)  # For adding significance bars

# Check the distribution of sample types
table(host_data_clean$sample_type)

# Remove rows where cross_species_count is missing or zero before log transformation
host_data_clean <- host_data_clean %>%
  filter(!is.na(cross_species_count) & cross_species_count > 0)

# Perform the Kruskal-Wallis test for log-transformed virus abundance across sample_type
kruskal_test <- kruskal.test(log10(cross_species_count) ~ sample_type, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Create the boxplot to visualize log-transformed virus abundance across sample types
plot <- ggplot(host_data_clean, aes(x = sample_type, y = log10(cross_species_count), fill = sample_type)) +
  geom_boxplot() +  # Create the boxplot
  labs(
    title = "Virus Abundance by Sample Type (Log Scale)", 
    x = "Sample Type", 
    y = "Virus Abundance (Log Scale)"
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(angle = 0, size = 20, face = "bold", hjust = 1),  # Increase size and adjust angle for x-axis labels
    axis.text.y = element_text(size = 20, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",  # Remove legend for this case
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border to the panel (axis frame)
  ) + 
  scale_fill_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936', '#FFB6C1')) +  # Customize colors for different sample types
  geom_signif(
    comparisons = list(c("Oral", "Liver"), c("Oral", "Rectal"), c("Oral", "Swab"), c("Oral", "Lung")),  # Specify groups to compare
    map_signif_level = TRUE,  # Add significance levels
    textsize = 8,  # Adjust text size of the p-value
    size = 1,  # Size of the significance bar
    y_position = max(log10(host_data_clean$cross_species_count)) + 2  # Position the significance bar above the boxes
  )

# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_abundance_by_Sample_type.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/Virus_abundance_by_Sample_type.pdf", plot = plot, width = 12, height = 8)






# Load necessary libraries
library(ggplot2)
library(ggsignif)  # For adding significance bars

# Check the distribution of host orders
table(host_data_clean$host_order)

# Correct typo in the comparisons, if necessary
# Perform the Kruskal-Wallis test for log-transformed virus abundance across host_order
kruskal_test <- kruskal.test(log10(cross_species_count) ~ host_order, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Create the boxplot to visualize log-transformed virus abundance across host orders
plot <- ggplot(host_data_clean, aes(x = host_order, y = log10(cross_species_count), fill = host_order)) +
  geom_boxplot() +  # Create the boxplot
  labs(
    title = "Virus Abundance by Host Order (Log Scale)", 
    x = "Host Order", 
    y = "Virus Abundance (Log Scale)"
  ) + 
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(angle = 0, size = 20, face = "bold", hjust = 1),  # Increase size and adjust angle for x-axis labels
    axis.text.y = element_text(size = 20, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(),  # Remove grid lines
    legend.position = "none",  # Remove legend for this case
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border to the panel (axis frame)
  ) + 
  scale_fill_manual(values = c('#FF33A6', '#7B68EE', '#4D7D73', '#884936')) +  # Customize colors for different host orders
  geom_signif(
    comparisons = list(c("Carnivores", "Herbivores"), c("Insectivores", "Carnivores"), c("Omnivores", "Herbivores"), c("Herbivores", "Carnivores")),  # Specify correct groups to compare
    map_signif_level = TRUE,  # Add significance levels
    textsize = 8,  # Adjust text size of the p-value
    size = 1,  # Size of the significance bar
    y_position = max(log10(host_data_clean$cross_species_count)) + 2  # Position the significance bar above the boxes
  )

# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_cross_by_feeding_type.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/Virus_cross_by_feeding_type.pdf", plot = plot, width = 12, height = 8)



# Load necessary libraries
library(ggplot2)
library(ggsignif)  # For adding significance bars
library(MASS)      # For negative binomial regression

# Check the distribution of virus richness and cross-species transmission
table(host_data_clean$virus_species_richness)

# Perform the Kruskal-Wallis test for log-transformed virus richness and cross-species transmission
kruskal_test <- kruskal.test(log10(cross_species_count) ~ virus_species_richness, data = host_data_clean)

# Print the result of the test
print(kruskal_test)

# Fit a negative binomial regression model to test the relationship
nb_model <- glm.nb(cross_species_count ~ virus_species_richness + 1, data = host_data_clean, na.action = na.fail)

# Print summary of the model
summary(nb_model)

# Create the plot to visualize the relationship between virus richness and cross-species transmission
plot <- ggplot(host_data_clean, aes(x = virus_species_richness, y = log10(cross_species_count))) +
  geom_point(aes(color = virus_species_richness), size = 3, alpha = 0.6) +  # Scatter plot for individual points with transparency
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid", size = 1) +  # Add regression line with confidence bands
  labs(
    title = "Virus Richness vs Cross-Species Transmission (Log Scale)",
    x = "Virus Species Richness",
    y = "Cross-Species Transmission (Log Scale)"
  ) +
  theme_minimal() +  # Minimal theme for a cleaner look
  theme(
    axis.text.x = element_text(size = 14, face = "bold"),  # Increase size of x-axis labels
    axis.text.y = element_text(size = 14, face = "bold"),  # Increase size of y-axis labels
    axis.title = element_text(size = 16),  # Increase axis title size
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size and make it bold
    panel.grid = element_blank(),  # Remove grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # Add border to the panel (axis frame)
  )

# Save the plot as a PNG file with increased dimensions
ggsave("./Desktop/Virus_fasta/Virus_richness_vs_cross_species_transmission.png", plot = plot, width = 12, height = 8, dpi = 150)

# Optionally, save it in PDF format
ggsave("./Desktop/Virus_fasta/Virus_richness_vs_cross_species_transmission.pdf", plot = plot, width = 12, height = 8)

























# =============================
# PART 3: Virus-Level Data Aggregation and GAMs
# =============================

# Derive multi_organ_distribution from sample_type.
# For each virus (sci_name) in a given host, flag if it is detected in >1 sample type.
virus_host_data <- data %>%
  group_by(sci_name, hosts) %>%
  summarise(unique_sample_types      = n_distinct(sample_type),
            virus_abundance          = sum(read_count),
            multi_organ_distribution = ifelse(unique_sample_types > 1, 1, 0),
            virus_family             = first(virus_family),
            .groups = "drop")

# Aggregate to a virus-level summary across hosts.
# A virus is considered cross-species transmitted if it is found in >1 host.
virus_summary <- virus_host_data %>%
  group_by(sci_name) %>%
  summarise(total_abundance          = sum(virus_abundance),
            host_count               = n_distinct(hosts),
            multi_organ_distribution = as.integer(any(multi_organ_distribution == 1)),
            virus_family             = first(virus_family),
            .groups = "drop") %>%
  mutate(cross_species_transmission_status = ifelse(host_count > 1, 1, 0))

# (G) GAM for Cross-Species Virus Transmission Potential --------------------
gam_model <- gam(cross_species_transmission_status ~ host_count + total_abundance +
                   multi_organ_distribution + s(virus_family, bs = "re"),
                 family = binomial(link = "logit"),
                 data = virus_summary,
                 method = "REML")
summary(gam_model)

# Partial effect plots for the GAM model:
visreg(gam_model, "host_count", 
       main = "Partial Effect: Host Count on Transmission Potential")
visreg(gam_model, "total_abundance", 
       main = "Partial Effect: Total Abundance on Transmission Potential")
visreg(gam_model, "multi_organ_distribution", 
       main = "Partial Effect: Multi-Organ Distribution on Transmission Potential")

# (H) GAM with Nested Random Effects -----------------------------------------
# If you have a grouping variable for host family, you can nest virus_family within it.
# Here, if host_family is not directly available, you might use hosts as a proxy.
# (This example uses a placeholder; replace with your actual host family variable if available.)
virus_summary$host_family <- virus_summary$host_count  # Placeholder

# (H) GAM with Nested Random Effects -----------------------------------------
# If you do not have an actual host_family variable, create a placeholder factor.
# Ensure that virus_summary has a proper 'host_family' column.
# If you don't have one, create a placeholder factor.
if(!"host_family" %in% names(virus_summary)) {
  virus_summary$host_family <- factor(rep("unknown", nrow(virus_summary)))
}

# Create a new combined grouping variable for the nested random effect.
virus_summary$group_nested <- with(virus_summary, interaction(virus_family, host_family, drop = TRUE))

# Check that the new variable has the correct length:
stopifnot(nrow(virus_summary) == length(virus_summary$group_nested))

# Now, fit the nested GAM using this combined factor:
gam_model_nested <- gam(cross_species_transmission_status ~ host_count + total_abundance +
                          multi_organ_distribution + s(group_nested, bs = "re"),
                        family = binomial(link = "logit"),
                        data = virus_summary,
                        method = "REML")
summary(gam_model_nested)

# Partial effect plots for the nested GAM model:
visreg(gam_model_nested, "host_count", 
       main = "Nested GAM: Partial Effect of Host Count")
visreg(gam_model_nested, "total_abundance", 
       main = "Nested GAM: Partial Effect of Total Abundance")
visreg(gam_model_nested, "multi_organ_distribution", 
       main = "Nested GAM: Partial Effect of Multi-Organ Distribution")

