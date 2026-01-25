# =============================================================================
# Dormitory Culture Analysis - Main Analysis Script
# =============================================================================
# Authors: Tsubasa Sato, Natsumi Negoro, Misa LoPresti, Misaki Iio
# Affiliation: Keio University
# Date: 2025
# 
# This script performs the complete analysis for the dormitory culture study.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

# Load required packages
required_packages <- c(
  "lme4",      # Linear mixed-effects models

"vegan",     # PERMANOVA analysis
  "boot",      # Bootstrap resampling
  "ggplot2",   # Visualization
  "dplyr",     # Data manipulation
  "tidyr",     # Data tidying
  "readr",     # Data import
  "patchwork", # Combining plots
  "ggpubr",    # Publication-ready plots
  "RColorBrewer" # Color palettes
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load packages
lapply(required_packages, library, character.only = TRUE)

# Set seed for reproducibility
set.seed(42)

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# -----------------------------------------------------------------------------

# Note: Due to privacy considerations, raw data is not included in this repository.
# The following code assumes data is loaded from a CSV file with the structure:
# - id: Participant ID
# - building: Building name (Basil, Turmeric, Rosemary, Paprika)
# - sex: Participant sex (Male, Female)
# - resident_type: New or Incumbent
# - survey_month: Survey month (May, June, July, October, November)
# - Q1-Q8: Response scores (1-5 Likert scale)

# Example data structure (replace with actual data loading)
# data <- read_csv("data/survey_data.csv")

# For demonstration, create simulated data structure
create_demo_data <- function() {
  set.seed(42)
  
  buildings <- c("Basil", "Turmeric", "Rosemary", "Paprika")
  building_types <- c("Male-only", "Female-only", "Co-ed", "Co-ed")
  n_per_building <- c(53, 29, 29, 42)
  
  data_list <- list()
  
  for (i in seq_along(buildings)) {
    building <- buildings[i]
    n <- n_per_building[i]
    
    # Determine sex based on building type
    if (building == "Basil") {
      sex <- rep("Male", n)
    } else if (building == "Turmeric") {
      sex <- rep("Female", n)
    } else {
      sex <- sample(c("Male", "Female"), n, replace = TRUE)
    }
    
    # Create participant data
    df <- data.frame(
      id = paste0(building, "_", 1:n),
      building = building,
      building_type = building_types[i],
      sex = sex,
      resident_type = sample(c("New", "Incumbent"), n, replace = TRUE, prob = c(0.4, 0.6))
    )
    
    # Add Q1-Q8 scores with building effects
    building_effects <- list(
      Basil = c(0.1, 0.3, 0.1, 0, 0.5, 0.1, 0, 0.3),
      Turmeric = c(0.2, -0.2, 0.2, 0, 0.6, 0.1, 0, 0.4),
      Rosemary = c(0, 0.1, -0.1, 0, -0.3, -0.2, 0, -0.2),
      Paprika = c(0.1, 0, 0.2, 0, 0.1, 0.1, 0, 0)
    )
    
    for (q in 1:8) {
      base_mean <- 3.5
      effect <- building_effects[[building]][q]
      df[[paste0("Q", q)]] <- pmin(5, pmax(1, round(rnorm(n, base_mean + effect, 0.8))))
    }
    
    data_list[[i]] <- df
  }
  
  bind_rows(data_list)
}

# Load or create data
# data <- read_csv("data/survey_data.csv")  # Uncomment for actual data
data <- create_demo_data()  # Demo data for testing

# Data summary
cat("=== Data Summary ===\n")
cat("Total observations:", nrow(data), "\n")
cat("Buildings:\n")
print(table(data$building))
cat("\nResident types:\n
")
print(table(data$resident_type))

# -----------------------------------------------------------------------------
# 2. ICC Analysis (RQ1 & RQ2)
# -----------------------------------------------------------------------------

cat("\n=== ICC Analysis ===\n")

# Function to calculate ICC
calculate_icc <- function(data, outcome_var, group_var = "building") {
  formula <- as.formula(paste(outcome_var, "~ 1 + (1 |", group_var, ")"))
  model <- lmer(formula, data = data, REML = TRUE)
  
  var_components <- as.data.frame(VarCorr(model))
  var_between <- var_components$vcov[var_components$grp == group_var]
  var_within <- var_components$vcov[var_components$grp == "Residual"]
  
  icc <- var_between / (var_between + var_within)
  
  # Likelihood ratio test
  model_null <- lm(as.formula(paste(outcome_var, "~ 1")), data = data)
  lrt <- anova(model, model_null)
  
  list(
    icc = icc,
    var_between = var_between,
    var_within = var_within,
    model = model,
    lrt_chisq = ifelse(nrow(lrt) > 1, lrt$Chisq[2], NA),
    lrt_p = ifelse(nrow(lrt) > 1, lrt$`Pr(>Chisq)`[2], NA)
  )
}

# Calculate ICC for each item
items <- paste0("Q", 1:8)
item_names <- c(
  "Communication ease",
  "Barriers with opposite sex",
  "Unit friendships",
  "Room cleanliness",
  "Event participation willingness",
  "Interest in presentations",
  "Self-disclosure",
  "Building attachment"
)

icc_results <- data.frame(
  item = items,
  name = item_names,
  icc = NA,
  lrt_chisq = NA,
  lrt_p = NA
)

for (i in seq_along(items)) {
  result <- calculate_icc(data, items[i])
  icc_results$icc[i] <- result$icc
  icc_results$lrt_chisq[i] <- result$lrt_chisq
  icc_results$lrt_p[i] <- result$lrt_p
}

# Calculate overall ICC (mean of Q1-Q8)
data$mean_score <- rowMeans(data[, items])
overall_icc <- calculate_icc(data, "mean_score")

cat("Overall ICC:", round(overall_icc$icc, 3), "\n\n")
cat("Item-specific ICCs:\n")
print(icc_results[order(-icc_results$icc), ])

# Bootstrap CI for overall ICC
bootstrap_icc <- function(data, indices) {
  d <- data[indices, ]
  calculate_icc(d, "mean_score")$icc
}

boot_results <- boot(data, bootstrap_icc, R = 1000)
boot_ci <- boot.ci(boot_results, type = "perc")

cat("\nBootstrap 95% CI for overall ICC:", 
    round(boot_ci$percent[4], 3), "-", round(boot_ci$percent[5], 3), "\n")

# -----------------------------------------------------------------------------
# 3. PERMANOVA Analysis
# -----------------------------------------------------------------------------

cat("\n=== PERMANOVA Analysis ===\n")

# Prepare response matrix
response_matrix <- as.matrix(data[, items])

# Aitchison distance (compositional data)
# Add small constant to avoid zeros
response_matrix_adj <- response_matrix + 0.5
response_clr <- log(response_matrix_adj) - rowMeans(log(response_matrix_adj))
dist_aitchison <- dist(response_clr)

# PERMANOVA by building
permanova_building <- adonis2(dist_aitchison ~ building, data = data, permutations = 999)
cat("PERMANOVA by building:\n")
print(permanova_building)

# PERMANOVA by sex
permanova_sex <- adonis2(dist_aitchison ~ sex, data = data, permutations = 999)
cat("\nPERMANOVA by sex:\n")
print(permanova_sex)

# -----------------------------------------------------------------------------
# 4. Cultural Similarity Analysis (RQ3)
# -----------------------------------------------------------------------------

cat("\n=== Cultural Similarity Analysis ===\n")

# Function to calculate cosine similarity
cosine_similarity <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2))) * 100
}

# Calculate incumbent profiles for each building
incumbent_profiles <- data %>%
  filter(resident_type == "Incumbent") %>%
  group_by(building) %>%
  summarise(across(all_of(items), mean)) %>%
  as.data.frame()

# Calculate similarity for each resident
data$similarity <- NA

for (i in 1:nrow(data)) {
  building <- data$building[i]
  profile <- incumbent_profiles[incumbent_profiles$building == building, items]
  individual <- as.numeric(data[i, items])
  data$similarity[i] <- cosine_similarity(individual, as.numeric(profile))
}

# Summary by building and resident type
similarity_summary <- data %>%
  group_by(building, resident_type) %>%
  summarise(
    mean_similarity = mean(similarity),
    sd_similarity = sd(similarity),
    n = n(),
    .groups = "drop"
  )

cat("Similarity by building and resident type:\n")
print(similarity_summary)

# Wilcoxon test for new vs incumbent within each building
cat("\nWilcoxon tests (New vs Incumbent):\n")
for (b in unique(data$building)) {
  new_sim <- data$similarity[data$building == b & data$resident_type == "New"]
  inc_sim <- data$similarity[data$building == b & data$resident_type == "Incumbent"]
  
  if (length(new_sim) > 0 & length(inc_sim) > 0) {
    test <- wilcox.test(new_sim, inc_sim)
    cat(b, ": W =", test$statistic, ", p =", round(test$p.value, 3), "\n")
  }
}

# -----------------------------------------------------------------------------
# 5. Cross-Gender Cultural Influence Analysis (RQ4)
# -----------------------------------------------------------------------------

cat("\n=== Cross-Gender Cultural Influence Analysis ===\n")

# Filter co-ed buildings and new residents
coed_new <- data %>%
  filter(building %in% c("Rosemary", "Paprika"), resident_type == "New")

if (nrow(coed_new) > 0) {
  # Calculate same-sex and opposite-sex incumbent profiles
  coed_incumbents <- data %>%
    filter(building %in% c("Rosemary", "Paprika"), resident_type == "Incumbent")
  
  # Calculate profiles by building and sex
  profiles_by_sex <- coed_incumbents %>%
    group_by(building, sex) %>%
    summarise(across(all_of(items), mean), .groups = "drop")
  
  # Calculate similarity to same-sex and opposite-sex seniors
  coed_new$sim_same_sex <- NA
  coed_new$sim_opposite_sex <- NA
  
  for (i in 1:nrow(coed_new)) {
    building <- coed_new$building[i]
    sex <- coed_new$sex[i]
    opposite_sex <- ifelse(sex == "Male", "Female", "Male")
    
    same_profile <- profiles_by_sex %>%
      filter(building == !!building, sex == !!sex)
    
    opposite_profile <- profiles_by_sex %>%
      filter(building == !!building, sex == !!opposite_sex)
    
    individual <- as.numeric(coed_new[i, items])
    
    if (nrow(same_profile) > 0) {
      coed_new$sim_same_sex[i] <- cosine_similarity(individual, as.numeric(same_profile[, items]))
    }
    if (nrow(opposite_profile) > 0) {
      coed_new$sim_opposite_sex[i] <- cosine_similarity(individual, as.numeric(opposite_profile[, items]))
    }
  }
  
  # Paired t-test
  valid_data <- coed_new %>% filter(!is.na(sim_same_sex) & !is.na(sim_opposite_sex))
  
  if (nrow(valid_data) > 1) {
    t_test_result <- t.test(valid_data$sim_same_sex, valid_data$sim_opposite_sex, paired = TRUE)
    
    cat("Same-sex similarity: M =", round(mean(valid_data$sim_same_sex), 1), 
        "%, SD =", round(sd(valid_data$sim_same_sex), 1), "%\n")
    cat("Opposite-sex similarity: M =", round(mean(valid_data$sim_opposite_sex), 1), 
        "%, SD =", round(sd(valid_data$sim_opposite_sex), 1), "%\n")
    cat("Paired t-test: t =", round(t_test_result$statistic, 2), 
        ", df =", t_test_result$parameter,
        ", p =", round(t_test_result$p.value, 3), "\n")
  }
}

# -----------------------------------------------------------------------------
# 6. Building Type Comparison (RQ5)
# -----------------------------------------------------------------------------

cat("\n=== Building Type Comparison ===\n")

# One-way ANOVA for each item
anova_results <- data.frame(
  item = items,
  F_stat = NA,
  p_value = NA
)

for (i in seq_along(items)) {
  formula <- as.formula(paste(items[i], "~ building_type"))
  aov_result <- aov(formula, data = data)
  aov_summary <- summary(aov_result)
  
  anova_results$F_stat[i] <- aov_summary[[1]]$`F value`[1]
  anova_results$p_value[i] <- aov_summary[[1]]$`Pr(>F)`[1]
}

cat("ANOVA results by building type:\n")
print(anova_results)

# Post-hoc comparisons for significant items
significant_items <- anova_results$item[anova_results$p_value < 0.05]

if (length(significant_items) > 0) {
  cat("\nPost-hoc Tukey HSD for significant items:\n")
  for (item in significant_items) {
    cat("\n", item, ":\n")
    formula <- as.formula(paste(item, "~ building_type"))
    aov_result <- aov(formula, data = data)
    tukey <- TukeyHSD(aov_result)
    print(tukey)
  }
}

# Comparison between Rosemary and Paprika
cat("\nComparison between Rosemary and Paprika:\n")
coed_data <- data %>% filter(building %in% c("Rosemary", "Paprika"))

for (item in items) {
  test <- t.test(as.formula(paste(item, "~ building")), data = coed_data)
  if (test$p.value < 0.05) {
    cat(item, ": t =", round(test$statistic, 2), 
        ", p =", round(test$p.value, 3), "*\n")
  }
}

# -----------------------------------------------------------------------------
# 7. Visualization
# -----------------------------------------------------------------------------

cat("\n=== Generating Figures ===\n")

# Color palette
building_colors <- c(
  "Basil" = "#1f77b4",
  "Turmeric" = "#ff7f0e",
  "Rosemary" = "#2ca02c",
  "Paprika" = "#d62728"
)

building_type_colors <- c(
  "Male-only" = "#1f77b4",
  "Female-only" = "#ff7f0e",
  "Co-ed" = "#2ca02c"
)

# Figure 2: Cultural Effects
# (a) ICC by item
p2a <- ggplot(icc_results, aes(x = reorder(item, icc), y = icc)) +
  geom_bar(stat = "identity", fill = "#4292c6") +
  geom_hline(yintercept = overall_icc$icc, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(x = "Item", y = "ICC", title = "(a) ICC by Item") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

# (b) Bootstrap distribution
boot_df <- data.frame(icc = boot_results$t)
p2b <- ggplot(boot_df, aes(x = icc)) +
  geom_histogram(bins = 30, fill = "#4292c6", color = "white", alpha = 0.7) +
  geom_vline(xintercept = overall_icc$icc, linetype = "solid", color = "black") +
  geom_vline(xintercept = boot_ci$percent[4:5], linetype = "dashed", color = "red") +
  labs(x = "ICC", y = "Frequency", title = "(b) Bootstrap Distribution") +
  theme_minimal()

# (c) PCoA plot
pcoa <- cmdscale(dist_aitchison, k = 2)
pcoa_df <- data.frame(
  PC1 = pcoa[, 1],
  PC2 = pcoa[, 2],
  building = data$building,
  resident_type = data$resident_type
)

p2c <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = building, shape = resident_type)) +
  geom_point(alpha = 0.7, size = 2) +
  stat_ellipse(aes(group = building), level = 0.95) +
  scale_color_manual(values = building_colors) +
  labs(x = "PC1", y = "PC2", title = "(c) PCoA Plot",
       color = "Building", shape = "Resident Type") +
  theme_minimal() +
  theme(legend.position = "right")

# Combine Figure 2
fig2 <- (p2a | p2b) / p2c + plot_layout(heights = c(1, 1.5))
ggsave("figures/fig2_cultural_effects.png", fig2, width = 12, height = 10, dpi = 300)

# Figure 3: Convergence Analysis
# (a) Similarity distribution by building
p3a <- ggplot(data, aes(x = building, y = similarity, fill = building)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = building_colors) +
  labs(x = "Building", y = "Similarity (%)", title = "(a) Similarity Distribution") +
  theme_minimal() +
  theme(legend.position = "none")

# (b) New vs Incumbent comparison
p3b <- ggplot(data, aes(x = building, y = similarity, fill = resident_type)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Building", y = "Similarity (%)", title = "(b) New vs Incumbent",
       fill = "Resident Type") +
  theme_minimal()

# Combine Figure 3
fig3 <- p3a + p3b + plot_layout(widths = c(1, 1.5))
ggsave("figures/fig3_convergence.png", fig3, width = 12, height = 5, dpi = 300)

# Figure 5: Building Type Comparison
# Mean scores by building type
mean_scores <- data %>%
  group_by(building_type) %>%
  summarise(across(all_of(items), list(mean = mean, se = ~sd(.)/sqrt(n())))) %>%
  pivot_longer(-building_type, names_to = c("item", "stat"), names_sep = "_") %>%
  pivot_wider(names_from = stat, values_from = value)

p5b <- ggplot(mean_scores, aes(x = item, y = mean, fill = building_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.9), width = 0.2) +
  scale_fill_manual(values = building_type_colors) +
  labs(x = "Item", y = "Mean Score", fill = "Building Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/fig5_building_type_comparison.png", p5b, width = 10, height = 6, dpi = 300)

cat("Figures saved to 'figures/' directory\n")

# -----------------------------------------------------------------------------
# 8. Export Results
# -----------------------------------------------------------------------------

cat("\n=== Exporting Results ===\n")

# Save ICC results
write_csv(icc_results, "results/icc_results.csv")

# Save ANOVA results
write_csv(anova_results, "results/anova_results.csv")

# Save similarity summary
write_csv(similarity_summary, "results/similarity_summary.csv")

# Save session info
sink("results/session_info.txt")
sessionInfo()
sink()

cat("Results saved to 'results/' directory\n")
cat("\n=== Analysis Complete ===\n")
