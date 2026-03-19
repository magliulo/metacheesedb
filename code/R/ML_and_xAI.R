# ============================================================================ #
#              Cheese Microbiome Catalogue — ML & SHAP Analysis                #
#                               Raffaele Magliulo                              #
#   Department of Agricultural Sciences, University of Naples Federico II      #
#   Department of Food Biosciences, Tegasc                                     #

# ============================================================================ #
#
# This script contains the supervised machine learning and SHAP-based
# explainability analyses described in:
#
# "Creation of a catalogue containing 1,593 metagenomes unearths the 
#  microbiome role in the terroir of cheese"
#
# Part_of_cheese (core vs rind classification) is provided here as a
# representative example. The same approach was applied to all 18 metadata
# categories listed in the manuscript. Full results are reported in
# Supplementary File 7.
#
# Input:
# - s_Part_of_cheese.csv  (provided in this repository)
#
# Output:
# - AUC-ROC value and 95% confidence interval for Part_of_cheese classification
# - SHAP beeswarm plot for Part_of_cheese classification
#
# Repository: https://github.com/magliulo/metacheesedb
#
# ============================================================================ #



#=== Setup environment =========================================================

setwd("/path/to/folder")
dir.create("output")

# Load required libraries:

## Data manipulation
library(tidyverse)

## Machine learning and explainable artificial intelligence
library(tidymodels)
library(ranger)
library(randomForest)
library(kernelshap)
library(shapviz)



#=== 1. Data ===================================================================

# Import pre-processed input data frame (species log1p abundances + Part_of_cheese label)
s_Part_of_cheese <- read.csv(
  "s_Part_of_cheese.csv",
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

cat("Classes:", levels(as.factor(s_Part_of_cheese$Part_of_cheese)), "\n")



#=== 2. sML ====================================================================

# The same workflow was applied to all 18 metadata categories.
# AUC-ROC values and 95% CIs for all categories are reported in Supplementary File 7.

### ===> Part_of_cheese | One-vs-rest
# Convert the outcome variable to a factor
s_Part_of_cheese$Part_of_cheese <- as.factor(s_Part_of_cheese$Part_of_cheese)

# Set seed for reproducibility
set.seed(123)

# Split the data into training and testing sets
data_split <- initial_split(s_Part_of_cheese, prop = 0.7, strata = Part_of_cheese)
train_data <- training(data_split)
test_data  <- testing(data_split)

# Specify the Random Forest model
rf_model <- rand_forest(trees = 100) %>%
  set_engine("randomForest") %>%
  set_mode("classification")

# Create a workflow
rf_workflow <- workflow() %>%
  add_model(rf_model) %>%
  add_formula(Part_of_cheese ~ .)

# Fit the model
rf_fit <- fit(rf_workflow, data = train_data)

# Predict class probabilities
rf_predictions <- predict(rf_fit, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(Part_of_cheese))

# Compute ROC curves
roc_data <- rf_predictions %>%
  roc_curve(Part_of_cheese, .pred_Core, .pred_Core_AND_rind, .pred_Rind)

# Compute AUC
auc_data <- rf_predictions %>%
  roc_auc(Part_of_cheese, .pred_Core, .pred_Core_AND_rind, .pred_Rind)

# Compute original AUC
auc_value <- auc_data$.estimate

# Bootstrap AUC with CI
set.seed(123)
boot_auc <- rf_predictions %>%
  bootstraps(times = 1000) %>%
  mutate(auc = map(splits, ~ analysis(.x) %>% 
                     roc_auc(Part_of_cheese, .pred_Core, .pred_Core_AND_rind, .pred_Rind))) %>%
  unnest(auc)

# Calculate 95% CI
ci_auc <- quantile(boot_auc$.estimate, probs = c(0.025, 0.975))

roc_data <- rf_predictions %>%
  roc_curve(Part_of_cheese, .pred_Core:.pred_Rind)

roc_plot <- autoplot(roc_data) + ggtitle("Tidymodels One-vs-All ROC") +
  ggtitle(paste0("Part of cheese, One-vs-Rest ROC curves\nBootsrapping AUC = ", round(auc_value, 2),
                 "\n95% CI: [", round(ci_auc[1], 2), "-", round(ci_auc[2], 2), "]"))
roc_plot
ggsave("output/9_OvR_Part_of_cheese.svg", plot = roc_plot, device = "svg", width = 6, height = 6)



### ===> Part_of_cheese | One-vs-one
classes <- levels(test_data$Part_of_cheese)
class_pairs <- combn(classes, 2, simplify = FALSE)  # All unique pairs
# Compute one-vs-one AUC and CIs for all pairs
pairwise_results <- map(class_pairs, ~ {
  # Filter data to the current pair
  pair_data <- rf_predictions %>%
    filter(Part_of_cheese %in% .x) %>%
    mutate(Part_of_cheese = factor(Part_of_cheese))  # Drop unused levels
  
  # Extract class names and probability column
  class1 <- .x[1]
  class2 <- .x[2]
  prob_col <- paste0(".pred_", class2)
  
  # Compute AUC
  auc_value <- pair_data %>%
    roc_auc(Part_of_cheese, !!sym(prob_col), event_level = "second")
  
  # Bootstrap AUC
  set.seed(123)
  boot_auc <- pair_data %>%
    bootstraps(times = 1000) %>%
    mutate(auc = map(splits, ~ analysis(.x) %>% 
                       roc_auc(Part_of_cheese, !!sym(prob_col), event_level = "second"))) %>%
    unnest(auc)
  
  # Calculate 95% CI
  ci_auc <- quantile(boot_auc$.estimate, probs = c(0.025, 0.975))
  
  # Compute ROC curve
  roc_data <- pair_data %>%
    roc_curve(Part_of_cheese, !!sym(prob_col), event_level = "second")
  
  # Return results
  list(
    pair = paste(class1, "vs", class2),
    auc = auc_value$.estimate,
    ci_lower = ci_auc[1],
    ci_upper = ci_auc[2],
    roc_data = roc_data
  )
})
# Plot all one-vs-one ROC curves
walk(pairwise_results, ~ {
  p <- autoplot(.x$roc_data) +
    ggtitle(paste0("Part_of_cheese, One-vs-One ROC: ", .x$pair,
                   "\nAUC = ", round(.x$auc, 2),
                   " 95% CI [", round(.x$ci_lower, 2), "-", round(.x$ci_upper, 2), "]"))

  # Define filename (safe for filesystem)
  filename <- paste0("output/9_OvO_Part_of_cheese_", gsub(" ", "_", .x$pair), ".svg")
  
  # Save the plot
  ggsave(filename, plot = p, device = "svg", width = 6, height = 6)
})



#=== 3. xAI ====================================================================

# The same SHAP workflow was applied to all 18 metadata categories.

# 1a) Make sure the “Part_of_cheese” column is a factor
s_Part_of_cheese$Part_of_cheese <- as.factor(s_Part_of_cheese$Part_of_cheese)

# 1b) Train a probability ranger model:
set.seed(123)
fit_Part_of_cheese <- ranger(
  Part_of_cheese ~ .,
  data        = s_Part_of_cheese,
  probability = TRUE,
  num.trees   = 500
)

# 1c) Compute kernel SHAP explanations:
ks_Part_of_cheese <- kernelshap(
  fit_Part_of_cheese,
  X    = s_Part_of_cheese[,-1],  # all predictors except the target column
  bg_X = s_Part_of_cheese[,-1]
)

# 1d) Wrap into shapviz:
sv_Part_of_cheese <- shapviz(ks_Part_of_cheese)

# 1e) Plot (standard feature‐importance) and save:
sv_importance(sv_Part_of_cheese)
ggsave(
  "output/shap_Part_of_cheese_bar.svg",
  width  = 10,
  height = 5,
  dpi    = 600
)

# 1f) Plot (“bee”‐style) and save:
sv_importance(sv_Part_of_cheese, kind = "bee")
ggsave(
  "output/shap_Part_of_cheese_bee.svg",
  width  = 25,
  height = 5,
  dpi    = 600
)
