# Fig. 3


# 0. Libraries
library(tidyverse)
library(tidymodels)
library(purrr)
library(ggplot2)
library(kernelshap)
library(shapviz)
library(ranger)


# 1. Output directory
dir.create(
  "output/Fig3",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Load public MetaCheeseDB data
metadata <- read.csv(
  "data/MetaCheeseDB_metadata.tsv",
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

taxonomy <- read.csv(
  "data/MetaCheeseDB_taxonomy.tsv",
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

stopifnot(
  identical(
    rownames(metadata),
    rownames(taxonomy)
  )
)

metadata <- metadata %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(across(everything(), as.character))


# 3. Prepare species profiles
s_Cheeses <- taxonomy %>%
  rownames_to_column(var = "SampleID") %>%
  mutate(
    SampleID = as.character(SampleID),
    across(-SampleID, as.numeric)
  ) %>%
  mutate(across(-SampleID, ~ replace_na(.x, 0))) %>%
  select(where(~ !is.numeric(.) || sum(., na.rm = TRUE) > 0))


s_metaCheeses <- inner_join(metadata, s_Cheeses, by = "SampleID")
rownames(s_metaCheeses) <- s_metaCheeses$SampleID
s_metaCheeses$SampleID <- NULL


# 4. Mean-abundance filtering (mean > 0.08) and log1p transformation
species_mean <- s_metaCheeses %>%
  .[, !sapply(., is.character)] %>%
  colMeans() %>%
  as.data.frame()

species_prevalence <- species_mean %>%
  filter_all(all_vars(. > 0.08))

selected_features <- rownames(species_prevalence)

species_subset <- s_metaCheeses %>%
  select(all_of(selected_features))

species_subset_log1p <- log1p(species_subset[, !sapply(species_subset, is.character)])

s_metaCheeses$samples <- rownames(s_metaCheeses)
species_subset_log1p$samples <- rownames(species_subset_log1p)

s_cheeses_meta_heatmap_log1p <- semi_join(
  x = species_subset_log1p,
  y = s_metaCheeses,
  by = "samples"
)
rownames(s_cheeses_meta_heatmap_log1p) <- s_cheeses_meta_heatmap_log1p$samples
s_cheeses_meta_heatmap_log1p$samples <- NULL
s_cheeses_meta_heatmap_log1p <- as.data.frame(s_cheeses_meta_heatmap_log1p)


# 5. Build the dataset (target is Rheological_properties)
s_cheeses_meta_heatmap_log1p$SampleID <- rownames(s_cheeses_meta_heatmap_log1p)

s_target_log1p <-
  s_cheeses_meta_heatmap_log1p %>%
  right_join(
    metadata %>% select(SampleID, Rheological_properties),
    by = "SampleID"
  ) %>%
  filter(!is.na(Rheological_properties)) %>%
  relocate(SampleID, Rheological_properties, .before = 1)

rownames(s_target_log1p) <- s_target_log1p$SampleID
s_target_log1p$SampleID <- NULL

# Convert target to factor
s_target_log1p$Rheological_properties <-
  as.factor(s_target_log1p$Rheological_properties)


# 6. Train/test split and Random Forest model
set.seed(123)

data_split <- initial_split(
  s_target_log1p,
  prop = 0.7,
  strata = Rheological_properties
)
train_data <- training(data_split)
test_data  <- testing(data_split)

rf_model <- rand_forest(trees = 100) %>%
  set_engine("randomForest") %>%
  set_mode("classification")

rf_workflow <- workflow() %>%
  add_model(rf_model) %>%
  add_formula(Rheological_properties ~ .)

rf_fit <- fit(rf_workflow, data = train_data)

# Predict class probabilities on test set
rf_predictions <- predict(rf_fit, test_data, type = "prob") %>%
  bind_cols(test_data %>% select(Rheological_properties))


# 7. One-vs-Rest ROC curve (multiclass) with bootstrap CI
roc_data_ovr <- rf_predictions %>%
  roc_curve(
    Rheological_properties,
    .pred_Hard,
    .pred_Semihard,
    .pred_Soft,
    .pred_Variable,
    .pred_Very_hard
  )

auc_data_ovr <- rf_predictions %>%
  roc_auc(
    Rheological_properties,
    .pred_Hard,
    .pred_Semihard,
    .pred_Soft,
    .pred_Variable,
    .pred_Very_hard
  )
auc_value_ovr <- auc_data_ovr$.estimate

# Bootstrap AUC (1000 replicates)
set.seed(123)
boot_auc_ovr <- rf_predictions %>%
  bootstraps(times = 1000) %>%
  mutate(
    auc = map(
      splits,
      ~ analysis(.x) %>%
        roc_auc(
          Rheological_properties,
          .pred_Hard,
          .pred_Semihard,
          .pred_Soft,
          .pred_Variable,
          .pred_Very_hard
        )
    )
  ) %>%
  unnest(auc)

ci_auc_ovr <- quantile(boot_auc_ovr$.estimate, probs = c(0.025, 0.975))

# Plot and save – now using "Texture" in title and filename
roc_plot_ovr <- autoplot(roc_data_ovr) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle(
    paste0(
      "Texture, One-vs-Rest ROC curves\n",
      "Bootstrap AUC = ", round(auc_value_ovr, 2),
      "\n95% CI: [", round(ci_auc_ovr[1], 2), "-", round(ci_auc_ovr[2], 2), "]"
    )
  )

roc_plot_ovr$layers <- roc_plot_ovr$layers[
  !vapply(
    roc_plot_ovr$layers,
    function(layer) inherits(layer$geom, "GeomAbline"),
    logical(1)
  )
]

ggsave(
  filename = "output/Fig3/Fig3_OvR_Texture.svg",
  plot     = roc_plot_ovr,
  device   = "svg",
  width    = 6,
  height   = 6
)


# 8. One-vs-One pairwise ROC curves
classes <- levels(test_data$Rheological_properties)
class_pairs <- combn(classes, 2, simplify = FALSE)

pairwise_results <- map(class_pairs, ~ {
  
  # Filter to current pair
  pair_data <- rf_predictions %>%
    filter(Rheological_properties %in% .x) %>%
    mutate(Rheological_properties = factor(Rheological_properties))
  
  class1 <- .x[1]
  class2 <- .x[2]
  prob_col <- paste0(".pred_", class2)
  
  # Compute AUC
  auc_value <- pair_data %>%
    roc_auc(Rheological_properties, !!sym(prob_col), event_level = "second")
  
  # Bootstrap AUC (1000 replicates)
  set.seed(123)
  boot_auc <- pair_data %>%
    bootstraps(times = 1000) %>%
    mutate(
      auc = map(
        splits,
        ~ analysis(.x) %>%
          roc_auc(Rheological_properties, !!sym(prob_col), event_level = "second")
      )
    ) %>%
    unnest(auc)
  
  ci_auc <- quantile(boot_auc$.estimate, probs = c(0.025, 0.975))
  
  # Compute ROC curve
  roc_data <- pair_data %>%
    roc_curve(Rheological_properties, !!sym(prob_col), event_level = "second")
  
  list(
    pair = paste(class1, "vs", class2),
    auc = auc_value$.estimate,
    ci_lower = ci_auc[1],
    ci_upper = ci_auc[2],
    roc_data = roc_data
  )
})

walk(pairwise_results, ~ {
  
  p <- autoplot(.x$roc_data) +
    ggtitle(
      paste0(
        "Texture, One-vs-One ROC: ", .x$pair,
        "\nAUC = ", round(.x$auc, 2),
        " 95% CI [", round(.x$ci_lower, 2), "-", round(.x$ci_upper, 2), "]"
      )
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p$layers <- p$layers[
    !vapply(
      p$layers,
      function(layer) inherits(layer$geom, "GeomAbline"),
      logical(1)
    )
  ]
  
  # Output filename
  filename <- paste0(
    "output/Fig3/ROC_Texture_",
    gsub(" ", "_", .x$pair),
    ".svg"
  )
  
  # Save plot
  ggsave(
    filename = filename,
    plot = p,
    device = "svg",
    width = 6,
    height = 6
  )
  
})


# 9. SHAP analysis
s_target_log1p$Rheological_properties <-
  as.factor(s_target_log1p$Rheological_properties)

# 9a. Train probability ranger model
fit_Rheological_properties <- ranger(
  Rheological_properties ~ .,
  data = s_target_log1p,
  probability = TRUE,
  num.trees = 1000
)

# 9b. Compute kernel SHAP explanations
ks_Rheological_properties <- kernelshap(
  fit_Rheological_properties,
  X = s_target_log1p[, -1],
  bg_X = s_target_log1p[, -1]
)


# 9c. Wrap into shapviz
sv_Rheological_properties <- shapviz(
  ks_Rheological_properties
)


# 9d. Standard feature importance
p_shap_Rheological_properties <- sv_importance(
  sv_Rheological_properties
)

ggsave(
  filename = "output/Fig3/Fig3_SHAP_Texture.svg",
  plot = p_shap_Rheological_properties,
  width = 10,
  height = 5,
  dpi = 600
)


# 9e. Beeswarm feature importance
p_shap_Rheological_properties_bee <- sv_importance(
  sv_Rheological_properties,
  kind = "bee"
)

ggsave(
  filename = "output/Fig3/Fig3_SHAP_Texture_bee.svg",
  plot = p_shap_Rheological_properties_bee,
  width = 25,
  height = 5,
  dpi = 600
)
