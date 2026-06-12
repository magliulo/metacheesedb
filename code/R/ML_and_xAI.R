#!/usr/bin/env Rscript

# MetaCheeseDB revised machine-learning and xAI analysis
# -------------------------------------------------------
# Purpose:
#   Reproduce the reviewer-sensitive Random Forest analyses used to evaluate
#   country-associated microbiome patterns in core cheese samples.
#
# Main analyses:
#   1. Italy vs Spain classification under four settings:
#      - optimistic sample-level cross-validation;
#      - leakage-aware dataset-grouped cross-validation;
#      - metadata-only confounder model;
#      - confounder-adjusted microbiome model.
#   2. One-vs-one country comparisons for Italy, Spain, Ireland and Austria.
#   3. SHAP interpretation of the Italy vs Spain microbiome model.
#
# Expected input files:
#   data/MetaCheeseDB_metadata.tsv
#   data/MetaCheeseDB_taxonomy.tsv
#
# Usage from the repository root:
#   Rscript code/R/ML_and_xAI.R
#
# Optional arguments:
#   --metadata <path>   path to MetaCheeseDB_metadata.tsv
#   --taxonomy <path>   path to MetaCheeseDB_taxonomy.tsv
#   --outdir <path>     output directory
#   --run-shap <TRUE/FALSE> whether to run kernel SHAP, which can be slow

parse_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  pos <- match(flag, args)
  if (!is.na(pos) && length(args) >= pos + 1) {
    return(args[pos + 1])
  }
  default
}

metadata_file <- parse_arg("--metadata", file.path("data", "MetaCheeseDB_metadata.tsv"))
taxonomy_file <- parse_arg("--taxonomy", file.path("data", "MetaCheeseDB_taxonomy.tsv"))
out_dir       <- parse_arg("--outdir", file.path("results", "ML_reviewer_sensitivity"))
run_shap      <- tolower(parse_arg("--run-shap", "TRUE")) %in% c("true", "t", "1", "yes", "y")

required_packages <- c(
  "tidyverse", "svglite", "caret", "pROC", "vegan", "ranger",
  "writexl", "ggsci", "scico", "kernelshap", "shapviz"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(svglite)
  library(caret)
  library(pROC)
  library(vegan)
  library(ranger)
  library(writexl)
  library(ggsci)
  library(scico)
  library(kernelshap)
  library(shapviz)
  library(grid)
})

set.seed(123)

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file, call. = FALSE)
}
if (!file.exists(taxonomy_file)) {
  stop("Taxonomy file not found: ", taxonomy_file, call. = FALSE)
}

# 1. Load and align data


metadata <- read.csv(
  metadata_file,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "NaN", "Unknown", "unknown", "not available")
)

taxa_table <- read.csv(
  taxonomy_file,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# Keep only common sample IDs and force identical ordering
common_ids <- intersect(rownames(metadata), rownames(taxa_table))
metadata <- metadata[common_ids, , drop = FALSE]
taxa_table <- taxa_table[common_ids, , drop = FALSE]
taxa_table <- taxa_table[rownames(metadata), , drop = FALSE]

stopifnot(identical(rownames(metadata), rownames(taxa_table)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# 2. Define the exact ML task


# Reviewer example: Italy vs Spain, Core only.
# Change these only if you want another comparison.
target_var     <- "Country"
positive_class <- "Italy"
negative_class <- "Spain"

meta0 <- metadata %>%
  filter(
    Part_of_cheese == "Core",
    Country %in% c(positive_class, negative_class)
  ) %>%
  mutate(
    y = factor(
      if_else(.data[[target_var]] == positive_class, positive_class, negative_class),
      levels = c(positive_class, negative_class)
    )
  )

taxa0 <- taxa_table[rownames(meta0), , drop = FALSE]

cat("\nSamples used in Italy vs Spain Core analysis:\n")
print(table(meta0$y))


# 3. Diagnose possible technical/biological replicate structure


# Candidate grouping variables.
# Inspect this table and set preferred_group_var below to the best available
# biological-sample / cheese / batch / study identifier.
candidate_check <- c(
  "Dataset",
  "Accession",
  "Run",
  "Experiment",
  "Study",
  "Project"
)

for (v in candidate_check) {
  
  cat("\n====================\n")
  cat("VARIABLE:", v, "\n")
  
  print(table(meta0[[v]], useNA = "ifany"))
  
  cat("\nUnique values:", n_distinct(meta0[[v]], na.rm = TRUE), "\n")
  
  cat("\nRepeated groups:\n")
  
  print(
    sort(table(meta0[[v]]), decreasing = TRUE)[1:10]
  )
}

# IMPORTANT:
# Set this to the best variable identifying non-independent samples.
# Examples: "Cheese_ID", "Sample_ID", "Study_ID", "Source_DOI", "Batch", etc.
# If no such metadata column exists, keep NA; the script will use each sample as
# its own group, but that is less conservative.
preferred_group_var <- "Dataset"

meta0$cv_group <- interaction(
  meta0$Dataset,
  meta0$Country,
  drop = TRUE
)

if (!is.na(preferred_group_var) && preferred_group_var %in% names(meta0)) {
  meta0$cv_group <- as.character(meta0[[preferred_group_var]])
  meta0$cv_group[is.na(meta0$cv_group) | meta0$cv_group == ""] <- rownames(meta0)[
    is.na(meta0$cv_group) | meta0$cv_group == ""
  ]
} else {
  meta0$cv_group <- rownames(meta0)
}

cat("\nCross-validation grouping variable used:\n")
cat(ifelse(is.na(preferred_group_var), "Each sample used as its own group", preferred_group_var), "\n")
cat("Number of CV groups:", n_distinct(meta0$cv_group), "\n")


# 4. Near-duplicate taxonomic profiles


to_numeric_matrix <- function(x) {
  x <- as.data.frame(x)
  x[] <- lapply(x, function(z) as.numeric(as.character(z)))
  m <- as.matrix(x)
  m[is.na(m)] <- 0
  storage.mode(m) <- "numeric"
  m
}

relative_abundance <- function(x) {
  m <- to_numeric_matrix(x)
  rs <- rowSums(m, na.rm = TRUE)
  keep <- rs > 0
  m[keep, ] <- sweep(m[keep, , drop = FALSE], 1, rs[keep], "/")
  m[!is.finite(m)] <- 0
  m
}

taxa_rel <- relative_abundance(taxa0)

bray <- as.matrix(vegdist(taxa_rel, method = "bray"))
diag(bray) <- NA

nearest_idx <- apply(bray, 1, function(z) which.min(z))

nearest_report <- tibble(
  sample_id = rownames(taxa_rel),
  nearest_sample = rownames(taxa_rel)[nearest_idx],
  bray_curtis_distance = bray[cbind(seq_len(nrow(bray)), nearest_idx)],
  sample_class = as.character(meta0$y),
  nearest_class = as.character(meta0$y[nearest_idx]),
  cv_group = meta0$cv_group,
  nearest_cv_group = meta0$cv_group[nearest_idx],
  same_cv_group = cv_group == nearest_cv_group,
  same_class = sample_class == nearest_class
) %>%
  arrange(bray_curtis_distance)

cat("\nNearest-neighbour duplicate check, first 10 pairs:\n")
print(head(nearest_report, 10))

# Pairs with Bray-Curtis distance near zero but different CV groups are suspicious:
possible_unblocked_replicates <- nearest_report %>%
  filter(bray_curtis_distance <= 1e-6, !same_cv_group)


# 5. Confounder diagnostics: are metadata variables associated with Italy/Spain?

confounder_vars <- c(
  "Milk_source",
  "Animal_feeding",
  "Milk_processing",
  "Thermisation",
  "Pasteurization",
  "Skimming",
  "Inoculation_of_moulds",
  "Curd_cutting",
  "Presence_of_rind",
  "Backslopping",
  "Starter_culture",
  "Technology_of_production",
  "Coagulation_method",
  "Ripening_period",
  "Rheological_properties",
  "Temperature_of_curd_processing_Celsius_degree"
)

confounder_vars <- intersect(confounder_vars, names(meta0))

safe_assoc_test <- function(tab) {
  if (any(dim(tab) < 2)) return(NA_real_)
  suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 5000)$p.value)
}

cramers_v <- function(tab) {
  if (any(dim(tab) < 2)) return(NA_real_)
  suppressWarnings({
    chi <- chisq.test(tab, simulate.p.value = TRUE, B = 5000)$statistic
    n <- sum(tab)
    k <- min(dim(tab))
    as.numeric(sqrt(chi / (n * (k - 1))))
  })
}

confounder_screen <- map_dfr(confounder_vars, function(v) {
  z <- as.character(meta0[[v]])
  z[is.na(z) | z == "" | z == "NA"] <- NA
  
  tab <- table(meta0$y, z, useNA = "no")
  
  tibble(
    variable = v,
    n_complete = sum(!is.na(z)),
    n_levels = n_distinct(z, na.rm = TRUE),
    p_value = safe_assoc_test(tab),
    cramers_v = cramers_v(tab),
    levels_Italy = paste(head(names(sort(table(z[meta0$y == positive_class]), decreasing = TRUE)), 6), collapse = "; "),
    levels_Spain = paste(head(names(sort(table(z[meta0$y == negative_class]), decreasing = TRUE)), 6), collapse = "; ")
  )
}) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(q_value, desc(cramers_v))


# 6. Utility functions for leakage-aware Leakage-aware


make_group_folds <- function(y, group, k = 5, seed = 123) {
  set.seed(seed)
  
  y <- droplevels(y)
  group <- as.character(group)
  
  # Assign one class label per group.
  group_y <- tapply(as.character(y), group, function(z) {
    names(sort(table(z), decreasing = TRUE))[1]
  })
  
  group_y <- factor(group_y, levels = levels(y))
  class_counts <- table(group_y)
  
  k_eff <- min(k, length(group_y), min(class_counts))
  
  if (k_eff < 2) {
    stop("Not enough independent groups per class for Leakage-aware.")
  }
  
  fold_groups <- caret::createFolds(group_y, k = k_eff, list = TRUE, returnTrain = FALSE)
  
  lapply(fold_groups, function(gidx) {
    which(group %in% names(group_y)[gidx])
  })
}

prep_taxa_train <- function(x_train, min_prevalence = 0.05) {
  x_rel <- relative_abundance(x_train)
  
  prevalence <- colMeans(x_rel > 0)
  keep <- prevalence >= min_prevalence & is.finite(prevalence)
  
  # Fallback if prevalence filter is too strict
  if (sum(keep) < 2) {
    keep <- rank(-prevalence, ties.method = "first") <= min(100, ncol(x_rel))
  }
  
  x_rel <- x_rel[, keep, drop = FALSE]
  
  pseudo <- suppressWarnings(min(x_rel[x_rel > 0], na.rm = TRUE) / 2)
  if (!is.finite(pseudo) || pseudo <= 0) pseudo <- 1e-6
  
  x_clr <- log(x_rel + pseudo)
  x_clr <- sweep(x_clr, 1, rowMeans(x_clr), "-")
  
  original_features <- colnames(x_clr)
  safe_features <- make.names(original_features, unique = TRUE)
  colnames(x_clr) <- safe_features
  
  list(
    features_original = original_features,
    features_safe = safe_features,
    pseudo = pseudo,
    train = as.data.frame(x_clr)
  )
}

prep_taxa_test <- function(x_test, prep) {
  x_rel <- relative_abundance(x_test)
  
  missing_features <- setdiff(prep$features_original, colnames(x_rel))
  if (length(missing_features) > 0) {
    stop("Some training features are absent from the test table.")
  }
  
  x_rel <- x_rel[, prep$features_original, drop = FALSE]
  
  x_clr <- log(x_rel + prep$pseudo)
  x_clr <- sweep(x_clr, 1, rowMeans(x_clr), "-")
  colnames(x_clr) <- prep$features_safe
  
  as.data.frame(x_clr)
}

prep_covariate_mm <- function(meta_train, meta_test, covariates) {
  covariates <- intersect(covariates, names(meta_train))
  
  if (length(covariates) == 0) {
    return(list(
      train = matrix(nrow = nrow(meta_train), ncol = 0),
      test = matrix(nrow = nrow(meta_test), ncol = 0),
      covariates_used = character()
    ))
  }
  
  tr <- meta_train[, covariates, drop = FALSE]
  te <- meta_test[, covariates, drop = FALSE]
  
  use <- character()
  
  for (v in covariates) {
    if (is.numeric(tr[[v]])) {
      med <- median(tr[[v]], na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      tr[[v]][is.na(tr[[v]])] <- med
      te[[v]][is.na(te[[v]])] <- med
      
      if (sd(tr[[v]], na.rm = TRUE) > 0) use <- c(use, v)
      
    } else {
      tr[[v]] <- as.character(tr[[v]])
      te[[v]] <- as.character(te[[v]])
      
      tr[[v]][is.na(tr[[v]]) | tr[[v]] == "" | tr[[v]] == "NA"] <- "Missing"
      te[[v]][is.na(te[[v]]) | te[[v]] == "" | te[[v]] == "NA"] <- "Missing"
      
      train_levels <- sort(unique(tr[[v]]))
      te[[v]][!te[[v]] %in% train_levels] <- "Other"
      
      all_levels <- unique(c(train_levels, "Other"))
      tr[[v]] <- factor(tr[[v]], levels = all_levels)
      te[[v]] <- factor(te[[v]], levels = all_levels)
      
      if (n_distinct(tr[[v]]) > 1) use <- c(use, v)
    }
  }
  
  if (length(use) == 0) {
    return(list(
      train = matrix(nrow = nrow(meta_train), ncol = 0),
      test = matrix(nrow = nrow(meta_test), ncol = 0),
      covariates_used = character()
    ))
  }
  
  form <- as.formula(paste("~", paste(sprintf("`%s`", use), collapse = " + ")))
  
  mm_train <- model.matrix(form, data = tr)[, -1, drop = FALSE]
  mm_test  <- model.matrix(form, data = te)[, -1, drop = FALSE]
  
  # Align columns defensively
  all_cols <- union(colnames(mm_train), colnames(mm_test))
  
  add_missing_cols <- function(m, cols) {
    missing <- setdiff(cols, colnames(m))
    if (length(missing) > 0) {
      add <- matrix(0, nrow = nrow(m), ncol = length(missing))
      colnames(add) <- missing
      m <- cbind(m, add)
    }
    m[, cols, drop = FALSE]
  }
  
  mm_train <- add_missing_cols(mm_train, all_cols)
  mm_test  <- add_missing_cols(mm_test, all_cols)
  
  list(
    train = mm_train,
    test = mm_test,
    covariates_used = use
  )
}

residualize_train_test <- function(x_train, x_test, meta_train, meta_test, covariates) {
  mm <- prep_covariate_mm(meta_train, meta_test, covariates)
  
  if (ncol(mm$train) == 0) {
    return(list(train = x_train, test = x_test, covariates_used = character()))
  }
  
  Xtr <- cbind(Intercept = 1, mm$train)
  Xte <- cbind(Intercept = 1, mm$test)
  
  Ytr <- as.matrix(x_train)
  Yte <- as.matrix(x_test)
  
  # Ridge-stabilised least squares for high-dimensional features
  lambda <- 1e-6
  beta <- solve(crossprod(Xtr) + diag(lambda, ncol(Xtr)), crossprod(Xtr, Ytr))
  
  res_train <- Ytr - Xtr %*% beta
  res_test  <- Yte - Xte %*% beta
  
  colnames(res_train) <- colnames(x_train)
  colnames(res_test)  <- colnames(x_test)
  
  list(
    train = as.data.frame(res_train),
    test = as.data.frame(res_test),
    covariates_used = mm$covariates_used
  )
}

summarise_binary <- function(preds, positive = positive_class, negative = negative_class) {
  preds$obs <- factor(preds$obs, levels = c(positive, negative))
  
  roc_obj <- pROC::roc(
    response = preds$obs,
    predictor = preds[[positive]],
    levels = c(negative, positive),
    quiet = TRUE
  )
  
  pred_class <- factor(
    if_else(preds[[positive]] >= 0.5, positive, negative),
    levels = c(positive, negative)
  )
  
  cm <- confusionMatrix(pred_class, preds$obs, positive = positive)
  
  tibble(
    n = nrow(preds),
    n_folds = n_distinct(preds$fold),
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUC_low = as.numeric(pROC::ci.auc(roc_obj)[1]),
    AUC_high = as.numeric(pROC::ci.auc(roc_obj)[3]),
    Accuracy = unname(cm$overall["Accuracy"]),
    Sensitivity = unname(cm$byClass["Sensitivity"]),
    Specificity = unname(cm$byClass["Specificity"]),
    BalancedAccuracy = unname(cm$byClass["Balanced Accuracy"])
  )
}

run_nested_grouped_rf <- function(taxa, meta,
                                  outcome = "y",
                                  group = "cv_group",
                                  confounders = NULL,
                                  residualize = FALSE,
                                  k_outer = 5,
                                  k_inner = 3,
                                  seed = 123) {
  
  y <- droplevels(meta[[outcome]])
  outer_test_folds <- make_group_folds(y, meta[[group]], k = k_outer, seed = seed)
  
  all_preds <- vector("list", length(outer_test_folds))
  
  for (i in seq_along(outer_test_folds)) {
    test_id <- outer_test_folds[[i]]
    train_id <- setdiff(seq_along(y), test_id)
    
    # Check that no group appears in both train and test
    overlap <- intersect(meta[[group]][train_id], meta[[group]][test_id])
    stopifnot(length(overlap) == 0)
    
    taxa_prep <- prep_taxa_train(taxa[train_id, , drop = FALSE])
    x_train <- taxa_prep$train
    x_test  <- prep_taxa_test(taxa[test_id, , drop = FALSE], taxa_prep)
    
    if (residualize) {
      rz <- residualize_train_test(
        x_train = x_train,
        x_test = x_test,
        meta_train = meta[train_id, , drop = FALSE],
        meta_test = meta[test_id, , drop = FALSE],
        covariates = confounders
      )
      x_train <- rz$train
      x_test  <- rz$test
    }
    
    y_train <- droplevels(y[train_id])
    group_train <- meta[[group]][train_id]
    
    # Inner grouped folds for hyperparameter tuning, using only training samples
    inner_test_folds <- tryCatch(
      make_group_folds(y_train, group_train, k = k_inner, seed = seed + i),
      error = function(e) NULL
    )
    
    if (is.null(inner_test_folds)) {
      ctrl <- trainControl(
        method = "none",
        classProbs = TRUE,
        summaryFunction = twoClassSummary
      )
      
      tune_grid <- expand.grid(
        mtry = max(1, floor(sqrt(ncol(x_train)))),
        splitrule = "gini",
        min.node.size = 5
      )
      
    } else {
      inner_index <- lapply(inner_test_folds, function(te) {
        setdiff(seq_along(y_train), te)
      })
      
      ctrl <- trainControl(
        method = "cv",
        index = inner_index,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "final",
        allowParallel = TRUE
      )
      
      mtry_values <- unique(pmax(
        1,
        floor(c(
          sqrt(ncol(x_train)),
          ncol(x_train) / 10,
          ncol(x_train) / 3
        ))
      ))
      
      tune_grid <- expand.grid(
        mtry = mtry_values,
        splitrule = "gini",
        min.node.size = c(1, 5, 10)
      )
    }
    
    set.seed(seed + 1000 + i)
    
    fit <- train(
      x = x_train,
      y = y_train,
      method = "ranger",
      metric = "ROC",
      trControl = ctrl,
      tuneGrid = tune_grid,
      num.trees = 1000,
      importance = "permutation"
    )
    
    prob <- predict(fit, x_test, type = "prob")
    
    all_preds[[i]] <- bind_cols(
      tibble(
        sample_id = rownames(meta)[test_id],
        obs = y[test_id],
        fold = i,
        cv_group = as.character(meta[[group]][test_id])
      ),
      as_tibble(prob)
    )
  }
  
  bind_rows(all_preds)
}

run_nested_grouped_metadata_rf <- function(meta,
                                           covariates,
                                           outcome = "y",
                                           group = "cv_group",
                                           k_outer = 5,
                                           k_inner = 3,
                                           seed = 123) {
  
  y <- droplevels(meta[[outcome]])
  outer_test_folds <- make_group_folds(y, meta[[group]], k = k_outer, seed = seed)
  
  all_preds <- vector("list", length(outer_test_folds))
  
  for (i in seq_along(outer_test_folds)) {
    test_id <- outer_test_folds[[i]]
    train_id <- setdiff(seq_along(y), test_id)
    
    mm <- prep_covariate_mm(
      meta_train = meta[train_id, , drop = FALSE],
      meta_test = meta[test_id, , drop = FALSE],
      covariates = covariates
    )
    
    if (ncol(mm$train) < 1) {
      stop("No usable metadata covariates available for Confounder-only.")
    }
    
    x_train <- as.data.frame(mm$train)
    x_test  <- as.data.frame(mm$test)
    y_train <- droplevels(y[train_id])
    group_train <- meta[[group]][train_id]
    
    inner_test_folds <- tryCatch(
      make_group_folds(y_train, group_train, k = k_inner, seed = seed + i),
      error = function(e) NULL
    )
    
    if (is.null(inner_test_folds)) {
      ctrl <- trainControl(
        method = "none",
        classProbs = TRUE,
        summaryFunction = twoClassSummary
      )
      
      tune_grid <- expand.grid(
        mtry = max(1, floor(sqrt(ncol(x_train)))),
        splitrule = "gini",
        min.node.size = 5
      )
      
    } else {
      inner_index <- lapply(inner_test_folds, function(te) {
        setdiff(seq_along(y_train), te)
      })
      
      ctrl <- trainControl(
        method = "cv",
        index = inner_index,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "final",
        allowParallel = TRUE
      )
      
      tune_grid <- expand.grid(
        mtry = unique(pmax(1, floor(c(sqrt(ncol(x_train)), ncol(x_train) / 2)))),
        splitrule = "gini",
        min.node.size = c(1, 5, 10)
      )
    }
    
    set.seed(seed + 2000 + i)
    
    fit <- train(
      x = x_train,
      y = y_train,
      method = "ranger",
      metric = "ROC",
      trControl = ctrl,
      tuneGrid = tune_grid,
      num.trees = 1000,
      importance = "permutation"
    )
    
    prob <- predict(fit, x_test, type = "prob")
    
    all_preds[[i]] <- bind_cols(
      tibble(
        sample_id = rownames(meta)[test_id],
        obs = y[test_id],
        fold = i,
        cv_group = as.character(meta[[group]][test_id])
      ),
      as_tibble(prob)
    )
  }
  
  bind_rows(all_preds)
}


# 7. Main sensitivity analyses

# A. Optimistic: this approximates the optimistic analysis.
#    It is included only as a comparator.
meta_sample_cv <- meta0
meta_sample_cv$cv_group <- rownames(meta_sample_cv)

pred_sample_cv <- run_nested_grouped_rf(
  taxa = taxa0,
  meta = meta_sample_cv,
  outcome = "y",
  group = "cv_group",
  residualize = FALSE,
  k_outer = 5,
  k_inner = 3,
  seed = 123
)

summary_sample_cv <- summarise_binary(pred_sample_cv) %>%
  mutate(analysis = "A_naive_sample_level_CV")

# B. Leakage-aware: use cv_group to prevent related samples from being split.
pred_grouped_cv <- run_nested_grouped_rf(
  taxa = taxa0,
  meta = meta0,
  outcome = "y",
  group = "cv_group",
  residualize = FALSE,
  k_outer = 5,
  k_inner = 3,
  seed = 123
)

summary_grouped_cv <- summarise_binary(pred_grouped_cv) %>%
  mutate(analysis = "B_grouped_CV_no_covariate_adjustment")

# C. Confounder-only: can country be predicted from covariates alone?
#    A high AUC here means the microbial model is likely confounded.
pred_metadata_only <- run_nested_grouped_metadata_rf(
  meta = meta0,
  covariates = confounder_vars,
  outcome = "y",
  group = "cv_group",
  k_outer = 5,
  k_inner = 3,
  seed = 123
)

summary_metadata_only <- summarise_binary(pred_metadata_only) %>%
  mutate(analysis = "C_metadata_only_confounder_model")

# D. Confounder-adjusted model:
#    Taxa are transformed within each training fold, then residualised against
#    covariates using training-fold coefficients only.
pred_residualised <- run_nested_grouped_rf(
  taxa = taxa0,
  meta = meta0,
  outcome = "y",
  group = "cv_group",
  confounders = confounder_vars,
  residualize = TRUE,
  k_outer = 5,
  k_inner = 3,
  seed = 123
)

summary_residualised <- summarise_binary(pred_residualised) %>%
  mutate(analysis = "D_grouped_CV_covariate_residualised_microbiome")

# Save predictions and summary table
all_predictions <- bind_rows(
  pred_sample_cv %>% mutate(analysis = "A_naive_sample_level_CV"),
  pred_grouped_cv %>% mutate(analysis = "B_grouped_CV_no_covariate_adjustment"),
  pred_metadata_only %>% mutate(analysis = "C_metadata_only_confounder_model"),
  pred_residualised %>% mutate(analysis = "D_grouped_CV_covariate_residualised_microbiome")
)

all_summaries <- bind_rows(
  summary_sample_cv,
  summary_grouped_cv,
  summary_metadata_only,
  summary_residualised
) %>%
  dplyr::select(analysis, dplyr::everything())

write_tsv(
  all_summaries,
  file.path(out_dir, "06_model_performance_summary.tsv")
)


performance_plot <- all_summaries %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(
      analysis,
      "A_naive_sample_level_CV" =
        "Optimistic",
      "B_grouped_CV_no_covariate_adjustment" =
        "Leakage-aware",
      "C_metadata_only_confounder_model" =
        "Confounder-only",
      "D_grouped_CV_covariate_residualised_microbiome" =
        "Confounder-adjusted"
    ),
    analysis_label = factor(
      analysis_label,
      levels = c(
        "Optimistic",
        "Leakage-aware",
        "Confounder-only",
        "Confounder-adjusted"
      )
    )
  )

ggplot(performance_plot,
       aes(x = analysis_label, y = AUC)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = AUC_low, ymax = AUC_high),
    width = 0.15
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    x = NULL,
    y = "AUC-ROC"
  ) +
  theme_bw(base_size = 12)

roc_df <- all_predictions %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(
      analysis,
      "A_naive_sample_level_CV" =
        "Optimistic",
      "B_grouped_CV_no_covariate_adjustment" =
        "Leakage-aware",
      "C_metadata_only_confounder_model" =
        "Confounder-only",
      "D_grouped_CV_covariate_residualised_microbiome" =
        "Confounder-adjusted"
    ),
    analysis_label = factor(
      analysis_label,
      levels = c(
        "Optimistic",
        "Leakage-aware",
        "Confounder-only",
        "Confounder-adjusted"
      )
    )
  ) %>%
  dplyr::group_by(analysis_label) %>%
  dplyr::group_modify(~ {
    roc_obj <- pROC::roc(
      response = .x$obs,
      predictor = .x[[positive_class]],
      levels = c(negative_class, positive_class),
      direction = "<",
      quiet = TRUE
    )
    
    tibble::tibble(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      AUC = as.numeric(pROC::auc(roc_obj))
    )
  }) %>%
  dplyr::ungroup()

p_roc <- ggplot(
  roc_df,
  aes(x = 1 - specificity, y = sensitivity)
) +
  geom_path(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~ analysis_label) +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    title = NULL,
    subtitle = NULL
  ) +
  theme_bw(base_size = 12)

p_roc

roc_df <- all_predictions %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(
      analysis,
      "A_naive_sample_level_CV" =
        "Optimistic",
      "B_grouped_CV_no_covariate_adjustment" =
        "Leakage-aware",
      "C_metadata_only_confounder_model" =
        "Confounder-only",
      "D_grouped_CV_covariate_residualised_microbiome" =
        "Confounder-adjusted"
    ),
    analysis_label = factor(
      analysis_label,
      levels = c(
        "Optimistic",
        "Leakage-aware",
        "Confounder-only",
        "Confounder-adjusted"
      )
    )
  ) %>%
  dplyr::group_by(analysis_label) %>%
  dplyr::group_modify(~ {
    roc_obj <- pROC::roc(
      response = .x$obs,
      predictor = .x[[positive_class]],
      levels = c(negative_class, positive_class),
      direction = "<",
      quiet = TRUE
    )
    
    tibble::tibble(
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      AUC = as.numeric(pROC::auc(roc_obj))
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    auc_label = paste0(
      analysis_label,
      " (AUC=",
      sprintf("%.3f", AUC),
      ")"
    )
  )

palette_okabe <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#D55E00"  # vermillion
)

library(ggsci)

p_roc_multiple <- ggplot(
  roc_df,
  aes(
    x = 1 - specificity,
    y = sensitivity,
    colour = auc_label,
    group = auc_label
  )
) +
  geom_path(linewidth = 1.2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed"
  ) +
  coord_equal() +
  labs(
    x = "False positive rate",
    y = "True positive rate",
    colour = NULL
  ) +
  scale_colour_npg() +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10)
  )

p_roc_multiple

svglite(
  file = file.path(out_dir, "multiroc.svg"),
  width = 7.5,
  height = 7.5
)
print(p_roc_multiple)
dev.off()

performance_long <- all_summaries %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(
      analysis,
      "A_naive_sample_level_CV" =
        "Optimistic",
      "B_grouped_CV_no_covariate_adjustment" =
        "Leakage-aware",
      "C_metadata_only_confounder_model" =
        "Confounder-only",
      "D_grouped_CV_covariate_residualised_microbiome" =
        "Confounder-adjusted"
    ),
    analysis_label = factor(
      analysis_label,
      levels = c(
        "Optimistic",
        "Leakage-aware",
        "Confounder-only",
        "Confounder-adjusted"
      )
    )
  ) %>%
  dplyr::select(
    analysis_label,
    AUC,
    Accuracy,
    Sensitivity,
    Specificity,
    BalancedAccuracy
  ) %>%
  tidyr::pivot_longer(
    cols = c(AUC, Accuracy, Sensitivity, Specificity, BalancedAccuracy),
    names_to = "metric",
    values_to = "value"
  )

ggplot(performance_long,
       aes(x = analysis_label, y = value, shape = metric)) +
  geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
  coord_flip() +
  ylim(0.5, 1.0) +
  labs(
    x = NULL,
    y = "Performance",
    title = "Model performance across validation and confounding analyses"
  ) +
  theme_bw(base_size = 12)

library(scico)

p_perf_heatmap <- performance_long %>%
  ggplot(aes(x = metric, y = analysis_label, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.5) +
  coord_flip() +
  scale_fill_scico(
    palette = "lipari", 
    limits = c(0.73, 1)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

p_perf_heatmap

svglite(
  file = file.path(out_dir, "performance_heatmap.svg"),
  width = 4,
  height = 3
)
print(p_perf_heatmap)
dev.off()

p_perf_bubble <- performance_long %>%
  ggplot(
    aes(
      x = metric,
      y = analysis_label,
      size = value,
      fill = value
    )
  ) +
  geom_point(
    shape = 21,
    colour = "black",
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    size = 3,
    colour = "black"
  ) +
  scale_size(
    range = c(6, 18),
    limits = c(0.73, 1),
    name = "Performance"
  ) +
  scale_fill_viridis_c(
    option = "C",
    limits = c(0.73, 1),
    name = "Performance"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "right"
  )

p_perf_bubble

analysis_labels <- c(
  "A_naive_sample_level_CV" =
    "Naive sample-level CV",
  
  "B_grouped_CV_no_covariate_adjustment" =
    "Grouped CV",
  
  "C_metadata_only_confounder_model" =
    "Confounder-only",
  
  "D_grouped_CV_covariate_residualised_microbiome" =
    "Confounder-adjusted"
)

prob_df <- all_predictions %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(analysis, !!!analysis_labels),
    analysis_label = factor(
      analysis_label,
      levels = unname(analysis_labels)
    ),
    prob_Italy = .data[[positive_class]]
  )

p_prob <- ggplot(
  prob_df,
  aes(x = obs, y = prob_Italy)
) +
  #geom_violin(trim = TRUE) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 0.8) +
  geom_boxplot(width = 0.15, alpha = 0.25, outlier.shape = NA) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  facet_wrap(~ analysis_label) +
  labs(
    x = "Observed class",
    y = paste0("Predicted probability of ", positive_class)
    ) +
  theme_bw(base_size = 12)

p_prob

country_cols <- c(
  "Spain" = "firebrick",
  "Italy" = "cornflowerblue"
)

prob_df <- all_predictions %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(analysis, !!!analysis_labels),
    analysis_label = factor(
      analysis_label,
      levels = unname(analysis_labels)
    ),
    prob_Italy = .data[[positive_class]],
    obs = factor(obs, levels = c("Italy", "Spain"))
  )

p_prob <- ggplot(
  prob_df,
  aes(x = obs, y = prob_Italy, fill = obs)
) +
  geom_jitter(
    aes(colour = obs),
    width = 0.05,
    alpha = 0.5,
    size = 0.8
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey35"
  ) +
  scale_fill_manual(values = country_cols) +
  scale_colour_manual(values = country_cols) +
  facet_wrap(~ analysis_label, ncol = 2) +
  labs(
    x = NULL,
    y = paste0("Predicted probability of ", positive_class)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

p_prob

svglite(
  file = file.path(out_dir, "pred_prob.svg"),
  width = 3.5,
  height = 4.5
)
print(p_prob)
dev.off()


fold_metrics <- all_predictions %>%
  dplyr::group_by(analysis, fold) %>%
  dplyr::group_modify(~ {
    summarise_binary(.x) %>%
      dplyr::select(
        AUC,
        Accuracy,
        Sensitivity,
        Specificity,
        BalancedAccuracy
      )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    analysis_label = dplyr::recode(analysis, !!!analysis_labels),
    analysis_label = factor(
      analysis_label,
      levels = unname(analysis_labels)
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(AUC, Accuracy, Sensitivity, Specificity, BalancedAccuracy),
    names_to = "metric",
    values_to = "value"
  )

p_resamples <- ggplot(
  fold_metrics,
  aes(x = analysis_label, y = value)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  #coord_flip() +
  labs(
    x = NULL,
    y = "Fold-level performance",
    title = "Cross-validation fold stability",
    subtitle = "Boxplots show variability across held-out folds"
  ) +
  theme_bw(base_size = 12)

p_resamples




library(caret)
library(ranger)

# Prepare all microbiome features using the same transformation function
taxa_prep_all <- prep_taxa_train(taxa0, min_prevalence = 0.05)

x_all <- taxa_prep_all$train
y_all <- meta0$y

feature_map <- tibble(
  feature_safe = taxa_prep_all$features_safe,
  feature_original = taxa_prep_all$features_original
)

# Leakage-aware indices for model tuning
outer_test_folds <- make_group_folds(
  y = y_all,
  group = meta0$cv_group,
  k = 5,
  seed = 123
)

outer_train_index <- lapply(
  outer_test_folds,
  function(test_id) setdiff(seq_along(y_all), test_id)
)

rf_ctrl <- trainControl(
  method = "cv",
  index = outer_train_index,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)

rf_microbiome_final <- train(
  x = x_all,
  y = y_all,
  method = "ranger",
  metric = "ROC",
  trControl = rf_ctrl,
  tuneGrid = expand.grid(
    mtry = unique(pmax(1, floor(c(
      sqrt(ncol(x_all)),
      ncol(x_all) / 10,
      ncol(x_all) / 3
    )))),
    splitrule = "gini",
    min.node.size = c(1, 5, 10)
  ),
  num.trees = 1000,
  importance = "permutation"
)

rf_varimp <- varImp(rf_microbiome_final, scale = FALSE)

plot(rf_varimp, top = 10, main = "Random forest variable importance: microbiome model")

make_species_label <- function(x) {
  
  # Extract species name after s__
  sp <- stringr::str_extract(x, "s__[^;]+")
  
  # Remove s__
  sp <- stringr::str_remove(sp, "^s__")
  
  # Replace underscores with spaces
  sp <- stringr::str_replace_all(sp, "_", " ")
  
  # If no species found, keep original
  sp[is.na(sp)] <- x[is.na(sp)]
  
  # Convert to italic expression
  paste0("italic('", sp, "')")
}


imp_df <- rf_varimp$importance %>%
  tibble::rownames_to_column("feature") %>%
  dplyr::rename(importance = Overall) %>%
  dplyr::arrange(desc(importance)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(
    label = make_species_label(feature)
  )

p_varimp_italic <- ggplot(
  imp_df,
  aes(
    x = reorder(label, importance),
    y = importance
  )
) +
  geom_col() +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) parse(text = x)
  ) +
  labs(
    x = NULL,
    y = "Permutation importance",
    title = "Random forest variable importance: microbiome model"
  ) +
  theme_bw(base_size = 12)

p_varimp_italic

p_varimp_lollipop <- ggplot(
  imp_df,
  aes(
    x = reorder(label, importance),
    y = importance
  )
) +
  # Lollipop stems
  geom_segment(
    aes(
      xend = reorder(label, importance),
      y = 0,
      yend = importance
    ),
    linewidth = 0.8,
    colour = "grey50"
  ) +
  
  # Lollipop heads
  geom_point(
    size = 5,
    shape = 21,
    fill = "#2C7FB8",
    colour = "black",
    stroke = 0.4
  ) +
  
  coord_flip() +
  
  scale_x_discrete(
    labels = function(x) parse(text = x)
  ) +
  
  labs(
    x = NULL,
    y = "Variable importance"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p_varimp_lollipop

svglite(
  file = file.path(out_dir, "p_varimp_lollipop.svg"),
  width = 5,
  height = 3
)
print(p_varimp_lollipop)
dev.off()

mm_all <- prep_covariate_mm(
  meta_train = meta0,
  meta_test = meta0,
  covariates = confounder_vars
)

x_meta_all <- as.data.frame(mm_all$train)

set.seed(123)

rf_metadata_final <- train(
  x = x_meta_all,
  y = y_all,
  method = "ranger",
  metric = "ROC",
  trControl = rf_ctrl,
  tuneGrid = expand.grid(
    mtry = unique(pmax(1, floor(c(
      sqrt(ncol(x_meta_all)),
      ncol(x_meta_all) / 2
    )))),
    splitrule = "gini",
    min.node.size = c(1, 5, 10)
  ),
  num.trees = 1000,
  importance = "permutation"
)

rf_metadata_varimp <- varImp(rf_metadata_final, scale = FALSE)

imp_metadata <- rf_metadata_varimp$importance %>%
  tibble::rownames_to_column("feature") %>%
  dplyr::rename(importance = Overall) %>%
  dplyr::arrange(desc(importance)) %>%
  dplyr::slice_head(n = 30)

p_varimp_metadata <- ggplot(
  imp_metadata,
  aes(
    x = reorder(feature, importance),
    y = importance
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Permutation importance",
    title = "Top metadata predictors of country",
    subtitle = "High metadata-only performance indicates strong covariate structure"
  ) +
  theme_bw(base_size = 12)

p_varimp_metadata

short_taxon <- function(x) {
  x %>%
    stringr::str_replace("^.*s__", "s__") %>%
    stringr::str_replace("^.*g__", "g__") %>%
    stringr::str_replace("^.*f__", "f__") %>%
    stringr::str_replace_all("_", " ")
}

imp_microbiome <- rf_varimp$importance %>%
  tibble::rownames_to_column("feature_safe") %>%
  dplyr::rename(importance = Overall) %>%
  dplyr::left_join(feature_map, by = "feature_safe") %>%
  dplyr::mutate(
    feature_label = short_taxon(feature_original)
  ) %>%
  dplyr::arrange(desc(importance))

top3_features <- imp_microbiome %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::pull(feature_safe)



featurePlot(
  x = x_all[, top3_features, drop = FALSE],
  y = y_all,
  plot = "box",
  strip = strip.custom(
    par.strip.text = list(cex = 0.7)
  ),
  scales = list(
    x = list(relation = "free"),
    y = list(relation = "free")
  ),
  main = "Top microbiome features by country"
)


top3_map <- imp_microbiome %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::select(feature_safe, feature_original, feature_label)

feature_box_df <- x_all %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, dplyr::all_of(top3_map$feature_safe)) %>%
  tidyr::pivot_longer(
    cols = -sample_id,
    names_to = "feature_safe",
    values_to = "clr_abundance"
  ) %>%
  dplyr::left_join(top3_map, by = "feature_safe") %>%
  dplyr::left_join(
    meta0 %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::select(
        sample_id,
        Country,
        Milk_processing,
        Technology_of_production,
        Rheological_properties
      ),
    by = "sample_id"
  )

p_feature_box <- ggplot(
  feature_box_df,
  aes(x = Country, y = clr_abundance)
) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.7) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  labs(
    x = NULL,
    y = "CLR-transformed relative abundance",
    title = "Top microbiome predictors differ between Italy and Spain",
    subtitle = "Features selected from final RF model; visualisation is exploratory"
  ) +
  theme_bw(base_size = 11)

p_feature_box


country_cols <- c(
  "Spain" = "firebrick",
  "Italy" = "cornflowerblue"
)

p_feature_box <- ggplot(
  feature_box_df,
  aes(x = Country, y = clr_abundance, fill = Country)
) +
  geom_jitter(
    aes(colour = Country),
    width = 0.15,
    alpha = 0.5,
    size = 0.8
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.55
  ) +
  #ggpubr::stat_compare_means(
    #comparisons = list(c("Spain", "Italy")),
    #method = "wilcox.test",
    #label = "p.format",
    #size = 3
  #) +
  scale_fill_manual(values = country_cols) +
  scale_colour_manual(values = country_cols) +
  facet_wrap(~ feature_label, scales = "free_y", ncol = 3) +
  labs(
    x = NULL,
    y = "CLR-transformed relative abundance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic"),
    panel.grid.minor = element_blank()
  )

p_feature_box

svglite(
  file = file.path(out_dir, "p_feature_box.svg"),
  width = 6,
  height = 6
)
print(p_feature_box)
dev.off()




p_confounders <- confounder_screen %>%
  dplyr::mutate(
    variable = factor(variable, levels = rev(variable)),
    neg_log10_q = -log10(q_value)
  ) %>%
  ggplot(
    aes(x = variable, y = cramers_v)
  ) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Cramer's V with country",
    title = "Cheese metadata variables associated with country",
    subtitle = "Strong associations indicate confounding structure in Italy vs Spain comparison"
  ) +
  theme_bw(base_size = 12)

p_confounders


## ==== 7.1 One-vs-One: Italy, Spain, Ireland, and Austria ====

# Supplementary Excel: all one-vs-one country comparisons
# Place this AFTER run_nested_grouped_rf(), run_nested_grouped_metadata_rf(),
# summarise_binary(), and confounder_vars have been defined.

countries <- c("Italy", "Spain", "Ireland", "Austria")
country_pairs <- combn(countries, 2, simplify = FALSE)

all_pairwise_long <- list()
all_pairwise_wide <- list()

for (pair in country_pairs) {
  
  positive_class <- pair[1]
  negative_class <- pair[2]
  comparison_name <- paste0(positive_class, "_vs_", negative_class)
  
  cat("\n====================\n")
  cat("Processing:", positive_class, "vs", negative_class, "\n")
  
  meta_pair <- metadata %>%
    dplyr::filter(
      Part_of_cheese == "Core",
      .data[[target_var]] %in% c(positive_class, negative_class)
    ) %>%
    dplyr::mutate(
      y = factor(
        dplyr::if_else(
          .data[[target_var]] == positive_class,
          positive_class,
          negative_class
        ),
        levels = c(positive_class, negative_class)
      )
    )
  
  taxa_pair <- taxa_table[rownames(meta_pair), , drop = FALSE]
  
  cat("Samples used:\n")
  print(table(meta_pair$y))
  
  # Use the same grouping logic as the original script
  if (!is.na(preferred_group_var) && preferred_group_var %in% names(meta_pair)) {
    meta_pair$cv_group <- as.character(meta_pair[[preferred_group_var]])
    meta_pair$cv_group[
      is.na(meta_pair$cv_group) | meta_pair$cv_group == ""
    ] <- rownames(meta_pair)[
      is.na(meta_pair$cv_group) | meta_pair$cv_group == ""
    ]
  } else {
    meta_pair$cv_group <- rownames(meta_pair)
  }
  
  confounder_vars_pair <- intersect(confounder_vars, names(meta_pair))
  
  # A. Optimistic sample-level CV
  meta_sample_cv <- meta_pair
  meta_sample_cv$cv_group <- rownames(meta_sample_cv)
  
  pred_sample_cv <- run_nested_grouped_rf(
    taxa = taxa_pair,
    meta = meta_sample_cv,
    outcome = "y",
    group = "cv_group",
    residualize = FALSE,
    k_outer = 5,
    k_inner = 3,
    seed = 123
  )
  
  summary_sample_cv <- summarise_binary(
    pred_sample_cv,
    positive = positive_class,
    negative = negative_class
  ) %>%
    dplyr::mutate(analysis = "A_naive_sample_level_CV")
  
  # B. Leakage-aware grouped CV
  pred_grouped_cv <- run_nested_grouped_rf(
    taxa = taxa_pair,
    meta = meta_pair,
    outcome = "y",
    group = "cv_group",
    residualize = FALSE,
    k_outer = 5,
    k_inner = 3,
    seed = 123
  )
  
  summary_grouped_cv <- summarise_binary(
    pred_grouped_cv,
    positive = positive_class,
    negative = negative_class
  ) %>%
    dplyr::mutate(analysis = "B_grouped_CV_no_covariate_adjustment")
  
  # C. Metadata-only confounder model
  pred_metadata_only <- run_nested_grouped_metadata_rf(
    meta = meta_pair,
    covariates = confounder_vars_pair,
    outcome = "y",
    group = "cv_group",
    k_outer = 5,
    k_inner = 3,
    seed = 123
  )
  
  summary_metadata_only <- summarise_binary(
    pred_metadata_only,
    positive = positive_class,
    negative = negative_class
  ) %>%
    dplyr::mutate(analysis = "C_metadata_only_confounder_model")
  
  # D. Confounder-adjusted microbiome model
  pred_residualised <- run_nested_grouped_rf(
    taxa = taxa_pair,
    meta = meta_pair,
    outcome = "y",
    group = "cv_group",
    confounders = confounder_vars_pair,
    residualize = TRUE,
    k_outer = 5,
    k_inner = 3,
    seed = 123
  )
  
  summary_residualised <- summarise_binary(
    pred_residualised,
    positive = positive_class,
    negative = negative_class
  ) %>%
    dplyr::mutate(analysis = "D_grouped_CV_covariate_residualised_microbiome")
  
  all_summaries_pair <- dplyr::bind_rows(
    summary_sample_cv,
    summary_grouped_cv,
    summary_metadata_only,
    summary_residualised
  ) %>%
    dplyr::mutate(
      comparison = comparison_name,
      positive_class = positive_class,
      negative_class = negative_class,
      n_positive = sum(meta_pair$y == positive_class),
      n_negative = sum(meta_pair$y == negative_class)
    ) %>%
    dplyr::select(
      comparison,
      positive_class,
      negative_class,
      n_positive,
      n_negative,
      analysis,
      dplyr::everything()
    )
  
  performance_long_pair <- all_summaries_pair %>%
    dplyr::mutate(
      analysis_label = dplyr::recode(
        analysis,
        "A_naive_sample_level_CV" = "Optimistic",
        "B_grouped_CV_no_covariate_adjustment" = "Leakage-aware",
        "C_metadata_only_confounder_model" = "Confounder-only",
        "D_grouped_CV_covariate_residualised_microbiome" = "Confounder-adjusted"
      ),
      analysis_label = factor(
        analysis_label,
        levels = c(
          "Optimistic",
          "Leakage-aware",
          "Confounder-only",
          "Confounder-adjusted"
        )
      )
    ) %>%
    dplyr::select(
      comparison,
      positive_class,
      negative_class,
      n_positive,
      n_negative,
      analysis,
      analysis_label,
      AUC,
      Accuracy,
      Sensitivity,
      Specificity,
      BalancedAccuracy
    ) %>%
    tidyr::pivot_longer(
      cols = c(AUC, Accuracy, Sensitivity, Specificity, BalancedAccuracy),
      names_to = "metric",
      values_to = "value"
    )
  
  all_pairwise_wide[[comparison_name]] <- all_summaries_pair
  all_pairwise_long[[comparison_name]] <- performance_long_pair
}

performance_long_all_pairs <- dplyr::bind_rows(all_pairwise_long)
performance_wide_all_pairs <- dplyr::bind_rows(all_pairwise_wide)

writexl::write_xlsx(
  list(
    performance_long = performance_long_all_pairs,
    performance_wide = performance_wide_all_pairs
  ),
  path = file.path(out_dir, "performance_long_all_pairs.xlsx")
)

cat("\nSaved supplementary Excel file:\n")
cat(file.path(out_dir, "performance_long_all_pairs.xlsx"), "\n", sep = "")





# 8. SHAP =============================================================================

# SHAP interpretation is performed on a final RF microbiome model trained on all
# Italy-vs-Spain core samples. This model is used for feature interpretation only;
# performance estimates are taken from the held-out cross-validation analyses above.

if (run_shap) {

  positive_class <- "Italy"
  negative_class <- "Spain"

  cat("\nRunning kernel SHAP for Italy vs Spain microbiome model...\n")

  pred_fun_prob <- function(object, newdata) {
    as.matrix(
      predict(
        object,
        newdata = as.data.frame(newdata),
        type = "prob"
      )
    )
  }

  # Stratified background distribution. This keeps SHAP computationally tractable
  # while representing both classes.
  set.seed(123)
  bg_size <- min(100, nrow(x_all))
  bg_index <- unlist(
    tapply(seq_len(nrow(x_all)), y_all, function(idx) {
      sample(idx, size = min(length(idx), ceiling(bg_size / 2)))
    }),
    use.names = FALSE
  )
  bg_index <- unique(bg_index)

  bg_data <- x_all[bg_index, , drop = FALSE]
  X_explain <- x_all

  ks_IT_ES <- kernelshap::kernelshap(
    object = rf_microbiome_final,
    X = X_explain,
    bg_X = bg_data,
    pred_fun = pred_fun_prob,
    verbose = TRUE
  )

  saveRDS(
    ks_IT_ES,
    file = file.path(out_dir, "kernelshap_Italy_vs_Spain.rds")
  )

  sv_Country <- shapviz::shapviz(ks_IT_ES)

  # Overall importance plot.
  p_shap_importance <- shapviz::sv_importance(sv_Country)
  ggsave(
    filename = file.path(out_dir, "9_3_shap_Country.svg"),
    plot = p_shap_importance,
    width = 10,
    height = 5,
    dpi = 600
  )

  # Beeswarm plot for all classes/outputs.
  p_shap_bee <- shapviz::sv_importance(sv_Country, kind = "bee")
  ggsave(
    filename = file.path(out_dir, "9_3_shap_Country_bee.svg"),
    plot = p_shap_bee,
    width = 25,
    height = 5,
    dpi = 600
  )

  # Class-specific beeswarm plots, if the multi-output object exposes class names.
  if ("Italy" %in% names(sv_Country)) {
    p_bee_Italy <- shapviz::sv_importance(sv_Country$Italy, kind = "bee") +
      scale_colour_gradient(low = "grey90", high = "#6495edff", na.value = "grey90") +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        legend.direction = "horizontal",
        axis.text.y      = element_text(face = "italic")
      ) +
      guides(
        colour = guide_colourbar(
          title.position = "top",
          title.hjust    = 0.5,
          barwidth       = unit(10, "cm"),
          barheight      = unit(0.15, "cm")
        )
      )

    ggsave(
      filename = file.path(out_dir, "SHAP_Italy.svg"),
      plot = p_bee_Italy,
      width = 17.5,
      height = 15,
      units = "cm",
      dpi = 600
    )
  }

  if ("Spain" %in% names(sv_Country)) {
    p_bee_Spain <- shapviz::sv_importance(sv_Country$Spain, kind = "bee") +
      scale_colour_gradient(low = "grey90", high = "#b22222ff", na.value = "grey90") +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "top",
        legend.direction = "horizontal",
        axis.text.y      = element_text(face = "italic")
      ) +
      guides(
        colour = guide_colourbar(
          title.position = "top",
          title.hjust    = 0.5,
          barwidth       = unit(10, "cm"),
          barheight      = unit(0.15, "cm")
        )
      )

    ggsave(
      filename = file.path(out_dir, "SHAP_Spain.svg"),
      plot = p_bee_Spain,
      width = 17.5,
      height = 15,
      units = "cm",
      dpi = 600
    )
  }

} else {
  cat("\nSkipping SHAP because --run-shap was set to FALSE.\n")
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

cat("\nAnalysis complete. Outputs saved in: ", out_dir, "\n", sep = "")
