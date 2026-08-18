# Fig. 4


# 0. Libraries
library(tidyverse)
library(vegan)
library(ggstatsplot)
library(ggsignif)
library(rstatix)
library(svglite)


# 1. Output directory
dir.create(
  "output/Fig4",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Load public MetaCheeseDB data
metadata <- read.delim(
  "data/MetaCheeseDB_metadata.tsv",
  row.names = 1,
  check.names = FALSE
)

taxonomy <- read.delim(
  "data/MetaCheeseDB_taxonomy.tsv",
  row.names = 1,
  check.names = FALSE
)

stopifnot(
  identical(
    rownames(metadata),
    rownames(taxonomy)
  )
)


# 3. Load alpha-diversity results
alpha_div <- read.table(
  "data/alpha_SampleID.csv",
  header = TRUE,
  sep = ",",
  check.names = FALSE,
  dec = ".",
  comment.char = "",
  quote = ""
)


cat(
  "alpha_div: ",
  nrow(alpha_div),
  " rows, ",
  ncol(alpha_div),
  " columns\n",
  sep = ""
)

print(
  head(
    alpha_div,
    3
  )
)



# 4. Merge alpha diversity with metadata and retain core samples
metadata_for_merge <- metadata %>%
  rownames_to_column(
    var = "SampleID"
  )


common_ids <- intersect(
  alpha_div$SampleID,
  metadata_for_merge$SampleID
)


cat(
  "Sample IDs in alpha_div: ",
  length(
    unique(
      alpha_div$SampleID
    )
  ),
  "\n",
  "Sample IDs in metadata: ",
  length(
    unique(
      metadata_for_merge$SampleID
    )
  ),
  "\n",
  "Sample IDs in BOTH: ",
  length(common_ids),
  "\n",
  sep = ""
)


alpha_meta <- alpha_div %>%
  inner_join(
    metadata_for_merge,
    by = "SampleID"
  )


cat(
  "Merged: ",
  nrow(alpha_meta),
  " samples\n",
  sep = ""
)


alpha_core <- alpha_meta %>%
  mutate(
    Part_of_cheese = na_if(
      na_if(
        as.character(
          Part_of_cheese
        ),
        ""
      ),
      "NA"
    )
  ) %>%
  filter(
    Part_of_cheese == "Core"
  ) %>%
  mutate(
    Country = factor(
      Country
    )
  )


cat(
  "Core samples: ",
  nrow(alpha_core),
  "\n",
  sep = ""
)

cat(
  "Samples per country (Core only):\n"
)

print(
  table(
    alpha_core$Country
  )
)



# 5. Retain countries represented by at least 30 core samples
country_counts <- alpha_core %>%
  count(
    Country,
    name = "n_country"
  )


cat(
  "Sample counts per country (Core only):\n"
)

print(
  country_counts %>%
    arrange(
      desc(n_country)
    )
)


countries_keep <- country_counts %>%
  filter(
    n_country >= 30
  ) %>%
  pull(
    Country
  ) %>%
  as.character()


cat(
  "\nCountries retained (n >= 30): ",
  paste(
    countries_keep,
    collapse = ", "
  ),
  "\n",
  sep = ""
)


cat(
  "Countries dropped (n < 30): ",
  paste(
    setdiff(
      as.character(
        country_counts$Country
      ),
      countries_keep
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)


alpha_core_30 <- alpha_core %>%
  filter(
    Country %in% countries_keep
  ) %>%
  mutate(
    Country = factor(
      Country,
      levels = countries_keep
    )
  )


cat(
  "\nFinal sample count: ",
  nrow(alpha_core_30),
  "\n",
  sep = ""
)

print(
  table(
    alpha_core_30$Country
  )
)



# 6. Country colours
my_col_individual <- c(
  "Spain"        = "firebrick",
  "France"       = "cadetblue",
  "Ireland"      = "darkslategray",
  "Austria"      = "sandybrown",
  "USA"          = "coral3",
  "Italy"        = "cornflowerblue",
  "Canada"       = "forestgreen",
  "<_30_samples" = "grey80"
)


# 7. Fig. 4A | Simpson
pairwise_results_Simpson <- alpha_core_30 %>%
  rstatix::wilcox_test(
    Simpson ~ Country,
    p.adjust.method = "fdr"
  ) %>%
  rstatix::add_significance(
    "p.adj",
    cutpoints = c(
      0,
      0.0001,
      0.001,
      0.01,
      0.05,
      1
    ),
    symbols = c(
      "****",
      "***",
      "**",
      "*",
      "ns"
    )
  ) %>%
  filter(
    p.adj.signif != "ns"
  )


write.table(
  pairwise_results_Simpson,
  "output/Fig4/Fig4A_pairwise_Simpson.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


Simpson_Country <- ggbetweenstats(
  data = alpha_core_30,
  x = Country,
  y = Simpson,
  type = "nonparametric",

  pairwise.display = "none",
  results.subtitle = TRUE,
  sample.size.label = TRUE,
  centrality.plotting = FALSE,

  point.args = list(
    position = position_jitterdodge(
      dodge.width = 0.5
    ),
    alpha = 0.75,
    size = 2,
    stroke = 0.5,
    na.rm = TRUE
  ),

  boxplot.args = list(
    width = 0.5,
    alpha = 0.2,
    na.rm = TRUE
  ),

  violin.args = list(
    width = 0.5,
    alpha = 0.1,
    na.rm = TRUE
  ),

  ggtheme = theme_ggstatsplot()
) +

  ggsignif::geom_signif(
    comparisons = lapply(
      seq_len(
        nrow(
          pairwise_results_Simpson
        )
      ),
      function(i) {
        c(
          pairwise_results_Simpson$group1[i],
          pairwise_results_Simpson$group2[i]
        )
      }
    ),
    annotations = pairwise_results_Simpson$p.adj.signif,
    map_signif_level = FALSE,
    step_increase = 0.04,
    textsize = 4,
    tip_length = 0.01,
    vjust = 0.9
  ) +

  scale_color_manual(
    values = my_col_individual
  ) +

  scale_fill_manual(
    values = my_col_individual
  ) +

  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(
      colour = "grey50"
    ),
    panel.grid = element_line(
      color = "#b4aea9"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      linetype = "dashed"
    )
  ) +

  labs(
    title = NULL,
    subtitle = NULL
  )


Simpson_Country


svglite(
  file = "output/Fig4/Fig4A_alpha_Simpson_Country_Core.svg",
  width = 3.5,
  height = 6
)

print(
  Simpson_Country
)

dev.off()


# 8. Fig. 4B | Shannon
pairwise_results_Shannon <- alpha_core_30 %>%
  rstatix::wilcox_test(
    Shannon ~ Country,
    p.adjust.method = "fdr"
  ) %>%
  rstatix::add_significance(
    "p.adj",
    cutpoints = c(
      0,
      0.0001,
      0.001,
      0.01,
      0.05,
      1
    ),
    symbols = c(
      "****",
      "***",
      "**",
      "*",
      "ns"
    )
  ) %>%
  filter(
    p.adj.signif != "ns"
  )


write.table(
  pairwise_results_Shannon,
  "output/Fig4/Fig4B_pairwise_Shannon.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


shannon_Country <- ggbetweenstats(
  data = alpha_core_30,
  x = Country,
  y = Shannon,
  type = "nonparametric",

  pairwise.display = "none",
  results.subtitle = TRUE,
  sample.size.label = TRUE,
  centrality.plotting = FALSE,

  point.args = list(
    position = position_jitterdodge(
      dodge.width = 0.5
    ),
    alpha = 0.75,
    size = 2,
    stroke = 0.5,
    na.rm = TRUE
  ),

  boxplot.args = list(
    width = 0.5,
    alpha = 0.2,
    na.rm = TRUE
  ),

  violin.args = list(
    width = 0.5,
    alpha = 0.1,
    na.rm = TRUE
  ),

  ggtheme = theme_ggstatsplot()
) +

  ggsignif::geom_signif(
    comparisons = lapply(
      seq_len(
        nrow(
          pairwise_results_Shannon
        )
      ),
      function(i) {
        c(
          pairwise_results_Shannon$group1[i],
          pairwise_results_Shannon$group2[i]
        )
      }
    ),
    annotations = pairwise_results_Shannon$p.adj.signif,
    map_signif_level = FALSE,
    step_increase = 0.04,
    textsize = 4,
    tip_length = 0.01,
    vjust = 0.9
  ) +

  scale_color_manual(
    values = my_col_individual
  ) +

  scale_fill_manual(
    values = my_col_individual
  ) +

  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(
      colour = "grey50"
    ),
    panel.grid = element_line(
      color = "#b4aea9"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      linetype = "dashed"
    )
  ) +

  labs(
    title = NULL,
    subtitle = NULL
  )


shannon_Country


svglite(
  file = "output/Fig4/Fig4B_alpha_Shannon_Country_Core.svg",
  width = 3.5,
  height = 6
)

print(
  shannon_Country
)

dev.off()


# 9. Fig. 4C | Pielou
pairwise_results_Pielou <- alpha_core_30 %>%
  rstatix::wilcox_test(
    Pielou ~ Country,
    p.adjust.method = "fdr"
  ) %>%
  rstatix::add_significance(
    "p.adj",
    cutpoints = c(
      0,
      0.0001,
      0.001,
      0.01,
      0.05,
      1
    ),
    symbols = c(
      "****",
      "***",
      "**",
      "*",
      "ns"
    )
  ) %>%
  filter(
    p.adj.signif != "ns"
  )


write.table(
  pairwise_results_Pielou,
  "output/Fig4/Fig4C_pairwise_Pielou.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


Pielou_Country <- ggbetweenstats(
  data = alpha_core_30,
  x = Country,
  y = Pielou,
  type = "nonparametric",

  pairwise.display = "none",
  results.subtitle = TRUE,
  sample.size.label = TRUE,
  centrality.plotting = FALSE,

  point.args = list(
    position = position_jitterdodge(
      dodge.width = 0.5
    ),
    alpha = 0.75,
    size = 2,
    stroke = 0.5,
    na.rm = TRUE
  ),

  boxplot.args = list(
    width = 0.5,
    alpha = 0.2,
    na.rm = TRUE
  ),

  violin.args = list(
    width = 0.5,
    alpha = 0.1,
    na.rm = TRUE
  ),

  ggtheme = theme_ggstatsplot()
) +

  ggsignif::geom_signif(
    comparisons = lapply(
      seq_len(
        nrow(
          pairwise_results_Pielou
        )
      ),
      function(i) {
        c(
          pairwise_results_Pielou$group1[i],
          pairwise_results_Pielou$group2[i]
        )
      }
    ),
    annotations = pairwise_results_Pielou$p.adj.signif,
    map_signif_level = FALSE,
    step_increase = 0.04,
    textsize = 4,
    tip_length = 0.01,
    vjust = 0.9
  ) +

  scale_color_manual(
    values = my_col_individual
  ) +

  scale_fill_manual(
    values = my_col_individual
  ) +

  theme(
    axis.ticks = element_blank(),
    axis.line = element_line(
      colour = "grey50"
    ),
    panel.grid = element_line(
      color = "#b4aea9"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      linetype = "dashed"
    )
  ) +

  labs(
    title = NULL,
    subtitle = NULL
  )


Pielou_Country


svglite(
  file = "output/Fig4/Fig4C_alpha_Pielou_Country_Core.svg",
  width = 3.5,
  height = 6
)

print(
  Pielou_Country
)

dev.off()


# 10. Fig. 4D | Bray-Curtis UMAP of core samples by country
metadata_core <- metadata %>%
  filter(
    Part_of_cheese == "Core"
  )

taxa_core <- taxonomy[
  rownames(metadata_core),
  ,
  drop = FALSE
]

core_country <- sort(
  table(
    metadata_core$Country
  ),
  decreasing = TRUE
)

print(
  core_country
)

viable_countries <- names(
  core_country[
    core_country >= 30
  ]
)

print(
  viable_countries
)

metadata_core_balanced <- metadata_core %>%
  filter(
    Country %in% viable_countries
  ) %>%
  mutate(
    Country = factor(
      Country,
      levels = viable_countries
    )
  )

taxa_core_balanced <- taxonomy[
  rownames(metadata_core_balanced),
  ,
  drop = FALSE
]

cat(
  "\nCore samples before Bray-Curtis cleanup:\n"
)

print(
  table(
    metadata_core_balanced$Country
  )
)

# 10.2 Align taxa and metadata
common_samples <- intersect(
  rownames(taxa_core_balanced),
  rownames(metadata_core_balanced)
)

cat(
  "Samples in taxa: ",
  nrow(taxa_core_balanced),
  "\n",
  "Samples in metadata: ",
  nrow(metadata_core_balanced),
  "\n",
  "Samples in both: ",
  length(common_samples),
  "\n",
  sep = ""
)

taxa_pcoa <- taxa_core_balanced[
  common_samples,
  ,
  drop = FALSE
]

meta_pcoa <- metadata_core_balanced[
  common_samples,
  ,
  drop = FALSE
]

# 10.3 Bray-Curtis distance
if (any(is.na(taxa_pcoa))) {

  cat(
    "Warning: ",
    sum(is.na(taxa_pcoa)),
    " NA values in taxa table - replacing with 0\n",
    sep = ""
  )

  taxa_pcoa[
    is.na(taxa_pcoa)
  ] <- 0
}


if (any(taxa_pcoa < 0, na.rm = TRUE)) {

  stop(
    "Negative values detected in taxa table. ",
    "Bray-Curtis requires non-negative abundance data."
  )
}

nonzero_species <- colSums(
  taxa_pcoa
) > 0

cat(
  "Species retained (non-zero): ",
  sum(nonzero_species),
  " of ",
  length(nonzero_species),
  "\n",
  sep = ""
)

taxa_pcoa <- taxa_pcoa[
  ,
  nonzero_species,
  drop = FALSE
]

sample_sums <- rowSums(
  taxa_pcoa
)

cat(
  "Samples with row sum = 0: ",
  sum(sample_sums == 0),
  "\n",
  "Samples with row sum > 0: ",
  sum(sample_sums > 0),
  "\n",
  sep = ""
)

zero_sum_samples <- rownames(taxa_pcoa)[
  sample_sums == 0
]

if (length(zero_sum_samples) > 0) {

  cat(
    "\nZero-sum samples (",
    length(zero_sum_samples),
    ") removed before Bray-Curtis.\n",
    sep = ""
  )
}

keep_samples <- sample_sums > 0

taxa_pcoa <- taxa_pcoa[
  keep_samples,
  ,
  drop = FALSE
]

meta_pcoa <- meta_pcoa[
  keep_samples,
  ,
  drop = FALSE
]

species_keep <- colSums(
  taxa_pcoa
) > 0

taxa_pcoa <- taxa_pcoa[
  ,
  species_keep,
  drop = FALSE
]

cat(
  "After cleanup: ",
  nrow(taxa_pcoa),
  " samples x ",
  ncol(taxa_pcoa),
  " species\n",
  sep = ""
)

cat(
  "\nCountry counts after Bray-Curtis cleanup:\n"
)

print(
  table(
    meta_pcoa$Country
  )
)

bray_dist <- vegdist(
  taxa_pcoa,
  method = "bray"
)

stopifnot(
  all(
    is.finite(
      bray_dist
    )
  )
)

cat(
  "Bray-Curtis distance matrix is clean (all finite values).\n"
)


# 10.4 Countries retained for colouring
country_counts_bray <- meta_pcoa %>%
  count(
    Country,
    name = "n_country"
  )

print(
  country_counts_bray %>%
    arrange(
      desc(n_country)
    )
)

countries_keep_bray <- country_counts_bray %>%
  filter(
    n_country >= 30
  ) %>%
  pull(
    Country
  ) %>%
  as.character()

cat(
  "Countries retained for colouring: ",
  paste(
    countries_keep_bray,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

# 10.5 PERMANOVA
set.seed(123)
permanova_country_bray <- adonis2(
  bray_dist ~ Country,
  data = meta_pcoa,
  permutations = 999
)

print(
  permanova_country_bray
)

r2_val_bray <- round(
  permanova_country_bray$R2[1],
  3
)

p_val_bray <- permanova_country_bray$`Pr(>F)`[1]

pval_fmt_bray <- ifelse(
  p_val_bray < 0.001,
  "< 0.001",
  sprintf(
    "%.3f",
    p_val_bray
  )
)

cat(
  "PERMANOVA - Bray-Curtis: R2 = ",
  r2_val_bray,
  ", p = ",
  pval_fmt_bray,
  "\n",
  sep = ""
)

write.table(
  as.data.frame(
    permanova_country_bray
  ),
  "output/Fig4/Fig4D_PERMANOVA_BrayCurtis.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# 10.6 UMAP on Bray-Curtis distances
library(uwot)

n_samp <- attr(
  bray_dist,
  "Size"
)

n_neighbors_v <- round(
  sqrt(
    n_samp
  )
)

cat(
  "Samples for UMAP: ",
  n_samp,
  "\n",
  "n_neighbors: ",
  n_neighbors_v,
  "\n",
  sep = ""
)

set.seed(123)
umap_out_bray <- umap2(
  X = bray_dist,
  n_neighbors = n_neighbors_v,
  n_components = 2,
  min_dist = 0.3,
  init = "spectral",
  verbose = TRUE,
  n_threads = 4
)

cat(
  "UMAP complete. Output dimensions: ",
  paste(
    dim(umap_out_bray),
    collapse = " x "
  ),
  "\n",
  sep = ""
)

# 10.7 Build UMAP plotting data
umap_df_bray <- tibble(
  SampleID = rownames(
    as.matrix(
      bray_dist
    )
  ),
  UMAP1 = umap_out_bray[, 1],
  UMAP2 = umap_out_bray[, 2]
) %>%
  left_join(
    meta_pcoa %>%
      rownames_to_column(
        "SampleID"
      ),
    by = "SampleID"
  ) %>%
  mutate(
    color_id = ifelse(
      as.character(Country) %in% countries_keep_bray,
      as.character(Country),
      "<_40_samples"
    ),
    color_id = factor(
      color_id,
      levels = c(
        countries_keep_bray,
        "<_40_samples"
      )
    )
  )

write.table(
  umap_df_bray %>%
    dplyr::select(
      SampleID,
      UMAP1,
      UMAP2,
      Country,
      color_id
    ),
  "output/Fig4/Fig4D_UMAP_coordinates.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Same country palette used in the revised analysis
my_col_individual <- c(
  "Spain"        = "firebrick",
  "France"       = "cadetblue",
  "Ireland"      = "darkslategray",
  "Austria"      = "sandybrown",
  "USA"          = "coral3",
  "Italy"        = "cornflowerblue",
  "Canada"       = "forestgreen",
  "<_40_samples" = "grey80"
)

# Base UMAP
fig_umap_country_bray <- ggplot(
  umap_df_bray,
  aes(
    x = UMAP1,
    y = UMAP2,
    colour = color_id
  )
) +

  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +

  scale_color_manual(
    values = my_col_individual,
    name = "Country"
  ) +

  scale_fill_manual(
    values = my_col_individual,
    name = "Country"
  ) +

  coord_fixed() +

  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +

  theme_bw(
    base_size = 11
  ) +

  theme(
    panel.grid = element_line(
      colour = "grey94",
      linewidth = 0.3
    ),
    legend.position = "right"
  )

# 10.8 Weighted-average species scores
wa_scores_umap <- wascores(
  umap_out_bray[, 1:2],
  taxa_pcoa
) %>%
  as_tibble(
    rownames = "species"
  ) %>%
  rename(
    UMAP1 = V1,
    UMAP2 = V2
  )

species_of_interest <- c(
  "Lactococcus_lactis",
  "Lactococcus_cremoris",
  "Streptococcus_thermophilus",
  "Lactobacillus_delbrueckii",
  "Lactobacillus_helveticus",
  "Lacticaseibacillus_paracasei",
  "Lacticaseibacillus_rhamnosus",
  "Escherichia_coli"
)

wa_subset_umap <- wa_scores_umap %>%
  filter(
    species %in% species_of_interest
  )

cat(
  "Species found and labelled: ",
  nrow(wa_subset_umap),
  " of ",
  length(species_of_interest),
  "\n",
  sep = ""
)

if (nrow(wa_subset_umap) < length(species_of_interest)) {

  cat(
    "Not found: ",
    paste(
      setdiff(
        species_of_interest,
        wa_subset_umap$species
      ),
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

write.table(
  wa_subset_umap,
  "output/Fig4/Fig4D_weighted_species_scores.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# 10.9 Add species labels and save Fig. 4D
fig_umap_country_species_bray <- fig_umap_country_bray +

  geom_text(
    data = wa_subset_umap,
    aes(
      x = UMAP1,
      y = UMAP2,
      label = species
    ),
    inherit.aes = FALSE,
    colour = "black",
    fontface = "italic",
    size = 3.5,
    check_overlap = TRUE
  )

fig_umap_country_species_bray
ggsave(
  filename = "output/Fig4/Fig4D_UMAP_Country_Core_BrayCurtis.svg",
  plot = fig_umap_country_species_bray,
  width = 5,
  height = 5,
  device = svg
)

# 11. Fig. 4E

# 11.1 Libraries used in the revised ML analysis
library(caret)
library(pROC)
library(ranger)
library(janitor)
library(MASS)
library(ggsci)

metadata <- read.csv(
  "data/MetaCheeseDB_metadata.tsv",
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "NaN",
    "Unknown",
    "unknown",
    "not available"
  )
)

taxa_table <- read.csv(
  "data/MetaCheeseDB_taxonomy.tsv",
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

common_ids <- intersect(
  rownames(metadata),
  rownames(taxa_table)
)

metadata <- metadata[
  common_ids,
  ,
  drop = FALSE
]

taxa_table <- taxa_table[
  common_ids,
  ,
  drop = FALSE
]

taxa_table <- taxa_table[
  rownames(metadata),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(metadata),
    rownames(taxa_table)
  )
)

target_var     <- "Country"
positive_class <- "Italy"
negative_class <- "Spain"

meta0 <- metadata %>%
  filter(
    Part_of_cheese == "Core",
    Country %in% c(
      positive_class,
      negative_class
    )
  ) %>%
  mutate(
    y = factor(
      if_else(
        .data[[target_var]] == positive_class,
        positive_class,
        negative_class
      ),
      levels = c(
        positive_class,
        negative_class
      )
    )
  )

taxa0 <- taxa_table[
  rownames(meta0),
  ,
  drop = FALSE
]

cat(
  "\nSamples used in Italy vs Spain Core analysis:\n"
)

print(
  table(
    meta0$y
  )
)


preferred_group_var <- "Dataset"
meta0$cv_group <- interaction(
  meta0$Dataset,
  meta0$Country,
  drop = TRUE
)

if (
  !is.na(preferred_group_var) &&
  preferred_group_var %in% names(meta0)
) {

  meta0$cv_group <- as.character(
    meta0[[preferred_group_var]]
  )

  meta0$cv_group[
    is.na(meta0$cv_group) |
      meta0$cv_group == ""
  ] <- rownames(meta0)[
    is.na(meta0$cv_group) |
      meta0$cv_group == ""
  ]

} else {

  meta0$cv_group <- rownames(meta0)
}

cat(
  "\nCross-validation grouping variable used:\n"
)

cat(
  ifelse(
    is.na(preferred_group_var),
    "Each sample used as its own group",
    preferred_group_var
  ),
  "\n"
)

cat(
  "Number of CV groups:",
  n_distinct(
    meta0$cv_group
  ),
  "\n"
)

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

confounder_vars <- intersect(
  confounder_vars,
  names(meta0)
)

to_numeric_matrix <- function(x) {

  x <- as.data.frame(x)

  x[] <- lapply(
    x,
    function(z) {
      as.numeric(
        as.character(z)
      )
    }
  )

  m <- as.matrix(x)

  m[
    is.na(m)
  ] <- 0

  storage.mode(m) <- "numeric"

  m
}


relative_abundance <- function(x) {

  m <- to_numeric_matrix(x)

  rs <- rowSums(
    m,
    na.rm = TRUE
  )

  keep <- rs > 0

  m[
    keep,
  ] <- sweep(
    m[
      keep,
      ,
      drop = FALSE
    ],
    1,
    rs[keep],
    "/"
  )

  m[
    !is.finite(m)
  ] <- 0

  m
}


make_group_folds <- function(
    y,
    group,
    k = 5,
    seed = 123) {

  set.seed(seed)

  y <- droplevels(y)

  group <- as.character(group)

  group_y <- tapply(
    as.character(y),
    group,
    function(z) {
      names(
        sort(
          table(z),
          decreasing = TRUE
        )
      )[1]
    }
  )

  group_y <- factor(
    group_y,
    levels = levels(y)
  )

  class_counts <- table(
    group_y
  )

  k_eff <- min(
    k,
    length(group_y),
    min(class_counts)
  )

  if (k_eff < 2) {

    stop(
      "Not enough independent groups per class for Leakage-aware."
    )
  }

  fold_groups <- caret::createFolds(
    group_y,
    k = k_eff,
    list = TRUE,
    returnTrain = FALSE
  )

  lapply(
    fold_groups,
    function(gidx) {

      which(
        group %in% names(group_y)[gidx]
      )
    }
  )
}


prep_taxa_train <- function(
    x_train,
    min_prevalence = 0.05) {

  x_rel <- relative_abundance(
    x_train
  )

  prevalence <- colMeans(
    x_rel > 0
  )

  keep <- prevalence >= min_prevalence &
    is.finite(prevalence)

  if (sum(keep) < 2) {

    keep <- rank(
      -prevalence,
      ties.method = "first"
    ) <= min(
      100,
      ncol(x_rel)
    )
  }

  x_rel <- x_rel[
    ,
    keep,
    drop = FALSE
  ]

  pseudo <- suppressWarnings(
    min(
      x_rel[
        x_rel > 0
      ],
      na.rm = TRUE
    ) / 2
  )

  if (
    !is.finite(pseudo) ||
    pseudo <= 0
  ) {

    pseudo <- 1e-6
  }

  x_clr <- log(
    x_rel + pseudo
  )

  x_clr <- sweep(
    x_clr,
    1,
    rowMeans(x_clr),
    "-"
  )

  original_features <- colnames(
    x_clr
  )

  safe_features <- make.names(
    original_features,
    unique = TRUE
  )

  colnames(
    x_clr
  ) <- safe_features

  list(
    features_original = original_features,
    features_safe = safe_features,
    pseudo = pseudo,
    train = as.data.frame(x_clr)
  )
}


prep_taxa_test <- function(
    x_test,
    prep) {

  x_rel <- relative_abundance(
    x_test
  )

  missing_features <- setdiff(
    prep$features_original,
    colnames(x_rel)
  )

  if (
    length(missing_features) > 0
  ) {

    stop(
      "Some training features are absent from the test table."
    )
  }

  x_rel <- x_rel[
    ,
    prep$features_original,
    drop = FALSE
  ]

  x_clr <- log(
    x_rel + prep$pseudo
  )

  x_clr <- sweep(
    x_clr,
    1,
    rowMeans(x_clr),
    "-"
  )

  colnames(
    x_clr
  ) <- prep$features_safe

  as.data.frame(
    x_clr
  )
}


prep_covariate_mm <- function(
    meta_train,
    meta_test,
    covariates) {

  covariates <- intersect(
    covariates,
    names(meta_train)
  )

  if (
    length(covariates) == 0
  ) {

    return(
      list(
        train = matrix(
          nrow = nrow(meta_train),
          ncol = 0
        ),
        test = matrix(
          nrow = nrow(meta_test),
          ncol = 0
        ),
        covariates_used = character()
      )
    )
  }

  tr <- meta_train[
    ,
    covariates,
    drop = FALSE
  ]

  te <- meta_test[
    ,
    covariates,
    drop = FALSE
  ]

  use <- character()

  for (v in covariates) {

    if (
      is.numeric(
        tr[[v]]
      )
    ) {

      med <- median(
        tr[[v]],
        na.rm = TRUE
      )

      if (
        !is.finite(med)
      ) {

        med <- 0
      }

      tr[[v]][
        is.na(
          tr[[v]]
        )
      ] <- med

      te[[v]][
        is.na(
          te[[v]]
        )
      ] <- med

      if (
        sd(
          tr[[v]],
          na.rm = TRUE
        ) > 0
      ) {

        use <- c(
          use,
          v
        )
      }

    } else {

      tr[[v]] <- as.character(
        tr[[v]]
      )

      te[[v]] <- as.character(
        te[[v]]
      )

      tr[[v]][
        is.na(tr[[v]]) |
          tr[[v]] == "" |
          tr[[v]] == "NA"
      ] <- "Missing"

      te[[v]][
        is.na(te[[v]]) |
          te[[v]] == "" |
          te[[v]] == "NA"
      ] <- "Missing"

      train_levels <- sort(
        unique(
          tr[[v]]
        )
      )

      te[[v]][
        !te[[v]] %in% train_levels
      ] <- "Other"

      all_levels <- unique(
        c(
          train_levels,
          "Other"
        )
      )

      tr[[v]] <- factor(
        tr[[v]],
        levels = all_levels
      )

      te[[v]] <- factor(
        te[[v]],
        levels = all_levels
      )

      if (
        n_distinct(
          tr[[v]]
        ) > 1
      ) {

        use <- c(
          use,
          v
        )
      }
    }
  }

  if (
    length(use) == 0
  ) {

    return(
      list(
        train = matrix(
          nrow = nrow(meta_train),
          ncol = 0
        ),
        test = matrix(
          nrow = nrow(meta_test),
          ncol = 0
        ),
        covariates_used = character()
      )
    )
  }

  form <- as.formula(
    paste(
      "~",
      paste(
        sprintf(
          "`%s`",
          use
        ),
        collapse = " + "
      )
    )
  )

  mm_train <- model.matrix(
    form,
    data = tr
  )[
    ,
    -1,
    drop = FALSE
  ]

  mm_test <- model.matrix(
    form,
    data = te
  )[
    ,
    -1,
    drop = FALSE
  ]

  all_cols <- union(
    colnames(mm_train),
    colnames(mm_test)
  )

  add_missing_cols <- function(
      m,
      cols) {

    missing <- setdiff(
      cols,
      colnames(m)
    )

    if (
      length(missing) > 0
    ) {

      add <- matrix(
        0,
        nrow = nrow(m),
        ncol = length(missing)
      )

      colnames(
        add
      ) <- missing

      m <- cbind(
        m,
        add
      )
    }

    m[
      ,
      cols,
      drop = FALSE
    ]
  }

  mm_train <- add_missing_cols(
    mm_train,
    all_cols
  )

  mm_test <- add_missing_cols(
    mm_test,
    all_cols
  )

  list(
    train = mm_train,
    test = mm_test,
    covariates_used = use
  )
}


residualize_train_test <- function(
    x_train,
    x_test,
    meta_train,
    meta_test,
    covariates) {

  mm <- prep_covariate_mm(
    meta_train,
    meta_test,
    covariates
  )

  if (
    ncol(mm$train) == 0
  ) {

    return(
      list(
        train = x_train,
        test = x_test,
        covariates_used = character()
      )
    )
  }

  Xtr <- cbind(
    Intercept = 1,
    mm$train
  )

  Xte <- cbind(
    Intercept = 1,
    mm$test
  )

  Ytr <- as.matrix(
    x_train
  )

  Yte <- as.matrix(
    x_test
  )

  lambda <- 1e-6

  beta <- solve(
    crossprod(Xtr) +
      diag(
        lambda,
        ncol(Xtr)
      ),
    crossprod(
      Xtr,
      Ytr
    )
  )

  res_train <- Ytr -
    Xtr %*% beta

  res_test <- Yte -
    Xte %*% beta

  colnames(
    res_train
  ) <- colnames(
    x_train
  )

  colnames(
    res_test
  ) <- colnames(
    x_test
  )

  list(
    train = as.data.frame(
      res_train
    ),
    test = as.data.frame(
      res_test
    ),
    covariates_used = mm$covariates_used
  )
}

summarise_binary <- function(
    preds,
    positive = positive_class,
    negative = negative_class) {

  preds$obs <- factor(
    preds$obs,
    levels = c(
      positive,
      negative
    )
  )

  roc_obj <- pROC::roc(
    response = preds$obs,
    predictor = preds[[positive]],
    levels = c(
      negative,
      positive
    ),
    quiet = TRUE
  )

  pred_class <- factor(
    if_else(
      preds[[positive]] >= 0.5,
      positive,
      negative
    ),
    levels = c(
      positive,
      negative
    )
  )

  cm <- confusionMatrix(
    pred_class,
    preds$obs,
    positive = positive
  )

  tibble(
    n = nrow(preds),
    n_folds = n_distinct(
      preds$fold
    ),
    AUC = as.numeric(
      pROC::auc(
        roc_obj
      )
    ),
    AUC_low = as.numeric(
      pROC::ci.auc(
        roc_obj
      )[1]
    ),
    AUC_high = as.numeric(
      pROC::ci.auc(
        roc_obj
      )[3]
    ),
    Accuracy = unname(
      cm$overall[
        "Accuracy"
      ]
    ),
    Sensitivity = unname(
      cm$byClass[
        "Sensitivity"
      ]
    ),
    Specificity = unname(
      cm$byClass[
        "Specificity"
      ]
    ),
    BalancedAccuracy = unname(
      cm$byClass[
        "Balanced Accuracy"
      ]
    )
  )
}

run_nested_grouped_rf <- function(
    taxa,
    meta,
    outcome = "y",
    group = "cv_group",
    confounders = NULL,
    residualize = FALSE,
    k_outer = 5,
    k_inner = 3,
    seed = 123) {

  y <- droplevels(
    meta[[outcome]]
  )

  outer_test_folds <- make_group_folds(
    y,
    meta[[group]],
    k = k_outer,
    seed = seed
  )

  all_preds <- vector(
    "list",
    length(
      outer_test_folds
    )
  )

  for (
    i in seq_along(
      outer_test_folds
    )
  ) {

    test_id <- outer_test_folds[[i]]

    train_id <- setdiff(
      seq_along(y),
      test_id
    )

    overlap <- intersect(
      meta[[group]][train_id],
      meta[[group]][test_id]
    )

    stopifnot(
      length(overlap) == 0
    )

    taxa_prep <- prep_taxa_train(
      taxa[
        train_id,
        ,
        drop = FALSE
      ]
    )

    x_train <- taxa_prep$train

    x_test <- prep_taxa_test(
      taxa[
        test_id,
        ,
        drop = FALSE
      ],
      taxa_prep
    )

    if (
      residualize
    ) {

      rz <- residualize_train_test(
        x_train = x_train,
        x_test = x_test,
        meta_train = meta[
          train_id,
          ,
          drop = FALSE
        ],
        meta_test = meta[
          test_id,
          ,
          drop = FALSE
        ],
        covariates = confounders
      )

      x_train <- rz$train

      x_test <- rz$test
    }

    y_train <- droplevels(
      y[train_id]
    )

    group_train <- meta[[group]][
      train_id
    ]

    inner_test_folds <- tryCatch(
      make_group_folds(
        y_train,
        group_train,
        k = k_inner,
        seed = seed + i
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      is.null(
        inner_test_folds
      )
    ) {

      ctrl <- trainControl(
        method = "none",
        classProbs = TRUE,
        summaryFunction = twoClassSummary
      )

      tune_grid <- expand.grid(
        mtry = max(
          1,
          floor(
            sqrt(
              ncol(
                x_train
              )
            )
          )
        ),
        splitrule = "gini",
        min.node.size = 5
      )

    } else {

      inner_index <- lapply(
        inner_test_folds,
        function(te) {

          setdiff(
            seq_along(y_train),
            te
          )
        }
      )

      ctrl <- trainControl(
        method = "cv",
        index = inner_index,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "final",
        allowParallel = TRUE
      )

      mtry_values <- unique(
        pmax(
          1,
          floor(
            c(
              sqrt(
                ncol(
                  x_train
                )
              ),
              ncol(x_train) / 10,
              ncol(x_train) / 3
            )
          )
        )
      )

      tune_grid <- expand.grid(
        mtry = mtry_values,
        splitrule = "gini",
        min.node.size = c(
          1,
          5,
          10
        )
      )
    }

    set.seed(
      seed + 1000 + i
    )

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

    prob <- predict(
      fit,
      x_test,
      type = "prob"
    )

    all_preds[[i]] <- bind_cols(
      tibble(
        sample_id = rownames(meta)[
          test_id
        ],
        obs = y[
          test_id
        ],
        fold = i,
        cv_group = as.character(
          meta[[group]][
            test_id
          ]
        )
      ),
      as_tibble(
        prob
      )
    )
  }

  bind_rows(
    all_preds
  )
}


run_nested_grouped_metadata_rf <- function(
    meta,
    covariates,
    outcome = "y",
    group = "cv_group",
    k_outer = 5,
    k_inner = 3,
    seed = 123) {

  y <- droplevels(
    meta[[outcome]]
  )

  outer_test_folds <- make_group_folds(
    y,
    meta[[group]],
    k = k_outer,
    seed = seed
  )

  all_preds <- vector(
    "list",
    length(
      outer_test_folds
    )
  )

  for (
    i in seq_along(
      outer_test_folds
    )
  ) {

    test_id <- outer_test_folds[[i]]

    train_id <- setdiff(
      seq_along(y),
      test_id
    )

    mm <- prep_covariate_mm(
      meta_train = meta[
        train_id,
        ,
        drop = FALSE
      ],
      meta_test = meta[
        test_id,
        ,
        drop = FALSE
      ],
      covariates = covariates
    )

    if (
      ncol(
        mm$train
      ) < 1
    ) {

      stop(
        "No usable metadata covariates available for Confounder-only."
      )
    }

    x_train <- as.data.frame(
      mm$train
    )

    x_test <- as.data.frame(
      mm$test
    )

    y_train <- droplevels(
      y[train_id]
    )

    group_train <- meta[[group]][
      train_id
    ]

    inner_test_folds <- tryCatch(
      make_group_folds(
        y_train,
        group_train,
        k = k_inner,
        seed = seed + i
      ),
      error = function(e) {
        NULL
      }
    )

    if (
      is.null(
        inner_test_folds
      )
    ) {

      ctrl <- trainControl(
        method = "none",
        classProbs = TRUE,
        summaryFunction = twoClassSummary
      )

      tune_grid <- expand.grid(
        mtry = max(
          1,
          floor(
            sqrt(
              ncol(
                x_train
              )
            )
          )
        ),
        splitrule = "gini",
        min.node.size = 5
      )

    } else {

      inner_index <- lapply(
        inner_test_folds,
        function(te) {

          setdiff(
            seq_along(y_train),
            te
          )
        }
      )

      ctrl <- trainControl(
        method = "cv",
        index = inner_index,
        classProbs = TRUE,
        summaryFunction = twoClassSummary,
        savePredictions = "final",
        allowParallel = TRUE
      )

      tune_grid <- expand.grid(
        mtry = unique(
          pmax(
            1,
            floor(
              c(
                sqrt(
                  ncol(
                    x_train
                  )
                ),
                ncol(x_train) / 2
              )
            )
          )
        ),
        splitrule = "gini",
        min.node.size = c(
          1,
          5,
          10
        )
      )
    }

    set.seed(
      seed + 2000 + i
    )

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

    prob <- predict(
      fit,
      x_test,
      type = "prob"
    )

    all_preds[[i]] <- bind_cols(
      tibble(
        sample_id = rownames(meta)[
          test_id
        ],
        obs = y[
          test_id
        ],
        fold = i,
        cv_group = as.character(
          meta[[group]][
            test_id
          ]
        )
      ),
      as_tibble(
        prob
      )
    )
  }

  bind_rows(
    all_preds
  )
}

# 11.7 Run the four sensitivity analyses
# A. Optimistic: sample-level CV, retained as comparator.
meta_sample_cv <- meta0

meta_sample_cv$cv_group <- rownames(
  meta_sample_cv
)


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


summary_sample_cv <- summarise_binary(
  pred_sample_cv
) %>%
  mutate(
    analysis = "A_naive_sample_level_CV"
  )


# B. Dataset-grouped CV.
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


summary_grouped_cv <- summarise_binary(
  pred_grouped_cv
) %>%
  mutate(
    analysis = "B_grouped_CV_no_covariate_adjustment"
  )


# C. Metadata-only confounder model.
pred_metadata_only <- run_nested_grouped_metadata_rf(
  meta = meta0,
  covariates = confounder_vars,
  outcome = "y",
  group = "cv_group",
  k_outer = 5,
  k_inner = 3,
  seed = 123
)


summary_metadata_only <- summarise_binary(
  pred_metadata_only
) %>%
  mutate(
    analysis = "C_metadata_only_confounder_model"
  )


# D. Dataset-grouped, covariate-residualised microbiome model.
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


summary_residualised <- summarise_binary(
  pred_residualised
) %>%
  mutate(
    analysis = "D_grouped_CV_covariate_residualised_microbiome"
  )


# 11.8 Save held-out predictions and performance summary
all_predictions <- bind_rows(
  pred_sample_cv %>%
    mutate(
      analysis = "A_naive_sample_level_CV"
    ),

  pred_grouped_cv %>%
    mutate(
      analysis = "B_grouped_CV_no_covariate_adjustment"
    ),

  pred_metadata_only %>%
    mutate(
      analysis = "C_metadata_only_confounder_model"
    ),

  pred_residualised %>%
    mutate(
      analysis = "D_grouped_CV_covariate_residualised_microbiome"
    )
)


all_summaries <- bind_rows(
  summary_sample_cv,
  summary_grouped_cv,
  summary_metadata_only,
  summary_residualised
) %>%
  dplyr::select(
    analysis,
    dplyr::everything()
  )


cat(
  "\nFig. 4E model-performance summary:\n"
)


print(
  all_summaries
)


write_tsv(
  all_predictions,
  "output/Fig4/Fig4E_model_predictions.tsv"
)


write_tsv(
  all_summaries,
  "output/Fig4/Fig4E_model_performance.tsv"
)


# 11.9 Construct the four ROC curves
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
  dplyr::group_by(
    analysis_label
  ) %>%
  dplyr::group_modify(
    ~ {

      roc_obj <- pROC::roc(
        response = .x$obs,
        predictor = .x[[positive_class]],
        levels = c(
          negative_class,
          positive_class
        ),
        direction = "<",
        quiet = TRUE
      )

      tibble::tibble(
        specificity = roc_obj$specificities,
        sensitivity = roc_obj$sensitivities,
        AUC = as.numeric(
          pROC::auc(
            roc_obj
          )
        )
      )
    }
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    auc_label = paste0(
      analysis_label,
      " (AUC=",
      sprintf(
        "%.3f",
        AUC
      ),
      ")"
    )
  )


write_tsv(
  roc_df,
  "output/Fig4/Fig4E_ROC_coordinates.tsv"
)

palette_okabe <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#D55E00"
)


p_roc_multiple <- ggplot(
  roc_df,
  aes(
    x = 1 - specificity,
    y = sensitivity,
    colour = auc_label,
    group = auc_label
  )
) +

  geom_path(
    linewidth = 1.2
  ) +

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

  theme_bw(
    base_size = 10
  ) +

  theme(
    legend.position = "right",
    legend.text = element_text(
      size = 10
    )
  )


p_roc_multiple


# 11.10 Save Fig. 4E
svglite(
  file = "output/Fig4/Fig4E_multiroc.svg",
  width = 7.5,
  height = 7.5
)


print(
  p_roc_multiple
)


dev.off()


# 12. Fig. 4F-G
# 12.1 Libraries
library(kernelshap)
library(shapviz)

# 12.2 Prepare the full Italy-vs-Spain microbiome feature matrix

taxa_prep_all <- prep_taxa_train(
  taxa0,
  min_prevalence = 0.05
)


x_all <- taxa_prep_all$train

y_all <- meta0$y


cat(
  "\nFig. 4F-G SHAP input:\n",
  nrow(x_all),
  " samples x ",
  ncol(x_all),
  " microbial features\n",
  sep = ""
)


cat(
  "Class counts:\n"
)


print(
  table(
    y_all
  )
)

feature_map <- tibble(
  feature_safe = taxa_prep_all$features_safe,
  feature_original = taxa_prep_all$features_original
)


write_tsv(
  feature_map,
  "output/Fig4/Fig4FG_feature_map.tsv"
)



# 12.3 Fit the final RF model used for SHAP interpretation

outer_test_folds_shap <- make_group_folds(
  y = y_all,
  group = meta0$cv_group,
  k = 5,
  seed = 123
)


outer_train_index_shap <- lapply(
  outer_test_folds_shap,
  function(test_id) {

    setdiff(
      seq_along(y_all),
      test_id
    )
  }
)


rf_ctrl_shap <- trainControl(
  method = "cv",
  index = outer_train_index_shap,
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
  trControl = rf_ctrl_shap,

  tuneGrid = expand.grid(
    mtry = unique(
      pmax(
        1,
        floor(
          c(
            sqrt(
              ncol(x_all)
            ),
            ncol(x_all) / 10,
            ncol(x_all) / 3
          )
        )
      )
    ),
    splitrule = "gini",
    min.node.size = c(
      1,
      5,
      10
    )
  ),

  num.trees = 1000,
  importance = "permutation"
)


cat(
  "\nFinal RF tuning parameters used for SHAP:\n"
)


print(
  rf_microbiome_final$bestTune
)


write_tsv(
  as_tibble(
    rf_microbiome_final$results
  ),
  "output/Fig4/Fig4FG_RF_tuning_results.tsv"
)


# 12.4 Compute kernel SHAP
pred_fun_prob <- function(
    object,
    newdata) {

  as.matrix(
    predict(
      object,
      newdata = as.data.frame(
        newdata
      ),
      type = "prob"
    )
  )
}

set.seed(123)


bg_size <- min(
  100,
  nrow(x_all)
)


bg_index <- unlist(
  tapply(
    seq_len(
      nrow(x_all)
    ),
    y_all,
    function(idx) {

      sample(
        idx,
        size = min(
          length(idx),
          ceiling(
            bg_size / 2
          )
        )
      )
    }
  ),
  use.names = FALSE
)


bg_index <- unique(
  bg_index
)


bg_data <- x_all[
  bg_index,
  ,
  drop = FALSE
]


X_explain <- x_all


cat(
  "\nBackground distribution size: ",
  nrow(bg_data),
  "\n",
  "Samples explained: ",
  nrow(X_explain),
  "\n",
  "Features explained: ",
  ncol(X_explain),
  "\n",
  sep = ""
)


ks_IT_ES <- kernelshap(
  object = rf_microbiome_final,
  X = X_explain,
  bg_X = bg_data,
  pred_fun = pred_fun_prob,
  verbose = TRUE
)


sv_Country <- shapviz(
  ks_IT_ES
)


stopifnot(
  all(
    c(
      "Italy",
      "Spain"
    ) %in% names(
      sv_Country
    )
  )
)


# 12.6 Fig. 4F | Italy SHAP beeswarm

p_bee_Italy <- sv_importance(
  sv_Country$Italy,
  kind = "bee"
)

ggsave(
  filename = "output/Fig4/Fig4F_SHAP_Italy.svg",
  plot = p_bee_Italy,
  width = 17.5,
  height = 15,
  units = "cm",
  dpi = 600
)


# 12.7 Fig. 4G | Spain SHAP beeswarm

p_bee_Spain <- sv_importance(
  sv_Country$Spain,
  kind = "bee"
)

ggsave(
  filename = "output/Fig4/Fig4G_SHAP_Spain.svg",
  plot = p_bee_Spain,
  width = 17.5,
  height = 15,
  units = "cm",
  dpi = 600
)
