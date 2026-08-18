# Fig. 2


# 0. Packages
library(tidyverse)
library(vegan)
library(rcompanion)
library(patchwork)

set.seed(123)
# 1. Load data
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
  identical(rownames(metadata), rownames(taxonomy))
)


# Harmonised metadata used in the analysis

metadata_vars <- c(
  "Country",
  "Part_of_cheese",
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


pretty_name <- function(x) {
  
  x <- str_replace_all(x, "_", " ")
  
  x <- str_replace(
    x,
    "^Rheological properties$",
    "Texture"
  )
  
  x
}


# 2. PANEL A
# 2.1 Prepare taxonomic profiles
taxonomy_permanova <- taxonomy

taxonomy_permanova[is.na(taxonomy_permanova)] <- 0

taxonomy_permanova <- taxonomy_permanova[
  rowSums(taxonomy_permanova) > 0,
  ,
  drop = FALSE
]

metadata_permanova <- metadata[
  rownames(taxonomy_permanova),
  ,
  drop = FALSE
]

# 2.2 Bray-Curtis distances
bray <- vegdist(
  taxonomy_permanova,
  method = "bray"
)

# 2.3 Multivariable PERMANOVA
set.seed(123)
permanova <- adonis2(
  bray ~ Country +
    Part_of_cheese +
    Milk_source +
    Animal_feeding +
    Milk_processing +
    Thermisation +
    Pasteurization +
    Skimming +
    Inoculation_of_moulds +
    Curd_cutting +
    Presence_of_rind +
    Backslopping +
    Starter_culture +
    Technology_of_production +
    Coagulation_method +
    Ripening_period +
    Rheological_properties +
    Temperature_of_curd_processing_Celsius_degree,
  data = metadata_permanova,
  by = "margin",
  permutations = 999,
  na.action = na.omit
)

# 2.4 Prepare PERMANOVA results
permanova_df <- as.data.frame(permanova) %>%
  rownames_to_column("variable") %>%
  filter(
    !variable %in% c("Residual", "Total")
  ) %>%
  mutate(
    R2_pct = round(R2 * 100, 2),
    p_value = `Pr(>F)`,
    
    status = case_when(
      is.na(p_value) ~ "NA",
      p_value <= 0.05 ~ "Significant (p≤0.05)",
      TRUE ~ "Not significant"
    ),
    
    status = factor(
      status,
      levels = c(
        "Significant (p≤0.05)",
        "Not significant",
        "NA"
      )
    ),
    
    variable = pretty_name(variable),
    
    variable = fct_reorder(
      variable,
      R2
    )
  )

# 2.5 Panel A plot
fig_A <- ggplot(
  permanova_df,
  aes(
    x = R2_pct,
    y = variable,
    colour = status
  )
) +
  
  geom_segment(
    aes(
      x = 0,
      xend = R2_pct,
      yend = variable
    ),
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 3.5
  ) +
  
  geom_text(
    aes(
      label = paste0(
        format(
          R2_pct,
          trim = TRUE
        ),
        "%"
      )
    ),
    hjust = -0.3,
    size = 3.5,
    colour = "grey30"
  ) +
  
  scale_colour_manual(
    values = c(
      "Significant (p≤0.05)" = "#254E70",
      "Not significant" = "#B4B2A9",
      "NA" = "grey50"
    ),
    breaks = c(
      "Significant (p≤0.05)",
      "Not significant",
      "NA"
    ),
    name = NULL
  ) +
  
  scale_x_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(
      mult = c(0, 0.18)
    )
  ) +
  
  labs(
    x = "R² pct",
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      colour = "grey92",
      linewidth = 0.4
    ),
    axis.text = element_text(
      colour = "grey30"
    ),
    legend.position = "right"
  )

# 3. PANEL B
# 3.1 Prepare metadata
metadata_cv <- metadata %>%
  mutate(
    across(
      all_of(metadata_vars),
      ~ na_if(
        na_if(
          as.character(.x),
          ""
        ),
        "NA"
      )
    )
  )

pairs <- as.data.frame(
  t(
    combn(
      metadata_vars,
      2
    )
  ),
  stringsAsFactors = FALSE
)

colnames(pairs) <- c(
  "var1",
  "var2"
)

# 3.2 Calculate pairwise Cramér's V and significance
set.seed(123)
covariation <- map_dfr(
  seq_len(nrow(pairs)),
  function(i) {
    
    var1 <- pairs$var1[i]
    var2 <- pairs$var2[i]
    
    dat <- metadata_cv %>%
      select(
        all_of(
          c(
            var1,
            var2
          )
        )
      ) %>%
      drop_na()
    
    if (
      nrow(dat) < 10 ||
      n_distinct(dat[[var1]]) < 2 ||
      n_distinct(dat[[var2]]) < 2
    ) {
      
      return(
        tibble(
          var1 = var1,
          var2 = var2,
          n = nrow(dat),
          cramers_v = NA_real_,
          p_value = NA_real_
        )
      )
    }
    
    tab <- table(
      dat[[var1]],
      dat[[var2]]
    )
    
    cramers_v <- rcompanion::cramerV(
      tab,
      bias.correct = TRUE
    )
    
    chi <- suppressWarnings(
      chisq.test(tab)
    )
    
    if (min(chi$expected) < 5) {
      
      test <- fisher.test(
        tab,
        simulate.p.value = TRUE,
        B = 10000
      )
      
    } else {
      
      test <- chi
    }
    
    tibble(
      var1 = var1,
      var2 = var2,
      n = nrow(dat),
      cramers_v = as.numeric(cramers_v),
      p_value = test$p.value
    )
  }
)

# 3.3 FDR correction
covariation <- covariation %>%
  mutate(
    q_value = p.adjust(
      p_value,
      method = "BH"
    ),
    
    sig = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01 ~ "**",
      q_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# 3.4 Order variables according to mean pairwise Cramér's V
variable_order <- covariation %>%
  filter(
    !is.na(cramers_v)
  ) %>%
  pivot_longer(
    cols = c(
      var1,
      var2
    ),
    values_to = "variable"
  ) %>%
  group_by(variable) %>%
  summarise(
    mean_v = mean(
      cramers_v,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(mean_v)
  ) %>%
  pull(variable)

# 3.5 Construct complete symmetric matrix
mat_long <- covariation %>%
  select(
    var1,
    var2,
    cramers_v,
    sig,
    q_value
  )

mat_mirror <- mat_long[
  ,
  c(
    "var2",
    "var1",
    "cramers_v",
    "sig",
    "q_value"
  )
]

colnames(mat_mirror)[1:2] <- c(
  "var1",
  "var2"
)

mat_diagonal <- tibble(
  var1 = metadata_vars,
  var2 = metadata_vars,
  cramers_v = 1,
  sig = "",
  q_value = 0
)

mat_full <- bind_rows(
  mat_long,
  mat_mirror,
  mat_diagonal
)

mat_full <- mat_full %>%
  mutate(
    var1 = pretty_name(var1),
    var2 = pretty_name(var2),
    
    var1 = factor(
      var1,
      levels = pretty_name(variable_order)
    ),
    
    var2 = factor(
      var2,
      levels = pretty_name(variable_order)
    )
  )


# 3.6 Panel B plot
fig_B <- ggplot(
  mat_full,
  aes(
    x = var1,
    y = var2,
    fill = cramers_v
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.4
  ) +
  
  geom_text(
    aes(
      label = ifelse(
        is.na(cramers_v),
        "",
        paste0(
          sprintf(
            "%.2f",
            cramers_v
          ),
          "\n",
          sig
        )
      ),
      
      colour = ifelse(
        cramers_v > 0.5,
        "white",
        "black"
      )
    ),
    size = 2.8,
    fontface = "bold",
    lineheight = 0.85
  ) +
  
  scale_colour_identity() +
  
  scale_fill_distiller(
    palette = "RdPu",
    direction = 1,
    name = "Cramér's V",
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.5,
      0.7,
      1
    ),
    labels = c(
      "0",
      "0.25",
      "0.5",
      "0.7",
      "1"
    ),
    na.value = "grey92"
  ) +
  
  coord_fixed() +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(
    base_size = 10
  ) +
  
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    panel.grid = element_blank(),
    
    legend.position = "right"
  )


# 4. Combine panels
fig2 <- fig_A / fig_B +
  
  plot_layout(
    heights = c(
      0.75,
      1.25
    )
  ) +
  
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        size = 22,
        face = "bold"
      )
    )
  )


# 5. Save outputs
dir.create(
  "output/Fig2",
  recursive = TRUE,
  showWarnings = FALSE
)

# Final figure
ggsave(
  filename = "output/Fig2/Fig2_metadata_covariation_PERMANOVA.svg",
  plot = fig2,
  width = 11,
  height = 13
)

# PERMANOVA results
write.table(
  permanova_df,
  "output/Fig2/Fig2A_PERMANOVA_marginal_tests.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Pairwise Cramér's V results
write.table(
  covariation,
  "output/Fig2/Fig2B_CramersV_pairwise.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
