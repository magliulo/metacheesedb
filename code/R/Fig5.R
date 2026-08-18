# Fig. 5A-C


# 0. Libraries
library(tidyverse)
library(ggpubr)
library(forcats)
library(ggsankey)
library(viridis)


# 1. Output directory
dir.create(
  "output/Fig5",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Fig. 5A
Fig5A <- read.delim(
  "data/Fig5A_top30_species.tsv",
  check.names = FALSE
)

rheo_cols <- c(
  "Soft"      = "#00bf7dcc",
  "Semihard"  = "#a3a500cc",
  "Hard"      = "#f8766dcc",
  "Very_hard" = "#e76bf3cc"
)

topN_tbl <- Fig5A %>%
  dplyr::distinct(
    species,
    total
  ) %>%
  dplyr::arrange(
    dplyr::desc(total)
  )

plot_df <- Fig5A %>%
  dplyr::mutate(
    species = factor(
      species,
      levels = rev(
        topN_tbl$species
      )
    ),
    
    Rheo = factor(
      Rheo,
      levels = c(
        "Soft",
        "Semihard",
        "Hard",
        "Very_hard"
      )
    )
  )

p_top30_stack <- ggplot(
  plot_df,
  aes(
    x = species,
    y = n,
    fill = Rheo
  )
) +
  geom_col(
    width = 0.8,
    colour = "white"
  ) +
  coord_flip() +
  scale_fill_manual(
    values = rheo_cols,
    name = "Rheological properties"
  ) +
  labs(
    x = NULL,
    y = "Number of MAGs",
    title = ""
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = 20
    ),
    legend.position = "right"
  )

totals <- topN_tbl %>%
  dplyr::mutate(
    species = factor(
      species,
      levels = levels(
        plot_df$species
      )
    )
  )

pad <- max(
  totals$total
) * 0.03

p_top30_stack_labeled <- p_top30_stack +
  geom_text(
    data = totals,
    aes(
      x = species,
      y = total,
      label = scales::comma(
        total
      )
    ),
    inherit.aes = FALSE,
    nudge_y = pad,
    hjust = 0,
    size = 4
  ) +
  expand_limits(
    y = max(
      totals$total
    ) * 1.12
  ) +
  ggpubr::theme_classic2() +
  coord_flip() +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = 20
    ),
    legend.position = "right"
  )

p_top30_stack_labeled

ggsave(
  filename = "output/Fig5/Fig5A_top30_MAGs_rheology_stacked.svg",
  plot = p_top30_stack_labeled,
  device = "svg",
  width = 60,
  height = 30,
  units = "cm",
  dpi = 600
)


# 3. Fig. 5B
Fig5B <- read.delim(
  "data/Fig5B_percentage.tsv",
  check.names = FALSE
) %>%
  dplyr::mutate(
    Rheo = factor(
      Rheo,
      levels = c(
        "Soft",
        "Semihard",
        "Hard",
        "Very_hard"
      )
    ),
    
    SGB = factor(
      SGB,
      levels = c(
        "kSGB",
        "uSGB"
      )
    )
  )


sgb_cols <- c(
  "kSGB" = "#0072B2",
  "uSGB" = "#E69F00"
)


centres <- Fig5B %>%
  dplyr::distinct(
    Rheo,
    total
  )


p_donut <- ggplot(
  Fig5B,
  aes(
    x = 1,
    y = pct,
    fill = SGB
  )
) +
  geom_col(
    width = 0.75,
    colour = "white"
  ) +
  geom_text(
    aes(
      label = lbl
    ),
    position = position_stack(
      vjust = 0.5
    ),
    colour = "white",
    size = 3
  ) +
  coord_polar(
    theta = "y"
  ) +
  xlim(
    0.5,
    1.5
  ) +
  facet_wrap(
    ~ Rheo
  ) +
  scale_fill_manual(
    values = sgb_cols,
    name = "SGB type"
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "right"
  )


p_donut


ggsave(
  filename = "output/Fig5/Fig5B_SGB_unknownness_by_texture.svg",
  plot = p_donut,
  device = "svg",
  height = 10,
  width = 10,
  units = "cm",
  dpi = 600
)



# 4. Fig. 5C
Fig5C <- read.delim(
  "data/Fig5C_sankey.tsv",
  check.names = FALSE
)

phy_order <- Fig5C %>%
  dplyr::count(
    phylum,
    name = "n"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n)
  ) %>%
  dplyr::pull(
    phylum
  )


rheo_order <- Fig5C %>%
  dplyr::count(
    Rheological_properties,
    name = "n",
    sort = TRUE
  ) %>%
  dplyr::pull(
    Rheological_properties
  )


unk_counts <- Fig5C %>%
  dplyr::count(
    unknowness,
    name = "n",
    sort = TRUE
  )


unk_order <- unk_counts %>%
  dplyr::pull(
    unknowness
  )

if (
  "uSGB" %in% unk_order
) {
  
  unk_order <- c(
    setdiff(
      unk_order,
      "uSGB"
    ),
    "uSGB"
  )
}


all_levels <- c(
  phy_order,
  rheo_order,
  unk_order
)


sankey_df <- Fig5C %>%
  ggsankey::make_long(
    phylum,
    Rheological_properties,
    unknowness
  ) %>%
  dplyr::mutate(
    node = ifelse(
      x == "unknowness" &
        node == "uSGB",
      "zzz_uSGB",
      node
    ),
    
    next_node = ifelse(
      next_x == "unknowness" &
        next_node == "uSGB",
      "zzz_uSGB",
      next_node
    ),
    
    x = factor(
      x,
      levels = c(
        "phylum",
        "Rheological_properties",
        "unknowness"
      )
    ),
    
    next_x = factor(
      next_x,
      levels = levels(x)
    ),
    
    node = factor(
      node,
      levels = c(
        phy_order,
        rheo_order,
        c(
          setdiff(
            unk_order,
            "uSGB"
          ),
          "zzz_uSGB"
        )
      )
    ),
    
    next_node = factor(
      next_node,
      levels = levels(node)
    ),
    
    node = forcats::fct_recode(
      node,
      "uSGB" = "zzz_uSGB"
    ),
    
    next_node = forcats::fct_recode(
      next_node,
      "uSGB" = "zzz_uSGB"
    )
  )


p_sankey <- ggplot(
  sankey_df,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = node,
    label = node
  )
) +
  ggsankey::geom_sankey(
    flow.alpha = 0.5,
    node.color = "gray30",
    show.legend = FALSE
  ) +
  ggsankey::geom_sankey_label(
    size = 3.2,
    color = "black",
    fill = "white",
    hjust = 0.5
  ) +
  ggsankey::theme_sankey(
    base_size = 15
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.margin = margin(
      15,
      15,
      15,
      15,
      "pt"
    )
  ) +
  scale_y_sqrt() +
  scale_x_discrete(
    labels = c(
      "phylum" = "Phylum",
      "Rheological_properties" = "Texture",
      "unknowness" = "Unknowness"
    )
  ) +
  labs(
    title = ""
  ) +
  scale_fill_viridis_d(
    option = "plasma",
    alpha = 1,
    begin = 0.1,
    end = 0.9,
    direction = 1
  )


p_sankey


ggsave(
  filename = "output/Fig5/Fig5C_sankey.svg",
  plot = p_sankey,
  device = "svg",
  height = 15,
  width = 20,
  units = "cm",
  dpi = 600
)
