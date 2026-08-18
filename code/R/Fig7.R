# Fig. 7


# 0. Libraries
library(tidyverse)
library(patchwork)
library(viridis)
library(scales)


# 1. Output directory
dir.create(
  "output/Fig7",
  recursive = TRUE,
  showWarnings = FALSE
)


make_fig7_panel <- function(
    bubble_file,
    bars_file,
    feature_col,
    italic_y = FALSE
) {

  


  bubble_df <- read.csv(
    bubble_file,
    check.names = FALSE
  )

  bars_df <- read.csv(
    bars_file,
    check.names = FALSE
  )


  required_bubble <- c(
    feature_col,
    "Rheology",
    "pct"
  )

  required_bars <- c(
    feature_col,
    "order",
    "n_hits"
  )

  stopifnot(
    all(
      required_bubble %in%
        names(bubble_df)
    ),
    all(
      required_bars %in%
        names(bars_df)
    )
  )


  

  feature_order_desc <- bars_df %>%
    dplyr::group_by(
      .data[[feature_col]]
    ) %>%
    dplyr::summarise(
      total = sum(
        n_hits,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(
      dplyr::desc(total)
    ) %>%
    dplyr::pull(
      .data[[feature_col]]
    )


  feature_order_plot <- rev(
    feature_order_desc
  )


  bubble_df[[feature_col]] <- factor(
    bubble_df[[feature_col]],
    levels = feature_order_plot
  )

  bars_df[[feature_col]] <- factor(
    bars_df[[feature_col]],
    levels = feature_order_plot
  )

  bubble_df <- bubble_df %>%
    dplyr::mutate(
      Rheology = factor(
        Rheology,
        levels = c(
          "Soft",
          "Semihard",
          "Hard",
          "Very_hard"
        )
      )
    )


  bars_df <- bars_df %>%
    dplyr::mutate(
      order = forcats::fct_reorder(
        order,
        n_hits,
        .fun = sum,
        .desc = TRUE
      )
    )


  


  y_text <- if (italic_y) {
    element_text(
      face = "italic",
      size = 12
    )
  } else {
    element_text(
      size = 12
    )
  }


  p_bubble <- ggplot(
    bubble_df,
    aes(
      x = Rheology,
      y = .data[[feature_col]]
    )
  ) +
    geom_point(
      aes(
        size = pct,
        fill = pct
      ),
      shape = 21,
      colour = "grey25",
      alpha = 0.95
    ) +
    scale_fill_viridis_c(
      option = "E",
      end = 0.95,
      limits = c(
        0,
        max(
          bubble_df$pct,
          na.rm = TRUE
        )
      ),
      name = "Percentage (%)"
    ) +
    scale_size_area(
      max_size = 10,
      name = "Percentage (%)"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = ""
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.y = y_text,
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 12
      ),
      legend.box = "vertical",
      legend.spacing.y = grid::unit(
        2,
        "mm"
      )
    )


  


  p_bars <- ggplot() +
    geom_col(
      data = bars_df,
      aes(
        y = .data[[feature_col]],
        x = n_hits,
        fill = order
      ),
      width = 0.75,
      colour = "black",
      size = 0.25
    ) +
    scale_x_continuous(
      labels = scales::label_number(
        big.mark = ","
      ),
      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )
    ) +
    scale_fill_brewer(
      palette = "Set3",
      name = "Order"
    ) +
    labs(
      x = NULL,
      y = NULL,
      title = ""
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.text = element_text(
        face = "italic"
      )
    ) +
    guides(
      fill = guide_legend(
        ncol = 1,
        bycol = TRUE
      )
    )


  


  p_combo <- p_bubble +
    p_bars +
    patchwork::plot_layout(
      widths = c(
        25,
        75
      ),
      guides = "collect"
    )


  return(
    list(
      bubble = p_bubble,
      bars = p_bars,
      combo = p_combo,
      bubble_data = bubble_df,
      bars_data = bars_df
    )
  )
}



# 3. Fig. 7A
Fig7A <- make_fig7_panel(
  bubble_file = "data/Fig7A_bubble.csv",
  bars_file = "data/Fig7A_bars.csv",
  feature_col = "abx_class",
  italic_y = FALSE
)

Fig7A$combo

ggsave(
  filename = "output/Fig7/Fig7A_AMR.svg",
  plot = Fig7A$combo,
  device = "svg",
  width = 50,
  height = 45,
  units = "cm",
  dpi = 600
)


# 4. Fig. 7B
Fig7B <- make_fig7_panel(
  bubble_file = "data/Fig7B_bubble.csv",
  bars_file = "data/Fig7B_bars.csv",
  feature_col = "gene",
  italic_y = TRUE
)

Fig7B$combo

ggsave(
  filename = "output/Fig7/Fig7B_vitamins.svg",
  plot = Fig7B$combo,
  device = "svg",
  width = 50,
  height = 25,
  units = "cm",
  dpi = 600
)


# 5. Fig. 7C

Fig7C <- make_fig7_panel(
  bubble_file = "data/Fig7C_bubble.csv",
  bars_file = "data/Fig7C_bars.csv",
  feature_col = "gene",
  italic_y = TRUE
)

Fig7C$combo

ggsave(
  filename = "output/Fig7/Fig7C_neurotransmitters.svg",
  plot = Fig7C$combo,
  device = "svg",
  width = 50,
  height = 25,
  units = "cm",
  dpi = 600
)
