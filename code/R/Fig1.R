# Fig. 1B


# 0. Libraries
library(tidyverse)
library(forcats)
library(ggplot2)
library(giscoR)
library(viridis)
library(scales)


# 1. Load metadata
metadata <- read.delim(
  "data/MetaCheeseDB_metadata.tsv",
  row.names = 1,
  check.names = FALSE
)

dir.create(
  "output/Fig1/pies",
  recursive = TRUE,
  showWarnings = FALSE
)


# 2. Maps

# 2.1 Obtain country boundaries
world <- gisco_get_countries(
  resolution = "10"
)

# 2.2 Count samples per country
metadata_map <- metadata %>%
  mutate(
    Origin = recode(
      Country,
      "USA" = "United States"
    )
  )

meta_counts <- metadata_map %>%
  group_by(Origin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )


# Join sample counts to map
world_meta <- world %>%
  left_join(
    meta_counts,
    by = c("NAME_ENGL" = "Origin")
  )


# 2.3 World map
map_world <- ggplot(world_meta) +
  
  geom_sf(
    aes(fill = n)
  ) +
  
  scale_fill_viridis_c(
    na.value = "grey90",
    name = "N° of collected cheese samples"
  ) +
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom"
  )

map_world

ggsave(
  filename = "output/Fig1/Fig1_world_map.svg",
  plot = map_world,
  width = 20,
  height = 12,
  units = "cm"
)

# 2.4 Europe inset
map_europe <- ggplot(world_meta) +
  
  geom_sf(
    aes(fill = n)
  ) +
  
  scale_fill_viridis_c(
    na.value = "grey90",
    name = "N° of collected cheese samples"
  ) +
  
  coord_sf(
    xlim = c(-10, 45),
    ylim = c(30, 70)
  ) +
  
  theme_classic() +
  
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )

map_europe


ggsave(
  filename = "output/Fig1/Fig1_Europe_map.svg",
  plot = map_europe,
  width = 10,
  height = 10,
  units = "cm"
)


# 3. Metadata levels and palettes

# 3.1 Part of cheese
part_of_cheese_levels <- c(
  "Core",
  "Core_AND_rind",
  "Rind"
)

pal_part_of_cheese <- c(
  Core          = "#023743",
  Core_AND_rind = "#cd1010",
  Rind          = "#FED789"
)


# 3.2 Milk source
milk_source_levels <- c(
  "Cow",
  "Goat",
  "Mixed",
  "Sheep",
  "Water_buffalo",
  "Yak"
)

pal_milk_source <- c(
  Cow           = "#6497B1",
  Goat          = "#6A359C",
  Mixed         = "#FFB04F",
  Sheep         = "#679C35",
  Water_buffalo = "#CD1076",
  Yak           = "#cd1010"
)

# 3.3 Animal feeding
animal_feeding_levels <- c(
  "Forage",
  "Pasture"
)

pal_animal_feeding <- c(
  Forage  = "#FED789",
  Pasture = "#023743"
)

# 3.4 Milk processing
milk_processing_levels <- c(
  "Pasteurized",
  "Raw"
)

pal_milk_processing <- c(
  Raw         = "#BE3B2C",
  Pasteurized = "#1F386E"
)

# 3.5 Thermisation
thermisation_levels <- c(
  "Facoltative",
  "No",
  "Yes"
)

pal_thermisation <- c(
  Yes         = "#F7E790FF",
  No          = "#FAB57CFF",
  Facoltative = "#73652DFF"
)

# 3.6 Pasteurization
pasteurization_levels <- c(
  "Facoltative",
  "No",
  "Yes"
)

pal_pasteurization <- c(
  Yes         = "#558934",
  No          = "#F2A241",
  Facoltative = "#0E54B6"
)

# 3.7 Skimming
skimming_levels <- c(
  "Facoltative",
  "No",
  "Partially",
  "Yes"
)

pal_skimming <- c(
  Yes         = "#090109",
  No          = "#FFEC04",
  Facoltative = "#E7CFB7",
  Partially   = "#BC000E"
)

# 3.8 Backslopping
backslopping_levels <- c(
  "No",
  "Yes"
)

pal_backslopping <- c(
  Yes = "#3D578E",
  No  = "#BFAB68"
)

# 3.9 Starter culture
starter_culture_levels <- c(
  "None",
  "NWC",
  "Selected"
)

pal_starter_culture <- c(
  None     = "#388F30",
  NWC      = "#E04B28",
  Selected = "#004042"
)

# 3.10 Inoculation of moulds
inoculation_moulds_levels <- c(
  "None",
  "Yes"
)

pal_inoculation_moulds <- c(
  Yes  = "#FED789",
  None = "#A4BED5"
)

# 3.11 Coagulation method
coagulation_method_levels <- c(
  "Acid",
  "Rennet",
  "Vegetable_rennet_from_Cynara_cardunculus",
  "Vegetable_rennet_from_Cynara_cardunculus_OR_Cynara_scolimus"
)

pal_coagulation_method <- c(
  Acid =
    "#CB6BCE",
  
  Rennet =
    "#240E31",
  
  Vegetable_rennet_from_Cynara_cardunculus =
    "#74F3D3",
  
  Vegetable_rennet_from_Cynara_cardunculus_OR_Cynara_scolimus =
    "#468892"
)

# 3.12 Ripening period
ripening_period_levels <- c(
  "None",
  "1_to_3_months",
  "3_to_6_months",
  "6_to_9_months",
  "Min_9_months"
)

pal_ripening_period <- c(
  None          = "#C5DAF6",
  `1_to_3_months` = "#A1C2ED",
  `3_to_6_months` = "#6996E3",
  `6_to_9_months` = "#4060C8",
  Min_9_months  = "#1A318B"
)

# 3.13 Rheological properties / texture
rheological_properties_levels <- c(
  "Hard",
  "Semihard",
  "Soft",
  "Very_hard"
)

pal_rheological_properties <- c(
  Hard      = "#fbaba6",
  Semihard  = "#c7c864",
  Soft      = "#64d8b0",
  Very_hard = "#f1a5f8"
)

# 3.14 Curd cutting
curd_cutting_levels <- c(
  "Ladled",
  "Large",
  "Small"
)

pal_curd_cutting <- c(
  Ladled = "#C24841",
  Large  = "#8B5B45",
  Small  = "#FFFF33"
)

# 3.15 Technology of production
technology_of_production_levels <- c(
  "Uncooked",
  "Uncooked_OR_Semicooked",
  "Semicooked",
  "Pasta_filata",
  "Cooked_curd"
)

pal_technology_of_production <- c(
  Uncooked               = "#FFE099",
  Uncooked_OR_Semicooked = "#FFAD72",
  Semicooked             = "#F76D5E",
  Pasta_filata           = "#D82632",
  Cooked_curd            = "#A50021"
)

# 3.16 Temperature of curd processing
temperature_of_curd_processing_levels <- c(
  "No_heating_beyond_coagulation",
  "Within_48",
  "48_to_70",
  "70_to_90",
  "Variable"
)

pal_temperature_of_curd_processing <- c(
  No_heating_beyond_coagulation = "#FFE099",
  Within_48                     = "#FFAD72",
  `48_to_70`                    = "#F76D5E",
  `70_to_90`                    = "#D82632",
  Variable                      = "#A50021"
)

# 3.17 Presence of rind
presence_of_rind_levels <- c(
  "Facoltative",
  "No",
  "Yes"
)

pal_presence_of_rind <- c(
  Yes         = "#FED789",
  No          = "#023743",
  Facoltative = "#FFFFBF"
)

# 4. Pie charts
pie_config <- list(
  
  Part_of_cheese = list(
    levels = part_of_cheese_levels,
    palette = pal_part_of_cheese
  ),
  
  Milk_source = list(
    levels = milk_source_levels,
    palette = pal_milk_source
  ),
  
  Animal_feeding = list(
    levels = animal_feeding_levels,
    palette = pal_animal_feeding
  ),
  
  Milk_processing = list(
    levels = milk_processing_levels,
    palette = pal_milk_processing
  ),
  
  Thermisation = list(
    levels = thermisation_levels,
    palette = pal_thermisation
  ),
  
  Pasteurization = list(
    levels = pasteurization_levels,
    palette = pal_pasteurization
  ),
  
  Skimming = list(
    levels = skimming_levels,
    palette = pal_skimming
  ),
  
  Backslopping = list(
    levels = backslopping_levels,
    palette = pal_backslopping
  ),
  
  Starter_culture = list(
    levels = starter_culture_levels,
    palette = pal_starter_culture
  ),
  
  Inoculation_of_moulds = list(
    levels = inoculation_moulds_levels,
    palette = pal_inoculation_moulds
  ),
  
  Coagulation_method = list(
    levels = coagulation_method_levels,
    palette = pal_coagulation_method
  ),
  
  Ripening_period = list(
    levels = ripening_period_levels,
    palette = pal_ripening_period
  ),
  
  Rheological_properties = list(
    levels = rheological_properties_levels,
    palette = pal_rheological_properties
  ),
  
  Curd_cutting = list(
    levels = curd_cutting_levels,
    palette = pal_curd_cutting
  ),
  
  Technology_of_production = list(
    levels = technology_of_production_levels,
    palette = pal_technology_of_production
  ),
  
  Temperature_of_curd_processing_Celsius_degree = list(
    levels = temperature_of_curd_processing_levels,
    palette = pal_temperature_of_curd_processing
  ),
  
  Presence_of_rind = list(
    levels = presence_of_rind_levels,
    palette = pal_presence_of_rind
  )
)


# 5. Functions

# 5.1 Prepare percentages
prepare_pie_data <- function(data, variable, allowed_levels) {
  
  data %>%
    
    dplyr::select(
      Country,
      Value = dplyr::all_of(variable)
    ) %>%
    
    dplyr::mutate(
      Country = stringr::str_squish(
        as.character(Country)
      ),
      
      Value = as.character(Value)
    ) %>%
    
    dplyr::mutate(
      across(
        c(Country, Value),
        ~ dplyr::na_if(.x, "")
      )
    ) %>%
    
    dplyr::filter(
      !is.na(Country)
    ) %>%
    
    dplyr::filter(
      !is.na(Value)
    ) %>%

  dplyr::mutate(
    Country = forcats::fct_drop(
      factor(Country)
    ),
    
    Value = ifelse(
      Value %in% allowed_levels,
      Value,
      NA_character_
    ),
    
    Value = factor(
      Value,
      levels = allowed_levels
    )
  ) %>%
    
    dplyr::filter(
      !is.na(Value)
    ) %>%
    
    dplyr::count(
      Country,
      Value,
      name = "n"
    ) %>%
    
    dplyr::group_by(
      Country
    ) %>%
    
    dplyr::mutate(
      pct = n / sum(n)
    ) %>%
    
    dplyr::ungroup()
}

# 5.2 Create a single pie for one country
make_single_pie <- function(
    pie_data,
    country_name,
    palette,
    allowed_levels) {
  
  data_country <- pie_data %>%
    dplyr::filter(
      Country == country_name
    )
  
  ggplot(
    data_country,
    aes(
      x = "",
      y = pct,
      fill = Value
    )
  ) +
    
    geom_col(
      width = 1,
      colour = "white"
    ) +
    
    coord_polar(
      theta = "y"
    ) +
    
    scale_fill_manual(
      values = palette,
      breaks = allowed_levels,
      limits = allowed_levels,
      na.translate = TRUE,
      na.value = "grey80",
      drop = FALSE
    ) +
    
    theme_void() +
    
    theme(
      legend.position = "none",
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )
}

# 6. Individual pie charts
countries <- metadata %>%
  dplyr::filter(
    !is.na(Country)
  ) %>%
  dplyr::count(
    Country,
    name = "n"
  ) %>%
  dplyr::arrange(
    desc(n)
  ) %>%
  dplyr::pull(
    Country
  )


for (variable in names(pie_config)) {
  
  message(
    "Generating pies for: ",
    variable
  )
  
  cfg <- pie_config[[variable]]
  
  
  pie_data <- prepare_pie_data(
    data = metadata,
    variable = variable,
    allowed_levels = cfg$levels
  )
  
  
  
  variable_dir <- file.path(
    "output/Fig1/pies",
    variable
  )
  
  dir.create(
    variable_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  
  for (country_name in countries) {
    
    p <- make_single_pie(
      pie_data = pie_data,
      country_name = country_name,
      palette = cfg$palette,
      allowed_levels = cfg$levels
    )
    
    filename <- file.path(
      variable_dir,
      paste0(
        country_name,
        ".svg"
      )
    )
    
    ggsave(
      filename = filename,
      plot = p,
      width = 2,
      height = 2,
      units = "cm"
    )
  }
}

dir.create(
  "output/Fig1/pies_faceted",
  recursive = TRUE,
  showWarnings = FALSE
)


for (variable in names(pie_config)) {
  
  cfg <- pie_config[[variable]]
  
  pie_data <- prepare_pie_data(
    data = metadata,
    variable = variable,
    allowed_levels = cfg$levels
  )
  
  
  p_faceted <- ggplot(
    pie_data,
    aes(
      x = "",
      y = pct,
      fill = Value
    )
  ) +
    
    geom_col(
      width = 1,
      colour = "white"
    ) +
    
    coord_polar(
      theta = "y"
    ) +
    
    facet_wrap(
      ~ Country
    ) +
    
    scale_fill_manual(
      values = cfg$palette,
      breaks = cfg$levels,
      limits = cfg$levels,
      na.translate = TRUE,
      na.value = "grey80",
      drop = FALSE
    ) +
    
    labs(
      fill = variable,
      title = variable
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
  
  
  ggsave(
    filename = file.path(
      "output/Fig1/pies_faceted",
      paste0(
        variable,
        ".svg"
      )
    ),
    plot = p_faceted,
    width = 20,
    height = 20,
    units = "cm"
  )
}
