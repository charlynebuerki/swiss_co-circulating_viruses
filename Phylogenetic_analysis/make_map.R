#make simple map for swiss regions
library(sf)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(data.table)
library(patchwork)

make_map <- function(small_mode = FALSE, legend_right=TRUE) 
{
  
  load_regions <- function(){ 
    regions <- st_read("Data/data/map/Boundaries_G1_REGCH_20250406_fr.shp/Boundaries_G1_REGCH_20250406_fr.shp")
    regions <- st_transform(regions, crs = 25830) 
    
    regions <- regions %>%  
      mutate(regions= case_when(
        REGCHf == "Ticino" ~ "Ticino",
        REGCHf ==  "Région lémanique" ~ "Lake Geneva",
        REGCHf == "Espace Mittelland" ~ "Espace Mittelland",
        REGCHf ==  "Nordwestschweiz" ~ "Northwestern CH.",
        REGCHf == "Zürich" ~ "Zurich",
        REGCHf == "Ostschweiz" ~ "Eastern CH.",
        REGCHf == "Zentralschweiz" ~ "Central CH.",
        TRUE ~ "non-CH."
      ))
    
    return(regions)
  }
  
  regions <- load_regions()
  
  # Define sizes based on mode
  if(small_mode) {
    # Tiny sizes for embedded plots
    txt_size <- 5
    title_size <- 6
    key_size <- 0.2 # in cm
    leg_pos <- if(legend_right)"right" else "bottom" # "right" is often better for short/wide maps than "bottom"
  } else {
    # Normal sizes for full plots
    txt_size <- 8
    title_size <- 10
    key_size <- 0.5
    leg_pos <- "right"
  }
  
  simple_map <- ggplot() +
    geom_sf(data = regions, aes(fill = regions), size=0.1) + # Added size=0.1 for thinner borders
    scale_fill_manual(values = get_location_colors()) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA), 
      plot.background = element_rect(fill = "transparent", color = NA), 
      
      # --- DYNAMIC LEGEND SIZING ---
      legend.position = leg_pos,
      
      # Reduce text size
      legend.text = element_text(size = txt_size),
      legend.title = element_text(size = title_size),
      
      # Shrink the colored squares
      legend.key.size = unit(key_size, "cm"),
      
      # Tighten the spacing between items
      legend.spacing.x = unit(0.1, "cm"),
      legend.spacing.y = unit(0.1, "cm"),
      
      # Reduce white space around the legend box
      legend.margin = margin(t=0, r=0, b=0, l=0),
      legend.box.margin = margin(0,0,0,0),
      plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm")
      
    ) +
    guides(
           fill=guide_legend(ncol = 2))
  
  return(simple_map)
}