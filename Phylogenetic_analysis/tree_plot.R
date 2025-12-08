#plotting functions for tree visualizations 

#function to convert the date decimal format to YYYY-MM-dd, specifically to use from the NEXUS formatted trees outputed from Nextstrain 
decimal_to_date <- function(decimal_dates) {
  decimal_years <- as.numeric(decimal_dates)
  years <- floor(decimal_dates)  # Extract integer years
  fractions <- decimal_dates - years  # Extract fractional part
  days_in_year <- ifelse(leap_year(years), 366, 365)  # Vectorized leap year check
  date_vector <- as.Date(paste0(years, "-01-01")) + round(fractions * days_in_year)  # Compute dates
  return(date_vector)
}

#rescale the tree in date format with the given amount of major gridline breaks
make_date_tree_scale<- function(tree, breaks=8) {
  
  #dates <- tree$data$date
  x_values <- tree$data$x
  x_values <- x_values[!is.na(x_values)]
  
  root_depth <- max(x_values, na.rm=TRUE)
  
  # Define breaks in the x-axis (divergence scale)
  breaks_mapping <- seq(from=min(x_values, na.rm=TRUE), 
                        to=root_depth, 
                        length.out=breaks)
  
  #get date standard format
  dates <- decimal_to_date(as.numeric(tree$data$num_date))
  
  dates_labels <- seq(from=min(dates, na.rm=TRUE), 
                      to=max(dates, na.rm=TRUE), 
                      length.out=breaks)
  
  return(list("labels"=dates_labels, "breaks"=breaks_mapping))
}

#read locations mapping from file
get_location_colors <- function()
{
  colors <- readr::read_tsv("Data/resources/colors.tsv", col_names = c("object", "key", "color")) %>% 
    #filter(object=="location" | object=="country" ) %>% 
    select(key, color)
  
  df_colors <- setNames(colors$color,colors$key)
  
  return(df_colors)
}

#nice theme/constant theme for all tree plots with various options to customize
plot_nice_tree<-function(gg_tree, axis_breaks, legend_on, date_format="%Y", title="no title", legend_position=c(0.3,0.8), small_mode = FALSE, rotate=F)
{
  dates <- make_date_tree_scale(gg_tree, breaks=axis_breaks)
  
  # --- DEFINE SIZES BASED ON MODE ---
  if(small_mode) {
    # Tiny sizes for multipanels
    leg_key_size <- 0.2   # cm
    leg_txt_size <- 5     # pt
    leg_title_size <- 6   # pt
    plot_title_size <- 8  # pt
    margin_pad <- 1       # pt
  } else {
    # Standard sizes
    leg_key_size <- 0.3
    leg_txt_size <- 7
    leg_title_size <- 8
    plot_title_size <- 10
    margin_pad <- 3
  }
  
  gg_tree + theme_minimal() +
    theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=0.3),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title = element_text(size=plot_title_size),
         
          plot.margin = unit(c(0.1,0.0,0.1,0.0), "cm"),
          
          plot.title = element_text(size = 8),
          
          axis.text.x = if(rotate) element_text(angle = 20, hjust = 1, size=7) else element_text(size= 8) ,
          
          #legend
          legend.position = if(legend_on) "inside" else "none",
          legend.position.inside = legend_position,
          legend.direction = "horizontal",
          
          # Use variables defined above
          # Use variables defined above
          legend.key.size = unit(leg_key_size, "cm"),  
          legend.title = element_text(size=leg_title_size, face="bold"), 
          legend.text = element_text(size=leg_txt_size), 
          
          legend.spacing.x = unit(0.1, "cm"), 
          legend.spacing.y = unit(0.1, "cm"), 
          
          legend.margin = margin(t=margin_pad, r=margin_pad, b=margin_pad, l=margin_pad, unit="pt"),
          legend.box.margin = margin(0,0,0,0),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    ) +
    scale_x_continuous(
      breaks = dates$breaks,  # Replace with appropriate numeric mapping
      name = "Time",
      labels = as.character(format(dates$labels, date_format)),    # Format the date column for the labels
      expand=expansion(mult=0.1),
    ) +
    guides(size=FALSE,
           color=guide_legend(nrow=4)
    ) +
    ggtitle(title)
  
}

#plotting function for pretty tree
#tree is a tree format readfrom NEXUS format containing approrpiate metadata
#color_by: variable to color by the tree; must be included in dataset; also needs column called datasetfor plotting dots
#zoom in option: if true, must provide 2 tips to zoom in, it will find the MCRA, should also adapt the axis break
#tip_color is one of "database" or "location"
plot_tree <- function(tree, color_by, plot_title ,legend_position=c(0.3,0.8), 
                      zoomed=FALSE, tips_to_zoom = NULL, axis_breaks=8, legend_on=TRUE, 
                      date_format="%Y", save=TRUE, tip_color="database", small_mode=FALSE, rotate=FALSE)
{
  options(ignore.negative.edge = TRUE)
  
  basic_tree <- ggtree(tree, 
                       if(!is.na(color_by)) aes(color=.data[[color_by]]) else NULL, 
                           show.legend = if(tip_color=="database") TRUE else FALSE )   #color by dynamically setting coloring
  
  if (tip_color=="database"){ 
    # Generate default ggplot colors dynamically based on `color_by`
    if(!is.na(color_by))
    {
      unique_categories <- na.omit(unique(basic_tree$data[[color_by]]))
      
      # 1. Remove "other" from the automatic color generator (so it doesn't get a random color)
      unique_categories <- unique_categories[unique_categories != "other"]
      
      default_colors <- if(color_by == "country" | color_by== "region" | color_by=="Clade") get_location_colors() else setNames(scales::hue_pal()(length(unique_categories)), unique_categories) 
      
      # 2. Manually assign "other" to lightgrey
      default_colors <- c(default_colors, "other" = "lightgrey")  
      
      # 3. Apply the scale
      basic_tree <- basic_tree + scale_color_manual(values = default_colors, na.value = "lightgrey")
    }

    basic_tree <- basic_tree + 
      geom_tippoint(aes(size=.data[[tip_color]]), shape=16, color="black") +  #default visualization of revseq samples
      scale_size_discrete(range = c(1, 1))
    
  } else {
    basic_tree <- basic_tree + 
      #scale_color_manual(values = scales::hue_pal()(length(unique(tree@data[[color_by]]))), guide = "none") +  # Hide tree legend
      geom_tippoint(aes(color=.data[[tip_color]]), size=0.5) + #based on location
      scale_color_manual(values = get_location_colors(), na.translate=FALSE, expand=expansion(mult=0.1))
  }
    
  
  #nice plotting format
  basic_tree <- if (zoomed) viewClade(basic_tree, MRCA(tree, tips_to_zoom[1], tips_to_zoom[2])) else basic_tree #(leave the same) 
  
  whole_tree<-plot_nice_tree(basic_tree, axis_breaks, legend_on, date_format, plot_title, legend_position, small_mode= small_mode, rotate=rotate)
  
  if(save) ggsave("Figures/tree_plots/test.pdf", whole_tree,  dpi="retina", units="mm", width=174, height=123)
  
  return(whole_tree)
}

#plotting function for a vertical colored band next to the tree
plot_vertical_data <- function(tree_metadata, taxa_order, color_by = "country", legend_on=TRUE, legend_columns=2, small_mode=FALSE) {
  
  df <- tree_metadata %>% 
    #filter(accession %in% get_taxa_name(gg_tree)) %>% 
    mutate(accession = factor(accession, levels = rev(taxa_order))) %>% 
    arrange(accession)  # Arrange the data based on tree ordering
  
  # Generate default ggplot colors dynamically based on `color_by`
  
  unique_categories <- na.omit(unique(df[[color_by]]))
  
  # 1. Remove "other" from the automatic color generator (so it doesn't get a random color)
  unique_categories <- unique_categories[unique_categories != "other"]
  
  default_colors <- if(color_by == "country" | color_by== "region" | color_by=="Clade") get_location_colors() else setNames(scales::hue_pal()(length(unique_categories)), unique_categories) 
  
  # 2. Manually assign "other" to lightgrey
  default_colors <- c(default_colors, "other" = "#ebebeb")  
  
  # --- DEFINE SIZES ---
  if(small_mode) {
    leg_txt_size <- 5
    leg_title_size <- 6
    leg_key_size <- 0.2 # cm
  } else {
    leg_txt_size <- 6
    leg_title_size <- 8
    leg_key_size <- 0.4 # cm (default is usually bigger for side bars)
  }
  
  
  # Create the ggplot
  plt <- ggplot(data = df) +
    geom_tile_rast(
      aes(fill = .data[[color_by]], x = " ", y = accession),
      dpi = 3500,
      dev = "ragg" 
    ) +
    theme_minimal() +
    xlab(NULL) +
    ylab(NULL) +
    theme( 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      legend.position = if(legend_on) "right" else "none",
      
      legend.title = element_text(size=leg_title_size),
      legend.text = element_text(size=leg_txt_size),
      legend.key.size = unit(leg_key_size, "cm"),
      legend.spacing.y = unit(0.1, "cm"), # Tighten rows
      plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm")
    ) +
    scale_fill_manual(values = default_colors, na.value = "#ebebeb") +  # Dynamically set colors
    guides(fill=guide_legend(ncol=legend_columns))
  return(plt)
}



#given a zoomed in plot creaded by ViewClade, get the ordered list of 
#visible tips to plot 
get_visible_tips_zoomed_plot<-function(gg_tree)
{
  panel_x_range <- ggplot_build(gg_tree)$layout$panel_params[[1]]$x.range
  panel_y_range <- ggplot_build(gg_tree)$layout$panel_params[[1]]$y.range
  
  visible_tips <- gg_tree$data %>%
    filter(isTip, 
           x >= panel_x_range[1], x <= panel_x_range[2],
           y >= panel_y_range[1], y <= panel_y_range[2]) %>%
    arrange(rev(y)) %>%
    pull(label)
  
  return(visible_tips)  
}


#main function to plot tree and accompagnying vertical data bars
#color_by_vertical: string vector of variable names to color by for the vertical plots; options: region, clade_membership, subclade_membership, country, location
plot_tree_with_data <- function(tree, tree_meta, title, color_by, zoomed=FALSE, save=TRUE, 
                                tips_to_zoom=NULL, axis_breaks=8, tree_legend_on=TRUE, 
                                vertical_legend_on=c("clade_membership"=TRUE,"region"=TRUE), 
                                tree_legend_position=c(0.3,0.8), date_format= "%Y", tip_color="database", 
                                color_by_vertical=c("clade_membership","region"), 
                                vertical_plot_legend_columns=2, full_page=TRUE,
                                place_tree_legend_inside = TRUE,
                                small_mode=FALSE,
                                rotate=F
                                ) {
  
  # --- STEP 1: EXTRACT LEGENDS (Side Column) ---
  
  all_legends_list <- list()
  leg_rel_heights <- c()
  
  # 1. GENERATE DUMMY TREE 
  # We always need this to calculate taxa order, regardless of where the legend goes
  dummy_tree <- plot_tree(tree, color_by, title, tree_legend_position, zoomed, tips_to_zoom, axis_breaks, 
                          legend_on=TRUE, date_format, save=FALSE, tip_color, small_mode=small_mode, rotate=rotate)
  
  # 2. Extract Tree Legend (CONDITIONAL)
  # Only extract if requested AND if we are NOT placing it inside
  if(tree_legend_on && !place_tree_legend_inside) {
    # Create a version with the legend on the right for extraction
    dummy_tree_for_extract <- dummy_tree + theme(legend.position = "right", legend.margin = margin(b = 10))
    
    all_legends_list[["tree"]] <- get_legend(dummy_tree_for_extract)
    
    n_cats <- length(unique(na.omit(dummy_tree$data[[color_by]])))
    leg_rel_heights <- c(leg_rel_heights, n_cats) 
  }
  
  # 3. Extract Vertical Plots Legends
  taxa_order_dummy <- if(!zoomed) get_taxa_name(dummy_tree) else get_visible_tips_zoomed_plot(dummy_tree)
  
  for(v in color_by_vertical) {
    show_leg <- if(v %in% names(vertical_legend_on)) vertical_legend_on[[v]] else TRUE
    
    if(isTRUE(show_leg)) {
      dummy_vert <- plot_vertical_data(tree_meta, taxa_order_dummy, color_by=v, 
                                       legend_on=TRUE, legend_columns=vertical_plot_legend_columns,
                                       small_mode=small_mode) +
        theme(legend.margin = margin(b = 10))
      
      all_legends_list[[v]] <- get_legend(dummy_vert)
      
      n_cats <- length(unique(na.omit(tree_meta[[v]])))
      n_rows <- ceiling(n_cats / vertical_plot_legend_columns)
      leg_rel_heights <- c(leg_rel_heights, n_rows)
    }
  }
  
  # 4. Combine Side Legends
  if(length(all_legends_list) > 0) {
    combined_legend <- plot_grid(
      plotlist = all_legends_list, 
      ncol = 1, 
      align = "none", 
      rel_heights = leg_rel_heights
    )
  } else {
    combined_legend <- NULL
  }
  
  
  # --- STEP 2: GENERATE REAL PLOTS ---
  
  # 1. Real Tree
  # LOGIC CHANGE: If place_tree_legend_inside is TRUE, we set legend_on=TRUE here.
  # If it is FALSE (meaning we extracted it to the side), we set legend_on=FALSE here.
  show_tree_legend_in_plot <- (tree_legend_on && place_tree_legend_inside)
  
  real_tree <- plot_tree(tree, color_by, title, tree_legend_position, zoomed, tips_to_zoom, axis_breaks, 
                         legend_on=show_tree_legend_in_plot, 
                         date_format, save=FALSE, tip_color,
                         small_mode=small_mode,
                         rotate=rotate)
  
  taxa_order <- if(!zoomed) get_taxa_name(real_tree) else get_visible_tips_zoomed_plot(real_tree)
  
  # 2. Real Vertical Plots (Always no legend, as they are on the side)
  real_verticals <- list()
  for(vertical_type in color_by_vertical) {
    real_verticals[[vertical_type]] <- plot_vertical_data(tree_meta, taxa_order, color_by=vertical_type, 
                                                          legend_on=FALSE, 
                                                          legend_columns=vertical_plot_legend_columns,
                                                          small_mode=small_mode)
  }
  
  
  # --- STEP 3: ASSEMBLE FIGURE (Using cowplot for robustness) ---
  # Replaces 'insert_right' logic. 
  
  # 1. Combine Tree + Vertical Bars into a single list
  plot_list_step3 <- c(list(real_tree), real_verticals)
  
  # 2. Define Widths: Tree gets 10, Each bar gets 1
  #    Example: c(10, 1) or c(10, 1, 1)
  #    This ensures bars are always ~1/10th of the tree width and visible.
  width_ratios <- c(10, rep(1, length(real_verticals)))
  
  full_plot_nolegend <- plot_grid(
    plotlist = plot_list_step3,
    nrow = 1,
    align = "h",  # Align the panels horizontally (Top/Bottom edges match)
    axis = "tb",  # Strictly align the Top and Bottom axes
    rel_widths = width_ratios
  )
  
  full_plot_grob <- as.ggplot(full_plot_nolegend)
  
  
  # --- STEP 4: SAVE ---
  
  if(save) {
    if(!is.null(combined_legend)) {
      final_save_plot <- plot_grid(full_plot_grob, NULL, combined_legend, nrow=1, rel_widths=c(1, 0.05, 0.15))
    } else {
      final_save_plot <- full_plot_grob
    }
    
    ggsave(paste0("Figures/tree_plots/", title, ifelse(full_page, "_supp", "_modif"), "_tree.pdf"), 
           final_save_plot,  
           dpi=3000, 
           units="mm", 
           width= 180, 
           height= ifelse(full_page, 230, 120)
    )
  }
  
  # --- STEP 5: RETURN ---
  return(list(
    plot = full_plot_grob, 
    legend = combined_legend 
  ))
}

##### reworking some figures -- notably to collapse to a simplified tree:
#function to grab the TMCR of given clades -< returns list of nodes to collapse,
#must have clade_membership present in tree data
get_nodes_to_collapse_by_clade<-function(tree, list_of_groups, clade=TRUE)
{
  tree_phylo <- as.phylo(tree)
  tree_data_full <- as_tibble(tree)
  
  target_column <- if(clade) "clade_membership" else "country"
  
  #Calculate the MRCA node for each Clade
  nodes_to_collapse_df <- tree_data_full %>%
    filter(.data[[target_column]] %in% list_of_groups) %>%
    # getMRCA expects a list of tips, not internal nodes.
    filter(node <= Ntip(tree_phylo)) %>% 
    group_by(.data[[target_column]]) %>%
    summarise(
      # getMRCA finds the common ancestor node of all tips in the group
      mrca_node = if(n() > 1) getMRCA(tree_phylo, node) else NA_integer_
    ) %>%
    filter(!is.na(mrca_node)) # Remove groups that failed (e.g. singletons)
  
  # Extract the vector of internal nodes
  nodes_to_collapse <- nodes_to_collapse_df$mrca_node
  
  return(nodes_to_collapse)
}


#function to make the collapsed tree with list of clades to collapse, addiitional nodes is a list of nodes to collapse (numeric)
tree_plot_collapsed<-function(tree, tree_metadata, color_by, plot_title, list_of_clades, additional_nodes=NULL,
                              legend_position=c(0.3,0.5), 
                              zoomed=FALSE, tips_to_zoom = NULL, axis_breaks=8, legend_on=TRUE, 
                              date_format="%Y", save=TRUE, tip_color="database", 
                              collapse_by_clade=TRUE,
                              collapse_proportion=0.003
                              )
{
  
  basic_tree <-plot_tree(tree, color_by, plot_title ,legend_position=legend_position, 
                         zoomed = zoomed, 
                         tips_to_zoom=tips_to_zoom, 
                         axis_breaks=axis_breaks, 
                         legend_on= legend_on, 
                         date_format=date_format, 
                         save=save, 
                         tip_color=tip_color
  )
  
  #color scheme
  unique_categories <- na.omit(unique(tr@data[[color_by]]))
  
  # 1. Remove "other" from the automatic color generator (so it doesn't get a random color)
  unique_categories <- unique_categories[unique_categories != "other"]
  
  default_colors <- if(color_by == "country" | color_by== "region" | color_by=="Clade") 
    get_location_colors() 
  else setNames(scales::hue_pal()(length(unique_categories)), unique_categories) 
  
  # 2. Manually assign "other" to lightgrey
  default_colors <- c(default_colors, "other" = "lightgrey")  
  
  nodes_to_collapse <- get_nodes_to_collapse_by_clade(tree,list_of_clades,clade=collapse_by_clade)
  
  if(!is.null(additional_nodes))
  {
    nodes_to_collapse <- c(nodes_to_collapse, additional_nodes)
  }
  
  p <- basic_tree
  
  # --- STEP 1: SCALE CLADES (Has to be first) ---
  # We must apply the scaling to the Y-axis coordinates before we start 
  for(n in nodes_to_collapse) {
    # Check if node exists to prevent crash
    if(n %in% p$data$node) {
      p <- scaleClade(p, node = n, scale = collapse_proportion) 
    }
  }
  
  # --- STEP 2: COLLAPSE MAIN CLADES (Color Triangles) ---
  for(n in nodes_to_collapse) {
    if(n %in% p$data$node) {
      # Color Lookup
      node_val <- p$data %>% filter(node == n) %>% pull(!!sym(color_by))
      fill_col <- default_colors[as.character(node_val)]
      if(length(fill_col) == 0 || is.na(fill_col)) fill_col <- "lightgrey"
      
      # Collapse
      p <- collapse(p, node = n, mode = "max", fill = fill_col, color = "black")
    }
  }
  
  #add colored diamonds
  final_tree <- p +
    geom_rootedge(rootedge = 0.25, color="lightgrey") +
    theme(
      panel.grid.minor.y=element_blank()
    ) +
    scale_fill_manual(values = default_colors, 
                      na.value = "lightgrey",
                      guide = guide_legend(ncol=1)
                      ) +   # Dynamically set colors
  
    guides(fill=guide_legend(ncol=1))
    
  if(save) {
    ggsave(paste0("Figures/tree_plots/", plot_title, "_collapsed", "_tree.pdf"), 
           final_tree,  
           dpi=300, 
           units="mm", 
           width= 180, 
           height= 120
    )
  }
  
  return(final_tree)
}