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
plot_nice_tree<-function(gg_tree, axis_breaks, legend_on, date_format="%Y", title="no title", legend_position=c(0.3,0.8))
{
  dates <- make_date_tree_scale(gg_tree, breaks=axis_breaks)
  
  gg_tree + theme_minimal() +
    theme(panel.grid.major.x=element_line(color="grey20", linetype="dotted", size=0.3),
          panel.grid.major.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.title = element_text(size=8),
          legend.title.position = "top",
          #legend.background = element_rect(fill = alpha("white", 0.8)),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.8, "lines"),
          legend.position = if(legend_on) "inside" else "none",
          legend.position.inside = legend_position,
          legend.direction= "horizontal",
          plot.margin = unit(c(0.5,0.0,0.5,0.0), "cm"),
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
plot_tree <- function(tree, color_by, plot_title ,legend_position=c(0.3,0.8), zoomed=FALSE, tips_to_zoom = NULL, axis_breaks=8, legend_on=TRUE, date_format="%Y", save=TRUE, tip_color="database")
{
  options(ignore.negative.edge = TRUE)
  
  basic_tree <- ggtree(tree, 
                       if(!is.na(color_by)) aes(color=.data[[color_by]]) else NULL, 
                           show.legend = if(tip_color=="database") TRUE else FALSE )   #color by dynamically setting coloring
  
  if (tip_color=="database"){ 
    # Generate default ggplot colors dynamically based on `color_by`
    # Generate default ggplot colors dynamically based on `color_by`
    if(!is.na(color_by))
    {
      unique_categories <- na.omit(unique(basic_tree$data[[color_by]]))  # Extract unique non-NA values
      
      default_colors <- if(color_by == "country" | color_by== "region") get_location_colors() else setNames(scales::hue_pal()(length(unique_categories)), unique_categories) 
      
      # Add NA color
      default_colors <- c(default_colors, "NA" = "lightgrey")  # Assign light grey for NA
      
      basic_tree <- basic_tree + scale_color_manual(values = default_colors, na.value = "lightgrey")  # Dynamically set colors
    }

    basic_tree <- basic_tree + 
      geom_tippoint(aes(size=.data[[tip_color]]), shape=16, color="black") #default visualization of revseq samples
      
  } else {
    basic_tree <- basic_tree + 
      #scale_color_manual(values = scales::hue_pal()(length(unique(tree@data[[color_by]]))), guide = "none") +  # Hide tree legend
      geom_tippoint(aes(color=.data[[tip_color]]), size=1.5) + #based on location
      scale_color_manual(values = get_location_colors(), na.translate=FALSE, expand=expansion(mult=0.1))
  }
    
  

  
  #nice plotting format
  basic_tree <- if (zoomed) viewClade(basic_tree, MRCA(tree, tips_to_zoom[1], tips_to_zoom[2])) else basic_tree #(leave the same) 
  
  whole_tree<-plot_nice_tree(basic_tree, axis_breaks, legend_on, date_format, plot_title, legend_position)
  
  if(save) ggsave("Figures/tree_plots/test.pdf", whole_tree,  dpi="retina", units="mm", width=174, height=123)
  
  return(whole_tree)
}

#plotting function for a vertical colored band next to the tree
plot_vertical_data <- function(tree_metadata, taxa_order, color_by = "country", legend_on=TRUE, legend_columns=2) {
  
  df <- tree_metadata %>% 
    #filter(accession %in% get_taxa_name(gg_tree)) %>% 
    mutate(accession = factor(accession, levels = rev(taxa_order))) %>% 
    arrange(accession)  # Arrange the data based on tree ordering
  
  # Generate default ggplot colors dynamically based on `color_by`
  unique_categories <- na.omit(unique(df[[color_by]]))  # Extract unique non-NA values
  
  default_colors <- if(color_by != "location" & color_by != "country" & color_by != "region") setNames(scales::hue_pal()(length(unique_categories)), unique_categories) else get_location_colors()
  
  # Add NA color
  default_colors <- c(default_colors, "NA" = "lightgrey")  # Assign light grey for NA
  
  
  # Create the ggplot
  plt <- ggplot(data = df) +
    geom_tile(aes(fill = .data[[color_by]], x = " ", y = accession)) +  # Dynamic fill
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
      legend.title = element_text(size=8),
      legend.text = element_text(size=6),
      plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm")
    ) +
    scale_fill_manual(values = default_colors, na.value = "lightgrey") +  # Dynamically set colors
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
plot_tree_with_data<-function(tree, tree_meta, title ,color_by, zoomed=FALSE, save=TRUE, tips_to_zoom=NULL, axis_breaks=8, tree_legend_on=TRUE, 
                              vertical_legend_on=c("clade_membership"=TRUE,"region"=TRUE), tree_legend_position=c(0.3,0.8) ,date_format= "%Y", tip_color="database", 
                              color_by_vertical=c("clade_membership","region"), vertical_plot_legend_columns=2)
{
  tree_plot <-plot_tree(tree, color_by, title ,tree_legend_position,zoomed, tips_to_zoom, axis_breaks, tree_legend_on, date_format, save ,tip_color)
  
  #get order of tips based on whether it's a zoomed tree or not
  taxa_order <- if(!zoomed) get_taxa_name(tree_plot) else get_visible_tips_zoomed_plot(tree_plot)
  vertical_plots <- list()
  for(vertical_type in color_by_vertical)
  {
    vertical_plots[[vertical_type]] <- plot_vertical_data(tree_meta, taxa_order, color_by=vertical_type, legend_on=vertical_legend_on[[vertical_type]], legend_columns=vertical_plot_legend_columns)
  }
  
  
  full_plot <- vertical_plots[[1]] %>% insert_left(tree_plot, width=10) %>% insert_right(vertical_plots[[2]])
  
  if(save)
  {
    ggsave(paste0("Figures/tree_plots/", title, "_tree.pdf"), full_plot,  dpi="retina", units="mm", width=174, height=123)
  }
  
  return(full_plot)
  
}