#figure plotting for fig 3

#expects data containing strain_1, strain_2 and n for plotting
plot_matrix <- function(data, strain_order = NULL) { 
  virus_abbreviations <- get_virus_abbreviation()
  
  # 1. Define the grid dimensions
  #    If strain_order is provided, use it. Otherwise calculate it.
  if(!is.null(strain_order)) {
    all_strains <- strain_order
  } else {
    all_strains <- unique(c(as.character(data$strain_1), as.character(data$strain_2)))
  }
  
  # 2. Reconstitute the Grid with explicit levels
  data_full <- data %>% 
    mutate(strain_1 = factor(strain_1, levels = all_strains),
           strain_2 = factor(strain_2, levels = all_strains)) %>%
    complete(strain_1, strain_2, fill = list(n = 0)) 
  
  fig <- ggplot() + 
    geom_tile(data = data_full, 
              aes(x = strain_1, y = strain_2, fill = n), 
              color = "grey") +
    geom_tile(data = data_full %>% filter(strain_1 == strain_2),
              aes(x = strain_1, y = strain_2),
              fill = "darkgrey", 
              linewidth = 0.5) +
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "#403e9a", na.value = "white") +
    
    # Force X-axis to follow the specific order
    scale_x_discrete(labels = virus_abbreviations, limits = all_strains) + 
    scale_y_discrete(labels = virus_abbreviations, limits = all_strains) +
    
    theme(
      axis.text.x = element_text(angle = 90, size=8, vjust = 0.5, hjust=1, margin = margin(t=0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size=8),
      axis.title.y = element_text(size=10),
      legend.title = element_text(size=8),
      legend.title.position = "top",
      legend.text = element_text(size=8),
      legend.position = "", 
      legend.direction = "horizontal",
      plot.margin = unit(c(0.5, 0.5, 0.0, 0.0), "cm")
    ) +
    ylab("Co-infecting viruses") +
    xlab("") +
    labs(fill = "no. co-infections")
  
  return(fig)
}
# plot frequency of co-infections
plot_percent_coinfection <- function(data, strain_order = NULL) { 
  
  virus_abbreviations <- get_virus_abbreviation()
  
  freq_plot <- ggplot() +
    geom_bar(data=data, aes(x=substrain_name, y=percent_co), stat= "identity", fill="#403e9a") +
    
    # CRITICAL FIX: Use 'limits' to force the exact same x-axis items and order as the heatmap
    scale_x_discrete(labels = virus_abbreviations, limits = strain_order) +
    
    theme_bw() +
    theme( 
      axis.text.x = element_blank(), 
      axis.title.x = element_blank(),
      legend.position = "",
      panel.grid.major = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size=8),
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")
    ) +
    ylab("% co-infections")
  
  return(freq_plot)
}

#aligning plots to make figure 2
make_figure_3 <- function(matrix_data, percent_data, save=TRUE)
{
  # 1. Calculate the MASTER list of strains from the matrix data
  #    (The union of all rows and columns).
  all_strains <- sort(unique(c(as.character(matrix_data$strain_1), as.character(matrix_data$strain_2))))
  
  # 2. Pass this order to fcts
  f2b <- plot_matrix(matrix_data, strain_order = all_strains)
  f2a <- plot_percent_coinfection(percent_data, strain_order = all_strains)
  
  # 3. Align vertically ("v") and by axis ("lr")
  f2 <- plot_grid(f2a, f2b, nrow=2, rel_heights = c(1,2), align="v", axis="lr")
  
  if(save)
  {
    
    ggsave("Figures/f3_coinf.pdf", 
           f2,  
           dpi=300, 
           units="mm", 
           width= 180, 
           height= 120)
  }
  return(f2)
}