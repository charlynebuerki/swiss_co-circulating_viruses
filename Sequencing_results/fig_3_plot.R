#figure plotting for fig 3

#expects data containing strain_1, strain_2 and n for plotting
plot_matrix <- function(data)
{
  virus_abbreviations <-get_virus_abbreviation()
 
  fig <- ggplot() + 
    geom_tile(data=data, 
              aes(x=strain_1, y=strain_2, fill=n), color="grey") +
    geom_tile(data = data %>% filter(strain_1 == strain_2),
              aes(x = strain_1, y = strain_2),
              fill = "darkgrey", # Set the fill color for diagonal elements to grey
              linewidth = 0.5) +
    theme_minimal() +
    #scale_fill_gradient(low= "white",high="#403e9a") +
    scale_fill_gradient(low = "white", high="#403e9a") +
    scale_x_discrete(labels = virus_abbreviations) +
    scale_y_discrete(labels = virus_abbreviations ) +
    theme(
      axis.text.x = element_text(angle = 90, size=8, vjust = 1, hjust=1, margin = margin(t=0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size=8),
      axis.title.y = element_text(size=10),
      legend.title = element_text(size=8),
      legend.title.position = "top",
      #legend.background = element_rect(fill=NA),
      legend.text = element_text(size=8),
      #legend.position.inside = c(0.9,0.05),
      legend.position = "",
      legend.direction = "horizontal",
      plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm")
    ) +
    ylab("Co-infecting viruses") +
    xlab("")+
    labs(fill="no. co-infections")
  return(fig)
}

# plot frequency of co-infections
plot_percent_coinfection <- function(data) {
  
  virus_abbreviations <-get_virus_abbreviation()
  
  freq_plot <- ggplot() +
    #geom_bar(data=freq_data, aes(x=substrain_name, y=n), stat= "identity", fill="darkgrey") +
    geom_bar(data=data, aes(x=substrain_name, y=percent_co), stat= "identity", fill="#403e9a") +
    scale_x_discrete(labels = virus_abbreviations) +
    theme_bw() +
    theme( axis.text.x = element_blank(), #element_text(size = 8, angle = 90, vjust = 1, hjust=1),
           axis.title.x = element_blank(),
           legend.position = "",
           panel.grid.major = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_text(size=8),
           plot.margin = unit(c(0.5,0.5,0,0.5), "cm")
    ) +
    
    ylab("% co-infections") +
    scale_colour_manual(values= colors, name="virus", aesthetics = c("colour", "fill")) 
  
  return(freq_plot)
}

#aligning plots to make figure 2
make_figure_3<-function(matrix_data, percent_data, save=TRUE)
{
  f2b<-plot_matrix(matrix_data)
  f2a<- plot_percent_coinfection(percent_data)
  
  f2 <- plot_grid(f2a, f2b, nrow=2, rel_heights = c(1,2), align="v", axis="l")
  
  if(save)
  {
    ggsave("Figures/f3_coinf.pdf",f2 , dpi="retina", units="mm", width=174, height=123)
  }
  return(f2)
}