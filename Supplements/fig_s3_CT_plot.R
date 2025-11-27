#plotting functions for additional CT value comparison

#function to plot the correlation between aligned reads, DP10 and CT values
plot_ct_value_supplements<-function(df, save=TRUE)
{
  
  #plot aligned
  aligned_reads_plot <-ggplot(data=df, aes(x=values_PCR, y=aligned)) + 
    geom_smooth(method = "lm", color = "black", size = 0.5, se = TRUE) +
    
    geom_point(aes(color=substrain_name)) + 
    
    stat_cor(method = "pearson", 
             label.x.npc = 0.01, 
             label.y.npc = 0.06, # Adjust position (top/bottom) as needed
             label.sep = ",",
             digits = 2) +
    
    facet_grid(~strain_to_plot) +
    theme_bw() +
    scale_y_log10(
      name = "No. aligned reads",
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    annotation_logticks(sides = "l") +
    theme(
      legend.position = "none",
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.title.x = element_blank()
    )
  
  
  #plot coverage
  coverage_plot <-ggplot(data=df %>% rename(Strain=substrain_name), aes(x=values_PCR, y=DP10)) + 
    geom_point(aes(color=Strain)) + 
    
    facet_grid(~strain_to_plot) +
    
    theme_bw() +
    labs(
      x= "CT value",
      y= "DP10"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y=element_blank(),
    ) +
    ylim(c(0,1.0))
  
  
  plt<-plot_grid(aligned_reads_plot, coverage_plot, nrow = 2, align="v", axis="lr", labels="AUTO")
  
  if(save)
  {
    ggsave("Figures/s3_ct_corr.pdf", plt,  dpi="retina", units="mm", width=174, height=123)
  }
  
  return(plt)
  
}