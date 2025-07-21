#functions to plot the coverage plots 

#requires column to have Name, start, end from gff3 files
plot_genome_annotation<-function(data)
{
  x_min=min(data$start, na.rm=TRUE)
  x_max=max(data$end, na.rm=TRUE)
  
  annotation<-ggplot(data=data) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = Name), color = "black", size=0.2) +
    geom_text(aes(x = (start + end) / 2, y = 0.5, label = Name), size = 2, color = "white") +
    scale_y_continuous(expand = c(0, 0)) +  # Remove extra space on y-axis
    #scale_x_continuous(breaks = seq(x_min,x_max, by=1200))+
    theme_minimal() +
    theme(
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.title.y = element_blank(),
      panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      legend.position = "",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm")
    ) +
    labs(x = NULL, y = NULL, fill = "Label")
  
  return(annotation)
}

#calculates the mean coverage over the genomes, requires pseudo_id, idx, and depth columns
format_mean_coverage<-function(data)
{
  mean_coverage <- data %>% select(pseudonymized_id, idx, depth) %>% 
    group_by(idx) %>% 
    summarize(mean_depth= mean(depth), .groups="drop") %>% 
    mutate(idx=as.numeric(idx))
  
  return(mean_coverage)
}

#plots the coverage (mean) over the whole genome; requires mean depth and idx 
plot_mean_coverage_genome<-function(data)
{
  mean_cov_fig<- ggplot(data=data %>% 
                          mutate(mean_depth= ifelse(mean_depth==0, NA_real_, mean_depth)))+
    geom_line(aes(x=idx, y=mean_depth), color="black", size=0.2) +
    scale_y_log10(limits=c(1, NA))+
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(0.2,0.0,0.0,0.0), "cm"),
      axis.title.y = element_text(size=4),
      axis.text.y = element_text(size=4),
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size=0.1),
    ) +
    ylab("avg. coverage \n(log)")
  
  return(mean_cov_fig)
}

#calculates the median coverage over the genomes, requires pseudo_id, idx, and depth columns returns in log
format_median_coverage<-function(data)
{
  med_coverage <- data %>% select(pseudonymized_id, idx, depth) %>% 
    group_by(idx) %>% 
    summarize(
      med_depth = median(log10(depth)), 
      p10 = quantile(log10(depth), 0.1),
      p90 = quantile(log10(depth), 0.9),
      .groups="drop") %>% 
    mutate(idx=as.numeric(idx))
  
  return(med_coverage)
}

#plots the coverage (log median) over the whole genome; requires med depth and idx and quantiles 10 and 90; floors negative values to 1
plot_median_coverage_genome<-function(data)
{
  med_cov_fig<- ggplot(data=data %>% 
                         mutate(med_depth= ifelse(med_depth < 1, 1, med_depth),
                                p10= ifelse(p10 <1, 1, p10),
                                p90= ifelse(p90 <1, 1, p90),           
                         )
  )+
    geom_ribbon(aes(x=idx,ymin = p10, ymax = p90), fill="black", alpha = 0.2) +
    geom_line(aes(x=idx, y=med_depth), color="black", size=0.2) +
    scale_y_log10(limits=c(1, NA))+
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(0.2,0.0,0.0,0.0), "cm"),
      axis.title.y = element_text(size=4),
      axis.text.y = element_text(size=4),
      panel.grid.minor.x = element_blank(), 
      panel.grid.minor.y=element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size=0.1),
    ) +
    ylab("med. coverage \n(log)")
  
  return(med_cov_fig)
}


#plots heatmap of genome coverage; requires pseudo id, idx, and depth
plot_coverage_genome_heatmap<-function(data)
{
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  min_x <-min(data$idx, na.rm = TRUE)
  max_x <-max(data$idx, na.rm = TRUE)
  
  #ensure that it is plotted chronologically:
  data <- data %>%
    arrange(ent_date) %>% 
    mutate(pseudonymized_id= factor(pseudonymized_id, levels=rev(unique(pseudonymized_id))))
  
  fig <-ggplot(data=data) +
    geom_tile(aes(x=idx, y=pseudonymized_id, fill=depth))+ 
    #scale_fill_gradient(low="#EBA091", high="#403e9a", trans="log", labels=scaleFUN) +
    scale_x_continuous(breaks = seq(min_x,max_x, by=1200))+
    scale_fill_viridis_c(option = "magma", trans="log10", direction=-1) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size=5), 
      axis.title.x = element_text(size=5),
      axis.text.y = element_text(size=4),
      legend.title = element_text(size=6),
      legend.title.position = "top",
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.text = element_text(size=4),
      legend.key.size = unit(0.3, "lines"),
      axis.ticks.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.02,0.5),
      legend.direction="vertical",
      plot.margin = unit(c(0.0,0.0,0.0,0.0), "cm"),
      axis.title.y = element_blank(),
      legend.spacing.x = unit(0.1, "cm"),  # Reduce spacing between legend items
      legend.spacing.y = unit(0.1, "cm")   # Reduce vertical space
    ) +
    xlab("Position")
  
  
  return(fig)
}

#main function to launch figure-generating process depending on pathogen segmentation
make_figure_coverage_genome<-function(data, pathogen_name,save)
{
  if(pathogen_name=="2015(H1N1)" | pathogen_name == "2012(H3N2)" | pathogen_name == "B/Brisbane")
  {
    make_figure_coverage_genome_influenza(data, pathogen_name, save)
  }
  else
  {
    make_figure_coverage_genome_single(data, pathogen_name, save)
  }
}


#main function to assemble the coverage sub figures of single segment viruses
make_figure_coverage_genome_single<-function(data, pathogen_name, save)
{
  virus_reference <-data$virus_ref %>% unique() #extract virus reference
  
  
  #gff annotation
  annotation_genome<-format_gff3_data(virus_reference) 
  annotation_plot <-plot_genome_annotation(annotation_genome)
  
  #mean coverage
  mean_coverage<-format_median_coverage(data)
  mean_coverage_plot<-plot_median_coverage_genome(mean_coverage)
  
  # heat map
  coverage_plot <- plot_coverage_genome_heatmap(data)
  
  #put it all together
  title <- ggdraw() + draw_label(pathogen_name, size = 4, x=0)
  
  fig<-plot_grid(title, mean_coverage_plot, annotation_plot, coverage_plot, 
                 nrow=4, rel_heights = c(0.03,0.1,0.04,1), align="v", axis="l")
  
  if(save)
  {
    ggsave(paste0("Figures/coverage_plots/", str_replace_all(pathogen_name, " ", "_"), "_coverage_plot.pdf"),fig , dpi="screen", units="mm", width=174, height=123)
    
  }
  
  return(fig)
}

#handle influenza cases to map to proper reference 
map_accession <- function(x) {
  switch(x,
         "KU509707.1" = "h1n1_NS", #h1n1
         "KU509706.1" = "h1n1_MP",
         "KU509705.1" = "h1n1_NA",
         "KU509701.1" = "h1n1_PB1",
         "KU509700.1" = "h1n1_PB2",
         "KU509702.1" = "h1n1_PA",
         "KU509703.1" = "h1n1_HA",
         "KU509704.1" = "h1n1_NP", 
         "KJ942620.1" = "h3n2_NS",#h3n2
         "KJ942617.1" = "h3n2_MP",
         "KJ942618.1" = "h3n2_NA",
         "KJ942622.1" = "h3n2_PB1",
         "KJ942623.1" = "h3n2_PB2",
         "KJ942621.1" = "h3n2_PA",
         "KJ942616.1" = "h3n2_HA",
         "KJ942619.1" = "h3n2_NP",
         "KC866606.1" = "vic_NS",#B/vic
         "KC866607.1" = "vic_MP",
         "FJ766839.1" = "vic_NA",
         "KC866603.1" = "vic_PB1",
         "KC866604.1" = "vic_PB2", 
         #missing PA gene
         "KX058884.1" = "vic_HA",
         "KC866605.1" = "vic_NP",
         "Unknown")  # Default case if not matched
}

map_segment_names <- function(segment) {
  switch(segment,
         "h1n1_NS"  = "KU509707.1", #h1n1
         "h1n1_MP"  = "KU509706.1",
         "h1n1_NA"  = "KU509705.1",
         "h1n1_PB1" = "KU509701.1",
         "h1n1_PB2" = "KU509700.1",
         "h1n1_PA"  = "KU509702.1",
         "h1n1_HA"  = "KU509703.1",
         "h1n1_NP"  = "KU509704.1",
         "h3n2_NS"  = "KJ942620.1", # h3n2
         "h3n2_MP"  = "KJ942617.1",
         "h3n2_NA"  = "KJ942618.1",
         "h3n2_PB1" = "KJ942622.1",
         "h3n2_PB2" = "KJ942623.1",
         "h3n2_PA"  = "KJ942621.1",
         "h3n2_HA"  = "KJ942616.1",
         "h3n2_NP"  = "KJ942619.1",
         "vic_NS"  = "KC866606.1", # B/vic
         "vic_MP"  = "KC866607.1",
         "vic_NA"  = "FJ766839.1",
         "vic_PB1" = "KC866603.1",
         "vic_PB2" = "KC866604.1",
         # Missing PA gene
         "vic_HA"  = "KX058884.1",
         "vic_NP"  = "KC866605.1",
         "Unknown")  # Default case if not matched
}


#main function to assemble the coverage sub figures containing segmented genomes
#contains 8 segments
make_figure_coverage_genome_influenza<-function(data, pathogen_name, save)
{
  #creating list of figures 
  figs <- list()
  segment_names <-NULL
  virus_reference <-data$virus_ref %>% unique() #extract virus references
  
  if(length(virus_reference) >1) { #case of influenza -- 8 segments (or 7 for victoria)
    # get names of segments
    segment_names <- sapply(virus_reference, map_accession)
  }
  
  for(segment_name in segment_names)
  {
    #gff annotation
    annotation_genome<-format_gff3_data(paste0("influenza/",segment_name))
    annotation_plot <-plot_genome_annotation(annotation_genome)
    
    accession<- map_segment_names(segment_name)
    
    #mean coverage
    mean_coverage<-format_median_coverage(data %>% filter(virus_ref==accession) )
    mean_coverage_plot<-plot_median_coverage_genome(mean_coverage)
    
    # heat map
    coverage_plot <- plot_coverage_genome_heatmap(data %>% filter(virus_ref==accession))
    
    #put it all together
    title <- ggdraw() + draw_label(segment_name, size = 4, x=0)
    
    fig<-plot_grid(title, mean_coverage_plot, annotation_plot, coverage_plot, 
                   nrow=4, rel_heights = c(0.03,0.1,0.04,1), align="v", axis="l")
    
    if(save)
    {
      ggsave(paste0("Figures/coverage_plots/influenza/", str_replace_all(segment_name, " ", "_"), "_coverage_plot.pdf"),fig , dpi="retina", units="mm", width=174, height=123)
    }
    
    #add semgent to list of plots
    figs[[segment_name]] <- fig
  }
  
  #return named list of segment figures
  return(figs)
  
  
}

#------------------------------------------
#example run:
#get data
#virus_strain <- "Human coronavirus 229E"
#test_data <-  seq_data %>% filter(grepl(virus_strain, substrain_name, fixed=TRUE))

#df <-format_coverage_plot_data(test_data, virus_strain)

#make_figure_coverage_genome(df, virus_strain, save=TRUE)