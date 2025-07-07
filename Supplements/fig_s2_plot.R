#script to generate s2 plots (pathogen correlations)


#plots the correlation figures given a pathogen and the data counts in parameter, will return a separate legend and figure
plot_correlations_pathogen<-function(data, date_limits=c("2023-03-15","2024-07-15"), pathogen)
{
  
  cor_matrix <-cor(data %>% select(n_seq, n_pcr, n_detected ,percentage), use="pairwise.complete.obs")^2
  
  cor_labels <- data.frame(
    x= c(as.Date(date_limits[2])-15,as.Date(date_limits[2])-15,as.Date(date_limits[2])-15, as.Date(date_limits[2])-15),
    y=c(max(data$percentage, na.rm=TRUE)-12,max(data$percentage, na.rm=TRUE)-8,max(data$percentage, na.rm=TRUE)-4,max(data$percentage, na.rm=TRUE)),
    type=c("n_pcr/n_detected", "n_seq/n_detected","n_detected/percentage","n_pcr/percentage"),
    label= c(
      paste0("R2 = ", round(cor_matrix["n_pcr", "n_detected"], 2)),
      paste0("R2 = ", round(cor_matrix["n_seq", "n_detected"], 2)),
      paste0("R2 = ", round(cor_matrix["n_detected", "percentage"], 2)),
      paste0("R2 = ", round(cor_matrix["n_pcr", "percentage"], 2))
    ))
  
  write.csv(cor_labels, file=paste0("Data/data/correlations/corr_sentinella_",pathogen, ".csv"), row.names=FALSE)
  
  fig<-ggplot(
    data=data %>% complete(date = seq(min(date), max(date), by = "week")), 
    aes(x=date), size=0.5) +
    #sequencing
    geom_line(aes(y=n_seq, color= "seq. hq"), size=0.2)+ 
    geom_point(aes(y=n_seq, color= "seq. hq"),  shape=".", size=1) +
    #detected
    geom_line(aes(y=n_detected, color= "seq. detected"), size=0.2)+ 
    geom_point(aes(y=n_detected, color= "seq. detected"),  shape=".", size=1) +
    #pcr
    geom_line( aes(y=n_pcr, color= "PCR"), size=0.2)+ 
    geom_point(aes(y=n_pcr, color="PCR"), shape=".", size=1) +
    #sentinella
    geom_line( aes(y=percentage, color="sentinella")) +
    #geom_text(data=cor_labels, aes(x=x, y=y, label=label), size=2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b-%y", limits=c(as.Date(date_limits[1]), as.Date(date_limits[2]))) +
    scale_y_continuous(name= "no. samples", 
                       sec.axis = sec_axis(~. , name="positivity rate")
    ) +
    # Custom color legend
    scale_color_manual(
      name = "data type",  # Legend title
      values = c("seq. hq" = "blue", "PCR" = "green", "seq. detected"="orange" ,"sentinella" = "black")
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=5),
      axis.ticks.y = element_line(color="grey"),
      axis.ticks.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(linetype = 2, linewidth = 0.1),
      plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
      axis.text.x = element_text(angle = 45, size=4, vjust = 1, hjust=1, margin = margin(t=0)),
      #legend.position = "bottom",
      #legend.direction = "horizontal",
      plot.title = element_text(size = 6),
      axis.title.x = element_text(size=5),
      axis.title.y = element_text(size=5)
    ) +  #removed guides override
    annotate("rect", xmin=as.Date("2023-06-01"), xmax=as.Date("2023-10-15"), ymin=0, ymax=max(data$percentage, na.rm=TRUE)+3, alpha=0.2, fill="grey") +
    annotate("text", x=as.Date("2023-07-24"), y=max(data$percentage, na.rm=TRUE), label= "no collection period", size=4/.pt ) +
    labs(title=str_replace_all(pathogen, "_", " ")) + theme(legend.direction = "horizontal") +
    guides(color = guide_legend(nrow = 2))
  
  legend <- get_legend(fig)
  fig <- fig + theme( legend.position = "")
  
  return(list("plot"=fig, "legend"=plot_grid(legend)))
}

#WASTEWATER plots the correlation figures given a pathogen and the data counts in parameter, will return a separate legend and figure
plot_correlations_pathogen_wastewater<-function(data, date_limits=c("2023-03-15","2024-07-15"), pathogen)
{
  
  cor_matrix <-cor(data %>% select(n_seq, n_pcr, n_detected ,seven_day_median_viral_load_avg), use="pairwise.complete.obs")^2

  cor_labels <- data.frame(
    x= c(as.Date(date_limits[2])-15,as.Date(date_limits[2])-15,as.Date(date_limits[2])-15, as.Date(date_limits[2])-15),
    y=c(max(data$n_pcr, na.rm=TRUE)-12,max(data$n_pcr, na.rm=TRUE)-8,max(data$n_pcr, na.rm=TRUE)-4,max(data$n_pcr, na.rm=TRUE)),
    type=c("n_pcr/n_detected", "n_seq/n_detected","n_detected/viralload","n_pcr/viralload"),
    label= c(
      paste0("R2 = ", round(cor_matrix["n_pcr", "n_detected"], 2)),
      paste0("R2 = ", round(cor_matrix["n_seq", "n_detected"], 2)),
      paste0("R2 = ", round(cor_matrix["n_detected", "seven_day_median_viral_load_avg"], 2)),
      paste0("R2 = ", round(cor_matrix["n_pcr", "seven_day_median_viral_load_avg"], 2))
    ))
  
  scale_factor <- min(data$seven_day_median_viral_load_avg, na.rm=TRUE) + 
    mean(data$seven_day_median_viral_load_avg, na.rm=TRUE)/5
  

  write.csv(cor_labels, file=paste0("Data/data/correlations/corr_wastewater_",pathogen, ".csv"), row.names=FALSE)
  

  fig<-ggplot(
    data=data %>% complete(date = seq(min(date), max(date), by = "week")), 
    aes(x=date), size=0.5) +
    #sequencing
    geom_line(aes(y=n_seq, color= "seq. hq"), size=0.2)+ 
    geom_point(aes(y=n_seq, color= "seq. hq"),  shape=".", size=1) +
    #detected
    geom_line(aes(y=n_detected, color= "seq. detected"), size=0.2)+ 
    geom_point(aes(y=n_detected, color= "seq. detected"),  shape=".", size=1) +
    #pcr
    geom_line( aes(y=n_pcr, color= "PCR"), size=0.2)+ 
    geom_point(aes(y=n_pcr, color="PCR"), shape=".", size=1) +
    #wastewater
    geom_line( aes(y=seven_day_median_viral_load_avg/scale_factor, color="viral load")) +
    #geom_text(data=cor_labels, aes(x=x, y=y, label=label), size=2) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b-%y", limits=c(as.Date(date_limits[1]), as.Date(date_limits[2]))) +
    scale_y_continuous(name= "no. samples", 
                       sec.axis = sec_axis(~.*scale_factor , name="7-day med. Viral Load (gc/person/day)"),
                                           limits=c(0, max(data$n_pcr, na.rm = TRUE)+10)
    ) +
    # Custom color legend
    scale_color_manual(
      name = "data type",  # Legend title
      values = c("seq. hq" = "blue", "PCR" = "green", "seq. detected"="orange" ,"viral load" = "black")
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=5),
      axis.ticks.y = element_line(color="grey"),
      axis.ticks.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(linetype = 2, linewidth = 0.1),
      plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
      axis.text.x = element_text(angle = 45, size=4, vjust = 1, hjust=1, margin = margin(t=0)),
      #legend.position = "bottom",
      #legend.direction = "horizontal",
      plot.title = element_text(size = 6),
      axis.title.x = element_text(size=5),
      axis.title.y = element_text(size=5)
    ) +  #removed guides override
    annotate("rect", xmin=as.Date("2023-06-01"), xmax=as.Date("2023-10-15"), ymin=0, ymax=max(data$n_pcr, na.rm=TRUE)+10, alpha=0.2, fill="grey") +
    annotate("text", x=as.Date("2023-07-24"), y=max(data$n_pcr, na.rm=TRUE)+5, label= "no collection period", size=4/.pt ) +
    labs(title=str_replace_all(pathogen, "_", " ")) + 
    theme(legend.direction = "horizontal") +
    guides(color = guide_legend(nrow = 2)) 
  
  legend <- get_legend(fig)
  fig <- fig + theme( legend.position = "")
  
  return(list("plot"=fig, "legend"=plot_grid(legend)))
}



#assembles and formats the data for all pathogens, generates individual pathogen plots and plots together in grid. based on pathogens give
make_supplementary_figure_2<-function(sentinella_data, pcr_data, sequencing_data, detected_data, wastewater_data, pathogens, save, date_limits=c("2023-03-15","2024-07-15"))
{
  #first assemble the data 
  
  corr_plots <- list()
  legend <- NULL
  
  corr_plots_ww <- list()
  legend_ww <- NULL
  
  #dynamically create all plots
  pathogens_wastewater <- c("influenza_A","influenza_B","respiratory_syncytial_virus","sars-cov-2")
  #plt<-plot_correlations_pathogen_wastewater(dat, pathogen = p, date_limits = date_limits)
  
  for(p in pathogens)
  {
    #format & plot data for individual pathogen
    print(p)
    dat<-format_data_correlation_pathogen(sentinella_data, pcr_data, sequencing_data, detected_data, wastewater_data, pathogen=p) 
    
    plt<-plot_correlations_pathogen(dat, pathogen=p, date_limits = date_limits)

    corr_plots[[length(corr_plots) + 1]] <- plt$plot
    
    legend <- plt$legend
    
    #wastewater
    if(p %in% pathogens_wastewater){
      plt_2<-plot_correlations_pathogen_wastewater(dat, pathogen = p, date_limits = date_limits)
      corr_plots_ww[[length(corr_plots_ww) + 1]] <- plt_2$plot
      
      legend_ww <- plt_2$legend
    }
    
    
  }
  
  corr_plots[[length(corr_plots) +1]] <- legend
  corr_plots_ww[[length(corr_plots_ww) +1]] <- legend_ww
  #put all the plots together
  # s2<-plot_grid(
  #   plot_grid(plotlist = corr_plots, nrow = 3, ncol = 2, align = "v", axis="l"),
  #   legend,  # Ensure the legend is properly stored in plt$legend
  #   ncol = 1, 
  #   rel_heights = c(1, 0.1)  # Adjust legend height
  # )
  
  s2<-plot_grid(plotlist = corr_plots, nrow = 4, ncol = 2, align = "v", axis="l")
  s2_ww<-plot_grid(plotlist = corr_plots_ww, nrow = 3, ncol = 2, rel_heights = c(1,1,0.2), align = "v", axis="l")
  
  if(save) { 
    ggsave("Figures/s1_corr_sentinella.pdf", s2,  dpi="retina", units="mm", width=174, height=123)
    ggsave("Figures/s2_corr_ww.pdf", s2_ww,  dpi="retina", units="mm", width=174, height=123)
    
    }
  
  return(list(corr=s2, corr_ww=s2_ww))
  
}