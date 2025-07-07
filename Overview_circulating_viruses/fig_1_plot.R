#plotting  for figure 1

# returns virus abbreviations for plot
get_virus_abbreviation<-function()
{
  abbreviations <- read.csv("Data/resources/virus_abbreviations.csv", header = TRUE, stringsAsFactors = FALSE)
  abbreviation_scheme <- as_vector(abbreviations$virus_abbr)
  names(abbreviation_scheme) <- abbreviations$virus
  return(abbreviation_scheme)
}

#controls virus order appearance on figure panel
get_virus_order<-function(PCR=FALSE)
{
  if(!PCR)
  {
    virus_order <-c("Human coronavirus OC43", "Human coronavirus HKU1", "Human coronavirus 229E","Human coronavirus NL63", #seasonal coronaviruses
                    "SARS-CoV-2",
                    "Human parainfluenza virus 1" , #parainfluenzas
                    "Human parainfluenza virus 2", 
                    "Human parainfluenza virus 3",
                    "Human parainfluenza virus 4a",
                    "Influenza B virus (B/Brisbane/60/2008)", #influenza B
                    "Influenza A virus (A/Michigan/45/2015(H1N1))", #influenza A
                    "Influenza A virus (A/Texas/50/2012(H3N2))",
                    "Respiratory syncytial virus (type A)", #RSVs
                    "Human Respiratory syncytial virus 9320 (type B)",
                    "Human metapneumovirus A",   
                    "Human metapneumovirus B",  #metapneumo
                    "Human adenovirus B1" , #adenoviruses
                    "Human adenovirus C2",
                    "Human parechovirus type 1 PicoBank/HPeV1/a", #picornaviridae family
                    "Human enterovirus C109 isolate NICA08-4327" ,
                    "Human rhinovirus A89",
                    "Human rhinovirus B14",
                    "WU polyomavirus",
                    "KI polyomavirus Stockholm 60",
                    "Human bocavirus 1 (Primate bocaparvovirus 1 isolate st2)"
    ) 
  } else 
  {
    virus_order <-c("coronavirus OC43", "coronavirus HKU1", "coronavirus 229E","coronavirus NL63", #seasonal coronaviruses
                    "SARS-CoV-2",
                    "Parainfluenza 1" , #parainfluenzas
                    "Parainfluenza 2", 
                    "Parainfluenza 3",
                    "Parainfluenza 4",
                    "Influenza B", #influenza B
                    "Influenza A", #influenza A
                    "RSV - A/B", #RSVs
                    "Metapneumovirus",  #metapneumo
                    "Adenovirus", #adenoviruses
                    "Rhino- / Enterovirus"
                    
    ) 
  }
  return (virus_order)
}

## MAIN HEATMAP FIGURE 
#assumes that data contains a column line_positions, date; pathogen_type_to_plot is the the name of the column containing pathogen names
#value to fill, name of column to use when filling; legend name is the legend name, whether to log transform or not, date_limits is the scale
#of limits for the plot;legend position is a coordinate 2d for inside positionin of legend
plot_grid_pathogen <- function(data, pathogen_type_to_plot, value_to_fill, legend_name, log_transform=TRUE, date_limits=c("2023-03-15","2024-07-15"), legend_position)
{
  #ad-hoc functions
  scaleFUN <- function(x) sprintf("%.0f", x)
  virus_abbreviations <-get_virus_abbreviation()
  
  
  bk <- ggplot()+ geom_hline(data=data, aes(yintercept = line_positions), color="grey")
  
  if(log_transform) {
    bk <-bk + scale_fill_gradient(low="#EBA091", high="#403e9a", trans="log", labels=scaleFUN) 
  } else {
    bk <- bk + scale_fill_gradient(low="#EBA091", high="#403e9a", labels=scaleFUN) 
  }
  
  fig<-bk + geom_tile(data=data , aes(date, !!sym(pathogen_type_to_plot), fill=!!sym(value_to_fill))) +
    labs(fill=legend_name) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b-%y", limits=c(as.Date(date_limits[1]), as.Date(date_limits[2]))) +
    scale_y_discrete(labels = virus_abbreviations) +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_text(size=6),
          axis.ticks.y = element_line(color="grey"),
          legend.title = element_text(size=6),
          legend.title.position = "top",
          legend.background = element_rect(fill = alpha("white", 0.8)),
          legend.text = element_text(size=4),
          legend.key.size = unit(0.8, "lines"),
          axis.ticks.x = element_blank(),
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(linetype = 2),
          legend.position = "inside",
          legend.position.inside = legend_position,
          #legend.position.inside = c(0.94,0.06),
          legend.direction="horizontal",
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          axis.title.y = element_blank()
    )  + #removed guides override
    annotate("rect", xmin=as.Date("2023-06-01"), xmax=as.Date("2023-10-15"), ymin=0.5, ymax=max(data$line_positions+1, na.rm=TRUE), alpha=0.2, fill="grey")
  return(fig)
}


#PLOTS THE FREQUENCIES BARCHART 
#requires plot of type count made from function format frequencies plot
plot_frequencies<-function(freq_data)
{
  fig <- ggplot() +
    geom_bar(data = freq_data, aes(x = date, y = n_all, fill = "PCR pos."), stat = "identity") +
    geom_bar(data = freq_data, aes(x = date, y = n_detect, fill = "seq. detect"), stat = "identity", alpha = 0.8) +
    geom_bar(data = freq_data, aes(x = date, y = n_hq, fill = "seq. HQ"), stat = "identity", alpha = 0.8) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b-%y", limits = c(as.Date("2023-03-15"), as.Date("2024-07-15"))) +
    scale_fill_manual(
      name = "Sample Type",
      values = c("PCR pos." = "lightgrey", "seq. detect" = "darkgrey", "seq. HQ" = "#636363"),
      labels = c("PCR pos." = "PCR pos.", "seq. detect" = "seq. detected", "seq. HQ" = "seq. HQ")
    ) +
    theme_bw() +
    theme(axis.text.x = element_text( size=6), 
          axis.text.y = element_text(size=8), 
          axis.title.y = element_text(size=6),
          axis.title.x = element_text(size=7),
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.x=element_line(linetype = 2),
          plot.margin = unit(c(0.5,0.5,0.0,0.0), "cm"),
          #plot.margin = unit(c(0.5,0.5,0.0,0.5), "cm"),
          legend.text= element_text(size=5),
          legend.title = element_text(size=6),
          legend.title.position = "top",
          legend.background = element_rect(fill = alpha("white",1.0)),
          legend.key.size = unit(0.8, "lines"),
          legend.position = "inside",
          legend.position.inside =c(0.32,0.5),#c(.93,.89),
          legend.direction="vertical"
    ) +
    xlab("Time")+
    ylab("no. \nsamples") +
    annotate("rect", xmin=as.Date("2023-06-01"), xmax=as.Date("2023-10-15"), ymin=0, ymax=max(freq_data$n_all, na.rm=TRUE)+1, alpha=0.2, fill="grey")+
    annotate("text", x=as.Date("2023-07-24"), y=min(freq_data$n_all, na.rm=TRUE)+10, label= "no collection period", size=8/.pt ) 
  
  return(fig)
  
}

#Function that puts all the plot together and aligns them with the option of highlighting a strain family. 
make_figure_one<-function(data, frequencies_data ,sentinella_data, highlight, substrain_to_highlight="", save, pcr=FALSE)
{
  f1a <- plot_grid_pathogen(df_plt, pathogen_type_to_plot = "substrain_name", value_to_fill = "n", 
                            legend_name = "no. samples", log_transform = TRUE, 
                            legend_position=c(0.32,0.08) )
  #frequency subfigure 
  f1b<- plot_frequencies(frequencies_data)
  
  #if we're highlighting 
  if(highlight)
  {
    df_plt <- data %>% 
      rowwise() %>% 
      mutate(
        highlight_max = ifelse(substrain_name %in% substrain_to_highlight, highlight_max, NA),
        highlight_min =ifelse(substrain_name %in% substrain_to_highlight, highlight_min, NA))
    
    f1a <- f1a+geom_rect(data=df_plt %>% filter(!is.na(highlight_max)) %>% group_by(substrain_name) %>%  sample_n(1), 
                         aes(xmin=as.Date("2023-03-15"), xmax = as.Date("2024-07-15") , 
                             ymin = highlight_max, ymax = highlight_min), 
                         alpha=0.5, fill=alpha("lightgreen", 0.5))
    f1b <- f1b + geom_bar(data=frequencies_data, aes(x=date, y=n_highlight), stat="identity", fill="lightgreen", alpha=0.5)
  }
  
  #sentinella data
  f1c <- plot_grid_pathogen(sentinella_data, pathogen_type_to_plot = "pathogen_name", 
                            value_to_fill = "percentage", legend_name= "% positivity", 
                            log_transform=FALSE,legend_position=c(0.32,0.08))
  
  
  f1<-plot_grid(f1a, f1b, f1c, nrow=3, rel_heights = c(2.5,1,1), align="v", axis="l")
  
  if(save)
  {
    if(!pcr){
      substrain_to_highlight <- substrain_to_highlight[1]
      ggsave(paste0("Figures/f1_DP10", ifelse(highlight, paste0("_highlight_", getElement(virus_abbreviations, substrain_to_highlight)), ""), ".pdf"),
             f1 , dpi="retina", units="mm", width=174, height=123)
      
    } else
    {
      ggsave("Figures/f1_pcr.pdf", f1,  dpi="retina", units="mm", width=174, height=123)
    }
    
  }
  return(f1)
}