#script to generate s2 data formats for pathogen correlations figure


#formats the 3 types of data: sentinella, PCR and sequencing counts for the particular pathogen given in parameter
#outputs a format with a date column, pathogen type, and counts of each by week
format_data_correlation_pathogen<-function(sentinella, counts_pcr, counts_sequencing, counts_detected, wastewater, pathogen)
{
  
  wastewater_pathogens <-c("influenza_A", "influenza_B", "sars-cov-2", "respiratory_syncytial_virus")
  
  
  #case the pathogen is "other"
  if(pathogen == "other")
  {
    pathogens_to_exclude <- c("influenza_A", "influenza_B", "sars-cov-2", "adenovirus", "rhinovirus", "respiratory_syncytial_virus")
    
    sent_data <- sentinella %>% 
      filter(pathogen_name == pathogen) %>% 
      mutate(pathogen_name =  tolower(str_replace_all(pathogen, "_", " "))) %>% 
      select(date, pathogen_name, percentage) %>% 
      mutate(week=lubridate::floor_date(date, "week"))
    
    pathogen_data <- counts_sequencing %>% 
      filter(!grepl(paste(str_replace_all(pathogens_to_exclude, "_", " "), collapse = "|"), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(pathogen_name = tolower(str_replace_all(pathogen, "_", " "))) %>% 
      rename(n_seq = n)
    
    detected_data <- counts_detected %>% 
      filter(!grepl(paste(str_replace_all(pathogens_to_exclude, "_", " "), collapse = "|"), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(pathogen_name = tolower(str_replace_all(pathogen, "_", " "))) %>% 
      rename(n_detected = n)   
    
    
    
  } else #other defined scenarios
  {
    sent_data <- sentinella %>% 
      filter(pathogen_name == pathogen) %>% 
      mutate(pathogen_name =  tolower(str_replace_all(pathogen, "_", " "))) %>% 
      select(date, pathogen_name, percentage) %>% 
      mutate(week=lubridate::floor_date(date, "week"))
    
    pathogen_data <- counts_sequencing %>% 
      filter(grepl(str_replace_all(pathogen, "_", " "), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(pathogen_name = tolower(str_replace_all(pathogen, "_", " "))) %>% 
      rename(n_seq = n)
    
    detected_data <- counts_detected %>% 
      filter(grepl(str_replace_all(pathogen, "_", " "), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(pathogen_name = tolower(str_replace_all(pathogen, "_", " "))) %>% 
      rename(n_detected = n)
    
    if(pathogen %in% wastewater_pathogens)
    {
      wastewater_data <- wastewater %>% 
        filter(pathogen_name == pathogen) %>% 
        mutate(pathogen_name =  tolower(str_replace_all(pathogen, "_", " "))) %>% 
        select(date, pathogen_name, seven_day_median_viral_load_avg)%>% 
        mutate(week=lubridate::floor_date(date, "week"))
    }
    
  }
  
  
  
  #if we have multi-substrain, we should consider them together:
  if(length(unique(pathogen_data$substrain_name, na.rm=TRUE)) > 1) {
    pathogen_data <- pathogen_data %>% 
      group_by(date, pathogen_name) %>% 
      summarise(n_seq=sum(n_seq), .groups="drop")
    
    detected_data <- detected_data %>% 
      group_by(date, pathogen_name) %>% 
      summarise(n_detected=sum(n_detected), .groups="drop")
  }
  
  if(pathogen == "rhinovirus") {pathogen_pcr = "rhino-"} 
  else if(pathogen == "respiratory_syncytial_virus") {pathogen_pcr= "RSV"}
  else { pathogen_pcr = pathogen}
  
  
  #taking care of PCR in other case
  if(pathogen == "other")
  {
    pathogens_to_exclude <- c("influenza_A", "influenza_B", "sars-cov-2", 
                              "adenovirus", "rhino-", "RSV")
    
    pathogen_pcr <- counts_pcr %>% 
      filter(!grepl(paste(str_replace_all(pathogens_to_exclude, "_", " "), collapse = "|"), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(substrain_name = tolower((str_replace_all(pathogen, "_", " ")))) %>% 
      rename(n_pcr = n, pathogen_name = substrain_name) %>% 
      group_by(date, pathogen_name) %>% 
      summarise(n_pcr=sum(n_pcr), .groups="drop")
    
    
  } else {
     pathogen_pcr <- counts_pcr %>% 
      filter(grepl(str_replace_all(pathogen_pcr, "_", " "), substrain_name, ignore.case=TRUE)) %>% 
      select(date, substrain_name, n) %>% 
      mutate(substrain_name = tolower((str_replace_all(pathogen, "_", " ")))) %>% 
      rename(n_pcr = n, pathogen_name = substrain_name)
  }
  path_data<-full_join(pathogen_data, pathogen_pcr, by = c("date", "pathogen_name")) %>% 
    full_join(., detected_data, by=c("date", "pathogen_name")) %>% 
    mutate(week= lubridate::floor_date(date, "week")) %>% 
    full_join(., sent_data %>% select(-date), by=c("week", "pathogen_name")) %>% 
    {if(pathogen %in% wastewater_pathogens) left_join(., wastewater_data %>% select(-date), by=c("week", "pathogen_name")) %>%  
        select(week, pathogen_name, n_seq, n_pcr, n_detected, percentage, seven_day_median_viral_load_avg) else select(.,week, pathogen_name, n_seq, n_pcr, n_detected, percentage)} %>% 
    rename(date=week)
  
  return(path_data)
  
}