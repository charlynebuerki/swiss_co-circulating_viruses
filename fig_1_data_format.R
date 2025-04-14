#file that formats the data for figure 1 


#format for the b plot for frequencies 
format_frequencies_plot <- function(meta_data, all_data, detected_data ,substrain_to_highlight="")
{
  
  df_sampling_freq <- meta_data %>% 
    ungroup() %>% 
    filter(!is.na(strains_PCR)) %>% 
    select(pseudonymized_id, ent_date) %>% 
    mutate(date= as.Date(ent_date)) %>% 
    dplyr::count(date) %>% 
    rename(n_all=n)
  
  #what sequencing detects without hq threshold
  df_sampling_freq_detect <- detected_data %>% 
    mutate(
      pseudonymized_id = make.unique(as.character(pseudonymized_id)),
      date= as.Date(ent_date)
    ) %>% 
    filter(!grepl(".1", pseudonymized_id)) %>% #just get all samples not doubles
    select(pseudonymized_id,date) %>% 
    dplyr::count(date) %>% 
    rename(n_detect=n)
  
  #hq data 
  df_sampling_freq_sub <- all_data %>% 
    mutate(
      pseudonymized_id = make.unique(as.character(pseudonymized_id)),
      date= as.Date(ent_date)
    ) %>% 
    filter(!grepl(".1", pseudonymized_id)) %>% #just get all samples not doubles
    select(pseudonymized_id,date) %>% 
    dplyr::count(date) %>% 
    rename(n_hq=n)
  
  #highlight hq data
  df_sampling_freq_highlight <- all_data %>% 
    mutate(
      pseudonymized_id = make.unique(as.character(pseudonymized_id)),
      date= as.Date(ent_date)
    ) %>% 
    filter(!grepl(".1", pseudonymized_id)) %>%
    filter(substrain_name %in% substrain_to_highlight) %>% #just get highlighting samples
    select(pseudonymized_id,date) %>% 
    dplyr::count(date) %>% 
    rename(n_highlight=n)
  
  df_sampling_freq_all <- df_sampling_freq %>% 
    left_join(df_sampling_freq_sub, by="date") %>% 
    left_join(df_sampling_freq_detect, by="date") %>% 
    left_join(df_sampling_freq_highlight, by="date") %>% 
    mutate(n_hq= ifelse(is.na(n_hq), 0, n_hq),
           n_highlight= ifelse(is.na(n_highlight), 0, n_highlight),
           n_detect = ifelse(is.na(n_detect),0, n_detect)
    ) %>% 
    group_by(date=lubridate::ceiling_date(date, "week", week_start=1)) %>% 
    summarise(across(n_all:n_highlight, sum))
  
  return(df_sampling_freq_all)
  
}

#formats for the main heatmap grid 
format_grid_pathogen_plot <- function(data, pcr=FALSE)
{
  res <-data %>%
    select(pseudonymized_id, ent_date, substrain_name) %>% 
    mutate(date= as.Date(ent_date)) %>% 
    arrange(date) %>% 
    mutate(
      pseudonymized_id = make.unique(as.character(pseudonymized_id)),  # Make pseudonymized_id unique
      pseudonymized_id = factor(pseudonymized_id, levels = pseudonymized_id[order(date)]),
      substrain_name = factor(substrain_name, levels = get_virus_order(pcr)),
      line_positions= as.numeric(factor(substrain_name,levels = get_virus_order(pcr))), #graphing space
      line_positions= line_positions + 0.5,
      line_max = max(line_positions),
      line_second_max = sort(unique(line_positions), decreasing = TRUE)[2],
      line_positions = ifelse(line_positions %in% c(line_max), NA, line_positions),
      highlight_max = line_positions, #highlighting rectangles
      highlight_min = line_positions -1
    ) %>% 
    dplyr::count(date, substrain_name, line_positions, highlight_max, highlight_min)%>% 
    group_by(date=lubridate::ceiling_date(date, "week", week_start=1), substrain_name,line_positions, highlight_max, highlight_min) %>% 
    summarise(across(n, sum)) %>% 
    ungroup()
  
  return(res)
}

#gets sentinella data to format
get_sentinella_data<- function()
{
  sentinella <- read.csv("data/final_analysis/sentinella/data.csv", stringsAsFactors = FALSE)
  
  sentinella_df <- sentinella %>% 
    mutate(
      year = str_extract(temporal, "\\d{4}") %>% as.integer(),
      week = str_extract(temporal, "(?<=W)\\d{2}") %>% as.integer(),
      date = ISOweek::ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1")),
      
      pathogen_name = ifelse(pathogen == "influenza", paste0(pathogen, "_", type), pathogen) #renaming pathogens 
      
    ) %>% 
    filter(date >= as.Date("2023-03-15") & date <= as.Date("2024-07-15")) %>% #only period we're intersted in
    filter(testResult_type == "pcr") %>%  #for now we only care about PCR result
    group_by(date) %>% 
    mutate(
      divisor= value[valueCategory == "samples" & testResult == "all"], #extract total samples for that week
      percentage = ifelse(valueCategory == "detections", 100*(round(value/divisor,3)), NA)
    ) %>% 
    ungroup() %>% 
    filter(pathogen_name != "influenza_all" & valueCategory != "samples") %>%  #removing unspecific influenza
    select(date, pathogen_name, percentage)
  
  #sentinella plot 
  sentinella_df_plt <- sentinella_df %>% 
    filter(pathogen_name !="all") %>%
    mutate(
      line_positions= as.numeric(factor(pathogen_name)), #graphing space #missing level determination 
      line_positions= line_positions + 0.5,
      line_positions = ifelse(line_positions==max(line_positions), NA, line_positions),
      pathogen_name = factor(pathogen_name, levels=c("other", "sars-cov-2", "influenza_B", "influenza_A", 
                                                     "respiratory_syncytial_virus", "adenovirus", "rhinovirus"))
    ) %>% 
    filter(date <= as.Date("2023-06-01") | date >= as.Date("2023-10-15"))
  
  return(sentinella_df_plt)
}
