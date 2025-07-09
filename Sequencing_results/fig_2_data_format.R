#formatting for fig 2 comparison between PCR & sequencing

#ad hoc function
#function to check the unordered match
check_string_match_unordered <- function(string1, string2) {
  if (is.na(string1) || is.na(string2)) {
    return(FALSE) # Return FALSE if either string is NA
  }
  
  # Split the strings by commas and trim whitespace
  elements1 <- str_split(string1, ",")[[1]] %>% str_trim() %>% sort()
  elements2 <- str_split(string2, ",")[[1]] %>% str_trim() %>% sort()
  
  # Check if the sorted vectors of elements are identical
  identical(elements1, elements2)
}

#format the detect data 
format_detect_data<-function(detect_data, pcr_strains)
{
  
  detected_df <-detect_data %>% select(pseudonymized_id, strain_name, substrain_name, is_unique) %>%
    mutate(
      has_additional_virus = case_when( 
        !(strain_name %in% pcr_strains) ~ TRUE
      ),
      result_detect = ifelse(is_unique, "single positive", "multi-positive") 
    ) %>% 
    group_by(pseudonymized_id) %>% 
    summarize(all_detected_strains = paste(unique(strain_name[!is.na(strain_name)]), collapse = ", "), 
              all_detected_substrains = paste(unique(substrain_name[!is.na(substrain_name)]), collapse = ", "),
              num_detected_strains = 1+str_count(all_detected_strains, ","),
              has_additional_virus = ifelse(any(has_additional_virus), TRUE,FALSE),
              result_detect = first(result_detect)
    ) 
  
  return(detected_df)
}

#format the pcr data
format_pcr_data<-function(pcr_data)
{
  pcr_df <- pcr_data %>%
    select(pseudonymized_id, strains_PCR, panels_PCR, values_PCR, is_unique_PCR) %>%
    mutate(is_unique_PCR = ifelse(pseudonymized_id %in% c("2LfHfe", "YRBBfu", "t9UZtV", "CJU9tE", "UQY8X6"), TRUE, is_unique_PCR), #tested on multiple panels 
           values_PCR = ifelse(is.na(values_PCR), "negative", values_PCR),
           panel_type = ifelse(grepl("resp", panels_PCR, ignore.case = TRUE), "RESP", "other"),
           strains_PCR = ifelse(is.na(strains_PCR), "negative", strains_PCR),
           panels_PCR = ifelse(is.na(panels_PCR), "resp0", panels_PCR),
           result_PCR = ifelse(values_PCR == "negative", "negative",
                               ifelse(is_unique_PCR, "single positive", "multi-positive")),
           num_pcr_strains = 1+str_count(strains_PCR, ",")
    ) 
  
  return(pcr_df)
}

#format the hq data
format_hq_data<-function(hq_data, pcr_strains)
{
  
  hq_df <- hq_data %>% select(pseudonymized_id, strain_name, substrain_name, is_unique) %>% 
    mutate(
      has_additional_virus = case_when( 
        !(strain_name %in% pcr_strains) ~ TRUE
      ),
      result_hq = ifelse(is_unique, "single positive", "multi-positive") 
    ) %>% 
    group_by(pseudonymized_id) %>% 
    summarize(all_hq_strains = paste(unique(strain_name[!is.na(strain_name)]), collapse = ", "), 
              all_hq_substrains = paste(unique(substrain_name[!is.na(substrain_name)]), collapse = ", "),
              #should count the number of strains
              num_hq_strains = 1+str_count(all_hq_strains, ","),
              has_additional_virus_hq = ifelse(any(has_additional_virus), TRUE,FALSE),
              result_hq = first(result_hq)
              
    ) 
  
  return(hq_df)
}


#format data for figure 2
#requires the detected and pcr data, returns a dataframe qualifying the matches between the strains detected in each sample
format_pcr_to_detect_data<-function(pcr_data, detect_data, pcr_panel_strains)
{
  pcr_df <- format_pcr_data(pcr_data)
  
  detected_df <-format_detect_data(detect_data, pcr_panel_strains)
  
  all_df <-left_join(pcr_df, detected_df, by = "pseudonymized_id") %>% 
    mutate(all_detected_strains = ifelse(is.na(all_detected_strains), "negative", all_detected_strains))
  
  
  ## second part, annotate the matches to characterize them
  df <-all_df %>% 
    mutate(
      fuzzy_match = 1- stringdist(strains_PCR, all_detected_strains, method="jw") >= 0.8,
      perfect_match =  map2_lgl(strains_PCR, all_detected_strains, check_string_match_unordered),
      detected_match_PCR_in_detect = map2_lgl(strains_PCR, all_detected_strains, ~grepl(.x, .y, fixed = TRUE)),
      detected_match_detect_in_PCR = map2_lgl(all_detected_strains, strains_PCR, ~grepl(.x, .y, fixed = TRUE))
    ) %>% 
    mutate(
      is_consistent = case_when(
        perfect_match ~ "Yes", #no worries there
        !fuzzy_match & !detected_match_detect_in_PCR & !detected_match_PCR_in_detect ~ "No", # total inconsistency
        num_pcr_strains== 1 & detected_match_PCR_in_detect ~ "Yes", # we find a match among multiple 
        num_detected_strains == 1 & detected_match_detect_in_PCR ~ "Partially", # we find one strain among many, lose detection
        fuzzy_match & num_pcr_strains == 1 & num_detected_strains == 1 & !detected_match_detect_in_PCR & !detected_match_PCR_in_detect ~ "Partially", # special case of Inf A to Inf B matching
        num_pcr_strains > num_detected_strains & num_pcr_strains > 1 & detected_match_detect_in_PCR ~ "Partially", # partial detection of the PCR strains
        (num_detected_strains > num_pcr_strains) & (num_pcr_strains > 1) & is.na(has_additional_virus) & (fuzzy_match) ~ "Yes", #fully consistent and we detect additional strains that are included in panel
        (num_detected_strains > num_pcr_strains) & (num_pcr_strains >1) & detected_match_PCR_in_detect & !fuzzy_match ~ "Yes", # fully consistent and we detect additional strains
        num_detected_strains == num_pcr_strains & num_pcr_strains>1  & has_additional_virus ~ "Partially",  #only partial match but detect additional virus
        TRUE ~ NA_character_
      ),
      
      notes = case_when(
        num_pcr_strains== 1 & detected_match_PCR_in_detect ~ "additional strains detected",
        has_additional_virus ~ "non-PCR-panel strain detected",
        num_detected_strains == 1 & detected_match_detect_in_PCR & is_consistent == "Partially" ~ "identify at least one of co-infection",
        fuzzy_match & num_pcr_strains == 1 & num_detected_strains == 1 & !detected_match_detect_in_PCR & !detected_match_PCR_in_detect ~ "consistent in family; inconsistent in strain", #inf A to inf B matching
        num_pcr_strains > num_detected_strains & num_pcr_strains > 1 & detected_match_detect_in_PCR & is_consistent == "Partially"  ~ "identify more than one co-infection", 
        num_detected_strains > num_pcr_strains & num_pcr_strains > 1 & is.na(has_additional_virus) & (fuzzy_match) ~ "consistent and identify additional panel-based viruses", 
        TRUE ~ NA_character_)
    )
  return(df)
}

#format data for figure 2
#requires the detected and hq data, returns a dataframe qualifying the matches between the strains detected in each sample and hq samples
format_detect_to_hq_data<-function(detect_data, hq_data, pcr_strains)
{
  #format detected data
  detected_df <- format_detect_data(detect_data, pcr_strains)
  
  hq_df <- format_hq_data(hq_data, pcr_strains)
  
  
  all_df_hq<-left_join(detected_df, hq_df, by="pseudonymized_id")  %>% 
    mutate(all_hq_strains = ifelse(is.na(all_hq_strains), "negative", all_hq_strains))
  
  
  #second part, characterize the match
  df_hq <-all_df_hq %>% 
    mutate(
      fuzzy_match = 1- stringdist(all_detected_strains, all_hq_strains, method="jw") >= 0.8,
      perfect_match =  map2_lgl(all_detected_strains, all_hq_strains, check_string_match_unordered),
      detected_match_hq_in_detected = map2_lgl(all_hq_strains, all_detected_strains, ~grepl(.x, .y, fixed = TRUE)),
    ) %>% 
    mutate(is_consistent = case_when(
      
      perfect_match ~ "Yes", #no worries there
      
      all_hq_strains=="negative" & !is.na(all_detected_strains) ~ "No", #we drop out all strains
      
      num_hq_strains < num_detected_strains ~ "Partially", #we detect some strains
      
      TRUE ~ NA_character_
    ),
    
    notes = case_when(
      
      TRUE ~ NA_character_)
    )
  
  return(df_hq)
}