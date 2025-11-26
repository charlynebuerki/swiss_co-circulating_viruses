#data formatting for the supplementary CT / qc on mapped reads formatting 

#format the data needed for the correlation plots between aligned reads, DP10 and CT values; requires df which is the hq data
format_ct_value_data<-function(df)
{
  ct_data_hq <-df %>% filter(panels_PCR %in% c("flub4p", "flua4p", "rsv4p") & coverage_status=="PASSED") %>% 
    select(pseudonymized_id, substrain_name, strain_name, strains_PCR ,aligned, DP10, values_PCR, is_unique_PCR, is_unique, match_PCR_sequencing)
  
  # separate by virus detected by sequencing, classify others 
  strains_to_classify <- c("Influenza B virus (B/Brisbane/60/2008)",  "Influenza B virus (B/Washington/02/2019)", 
                           "Influenza A virus (A/Michigan/45/2015(H1N1))", "Influenza A virus (A/Puerto Rico/8/1934(H1N1))",
                           "Influenza A virus (A/New York/392/2004(H3N2))",
                           "Human Respiratory syncytial virus 9320 (type B)", "Respiratory syncytial virus (type A)"
  )
  desired_order <- c(
    "Influenza A/H1N1", 
    "Influenza A/H3N2", 
    "Influenza B", 
    "RSV-A", 
    "RSV-B", 
    "other")
  
  final <- ct_data_hq %>% mutate( 
    strain_to_plot = ifelse(substrain_name %in% strains_to_classify, substrain_name, "other"),
    strain_to_plot = case_when(
      grepl("H1N1", strain_to_plot)   ~ "Influenza A/H1N1",
      grepl("H3N2", strain_to_plot)   ~ "Influenza A/H3N2",
      grepl("type B", strain_to_plot) ~ "RSV-B",
      grepl("type A", strain_to_plot) ~ "RSV-A",
      grepl("B/", strain_to_plot)     ~ "Influenza B",
      TRUE                            ~ strain_to_plot  # Keeps "other" as "other"
    ),
    strain_to_plot = factor(strain_to_plot, levels = desired_order),
    aligned = as.numeric(aligned),
    values_PCR = as.numeric(values_PCR)
  )
  
  return(final)
}