#script to get the formatting for supp table x for reporting quantifications of mapped reads and DP10 in supp table. 


#format the data correctly
format_data_table_quantification<-function(df)
{
  df_hq <-df %>% 
    separate_rows(strains_PCR, panels_PCR, values_PCR, sep = ",\\s*") %>%
    distinct(pseudonymized_id, strain_name, .keep_all = TRUE) %>% 
    mutate(
      panel_type = ifelse(grepl("resp", panels_PCR), "resp", ifelse(panels_PCR == "flua" | panels_PCR == "flub" | panels_PCR == "rspc" | panels_PCR== "wuhan0", "p2", "p1"))
    ) %>% 
    group_by(strain_name, panel_type) %>% 
    summarise(
      n = n(), # Calculate n first so we can use it below
      
      mean_reads = mean(aligned),
      # Calculate 95% Margin of Error for Reads
      # qt(0.975) gives the critical value for a two-tailed 95% CI
      me_reads = qt(0.975, df = n() - 1) * (sd(aligned) / sqrt(n())),
      
      mean_dp = mean(DP10),
      # Calculate 95% Margin of Error for DP
      me_dp = qt(0.975, df = n() - 1) * (sd(DP10) / sqrt(n())),
      .groups = "drop"
    )
  return(df_hq)
}


# Helper function to format scientific notation consistently
fmt_sci <- function(x) {
  s <- formatC(x, format = "e", digits = 2)
  gsub("e\\+?0?", "x10^", s)
}

format_data_table_quantification_total<-function(df)
{
  df_hq <-df %>% 
    separate_rows(strains_PCR, panels_PCR, values_PCR, sep = ",\\s*") %>%
    distinct(pseudonymized_id, strain_name, .keep_all = TRUE) %>% 
    group_by(strain_name) %>% 
    summarise(
      n = n(), # Calculate n first so we can use it below
      
      mean_reads = mean(aligned),
      # Calculate 95% Margin of Error for Reads
      # qt(0.975) gives the critical value for a two-tailed 95% CI
      me_reads = qt(0.975, df = n() - 1) * (sd(aligned) / sqrt(n())),
      
      mean_dp = mean(DP10),
      # Calculate 95% Margin of Error for DP
      me_dp = qt(0.975, df = n() - 1) * (sd(DP10) / sqrt(n())),
      .groups = "drop"
    )
  return(df_hq)
  
}


#function that returns the quantifiation (mapped reads and genome coverage by pcr panel employed and strain in reportable format for paper)
make_supplementary_table_quantification_sequencing<-function(df)
{
  
  formatted_df<-format_data_table_quantification(df)
  total_df <- format_data_table_quantification_total(df)
  
  table <- formatted_df %>%
    
    # --- Calculate Total N for the final column ---
    group_by(strain_name) %>%
    mutate(total_n = sum(n)) %>%
    ungroup() %>%
    
    mutate(
      # --- 1. FORMAT READS (Scientific Notation) ---
      # Calculate Lower and Upper bounds
      r_lower = mean_reads - me_reads,
      r_upper = mean_reads + me_reads,
      
      # Format: "Mean (Lower-Upper)"
      formatted_reads = ifelse(
        n > 1,
        paste0(fmt_sci(mean_reads), " (", fmt_sci(r_lower), "-", fmt_sci(r_upper), ")"),
        paste0(fmt_sci(mean_reads), " (NA)") # Handle n=1 case
      ),
      
      # --- 2. FORMAT DP (Standard Rounding 2 decimals) ---
      d_lower = mean_dp - me_dp,
      d_upper = mean_dp + me_dp,
      
      formatted_dp = ifelse(
        n > 1,
        paste0(format(round(mean_dp, 2), nsmall=2), 
               " (", format(round(d_lower, 2), nsmall=2), "-", 
               format(round(d_upper, 2), nsmall=2), ")"),
        paste0(format(round(mean_dp, 2), nsmall=2), " (NA)")
      )
    ) %>%
    
    # --- Pivot and Select ---
    pivot_wider(
      id_cols = c(strain_name, total_n), 
      names_from = panel_type, 
      values_from = c(formatted_reads, formatted_dp, n),
      names_glue = "{.value}_{panel_type}" 
    ) %>%
    
    select(
      strain_name, 
      starts_with("formatted_reads"), 
      starts_with("formatted_dp"), 
      starts_with("n_"), 
      Total_Samples = total_n
    )
  
  final_table <- table %>% 
    left_join(total_df %>% select(strain_name, mean_reads, me_reads, mean_dp, me_dp), by= "strain_name") %>% 
    mutate(
      # --- 1. FORMAT READS (Scientific Notation) ---
      # Calculate Lower and Upper bounds
      r_lower = mean_reads - me_reads,
      r_upper = mean_reads + me_reads,
      
      # Format: "Mean (Lower-Upper)"
      formatted_reads = ifelse(
        Total_Samples > 1,
        paste0(fmt_sci(mean_reads), " (", fmt_sci(r_lower), "-", fmt_sci(r_upper), ")"),
        paste0(fmt_sci(mean_reads), " (NA)") # Handle n=1 case
      ),
      
      # --- 2. FORMAT DP (Standard Rounding 2 decimals) ---
      d_lower = mean_dp - me_dp,
      d_upper = mean_dp + me_dp,
      
      formatted_dp = ifelse(
        Total_Samples > 1,
        paste0(format(round(mean_dp, 2), nsmall=2), 
               " (", format(round(d_lower, 2), nsmall=2), "-", 
               format(round(d_upper, 2), nsmall=2), ")"),
        paste0(format(round(mean_dp, 2), nsmall=2), " (NA)")
      )
    ) %>% select(-c(mean_reads, me_reads, mean_dp, me_dp, r_lower, r_upper,d_lower, d_upper))
  
  return(final_table)
}