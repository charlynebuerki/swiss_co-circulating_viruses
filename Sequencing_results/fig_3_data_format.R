#formatting for data into figure 2 co-infections
#formats the data in 3 columns: strain_1, strain_2 and a count to detect all combinations of co-infections 
format_matrix_co_infection_plot<-function(data)
{
  co <- data  %>% 
    select(pseudonymized_id, substrain_name, is_unique) %>% 
    group_by(pseudonymized_id) %>% 
    filter(!is_unique) %>% 
    summarise(pairs = list(combn(substrain_name, 2, simplify =FALSE)), .groups="drop") %>% 
    unnest(pairs) %>% 
    mutate(strain_1 = map_chr(pairs,1),
           strain_2 = map_chr(pairs,2),
           n=1) %>% 
    select(-pairs) %>% 
    bind_rows(., rename(., strain_1=strain_2, strain_2=strain_1))
  
  single <- expand.grid(strain_1 = data$substrain_name %>% unique(), strain_2 = data$substrain_name %>% unique()) %>% 
    #filter(strain_1 != strain_2) %>% 
    mutate(n=0)
  
  matrix_co <- left_join(single, co, by=c("strain_1", "strain_2")) %>% 
    mutate(n=coalesce(n.y, n.x)) %>% 
    select(strain_1, strain_2, n)
  
  return(matrix_co)
}

#returns dataframe formatted with substrain_name, n_co, n, percent_co, where n_co is the count of co-infections, n is total samples, and percent_co is the co_infection percentage
#expects data to have at the minimum: is_unique column_substrain_name
format_percent_count_co_infection_plot<-function(data)
{
  freq_data <- data %>% 
    group_by(is_unique) %>% 
    dplyr::count(substrain_name) %>% 
    pivot_wider(
      names_from = is_unique,
      values_from = n,
      names_prefix = "n_",
      values_fill = list(n=0)
    ) %>% 
    rename(n_co=n_FALSE, n=n_TRUE) %>% 
    ungroup() %>% 
    mutate(percent_co = (n_co/(n+n_co))*100)
  return(freq_data)
}