#functions to format the data for coverage plots by pathogen

#get a full folder path based on the barcode, pseudo_id and base filepath
get_full_folder_path<-function(barcode, pseudonymized_id, base_filepath="Data/data/depth_files/") 
{
  return(paste0(base_filepath, barcode, "/", pseudonymized_id, "/depth/", pseudonymized_id, "_depth.tsv"))
}

#read the tsv files based on the given filepath 
read_tsv_file <- function(path) {
  if (file.exists(path)) {
    readr::read_tsv(path, col_names = c("virus_ref", "idx", "depth"), show_col_types = FALSE) %>% 
      mutate(source=path,
             idx = as.numeric(idx),
             depth= as.numeric(depth)) # Read file if it exists
  } else {
    tibble(source = path, idx = NA_real_, depth= NA_real_, message = NA_character_)  # Return NA if file is missing
  }
}


#given a virus reference, read from the lookup table and translate to substrain name
get_virus_reference_conversion_table<-function()
{
  virus_lookup<-read.table("Data/resources/virus_lookup_table.bed", sep="\t", header = FALSE)
  names(virus_lookup) <- c("accession_number", "start", "end", "substrain_name")
  
  return(virus_lookup)
}


#convert to virus name from the given virus reference according to the given conversion table 
convert_virus_reference_to_virus_name<-function(virus_reference, conversion_table)
{
  result <-conversion_table %>% 
    filter(grepl(virus_reference, accession_number)) %>%
    pull(substrain_name)
  
  if (length(result) == 0) return(NA_character_)  # Return NA if no match is found
  
  return(result[1])
}


#fetches gff3 files and reads them in correct format
format_gff3_data<-function(virus_reference, filepath= "../Nextstrain/z_gff_files/")
{
  full_address <- paste0(filepath, virus_reference, "_annotation.gff3")
  
  gff_df <-readr::read_tsv(full_address, comment= '#', col_names=FALSE)
  colnames(gff_df) <- c("ref", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  
  gff_df <- gff_df %>% 
    separate_rows(attributes, sep= ";") %>% 
    separate(attributes, into = c("key", "value"), sep = "=", extra = "merge") %>% # Split into key, value
    pivot_wider(names_from = key, values_from = value)  # Reshape into columns
  
  if(virus_reference == "MN908947.3") {
    gff_df <- gff_df %>% rename(Name=gene_name)
  }
  
  return(gff_df)
}


#coverage plots; requires data to have pseuodnymized column and a barcode (plate number) for file lookup purposes as well as an "ent_date" for date of sampling; requires the name of the pathogen to seek
format_coverage_plot_data<- function(data, virus_strain)
{
  #fetch the lookup table
  virus_reference_table <- get_virus_reference_conversion_table()
  
  df<-data %>% 
    select(pseudonymized_id, barcode) %>% 
    mutate(
      pseudonymized_id_bis = substr(pseudonymized_id, start=4, stop=nchar(pseudonymized_id)),
      full_coverage_address= get_full_folder_path(substr(barcode, start=1, stop=nchar(barcode)-3), pseudonymized_id_bis)
           ) 
  
  result <-df %>% pull(full_coverage_address) %>% 
    map_df(read_tsv_file) %>% 
    mutate(pseudonymized_id = str_sub(source, start=-16, end=-11)) %>% 
    group_by(virus_ref) %>% 
    mutate(virus_name = convert_virus_reference_to_virus_name(virus_ref, virus_reference_table)) %>% 
    filter(grepl(virus_strain, virus_name, fixed = TRUE)) %>% 
    ungroup() %>%  #no we reorder based on date 
    left_join(., data %>% select(pseudonymized_id, ent_date), by="pseudonymized_id") %>% 
    mutate(ent_date = as.Date(ent_date)) %>% 
    arrange(ent_date)
  
  return(result)
}


#--------------------------------
#can run: #get data
#virus_strain <- "SARS-CoV-2"
#test_data <-  seq_data %>% filter(grepl(virus_strain, substrain_name))

#df <-format_coverage_plot_data(test_data, virus_strain)