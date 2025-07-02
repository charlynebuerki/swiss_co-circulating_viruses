#figure plotting for fig 2

#plot alluvial function for plotting; requires number of parameters to distinguish between PCR/seq to seq/seq
plot_alluvial<-function(data, stratum_labels, fill_labels, legend_title, plot_title)
{
  
  fig_all<-ggplot(data %>% 
                    mutate(result= factor(result, levels = c("multi-positive", "single positive", "negative", "no detection"))),
                  aes(x = test, stratum = result, alluvium = alluvium, y = n, fill = is_consistent)) +
    geom_alluvium(width = 1/8, aes(fill = is_consistent),  knot.pos = 0.2) +
    geom_stratum(width = 1/8, alpha = 1, aes(fill=result), color="black", linewidth=0.2) +
    geom_text(stat = "stratum",
              aes(label = paste0(after_stat(stratum), " \n(n=", after_stat(count), ")")),
              size = 3, vjust = 0.5, hjust = 0.5, nudge_x = -0.2) +
    scale_x_discrete(labels = stratum_labels) +
    scale_fill_manual(
      values = c(
        "Partially" = "orange",
        "Yes" = "#609078",
        "No" = "red"
      ),
      labels = fill_labels, # Use the formatted labels for the legend
      na.value = "lightgrey",
      name = legend_title
    ) +
    labs(
      title = plot_title,
      y = "No. Samples",
      x =""
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          panel.grid.major.x = element_blank())
  
  legend <- get_legend(fig_all)
  fig<-fig_all + theme(legend.position = "none")
  
  #return legend and figure separately 
  return(list(plot=fig, legend=legend))
}


# plot pcr to detected
plot_alluvial_pcr_to_detected<-function(pcr_detect_df)
{
  all_alluvial<- pcr_detect_df %>% ungroup() %>% 
    select(pseudonymized_id, result_PCR, result_detect, is_consistent) %>% 
    dplyr::count(result_PCR, result_detect, is_consistent) %>% 
    mutate(result_detect = ifelse(is.na(result_detect), "no detection", result_detect))
  
  
  # 1. Convert to lodes form using to_lodes_form()
  alluvial_lodes <- all_alluvial %>%
    to_lodes_form(
      key = "test",
      axes = c("result_PCR", "result_detect"),
      value = "result",
      key.sort = TRUE # Optional: Sort the axes alphabetically
    ) %>%
    filter(!is.na(result))
  
  # 2. Calculate stratum counts
  stratum_counts <- alluvial_lodes %>%
    group_by(test, result) %>%
    summarise(count = sum(n), .groups = "drop")
  
  # 3. Join stratum counts back to the lodes data 
  alluvial_lodes_with_counts <- alluvial_lodes %>%
    left_join(stratum_counts, by = c("test", "result"))
  
  # 4. Calculate counts per match_type for the legend
  legend_counts <- all_alluvial %>%
    group_by(is_consistent) %>%
    summarise(n_samples = sum(n), .groups = "drop") %>%
    mutate(legend_label = paste0(is_consistent, "\n(n=", n_samples, ")"))
  
  # 5. Create a named vector for the fill scale with the formatted labels
  fill_labels <- setNames(legend_counts$legend_label, legend_counts$is_consistent)
  
  plots<-plot_alluvial(alluvial_lodes_with_counts, stratum_labels = c("PCR", "seq. detected"), fill_labels=fill_labels, legend_title= "Consistent PCR/seq. detection?", plot_title="PCR/seq. detection")
  
  return(plots)
}

#plot detected to hq
plot_alluvial_detected_to_hq<-function(detected_hq_df)
{
  all_alluvial<- detected_hq_df %>% ungroup() %>% 
    select(pseudonymized_id, result_hq, result_detect, is_consistent) %>% 
    dplyr::count(result_hq, result_detect, is_consistent) %>% 
    mutate(result_hq = ifelse(is.na(result_hq), "no detection", result_hq))
  
  # 1. Convert to lodes form using to_lodes_form()
  alluvial_lodes <- all_alluvial %>%
    to_lodes_form(
      key = "test",
      axes = c("result_detect", "result_hq"),
      value = "result",
      key.sort = TRUE # Optional: Sort the axes alphabetically
    ) %>%
    filter(!is.na(result))
  
  # 2. Calculate stratum counts (already done in your provided code)
  stratum_counts <- alluvial_lodes %>%
    group_by(test, result) %>%
    summarise(count = sum(n), .groups = "drop")
  
  # 3. Join stratum counts back to the lodes data (already done in your provided code)
  alluvial_lodes_with_counts <- alluvial_lodes %>%
    left_join(stratum_counts, by = c("test", "result"))
  
  # 4. Calculate counts per match_type for the legend
  legend_counts <- all_alluvial %>%
    group_by(is_consistent) %>%
    summarise(n_samples = sum(n), .groups = "drop") %>%
    mutate(legend_label = paste0(is_consistent, "\n(n=", n_samples, ")"))
  
  # 5. Create a named vector for the fill scale with the formatted labels
  fill_labels <- setNames(legend_counts$legend_label, legend_counts$is_consistent)
  
  plots<-plot_alluvial(alluvial_lodes_with_counts, stratum_labels = c("seq. detected", "seq. hq"), fill_labels=fill_labels, legend_title= "Consistent seq. detect/hq?", plot_title="detection/hq seq.")
  
  return(plots)
  
}

#overall figure to make the alluvial plot
make_figure_comparison<-function(formated_pcr_detected, formated_detected_hq, save=TRUE)
{
  plot_a<-plot_alluvial_pcr_to_detected(formated_pcr_detected)
  plot_b<-plot_alluvial_detected_to_hq(formated_detected_hq)
  
  fig <-plot_grid(plot_a$plot, plot_b$plot, plot_a$legend,  plot_b$legend, nrow=2, ncol = 2 ,rel_heights = c(1,0.15), align="v", axis="l")
  
  if(save) {ggsave("images/final_analysis/fig_pcr_comparison.pdf",fig , dpi="retina", units="mm", width=174, height=123)}
  
  return(fig)
}