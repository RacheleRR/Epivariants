ANNOTATION PLOTS 

Annotation plots ----
plot_package_annotation <- function(epi_df, label = "epivariants") {
  library(dplyr)
library(stringr)

  simplify_direction <- function(val) {
    # If no "|" assume it's a single direction
    if (!str_detect(val, "\\|")) {
      return(val)
    }
    
    # Split by "|" and get part after ":"
    directions <- str_split(val, "\\|")[[1]] %>%
      str_trim() %>%
      str_extract("(?<=:).*") %>%
      str_to_lower()
    
    # Normalize to 'hypo' or 'hyper'
    normalized <- ifelse(str_detect(directions, "hypo"), "hypo",
                        ifelse(str_detect(directions, "hyper"), "hyper", directions))
    
    # If all same → return that direction
    if (length(unique(normalized)) == 1) {
      if (unique(normalized) == "hypo") return("hypomethylation")
      if (unique(normalized) == "hyper") return("hypermethylation")
    } else {
      return("bidirectional")
    }
  }

# Apply to your dataframe
epi_df <- epi_df %>%
  mutate(simplified_direction = sapply(consensus_direction, simplify_direction))


  # Check annotation columns exist
  if (!"annotation_broad" %in% names(epi_df)) {
    stop("No annotation columns found — run get_epivariants() with master_full")
  }
  
  plots <- list()
  
  # 2a. Genomic feature pie/bar
  plots$feature_bar <- epi_df %>%
    dplyr::count(epivariant_class, annotation_broad) %>%
    group_by(epivariant_class) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    ggplot(aes(x = epivariant_class, y = pct, fill = annotation_broad)) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste("Genomic feature —", label),
         x = NULL, y = "%", fill = "Feature") +
    theme_minimal(base_size = 11)
  
  # 2b. CpG island context
  plots$cpgi_bar <- epi_df %>%
    dplyr::count(epivariant_class, cpgi_context) %>%
    group_by(epivariant_class) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup() %>%
    ggplot(aes(x = epivariant_class, y = pct, fill = cpgi_context)) +
    geom_col(position = "stack", color = "white", linewidth = 0.3) +
    scale_fill_viridis_d() +
    labs(title = paste("CpG island context —", label),
         x = NULL, y = "%", fill = "Context") +
    theme_minimal(base_size = 11)
  
  # 2c. CRE overlap
  plots$cre_bar <- epi_df %>%
    dplyr::count(epivariant_class, in_cre) %>%
    ggplot(aes(x = epivariant_class, y = n, fill = in_cre)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("TRUE"="#E31A1C","FALSE"="#BDBDBD"),
                      labels = c("TRUE"="In CRE","FALSE"="Not in CRE")) +
    labs(title = paste("ENCODE cCRE overlap —", label),
         x = NULL, y = "Count", fill = NULL) +
    theme_minimal(base_size = 11)
  
  # 2d. Feature × direction
  plots$feature_direction <- epi_df %>%
    dplyr::count(annotation_broad, simplified_direction) %>%
    ggplot(aes(x = reorder(annotation_broad, -n),
               y = n, fill = simplified_direction)) +
    geom_col(position = "dodge") +
    geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
    scale_fill_manual(values = c("hypermethylation"="#D73027",
                                 "hypomethylation" ="#4575B4",
                                 "bidirectional" = "#984EA3")) +
    labs(title = paste("Feature × Direction —", label),
         x = NULL, y = "Count", fill = "Direction") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
  combined <- wrap_plots(plots, ncol = 2)
  ggsave(paste0("annotation_", gsub(" ", "_", label), ".png"),
         combined, width = 16, height = 10, dpi = 300)
  
  cat("Saved: annotation_", label, ".png\n")
  invisible(plots)
}
