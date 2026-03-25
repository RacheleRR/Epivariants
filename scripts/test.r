
epi_beta_filtered_1_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.03)
epi_beta_filtered_2_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.05)
epi_beta_filtered_3_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.08)
epi_beta_filtered_4_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.09)
epi_beta_filtered_5_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.1)
epi_beta_filtered_6_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.13)

epi_beta_filtered_1_bg <-epi_beta_filtered_1_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  
epi_beta_filtered_2_bg <-epi_beta_filtered_2_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  
epi_beta_filtered_3_bg <- epi_beta_filtered_3_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  
epi_beta_filtered_4_bg <- epi_beta_filtered_4_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  
epi_beta_filtered_5_bg <-epi_beta_filtered_5_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  
epi_beta_filtered_6_bg <-epi_beta_filtered_6_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(canonical_cluster_id = TRUE)  


 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_1_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_1"
)

 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_2_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_2"
)
 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_3_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_3"
)

 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_4_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_4_bg"
)

 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_5_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_5_bg"
)

 plot_cluster_landscape(
  epi_clustered    = epi_beta_filtered_6_bg,    # hybrid_clustering() output
  master_annotated = master_annotated,
  beta_long_base   = beta_long_base,    # already built -- passed in, not rebuilt
  lookup_mega      = lookup_mega,       # already built -- passed in, not rebuilt
  cpg_positions    = cpg_positions,     # already built -- passed in, not rebuilt
  n_clusters       = 2,
  window_cpgs      = 200,
  cluster_type     = "any",             # "any" | "proximity" | "overlap" | "dbscan"
   min_epi_pts   = 5,
   title = "epi_beta_filtered_6_bg"
)

