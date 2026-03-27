
# test to see which cluster beta to choose 

# test 1
epi_beta_filtered_1 <- epi_1_total %>% filter(mean_abs_delta_detail >= 0.15)  #
epi_beta_filtered_2 <- epi_1_total %>% filter(pct_cpgs_02_detail >=50 )  #
epi_beta_filtered_3 <- epi_1_total %>% filter(abs(median_nSD_detail) >= 10) # needs to be made absolute 
epi_beta_filtered_4_bg <- epi_1_total %>% filter(mean_nSD_detail_bg >= 0.02)
epi_beta_filtered_5_bg <- epi_1_total %>% filter(mean_nSD_detail_bg  < 0.15)
epi_beta_filtered_6_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.03)

# ============================================================
epi_beta_filtered_1 <-epi_beta_filtered_1 %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)
epi_beta_filtered_2 <-epi_beta_filtered_2 %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)
epi_beta_filtered_3 <- epi_beta_filtered_3 %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)
epi_beta_filtered_4_bg <- epi_beta_filtered_4_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)
epi_beta_filtered_5_bg <-epi_beta_filtered_5_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)
epi_beta_filtered_6_bg <-epi_beta_filtered_6_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)

plot_cluster_landscape(
epi_clustered    = epi_beta_filtered_1,    # hybrid_clustering() output
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
epi_clustered    = epi_beta_filtered_2,    # hybrid_clustering() output
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
epi_clustered    = epi_beta_filtered_3,    # hybrid_clustering() output
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




# test 2 
#! I THINK WE SHOULD EITHER CONCENTRATE ON ONLY unvarible in bg or use the others in combination 
# OR MAYBE WE SHOULD JUST USE all 
# for first filtering we will conentrate on unvaribale (kkep only those )
# wewill try all the diffrent mean_nSD_detail_bg by themselves first 
# then by themselfs but with enrich_clusters and filtered for cluster mean delta detail and maybe as well cluster pct first separate and then together 
 # finally we will try 0.03, 0.06. 0.1 and 0.15 in mean_nSD_detail_bg with mean_abs_delta_detail 0.15 and or pct_cpgs_02_detail >=50 and or abs(median_nSD_detail >= 10) before doing the clustering 


epi_beta_filtered_1_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.03)
epi_beta_filtered_2_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.05)
epi_beta_filtered_3_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.08)
epi_beta_filtered_4_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.09)
epi_beta_filtered_5_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.1)
epi_beta_filtered_6_bg <- epi_1_total %>% filter(mean_nSD_detail_bg < 0.13)

epi_beta_filtered_1_bg <-epi_beta_filtered_1_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  
epi_beta_filtered_2_bg <-epi_beta_filtered_2_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  
epi_beta_filtered_3_bg <- epi_beta_filtered_3_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  
epi_beta_filtered_4_bg <- epi_beta_filtered_4_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  
epi_beta_filtered_5_bg <-epi_beta_filtered_5_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  
epi_beta_filtered_6_bg <-epi_beta_filtered_6_bg %>% hybrid_clustering() %>%assign_canonical_cluster_ids() %>% filter(in_canonical_cluster = TRUE)  


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

# ============================================================
# HELPER: reduces the repetitive plot calls
# NOTE: plot_cluster_landscape() doesn't have a 'title' param
# in its current definition — either add one or remove below
# ============================================================

make_plot <- function(epi_df, label, n_clusters = 2, min_epi_pts = 5) {
  cat("\n=== Plotting:", label, "| n regions:", nrow(epi_df), "===\n")
  plot_cluster_landscape(
    epi_clustered    = epi_df,
    master_annotated = master_annotated,
    beta_long_base   = beta_long_base,
    lookup_mega      = lookup_mega,
    cpg_positions    = cpg_positions,
    n_clusters       = n_clusters,
    window_cpgs      = 200,
    cluster_type     = "any",
    min_epi_pts      = min_epi_pts
  )
}

# ============================================================
# HELPER: standard clustering pipeline
# ============================================================

cluster_pipeline <- function(epi_df) {
  epi_df %>%
    hybrid_clustering() %>%
    assign_canonical_cluster_ids() %>%
    filter(in_canonical_cluster == TRUE)   # == not =
}

# ============================================================
# PART A: BG SD THRESHOLDS ALONE
# Logic: keep only near-invariant CpG regions (low background SD)
# so nSD values are not artificially inflated
# ============================================================

bg_thresholds <- c(0.03, 0.05, 0.08, 0.09, 0.10, 0.13)
bg_labels     <- paste0("bg_lt_", gsub("\\.", "", sprintf("%.2f", bg_thresholds)))

partA_results <- list()

for (i in seq_along(bg_thresholds)) {
  
  thresh <- bg_thresholds[i]
  label  <- bg_labels[i]
  
  cat("\n=== PART A |", label, "===\n")
  
  filtered <- epi_1_total %>%
    filter(mean_nSD_detail_bg < thresh)
  
  cat("Regions after bg filter:", nrow(filtered), "\n")
  
  clustered <- cluster_pipeline(filtered)
  
  partA_results[[label]] <- clustered
  
  make_plot(clustered, label)
}


# ============================================================
# PART B: BG SD THRESHOLDS + ENRICH CLUSTERS
#         then filter on post-cluster metrics
#
# B1: filter by cluster_mean_abs_delta_detail alone
# B2: filter by cluster_pct_cpgs_02_detail alone
# B3: filter by both together
# ============================================================

# Post-cluster thresholds to try
cluster_delta_thresh <- 0.15   # cluster_mean_abs_delta_detail
cluster_pct_thresh   <- 50     # cluster_pct_cpgs_02_detail

partB_results <- list()

for (i in seq_along(bg_thresholds)) {
  
  thresh <- bg_thresholds[i]
  label  <- bg_labels[i]
  
  cat("\n=== PART B |", label, "===\n")
  
  # Step 1: pre-cluster bg filter
  filtered <- epi_1_total %>%
    filter(mean_nSD_detail_bg < thresh)
  
  # Step 2: cluster + enrich (enrich adds cluster_mean_abs_delta_detail etc.)
  clustered_enriched <- filtered %>%
    hybrid_clustering() %>%
    assign_canonical_cluster_ids() %>%
    filter(in_canonical_cluster == TRUE) %>%
    enrich_clusters(
      beta_context = beta_context,
      sample_info  = sample_info,
      cluster_col  = "canonical_cluster_id",
      sample_col   = "samples_final"
    )
  
  # B1: cluster effect size only
  b1 <- clustered_enriched %>%
    filter(cluster_mean_abs_delta_detail >= cluster_delta_thresh)
  lbl_b1 <- paste0(label, "_B1_delta")
  partB_results[[lbl_b1]] <- b1
  cat("B1 (delta filter) regions:", nrow(b1), "\n")
  make_plot(b1, lbl_b1)
  
  # B2: cluster % CpGs above 0.2 only
  b2 <- clustered_enriched %>%
    filter(cluster_pct_cpgs_02_detail >= cluster_pct_thresh)
  lbl_b2 <- paste0(label, "_B2_pct")
  partB_results[[lbl_b2]] <- b2
  cat("B2 (pct filter) regions:", nrow(b2), "\n")
  make_plot(b2, lbl_b2)
  
  # B3: both together
  b3 <- clustered_enriched %>%
    filter(
      cluster_mean_abs_delta_detail >= cluster_delta_thresh,
      cluster_pct_cpgs_02_detail    >= cluster_pct_thresh
    )
  lbl_b3 <- paste0(label, "_B3_delta_AND_pct")
  partB_results[[lbl_b3]] <- b3
  cat("B3 (delta AND pct filter) regions:", nrow(b3), "\n")
  make_plot(b3, lbl_b3)
}


# ============================================================
# PART C: COMBINED PRE-CLUSTER FILTERS
# bg threshold × effect size filter combinations
#
# bg thresholds:     0.03, 0.06, 0.10, 0.15
# effect filters:
#   C1  mean_abs_delta_detail >= 0.15
#   C2  pct_cpgs_02_detail    >= 50
#   C3  abs(median_nSD_detail) >= 10   <-- fixed: abs() wraps the column, not the comparison
#   C4  C1 + C2
#   C5  C1 + C3
#   C6  C2 + C3
#   C7  C1 + C2 + C3
# ============================================================

bg_thresholds_c <- c(0.03, 0.06, 0.10, 0.15)
bg_labels_c     <- paste0("bg_lt_", gsub("\\.", "", sprintf("%.2f", bg_thresholds_c)))

partC_results <- list()

for (i in seq_along(bg_thresholds_c)) {
  
  thresh <- bg_thresholds_c[i]
  bg_lbl <- bg_labels_c[i]
  
  # Base: bg filter only
  base <- epi_1_total %>%
    filter(mean_nSD_detail_bg < thresh)
  
  cat("\n=== PART C |", bg_lbl, "| base regions:", nrow(base), "===\n")
  
  # Define each effect filter as a named list of dplyr expressions
  effect_filters <- list(
    C1_delta       = quote(mean_abs_delta_detail >= 0.15),
    C2_pct         = quote(pct_cpgs_02_detail    >= 50),
    C3_nsd         = quote(abs(median_nSD_detail) >= 10),    # fixed
    C4_delta_pct   = quote(mean_abs_delta_detail >= 0.15 & pct_cpgs_02_detail    >= 50),
    C5_delta_nsd   = quote(mean_abs_delta_detail >= 0.15 & abs(median_nSD_detail) >= 10),
    C6_pct_nsd     = quote(pct_cpgs_02_detail    >= 50   & abs(median_nSD_detail) >= 10),
    C7_all_three   = quote(mean_abs_delta_detail >= 0.15 & pct_cpgs_02_detail >= 50 & abs(median_nSD_detail) >= 10)
  )
  
  for (ef_name in names(effect_filters)) {
    
    label <- paste0(bg_lbl, "_", ef_name)
    
    filtered <- base %>% filter(!!effect_filters[[ef_name]])
    
    cat("  ", label, "| regions after effect filter:", nrow(filtered), "\n")
    
    if (nrow(filtered) < 3) {
      cat("  Skipping", label, "— too few regions to cluster\n")
      next
    }
    
    clustered <- cluster_pipeline(filtered)
    
    if (nrow(clustered) == 0) {
      cat("  Skipping", label, "— no canonical clusters formed\n")
      next
    }
    
    cat("  ", label, "| regions in canonical clusters:", nrow(clustered), "\n")
    
    partC_results[[label]] <- clustered
    
    make_plot(clustered, label)
  }
}


# ============================================================
# SUMMARY TABLE: how many regions survive each combination
# ============================================================

summary_rows <- bind_rows(
  
  imap_dfr(partA_results, ~ tibble(
    part    = "A",
    label   = .y,
    n_regions = nrow(.x)
  )),
  
  imap_dfr(partB_results, ~ tibble(
    part    = "B",
    label   = .y,
    n_regions = nrow(.x)
  )),
  
  imap_dfr(partC_results, ~ tibble(
    part    = "C",
    label   = .y,
    n_regions = nrow(.x)
  ))
)

print(summary_rows, n = Inf)
write_tsv(summary_rows, "filter_comparison_summary.tsv")