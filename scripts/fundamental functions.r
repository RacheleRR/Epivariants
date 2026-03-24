# FUNDAMENTAL FUNCTIONS

  # ============================================================
  # PREREQUISITE: rebuild per-method, per-sample raw data
  # (one row per region x method x sample)
  build_raw_table <- function(df_names, comparison_label) {
    
    map_dfr(df_names, function(df_name) {
      
      df <- get(df_name)
      
      # Create region_id and extract method name
      df %>%
        mutate(
          region_id  = paste(chromosome, start, end, sep = "_"),
          method     = gsub(paste0("epimutations_", comparison_label, "_"), "", df_name) %>%
                      gsub("_allSamples", "", .),
          comparison = comparison_label
        ) %>%
        # Keep only rows with valid coordinates
        filter(
          !is.na(chromosome), chromosome != 0,
          !is.na(start), start != 0,
          !is.na(end),   end   != 0
        ) %>%
        dplyr::select(region_id, chromosome, start, end,
              sample, method, comparison,
              outlier_direction, cpg_ids)
    })
  }


# ============================================================
# CORE FILTERING FUNCTION
# Applies all 3 filters to one comparison at a time
# ============================================================

  filter_epivariants <- function(raw_data, comparison_label, min_methods = 3) {
    
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("Filtering:", comparison_label, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    
    # ---- Step A: Summarise per region ----
    # For each region, collect per-method sample lists and directions
    
    region_summary <- raw_data %>%
      group_by(region_id, chromosome, start, end) %>%
      summarise(
        # How many methods detected this region?
        n_methods       = n_distinct(method),
        methods_list    = list(sort(unique(method))),
        
        # Per-method sample lists (named list: method -> samples)
        per_method_samples = list(
          split(sample, method) %>% lapply(function(s) sort(unique(s)))
        ),
        
        # Per-method directions (named list: method -> directions)
        per_method_dirs = list(
          split(outlier_direction, method) %>%
            lapply(function(d) sort(unique(trimws(d[!is.na(d) & d != ""]))))
        ),
        
        # All unique samples across all methods
        all_samples     = list(sort(unique(sample))),
        n_samples_total = n_distinct(sample),
        
        # All unique cpg_ids
        cpg_ids_all     = paste(unique(unlist(str_split(cpg_ids, ","))), collapse = ","),
        n_cpgs          = n_distinct(unlist(str_split(cpg_ids, ","))),
        
        .groups = "drop"
      )
    
    cat("Regions before filtering:", nrow(region_summary), "\n")
    
    # ---- Filter A: >= min_methods ----
    
    after_A <- region_summary %>%
      filter(n_methods >= min_methods)
    
    cat("After filter A (>=", min_methods, "methods):", nrow(after_A), "\n")
    
    # ---- Filter B: Same samples across all detecting methods ----
    # "Same samples" = intersection of per-method samples equals union
    # i.e., every method found exactly the same set of samples
    
    after_B <- after_A %>%
      rowwise() %>%
      mutate(
        samples_intersect     = list(Reduce(intersect, per_method_samples)),
        samples_union         = list(Reduce(union,     per_method_samples)),
        n_samples_intersect   = length(samples_intersect),
        n_samples_union       = length(samples_union),
        
        # TRUE only if every method agreed on exactly the same samples
        same_samples_strict   = identical(sort(unlist(samples_intersect)),
                                          sort(unlist(samples_union)))
      ) %>%
      ungroup() %>%
      filter(same_samples_strict)
    
    cat("After filter B (same samples across methods):", nrow(after_B), "\n")
    
    # ---- Filter C: All methods agree on direction ----
    # Direction must be uniformly hypermethylation OR uniformly hypomethylation
    # across ALL methods. Mixed = excluded.
    
    after_C <- after_B %>%
      rowwise() %>%
      mutate(
        # Flatten all directions across all methods into one set
        all_directions      = list(unique(unlist(per_method_dirs))),
        n_directions        = length(all_directions),
        
        # Consensus direction: hyper, hypo, or conflicting
        consensus_direction = case_when(
          all(grepl("hyper", unlist(all_directions))) &
            !any(grepl("hypo", unlist(all_directions))) ~ "hypermethylation",
          all(grepl("hypo", unlist(all_directions))) &
            !any(grepl("hyper", unlist(all_directions))) ~ "hypomethylation",
          TRUE ~ "CONFLICTING"
        ),
        
        same_direction = (consensus_direction != "CONFLICTING")
      ) %>%
      ungroup() %>%
      filter(same_direction)
    
    cat("After filter C (same direction across methods):", nrow(after_C), "\n")
    
    # ---- Build clean output table ----
    
    clean <- after_C %>%
      mutate(
        comparison        = comparison_label,
        # Collapse samples and methods to strings for saving
        samples_final     = map_chr(samples_intersect, paste, collapse = ","),
        methods_final     = map_chr(methods_list,      paste, collapse = ","),
        n_samples_final   = map_int(samples_intersect, length),
        
        # Tier by method count
        tier = case_when(
          n_methods == 5 ~ "TIER 1: All 5 Methods",
          n_methods == 4 ~ "TIER 2: 4 Methods",
          n_methods == 3 ~ "TIER 3: 3 Methods",
          TRUE ~ "OTHER"
        )
      ) %>%
      dplyr::select(
        region_id, chromosome, start, end,
        comparison, tier,
        n_methods, methods_final,
        n_samples_final, samples_final,
        n_cpgs, cpg_ids_all,
        consensus_direction
      ) %>%
      arrange(desc(n_methods), desc(n_samples_final), chromosome, start)
    
    # ---- Summary ----
    cat("\n--- CLEAN SET SUMMARY:", comparison_label, "---\n")
    cat("Total clean regions:", nrow(clean), "\n")
    cat("By tier:\n"); print(table(clean$tier))
    cat("By direction:\n"); print(table(clean$consensus_direction))
    cat("Samples per region (median):", median(clean$n_samples_final), "\n")
    cat("CpGs per region (median):",    median(clean$n_cpgs), "\n")
    
    return(clean)
  }

#Unified table 

  # UNIFIED ANALYSIS - One dataframe with all information
  unified_epivariant_analysis <- function(r_dataframes, mr_dataframes) {
      
      # 1. Collect all data
      all_epi_data <- data.frame()
      
      for(df_name in c(r_dataframes, mr_dataframes)) {
          df <- get(df_name)
          df$region_id <- paste(df$chromosome, df$start, df$end, sep = "_")
          df$comparison <- ifelse(grepl("MR_vs_NR", df_name), "MR_vs_NR", 
                                  ifelse(grepl("R_vs_NR", df_name), "R_vs_NR", "Unknown"))
          df$method <- gsub(".*_(beta|manova|mlm|quantile|iForest)_.*", "\\1", df_name)
          
          all_epi_data <- rbind(all_epi_data, 
                              df[, c("region_id", "sample", "comparison", "method", 
                                    "chromosome", "start", "end","outlier_direction","cpg_ids")])
      }
      
      # 2. Analyze each region comprehensively
      region_summary <- all_epi_data %>%
          group_by(region_id, chromosome, start, end) %>%
          summarise(
              # === OVERALL STATS ===
              n_samples_total = n_distinct(sample),
              samples_all = paste(unique(sample), collapse = ","),
              n_methods_total = n_distinct(method),
              methods_all = paste(unique(method), collapse = ","),
              n_comparisons = n_distinct(comparison),
              comparisons_all = paste(unique(comparison), collapse = ","),
              outlier_direction_all = paste(unique(outlier_direction), collapse = ","),
                  cpg_ids_all = paste(unique(cpg_ids), collapse = ","),
                  n_cpgs = length(unique(unlist(strsplit(cpg_ids_all, ",")))),
              
              # === R_vs_NR SPECIFIC ===
              n_samples_R_vs_NR = n_distinct(sample[comparison == "R_vs_NR"]),
              samples_R_vs_NR = paste(unique(sample[comparison == "R_vs_NR"]), collapse = ","),
              n_methods_R_vs_NR = n_distinct(method[comparison == "R_vs_NR"]),
              methods_R_vs_NR = paste(unique(method[comparison == "R_vs_NR"]), collapse = ","),
              outlier_direction_R_vs_NR = paste(unique(outlier_direction[comparison == "R_vs_NR"]), collapse = ","),
              # === MR_vs_NR SPECIFIC ===
              n_samples_MR_vs_NR = n_distinct(sample[comparison == "MR_vs_NR"]),
              samples_MR_vs_NR = paste(unique(sample[comparison == "MR_vs_NR"]), collapse = ","),
              n_methods_MR_vs_NR = n_distinct(method[comparison == "MR_vs_NR"]),
              methods_MR_vs_NR = paste(unique(method[comparison == "MR_vs_NR"]), collapse = ","),
              outlier_direction_MR_vs_NR =  paste(unique(outlier_direction[comparison == "MR_vs_NR"]), collapse = ","),
              .groups = 'drop'
          ) %>%
          mutate(
              # Replace empty strings with NA
              across(starts_with("samples_"), ~ifelse(. == "", NA_character_, .)),
              across(starts_with("methods_"), ~ifelse(. == "", NA_character_, .)),
              
              # === TIER (based on number of methods) ===
              tier = case_when(
                  n_methods_total == 5 ~ "TIER 1: All 5 Methods",
                  n_methods_total == 4 ~ "TIER 2: 4 Methods",
                  n_methods_total == 3 ~ "TIER 3: 3 Methods",
                  n_methods_total == 2 ~ "TIER 4: 2 Methods",
                  n_methods_total == 1 ~ "TIER 5: 1 Method",
                  TRUE ~ "TIER 6: Unknown"
              ),
              
              # === EPIVARIANT TYPE (based on sample count) ===
              epivariant_type = case_when(
                  n_samples_total == 1 ~ "Individual (1 sample)",
                  n_samples_total == 2 ~ "Bi-dividual (2 samples)",
                  n_samples_total >= 3 & n_samples_total <= 10 ~ "Recurring (3-10 samples)",
                  n_samples_total >= 11 & n_samples_total <= 20 ~ "Multi-recurring (11-20 samples)",
                  n_samples_total > 20 ~ "Max-recurring (>20 samples)",
                  TRUE ~ "Unknown"
              ),
              
              # === COMPARISON TYPE ===
              comparison_type = case_when(
                  n_comparisons == 2 ~ "Shared (Both R_vs_NR & MR_vs_NR)",
                  n_comparisons == 1 & grepl("R_vs_NR", comparisons_all) & 
                      !grepl("MR", comparisons_all) ~ "Unique to R_vs_NR",
                  n_comparisons == 1 & grepl("MR_vs_NR", comparisons_all) ~ "Unique to MR_vs_NR",
                  TRUE ~ "Other"
              ),
              
              # === CONFIDENCE SCORE (composite) ===
              confidence_score = n_methods_total * 10 + 
                                ifelse(n_comparisons == 2, 20, 0) +
                                pmin(n_samples_total, 5) * 2
          ) %>%
          arrange(desc(confidence_score), desc(n_methods_total), desc(n_comparisons))
      
      return(region_summary)
  }

# ANNOTATION 


  annotate_regions <- function(region_df) {

    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    library(ChIPseeker)
    library(annotatr)

    gr <- makeGRangesFromDataFrame(
      region_df %>%
        mutate(chr = chromosome,
              start = as.numeric(start),
              end   = as.numeric(end)),
      keep.extra.columns = TRUE,
      seqnames.field     = "chr"
    )
    names(gr) <- region_df$region_id

    # --- 1A: ChIPseeker ---
    txdb       <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peak_anno  <- annotatePeak(
      gr,
      tssRegion  = c(-2000, 200),
      TxDb       = txdb,
      annoDb     = "org.Hs.eg.db",
      overlap    = "all"
    )

    anno_df <- as.data.frame(peak_anno) %>%
      mutate(
        region_id = paste(seqnames, start, end, sep = "_"),
        annotation_broad = case_when(
          grepl("Promoter",   annotation) ~ "Promoter",
          grepl("5' UTR",     annotation) ~ "5'UTR",
          grepl("3' UTR",     annotation) ~ "3'UTR",
          grepl("Exon",       annotation) ~ "Exon",
          grepl("Intron",     annotation) ~ "Intron",
          grepl("Downstream", annotation) ~ "Downstream",
          grepl("Intergenic", annotation) ~ "Intergenic",
          TRUE                            ~ "Other"
        )
      ) %>%
      dplyr::select(region_id, annotation_broad,
            annotation_detail = annotation,
            distanceToTSS,
            gene_id     = geneId,
                  ENSEMBL,
            gene_symbol = SYMBOL,
            gene_name   = GENENAME)

    # --- 1B: CpG island context ---
    annots_cpgi  <- build_annotations(genome = "hg38", annotations = "hg38_cpgs")
    cpgi_hits    <- findOverlaps(gr, annots_cpgi)

    cpgi_df <- tibble(
      region_id    = names(gr)[queryHits(cpgi_hits)],
      cpgi_context = annots_cpgi$type[subjectHits(cpgi_hits)]
    ) %>%
      mutate(cpgi_context = gsub("hg38_cpg_", "", cpgi_context)) %>%
      group_by(region_id) %>%
      # If a region spans multiple contexts, take the most specific
      summarise(
        cpgi_context = paste(sort(unique(cpgi_context)), collapse = ","),
        .groups = "drop"
      )

    # --- 1C: ENCODE cCRE overlap ---
    cre_hits <- findOverlaps(gr, cre_screen)

    cre_df <- tibble(
      region_id = names(gr)[queryHits(cre_hits)],
      cre_type  = mcols(cre_screen)$cCRE_type[subjectHits(cre_hits)],
      cre_id    = mcols(cre_screen)$cCRE_id[subjectHits(cre_hits)]
    ) %>%
      group_by(region_id) %>%
      summarise(
        cre_types = paste(sort(unique(cre_type)), collapse = ","),
        n_cres    = n_distinct(cre_id),
        .groups   = "drop"
      )

    # --- Combine ---
    region_df %>%
      left_join(anno_df, by = "region_id") %>%
      left_join(cpgi_df, by = "region_id") %>%
      left_join(cre_df,  by = "region_id") %>%
      mutate(
        cpgi_context = replace_na(cpgi_context, "open_sea"),
        in_cre       = !is.na(cre_types)
      )
  }

#Epivariant  beta value context 

  build_beta_context <- function(master_annotated) {

    results <- map_dfr(seq_len(nrow(master_annotated)), function(i) {

      row        <- master_annotated[i, ]
      cpgs       <- unique(trimws(str_split(row$cpg_ids_all, ",")[[1]]))
      cpgs_valid <- intersect(cpgs, rownames(beta_mat))
      if (length(cpgs_valid) == 0) return(NULL)

      epi_class <- row$epivariant_class

      # Parse epivariant samples
      if (epi_class == "R_only") {
        epi_samples <- trimws(str_split(row$samples_final, ",")[[1]])
        # Background: NR reference + non-epivariant R samples
        background_samples_detail <- c(nr_samples,
                                setdiff(r_samples, epi_samples))
        background_samples_general   <- c(nr_samples,
                                setdiff(r_samples,  epi_samples),
                                setdiff(mr_samples, epi_samples))                        

      } else if (epi_class == "MR_only") {
        epi_samples <- trimws(str_split(row$samples_final, ",")[[1]])
        # Background: NR reference + non-epivariant MR samples
        background_samples_detail <- c(nr_samples,
                                setdiff(mr_samples, epi_samples))
        background_samples_general   <- c(nr_samples,
                                setdiff(r_samples,  epi_samples),
                                setdiff(mr_samples, epi_samples))                    


      } else {
        # Shared: has both R and MR epivariant samples
        epi_samples_R  <- str_extract(row$samples_final, "(?<=R\\[)[^\\]]+") %>%
                          str_split(",") %>% `[[`(1) %>% trimws()
        epi_samples_MR <- str_extract(row$samples_final, "(?<=MR\\[)[^\\]]+") %>%
                          str_split(",") %>% `[[`(1) %>% trimws()
        epi_samples    <- unique(c(epi_samples_R, epi_samples_MR))
        # Background: NR + non-epivariant R + non-epivariant MR
        background_samples_detail <- c(nr_samples,
                                setdiff(r_samples,  epi_samples),
                                setdiff(mr_samples, epi_samples))
        background_samples_general   <- c(nr_samples,
                                setdiff(r_samples,  epi_samples),
                                setdiff(mr_samples, epi_samples))                  
      }

      epi_samples_valid <- intersect(epi_samples,    colnames(beta_mat))
      bg_samples_detail_valid  <- intersect(background_samples_detail, colnames(beta_mat))
      bg_samples_general_valid  <- intersect(background_samples_general, colnames(beta_mat))

      # Vectorised per-CpG background stats
      beta_bg_dt      <- beta_mat[cpgs_valid, bg_samples_detail_valid,  drop = FALSE]
      beta_bg_general <- beta_mat[cpgs_valid, bg_samples_general_valid, drop = FALSE]

      beta_epi_mat <- beta_mat[cpgs_valid, epi_samples_valid, drop = FALSE]

      bg_mean_dt <- rowMeans(beta_bg_dt, na.rm = TRUE)
      bg_sd_dt   <- apply(beta_bg_dt, 1, sd, na.rm = TRUE)
      bg_n_dt    <- rowSums(!is.na(beta_bg_dt))

      bg_mean_general <- rowMeans(beta_bg_general, na.rm = TRUE)
      bg_sd_general   <- apply(beta_bg_general, 1, sd, na.rm = TRUE)
      bg_n_general    <- rowSums(!is.na(beta_bg_general)) 


      # One row per CpG x epivariant sample
      expand.grid(
        cpg_id    = cpgs_valid,
        sample_id = epi_samples_valid,
        stringsAsFactors = FALSE
      ) %>%
        as_tibble() %>%
        mutate(
          region_id           = row$region_id,
          epivariant_class    = epi_class,
          tier                = row$tier,
          consensus_direction = row$consensus_direction,
          gene_symbol         = row$gene_symbol,
          annotation_broad    = row$annotation_broad,
          cpgi_context        = row$cpgi_context,
          in_cre              = row$in_cre,
          response            = sample_info$response[
                                  match(sample_id, sample_info$sample_id)],

          beta_epivariant     = map2_dbl(cpg_id, sample_id,
                                        ~ beta_mat[.x, .y]),
          beta_bg_mean_detail        = bg_mean_dt[cpg_id],
          beta_bg_sd_detail          = bg_sd_dt[cpg_id],
          n_bg_samples_detail        = bg_n_dt[cpg_id],
          beta_bg_mean_general       = bg_mean_general[cpg_id],
          beta_bg_sd_general         = bg_sd_general[cpg_id],
          n_bg_samples_general       = bg_n_general[cpg_id],
          delta_beta_detail          = beta_epivariant - beta_bg_mean_detail,
          delta_beta_general         = beta_epivariant - beta_bg_mean_general,
          abs_delta_beta_detail      = abs(delta_beta_detail),
          abs_delta_beta_general     = abs(delta_beta_general),
          # How many SDs from background mean?
          n_sd_from_bg_detail        = delta_beta_detail / beta_bg_sd_detail,
          n_sd_from_bg_general       = delta_beta_general / beta_bg_sd_general
        )
    })

    return(results)
  }


#hotspots
# ============================================================
# DBSCAN CLUSTERING
# Finds clusters of any shape based on local density
# ============================================================

library(dbscan)

dbscan_clustering <- function(epi_df, 
                              eps_kb = 10,      # neighborhood radius
                              min_points = 3) {  # min points to form cluster
  
  # For each chromosome separately
  chr_clusters <- map_dfr(unique(epi_df$chromosome), function(chr) {
    
    chr_data <- epi_df %>% filter(chromosome == chr)
    
    if (nrow(chr_data) < min_points) return(NULL)
    
    # Use midpoint positions for clustering
    positions <- (chr_data$start + chr_data$end) / 2
    
    # DBSCAN expects matrix
    pos_matrix <- matrix(positions / 1000, ncol = 1)  # in kb
    
    # Run DBSCAN
    db <- dbscan(pos_matrix, eps = eps_kb, minPts = min_points)
    
    # Add cluster assignments (0 = noise)
    chr_data %>%
      mutate(
        dbscan_cluster = db$cluster,
        in_dbscan_cluster = dbscan_cluster > 0,
        dbscan_cluster_id = ifelse(
          dbscan_cluster > 0,
          paste0("DBS_", chr, "_", str_pad(dbscan_cluster, 3, pad = "0")),
          NA
        )
      )
  })
  
  cat("DBSCAN clustering (eps =", eps_kb, "kb, minPts =", min_points, "):\n")
  cat("  Total epivariants:     ", nrow(chr_clusters), "\n")
  cat("  In DBSCAN clusters:    ", sum(chr_clusters$in_dbscan_cluster), "\n")
  cat("  Noise (not clustered): ", sum(!chr_clusters$in_dbscan_cluster), "\n")
  cat("  Unique clusters:       ", n_distinct(na.omit(chr_clusters$dbscan_cluster_id)), "\n")
  
  list(
    epi_df_clustered = chr_clusters,
    cluster_summary = chr_clusters %>%
      filter(in_dbscan_cluster) %>%
      group_by(dbscan_cluster_id, chromosome) %>%
      summarise(
        n_epivariants = n(),
        span_start = min(start),
        span_end = max(end),
        span_kb = round((span_end - span_start) / 1000, 2),
        .groups = "drop"
      )
  )
}

# ============================================================
# OVERLAP-BASED CLUSTERING
# Cluster epivariants that overlap by ≥ X%
# ============================================================

overlap_clustering <- function(epi_df, 
                               min_overlap_pct = 50,
                               min_cluster_size = 2) {
  
  library(igraph)
  
  # Convert to GRanges
  epi_gr <- epi_df %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = TRUE
    )
  
  # Find all pairwise overlaps
  hits <- findOverlaps(epi_gr, epi_gr)
  
  # Remove self-hits
  hits <- hits[queryHits(hits) != subjectHits(hits)]
  
  # Calculate overlap percentage
  query_idx <- queryHits(hits)
  subject_idx <- subjectHits(hits)
  
  overlap_width <- width(pintersect(epi_gr[query_idx], epi_gr[subject_idx]))
  query_width <- width(epi_gr[query_idx])
  subject_width <- width(epi_gr[subject_idx])
  
  # Overlap as % of smaller region
  overlap_pct <- overlap_width / pmin(query_width, subject_width) * 100
  
  # Keep only significant overlaps
  significant <- overlap_pct >= min_overlap_pct
  
  # Build graph of overlapping regions
  edges <- tibble(
    from = query_idx[significant],
    to = subject_idx[significant],
    overlap_pct = overlap_pct[significant]
  )
  
  if (nrow(edges) == 0) {
    cat("No overlapping regions found at", min_overlap_pct, "% threshold\n")
    return(list(epi_df_clustered = epi_df %>% mutate(overlap_cluster_id = NA)))
  }
  
  # Create graph and find connected components
  g <- graph_from_data_frame(edges, directed = FALSE)
  components <- components(g)
  
  # Map back to original regions
  cluster_membership <- rep(NA, length(epi_gr))
  cluster_membership[as.numeric(names(components$membership))] <- components$membership
  
  # Create cluster summary
  cluster_sizes <- table(cluster_membership)
  cluster_sizes <- cluster_sizes[cluster_sizes >= min_cluster_size]
  
  epi_df_clustered <- epi_df %>%
    mutate(
      overlap_cluster_id = ifelse(cluster_membership %in% names(cluster_sizes),
                                   paste0("OC", str_pad(cluster_membership, 4, pad = "0")),
                                   NA),
      in_overlap_cluster = !is.na(overlap_cluster_id)
    )
  
  cat("Overlap clustering (min_overlap ≥", min_overlap_pct, "%):\n")
  cat("  Total epivariants:     ", nrow(epi_df), "\n")
  cat("  In overlap clusters:   ", sum(epi_df_clustered$in_overlap_cluster, na.rm = TRUE), "\n")
  cat("  Unique clusters:       ", n_distinct(na.omit(epi_df_clustered$overlap_cluster_id)), "\n")
  
  list(
    epi_df_clustered = epi_df_clustered,
    edges = edges
  )
}

# ============================================================
# PROXIMITY-BASED CLUSTERING
# Cluster epivariants that are within X bp of each other
# ============================================================

library(GenomicRanges)

proximity_clustering <- function(epi_df, 
                                 max_distance = 1000,  # 1 kb
                                 min_cluster_size = 2) {
  
  # Convert to GRanges
  epi_gr <- epi_df %>%
    makeGRangesFromDataFrame(
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = TRUE
    )
  
  # Find overlaps within max_distance
  # Use reduce() with max.gap to merge nearby regions
  clusters_gr <- reduce(epi_gr, 
                        min.gapwidth = max_distance,
                        with.revmap = TRUE)
  
  # Get which epivariants belong to which cluster
  cluster_membership <- as.list(mcols(clusters_gr)$revmap)
  
  # Create cluster table
  cluster_df <- tibble(
    prox_cluster_id = paste0("PC", str_pad(1:length(clusters_gr), 4, pad = "0")),
    chr = as.character(seqnames(clusters_gr)),
    cluster_start = start(clusters_gr),
    cluster_end = end(clusters_gr),
    cluster_width = width(clusters_gr),
    n_epivariants = lengths(cluster_membership)
  ) %>%
    filter(n_epivariants >= min_cluster_size)
  
  # Add cluster membership to original data
  cluster_assignments <- map_dfr(1:nrow(cluster_df), function(i) {
    epi_indices <- cluster_membership[[i]]
    tibble(
      region_id = epi_df$region_id[epi_indices],
      prox_cluster_id = cluster_df$prox_cluster_id[i],
      cluster_n = cluster_df$n_epivariants[i],
      cluster_width_kb = round(cluster_df$cluster_width[i] / 1000, 2)
    )
  })
  
  epi_df_clustered <- epi_df %>%
    left_join(cluster_assignments, by = "region_id") %>%
    mutate(
      in_proximity_cluster = !is.na(prox_cluster_id)
    )
  
  cat("Proximity clustering (max_distance =", max_distance, "bp):\n")
  cat("  Total epivariants:        ", nrow(epi_df), "\n")
  cat("  In proximity clusters:    ", sum(epi_df_clustered$in_proximity_cluster), "\n")
  cat("  Unique clusters:          ", n_distinct(na.omit(epi_df_clustered$prox_cluster_id)), "\n")
  cat("  Median cluster size:      ", median(cluster_df$n_epivariants), "\n")
  cat("  Median cluster width (kb):", median(cluster_df$cluster_width / 1000), "\n")
  
  list(
    epi_df_clustered = epi_df_clustered,
    cluster_summary = cluster_df
  )
}

# ============================================================
# HYBRID CLUSTERING PIPELINE
# Combine multiple clustering methods
# ============================================================

hybrid_clustering <- function(epi_df) {
  
  cat("Running hybrid clustering pipeline...\n\n")
  
  # 1. Proximity clustering (physical closeness)
  cat("Step 1: Proximity clustering...\n")
  prox <- proximity_clustering(epi_df, max_distance = 1000, min_cluster_size = 2)
  epi_df <- prox$epi_df_clustered
  
  # 2. Overlap clustering (shared CpGs)
  cat("\nStep 2: Overlap clustering...\n")
  ovlp <- overlap_clustering(epi_df, min_overlap_pct = 10, min_cluster_size = 2)
  epi_df <- ovlp$epi_df_clustered
  

  
  # 4. DBSCAN (density-based)
  cat("\nStep 4: DBSCAN clustering...\n")
  dbs <- dbscan_clustering(epi_df, eps_kb = 10, min_points = 3)
  epi_df <- dbs$epi_df_clustered
  
  # 5. Create unified cluster score
  epi_df <- epi_df %>%
    mutate(
      n_clustering_methods = (
        as.integer(in_proximity_cluster) +
        as.integer(in_overlap_cluster) +
        as.integer(in_dbscan_cluster)
      ),
      is_truly_clustered = n_clustering_methods >= 2  # ← at least 2 methods agree
    )
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("HYBRID CLUSTERING SUMMARY\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Clustering method agreement:\n")
  print(table(epi_df$n_clustering_methods))
  
  cat("\nEpivariants clustered by ≥2 methods:", sum(epi_df$is_truly_clustered), "\n")
  
  return(epi_df)
}



# beta long 
  # ============================================================
  # STEP A: BUILD BETA_LONG_BASE (run once, ~5 min)
  # Raw β values only — no epivariant flags yet
  # ============================================================

  build_beta_long_base <- function(grset,
                                    chr_filter = NULL,
                                    max_cpgs   = NULL) {
    
    library(tidyverse)
    library(SummarizedExperiment)
    
    cat("=== build_beta_long_base() — run once ===\n")
    
    # CpG positions
    cpg_positions <- rowRanges(grset) %>%
      as.data.frame() %>%
      rownames_to_column("cpg_id") %>%
      arrange(seqnames, start)
    
    if (!is.null(chr_filter)) {
      cpg_positions <- cpg_positions %>%
        filter(seqnames %in% chr_filter)
      cat("  Chr filter:", paste(chr_filter, collapse = ", "), "\n")
    }
    
    if (!is.null(max_cpgs)) {
      cpg_positions <- cpg_positions %>%
        group_by(seqnames) %>%
        slice_head(n = max_cpgs) %>%
        ungroup()
      cat("  Capped at", max_cpgs, "CpGs per chromosome\n")
    }
    
    cpg_positions <- cpg_positions %>%
      mutate(cpg_index = row_number())
    
    # Sample info
    sample_info <- colData(grset) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      dplyr::select(sample_id, response = LiResponse)
    
    # Beta matrix → long format
    beta_mat   <- getBeta(grset)
    valid_cpgs <- intersect(cpg_positions$cpg_id, rownames(beta_mat))
    cpg_positions <- cpg_positions %>% filter(cpg_id %in% valid_cpgs)
    
    cat("  CpGs:", nrow(cpg_positions),
        "| Samples:", nrow(sample_info), "\n")
    
    beta_long_base <- beta_mat[cpg_positions$cpg_id, , drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("cpg_id") %>%
      pivot_longer(-cpg_id, names_to = "sample_id", values_to = "beta") %>%
      left_join(cpg_positions %>%
                  dplyr::select(cpg_id, cpg_index, seqnames, start),
                by = "cpg_id") %>%
      left_join(sample_info, by = "sample_id")
    
    cat("  beta_long_base rows:", format(nrow(beta_long_base), big.mark = ","), "\n")
    cat("  No epivariant flags — use get_beta_long() per combination\n")
    cat("=== Done ===\n")
    
    list(
      beta_long_base = beta_long_base,
      cpg_positions  = cpg_positions,
      sample_info    = sample_info
    )
  }


  # ============================================================
  # STEP B: BUILD EPIVARIANT_LOOKUP_MEGA (run once, ~1 min)
  # All region→cpg→sample pairs for ALL tiers and both comparisons.
  # This is the join key used by get_beta_long() later.
  # ============================================================

  build_epivariant_lookup_mega <- function(unified_summary,
                                            cpg_positions   # from build_beta_long_base()
                                          ) {
    
    library(tidyverse)
    
    cat("=== build_epivariant_lookup_mega() ===\n")
    
    
    all_epivariants <- unified_summary %>%
      filter(grepl("R_vs_NR|MR_vs_NR", comparison_type)) %>%
      mutate(
        comparison_simplified = case_when(
          comparison_type == "Shared (Both R_vs_NR & MR_vs_NR)" ~ "Shared",
          grepl("MR_vs_NR", comparison_type) & !grepl("Shared", comparison_type) ~ "MR_vs_NR only",
          grepl("R_vs_NR",  comparison_type) & !grepl("Shared", comparison_type) ~ "R_vs_NR only",
          TRUE ~ "Other"
        )
      ) %>%
      filter(comparison_simplified != "Other")
    
    cat("  Source epivariant regions:", nrow(all_epivariants), "\n")
    
    # ---- Expand R samples ----
    lookup_R <- all_epivariants %>%
      filter(comparison_simplified %in% c("R_vs_NR only", "Shared"),
            !is.na(samples_R_vs_NR), samples_R_vs_NR != "") %>%
      dplyr::select(region_id, cpg_ids_all, samples_R_vs_NR, comparison_simplified) %>%
      mutate(
        cpg_list    = str_split(cpg_ids_all,    ","),
        sample_list = str_split(samples_R_vs_NR, ",")
      ) %>%
      rowwise() %>%
      mutate(grid = list(expand.grid(
        cpg_id    = cpg_list,
        sample_id = sample_list,
        stringsAsFactors = FALSE
      ))) %>%
      ungroup() %>%
      unnest(grid) %>%
      dplyr::select(region_id, cpg_id, sample_id, comparison_simplified) %>%
      mutate(cpg_id    = trimws(cpg_id),
            sample_id = trimws(sample_id),
            response_group = "R") %>%
      distinct()
    
    # ---- Expand MR samples ----
    lookup_MR <- all_epivariants %>%
      filter(comparison_simplified %in% c("MR_vs_NR only", "Shared"),
            !is.na(samples_MR_vs_NR), samples_MR_vs_NR != "") %>%
      dplyr::select(region_id, cpg_ids_all, samples_MR_vs_NR, comparison_simplified) %>%
      mutate(
        cpg_list    = str_split(cpg_ids_all,    ","),
        sample_list = str_split(samples_MR_vs_NR, ",")
      ) %>%
      rowwise() %>%
      mutate(grid = list(expand.grid(
        cpg_id    = cpg_list,
        sample_id = sample_list,
        stringsAsFactors = FALSE
      ))) %>%
      ungroup() %>%
      unnest(grid) %>%
      dplyr::select(region_id, cpg_id, sample_id, comparison_simplified) %>%
      mutate(cpg_id    = trimws(cpg_id),
            sample_id = trimws(sample_id),
            response_group = "MR") %>%
      distinct()
    
    lookup_mega <- bind_rows(lookup_R, lookup_MR) %>%
      # Only keep CpGs that exist in our beta_long_base
      filter(cpg_id %in% cpg_positions$cpg_id) %>%
      distinct()
    
    cat("  Total lookup rows:", format(nrow(lookup_mega), big.mark = ","), "\n")
    cat("  Unique regions:   ", n_distinct(lookup_mega$region_id), "\n")
    cat("  Unique CpGs:      ", n_distinct(lookup_mega$cpg_id),    "\n")
    cat("=== Done ===\n")
    
    lookup_mega
  }


  # ============================================================
  # STEP C: get_beta_long()  — FAST, per combination
  # Attaches epivariant flags from the pre-built lookup slice
  # ============================================================
  # Takes:
  #   epi_df              — any output of get_epivariants()
  #   beta_long_base      — from build_beta_long_base()$beta_long_base
  #   lookup_mega         — from build_epivariant_lookup_mega()
  #
  # Returns beta_long with flags for this specific combination only.
  # Typical runtime: a few seconds (just a filter + left_join).
  # ============================================================

  get_beta_long <- function(epi_df,
                            beta_long_base,
                            lookup_mega) {
    
    epi_region_ids <- unique(epi_df$region_id)
    
    # Slice lookup to only this combination's regions
    lookup_slice <- lookup_mega %>%
      filter(region_id %in% epi_region_ids) %>%
      mutate(epivariant = TRUE)
    
    # Join onto base (only flag rows that match — everything else = FALSE)
    beta_long <- beta_long_base %>%
      left_join(
        lookup_slice %>%
          dplyr::select(cpg_id, sample_id, epivariant,
                        region_id, response_group, comparison_simplified),
        by = c("cpg_id", "sample_id")
      ) %>%
      mutate(
        epivariant            = replace_na(epivariant, FALSE),
        region_id             = ifelse(is.na(region_id),             NA, region_id),
        response_group        = ifelse(is.na(response_group),        NA, response_group),
        comparison_simplified = ifelse(is.na(comparison_simplified), NA, comparison_simplified),
        # plot_group and epivariant_type derived once here
        plot_group = case_when(
          response == "NR"                                         ~ "NR",
          response == "R"  & !epivariant                          ~ "R (no epivar)",
          response == "MR" & !epivariant                          ~ "MR (no epivar)",
          response == "R"  & epivariant & response_group == "R"   ~ "R - EPIVARIANT",
          response == "MR" & epivariant & response_group == "MR"  ~ "MR - EPIVARIANT",
          response == "R"  & epivariant & response_group == "MR"  ~ "R - MR EPIVARIANT",
          response == "MR" & epivariant & response_group == "R"   ~ "MR - R EPIVARIANT",
          TRUE ~ as.character(response)
        ),
        epivariant_type = case_when(
          epivariant & response_group == "R"  ~ "R Epivariant",
          epivariant & response_group == "MR" ~ "MR Epivariant",
          TRUE ~ NA_character_
        )
      )
    
    cat("get_beta_long(): regions =", length(epi_region_ids),
        "| epivariant points =",
        format(sum(beta_long$epivariant, na.rm = TRUE), big.mark = ","), "\n")
    
    beta_long
  }

# beta context filtering (NOT 100% sure if this is correct or not ) 
  # ============================================================
  # STEP D: get_beta_context() — FAST, per combination
  # Filters the pre-built mega beta_context by region_id
  # ============================================================

  get_beta_context <- function(epi_df,
                                beta_context_mega) {
    
    epi_region_ids <- unique(epi_df$region_id)
    
    bc <- beta_context_mega %>%
      filter(region_id %in% epi_region_ids)
    
    cat("get_beta_context(): regions =", n_distinct(bc$region_id),
        "| rows =", format(nrow(bc), big.mark = ","), "\n")
    
    bc
  }
