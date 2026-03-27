BASIC SCORE METRIC 





**1. Technical Robustness** [0-30 points]
- Detected by **4-5 methods** (not 3)
- Effect size **Δβ ≥ 0.15** (functional threshold)
- **Consistent direction** (all hyper OR all hypo)
- Why: 5 methods agreeing = NOT an artifact

**2. Genomic Location** [0-20 points]
- **Promoter** (20 pts) - direct gene regulation
- **Enhancer/cCRE** (15 pts) - distal regulation
- **CpG island shores** (10 pts) - dynamic regions
- Gene body (3 pts) - unclear mechanism
- near to TSS 
- Why: Promoter methylation has PROVEN mechanism

# **3. Response Specificity** [0-15 points]
# - **R-only** (15 pts) - resistance mechanism (BEST!)
# - **Shared R+MR** (10 pts) - common pathway
# - MR-only (5 pts) - unclear relevance
# - Why: R-specific = drug response marker, not disease marker

**4. Sample Penetrance** [0-10 points]
- **Recurrent** (≥5 samples, 10 pts) - subgroup marker
- Moderate (3-4 samples, 5 pts) - reproducible
- Individual (1-2 samples, 0 pts) - not generalizable
- Why: Single-patient events could be random



calculate_priority_score <- function(epi_df) {
  
  epi_df %>%
    mutate(
      # First, compute n_samples_numeric from samples_final
      n_samples_numeric = sapply(n_samples_final, function(x) {
        if (grepl("R:|MR:", x)) {
          # Format: "R:5 MR:2" - extract numbers and sum
          nums <- str_extract_all(x, "\\d+")[[1]]
          sum(as.numeric(nums))
        } else {
          # Format: comma-separated sample IDs
          length(str_split(x, ",")[[1]])
        }
      }), 

      # 1. Technical robustness (0-30 points)
      score_technical = case_when(
        tier == "TIER 1: All 5 Methods" ~ 30,
        tier == "TIER 2: 4 Methods" ~ 20,
        tier == "TIER 3: 3 Methods" ~ 10,
        TRUE ~ 0
      ),
      
      # 2. Effect size (0-20 points)
      score_effect = case_when(
        mean_abs_delta_detail >= 0.4 ~ 20,  # huge effect
        mean_abs_delta_detail >= 0.3 ~ 15,
        mean_abs_delta_detail >= 0.2 ~ 10,
        TRUE ~ 0
      ),
      # 3. Genomic location (0-20 points)
      score_location_promoter = ifelse(annotation_broad == "Promoter", 20, 0),
      score_location_cre = ifelse(in_cre == "Yes", 15, 0),
      score_location_cpgi_good = ifelse(
        grepl("shores|shelves|islands", cpgi_context) & !grepl("inter", cpgi_context), 
        10, 0
      ),
      score_location_cpgi_mixed = ifelse(
        grepl("shores|shelves|islands", cpgi_context) & grepl("inter", cpgi_context), 
        5, 0
      ),
      score_location_utr = ifelse(annotation_broad %in% c("5'UTR", "3'UTR"), 8, 0),
      score_location_tss = ifelse(abs(distanceToTSS) <= 1500, 5, 0),
      
      # Sum up location scores (but cap at 20 points maximum)
      score_location = pmin(
        score_location_promoter + score_location_cre + 
        score_location_cpgi_good + score_location_cpgi_mixed + 
        score_location_utr + score_location_tss,
        60 # Cap at 20 points max
      ),
      
      # 5. Sample penetrance (0-10 points)
      score_penetrance = case_when(
        n_samples_numeric >= 10 ~ 10,
        n_samples_numeric >= 5 ~ 8,
        n_samples_numeric >= 3 ~ 5,
        TRUE ~ 0
      ),
      
      # TOTAL SCORE (0-120 points)
      priority_score = score_technical + score_effect + score_location +
                      score_penetrance ,
      middle_score = score_technical + score_effect ,
      nd_score = score_technical + score_location,          
      
      # Priority tier
      priority_tier = case_when(
        priority_score >= 90 ~ "TIER 1: Top candidates (≥90)",
        priority_score >= 70 ~ "TIER 2: Strong candidates (70-89)",
        priority_score >= 50 ~ "TIER 3: Moderate candidates (50-69)",
        TRUE ~ "TIER 4: Low priority (<50)"
      )
    ) %>%
    arrange(desc(priority_score))
}

# Apply scoring
epi_scored <- epi_all %>%
  calculate_priority_score()

# Summary
cat("\n=== PRIORITY DISTRIBUTION ===\n")
print(table(epi_scored$priority_tier))

# Top candidates
top_candidates <- epi_scored %>%
  filter(priority_tier == "TIER 1: Top candidates (≥90)") %>%
  select(region_id, gene_symbol, priority_score, score_technical,
         score_effect, score_location, score_specificity,
         score_penetrance, score_clustering, score_biology,
         chromosome, start, end, epivariant_class,
         mean_abs_delta_detail, tier)