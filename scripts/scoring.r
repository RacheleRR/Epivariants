# FOR SCORING 
we will use genes that are direct target or transporter of lithium 

## DrugBank
# URL: https://go.drugbank.com/
# Search: "lithium" or "lithium carbonate"



WE will use the genes associated to lithium response found by gwas 
# URL: https://www.ebi.ac.uk/gwas/
# Search: "lithium response"


# URL: https://www.pharmgkb.org/
# Search: "lithium"


## DisGeNET
# URL: https://www.disgenet.org/
# Search: "bipolar disorder" or "lithium response"

## gene ontology
#search lithium 
# https://amigo.geneontology.org/amigo/term/GO:0010226


We will use genes asscoiated to neurodrugs in general or psychodrugs  in general 
# URL: https://www.pharmgkb.org/
# Search: "antipsychotics" "Psychoanaleptics" "antidepressants" "Psycholeptics"


We will use genes asscoiated to functional impact oflithium or disorder 
# we will use genes asscoiated to neurotransmitters
# we will use genes associated to neuroprotection and plasticity
# Or maybe we wil use genes asscoiated to neuro and synpases in general 


# we will use genes associated to circadian rhythm because lithium is known to affect circadian rhythm and this is thought to be one of the mechanisms of action of lithium
# GO:0007623

# we will use genes asscoiated to mitocondria function 
# GO:0140053

# we will use genes associated to oxidative stress 
# GO:0006979

# we will use genes asscoiated to signaling patways because second messengers are important for lithium response
# GO:0023052 signaling


We will use genes associated to bipolar disorder 
# found in bipex bd 


epi_truly_clustered <- epi_clustered %>%
  filter(n_clustering_methods >= 2)
```

This reduces false positives exponentially!

---

## **📋 PROBLEM 2: Prioritizing 50,000 Epivariants**

### **The Framework: 7 Lines of Evidence**

**Core Principle:**
> Epivariants most likely related to lithium response are those where **MULTIPLE independent lines of evidence converge**.

### **The 7 Criteria:**

**1. Technical Robustness** [0-30 points]
- Detected by **4-5 methods** (not 3)
- Effect size **Δβ ≥ 0.2** (functional threshold)
- **Consistent direction** (all hyper OR all hypo)
- Why: 5 methods agreeing = NOT an artifact

**2. Genomic Location** [0-20 points]
- **Promoter** (20 pts) - direct gene regulation
- **Enhancer/cCRE** (15 pts) - distal regulation
- **CpG island shores** (10 pts) - dynamic regions
- Gene body (3 pts) - unclear mechanism
- Why: Promoter methylation has PROVEN mechanism

**3. Response Specificity** [0-15 points]
- **R-only** (15 pts) - resistance mechanism (BEST!)
- **Shared R+MR** (10 pts) - common pathway
- MR-only (5 pts) - unclear relevance
- Why: R-specific = drug response marker, not disease marker

**4. Sample Penetrance** [0-10 points]
- **Recurrent** (≥5 samples, 10 pts) - subgroup marker
- Moderate (3-4 samples, 5 pts) - reproducible
- Individual (1-2 samples, 0 pts) - not generalizable
- Why: Single-patient events could be random

**5. True Clustering** [0-10 points]
- **≥3 methods** (10 pts) - strong cluster
- 2 methods (7 pts) - good cluster
- 1 method (3 pts) - weak
- Why: Clustered = coordinated regulation = functional

**6. Biological Relevance** [0-15 points]
- **Direct lithium target** (15 pts) - GSK3B, IMPA1
- **GWAS hit** (12 pts) - CACNA1C, ANK3
- **Neurotransmitter** (10 pts) - SLC6A4, DRD2
- Any gene (5 pts)
- Why: Known pathway = mechanistic link

**7. Pathway Convergence**
- Multiple genes in **same pathway** = stronger signal
- GSK3B + AKT1 + CREB1 = pathway disrupted
- Why: Independent validation

### **The Filtering Cascade:**
```
50,000 epivariants
   ↓ [Technical: 4-5 methods, Δβ≥0.2]
10,000 epivariants
   ↓ [Location: Promoter/enhancer/shore]
2,000 epivariants
   ↓ [Specificity: R-only or Shared]
800 epivariants
   ↓ [Penetrance: ≥2 samples]
400 epivariants
   ↓ [Clustering: ≥2 methods agree]
120 epivariants
   ↓ [Integrated Score ≥90]
20-30 TOP CANDIDATES for validation

To reviewers/skeptics:
Skeptic: "You're just cherry-picking candidates."
Response: "No, we applied 7 objective criteria, each with biological rationale. We can provide the full scoring breakdown for any region."

Skeptic: "Effect size alone should be enough."
Response: "Effect size can be inflated by technical artifacts or cell-type heterogeneity. We require convergence of multiple independent criteria."

Skeptic: "Why exclude gene body epivariants?"
Response: "The mechanism linking gene body methylation to expression is weak/unclear. We prioritize regions with proven regulatory mechanisms (promoters, enhancers)."

Skeptic: "Why is R-only more important than shared?"
Response: "R-only epivariants are specific to responders, suggesting a resistance mechanism. Shared epivariants could be disease-related rather than response-related."

Skeptic: "Why weight GSK3B so highly?"
Response: "GSK3B is the primary known target of lithium. Finding methylation changes in its promoter provides a direct mechanistic link."


# ============================================================
# INTEGRATED PRIORITY SCORING
# ============================================================

calculate_priority_score <- function(epi_df) {
  
  epi_df %>%
    mutate(
      # 1. Technical robustness (0-30 points)
      score_technical = case_when(
        tier == "TIER 1: All 5 Methods" ~ 30,
        tier == "TIER 2: 4 Methods" ~ 20,
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
      score_location = case_when(
        annotation_broad == "Promoter" ~ 20,
        in_cre == "Yes" ~ 15,
        cpgi_context %in% c("shores", "shelves") ~ 10,
        annotation_broad %in% c("5'UTR", "3'UTR") ~ 8,
        TRUE ~ 0
      ),
      
      # 4. Response specificity (0-15 points)
      score_specificity = case_when(
        epivariant_class == "R_only" ~ 15,
        epivariant_class == "Shared_R_and_MR" ~ 10,
        epivariant_class == "MR_only" ~ 5,
        TRUE ~ 0
      ),
      
      # 5. Sample penetrance (0-10 points)
      score_penetrance = case_when(
        n_samples_numeric >= 10 ~ 10,
        n_samples_numeric >= 5 ~ 8,
        n_samples_numeric >= 3 ~ 5,
        TRUE ~ 0
      ),
      
      # 6. Clustering (0-10 points)
      score_clustering = case_when(
        n_clustering_methods >= 3 ~ 10,
        n_clustering_methods == 2 ~ 7,
        n_clustering_methods == 1 ~ 3,
        TRUE ~ 0
      ),
      
      # 7. Biological relevance (0-15 points)
      score_biology = case_when(
        gene_symbol %in% lithium_direct ~ 15,
        gene_symbol %in% gwas_bipolar ~ 12,
        gene_symbol %in% c(neurotrans, neuroprotect) ~ 10,
        !is.na(gene_symbol) ~ 5,  # any gene better than intergenic
        TRUE ~ 0
      ),
      
      # TOTAL SCORE (0-120 points)
      priority_score = score_technical + score_effect + score_location +
                      score_specificity + score_penetrance + score_clustering +
                      score_biology,
      
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