epivariants v2 
#==========================================
#!create the candindate regions 
#=========================================

library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(GenomicRanges)
library(IRanges)

## Load probe annotation
data(Locations, package = "IlluminaHumanMethylationEPICv2anno.20a1.hg38")

## Select autosomal CpGs
locs.epicv2 <- subset(
  Locations,
  grepl("^cg", rownames(Locations)) &
    chr %in% paste0("chr", 1:22)
)

## Convert to GRanges
locs.epicv2GR <- makeGRangesFromDataFrame(
  locs.epicv2,
  start.field = "pos",
  end.field   = "pos",
  strand      = "*",
  keep.extra.columns = FALSE
)

locs.epicv2GR <- sort(locs.epicv2GR)

## Dummy methylation matrix
mat <- matrix(
  0,
  nrow = length(locs.epicv2GR),
  ncol = 2,
  dimnames = list(names(locs.epicv2GR), c("A", "B"))
)

## Set sample B to all 1
mat[, 2] <- 1

## Design matrix
pheno <- data.frame(var = c(0, 1))
model <- model.matrix(~ var, pheno)

## Run bumphunter
bumps <- bumphunter(
  mat,
  design = model,
  pos = start(locs.epicv2GR),
  chr = as.character(seqnames(locs.epicv2GR)),
  cutoff = 0.05
)$table

## Filter regions
bumps.fil <- subset(bumps, L >= 3)

#create GRanges object with candidate regions

# 1) Basic GRanges creation -------------------------------------------------
candRegsGR <- GRanges(
  seqnames = bumps.fil$chr,
  ranges   = IRanges(start = bumps.fil$start, end = bumps.fil$end),
  strand   = Rle("*", nrow(bumps.fil))
)

# 2) Set stable, unique names (mimic EH names like "chr6_32128101")
names(candRegsGR) <- paste0(as.character(seqnames(candRegsGR)), "_", start(candRegsGR))

# 3) Copy metadata columns from bumps.fil into mcols (keep order similar to epimutacions)
meta_cols <- c("value", "area", "cluster", "indexStart", "indexEnd", "L", "clusterL")
# keep only the ones actually present in bumps.fil to avoid errors
meta_cols <- meta_cols[meta_cols %in% colnames(bumps.fil)]
mcols(candRegsGR) <- as.data.frame(bumps.fil[, meta_cols, drop = FALSE], stringsAsFactors = FALSE)

# 4) Sort (optional) and inspect
candRegsGR <- sort(candRegsGR)
candRegsGR
str(mcols(candRegsGR))

#downloaded the cis-Regulatory Elements file from https://screen.encodeproject.org/
GRCh38.cCREs <- read.delim2("~/Downloads/GRCh38-cCREs.bed", header=FALSE)


colnames(GRCh38.cCREs)[1:3] <- c("chr", "start", "end")

cre_screen <- makeGRangesFromDataFrame(
  GRCh38.cCREs,
  seqnames.field = "chr",
  start.field    = "start",
  end.field      = "end",
  keep.extra.columns = TRUE
)

colnames(mcols(cre_screen))[1:3] <- c("DNase_CREs","cCRE_id","cCRE_type" )


# overlap
ov <- findOverlaps(candRegsGR, cre_screen)

if (length(ov) > 0) {
  # paste overlapping CRE IDs and types for each candidate region
  cre_ids_by_query <- tapply(mcols(cre_screen)$cCRE_id[subjectHits(ov)],
                             queryHits(ov), paste, collapse = ",")
  cre_types_by_query <- tapply(mcols(cre_screen)$cCRE_type[subjectHits(ov)],
                               queryHits(ov), paste, collapse = ",")

  # initialize columns
  mcols(candRegsGR)$CRE <- NA_character_
  mcols(candRegsGR)$CRE_type <- NA_character_

  # assign by index (names of tapply are queryHits as character)
  qidx <- as.integer(names(cre_ids_by_query))
  mcols(candRegsGR)$CRE[qidx] <- as.character(cre_ids_by_query)
  mcols(candRegsGR)$CRE_type[qidx] <- as.character(cre_types_by_query)
} else {
  message("No overlaps found between candRegsGR and cre_screen.")
}


#! canidateregion file candRegsGR !!!

#=========================================
# Load and prepare data 
#=========================================

library(SummarizedExperiment)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

load("~/epivariants/data/beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi.RData")
pd_Li_617 <- read.csv("~/epivariants/data/pd_Li_617.csv")

dim(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi)
dim(pd_Li_617)

# 1. Get probe annotation
data(Locations)
locs <- Locations[rownames(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi), ]
locs <- locs[!is.na(locs$pos) & !is.na(locs$chr), ]
#529 

cpg_gr <- makeGRangesFromDataFrame(
  locs,
  seqnames.field = "chr",
  start.field = "pos",
  end.field = "pos",
  keep.extra.columns = TRUE
)

idx <- match(names(cpg_gr),
             rownames(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi))

# Sanity check
stopifnot(!anyNA(idx))

beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi <-
    beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi[idx, ]


dim(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi)
dim(pd_Li_617)


#CHECK
identical(
  rownames(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi),
  names(cpg_gr)
)
# MUST be TRUE

rownames(pd_Li_617) <- pd_Li_617$Basename
pd_Li_617 <- pd_Li_617[
    colnames(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi),
]

identical(
  colnames(beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi),
  rownames(pd_Li_617)
)
# MUST be TRUE

grSet <- GenomicRatioSet(
    gr = cpg_gr,
    Beta = beta_Noob_ChampFiltered_noXY_BMIQ_filteredLi,
    colData = pd_Li_617,
    annotation = c(
        array = "IlluminaHumanMethylationEPICv2",
        annotation = "20a1.hg38"
    )
)

R_grSet <- grSet[, colData(grSet)$LiResponse == "R"]
NR_grSet <- grSet[, colData(grSet)$LiResponse == "NR"]
MR_grSet <- grSet[, colData(grSet)$LiResponse == "MR"]


#=========================================
# DO EPIMUTAIONS 
#========================================
library(epimutacions)

#! first lets try using just a small amount of data to see if it works
set.seed(123)


# Function to safely sample
safe_sample <- function(grSet, n) {
    if (ncol(grSet) < n) {
        warning("Requested ", n, " samples but only ", ncol(grSet), 
                " available. Using all.")
        return(grSet)
    }
    sampled_idx <- sample(1:ncol(grSet), n, replace = FALSE)
    return(grSet[, sampled_idx])
}

R_random <- safe_sample(R_grSet, 20)
NR_random <- safe_sample(NR_grSet, 20)

# Show sample names
cat("Randomly selected R samples:\n")
print(colnames(R_random))


#ADD your function so that it will use your candRegsGR as candidate regions

get_candRegsGR <- function() {
    # just return your precomputed candidate regions
    return(candRegsGR)
}

#! FIND THE ACTUAL EPIMUTATIONS 

#using non responders as reference
epi_mvo <- epimutations(R_random, 
                        NR_random, 
                        method = "manova")


saveRDS(grSet, file = "grset.rds")

# ACTUAL EPIMUTATIONS
#! R vs NR using all samples
epi_mvo_true <- epimutations(
  R_grSet, 
  NR_grSet, 
  method = "manova")


write.table(epi_mvo_true,"epimutations_R_vs_NR_manova_allSamples.tsv",sep="\t",row.names=FALSE)


# epi_mahdist <- epimutations(
#   R_grSet, 
#   NR_grSet, 
#   method = "mahdist")

# write.table(epi_mahdist,"epimutations_R_vs_NR_mahdist_allSamples.tsv",sep="\t",row.names=FALSE)

epi_mlm <- epimutations(
  R_grSet, 
  NR_grSet, 
  method = "mlm")
write.table(epi_mlm,"epimutations_R_vs_NR_mlm_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_quantile <- epimutations(
  R_grSet, 
  NR_grSet, 
  method = "quantile")
write.table(epi_quantile,"epimutations_R_vs_NR_quantile_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_beta <- epimutations(
  R_grSet, 
  NR_grSet, 
  method = "beta")
write.table(epi_beta,"epimutations_R_vs_NR_beta_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_iForest <- epimutations(
  R_grSet, 
  NR_grSet, 
  method = "iForest")
write.table(epi_iForest,"epimutations_R_vs_NR_iForest_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

#! MR vs NR using all samples
epi_mvo_MR_NR <- epimutations(
  MR_grSet, 
  NR_grSet, 
  method = "manova")
write.table(epi_mvo_MR_NR,"epimutations_MR_vs_NR_manova_allSamples.tsv",sep="\t",row.names=FALSE)
gc()


# epi_mahdist_MR_NR <- epimutations(
#   MR_grSet, 
#   NR_grSet, 
#   method = "mahdist")

# write.table(epi_mahdist_MR_NR,"epimutations_MR_vs_NR_mahdist_allSamples.tsv",sep="\t",row.names=FALSE)

epi_mlm_MR_NR <- epimutations(
  MR_grSet, 
  NR_grSet, 
  method = "mlm")
write.table(epi_mlm_MR_NR,"epimutations_MR_vs_NR_mlm_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_quantile_MR_NR <- epimutations(
  MR_grSet, 
  NR_grSet, 
  method = "quantile")
write.table(epi_quantile_MR_NR,"epimutations_MR_vs_NR_quantile_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_beta_MR_NR <- epimutations(
  MR_grSet, 
  NR_grSet, 
  method = "beta")
write.table(epi_beta_MR_NR,"epimutations_MR_vs_NR_beta_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_iForest_MR_NR <- epimutations(
  MR_grSet, 
  NR_grSet, 
  method = "iForest")
write.table(epi_iForest_MR_NR,"epimutations_MR_vs_NR_iForest_allSamples.tsv",sep="\t",row.names=FALSE)
gc()


#TEST
setwd("/home/rachele/epivariants/results")
#! R vs MR using all samples
epi_mvo_R_MR <- epimutations(
  R_grSet, 
  MR_grSet, 
  method = "manova")
write.table(epi_mvo_R_MR,"epimutations_R_vs_MR_manova_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_mlm_R_MR <- epimutations(
  R_grSet, 
  MR_grSet, 
  method = "mlm")
write.table(epi_mlm_R_MR,"epimutations_R_vs_MR_mlm_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_quantile_R_MR <- epimutations(
  R_grSet, 
  MR_grSet, 
  method = "quantile")
write.table(epi_quantile_R_MR,"epimutations_R_vs_MR_quantile_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_beta_R_MR <- epimutations(
  R_grSet, 
  MR_grSet, 
  method = "beta")
write.table(epi_beta_R_MR,"epimutations_R_vs_MR_beta_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_iForest_R_MR <- epimutations(
  R_grSet, 
  MR_grSet, 
  method = "iForest")
write.table(epi_iForest_R_MR,"epimutations_R_vs_MR_iForest_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

#! MR vs R using all samples
epi_mvo_MR_R <- epimutations(
  MR_grSet, 
  R_grSet, 
  method = "manova")
write.table(epi_mvo_MR_R,"epimutations_MR_vs_R_manova_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_mlm_MR_R <- epimutations(
  MR_grSet, 
  R_grSet, 
  method = "mlm")
write.table(epi_mlm_MR_R,"epimutations_MR_vs_R_mlm_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_quantile_MR_R <- epimutations(
  MR_grSet, 
  R_grSet, 
  method = "quantile")
write.table(epi_quantile_MR_R,"epimutations_MR_vs_R_quantile_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_beta_MR_R <- epimutations(
  MR_grSet, 
  R_grSet, 
  method = "beta")
write.table(epi_beta_MR_R,"epimutations_MR_vs_R_beta_allSamples.tsv",sep="\t",row.names=FALSE)
gc()

epi_iForest_MR_R <- epimutations(
  MR_grSet, 
  R_grSet, 
  method = "iForest")
write.table(epi_iForest_MR_R,"epimutations_MR_vs_R_iForest_allSamples.tsv",sep="\t",row.names=FALSE)
gc()


#LOAD THEM 
library(readr)
library(data.table)
library(dplyr)

# ! load data 
# Specify the folder containing your TSV files
folder_path <- "/home/rachele/epivariants/results"

# List all the TSV files in the folder
tsv_files <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through each file and assign it to an object
for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  assign(file_name, read_tsv(file))
}

# Get the names of all dataframes in the environment
dataframes <- ls(pattern = ".*")

#!NOT SURE IF I SHOULD ANNOTATE OR NOT BEFORE WILL DECIDE AFTER WARDS 

# remove the once with 0 0 0 
dataframe_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
dataframe_names <- dataframe_names[dataframe_names != "tsv_files"]
dataframe_names <- dataframe_names[startsWith(dataframe_names, "epimutations")]

for (df_name in dataframe_names) {
  df <- get(df_name)

  df <- df[
    !is.na(df$chromosome) & !is.na(df$start) & !is.na(df$end) &
    df$chromosome != 0 & df$start != 0 & df$end != 0,
  ]

  assign(df_name, df, envir = .GlobalEnv)
}



#! FIND EPIVARIANTS SHARED BY MUTLIPLE SAMPLE IN SAME DF 

#OPTION!1: just look at the epimutations tables and see if there are any regions that appear in multiple samples (e.g. same chr, start, end)
library(dplyr)
library(tidyr)
library(purrr)

# List all your dataframes (excluding the tsv_files variable)
dataframe_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
dataframe_names <- dataframe_names[dataframe_names != "tsv_files"]

# Function to analyze epivariants within a dataframe
analyze_within_dataframe <- function(df, df_name) {
  cat("\n=== Analyzing:", df_name, "===\n")
  
  # Create a region identifier
  df$region_id <- paste(df$chromosome, df$start, df$end, sep = "_")
  
  # Check if region appears in multiple samples
  region_counts <- df %>%
    group_by(region_id, chromosome, start, end) %>%
    summarise(
      sample_count = n_distinct(sample),
      samples = paste(unique(sample), collapse = ", "),
      .groups = 'drop'
    )
  
  # Find regions in multiple samples
  multi_sample_regions <- region_counts %>%
    filter(sample_count > 1) %>%
    arrange(desc(sample_count))
  
  cat("Total unique regions:", nrow(region_counts), "\n")
  cat("Regions found in multiple samples:", nrow(multi_sample_regions), "\n")
  
  if(nrow(multi_sample_regions) > 0) {
    cat("\nTop regions in multiple samples:\n")
    print(head(multi_sample_regions, 10))
  }
  
  # Add dataframe name to results
  multi_sample_regions$source_df <- df_name
  
  return(list(
    all_regions = region_counts,
    multi_sample_regions = multi_sample_regions,
    df_name = df_name
  ))
}

# Run analysis on all dataframes
within_results <- list()
for(df_name in dataframe_names) {
  df <- get(df_name)
  within_results[[df_name]] <- analyze_within_dataframe(df, df_name)
}

# Combine all multi-sample regions
all_multi_sample <- map_dfr(within_results, ~ .x$multi_sample_regions)

cat("\n\n=== SUMMARY: Regions in Multiple Samples ===\n")
print(all_multi_sample %>%
  group_by(source_df) %>%
  summarise(
    n_multi_sample_regions = n(),
    max_samples_per_region = max(sample_count),
    .groups = 'drop'
  ))

#! OPTION 2
library(dplyr)
library(tidyr)

pd_Li_617 <- read.csv("~/epivariants/data/pd_Li_617.csv")

# First, let's see the structure of pd_Li_617 to understand the status mapping
cat("=== Checking pd_Li_617 structure ===\n")
print(head(pd_Li_617))
cat("\nUnique values in LiResponse:\n")
print(table(pd_Li_617$LiResponse))

# Create a mapping from Basename to LiResponse
sample_status <- pd_Li_617 %>%
  select(Basename, LiResponse) %>%
  distinct()

# List all epimutation dataframes
# dataframe_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
# dataframe_names <- dataframe_names[dataframe_names != "tsv_files"]
# dataframe_names <- dataframe_names[startsWith(dataframe_names, "epimutations")]
# dataframe_names <- dataframe_names[!dataframe_names %in% c("tsv_files", "pd_Li_617", "sample_status")]

cat("\n=== Found", length(dataframe_names), "epimutation dataframes ===\n")
print(dataframe_names)

# Function to analyze a SINGLE dataframe
analyze_single_dataframe <- function(df, df_name) {
  cat("\n" , rep("=", 50), "\n", sep="")
  cat("ANALYZING:", df_name, "\n")
  cat(rep("=", 50), "\n\n", sep="")
  
  # STEP 1: Add Status column by matching sample with pd_Li_617
  df_with_status <- df %>%
    left_join(sample_status, by = c("sample" = "Basename")) %>%
    dplyr::rename(Status = LiResponse)
  
  # Check how many samples matched
  matched_count <- sum(!is.na(df_with_status$Status))
  total_count <- nrow(df_with_status)
  cat("Added Status column: matched", matched_count, "of", total_count, "rows\n")
  
  # STEP 2: Create region identifier
  df_with_status$region_id <- paste(df_with_status$chromosome, 
                                    df_with_status$start, 
                                    df_with_status$end, sep = "_")
  
  # STEP 3: Group by region and analyze
  region_analysis <- df_with_status %>%
    group_by(region_id, chromosome, start, end) %>%
    summarise(
      # Count information
      sample_count = n_distinct(sample),
      
      # Sample information
      samples = paste(unique(sample), collapse = ", "),
      
      # Status information - comma separated
      Status = paste(unique(Status), collapse = ", "),
      
      # Sample + Status pairs
      sample_status_pairs = paste(paste(sample, "(", Status, ")", sep = ""), collapse = ", "),
      
      .groups = 'drop'
    ) %>%
    arrange(desc(sample_count))
  
  # STEP 4: Separate regions found in single vs multiple samples
  single_sample_regions <- region_analysis %>%
    filter(sample_count == 1)
  
  multi_sample_regions <- region_analysis %>%
    filter(sample_count > 1)
  
  # STEP 5: Print results for this dataframe
  cat("\n--- Summary for", df_name, "---\n")
  cat("Total unique regions:", nrow(region_analysis), "\n")
  cat("Regions in single sample:", nrow(single_sample_regions), "\n")
  cat("Regions in multiple samples:", nrow(multi_sample_regions), "\n")
  
  if(nrow(multi_sample_regions) > 0) {
    cat("\n=== REGIONS IN MULTIPLE SAMPLES ===\n")
    
    # Show detailed info for each multi-sample region
    for(i in 1:min(10, nrow(multi_sample_regions))) {
      region <- multi_sample_regions[i, ]
      cat("\nRegion", i, ":\n")
      cat("  Location:", region$chromosome, ":", region$start, "-", region$end, "\n")
      cat("  Found in", region$sample_count, "samples\n")
      cat("  Samples:", region$samples, "\n")
      cat("  Statuses:", region$Status, "\n")
      cat("  Sample(Status) pairs:", region$sample_status_pairs, "\n")
    }
    
    # If there are more than 10, mention it
    if(nrow(multi_sample_regions) > 10) {
      cat("\n... and", nrow(multi_sample_regions) - 10, "more regions\n")
    }
    
    # Count regions by number of samples
    cat("\n--- Breakdown by sample count ---\n")
    sample_count_summary <- multi_sample_regions %>%
      group_by(sample_count) %>%
      summarise(n_regions = n(), .groups = 'drop') %>%
      arrange(desc(sample_count))
    
    print(sample_count_summary)
    
    # Analyze status patterns in multi-sample regions
    cat("\n--- Status patterns in multi-sample regions ---\n")
    status_patterns <- multi_sample_regions %>%
      group_by(Status) %>%
      summarise(
        n_regions = n(),
        example_region = paste(head(region_id, 2), collapse = "; "),
        .groups = 'drop'
      ) %>%
      arrange(desc(n_regions))
    
    print(status_patterns)
  }
  
  # STEP 6: Return the results
  return(list(
    df_name = df_name,
    df_with_status = df_with_status,
    region_analysis = region_analysis,
    multi_sample_regions = multi_sample_regions,
    single_sample_regions = single_sample_regions
  ))
}

# STEP 7: Analyze each dataframe individually
all_results <- list()
for(df_name in dataframe_names) {
  df <- get(df_name)
  result <- analyze_single_dataframe(df, df_name)
  all_results[[df_name]] <- result
  
  # Optionally save the dataframe with status added
  assign(paste0(df_name, "_with_status"), result$df_with_status)
}




# STEP 8: Create a summary across all dataframes
cat("\n" , rep("=", 60), "\n", sep="")
cat("OVERALL SUMMARY ACROSS ALL DATAFRAMES\n")
cat(rep("=", 60), "\n\n", sep="")

# Create summary table
summary_table <- data.frame(
  Dataframe = character(),
  Total_Regions = numeric(),
  Single_Sample_Regions = numeric(),
  Multi_Sample_Regions = numeric(),
  Max_Samples_Per_Region = numeric(),
  stringsAsFactors = FALSE
)

for(df_name in dataframe_names) {
  result <- all_results[[df_name]]
  
  max_samples <- if(nrow(result$multi_sample_regions) > 0) {
    max(result$multi_sample_regions$sample_count)
  } else {
    0
  }
  
  summary_table <- rbind(summary_table, data.frame(
    Dataframe = df_name,
    Total_Regions = nrow(result$region_analysis),
    Single_Sample_Regions = nrow(result$single_sample_regions),
    Multi_Sample_Regions = nrow(result$multi_sample_regions),
    Max_Samples_Per_Region = max_samples
  ))
}

print(summary_table)

#!ITS USEFUL BUT I DONT NEED TO SAVE IT AT THE MOMENT 
# # STEP 9: Save results to files
# cat("\n=== Saving results to files ===\n")

# # Create output directory if it doesn't exist
# output_dir <- "within_df_analysis"
# if(!dir.exists(output_dir)) {
#   dir.create(output_dir)
# }

# # Save individual results for each dataframe
# for(df_name in dataframe_names) {
#   result <- all_results[[df_name]]
  
#   # Save multi-sample regions for this dataframe
#   if(nrow(result$multi_sample_regions) > 0) {
#     filename <- paste0(output_dir, "/", df_name, "_multi_sample_regions.tsv")
#     write.table(result$multi_sample_regions, 
#                 filename, 
#                 sep = "\t", 
#                 row.names = FALSE, 
#                 quote = FALSE)
#     cat("Saved:", filename, "\n")
#   }
  
#   # Save full region analysis for this dataframe
#   filename <- paste0(output_dir, "/", df_name, "_all_regions.tsv")
#   write.table(result$region_analysis, 
#               filename, 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
#   cat("Saved:", filename, "\n")
  
#   # Save dataframe with status added
#   filename <- paste0(output_dir, "/", df_name, "_with_status.tsv")
#   write.table(result$df_with_status, 
#               filename, 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
#   cat("Saved:", filename, "\n")
# }

# # Save overall summary
# filename <- paste0(output_dir, "/overall_summary.tsv")
# write.table(summary_table, 
#             filename, 
#             sep = "\t", 
#             row.names = FALSE, 
#             quote = FALSE)
# cat("Saved:", filename, "\n")

# # STEP 10: Find most interesting cases
# cat("\n=== MOST INTERESTING CASES ===\n")

# # Find regions with the most samples
# all_multi_sample_regions <- bind_rows(
#   lapply(all_results, function(x) {
#     if(nrow(x$multi_sample_regions) > 0) {
#       x$multi_sample_regions %>% mutate(df_name = x$df_name)
#     }
#   }),
#   .id = "origin"
# )

# if(nrow(all_multi_sample_regions) > 0) {
#   # Top 5 regions with most samples
#   top_regions <- all_multi_sample_regions %>%
#     arrange(desc(sample_count)) %>%
#     head(5)
  
#   cat("\nTop 5 regions with most samples:\n")
#   for(i in 1:nrow(top_regions)) {
#     region <- top_regions[i, ]
#     cat("\n", i, ". ", region$df_name, ":\n", sep="")
#     cat("   Region:", region$chromosome, ":", region$start, "-", region$end, "\n")
#     cat("   Samples:", region$sample_count, "\n")
#     cat("   Statuses:", region$Status, "\n")
#     cat("   Samples:", region$samples, "\n")
#   }
  
#   # Save top regions
#   filename <- paste0(output_dir, "/top_regions_across_dataframes.tsv")
#   write.table(top_regions, 
#               filename, 
#               sep = "\t", 
#               row.names = FALSE, 
#               quote = FALSE)
#   cat("\nSaved top regions to:", filename, "\n")
# }

# cat("\n" , rep("=", 60), "\n", sep="")
# cat("ANALYSIS COMPLETE\n")
# cat(rep("=", 60), "\n", sep="")
# cat("\nAll results saved in the '", output_dir, "' directory\n", sep="")





#FIND EPIVARIANTS SHARED BY MULTIPLE SAMPLE and or  IN DIFFERENT DF (e.g. R vs NR and MR vs NR)










































#=========================================
# ANNOTATE EPIMUTATIONS
#=========================================
library(epimutacions)
#TRIAL 1
epi_results_annotated <- annotate_epimutations(
    epi_results,
    db = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    build = "38",
    gene_col = "GencodeV41_Name",       # override gene column
    feat_col = "Regulatory_Feature_Group",
    relat_col = "Relation_to_Island"
)

# #TRIAL 2
# result <- annotate_cpg(epi_results,
#              db = IlluminaHumanMethylationEPICv2anno.20a1.hg38,
#              gene_col = "GencodeV41_Name",
#              feat_col = "Regulatory_Feature_Group",
#              relat_col = "Relation_to_Island",
#              build = "38")


# #TRAIL:3 in case the others didnt work 
# annotate_epimutations_hg38 <- function(
#     epi_results,
#     db = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
#     gene_col = "GencodeV41_Name",
#     feat_col = "Regulatory_Feature_Group",
#     relat_col = "Relation_to_Island",
#     build = "38",
#     ...
# ) {
#     # Annotate CpGs
#     epi_results <- annotate_cpg(
#         epi_results,
#         db = db,
#         gene_col = gene_col,
#         feat_col = feat_col,
#         relat_col = relat_col,
#         ...
#     )

#     # Add ENSEMBL regulatory annotation
#     epi_results <- add_ensemble_regulatory(epi_results, build = build)

#     return(epi_results)
# }

# epi_results_annotated <- annotate_epimutations_hg38(epi_results)

#=========================================
# PLOT EPIMUTATIONS
#=========================================
plot_epimutations(
    dmr = epi_results_annotated[1, ], 
    methy = grSet, 
    genome = "hg38",
    genes_annot = TRUE,
    regulation = TRUE
)


