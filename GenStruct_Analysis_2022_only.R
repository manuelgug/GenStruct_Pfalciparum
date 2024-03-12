
library(haven)
library(progress)
library(fs)
library(moire)
library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)
library(RColorBrewer)
library(Rtsne)
library(ggpubr)

# Objective 3: Describe genetic structure of P. falciparum population in Mozambique in 2021 and 2022 and investigate the origin of parasites with variants of concern

### INDEX

# a) Determine genetic diversity of Plasmodium falciparum in Mozambique
# • Intra-host diversity
#   ◦ Calculate MOI overall and means per province and region for each year (Naïve and eMOI, which incorporates intra-host relatedness between clones) ✔
#   ◦ Calculate % of polyclonal infections per province and region for each year ✔
#   ◦ Fixation index Fws?? ✔
# • Population genetic diversity: 
#   ◦ Calculate He per locus per province and region ✔
# Refs: Brokhattingen et al, preprint; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4 ; Kattenberg et al https://doi.org/10.1128/spectrum.00960-22 

# b) Determine population structure and connectivity between parasite populations in Mozambique
# • PCA to determine population structure ✔
# • Genetic differentiation (Fst) between provinces and regions ✔
# • Proportion of related pairwise infections using IBD between provinces and regions ✔
# Refs: Arnau; Brokhattingen et al, preprint; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4; Kattenberg et al https://doi.org/10.1128/spectrum.00960-22; Gerlovina et al, Genetics 2022

# c) Determine clonal emergence vs spread from other regions for parasites carrying VOC (dhps581 and dhps436)
# • Map parasites with VOC on PCA (color) ✔
# • Determine pairwise IBD between samples with variants of concern and wildtype parasites from the same or different areas.
# Refs: daSilva et al, https://doi.org/10.1038/s42003-023-04997-7; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4

# d) Kruskal–Wallis rank sum test will be used for the comparison of the distribution between populations, with Dunn’s test and Bonferroni correction for multiple testing in pairwise comparisons.


######################################################################
#-----------------------INPUT AND PREP DATA--------------------------#
######################################################################

#######################################################
# 1.- metadata (provided by Simone)
#######################################################

db <- read_dta('DrugRes_2021_2022_DB_ALLDATA_1Jan2024.dta')
colnames(db)
unique(db$run_id) #incomplete... must complete

# gonna need nida, year, region, province, run_id (ok)
# gonna need run_id of all samples (i can complete the ones that doesn't have it by grepping, np)
# which NIDA is the correct one? (according to a grep, it's NIDA2)
# MOI, calculate all together or per run? (ask simone/nanna!)
# entire runs or only samples in the DB?
# i don't have IMP_NEXTSEQ02_300522...
# no SMC2_run3?

# 1a.- complete run_id column with data from the actual runs

# Initialize an empty data frame to store the results
result_df <- data.frame(NIDA2 = character(), NEW_run_id = character(), stringsAsFactors = FALSE)

directory_path <- "../results_v0.1.8_RESMARKERS_FIX/"

# Create a progress bar
pb <- progress_bar$new(
  format = "[:bar] :percent ETA: :eta",
  total = length(dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$"))
)

# Iterate through the folders ending with _RESULTS_v0.1.8
for (folder_path in dir_ls(path = directory_path, regexp = "_RESULTS_v0.1.8$")) {
  pb$tick()  # Update progress bar
  
  folder_name <- path_file(folder_path)
  file_path <- file.path(folder_path, "amplicon_coverage.txt")
  
  cat("\n")
  print(folder_name)
  
  # Read the contents of /quality_report/amplicon_stats.txt
  if (file.exists(file_path)) {
    sample_coverage_content <- readLines(file_path)
    
    # Truncate each line after the first tab (\t) and return unique values and edit nida format to match the one from the db
    truncated_values <- unique(sapply(strsplit(sample_coverage_content, "\t"), function(x) x[1]))
    truncated_values <- gsub("_S.*$", "", truncated_values)
    truncated_values <- gsub("_", ".", truncated_values)
    truncated_values <- gsub("-", ".", truncated_values)
    truncated_values <- gsub("N", "", truncated_values)
    
    # Extract NIDA2 from runs using grep
    nida_values <- grep(paste(db$NIDA2, collapse = "|"), truncated_values, value = TRUE) #not all samples from runs are needed, only those in the db, hence this
    
    # Create a data frame with NIDA2 and folder_name
    if (length(nida_values) > 0) {
      temp_df <- data.frame(NIDA2 = nida_values, NEW_run_id = folder_name, stringsAsFactors = FALSE)
      
      # Append the results to the main data frame
      result_df <- rbind(result_df, temp_df)
    }
  }
}

result_df$NEW_run_id <- sub("_RESULTS_v0.1.8$", "", result_df$NEW_run_id)

#this should be TRUE if the same number of nidas is in the db and the results from the grep
length(result_df$NIDA2) == length(db$NIDA2) #it's not....
length(result_df$NIDA2) - length(db$NIDA2) # 113 repeated nidas...

#check for repeated nidas 
repeated_nidas <- names(table(result_df$NIDA2)[table(result_df$NIDA2) > 1]) #there are duplicate nidas OOF
repeated_nidas_df<- result_df[result_df$NIDA2 %in% repeated_nidas,]
repeated_nidas_df <- repeated_nidas_df[order(repeated_nidas_df$NIDA2), ]
length(repeated_nidas_df$NIDA2)
length(unique(repeated_nidas_df$NIDA2))

#ask team about these nidas
write.csv(repeated_nidas_df, "repeated_nidas.csv", row.names = F)

#######################################################
# 2.- genomic data from all runs (allele data, resmarkers, haplos?)
#######################################################
runs <- unique(paste0(result_df$NEW_run_id, "_RESULTS_v0.1.8_FILTERED"))

folder_path <- paste0("../results_v0.1.8_RESMARKERS_FIX/", runs)
allele_data_files <- file.path(folder_path, "allele_data_global_max_0_filtered.csv")

# Import filtered allele data tables to the list
allele_data_list <- list()
for (file in allele_data_files) {
  allele_data <- read.csv(file)
  allele_data_list <- append(allele_data_list, list(allele_data))
}

#format the imported dfs
for (i in seq_along(allele_data_list)) {
  df <- allele_data_list[[i]]
  
  df$sampleID <- gsub("_S.*$", "", df$sampleID)
  df$sampleID <- gsub("_", ".", df$sampleID)
  df$sampleID <- gsub("-", ".", df$sampleID)
  df$sampleID <- gsub("N", "", df$sampleID)
  df$sampleID <- gsub("\\.0", "", df$sampleID)
  
  colnames(df)[1] <- "sample_id"
  
    df <- df %>% ### TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
      group_by(sample_id, locus, pseudo_cigar) %>%
      summarize(reads = sum(reads),
                norm.reads.locus = sum(norm.reads.locus),
                Category = first(Category)) %>%
      mutate(allele = paste(locus, ".", row_number(), sep = ""))
    
    
  # Update the modified data frame back in the list
  df<- as.data.frame(df)
  allele_data_list[[i]] <- df
}

names(allele_data_list) <- runs

#get rid of replicate nidas keep the sample with the most reads across runs for each one of the replicates
sum_reads <- function(df) {
  aggregate(reads ~ sample_id, data = df, sum)
}

summed_reads <- lapply(allele_data_list, sum_reads)

# Change colnames of reads for the name of the respective df
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})

# Merge all data frames by sample_id
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)

#keep rows with replicates
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

# Exclude the first column and find the column with the maximum value for each row
merged_df_dups$BEST_RUN <- colnames(merged_df_dups)[apply(merged_df_dups[-1], 1, which.max) + 1]

# #did i found all replicates?
# if (all(repeated_nidas_df$NIDA2 %in% merged_df_dups$sample_id)) {
#   print("You found all replicates. Proceed with removal")
# }else{
#   "grab a coffee"
# }

#remove replicates, keep the best
for (i in 1:nrow(merged_df_dups)) {
  best_run <- merged_df_dups$BEST_RUN[i]
  sample_id <- merged_df_dups$sample_id[i]
  
  # Loop through allele_data_list
  for (j in seq_along(allele_data_list)) {
    df <- allele_data_list[[j]]
    
    # Exclude the df named after BEST_RUN
    if (names(allele_data_list)[j] != best_run) {
      allele_data_list[[j]] <- df[!(df$sample_id %in% sample_id), ]
    }
  }
}

#final check:
summed_reads <- lapply(allele_data_list, sum_reads)
summed_reads <- lapply(names(summed_reads), function(df_name) {
  df <- summed_reads[[df_name]]
  names(df)[2] <- df_name
  return(df)
})
merged_df_dups <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), summed_reads)
merged_df_dups <- merged_df_dups[rowSums(is.na(merged_df_dups)) < length(colnames(merged_df_dups)) - 2, ] #keep rows with more than 14 NAs (that is, that have reads in more than one run)

if (dim(merged_df_dups)[1] == 0){
  print("NO MORE REPLICATES.")
}else{
  "grab a coffee"
}

# Define a function to filter rows based on criteria: remove controls, basically
filter_rows <- function(df) {
  # Filter rows based on criteria
  filtered_df <- df[!grepl("^[A-Za-z]|3D", df$sample_id), ]
  return(filtered_df)
}

# Apply the filter_rows function to each dataframe in allele_data_list
allele_data_list <- lapply(allele_data_list, filter_rows)

#visual check:
for (df in allele_data_list) {
  cat("sample size:", as.character(length(unique(df$sample_id))))
  cat("\n")
  print(unique(df$sample_id))
  
}

#save allele_data_list
saveRDS(allele_data_list, "allele_data_list.RDS")

#######################################################
# 3.- GENOMIC + DB MERGING, FILTERING ETC. (DATA PREP)
#######################################################

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
if (!("n.alleles" %in% colnames(combined_df))){
  combined_df <- combined_df %>%
    group_by(sample_id, locus) %>%
    mutate(n.alleles = n_distinct(pseudo_cigar))
}

# merge with metadata
colnames(combined_df)[1]<- c("NIDA2")
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region", "seasonality", "run_id")], by="NIDA2", all.y =T) #forcing adding all db nidas

#check for columns with NAs. THESE SAMPLES WERE BAD QUALITY AND THUS FILTERED OUT DURING CONTAMINANTS FILTERING
removed_samples <- combined_df_merged[is.na(combined_df_merged$Category),]
removed_samples$NIDA2

# Remove rows with NIDA2 matching values in removed_samples
combined_df_merged <- combined_df_merged[!combined_df_merged$NIDA2 %in% removed_samples$NIDA2, ]


#sanity check
if (sum(is.na(combined_df_merged[, !colnames(combined_df_merged) %in% "seasonality"])) == 0) {
  print("No NAs ✔")
} else {
  print("grab a coffee.")
}

if( sum(!(combined_df_merged$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also the metadata db. No weird samples ✔")
}else{
  print("grab another coffee.")
}


#KEEP ONLY DIVERTSITY LOCI!
# are there samples with no Diversity loci sequenced?
unique_categories <- combined_df_merged %>%
  group_by(NIDA2) %>%
  summarize(unique_categories = toString(unique(Category)))

more_removed_sample <- unique_categories[!grepl("Diversity", unique_categories$unique_categories), ] # 1 sample didn't have diversity loci, probably due to a missing pool or something... no problem.
more_removed_sample

combined_df_merged <- combined_df_merged[combined_df_merged$Category == "Diversity",] 

#KEEP ONLY 2022 samples!
combined_df_merged <- combined_df_merged[combined_df_merged$year == 2022,] 


## FURTHER FILTERING

# 1) MAF filering (< 0.01)
combined_df_merged <- combined_df_merged[combined_df_merged$norm.reads.locus  > 0.01, ]

# 2) check for coverage (>100 loci with >threshold reads)
# Define thresholds for read depth
thresholds <- c(25, 50, 100, 200)
count_list <- list()

# Loop over each threshold
for (threshold in thresholds) {
  # Calculate unique loci counts for the current threshold
  count <- combined_df_merged %>%
    group_by(NIDA2, locus) %>%
    summarize(total_reads = sum(reads)) %>%
    group_by(NIDA2) %>%
    filter(total_reads > threshold) %>%
    summarize(!!paste("unique_loci_", threshold, sep = "") := n_distinct(locus))

  count_list[[paste("count_", threshold, sep = "")]] <- count
}

result_df <- Reduce(function(x, y) left_join(x, y, by = "NIDA2"), count_list)

#count cells above 100 for each column with the "unique" substring: here, i'm calculating the sample size for each reads count threshold using 100 loci as cutoff: samples with <100 read count per loci below the threshold should be removed
count_above_100 <- function(x) sum(x >= 100, na.rm = TRUE)
unique_columns <- grep("unique", colnames(result_df), value = TRUE)

count_results <- result_df %>%
  summarise_at(vars(unique_columns), count_above_100)

count_results

# decided going with a threshold of >= 100 read depth,
samples_to_keep <- result_df[result_df$unique_loci_100 >= 100, ]$NIDA2

combined_df_merged <- combined_df_merged[combined_df_merged$NIDA2 %in% samples_to_keep, ]


# 3) remove bad diversity loci: those that are not in at least 100 samples with a read depth of 100. 
# Group by NIDA2 and locus, then summarize the total reads
locus_read_depth <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  summarize(read_depth = sum(reads))

# Count the number of samples with read depth greater than 100 for each locus
locus_counts <- locus_read_depth %>%
  group_by(locus) %>%
  summarize(samples_above_100 = sum(read_depth > 100))

loci_to_keep <- locus_counts$locus

combined_df_merged <- combined_df_merged[combined_df_merged$locus %in% loci_to_keep, ]

library(forcats)

# MANAGE SEASONALITY
seasonality_factor <- as_factor(combined_df_merged$seasonality)
seasonality_recode <- recode(seasonality_factor, "1" = "Rainy", "2" = "Dry")
seasonality_char <- as.character(seasonality_recode)
combined_df_merged$seasonality <- seasonality_char

combined_df_merged$province <- ifelse(
  is.na(combined_df_merged$seasonality),
  combined_df_merged$province,
  paste0(combined_df_merged$province, "_", combined_df_merged$seasonality)
)


# # JUST IN CASE... recalculate n.alleles for each locus of each sample
# combined_df_merged <- combined_df_merged %>%
#     group_by(NIDA2, locus) %>%
#     mutate(n.alleles = n_distinct(pseudo_cigar))

#recount n.alleles after filtering
combined_df_merged <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(n.alleles = n_distinct(pseudo_cigar))

combined_df_merged <- as.data.frame(combined_df_merged)


# FILTERING RESULTS
SS <- length(unique(combined_df_merged$NIDA2))
cat("Final sample size is", as.character(SS))

LC <- length(unique(combined_df_merged$locus))
cat("Final Diversity loci count is", as.character(LC))

#save allele_data_list
saveRDS(combined_df_merged, "combined_df_merged_2022_only.RDS")

#######################################################
# 4.- CHECK SAMPLE SIZES FOR EACH PAIR OF VARIABLES: sample size affects He calculation, probably will need rarefactions or something similar
#######################################################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


sample_size_provinces <- combined_df_merged %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

combined_df_merged_nodry <- combined_df_merged %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions


######################################################################
#----------------------------ANALYZE DATA----------------------------#
######################################################################

#######################################################
# 5.- SUFFICIENCY OF SAMPLINGS
#######################################################

#ACCUMULATION CURVES (read this https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885658/; ver genotype_curve de paquete poppr?)

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
                                   year = combined_df_merged$year, 
                                   province = combined_df_merged$province,
                                   region = combined_df_merged$region,
                                   locus = combined_df_merged$locus,
                                   n.alleles = combined_df_merged$n.alleles,
                                   norm.reads.locus = combined_df_merged$norm.reads.locus,
                                   allele = paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar),
                                   run_id = combined_df_merged$run_id))

# PROVINCE

# Initialize a list to store the rarefaction curves for each year
accum_curves_2022 <- list()

# Get unique years and provinces
unique_provinces <- unique(raref_input$province)

# Iterate over each province
for (province in unique_provinces) {
  
  print(province)
  
  # Subsetting the data for 2022
  sub_2022 <- raref_input[raref_input$year == 2022 & raref_input$province == province, ]
  
  # Check if there are unique NIDA2s for 2022
  if (length(unique(sub_2022$NIDA2)) > 0) {
    # Initialize a list to store unique alleles for each NIDA2 for 2022
    unique_alleles_2022 <- list()
    
    # Iterate over each unique NIDA2 for 2022
    for (nida in unique(sub_2022$NIDA2)) {
      subset_data <- sub_2022[sub_2022$NIDA2 == nida, ]
      unique_alleles_it <- unique(subset_data$allele)
      unique_alleles_2022[[as.character(nida)]] <- unique_alleles_it
    }
    
    # Get unique alleles across all elements for 2022
    all_unique_alleles_2022 <- unique(unlist(unique_alleles_2022))
    
    # Create a matrix to store presence/absence of unique alleles for each element for 2022
    presence_matrix_2022 <- sapply(unique_alleles_2022, function(x) {
      as.integer(all_unique_alleles_2022 %in% x)
    })
    
    # Convert the matrix to a dataframe for 2022
    presence_df_2022 <- as.data.frame(presence_matrix_2022)
    presence_df_2022 <- t(presence_df_2022)
    rownames(presence_df_2022) <- names(unique_alleles_2022)
    colnames(presence_df_2022) <- all_unique_alleles_2022
    
    # CALCULATE CURVE for 2022
    accum_curve_2022 <- specaccum(presence_df_2022, 'random', permutations = 100)
    accum_curves_2022[[province]] <- accum_curve_2022
  }
}

# Select 9 colors from the Paired palette
colors <- brewer.pal(11, "Paired")

pdf("accumulation_curves_provinces.pdf", width = 12, height = 8)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", ylab ="Alleles", main = "Allele Accumulation Curves per Province", xlim = c(0,130), ylim = c(0,1800))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 110, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.7)

dev.off()
#conclusion....

raref_input_region <- raref_input[!(raref_input$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove DRY season pops from region analysis

# REGION

# Initialize a list to store the rarefaction curves for each year
accum_curves_2022 <- list()

# Get unique years and provinces
unique_provinces <- unique(raref_input_region$region)

# Iterate over each region
for (region in unique_provinces) {
  
  print(region)
  
  # Subsetting the data for 2022
  sub_2022 <- raref_input_region[raref_input_region$year == 2022 & raref_input_region$region == region, ]
  
  # Check if there are unique NIDA2s for 2022
  if (length(unique(sub_2022$NIDA2)) > 0) {
    # Initialize a list to store unique alleles for each NIDA2 for 2022
    unique_alleles_2022 <- list()
    
    # Iterate over each unique NIDA2 for 2022
    for (nida in unique(sub_2022$NIDA2)) {
      subset_data <- sub_2022[sub_2022$NIDA2 == nida, ]
      unique_alleles_it <- unique(subset_data$allele)
      unique_alleles_2022[[as.character(nida)]] <- unique_alleles_it
    }
    
    # Get unique alleles across all elements for 2022
    all_unique_alleles_2022 <- unique(unlist(unique_alleles_2022))
    
    # Create a matrix to store presence/absence of unique alleles for each element for 2022
    presence_matrix_2022 <- sapply(unique_alleles_2022, function(x) {
      as.integer(all_unique_alleles_2022 %in% x)
    })
    
    # Convert the matrix to a dataframe for 2022
    presence_df_2022 <- as.data.frame(presence_matrix_2022)
    presence_df_2022 <- t(presence_df_2022)
    rownames(presence_df_2022) <- names(unique_alleles_2022)
    colnames(presence_df_2022) <- all_unique_alleles_2022
    
    # CALCULATE CURVE for 2022
    accum_curve_2022 <- specaccum(presence_df_2022, 'random', permutations = 100)
    accum_curves_2022[[region]] <- accum_curve_2022
  }
}

# Select 3 colors from the Paired palette
colors <- brewer.pal(3, "Paired")

pdf("accumulation_curves_regions.pdf", width = 12, height = 8)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", ylab ="Alleles", main = "Allele Accumulation Curves per Region", xlim = c(0,350), ylim = c(0,2700))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 310, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.7)

dev.off()


#######################################################
# 6.- calculate MOI and eMOI for each run
#######################################################

combined_df_merged <- readRDS("combined_df_merged_2022_only.RDS") 

combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

colnames(combined_df_merged)[1]<- c("sample_id")

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(combined_df_merged)
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 20)
pt_num_threads <- 20

# run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10); saveRDS(mcmc_results, "all_samples_complete_filtered_MOIRE-RESULTS_2022_only_FOR_MOI.RDS")


#######################################################
# 7.- Present MOI/eMOI results overall and means per province and region for each year
#######################################################

mcmc_results <- readRDS("FINAL_MOIRE_RESULTS/all_samples_complete_filtered_MOIRE-RESULTS_2022_only_FOR_MOI.RDS")
combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


eff_coi <- moire::summarize_effective_coi(mcmc_results)
naive_coi <- moire::summarize_coi(mcmc_results)

# Merge the summaries by sample_id
coi_results <- merge(eff_coi, naive_coi, by = "sample_id")

# label mono and poly infections. NOTE: "proportion of polyclonal infections (eMOI>1.1)" from Nanna's manuscript
coi_results$polyclonal_from_ecoi_med <- ifelse(coi_results$post_effective_coi_med > 1.1, "polyclonal", "monoclonal")

#merge with categorical variables
colnames(coi_results)[1] <- "NIDA2"
coi_results <- merge(coi_results, db[c("NIDA2", "province", "region", "seasonality")], by="NIDA2")


# MANAGE SEASONALITY
seasonality_factor <- as_factor(coi_results$seasonality)
seasonality_recode <- recode(seasonality_factor, "1" = "Rainy", "2" = "Dry")
seasonality_char <- as.character(seasonality_recode)
coi_results$seasonality <- seasonality_char

coi_results$province <- ifelse(
  is.na(coi_results$seasonality),
  coi_results$province,
  paste0(coi_results$province, "_", coi_results$seasonality)
)

#coi_results <- coi_results[coi_results$year == 2022,]  #TEMPORAL FIX UNTIL MOIRE DATA FOR 2022 ONLY COMES OOUT


coi_results_region <- coi_results[!(coi_results$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove DRY season pops from region analysis

# % polyclonal infections on each province and region
polyclonal_percentage_region <- coi_results_region %>%
  group_by(region) %>%
  summarise(polyclonal_percentage_region = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

polyclonal_percentage_province <- coi_results %>%
  group_by(province) %>%
  summarise(polyclonal_percentage_province = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()


provinces <- c("Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")

coi_results$province <- factor(coi_results$province, levels = provinces)
coi_results$region <- factor(coi_results$region, levels = regions)
coi_results_region$region <- factor(coi_results_region$region, levels = regions)
polyclonal_percentage_region$region <- factor(polyclonal_percentage_region$region, levels = regions)
polyclonal_percentage_province$province <- factor(polyclonal_percentage_province$province, levels = provinces)


a <- ggplot(coi_results, aes(x = naive_coi, fill = region)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_wrap(~ province , scales = "fixed", nrow = 1) + 
  labs(title = "",
       x = "Naive COI",
       y = "Frequency",
       fill = "Province") +
  theme_minimal()

a

ggsave("naive_coi_provinces_ditros.png", a, width = 14, height = 6, bg = "white")


#naive coi pairwise kruskal wallis provinces
pairwise_province_naive_coi <- pairwise.wilcox.test(coi_results$naive_coi, 
                                                    coi_results$province, p.adjust.method = "bonferroni")

pairwise_province_naive_coi <- melt(pairwise_province_naive_coi[[3]])
signif_p.pairwise_province_naive_coi<- pairwise_province_naive_coi[pairwise_province_naive_coi$value <0.05 & !is.na(pairwise_province_naive_coi$value),]


pairwise_province_combinations <- lapply(1:nrow(signif_p.pairwise_province_naive_coi), function(i) {
  as.character(c(signif_p.pairwise_province_naive_coi[i, "Var1"], 
    signif_p.pairwise_province_naive_coi[i, "Var2"]))
})

a1 <- ggplot(coi_results, aes(x = province, y = naive_coi, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "Naive COI") +
  guides(color = FALSE)+
  stat_compare_means(comparisons = pairwise_province_combinations, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

a1

ggsave("naive_coi_provinces_violin.png", a1, width = 8, height = 6, bg = "white")

b <- ggplot(coi_results_region, aes(x = naive_coi, fill = region)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_grid( ~ region, scales = "free") +
  labs(title = "",
       x = "Naive COI",
       y = "Frequency",
       fill = "Region") +
  theme_minimal() +
  guides(fill = FALSE)

b

ggsave("naive_coi_regions_ditros.png", b, width = 14, height = 6, bg = "white")


#naive coi pairwise kruskal wallis regions
pairwise_region_naive_coi <- pairwise.wilcox.test(coi_results_region$naive_coi, 
                                                  coi_results_region$region, p.adjust.method = "bonferroni")

pairwise_region_naive_coi <- melt(pairwise_region_naive_coi[[3]])
signif_p.pairwise_region_naive_coi<- pairwise_region_naive_coi[pairwise_region_naive_coi$value <0.05 & !is.na(pairwise_region_naive_coi$value),]


pairwise_region_combinations <- lapply(1:nrow(signif_p.pairwise_region_naive_coi), function(i) {
  as.character(c(signif_p.pairwise_region_naive_coi[i, "Var1"], 
                 signif_p.pairwise_region_naive_coi[i, "Var2"]))
})


b1 <- ggplot(coi_results_region, aes(x = region, y = naive_coi, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "Naive COI") +
  guides(color = FALSE) +
  stat_compare_means(comparisons = pairwise_region_combinations, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

b1

ggsave("naive_coi_regions_violin.png", b1, width = 8, height = 6, bg = "white")


c <- ggplot(coi_results, aes(x = post_effective_coi_med, fill = region)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_wrap(~ province , scales = "fixed", nrow = 1) + 
  labs(title = "",
       x = "Post Effective COI Median",
       y = "Frequency",
       fill = "Province") +
  theme_minimal() +
  guides(fill = FALSE) 

c

ggsave("ecoi_provinces_ditros.png", c, width = 14, height = 6, bg = "white")


#ecoi pairwise kruskal wallis provinces
pairwise_province_ecoi <- pairwise.wilcox.test(coi_results$post_effective_coi_med, 
                                               coi_results$province, p.adjust.method = "bonferroni")

pairwise_province_ecoi <- melt(pairwise_province_ecoi[[3]])
signif_p.pairwise_province_ecoi<- pairwise_province_ecoi[pairwise_province_ecoi$value <0.05 & !is.na(pairwise_province_ecoi$value),]

pairwise_province_combinations <- lapply(1:nrow(signif_p.pairwise_province_ecoi), function(i) {
  as.character(c(signif_p.pairwise_province_ecoi[i, "Var1"], 
                 signif_p.pairwise_province_ecoi[i, "Var2"]))
})

c1 <- ggplot(coi_results, aes(x = province, y = post_effective_coi_med, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "eCOI") +
  guides(color = FALSE) +
  stat_compare_means(comparisons = pairwise_province_combinations, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

c1

ggsave("ecoi_provinces_violin.png", c1, width = 8, height = 6, bg = "white")


d <- ggplot(coi_results_region, aes(x = post_effective_coi_med, fill = region)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_grid( ~ region, scales = "free") +
  labs(title = "",
       x = "Post Effective COI Median",
       y = "Frequency",
       fill = "Region") +
  theme_minimal() +
  guides(fill = FALSE) 

d

ggsave("ecoi_regions_ditros.png", d, width = 10, height = 6, bg = "white")


#ecoi pairwise kruskal wallis regions
pairwise_region_ecoi <- pairwise.wilcox.test(coi_results_region$post_effective_coi_med, 
                                             coi_results_region$region, p.adjust.method = "bonferroni")

pairwise_region_ecoi <- melt(pairwise_region_ecoi[[3]])
signif_p.pairwise_region_ecoi<- pairwise_region_ecoi[pairwise_region_ecoi$value <0.05 & !is.na(pairwise_region_ecoi$value),]

pairwise_region_combinations <- lapply(1:nrow(signif_p.pairwise_region_ecoi), function(i) {
  as.character(c(signif_p.pairwise_region_ecoi[i, "Var1"], 
                 signif_p.pairwise_region_ecoi[i, "Var2"]))
})


d1 <- ggplot(coi_results_region, aes(x = region, y = post_effective_coi_med, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +
  labs(x = "", y = "eCOI") +
  guides(color = FALSE) +
  stat_compare_means(comparisons = pairwise_region_combinations, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

d1

ggsave("ecoi_region_violin.png", d1, width = 8, height = 6, bg = "white")




##polyclonal percentage
e <- ggplot(polyclonal_percentage_region, aes(x = region, y = polyclonal_percentage_region, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "% Polyclonal Infections") +
  #facet_wrap(~region, scales = "fixed", ncol = 3) +
  #scale_fill_manual(values = c("2022" = "orange")) + 
  theme_minimal()+
  guides(fill = FALSE) 
e

ggsave("perc_polyclonal_regions.png", e, width = 6, height = 4, bg = "white")


polyclonal_percentage_province <- merge(polyclonal_percentage_province, unique(coi_results[c("province", "region")]), by="province")

f <- ggplot(polyclonal_percentage_province, aes(x = province, y = polyclonal_percentage_province,  fill = region))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "%Polyclonal Infections") +
  #facet_wrap(~province, scales = "fixed", nrow = 1) +
  #scale_fill_manual(values = c("2022" = "orange")) + 
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
f

ggsave("perc_polyclonal_provinces.png", f, width = 8, height = 6, bg = "white")


###
# coi results for Simone
coi_for_db <- coi_results[c("NIDA2", "naive_coi", "post_effective_coi_med")]

write.csv(coi_for_db, "coi_for_db.csv", row.names = F)

plot(coi_for_db$naive_coi, coi_for_db$post_effective_coi_med)
###

#######################################################
# 8.- calculate He for each population (per region/province)
#######################################################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

# RUN IN CLUSTER:
# Define function to run MOIRE and save results
run_moire <- function(df, output_name) {
  colnames(df)[1] <- "sample_id"
  
  # set MOIRE parameters
  dat_filter <- moire::load_long_form_data(df)
  burnin <- 1e4
  num_samples <- 1e4
  pt_chains <- seq(1, .5, length.out = 20)
  
  # run moire
  mcmc_results <- moire::run_mcmc(
    dat_filter, is_missing = dat_filter$is_missing,
    verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
    pt_chains = pt_chains, pt_num_threads = length(pt_chains),
    thin = 10
  )
  
  # checkpoint
  saveRDS(mcmc_results, paste0(output_name, "_MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS"))
}

# Create a list of data frames and corresponding years
data_frames <- list(combined_df_merged)
years <- list("2022")

# Loop over each province
for (i in seq_along(data_frames)) {
  year <- years[[i]]
  df <- data_frames[[i]]
  
  for (province in unique(df$province)) {
    province_df <- df[df$province == province, ]
    
    # Run MOIRE
    it_pr <-  paste0(province, "_", year)
    print(it_pr)
    
    run_moire(province_df, it_pr)
  }
}

data_frames[[1]] <- data_frames[[1]][!(data_frames[[1]]$province %in% c("Maputo_Dry", "Manica_Dry")), ] # remove DRY season pops from region analysis

# Loop over each region
for (i in seq_along(data_frames)) {
  year <- years[[i]]
  df <- data_frames[[i]]
  
  for (region in unique(df$region)) {
    province_df <- df[df$region == region, ]
    
    # Run MOIRE
    it_re <-  paste0(region, "_", year)
    print(it_re)
    
    run_moire(province_df, it_re)
  }
}

#######################################################
# 9.- He and Fws results 
#######################################################

### no need to remove DRY season pops from region analysis because it already was removed qhen running moire by population 


combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = "FINAL_MOIRE_RESULTS", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)
rds_files <- rds_files[!rds_files %in% "FINAL_MOIRE_RESULTS/all_samples_complete_filtered_MOIRE-RESULTS_2022_only_FOR_MOI.RDS"] 


# Load each RDS file into the list with the file name as the list name
He_results_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  He_results_list[[file_name]] <- readRDS(file)
}

# Loop through each element in He_results_list
processed_He_results <- data.frame()

for (i in seq_along(He_results_list)) {
  
  # Summarize He
  He_results <- moire::summarize_he(He_results_list[[i]])
  He_results$population <- names(He_results_list[i])
  
  # Add the processed He results to the list
  processed_He_results <- rbind(processed_He_results, He_results)
}

#formatting categories
processed_He_results$population <- gsub("_2022_MOIRE-RESULTS_FOR_ALLELE_FREQS", "", processed_He_results$population)

######
# keep amplicons with high He:

he_amps <- processed_He_results %>%
  group_by(locus) %>%
  summarize(mean = mean(post_stat_mean)) %>%
  arrange(desc(mean))

#keep top n% amplicons with highest He (or use a numeric threshold like > 0.1 of He or whatever... test.)
perc_25<- round(length(unique(he_amps$locus))*0.95)
he_amps_top50 <- he_amps[1:perc_25,]

# FILTER
processed_He_results <- processed_He_results[processed_He_results$locus %in% he_amps_top50$locus,]
#####


library(stringr)
processed_He_results <- processed_He_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))

#processed_He_results$population <- gsub("TEST_|_202.*", "", processed_He_results$population)
#processed_He_results$year <- as.numeric(processed_He_results$year)


# 2) calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
heterozygosity_data <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))


# 3) calculate 1-Fws: 1-Fws = Hw/He
#merge He from provinces
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus")) %>%
  filter(population == province) %>%
  mutate(He_province = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, He_province)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus"))

heterozygosity_data <- distinct(heterozygosity_data)
heterozygosity_data

#merge He from regions
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus")) %>%
  filter(population == region) %>%
  mutate(He_region = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, He_region)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus"))

heterozygosity_data <- distinct(heterozygosity_data)

#calculate 1-Fws for province and region as populations, for each year
heterozygosity_data$fws_province <- heterozygosity_data$Hw/heterozygosity_data$He_province
heterozygosity_data$fws_region <- heterozygosity_data$Hw/heterozygosity_data$He_region

 # Columns to keep
columns_to_keep <- c("NIDA2", "locus", "province", "region", "n.alleles",
                     "Hw", "He_province", "He_region", "fws_province", "fws_region")
# Filter columns
heterozygosity_data_filtered <- heterozygosity_data %>%
  select(all_of(columns_to_keep))

# Keep unique rows
heterozygosity_data_filtered <- distinct(heterozygosity_data_filtered)

# convert NA to 0 (monoallelic loci that doesn't have heterozygosity)
heterozygosity_data_filtered[is.na(heterozygosity_data_filtered)] <- 0

#sanity check
if ((length(unique(heterozygosity_data_filtered$NIDA2)) == length(unique(combined_df_merged$NIDA2))) & 
    (length(unique(heterozygosity_data_filtered$locus)) == length(unique(combined_df_merged$locus))) & 
    (length(unique(heterozygosity_data_filtered$province)) == length(unique(combined_df_merged$province))) & 
    (length(unique(heterozygosity_data_filtered$region)) == length(unique(combined_df_merged$region)))) {
  print("All looks good.")
}else{
  print("grab a coffee.")
}

#WHY IS IT NOT BETWEEN 0 AND 1??? (should it be?)
mean_Fws_per_individual<- heterozygosity_data_filtered %>%
  group_by(NIDA2, province, region) %>%
  summarize(mean_indiv_fws_province = mean(fws_province),
            mean_indiv_fws_region = mean(fws_region))


provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")

mean_Fws_per_individual$province <- factor(mean_Fws_per_individual$province, levels = provinces)
mean_Fws_per_individual$region <- factor(mean_Fws_per_individual$region, levels = regions)


#naive coi pairwise kruskal wallis provinces
pairwise_province_fws <- pairwise.wilcox.test(mean_Fws_per_individual$mean_indiv_fws_province, 
                                              mean_Fws_per_individual$province, p.adjust.method = "bonferroni")

pairwise_province_fws <- melt(pairwise_province_fws[[3]])
signif_p.pairwise_province_fws<- pairwise_province_fws[pairwise_province_fws$value <0.05 & !is.na(pairwise_province_fws$value),]

combos_pairwise_province_fws_province <- lapply(1:nrow(signif_p.pairwise_province_fws), function(i) {
  as.character(c(signif_p.pairwise_province_fws[i, "Var1"], 
                 signif_p.pairwise_province_fws[i, "Var2"]))
})

prov_fws <- ggplot(mean_Fws_per_individual, aes(x = province, y = mean_indiv_fws_province, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_discrete(name = "Region") +  # Customize legend title
  #ggtitle("Province Connectivity") +
  labs(x = "Province", y = "Genome-wide 1-Fws") +
  guides(color = FALSE) +
  stat_compare_means(comparisons = combos_pairwise_province_fws_province, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

prov_fws

ggsave("province_fws.png", prov_fws, width = 8, height = 6, bg = "white")

mean_Fws_per_individual_nodry <- mean_Fws_per_individual %>%
  filter(!grepl("Dry", province))

pairwise_region_fws <- pairwise.wilcox.test(mean_Fws_per_individual_nodry$mean_indiv_fws_region, 
                                            mean_Fws_per_individual_nodry$region, p.adjust.method = "bonferroni")

pairwise_region_fws <- melt(pairwise_region_fws[[3]])
signif_p.pairwise_region_fws<- pairwise_region_fws[pairwise_region_fws$value <0.05 & !is.na(pairwise_region_fws$value),]

combos_pairwise_region_fws_province <- lapply(1:nrow(signif_p.pairwise_region_fws), function(i) {
  as.character(c(signif_p.pairwise_region_fws[i, "Var1"], 
                 signif_p.pairwise_region_fws[i, "Var2"]))
})

reg_fws <- ggplot(mean_Fws_per_individual_nodry, aes(x = region, y = mean_indiv_fws_region, fill = region)) +
  geom_violin(width = 1, aes(color = region), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "", y = "Genome-wide 1-Fws") +
  guides(fill = FALSE, color = FALSE) +
  stat_compare_means(comparisons = combos_pairwise_region_fws_province, aes(label = after_stat(p.signif)),
                     method = "wilcox.test")

reg_fws

ggsave("region_fws.png", reg_fws, width = 8, height = 6, bg = "white")


# # STATISTICAL ANALYSES
# Filter data by province and region
combined_data_province <- rbind(data.frame(He_province = processed_He_results[processed_He_results$geo == "province", ]$post_stat_mean,
                                           province = processed_He_results[processed_He_results$geo == "province", ]$population))

combined_data_region <- rbind(data.frame(He_region = processed_He_results[processed_He_results$geo == "region", ]$post_stat_mean,
                                         region = processed_He_results[processed_He_results$geo == "region", ]$population))

# 
# # CHECK FOR NORMALITY IN THE He DATA
# # PROVINCE
# # Create Q-Q plots
# ggplot(combined_data_province, aes(sample = He_province)) +
#   geom_qq() +
#   geom_qq_line() +
#   labs(title = "Q-Q Plot of He_province by Province",
#        x = "Theoretical Quantiles",
#        y = "Sample Quantiles") +
#   facet_grid(~ province, scales = "free")
# 
# # Perform Shapiro-Wilk
# shapiro_test_results_province <- combined_data_province %>%
#   group_by(province) %>%
#   summarize(
#     p_value = shapiro.test(He_province)$p.value
#   )
# 
# shapiro_test_results_province
# 
# #REGION
# # Create Q-Q plots
# ggplot(combined_data_region, aes(sample = He_region)) +
#   geom_qq() +
#   geom_qq_line() +
#   labs(title = "Q-Q Plot of He_region by Province and Year",
#        x = "Theoretical Quantiles",
#        y = "Sample Quantiles") +
#   facet_grid(~region , scales = "free")
# 
# # Perform Shapiro-Wilk test
# shapiro_test_results_region <- combined_data_region %>%
#   group_by(region) %>%
#   summarize(
#     p_value = shapiro.test(He_region)$p.value
#   )
# 
# shapiro_test_results_region
# 
# 
# ## He DATA IS NOT NORMAL, SO KRUSKAL-WALLIS (wilcox): 
# pairwise_province_He_2022 <- pairwise.wilcox.test(combined_data_province$He_province, 
#                                                   combined_data_province$province, p.adjust.method = "bonferroni")
# 
# p.val_province_He_2022 <- melt(pairwise_province_He_2022[[3]])
# signif_p.val_province_He_2022 <- p.val_province_He_2022[p.val_province_He_2022$value <0.05 & !is.na(p.val_province_He_2022$value),]
# 
# pairwise_region_He_2022 <- pairwise.wilcox.test(combined_data_region$He_region, 
#                                                 combined_data_region$region, p.adjust.method = "bonferroni")
# 
# p.val_region_He_2022 <- melt(pairwise_region_He_2022[[3]])
# signif_p.val_region_He_2022 <- p.val_region_He_2022[p.val_region_He_2022$value <0.05 & !is.na(p.val_region_He_2022$value),]


# # Changes in He distributions and means by year
# 
# provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
# regions <- c("North", "Centre", "South")
# 
# combined_data_province$province <- factor(combined_data_province$province, levels = provinces)
# combined_data_region$region <- factor(combined_data_region$region, levels = regions)
# 
# #provinces
# mean_data <- combined_data_province %>%
#   group_by(province) %>%
#   summarize(mean_He = mean(He_province, na.rm = TRUE),
#             median_He = median(He_province, na.rm = TRUE))
# 
# p1 <- ggplot(combined_data_province, aes(x = He_province)) +
#   geom_histogram(alpha = 0.7, bins = 30) +
#   geom_vline(data = mean_data, aes(xintercept = mean_He), linetype = "solid", color = "limegreen") +
#   geom_vline(data = mean_data, aes(xintercept = median_He), linetype = "solid", color = "orange") +
#   labs(title = "Province Heterozygosity",
#        x = "Genome-wide He",
#        y = "Frequency") +
#   facet_wrap(~ province) +
#   theme_minimal()
# 
# p1
# 
# p2 <- ggplot(combined_data_province, aes(x = province, y = He_province, fill = province)) +
#   geom_violin(width = 1, aes(color = province), alpha = 0.4) +
#   geom_boxplot(width = 0.1, aes(color = province), fill = "white", alpha = 0.4) +
#   labs(x = "", y = "Genome-wide He") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Adjust x-axis label angle
#   guides(fill = FALSE, color = FALSE)
# 
# p2
# 
# library(cowplot)
# combined_plot_prov <- plot_grid(p1, p2, ncol = 2)
# 
# ggsave("province_He.png", combined_plot_prov, width = 16, height = 8, bg = "white")
# 
# #regions
# mean_data_region <- combined_data_region %>%
#   group_by(region) %>%
#   summarize(mean_He = mean(He_region, na.rm = TRUE),
#             median_He = median(He_region, na.rm = TRUE))
# 
# p3 <- ggplot(combined_data_region, aes(x = He_region)) +
#   geom_histogram(alpha = 0.7, bins = 30) +
#   geom_vline(data = mean_data_region, aes(xintercept = mean_He), linetype = "solid", color = "limegreen") +
#   geom_vline(data = mean_data_region, aes(xintercept = median_He), linetype = "solid", color = "orange") +
#   labs(title = "Region Heterozygosity",
#        x = "He",
#        y = "Frequency") +
#   facet_wrap(~ region) +
#   theme_minimal()
# 
# p3
# 
# p4 <- ggplot(combined_data_region, aes(x = region, y = He_region, fill = region)) +
#   geom_violin(width = 1, aes(color = region), alpha = 0.4) +
#   geom_boxplot(width = 0.1, aes(color = region), fill = "white", alpha = 0.4) +
#   labs(x = "", y = "Mean Genome-wide He") +
#   theme_minimal()+
#   guides(fill = FALSE, color =FALSE)
# 
# p4
# 
# combined_plot_region <- plot_grid(p3, p4, ncol = 2)
# 
# ggsave("region_He.png", combined_plot_region, width = 16, height = 8, bg = "white")
# 


## linear mixed models to assess difference in He (from Nanna's scripts)
# LLM (EACH SITE AS REFERENCE)
library(nlme)

# Define the population levels
provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")


###
# Function to fit LLM and collect results for each reference level
fit_and_collect_results <- function(reference_population, pop = "region") {
  
  # Subset data and set reference population
  he_province <- processed_He_results[processed_He_results$geo == pop,] %>%
    mutate(population = factor(population, levels = population_levels)) %>%
    mutate(population = relevel(population, reference_population))
  
  # Fit LLM
  he.model.province <- lme(post_stat_med ~ population,
                           random = ~ 1 | locus,
                           data = he_province,
                           na.action = na.omit)

  # Extract and return summary statistics
  summary_data <- summary(he.model.province)
  aic <- AIC(logLik(he.model.province))
  t_table <- summary_data$tTable
  ci_mod <- intervals(he.model.province, which = "fixed")
  p_values <- anova(he.model.province, type = "marginal")$'p-value'[2]
  
  # Extract CI for reference population
  ci_reference <- ci_mod$fixed[1,]
  ci_reference <- t(as.data.frame(ci_reference))
  rownames(ci_reference) <- reference_population
  
  # Extract CIs for other populations
  ci_others <- t(sapply(2:length(population_levels), function(i) ci_reference + ci_mod$fixed[i,]))
  rownames(ci_others) <- population_levels[!population_levels %in% reference_population]
  
  ci_all <- rbind(ci_reference, ci_others)
  
  #merge tables
  t_table <- cbind(t_table, ci_all)
  rownames(t_table) <- rownames(ci_all)
  t_table <- as.data.frame(t_table)
  t_table$reference_pop <- reference_population
  t_table$compared_pop <- rownames(t_table)
  
  # Combine results
  results <- list(reference_population = reference_population,
                  t_table = t_table)
  
  return(results)
}

###

analyze_results <- function(results_list, pop = "region") {
  # edit list
  names(results_list) <- population_levels
  
  results_list <- lapply(results_list, function(sublist) {
    # Remove the reference_population element from each sublist
    sublist$reference_population <- NULL
    return(sublist)
  })
  
  
  llm_tables_list <- list()
  
  # Extract t_table from each element in results_list
  for (population_name in names(results_list)) {
    t_table <- results_list[[population_name]]$t_table
    llm_tables_list[[population_name]] <- t_table
  }
  
  # Combine all t_tables into a single dataframe
  llm_tables_all <- do.call(rbind, llm_tables_list)
  
  llm_tables_all$lower <- round(llm_tables_all$lower, 3)
  llm_tables_all$est. <- round(llm_tables_all$est., 3)
  llm_tables_all$upper <- round(llm_tables_all$upper, 3)
  
  llm_tables_all <- llm_tables_all[llm_tables_all$reference_pop != llm_tables_all$compared_pop,]
  
  #plot He and 95% CI for each pop
  CIs <- llm_tables_all[c("compared_pop", "lower", "est.", "upper")]
  rownames(CIs) <- NULL
  
  CIs <- unique(CIs)
  
  CIs$compared_pop <- factor(CIs$compared_pop, levels = rev(population_levels))
  
  # Create a ggplot object
  plotci <- ggplot(CIs, aes(x = compared_pop, y = est.)) +
    # Add the center point
    geom_point(color = "black") +
    # Add the error bars for confidence intervals
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "black") +
    # Adjust the appearance
    labs(title = "",
         x = "",
         y = "He Estimate") +
    theme_minimal()+
    coord_flip()
  
  plotci
  
  #put zero to self comparisons
  #llm_tables_all[llm_tables_all$`p-value` < 0.0000001, ]$`p-value` <- 0
  
  #significant pairwise comparions
  llm_tables_significant_pariwise<- llm_tables_all[llm_tables_all$reference_pop != llm_tables_all$compared_pop &  llm_tables_all$`p-value`< 0.05,]
  llm_tables_significant_pariwise <- subset(llm_tables_significant_pariwise, !duplicated(t(apply(llm_tables_significant_pariwise[c("compared_pop", "reference_pop")], 1, sort))))
  
  alpha <- 0.05
  sig_levels <- c("***", "**", "*", "")
  
  # Define significance labels based on p-values
  llm_tables_significant_pariwise$significance <- cut(llm_tables_significant_pariwise$`p-value`, 
                                                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                                      labels = sig_levels,
                                                      include.lowest = TRUE)
  
  #AIC AND ANOVA'S P-VALUE (it's the same for every comparison)
  
  # Subset data and set reference population
  he_province <- processed_He_results[processed_He_results$geo == pop ,] %>%
    mutate(population = factor(population, levels = population_levels))
  
  he.model.province <- lme(post_stat_med ~ population,
                           random = ~ 1 | locus,
                           data = he_province,
                           na.action = na.omit)
  
  anovap <- anova(he.model.province, type = "marginal")
  aicval <- AIC(logLik(he.model.province))
  
  
  # Combine results
  results <- list(plotCI = plotci,
                  llm_tables_significant_pariwise = llm_tables_significant_pariwise,
                  anovap = anovap,
                  AIC = aicval)
 
  return(results)
}
###

# PROVINCES
population_levels <- provinces
results_list_provinces <- lapply(provinces, fit_and_collect_results, pop = "province")
results_provinces <- analyze_results(results_list_provinces, pop = "province")

#PREGIONS
population_levels <- regions
results_list_regions <- lapply(regions, fit_and_collect_results, pop = "region")
results_regions <- analyze_results(results_list_regions, pop = "region")


#######################################3
# 10.- pairwise FST  (https://biology.stackexchange.com/questions/40756/calculating-pairwise-fst-from-allele-frequencies)
#########################################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

### no need to remove DRY season pops from region analysis because it already was removed qhen running moire by population 


# pairwise FST comparion = 1 - Hs/HT.... Hs: heterozygosity of each loci per population; HT: average heterozygosity between pop1 and 2 for each loci. heatmap AND distributions (mean tests also)

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(path = "FINAL_MOIRE_RESULTS", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)
rds_files <- rds_files[!rds_files %in% "FINAL_MOIRE_RESULTS/all_samples_complete_filtered_MOIRE-RESULTS_2022_only_FOR_MOI.RDS"] 


# Load each RDS file into the list with the file name as the list name
He_results_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  He_results_list[[file_name]] <- readRDS(file)
}

# Loop through each element in He_results_list
processed_He_results <- data.frame()

for (i in seq_along(He_results_list)) {
  
  # Summarize He
  He_results <- moire::summarize_he(He_results_list[[i]])
  He_results$population <- names(He_results_list[i])
  
  # Add the processed He results to the list
  processed_He_results <- rbind(processed_He_results, He_results)
}

#formatting categories
processed_He_results$population <- gsub("_2022_MOIRE-RESULTS_FOR_ALLELE_FREQS", "", processed_He_results$population)

library(stringr)
processed_He_results <- processed_He_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))

######
# keep amplicons with high He:

he_amps <- processed_He_results %>%
  group_by(locus) %>%
  summarize(mean = mean(post_stat_mean)) %>%
  arrange(desc(mean))

#keep top n% amplicons with highest He (or use a numeric threshold like > 0.1 of He or whatever... test.)
perc_25<- round(length(unique(he_amps$locus))*0.95)
he_amps_top50 <- he_amps[1:perc_25,]

# FILTER
processed_He_results <- processed_He_results[processed_He_results$locus %in% he_amps_top50$locus,]
#####


#separate provinces and regions
processed_He_results_provinces <- processed_He_results[processed_He_results$geo == "province", ]
processed_He_results_provinces$pop <- processed_He_results_provinces$population
processed_He_results_regions<- processed_He_results[processed_He_results$geo == "region", ]
processed_He_results_regions$pop <- processed_He_results_regions$population


#sample sizes for each population

sample_size_provinces <- combined_df_merged %>%
  group_by(year, province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces$pop <- sample_size_provinces$province
sample_size_provinces

# sample_size_regions <- combined_df_merged %>%
#   group_by(year, region) %>%
#   summarise(unique_NIDA2_count = n_distinct(NIDA2))
# 
# sample_size_regions$pop <- sample_size_regions$region
# sample_size_regions

combined_df_merged_nodry <- combined_df_merged %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions$pop <- sample_size_regions$region
sample_size_regions

# CALCULATE pairwise Fst ## STILL NOT COMPLETELY SURE ABOUT THIS, BUT HEY...
# Fst = (Ht - Hs) / Ht [same as 1 - Hs / Ht]*;
# Ht is (n1*Hs1 + n2*Hs2) / n1+n2 [n is individuals of each pop]; 
# Hs = Hs1 + Hs2 / 2 [Hs1 and Hs2 are post_stat_mean (heterozygosity for each locus) calculated by moire for each pop]. calculated genome-wide average, just avergaae Hs1 and Hs2 of all loci
# THIS IS DONE FOR EACH LOCUS

# 1.- calculate Hs and Ht for each locus
# 2.- calculate Fst [(Ht - Hs) / Ht] through a Linear Mixed Model


### 1) FOR REGIONS

fts_input_regions <- merge(processed_He_results_regions[c("locus", "post_stat_mean", "pop")], sample_size_regions[c("unique_NIDA2_count", "pop")], by = c("pop"))

# Generate all possible combinations of unique values
unique_pops <- unique(fts_input_regions$pop)
combinations <- expand.grid(pop1 = unique_pops, pop2 = unique_pops)

# Create an empty dataframe to store the results
fst_results_df <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(fst_results_df) <- c("pop1", "pop2", "Hs1", "Hs2", "n1", "n2", "locus")

# Loop through each locus
for (locus in unique(fts_input_regions$locus)) {
  locus_subset <- fts_input_regions[fts_input_regions$locus == locus,]
  
  # Loop through each combination of populations
  for (i in 1:nrow(combinations)) {
    pop1 <- as.character(combinations$pop1[i])
    pop2 <- as.character(combinations$pop2[i])
    
    Hs1 <- locus_subset[locus_subset$pop == pop1, ]$post_stat_mean
    n1 <- locus_subset[locus_subset$pop == pop1, ]$unique_NIDA2_count
    
    Hs2 <- locus_subset[locus_subset$pop == pop2, ]$post_stat_mean
    n2 <- locus_subset[locus_subset$pop == pop2, ]$unique_NIDA2_count
    
    # Check if both populations have data for this locus
    if (length(Hs1) > 0 && length(Hs2) > 0) {
      row <- c(pop1, pop2, Hs1, Hs2, n1, n2, locus)
    } else {
      # If one population does not have data, assign NA to the corresponding entry
      row <- c(pop1, pop2, NA, NA, NA, NA, locus)
    }
    
    fst_results_df <- rbind(fst_results_df, row)
  }
}

#formating
fst_results_df <- fst_results_df[-1,]
fst_results_df$Hs1 <- round(as.numeric(fst_results_df$Hs1), 3)
fst_results_df$Hs2 <- round(as.numeric(fst_results_df$Hs2),3)
fst_results_df$n1 <- round(as.numeric(fst_results_df$n1),3)
fst_results_df$n2 <- round(as.numeric(fst_results_df$n2),3)

if (nrow(fst_results_df) == nrow(combinations)*length(unique(fts_input_regions$locus))){
  print("Got all expected combinations for per locus pairwise Fst calculations")
}else{
  print("grab a coffee.")
}

fst_results_df <- fst_results_df[complete.cases(fst_results_df), ] # remove NA rows (those on which the amplicon wasn't present in both populations, hence Fst can't be calculated)

#calculate Hs for populations PER LOCUS
fst_results_df$Hs <- (fst_results_df$Hs1 + fst_results_df$Hs2) / 2
fst_results_df$Hs <- round(as.numeric(fst_results_df$Hs),3)

#Calculate Ht (uses sample sizes for each pop)
fst_results_df$Ht <- ((fst_results_df$n1 * fst_results_df$Hs1)+ (fst_results_df$n2 * fst_results_df$Hs2))/(fst_results_df$n1 + fst_results_df$n2)
fst_results_df$Ht <- round(as.numeric(fst_results_df$Ht),3)

# Define a function to calculate FST
calculate_FST <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # Extract pre-calculated Hs and Ht values
  Hs <- sampled_data$Hs
  Ht <- sampled_data$Ht
  
  # Calculate FST
  FST <- ((Ht - Hs) / Ht)
  
  return(FST)
}

#for each locus, just in case it's needed
fst_results_df$Fst <- (fst_results_df$Ht - fst_results_df$Hs)/fst_results_df$Ht


#llm (interchangeable with boostrat analysis)
library(nlme)

FST_LLM <- as.data.frame(cbind(pop1 =fst_results_df$pop1,
                               pop2 = fst_results_df$pop2,
                               comparison = paste0(fst_results_df$pop1, "_",fst_results_df$pop2),
                               locus = fst_results_df$locus,
                               fst = as.numeric(fst_results_df$Fst)))


FST_LLM$fst <- as.numeric(FST_LLM$fst)

fst.model.region <- lme(fst ~ comparison,
                        random = ~ 1 | locus,
                        data = FST_LLM,
                        na.action = na.omit)

summary_data <- summary(fst.model.region)
summary_table <- as.data.frame(summary_data$tTable)
summary_table <- unique(summary_table)
summary_table$comparison <-  rownames(summary_table)

dim(summary_table)

cis <- intervals(fst.model.region, which = "fixed")
cis <- as.data.frame(cis$fixed)
cis <- unique(cis)
cis$comparison <-  rownames(cis)

dim(cis)

#merge estimates with CIs
final_table<- merge(cis, summary_table, by = c("comparison"))
final_table$comparison <- gsub("comparison", "", final_table$comparison)
final_table <- unique(merge(final_table, FST_LLM[c("pop1", "pop2", "comparison")], by = c("comparison")))

#heatmap
final_table$He_estimate<- round(final_table$est., 3)
final_table <- complete(final_table, pop1, pop2, fill = list(He_estimate = 0, `p-value` = 1))
final_table$He_estimate<- ifelse(final_table$He_estimate < 0, 0, final_table$He_estimate) #TURN NEGATIVE VALUES TO 0
final_table$label <- ifelse(final_table$`p-value` < 0.05 & final_table$He_estimate > 0 , paste0(final_table$He_estimate, "*"), as.character(final_table$He_estimate))

regions <- c("North", "Centre", "South") #ordered from north to south
regions <- rev(regions)

final_table$pop1 <- factor(final_table$pop1, levels = regions)
final_table$pop2 <- factor(final_table$pop2, levels = regions)

heatmap_regions <- ggplot(final_table, aes(x = pop2, y = pop1, fill = He_estimate, label = label)) +
  geom_tile() +
  geom_text(color = "black") +
  scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(final_table$He_estimate), max(final_table$He_estimate))) +  # Adjust scale limits
  theme_minimal() +
  labs(x = "", y = "")

heatmap_regions

#significance
final_table$significance <- ifelse(final_table$`p-value` < 0.05, "p < 0.05", "not_signiff.")

final_table <- final_table %>%
  arrange(`est.`)

#keep fst > 0
final_table <- final_table[final_table$est. > 0.000001,]

# Create logical index to keep every other row
keep_rows <- seq(nrow(final_table)) %% 2 == 1

# Subset final_table to keep every other row
final_table <- final_table[keep_rows, ]

final_table <- final_table %>%
  mutate(comparison = factor(comparison, levels = comparison[order(est.)]))

fst_regions <- ggplot(na.omit(final_table), aes(x = comparison, y = est., color = significance)) +
  # Add the center point
  geom_point() +
  # Add the error bars for confidence intervals
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # Adjust the appearance
  labs(title = "",
       x = "Pairwise comparisons",
       y = "Fst Estimate") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()+
  coord_flip()

fst_regions

anovap <- anova(fst.model.region, type = "marginal")
aicval <- AIC(logLik(fst.model.region))

ggsave("fst_heatmap_regions.png", heatmap_regions, width = 8, height = 6, bg = "white")
ggsave("fst_CI_regions.png", fst_regions, width = 8, height = 6, bg = "white")


### 2) FOR PROVINCES

fts_input_provinces <- merge(processed_He_results_provinces[c("locus", "post_stat_mean", "pop")], sample_size_provinces[c("unique_NIDA2_count", "pop")], by = c("pop"))

# Generate all possible combinations of unique values
unique_pops <- unique(fts_input_provinces$pop)
combinations <- expand.grid(pop1 = unique_pops, pop2 = unique_pops)

# Create an empty dataframe to store the results
fst_results_df <- data.frame(matrix(ncol = 7, nrow = 1))
colnames(fst_results_df) <- c("pop1", "pop2", "Hs1", "Hs2", "n1", "n2", "locus")

# Loop through each locus
for (locus in unique(fts_input_provinces$locus)) {
  locus_subset <- fts_input_provinces[fts_input_provinces$locus == locus,]
  
  # Loop through each combination of populations
  for (i in 1:nrow(combinations)) {
    pop1 <- as.character(combinations$pop1[i])
    pop2 <- as.character(combinations$pop2[i])
    
    Hs1 <- locus_subset[locus_subset$pop == pop1, ]$post_stat_mean
    n1 <- locus_subset[locus_subset$pop == pop1, ]$unique_NIDA2_count
    
    Hs2 <- locus_subset[locus_subset$pop == pop2, ]$post_stat_mean
    n2 <- locus_subset[locus_subset$pop == pop2, ]$unique_NIDA2_count
    
    # Check if both populations have data for this locus
    if (length(Hs1) > 0 && length(Hs2) > 0) {
      row <- c(pop1, pop2, Hs1, Hs2, n1, n2, locus)
    } else {
      # If one population does not have data, assign NA to the corresponding entry
      row <- c(pop1, pop2, NA, NA, NA, NA, locus)
    }
    
    fst_results_df <- rbind(fst_results_df, row)
  }
}

#formating
fst_results_df <- fst_results_df[-1,]
fst_results_df$Hs1 <- round(as.numeric(fst_results_df$Hs1), 3)
fst_results_df$Hs2 <- round(as.numeric(fst_results_df$Hs2),3)
fst_results_df$n1 <- round(as.numeric(fst_results_df$n1),3)
fst_results_df$n2 <- round(as.numeric(fst_results_df$n2),3)

if (nrow(fst_results_df) == nrow(combinations)*length(unique(fts_input_provinces$locus))){
  print("Got all expected combinations for per locus pairwise Fst calculations")
}else{
  print("grab a coffee.")
}


fst_results_df <- fst_results_df[complete.cases(fst_results_df), ] # remove NA rows (those on which the amplicon wasn't present in both populations, hence Fst can't be calculated)

#calculate Hs for populations PER LOCUS
fst_results_df$Hs <- (fst_results_df$Hs1 + fst_results_df$Hs2) / 2
fst_results_df$Hs <- round(as.numeric(fst_results_df$Hs),3)

#Calculate Ht (uses sample sizes for each pop)
fst_results_df$Ht <- ((fst_results_df$n1 * fst_results_df$Hs1)+ (fst_results_df$n2 * fst_results_df$Hs2))/(fst_results_df$n1 + fst_results_df$n2)
fst_results_df$Ht <- round(as.numeric(fst_results_df$Ht),3)

# Define a function to calculate FST
calculate_FST <- function(data, indices) {
  sampled_data <- data[indices, ]
  
  # Extract pre-calculated Hs and Ht values
  Hs <- sampled_data$Hs
  Ht <- sampled_data$Ht
  
  # Calculate FST
  FST <- ((Ht - Hs) / Ht)
  
  return(FST)
}

#for each locus, just in case it's needed
fst_results_df$Fst <- (fst_results_df$Ht - fst_results_df$Hs)/fst_results_df$Ht


#llm (interchangeable with boostrat analysis)
library(nlme)

FST_LLM <- as.data.frame(cbind(pop1 =fst_results_df$pop1,
                               pop2 = fst_results_df$pop2,
                               comparison = paste0(fst_results_df$pop1, "_",fst_results_df$pop2),
                               locus = fst_results_df$locus,
                               fst = as.numeric(fst_results_df$Fst)))

FST_LLM$fst <- as.numeric(FST_LLM$fst)

fst.model.region <- lme(fst ~ comparison,
                        random = ~ 1 | locus,
                        data = FST_LLM,
                        na.action = na.omit)

summary_data <- summary(fst.model.region)
summary_table <- as.data.frame(summary_data$tTable)
summary_table <- unique(summary_table)
summary_table$comparison <-  rownames(summary_table)

dim(summary_table)

cis <- intervals(fst.model.region, which = "fixed")
cis <- as.data.frame(cis$fixed)
cis <- unique(cis)
cis$comparison <-  rownames(cis)

dim(cis)

#merge estimates with CIs
final_table<- merge(cis, summary_table, by = c("comparison"))
final_table$comparison <- gsub("comparison", "", final_table$comparison)
final_table <- unique(merge(final_table, FST_LLM[c("pop1", "pop2", "comparison")], by = c("comparison")))

#heatmap
final_table$He_estimate<- round(final_table$est., 3)
final_table <- complete(final_table, pop1, pop2, fill = list(He_estimate = 0, `p-value` = 1))
final_table$He_estimate<- ifelse(final_table$He_estimate < 0, 0, final_table$He_estimate) #TURN NEGATIVE VALUES TO 0
final_table$label <- ifelse(final_table$`p-value` < 0.05 & final_table$He_estimate > 0 , paste0(final_table$He_estimate, "*"), as.character(final_table$He_estimate))

provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
provinces <- rev(provinces)


final_table$pop1 <- factor(final_table$pop1, levels = provinces)
final_table$pop2 <- factor(final_table$pop2, levels = provinces)

heatmap_provinces <- ggplot(final_table, aes(x = pop2, y = pop1, fill = He_estimate, label = label)) +
  geom_tile() +
  geom_text(color = "black") +
  scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(final_table$He_estimate), max(final_table$He_estimate))) +  # Adjust scale limits
  theme_minimal() +
  labs(x = "", y = "")

heatmap_provinces

#significance
final_table$significance <- ifelse(final_table$`p-value` < 0.05, "p < 0.05", "not_signiff.")

final_table <- final_table %>%
  arrange(`est.`)

#keep fst > 0 
final_table <- final_table[final_table$est. > 0.000001,]

# Create logical index to keep every other row
keep_rows <- seq(nrow(final_table)) %% 2 == 1

# Subset final_table to keep every other row
final_table <- final_table[keep_rows, ]

final_table <- final_table %>%
  mutate(comparison = factor(comparison, levels = comparison[order(est.)]))

fst_provinces <- ggplot(na.omit(final_table), aes(x = comparison, y = est., color = significance)) +
  # Add the center point
  geom_point() +
  # Add the error bars for confidence intervals
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  # Adjust the appearance
  labs(title = "",
       x = "Pairwise comparisons",
       y = "Fst Estimate") +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal()+
  coord_flip()

fst_provinces

anovap <- anova(fst.model.region, type = "marginal")
aicval <- AIC(logLik(fst.model.region))

ggsave("fst_heatmap_provinces.png", heatmap_provinces, width = 12, height = 10, bg = "white")
ggsave("fst_CI_provinces.png", fst_provinces, width = 8, height = 6, bg = "white")



# #bootstraping analysis (interchangeable with llm) USING THIS
# library(boot)
# library(boot.pval)
# 
# # Create a list to store bootstrap results
# bootstrap_results <- list()
# 
# # Loop through each pairwise comparison
# for (i in 1:nrow(combinations)) {
#   pop1 <- as.character(combinations$pop1[i])
#   pop2 <- as.character(combinations$pop2[i])
#   
#   # Subset data for the current pairwise comparison
#   pairwise_data <- fst_results_df[fst_results_df$pop1 == pop1 & fst_results_df$pop2 == pop2, ]
#   
#   # Define a unique name for the pairwise comparison
#   comparison_name <- paste(pop1, pop2, sep = "_vs_")
#   
#   # Perform bootstrap resampling
#   bootstrap_results[[comparison_name]] <- boot(pairwise_data, calculate_FST, R = 10000)
# }
# 
# 
# # Calculate the mean of bootstrap replicates for each pairwise comparison
# mean_FST <- sapply(bootstrap_results, function(result) mean(result$t0))
# 
# p_value <- data.frame()
# 
# for (i in seq(1:length(bootstrap_results))){
#   p <- boot.pval(bootstrap_results[[i]], type = "perc", theta_null = 0)
#   p_value <- rbind(p_value, p)
# }
# 
# colnames(p_value)<- "p_value"
# 
# # Create a function to extract confidence intervals for the mean
# get_mean_confidence_intervals <- function(bootstrap_result) {
#   ci <- boot.ci(bootstrap_result, type = "perc")
#   lower_limit <- ci$percent[4]
#   upper_limit <- ci$percent[5]
#   return(c(lower_limit, upper_limit))
# }
# 
# # Extract confidence intervals for the mean of FST
# mean_confidence_intervals <- lapply(bootstrap_results, get_mean_confidence_intervals)
# 
# # Create a data frame to store the results
# mean_FST_df <- data.frame(
#   pairwise_comparisons = names(mean_FST),
#   mean_FST = ifelse(is.na(mean_FST), 0, mean_FST),
#   lower_limit = sapply(mean_confidence_intervals, `[`, 1),
#   upper_limit = sapply(mean_confidence_intervals, `[`, 2),
#   p_val = p_value
# )
# 
# # Remove "_2022" from population names and reorder according to the specified order
# #mean_FST_df$pairwise_comparisons <- gsub("_2022", "", mean_FST_df$pairwise_comparisons)
# 
# # Exclude rows with 0 in the last three columns
# mean_FST_df_filtered <- subset(mean_FST_df, lower_limit != 0 | upper_limit != 0 | mean_FST != 0)
# 
# #exclude repeated comparisons: CAREFUL!! IF MEANS FOR DIFFERENT COMPARISONS ARE THE SAME IT WILL ELIMINATE ONE. not the best method.
# mean_FST_df_filtered<- mean_FST_df_filtered %>%
#   distinct(mean_FST, .keep_all = TRUE)
# 
# # Reorder pairwise_comparisons based on mean_FST
# mean_FST_df_filtered$pairwise_comparisons <- factor(mean_FST_df_filtered$pairwise_comparisons, 
#                                                     levels = mean_FST_df_filtered$pairwise_comparisons[order(mean_FST_df_filtered$mean_FST)])
# 
# # Add a column indicating significance
# mean_FST_df_filtered$significance <- ifelse(mean_FST_df_filtered$p_value < 0.05, "p < 0.05", "not_signiff.")
# 
# # Plot with significance indication
# fst_regions <- ggplot(mean_FST_df_filtered, aes(x = pairwise_comparisons, y = mean_FST, color = significance)) +
#   geom_boxplot() +
#   geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = 0.2) +
#   scale_color_manual(values = c("red")) +  # Set colors for significance levels
#   labs(x = "Pairwise Comparisons", y = "Genome-wide Fst") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme_minimal() +
#   coord_flip()
# 
# fst_regions
# 
# ggsave("fst_CI_regions.png", fst_regions, width = 8, height = 6, bg = "white")

# # Create heatmap for regions comparisons
# create_heatmap <- function(data, title) {
#   
#   data$mean_FST<- round(data$mean_FST, 3)
#   data$mean_FST<- ifelse(data$mean_FST < 0, 0, data$mean_FST) #TURN NEGATIVE VALUES TO 0
#   data$label <- ifelse(data$p_val < 0.05 & data$mean_FST > 0 , paste0(data$mean_FST, "*"), as.character(data$mean_FST))
#   
#   ggplot(data, aes(x = pop2, y = pop1, fill = mean_FST, label = label)) +
#     geom_tile() +
#     geom_text(color = "black") +
#     scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(data$mean_FST), max(data$mean_FST))) +  # Adjust scale limits
#     theme_minimal() +
#     labs(x = "", y = "")
# }
# 
# # Split pairwise_comparisons into pop1 and pop2
# library(tidyr)
# mean_FST_df <- separate(mean_FST_df, pairwise_comparisons, into = c("pop1", "pop2"), sep = "_vs_")
# 
# # Reorder pop1 and pop2 columns based on regions order
# regions <- c("North", "Centre", "South") #ordered from north to south
# regions <- rev(regions)
# 
# mean_FST_df$pop1 <- factor(mean_FST_df$pop1, levels = regions)
# mean_FST_df$pop2 <- factor(mean_FST_df$pop2, levels = regions)
# 
# heatmap_2022_regions <- create_heatmap(mean_FST_df)
# 
# heatmap_2022_regions
# 
# ggsave("fst_heatmap_regions.png", heatmap_2022_regions, width = 8, height = 6, bg = "white")


# #bootstraping analysis (interchangeable with llm) USING THIS
# library(boot)
# library(boot.pval)
# 
# # Create a list to store bootstrap results
# bootstrap_results <- list()
# 
# # Loop through each pairwise comparison
# for (i in 1:nrow(combinations)) {
#   pop1 <- as.character(combinations$pop1[i])
#   pop2 <- as.character(combinations$pop2[i])
#   
#   # Subset data for the current pairwise comparison
#   pairwise_data <- fst_results_df[fst_results_df$pop1 == pop1 & fst_results_df$pop2 == pop2, ]
#   
#   # Define a unique name for the pairwise comparison
#   comparison_name <- paste(pop1, pop2, sep = "_vs_")
#   
#   # Perform bootstrap resampling
#   bootstrap_results[[comparison_name]] <- boot(pairwise_data, calculate_FST, R = 10000)
# }
# 
# 
# # Calculate the mean of bootstrap replicates for each pairwise comparison
# mean_FST <- sapply(bootstrap_results, function(result) mean(result$t0))
# 
# p_value <- data.frame()
# 
# for (i in seq(1:length(bootstrap_results))){
#   p <- boot.pval(bootstrap_results[[i]], type = "perc", theta_null = 0)
#   p_value <- rbind(p_value, p)
# }
# 
# colnames(p_value)<- "p_value"
# 
# # Create a function to extract confidence intervals for the mean
# get_mean_confidence_intervals <- function(bootstrap_result) {
#   ci <- boot.ci(bootstrap_result, type = "perc")
#   lower_limit <- ci$percent[4]
#   upper_limit <- ci$percent[5]
#   return(c(lower_limit, upper_limit))
# }
# 
# # Extract confidence intervals for the mean of FST
# mean_confidence_intervals <- lapply(bootstrap_results, get_mean_confidence_intervals)
# 
# # Create a data frame to store the results
# mean_FST_df <- data.frame(
#   pairwise_comparisons = names(mean_FST),
#   mean_FST = ifelse(is.na(mean_FST), 0, mean_FST),
#   lower_limit = sapply(mean_confidence_intervals, `[`, 1),
#   upper_limit = sapply(mean_confidence_intervals, `[`, 2),
#   p_val = p_value
# )
# 
# # Remove "_2022" from population names and reorder according to the specified order
# #mean_FST_df$pairwise_comparisons <- gsub("_2022", "", mean_FST_df$pairwise_comparisons)
# 
# # Exclude rows with 0 in the last three columns
# mean_FST_df_filtered <- subset(mean_FST_df, lower_limit != 0 | upper_limit != 0 | mean_FST != 0)
# 
# #exclude repeated comparisons: CAREFUL!! IF MEANS FOR DIFFERENT COMPARISONS ARE THE SAME IT WILL ELIMINATE ONE. not the best method.
# mean_FST_df_filtered<- mean_FST_df_filtered %>%
#   distinct(mean_FST, .keep_all = TRUE)
# 
# # Reorder pairwise_comparisons based on mean_FST
# mean_FST_df_filtered$pairwise_comparisons <- factor(mean_FST_df_filtered$pairwise_comparisons, 
#                                                     levels = mean_FST_df_filtered$pairwise_comparisons[order(mean_FST_df_filtered$mean_FST)])
# # Add a column indicating significance
# mean_FST_df_filtered$significance <- ifelse(mean_FST_df_filtered$p_value < 0.05, "p < 0.05", "not_signiff.")

# Plot with significance indication
# fst_provinces <- ggplot(mean_FST_df_filtered, aes(x = pairwise_comparisons, y = mean_FST, color = significance)) +
#   geom_boxplot() +
#   geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = 0.2) +
#   scale_color_manual(values = c("black", "red")) +  # Set colors for significance levels
#   labs(x = "Pairwise Comparisons", y = "Genome-wide Fst") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme_minimal() +
#   coord_flip()
# 
# fst_provinces
# 
# ggsave("fst_CI_provinces.png", fst_provinces, width = 8, height = 10, bg = "white")

# # Create heatmap for regions comparisons
# create_heatmap <- function(data, title) {
#   
#   data$mean_FST<- round(data$mean_FST, 3)
#   data$mean_FST<- ifelse(data$mean_FST < 0, 0, data$mean_FST) #TURN NEGATIVE VALUES TO 0
#   data$label <- ifelse(data$p_val < 0.05 & data$mean_FST > 0 , paste0(data$mean_FST, "*"), as.character(data$mean_FST))
#   
#   ggplot(data, aes(x = pop2, y = pop1, fill = mean_FST, label = label)) +
#     geom_tile() +
#     geom_text(color = "black") +
#     scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(data$mean_FST), max(data$mean_FST))) +  # Adjust scale limits
#     theme_minimal() +
#     labs(x = "", y = "")
# }
# 
# # Split pairwise_comparisons into pop1 and pop2
# library(tidyr)
# mean_FST_df <- separate(mean_FST_df, pairwise_comparisons, into = c("pop1", "pop2"), sep = "_vs_")
# 
# # Reorder pop1 and pop2 columns based on provinces order
# provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
# provinces <- rev(provinces)
# 
# mean_FST_df$pop1 <- factor(mean_FST_df$pop1, levels = provinces)
# mean_FST_df$pop2 <- factor(mean_FST_df$pop2, levels = provinces)
# 
# heatmap_2022_provinces <- create_heatmap(mean_FST_df)
# heatmap_2022_provinces
# 
# ggsave("heatmap_fst_provinces.png", heatmap_2022_provinces, width = 12, height = 10, bg = "white")


########################
# 11. Ordination, population structure
########################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

#add VOC
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_doub_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhfr_tr_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_436_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_581_b")], by = "NIDA2")

  
combined_df_merged$VOC <- ifelse(combined_df_merged$dhps_doub_95_b == 1 & combined_df_merged$dhfr_tr_95_b == 1, "dhps_d_dhfr_tr",
                                 ifelse(combined_df_merged$dhfr_tr_95_b == 1 & combined_df_merged$dhps_doub_95_b == 0, "dhfr_tr",
                                        ifelse(combined_df_merged$dhfr_tr_95_b == 0 & combined_df_merged$dhps_doub_95_b == 1, "dhps_d", "WT")))

#no_genotype
combined_df_merged$VOC <- ifelse(is.na(combined_df_merged$VOC), "no_genotype", combined_df_merged$VOC)

#436 581 genotypes
# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_436
combined_df_merged$dhps_436 <- ifelse(combined_df_merged$dhps_436 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_436 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_436 == 2, "mix", "no_genotype")))

# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_581
combined_df_merged$dhps_581 <- ifelse(combined_df_merged$dhps_581 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_581 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_581 == 2, "mix", "no_genotype")))
                                             
combined_df_merged$VOC_436_581 <- paste0("dhps_436", "-", combined_df_merged$dhps_436, "_", "581", "-", combined_df_merged$dhps_581)

combined_df_merged$VOC_436_581 <- gsub("dhps_436-NA_581-NA", "no_genotype", combined_df_merged$VOC_436_581 )
combined_df_merged$VOC_436_581 <- gsub("dhps_436-WT_581-WT", "WT", combined_df_merged$VOC_436_581 )

unique(combined_df_merged$VOC_436_581)


#input for multivariate analyses
raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
                                   year = combined_df_merged$year, 
                                   province = combined_df_merged$province,
                                   region = combined_df_merged$region,
                                   locus = combined_df_merged$locus,
                                   n.alleles = combined_df_merged$n.alleles,
                                   norm.reads.locus = combined_df_merged$norm.reads.locus,
                                   allele = paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar),
                                   run_id = combined_df_merged$run_id,
                                   VOC = combined_df_merged$VOC,
                                   VOC_436_581 = combined_df_merged$VOC_436_581))


#remove Dry
raref_input <- raref_input %>%
  filter(!grepl("Dry", province))


# Melt the data frame to convert it from wide to long format
melted <- melt(raref_input, id.vars = c("NIDA2", "norm.reads.locus"), measure.vars = "allele")
melted<-melted[,-3]

library(tidyr)

rearranged <- melted %>%
  pivot_wider(
    names_from = value,
    values_from = norm.reads.locus
  )

# format
rearranged <- as.data.frame(rearranged)
rownames(rearranged) <- rearranged$NIDA2
rearranged <- rearranged[, -1]
rearranged <- replace(rearranged, is.na(rearranged), 0)

#pca labels:
NIDA2 <-data.frame(NIDA2 = rownames(rearranged))
pca_labels<- combined_df_merged %>% distinct(NIDA2, year, province, region, VOC_436_581)

pca_labels <- pca_labels %>%
  filter(!grepl("Dry", province))

if (all(NIDA2$NIDA2 == pca_labels$NIDA2)){
  print("Order of categorical variables is ok.")
}else{
  "grab a coffee."
}

# format freq df
rearranged <- rearranged %>%
  mutate_all(as.numeric)

zero_cols <- sapply(rearranged, function(x) all(x == 0))
rearranged_filtered <- rearranged[, !zero_cols]

# Replace values greater than 0 with 1: MAKE IT PRESENCE/ABSENCE
rearranged_pres_abs <- rearranged %>%
    mutate_all(~ ifelse(. > 0, 1, .))

rearranged_pres_abs <- rearranged_pres_abs %>%
  mutate_all(as.numeric)
 

provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")

#PCA
pca_result <- prcomp(rearranged, scale. = F)

# Extract the principal component scores
pcs <- as.data.frame(pca_result$x)

# Combine principal component scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 =pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# Calculate percentage variance explained by each principal component
pc_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create the plot
af_pca<- ggplot(pcs_with_labels, aes(x = PC1, y = PC2, color = factor(pcs_with_labels$province), shape = factor(pcs_with_labels$VOC_436_581))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-sample Allele Frequencies",
       x = paste0("PC1: ", round(pc_variance[1], 2), "%\n"),
       y = paste0("PC2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal() +guides(fill = FALSE, color = FALSE, shape = FALSE)

af_pca



# Perform PCA prsence/absence
pca_result <- prcomp(rearranged_pres_abs, scale. = F)

# Extract the principal component scores
pcs <- as.data.frame(pca_result$x)

# Combine principal component scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 =pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# Calculate percentage variance explained by each principal component
pc_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create the plot
pa_pca<- ggplot(pcs_with_labels, aes(x = PC1, y = PC2, color = factor(pcs_with_labels$province), shape = factor(pca_labels$VOC_436_581))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-Sample Allele Presence/Absence",
       x = paste0("PC1: ", round(pc_variance[1], 2), "%\n"),
       y = paste0("PC2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal()

pa_pca

library(cowplot)
combined_plot_pca <- plot_grid(af_pca, pa_pca, ncol = 2)

ggsave("PCA_regions.png", combined_plot_pca, width = 16, height = 10, bg = "white")

#PCoA
library(vegan)  # For the vegdist() function
library(ape)    # For the pcoa() function
library(ggplot2)

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values / sum(pcoa_result$values) * 100, 2)
variance_explained_axis1 <- variance_explained$Eigenvalues[1]
variance_explained_axis2 <- variance_explained$Eigenvalues[2]

# Plot PCoA with variance explained in title
af_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-Sample Allele Frequecnies",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()+
  guides(fill = FALSE, color = FALSE, shape = FALSE)

af_pcoa

#PCoA presence/absence

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged_pres_abs, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
variance_explained <- round(pcoa_result$values / sum(pcoa_result$values) * 100, 2)
variance_explained_axis1 <- variance_explained$Eigenvalues[1]
variance_explained_axis2 <- variance_explained$Eigenvalues[2]

# Plot PCoA with variance explained in title
pa_pcoa <- ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-Sample Allele Presence/Absence",
       x = paste0("PCo 1: ", variance_explained_axis1, "%\n"),
       y = paste0("PCo 2: ", variance_explained_axis2, "%")) +
  theme_minimal()

pa_pcoa


combined_plot_pcoa <- plot_grid(af_pcoa, pa_pcoa, ncol = 2)

ggsave("PCoA_regions.png", combined_plot_pcoa, width = 16, height = 10, bg = "white")


#####TSNE
set.seed(69)
perplexity <- floor((nrow(rearranged_filtered) - 1) / 3) #highest possible, if needed
tsne_result_freqs <- Rtsne(as.matrix(rearranged_filtered), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = F, max_iter = 1e4, num_threads = 0, perplexity = 100)
set.seed(69)
tsne_result_pres_abs <- Rtsne(as.matrix(rearranged_pres_abs), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = F, max_iter = 1e4, num_threads = 0, perplexity = 100)

# Convert t-SNE results to data frame
tsne_data_freqs <- as.data.frame(tsne_result_freqs$Y)
tsne_data_freqs <- cbind(tsne_data_freqs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)
tsne_data_pres_abs <- as.data.frame(tsne_result_pres_abs$Y)
tsne_data_pres_abs <- cbind(tsne_data_pres_abs, province = pca_labels$province, region = pca_labels$region, VOC_436_581 = pca_labels$VOC_436_581)

# Order factors
tsne_data_freqs$province <- factor(tsne_data_freqs$province, levels = provinces)
tsne_data_freqs$region <- factor(tsne_data_freqs$region, levels = regions)
tsne_data_pres_abs$province <- factor(tsne_data_pres_abs$province, levels = provinces)
tsne_data_pres_abs$region <- factor(tsne_data_pres_abs$region, levels = regions)

# Plot t-SNE of freqs
af_tsne <- ggplot(tsne_data_freqs, aes(V1, V2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-Sample Allele Frequencies",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()+
  guides(fill = FALSE, color = FALSE, shape = FALSE)

af_tsne

# Plot t-SNE of presence/absence
pa_tsne<- ggplot(tsne_data_pres_abs, aes(V1, V2, color = province, shape = VOC_436_581)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "In-Sample Allele Presence/Absence",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

pa_tsne

combined_plot_tsne <- plot_grid(af_tsne, pa_tsne, ncol = 2)

ggsave("tsne_regions.png", combined_plot_tsne, width = 16, height = 10, bg = "white")



# POPULATION ALLELE FREQ TSNE

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

# 1) extract allele freqs
#import everything into lists
rds_files <- list.files(path = "FINAL_MOIRE_RESULTS", pattern = "\\MOIRE-RESULTS_FOR_ALLELE_FREQS.RDS$", full.names = TRUE)
rds_files <- rds_files[!rds_files %in% "FINAL_MOIRE_RESULTS/all_samples_complete_filtered_MOIRE-RESULTS_2022_only_FOR_MOI.RDS"] 

# Load each RDS file into the list with the file name as the list name
allele_freqs_list <- list()

for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  allele_freqs_list[[file_name]] <- readRDS(file)
}

# Loop through each element in allele_freqs_list
processed_allele_freq_results <- data.frame()

for (i in seq_along(allele_freqs_list)) {
  
  # Summarize He
  allele_freq_results <- moire::summarize_allele_freqs(allele_freqs_list[[i]])
  allele_freq_results$population <- names(allele_freqs_list[i])
  
  # Add the processed He results to the list
  processed_allele_freq_results <- rbind(processed_allele_freq_results, allele_freq_results)
}

#formatting categories
processed_allele_freq_results$population <- gsub("_2022_MOIRE-RESULTS_FOR_ALLELE_FREQS", "", processed_allele_freq_results$population)

library(stringr)
processed_allele_freq_results <- processed_allele_freq_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"))


# Melt the data frame to convert it from wide to long format
melted <- melt(processed_allele_freq_results, id.vars = c("population", "post_allele_freqs_mean"), measure.vars = "allele")
melted<-melted[,-3]

library(tidyr)

rearranged_processed_allele_freq_results <- melted %>%
  pivot_wider(
    names_from = value,
    values_from = post_allele_freqs_mean
  )

# format
rearranged_processed_allele_freq_results <- as.data.frame(rearranged_processed_allele_freq_results)
rownames(rearranged_processed_allele_freq_results) <- rearranged_processed_allele_freq_results$population
rearranged_processed_allele_freq_results <- rearranged_processed_allele_freq_results[, -1]
rearranged_processed_allele_freq_results <- replace(rearranged_processed_allele_freq_results, is.na(rearranged_processed_allele_freq_results), 0)


# Find rows with "Centre", "North", or "South" in their names
region_columns <- grepl("Centre|North|South", rownames(rearranged_processed_allele_freq_results))
rearranged_processed_allele_freq_results_region <- rearranged_processed_allele_freq_results[region_columns, ]
rearranged_processed_allele_freq_results_province <- rearranged_processed_allele_freq_results[!region_columns, ]

library(stringr)
# Split row names by the last "_"
metadata_province <- rownames(rearranged_processed_allele_freq_results_province)

provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south

metadata_province <- factor(metadata_province, levels = provinces)
metadata_region <- c("North", "South", "Centre", "Centre", "South", "South", "North", "North", "Centre", "Centre", "Centre")

# tsne of provinces
perplexity <- floor((length(metadata_province) - 1) / 3) #highest possible
set.seed(69)
tsne_result_freqs <- Rtsne(as.matrix(rearranged_processed_allele_freq_results_province), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = T, max_iter = 2e4, num_threads = 0, perplexity = perplexity)

# Convert t-SNE results to data frame
tsne_data_freqs <- as.data.frame(tsne_result_freqs$Y)

tsne_data_freqs<- cbind(tsne_data_freqs, metadata_region)

# Plot t-SNE of freqs
pop_allele_freq_tsne <- ggplot(tsne_data_freqs, aes(V1, V2, color = metadata_province, shape = metadata_region)) + # shape = factor(metadata_province$site)
  geom_point(size = 8, alpha = 0.7) +
  labs(title = "",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

pop_allele_freq_tsne

ggsave("pop_allele_freqs_tsne.png", pop_allele_freq_tsne, width = 10, height = 8, bg = "white")


### MANTEL TEST
library(vegan)
library(geosphere)

# Remove rows based on the indices and reorder rows
rows_to_remove <- grepl("Dry", rownames(rearranged_processed_allele_freq_results_province))
rearranged_processed_allele_freq_results_province_nodry <- rearranged_processed_allele_freq_results_province[!rows_to_remove, ]

#reorder rows
desired_order <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy")
indices <- match(desired_order, rownames(rearranged_processed_allele_freq_results_province_nodry))
rearranged_processed_allele_freq_results_province_nodry <- rearranged_processed_allele_freq_results_province_nodry[indices, ]

#calculate distances
bray_curtis_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "bray", diag = T, upper = T)
# manhattan_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "manhattan", diag = T, upper = T)
# euclidean_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "euclidean", diag = T, upper = T)
# gower_dist <- vegdist(rearranged_processed_allele_freq_results_province_nodry, method = "gower", diag = T, upper = T)

# not actual coordinates??
coordinates <- data.frame(
  longitude = c(36.4931, 39.6123, 38.3473, 36.8663, 33.5867, 33.898, 34.8444, 35.3837, 32.5716),
  latitude = c(-12.8779, -11.6455, -15.1056, -16.2828, -16.1742, -19.0834, -19.1548, -23.865, -25.9692)
)

rownames(coordinates) <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy")

#calculate harvesine distance
geo_dist <- distm(coordinates, fun = distHaversine)
rownames(geo_dist) <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy")
colnames(geo_dist) <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Rainy")


#mantel test fucntion to test multiple distance metrics
perform_mantel_test <- function(distance, geo_dist) {
  
  # Perform Mantel test
  mantel_result <- mantel(distance, geo_dist, method = "spearman", permutations = 10000)
  
  # Extract distance vectors
  distance_vec <- as.vector(as.matrix(distance))
  geo_vec <- as.vector(geo_dist)
  
  # Create data frame for plotting
  mat <- data.frame(deest = distance_vec, geo = geo_vec)
  
  # Filter out zero values for Bray-Curtis distance
  mat <- mat[mat$deest > 0, ]
  
  # Create scatter plot
  plot <- ggplot(mat, aes(y = deest, x = geo / 1000)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, aes(fill = geo / 1000)) + 
    geom_smooth(method = "lm", colour = "black", alpha = 0.2) + 
    labs(x = "Harvesine Distance", y = "Bray-Curtis Dissimilarity", fill = "Kilometers") + 
    theme(axis.text.x = element_text(colour = "black", size = 12), 
          axis.text.y = element_text(size = 11, colour = "black"), 
          axis.title = element_text(size = 14, colour = "black"), 
          panel.background = element_blank(), 
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.position = "right",
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11)) +
    scale_fill_continuous(high = "navy", low = "skyblue")
  
  # Return Mantel test result and the scatter plot
  return(list(mantel_result = mantel_result, plot = plot))
}

res <- perform_mantel_test(bray_curtis_dist, geo_dist)
# perform_mantel_test(manhattan_dist, geo_dist)
# perform_mantel_test(euclidean_dist, geo_dist)
# perform_mantel_test(gower_dist, geo_dist)

p <- res$plot
p

ggsave("mantel_pop_allele_Freqs.png", p, width = 8, height = 6, bg = "white")

#######################################################
# 12.- IBD: Proportion of related pairwise infections using IBD between provinces and regions
#######################################################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)

#add VOC
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_doub_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhfr_tr_95_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_436_b")], by = "NIDA2")
combined_df_merged <- merge(combined_df_merged, db[c("NIDA2", "dhps_581_b")], by = "NIDA2")


combined_df_merged$VOC <- ifelse(combined_df_merged$dhps_doub_95_b == 1 & combined_df_merged$dhfr_tr_95_b == 1, "dhps_d_dhfr_tr",
                                 ifelse(combined_df_merged$dhfr_tr_95_b == 1 & combined_df_merged$dhps_doub_95_b == 0, "dhfr_tr",
                                        ifelse(combined_df_merged$dhfr_tr_95_b == 0 & combined_df_merged$dhps_doub_95_b == 1, "dhps_d", "WT")))

#no_genotype
combined_df_merged$VOC <- ifelse(is.na(combined_df_merged$VOC), "no_genotype", combined_df_merged$VOC)

#436 581 genotypes
# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_436
combined_df_merged$dhps_436 <- ifelse(combined_df_merged$dhps_436 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_436 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_436 == 2, "mix", "no_genotype")))

# Replace 0 with WT, 1 with mut, 2 with mix, and NA with no_Genotype in dhps_581
combined_df_merged$dhps_581 <- ifelse(combined_df_merged$dhps_581 == 0, "WT",
                                      ifelse(combined_df_merged$dhps_581 == 1, "mut",
                                             ifelse(combined_df_merged$dhps_581 == 2, "mix", "no_genotype")))

combined_df_merged$VOC_436_581 <- paste0("dhps_436", "-", combined_df_merged$dhps_436, "_", "581", "-", combined_df_merged$dhps_581)

combined_df_merged$VOC_436_581 <- gsub("dhps_436-NA_581-NA", "no_genotype", combined_df_merged$VOC_436_581 )
combined_df_merged$VOC_436_581 <- gsub("dhps_436-WT_581-WT", "WT", combined_df_merged$VOC_436_581 )

unique(combined_df_merged$VOC_436_581)

library(dcifer)
pardef <- par(no.readonly = TRUE)

## 2022 samples ##
#format data
dsmp <- formatDat(combined_df_merged, svar = "NIDA2", lvar = "locus", avar = "pseudo_cigar")
str(dsmp, list.len = 2)

# format metadata
meta <- unique(combined_df_merged[c("NIDA2", "region", "province")])
meta <- meta[match(names(dsmp), meta$NIDA2), ]  # order samples as in dsmp

#estimate naive coi
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
min(coi)

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#order provinces from north to wouth
provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

#calculate ibd
dres0_2022 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                     alpha = 0.05, nr = 1e3)  

# saveRDS(dres0_2022, "dres0_2022.RDS")
dres0_2022 <- readRDS("dres0_2022.RDS")

pdf("dres0_2022_plot.pdf", width = 15, height = 15) 

layout(matrix(1:2, 1), width = c(15, 1))
par(mar = c(1, 1, 2, 1))
alpha <- 0.01         
nsmp  <- length(dsmp)
atsep <- cumsum(nsite)[-length(nsite)]
isig  <- which(dres0_2022[, , "p_value"] <= alpha, arr.ind = TRUE)
dmat  <- dres0_2022[, , "estimate"]
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))] 

plotRel(dmat, isig = isig, draw_diag = TRUE, alpha = alpha, idlab = FALSE, side_id = c(2, 3), srt_id = c(25, 65), lwd_diag = 0.5, border_sig = "darkviolet")

abline(v = atsep, h = atsep, col = "gray45", lty = 5)
atclin <- cumsum(nsite) - nsite/2
mtext(provinces, side = 3, at = atclin, line = 0.2)
mtext(provinces, side = 2, at = atclin, line = 0.2)

par(mar = c(1, 0, 2, 3))
plotColorbar()

dev.off()


### EXAMINE PAIRS OF SAMPLES!
str(dres0_2022)

#extract info
estimates_df <-as.data.frame(dres0_2022[1:924, 1:924, "estimate"])
estimates_p <-as.data.frame(dres0_2022[1:924, 1:924, "p_value"])
estimates_CI_lower <- as.data.frame(dres0_2022[1:924, 1:924, "CI_lower"])
estimates_CI_upper <- as.data.frame(dres0_2022[1:924, 1:924, "CI_upper"])

library(reshape2)
#sort info into a single df

estimates_df_long <- reshape2::melt(as.matrix(estimates_df))
names(estimates_df_long) <- c("sample1", "sample2", "estimate")
estimates_p_long <- reshape2::melt(as.matrix(estimates_p))
names(estimates_p_long) <- c("sample1", "sample2", "p_value")
estimates_CI_lower_long <- reshape2::melt(as.matrix(estimates_CI_lower))
names(estimates_CI_lower_long) <- c("sample1", "sample2", "CI_lower")
estimates_CI_upper_long <- reshape2::melt(as.matrix(estimates_CI_upper))
names(estimates_CI_upper_long) <- c("sample1", "sample2", "CI_upper")

merged_df <- merge(estimates_df_long, estimates_p_long, by = c("sample1", "sample2"))
merged_df <- merge(merged_df, estimates_CI_lower_long, by = c("sample1", "sample2"))
merged_df <- merge(merged_df, estimates_CI_upper_long, by = c("sample1", "sample2"))



#saveRDS(merged_df, "dres0_2022_TABLE.RDS")
merged_df <- readRDS("dres0_2022_TABLE.RDS")

#remove NA rows
merged_df <- merged_df[complete.cases(merged_df),]

#keep significant cases
merged_df_signif <- merged_df[merged_df$p_value < 0.01,]
merged_df_signif <- merged_df[merged_df$estimate > 0.25,]

#merge with geo locations for each pair of samples compared
merged_df_signif_geo <- merge(merged_df_signif, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample1", by.y = "NIDA2")
colnames(merged_df_signif_geo)[7:9] <- c("province_s1", "region_s1", "VOC_s1")
merged_df_signif_geo<- distinct(merged_df_signif_geo)

merged_df_signif_geo <- merge(merged_df_signif_geo, combined_df_merged[, c("NIDA2", "province", "region", "VOC_436_581")], by.x = "sample2", by.y = "NIDA2")
colnames(merged_df_signif_geo)[10:12] <- c("province_s2", "region_s2", "VOC_s2")
merged_df_signif_geo<- distinct(merged_df_signif_geo)

combos_iniciales<- paste0(merged_df_signif[,1], merged_df_signif[,2])
combos_finales<- paste0(merged_df_signif_geo[,1], merged_df_signif_geo[,2])

#sanity check
if (length(combos_finales %in% combos_finales) == dim(merged_df_signif_geo)[1]){
  print("merge was successful.")
} else{
  print("grab a coffee.")
}

#make connectivity columns
merged_df_signif_geo$conn_provinces <- paste0(merged_df_signif_geo$province_s1, "_", merged_df_signif_geo$province_s2)
merged_df_signif_geo$conn_regions <- paste0(merged_df_signif_geo$region_s1, "_", merged_df_signif_geo$region_s2)

# Sort the data frame by estimate in descending order
sorted_df <- merged_df_signif_geo %>% arrange(desc(estimate))


#remove Dry for connectivity comparisons
sorted_df<- sorted_df[!grepl("Dry", sorted_df$conn_provinces), ] #REMOVE TO INCLUDE DRY SEASON IN THE CONNECTIVITY COMPARISON


## PROVINCES CONNECTIVITY

# Calculate the median for each group
median_data <- aggregate(estimate ~ conn_provinces, sorted_df, median)

# Reorder conn_regions based on the median values in descending order
sorted_df$conn_provinces <- factor(sorted_df$conn_provinces, levels = median_data[order(-median_data$estimate), "conn_provinces"])

# Plot with legend and sorted x-axis
prov_conn<- ggplot(sorted_df, aes(x = conn_provinces, y = estimate, fill = conn_regions)) +
  geom_violin(width = 1, aes(color = conn_regions), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = conn_regions), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_fill_discrete(name = "Province") +  # Customize legend title
  ggtitle("Province Connectivity") +
  xlab("")+
  ylab("IBD")+
  guides(color =FALSE)

prov_conn

ggsave("province_connectivity.png", prov_conn, width = 14, height = 7, bg = "white")

## REGIONS CONNECTIVITY

# Calculate the median for each group
median_data <- aggregate(estimate ~ conn_regions, sorted_df, median)

# Reorder conn_regions based on the median values in descending order
sorted_df$conn_regions <- factor(sorted_df$conn_regions, levels = median_data[order(-median_data$estimate), "conn_regions"])

# Plot with legend and sorted x-axis
reg_conn <- ggplot(sorted_df, aes(x = conn_regions, y = estimate, fill = conn_regions)) +
  geom_violin(width = 1, aes(color = conn_regions), alpha = 0.4) +
  geom_boxplot(width = 0.1, aes(color = conn_regions), fill = "white", alpha = 0.4) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11), 
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Region Connectivity") +
  xlab("")+
  ylab("IBD")+
  guides(color = FALSE, fill =FALSE) 

reg_conn

ggsave("region_connectivity.png", reg_conn, width = 12, height = 8, bg = "white")

# pairwise proportions of related infections #

#number of pairwise significantly related (IBD) infections for provinces and regions
table(sorted_df$conn_provinces)
table(sorted_df$conn_regions)


#sample sizes
sample_size_provinces <- combined_df_merged %>%
  group_by(province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

combined_df_merged_nodry <- combined_df_merged %>%
  filter(!grepl("Dry", province))

sample_size_regions <- combined_df_merged_nodry %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions


#calculate percentages for provinces
ibd_samples_provinces <- sorted_df %>%
  group_by(paste0(province_s1, "_", province_s2)) %>%
  summarize(ibd_samples = length(unique(c(sample1, sample2))),
            province_s1 = province_s1,
            province_s2 = province_s2) %>%
  distinct()

ibd_samples_provinces <- ibd_samples_provinces[,-1]

ibd_samples_provinces <- merge(ibd_samples_provinces, sample_size_provinces, by.x = "province_s1", by.y = "province")
colnames(ibd_samples_provinces)[4] <- "province_s1_ss"
  
ibd_samples_provinces <- merge(ibd_samples_provinces, sample_size_provinces, by.x = "province_s2", by.y = "province")
colnames(ibd_samples_provinces)[5] <- "province_s2_ss"

ibd_samples_provinces$total_samples_pairwise <- ifelse(ibd_samples_provinces$province_s1 != ibd_samples_provinces$province_s2, rowSums(ibd_samples_provinces[c("province_s1_ss", "province_s2_ss")]), ibd_samples_provinces$province_s1_ss)

ibd_samples_provinces$perc_ibd_samples_pairwise <- ibd_samples_provinces$ibd_samples / ibd_samples_provinces$total_samples_pairwise

ibd_samples_provinces$pairwise_comparison <- paste0(ibd_samples_provinces$province_s1, "_", ibd_samples_provinces$province_s2)

# Assuming ibd_samples_regions is your data frame
ibd_samples_provinces <- ibd_samples_provinces %>%
  arrange(desc(perc_ibd_samples_pairwise))  # Sort the data frame by perc_ibd_samples_pairwise in descending order

# Reorder the levels of the pairwise_comparison factor
ibd_samples_provinces$pairwise_comparison <- factor(ibd_samples_provinces$pairwise_comparison, 
                                                  levels = ibd_samples_provinces$pairwise_comparison)

ibd_samples_provinces <- merge(ibd_samples_provinces, sorted_df[c("conn_provinces", "conn_regions")], by.x = "pairwise_comparison", by.y = "conn_provinces")

ibd_samples_provinces <- distinct(ibd_samples_provinces)

# Create the bar plot
prop_ibd_prov <- ggplot(ibd_samples_provinces, aes(x = pairwise_comparison, y = perc_ibd_samples_pairwise, fill = conn_regions)) +
  geom_bar(stat = "identity", alpha =0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "proportion of samples") +
  #ggtitle("proportion of significantly related IBD samples")+
  guides(color = FALSE) 

prop_ibd_prov

ggsave("prov_prop_IBD_samples.png", prop_ibd_prov, width = 16, height = 10, bg = "white")

#calculate percentages for regions
ibd_samples_regions <- sorted_df %>%
  group_by(paste0(region_s1, "_", region_s2)) %>%
  summarize(ibd_samples = length(unique(c(sample1, sample2))),
            region_s1 = region_s1,
            region_s2 = region_s2) %>%
  distinct()

ibd_samples_regions <- ibd_samples_regions[,-1]

ibd_samples_regions <- merge(ibd_samples_regions, sample_size_regions, by.x = "region_s1", by.y = "region")
ibd_samples_regions <- ibd_samples_regions[,-4]
colnames(ibd_samples_regions)[4] <- "region_s1_ss"

ibd_samples_regions <- merge(ibd_samples_regions, sample_size_regions, by.x = "region_s2", by.y = "region")
ibd_samples_regions <- ibd_samples_regions[,-5]
colnames(ibd_samples_regions)[5] <- "region_s2_ss"

ibd_samples_regions$total_samples_pairwise <- ifelse(ibd_samples_regions$region_s1 != ibd_samples_regions$region_s2, rowSums(ibd_samples_regions[c("region_s1_ss", "region_s2_ss")]), ibd_samples_regions$region_s1_ss)

ibd_samples_regions$perc_ibd_samples_pairwise <- ibd_samples_regions$ibd_samples / ibd_samples_regions$total_samples_pairwise

ibd_samples_regions$pairwise_comparison <- paste0(ibd_samples_regions$region_s1, "_", ibd_samples_regions$region_s2)

# Assuming ibd_samples_regions is your data frame
ibd_samples_regions <- ibd_samples_regions %>%
  arrange(desc(perc_ibd_samples_pairwise))  # Sort the data frame by perc_ibd_samples_pairwise in descending order

# Reorder the levels of the pairwise_comparison factor
ibd_samples_regions$pairwise_comparison <- factor(ibd_samples_regions$pairwise_comparison, 
                                                  levels = ibd_samples_regions$pairwise_comparison)

# Create the bar plot
prop_ibd_reg <- ggplot(ibd_samples_regions, aes(x = pairwise_comparison, y = perc_ibd_samples_pairwise, fill = pairwise_comparison)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Pairwise Comparison", y = "Percentage of IBD Samples") +
  ggtitle("Percentage of IBD Samples by Pairwise Comparison")+
  guides(color = FALSE, fill =FALSE) 

prop_ibd_reg

ggsave("region_prop_IBD_samples.png", prop_ibd_reg, width = 16, height = 10, bg = "white")


#pairwise IBD between samples with variants of concern and wildtype parasites from the same or different areas.
#remove "no_genotype" category
sorted_df_full_geno <- sorted_df[sorted_df$VOC_s1 != "no_genotype",]
sorted_df_full_geno <- sorted_df_full_geno[sorted_df_full_geno$VOC_s2 != "no_genotype",]

#mean IBD for each pair of genotypes for each pairwise comparison
sorted_df_full_geno_summary <- sorted_df_full_geno %>% 
  group_by(conn_provinces, pairwise_geno = paste0(VOC_s1, "_", VOC_s2)) %>% 
  summarize(IBD = estimate) %>%
  mutate(geno = ifelse(grepl("WT_WT", pairwise_geno), "WT_vs_WT", "WT_vs_res"))

# 
# # Calculate the median IBD for each conn_provinces value
# median_IBD <- aggregate(IBD ~ conn_provinces, data = sorted_df_full_geno_summary[sorted_df_full_geno_summary$geno == "WT_vs_res",], FUN = median)
# 
# # Reorder the levels of conn_provinces based on the median IBD of WT_vs_WT category
# sorted_df_full_geno_summary$conn_provinces <- factor(sorted_df_full_geno_summary$conn_provinces, 
#                                                      levels = median_IBD[order(-median_IBD$IBD), "conn_provinces"])

# Plot with conn_provinces sorted by median IBD of WT_vs_WT category
res_wt <- ggplot(sorted_df_full_geno_summary, aes(x = conn_provinces, y = IBD, fill = geno)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "IBD") +
  ggtitle("")

res_wt

ggsave("province_IBD_res_vs_WT.png", res_wt, width = 24, height = 12, bg = "white")


# sorted_df_full_geno_summary_WT_only<- sorted_df_full_geno_summary %>%
#   filter(grepl("Sofala", conn_provinces))

ggplot(sorted_df_full_geno_summary, aes(x = geno, y = IBD, fill = pairwise_geno)) +
  geom_boxplot() +
  facet_wrap(~ conn_provinces, scales = "free", nrow = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "IBD") +
  ggtitle("IBD by Connection Provinces and Pairwise Genotype")+
  guides(color = FALSE) 

