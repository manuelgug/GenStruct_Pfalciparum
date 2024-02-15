
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

# Objective 3: Describe genetic structure of P. falciparum population in Mozambique in 2021 and 2022 and investigate the origin of parasites with variants of concern

### INDEX

# a) Determine genetic diversity of Plasmodium falciparum in Mozambique
# • Intra-host diversity
#   ◦ Calculate MOI overall and means per province and region for each year (Naïve and eMOI, which incorporates intra-host relatedness between clones) ✔
#   ◦ Calculate % of polyclonal infections per province and region for each year ✔
#   ◦ Fixation index Fws??
# • Population genetic diversity: 
#   ◦ Calculate He per locus per province and region ✔
#   ◦ Correlates with transmission intensity?
# Refs: Brokhattingen et al, preprint; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4 ; Kattenberg et al https://doi.org/10.1128/spectrum.00960-22 

# b) Determine population structure and connectivity between parasite populations in Mozambique
# • PCA to determine population structure
# • Genetic differentiation (Fst) between provinces and regions
# • Proportion of related pairwise infections using IBD between provinces and regions 
# Refs: Arnau; Brokhattingen et al, preprint; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4; Kattenberg et al https://doi.org/10.1128/spectrum.00960-22; Gerlovina et al, Genetics 2022

# c) Determine clonal emergence vs spread from other regions for parasites carrying VOC (dhps581 and dhps436)
# • Map parasites with VOC on PCA (color)
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

#did i found all replicates?
if (all(repeated_nidas_df$NIDA2 %in% merged_df_dups$sample_id)) {
  print("You found all replicates. Proceed with removal")
}else{
  "grab a coffee"
}

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
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region", "run_id")], by="NIDA2", all.y =T) #forcing adding all db nidas

#check for columns with NAs. THESE SAMPLES WERE BAD QUALITY AND THUS FILTERED OUT DURING CONTAMINANTS FILTERING
removed_samples <- combined_df_merged[is.na(combined_df_merged$Category),]
removed_samples$NIDA2

# Remove rows with NIDA2 matching values in removed_samples
combined_df_merged <- combined_df_merged[!combined_df_merged$NIDA2 %in% removed_samples$NIDA2, ]


#sanity check
if (sum(is.na(combined_df_merged)) == 0){
  print("No NAs ✔")
}else{
  print("grab a coffee.")
}

if( sum(!(combined_df_merged$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also the metadata db. No weird samples ✔")
}else{
  print("grab another coffee.")
}

# if (length(unique(db$NIDA2)) == length(unique(combined_df_merged$NIDA2))){
#   print("same amount of samples in metadata db and combined_merged_df")
# }else{
#   print("cocaine maybe?")
# }

#KEEP ONLY DIVERTSITY LOCI!
# are there samples with no Diversity loci sequenced?
unique_categories <- combined_df_merged %>%
  group_by(NIDA2) %>%
  summarize(unique_categories = toString(unique(Category)))

more_removed_sample <- unique_categories[!grepl("Diversity", unique_categories$unique_categories), ] # 1 sample didn't have diversity loci, probably due to a missing pool or something... no problem.
more_removed_sample

combined_df_merged <- combined_df_merged[combined_df_merged$Category == "Diversity",] 


## FURTHER FILTERING

# 1) MAF filering (< 0.01)
combined_df_merged <- combined_df_merged[combined_df_merged$norm.reads.locus  > 0.01, ]

# 2) check for coverage (>50 loci with >threshold reads)
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

#count cells above 50 for each column with the "unique" substring: here, i'm calculating the sample size for each reads count threshold using 50 loci as cutoff: samples with <50 read count per loci below the threshold should be removed
count_above_50 <- function(x) sum(x >= 50, na.rm = TRUE)
unique_columns <- grep("unique", colnames(result_df), value = TRUE)

count_results <- result_df %>%
  summarise_at(vars(unique_columns), count_above_50)

count_results

# decided going with a threshold of >= 50 read depth, 100 is to astringnet for miseq: samples with a coverage of <50 diversity loci with a read depth <50 should be removed.
samples_to_keep <- result_df[result_df$unique_loci_50 >= 50, ]$NIDA2

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


# JUST IN CASE... recalculate n.alleles for each locus of each sample
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
saveRDS(combined_df_merged, "combined_df_merged.RDS")

######################################################################
#----------------------------ANALYZE DATA----------------------------#
######################################################################

#######################################################
# 4.- calculate MOI and eMOI for each run
#######################################################

combined_df_merged <- readRDS("combined_df_merged.RDS") 

colnames(combined_df_merged)[1]<- c("sample_id")

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(combined_df_merged)
burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, .5, length.out = 20)

# run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10)

saveRDS(mcmc_results, "all_samples_complete_filtered_MOIRE-RESULTS.RDS")

#######################################################
# 5.- Present MOI/eMOI results overall and means per province and region for each year
#######################################################

mcmc_results <- readRDS("all__samples_no_further_filtering_MOIRE-RESULTS.RDS") # change name to final file, THIS IS TEST!!

eff_coi <- moire::summarize_effective_coi(mcmc_results)
naive_coi <- moire::summarize_coi(mcmc_results)

# Merge the summaries by sample_id
coi_results <- merge(eff_coi, naive_coi, by = "sample_id")

# label mono and poly infections. NOTE: "proportion of polyclonal infections (eMOI>1.1)" from Nanna's manuscript
coi_results$polyclonal_from_ecoi_med <- ifelse(coi_results$post_effective_coi_med > 1.1, "polyclonal", "monoclonal")

#merge with categorical variables
colnames(coi_results)[1] <- "NIDA2"
coi_results <- merge(coi_results, db[c("NIDA2", "year", "province", "region")], by="NIDA2")

# % polyclonal infections on each province and region per year
polyclonal_percentage_region <- coi_results %>%
  group_by(region, year) %>%
  summarise(polyclonal_percentage_region = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

a <- ggplot(polyclonal_percentage_region, aes(x = region, y = polyclonal_percentage_region, fill = factor(year))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Region", y = "% Polyclonal Infections") +
  facet_wrap(~region, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()
a

polyclonal_percentage_province <- coi_results %>%
  group_by(province, year) %>%
  summarise(polyclonal_percentage_province = mean(polyclonal_from_ecoi_med == "polyclonal") * 100) %>%
  ungroup()

b <- ggplot(polyclonal_percentage_province, aes(x = province, y = polyclonal_percentage_province, fill = factor(year))) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Province", y = "%Polyclonal Infections") +
  facet_wrap(~province, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()
b

# post_effective_coi_med
coi_results$year <- factor(coi_results$year)

c <- ggplot(coi_results, aes(x = region, y = post_effective_coi_med, fill = year)) +
  geom_boxplot() +
  labs(x = "Region", y = "Post Effective COI Median") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
c

d <-ggplot(coi_results, aes(x = province, y = post_effective_coi_med, fill = year)) +
  geom_boxplot() +
  labs(x = "Province", y = "Post Effective COI Median") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
d

# naive coi
e <- ggplot(coi_results, aes(x = region, y = naive_coi, fill = year)) +
  geom_boxplot() +
  labs(x = "Region", y = "Naive COI") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
e

f <-ggplot(coi_results, aes(x = province, y = naive_coi, fill = year)) +
  geom_boxplot() +
  labs(x = "Province", y = "Naive COI") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
f

#######################################################
# 6.- CHECK SAMPLE SIZES FOR EACH PAIR OF VARIABLES: sample size affects He calculation, probably will need rarefactions or something similar
#######################################################

combined_df_merged <- readRDS("combined_df_merged.RDS") 

sample_size_provinces <- combined_df_merged %>%
  group_by(year, province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_provinces

sample_size_regions <- combined_df_merged %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions

################################### 
#ACCUMULATION CURVES (read this https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4885658/; ver genotype_curve de paquete poppr?)

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
accum_curves_2021 <- list()
accum_curves_2022 <- list()

# Get unique years and provinces
unique_provinces <- unique(raref_input$province)

# Iterate over each province
for (province in unique_provinces) {
  
  print(province)
  
  # Subsetting the data for 2021
  sub_2021 <- raref_input[raref_input$year == 2021 & raref_input$province == province, ]
  
  # Check if there are unique NIDA2s for 2021
  if (length(unique(sub_2021$NIDA2)) > 0) {
    # Initialize a list to store unique alleles for each NIDA2 for 2021
    unique_alleles_2021 <- list()
    
    # Iterate over each unique NIDA2 for 2021
    for (nida in unique(sub_2021$NIDA2)) {
      subset_data <- sub_2021[sub_2021$NIDA2 == nida, ]
      unique_alleles_it <- unique(subset_data$allele)
      unique_alleles_2021[[as.character(nida)]] <- unique_alleles_it
    }
    
    # Get unique alleles across all elements for 2021
    all_unique_alleles_2021 <- unique(unlist(unique_alleles_2021))
    
    # Create a matrix to store presence/absence of unique alleles for each element for 2021
    presence_matrix_2021 <- sapply(unique_alleles_2021, function(x) {
      as.integer(all_unique_alleles_2021 %in% x)
    })
    
    # Convert the matrix to a dataframe for 2021
    presence_df_2021 <- as.data.frame(presence_matrix_2021)
    presence_df_2021 <- t(presence_df_2021)
    rownames(presence_df_2021) <- names(unique_alleles_2021)
    colnames(presence_df_2021) <- all_unique_alleles_2021
    
    # CALCULATE CURVE for 2021
    accum_curve_2021 <- specaccum(presence_df_2021, 'random', permutations = 100)
    accum_curves_2021[[province]] <- accum_curve_2021
  }
  
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
colors <- brewer.pal(9, "Paired")

# Plot the curves for 2021
plot(accum_curves_2021[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2021 (per Province)", xlim = c(0,85), ylim = c(0,2000))
for (i in 2:length(accum_curves_2021)) {
  lines(accum_curves_2021[[i]], col = colors[i], lw = 1.5)
}
legend(x = 65, y = 550, legend = names(accum_curves_2021), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022 (per Province)", xlim = c(0,200), ylim = c(0,2200))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 160, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

#conclusion....

# REGION

# Initialize a list to store the rarefaction curves for each year
accum_curves_2021 <- list()
accum_curves_2022 <- list()

# Get unique years and provinces
unique_provinces <- unique(raref_input$region)

# Iterate over each region
for (region in unique_provinces) {
  
  print(region)
  
  # Subsetting the data for 2021
  sub_2021 <- raref_input[raref_input$year == 2021 & raref_input$region == region, ]
  
  # Check if there are unique NIDA2s for 2021
  if (length(unique(sub_2021$NIDA2)) > 0) {
    # Initialize a list to store unique alleles for each NIDA2 for 2021
    unique_alleles_2021 <- list()
    
    # Iterate over each unique NIDA2 for 2021
    for (nida in unique(sub_2021$NIDA2)) {
      subset_data <- sub_2021[sub_2021$NIDA2 == nida, ]
      unique_alleles_it <- unique(subset_data$allele)
      unique_alleles_2021[[as.character(nida)]] <- unique_alleles_it
    }
    
    # Get unique alleles across all elements for 2021
    all_unique_alleles_2021 <- unique(unlist(unique_alleles_2021))
    
    # Create a matrix to store presence/absence of unique alleles for each element for 2021
    presence_matrix_2021 <- sapply(unique_alleles_2021, function(x) {
      as.integer(all_unique_alleles_2021 %in% x)
    })
    
    # Convert the matrix to a dataframe for 2021
    presence_df_2021 <- as.data.frame(presence_matrix_2021)
    presence_df_2021 <- t(presence_df_2021)
    rownames(presence_df_2021) <- names(unique_alleles_2021)
    colnames(presence_df_2021) <- all_unique_alleles_2021
    
    # CALCULATE CURVE for 2021
    accum_curve_2021 <- specaccum(presence_df_2021, 'random', permutations = 100)
    accum_curves_2021[[region]] <- accum_curve_2021
  }
  
  # Subsetting the data for 2022
  sub_2022 <- raref_input[raref_input$year == 2022 & raref_input$region == region, ]
  
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

# Plot the curves for 2021
plot(accum_curves_2021[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2021 (per Region)", xlim = c(0,150), ylim = c(0,2000))
for (i in 2:length(accum_curves_2021)) {
  lines(accum_curves_2021[[i]], col = colors[i], lw = 1.5)
}
legend(x = 120, y = 550, legend = names(accum_curves_2021), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022 (per Region)", xlim = c(0,500), ylim = c(0,3000))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 400, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.5)



########################
###### PCA ########
########################

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
pca_labels<- combined_df_merged %>% distinct(NIDA2, year, province, region)

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
 
#TSNE
set.seed(420)
tsne_result_freqs <- Rtsne(as.matrix(rearranged_filtered), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = T, max_iter = 2e4, num_threads = 0)
set.seed(420)
tsne_result_pres_abs <- Rtsne(as.matrix(rearranged_pres_abs), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = T, max_iter = 2e4, num_threads = 0)

# Convert t-SNE results to data frame
tsne_data_freqs <- as.data.frame(tsne_result_freqs$Y)
tsne_data_pres_abs <- as.data.frame(tsne_result_pres_abs$Y)

# Plot t-SNE of freqs
ggplot(tsne_data_freqs, aes(V1, V2, color = factor(pca_labels$region), shape = factor(pca_labels$year))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "t-SNE of Genetic Content (allele frequency)",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# Plot t-SNE of presence/absence
ggplot(tsne_data_pres_abs, aes(V1, V2, color = factor(pca_labels$region), shape = factor(pca_labels$year))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "t-SNE of Genetic Content (presence/absence of alleles)",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

#######################################################
# 7.- Fst: Genetic differentiation (Fst) between provinces and regions
#######################################################


#######################################################
# 8.- IBD: Proportion of related pairwise infections using IBD between provinces and regions
#######################################################

combined_df_merged <- readRDS("combined_df_merged.RDS")

library(dcifer)
pardef <- par(no.readonly = TRUE)

#subset 2021 data
combined_df_merged_2021 <- combined_df_merged[combined_df_merged$year == "2021", ]
#subset 2022 data
combined_df_merged_2022 <- combined_df_merged[combined_df_merged$year == "2022", ]


## 2021 samples ##
#format data
dsmp <- formatDat(combined_df_merged_2021, svar = "NIDA2", lvar = "locus", avar = "pseudo_cigar")
str(dsmp, list.len = 2)

# format metadata
meta <- unique(combined_df_merged_2021[c("NIDA2", "region", "province")])
meta <- meta[match(names(dsmp), meta$NIDA2), ]  # order samples as in dsmp

#estimate naive coi
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
min(coi)

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#order provinces from north to south
provinces <- c("Niassa", "Zambezia", "Nampula", "Manica", "Inhambane", "Maputo") #ordered from north to south
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

#calculate ibd
dres0_2021 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                 alpha = 0.05, nr = 1e3)  

# saveRDS(dres0_2021, "dres0_2021.RDS")

pdf("dres0_2021_plot.pdf", width = 15, height = 15) 

layout(matrix(1:2, 1), width = c(15, 1))
alpha <- 0.05                        
par(mar = c(1, 1, 2, 1))
nsmp  <- length(dsmp)
atsep <- cumsum(nsite)[-length(nsite)]
isig  <- which(dres0_2021[, , "p_value"] <= alpha, arr.ind = TRUE)
dmat  <- dres0_2021[, , "estimate"]
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))] 

plotRel(dmat, isig = isig, draw_diag = TRUE, alpha = alpha, idlab = FALSE, side_id = c(2, 3), srt_id = c(25, 65), lwd_diag = 0.5, border_sig = "darkviolet")

abline(v = atsep, h = atsep, col = "gray45", lty = 5)
atclin <- cumsum(nsite) - nsite/2
mtext(provinces, side = 3, at = atclin, line = 0.2)
mtext(provinces, side = 2, at = atclin, line = 0.2)

par(mar = c(1, 0, 2, 3))
plotColorbar()

dev.off()


## 2022 samples ##
#format data
dsmp <- formatDat(combined_df_merged_2022, svar = "NIDA2", lvar = "locus", avar = "pseudo_cigar")
str(dsmp, list.len = 2)

# format metadata
meta <- unique(combined_df_merged_2022[c("NIDA2", "region", "province")])
meta <- meta[match(names(dsmp), meta$NIDA2), ]  # order samples as in dsmp

#estimate naive coi
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
min(coi)

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#order provinces from north to wouth
provinces <- c("Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica", "Sofala", "Inhambane", "Maputo") #ordered from north to south
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

#calculate ibd
dres0_2022 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                     alpha = 0.05, nr = 1e3)  

# saveRDS(dres0_2022, "dres0_2022.RDS")

pdf("dres0_2022_plot.pdf", width = 15, height = 15) 

layout(matrix(1:2, 1), width = c(15, 1))
par(mar = c(1, 1, 2, 1))
alpha <- 0.05         
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


#######################################################
# 7.- calculate He for each population (per year per region/province)
#######################################################

#subset 2021 data
combined_df_merged_2021 <- combined_df_merged[combined_df_merged$year == "2021", ]
#subset 2022 data
combined_df_merged_2022 <- combined_df_merged[combined_df_merged$year == "2022", ]

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
  saveRDS(mcmc_results, paste0(output_name, "_MOIRE-RESULTS.RDS"))
}

# Create a list of data frames and corresponding years
data_frames <- list(combined_df_merged_2021, combined_df_merged_2022)
years <- list("2021", "2022")

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
# 8.- He and Fws results
#######################################################

combined_df_merged <- readRDS("combined_df_merged.RDS")

# 1) calculate heterozygosity of the population (He); pop = province, region
#import everything into lists
rds_files <- list.files(pattern = "\\MOIRE-RESULTS.RDS$", full.names = TRUE)
rds_files <- rds_files[!rds_files %in% "./all__samples_no_further_filtering_MOIRE-RESULTS.RDS"] #check name for later, may not even be needed.

He_results_list <- list()

# Load each RDS file into the list with the file name as the list name
for (file in rds_files) {
  print(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  He_results_list[[file_name]] <- readRDS(file)
}

processed_He_results <- data.frame()

# Loop through each element in He_results_list
for (i in seq_along(He_results_list)) {
  
  # Summarize He
  He_results <- moire::summarize_he(He_results_list[[i]])
  He_results$population <- names(He_results_list[i])
  
  # Add the processed coi_results to the list
  processed_He_results <- rbind(processed_He_results, He_results)
}

#formatting categories
processed_He_results$population <- gsub("_MOIRE-RESULTS", "", processed_He_results$population)

library(stringr)
processed_He_results <- processed_He_results %>%
  mutate(geo = ifelse(str_detect(population, "North|South|Center"), "region", "province"),
         year = substr(population, nchar(population) - 3, nchar(population)))

processed_He_results$population <- gsub("TEST_|_202.*", "", processed_He_results$population)
processed_He_results$year <- as.numeric(processed_He_results$year)


# 2) calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
heterozygosity_data <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))


# 3) calculate 1-Fws: 1 - Fws = Hw/He
#add processed_He_results$post_stat_mean co heterozygosity_data$He when processed_He_results$population, processed_He_results$year and processed_He_results$lous match heterozygosity_data$province, heterozygosity_data$year and heterozygosity_data$locus, respectively

#merge He from provinces
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus", "year" = "year")) %>%
  filter(population == province) %>%
  mutate(He_province = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, year, He_province)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus", "year"))

heterozygosity_data <- distinct(heterozygosity_data)
heterozygosity_data

#merge He from regions
merged_data <- heterozygosity_data %>%
  left_join(processed_He_results, by = c("locus" = "locus", "year" = "year")) %>%
  filter(population == region) %>%
  mutate(He_region = ifelse(is.na(post_stat_mean), NA, post_stat_mean)) %>%
  select(NIDA2, locus, year, He_region)

heterozygosity_data <- heterozygosity_data %>%
  left_join(merged_data, by = c("NIDA2", "locus", "year"))

heterozygosity_data <- distinct(heterozygosity_data)

#calculate 1-Fws for province and region as populations
heterozygosity_data$fws_province <- heterozygosity_data$Hw/heterozygosity_data$He_province
heterozygosity_data$fws_region <- heterozygosity_data$Hw/heterozygosity_data$He_region

# Columns to keep
columns_to_keep <- c("NIDA2", "locus", "year", "province", "region", "run_id",
                     "Hw", "He_province", "He_region", "fws_province", "fws_region")
# Filter columns
heterozygosity_data_filtered <- heterozygosity_data %>%
  select(all_of(columns_to_keep))

# Keep unique rows
heterozygosity_data_filtered <- distinct(heterozygosity_data_filtered)

#sanity check
if ((length(unique(heterozygosity_data_filtered$NIDA2)) == length(unique(combined_df_merged$NIDA2))) & 
    (length(unique(heterozygosity_data_filtered$locus)) == length(unique(combined_df_merged$locus))) & 
    (length(unique(heterozygosity_data_filtered$year)) == length(unique(combined_df_merged$year))) & 
    (length(unique(heterozygosity_data_filtered$province)) == length(unique(combined_df_merged$province))) & 
    (length(unique(heterozygosity_data_filtered$region)) == length(unique(combined_df_merged$region))) & 
    (length(unique(heterozygosity_data_filtered$run_id)) == length(unique(combined_df_merged$run_id)))) {
  print("All looks good.")
}



# linear regression and correlation coeff of He vs eMOI
