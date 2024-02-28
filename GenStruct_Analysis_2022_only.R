
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
#   ◦ Fixation index Fws?? ✔
# • Population genetic diversity: 
#   ◦ Calculate He per locus per province and region ✔
#   ◦ Correlates with transmission intensity?
# Refs: Brokhattingen et al, preprint; Fola et al, Nature https://doi.org/10.1038/s41564-023-01461-4 ; Kattenberg et al https://doi.org/10.1128/spectrum.00960-22 

# b) Determine population structure and connectivity between parasite populations in Mozambique
# • PCA to determine population structure ✔
# • Genetic differentiation (Fst) between provinces and regions ✔
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

#count cells above 100 for each column with the "unique" substring: here, i'm calculating the sample size for each reads count threshold using 50 loci as cutoff: samples with <50 read count per loci below the threshold should be removed
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

sample_size_regions <- combined_df_merged %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions


######################################################################
#----------------------------ANALYZE DATA----------------------------#
######################################################################

#######################################################
# 5.- calculate MOI and eMOI for each run
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
# 5.- Present MOI/eMOI results overall and means per province and region for each year
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
coi_results_region$region <- factor(coi_results_region$region, levels = regions)
polyclonal_percentage_region$region <- factor(polyclonal_percentage_region$region, levels = regions)
polyclonal_percentage_province$province <- factor(polyclonal_percentage_province$province, levels = provinces)

a <- ggplot(coi_results, aes(x = naive_coi, fill = province)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_wrap(~ province , scales = "fixed", nrow = 1) + 
  labs(title = "",
       x = "Naive COI",
       y = "Frequency",
       fill = "Province") +
  theme_minimal() +
  guides(fill = FALSE) 

a
#ggsave("naive_coi_provinces_ditros.png", a, width = 14, height = 6, bg = "white")


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

c <- ggplot(coi_results, aes(x = post_effective_coi_med, fill = province)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.7) +
  facet_wrap(~ province , scales = "fixed", nrow = 1) + 
  labs(title = "",
       x = "Post Effective COI Median",
       y = "Frequency",
       fill = "Province") +
  theme_minimal() +
  guides(fill = FALSE) 

c

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


e <- ggplot(polyclonal_percentage_region, aes(x = region, y = polyclonal_percentage_region, fill = region)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "% Polyclonal Infections") +
  #facet_wrap(~region, scales = "fixed", ncol = 3) +
  #scale_fill_manual(values = c("2022" = "orange")) + 
  theme_minimal()+
  guides(fill = FALSE) 
e


f <- ggplot(polyclonal_percentage_province, aes(x = province, y = polyclonal_percentage_province,  fill = province))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "%Polyclonal Infections") +
  #facet_wrap(~province, scales = "fixed", nrow = 1) +
  #scale_fill_manual(values = c("2022" = "orange")) + 
  theme_minimal()+
  guides(fill = FALSE) 
f
#ggsave("perc_polyclonal_provinces.png", f, width = 14, height = 10, bg = "white")


################################### 
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

#pdf("accumulation_curves_provinces.pdf", width = 12, height = 8)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Allele Accumulation Curves per Province", xlim = c(0,130), ylim = c(0,1800))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 110, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.7)

#dev.off()
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

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022 (per Region)", xlim = c(0,350), ylim = c(0,2700))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend(x = 310, y = 950, legend = names(accum_curves_2022), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

#######################################################
# 9.- calculate He for each population (per region/province)
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
# 10.- He and Fws results 
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

ggplot(mean_Fws_per_individual, aes(x = province, y = mean_indiv_fws_province, fill = province)) +
  geom_violin() +  # Add violin plot
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.4) +  # Add swarm plot
  labs(x = "Province", y = "Mean Fws per Individual") +
  theme_minimal()+
  guides(fill = FALSE)

ggplot(mean_Fws_per_individual, aes(x = region, y = mean_indiv_fws_region, fill = region)) +
  geom_violin() +  # Add violin plot
  geom_point(position = position_jitter(width = 0.2), size = 2, alpha = 0.4) +  # Add swarm plot
  labs(x = "Region", y = "Mean Fws per Individual") +
  theme_minimal()+
  guides(fill = FALSE)


# STATISTICAL ANALYSES:::     CHECK CAREFULLY!!!!!
# Filter data by province and region
combined_data_province <- rbind(data.frame(He_province = processed_He_results[processed_He_results$geo == "province", ]$post_stat_mean, 
                                           province = processed_He_results[processed_He_results$geo == "province", ]$population))

combined_data_region <- rbind(data.frame(He_region = processed_He_results[processed_He_results$geo == "region", ]$post_stat_mean, 
                                         region = processed_He_results[processed_He_results$geo == "region", ]$population))


# CHECK FOR NORMALITY IN THE He DATA
# PROVINCE
# Create Q-Q plots
ggplot(combined_data_province, aes(sample = He_province)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Q-Q Plot of He_province by Province",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  facet_grid(~ province, scales = "free")

# Perform Shapiro-Wilk
shapiro_test_results_province <- combined_data_province %>%
  group_by(province) %>%
  summarize(
    p_value = shapiro.test(He_province)$p.value
  )

shapiro_test_results_province

#REGION
# Create Q-Q plots
ggplot(combined_data_region, aes(sample = He_region)) +
  geom_qq() +
  geom_qq_line() +
  labs(title = "Q-Q Plot of He_region by Province and Year",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  facet_grid(~region , scales = "free")

# Perform Shapiro-Wilk test
shapiro_test_results_region <- combined_data_region %>%
  group_by(region) %>%
  summarize(
    p_value = shapiro.test(He_region)$p.value
  )

shapiro_test_results_region


## He DATA IS NOT NORMAL, SO KRUSKAL-WALLIS (wilcox): 
# 1) SPATIAL COMPARISNS
pairwise_province_He_2022 <- pairwise.wilcox.test(combined_data_province$He_province, 
                                                  combined_data_province$province, p.adjust.method = "bonferroni")

p.val_province_He_2022 <- melt(pairwise_province_He_2022[[3]])
signif_p.val_province_He_2022 <- p.val_province_He_2022[p.val_province_He_2022$value <0.05 & !is.na(p.val_province_He_2022$value),]

pairwise_region_He_2022 <- pairwise.wilcox.test(combined_data_region$He_region, 
                                                combined_data_region$region, p.adjust.method = "bonferroni")

p.val_region_He_2022 <- melt(pairwise_region_He_2022[[3]])
signif_p.val_region_He_2022 <- p.val_region_He_2022[p.val_region_He_2022$value <0.05 & !is.na(p.val_region_He_2022$value),]


# Changes in He distributions and means by year

provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")

combined_data_province$province <- factor(combined_data_province$province, levels = provinces)
combined_data_region$region <- factor(combined_data_region$region, levels = regions)

#provinces
mean_data <- combined_data_province %>%
  group_by(province) %>%
  summarize(mean_He = mean(He_province, na.rm = TRUE),
            median_He = median(He_province, na.rm = TRUE))

ggplot(combined_data_province, aes(x = He_province)) +
  geom_histogram(alpha = 0.7, bins = 30) +
  geom_vline(data = mean_data, aes(xintercept = mean_He), linetype = "solid", color = "limegreen") +
  geom_vline(data = mean_data, aes(xintercept = median_He), linetype = "solid", color = "orange") +
  labs(title = "Province Heterozygosity",
       x = "He",
       y = "Frequency") +
  facet_wrap(~ province) +
  theme_minimal()


ggplot(combined_data_province, aes(x = province, y = He_province, fill = province)) +
  geom_violin(alpha =0.7) +  # Add violin plot
  ggbeeswarm::geom_beeswarm() +
  labs(x = "Region", y = "Mean Fws per Individual") +
  theme_minimal()+
  guides(fill = FALSE)


#ggsave("mean_he_province_distros.png", ph, width = 14, height = 10, bg = "white")

#regions
mean_data_region <- combined_data_region %>%
  group_by(region) %>%
  summarize(mean_He = mean(He_region, na.rm = TRUE),
            median_He = median(He_region, na.rm = TRUE))

ggplot(combined_data_region, aes(x = He_region)) +
  geom_histogram(alpha = 0.7, bins = 30) +
  geom_vline(data = mean_data_region, aes(xintercept = mean_He), linetype = "solid", color = "limegreen") +
  geom_vline(data = mean_data_region, aes(xintercept = median_He), linetype = "solid", color = "orange") +
  labs(title = "Region Heterozygosity",
       x = "He",
       y = "Frequency") +
  facet_wrap(~ region) +
  theme_minimal()

ggplot(combined_data_region, aes(x = region, y = He_region, fill = region)) +
  geom_violin() +  # Add violin plot
  ggbeeswarm::geom_beeswarm(dodge.width = 0.5) +
  labs(x = "Region", y = "Mean Fws per Individual") +
  theme_minimal()+
  guides(fill = FALSE)


########################
###### PCA ########
########################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


#input for multivariate analyses
raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
                                   year = combined_df_merged$year, 
                                   province = combined_df_merged$province,
                                   region = combined_df_merged$region,
                                   locus = combined_df_merged$locus,
                                   n.alleles = combined_df_merged$n.alleles,
                                   norm.reads.locus = combined_df_merged$norm.reads.locus,
                                   allele = paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar),
                                   run_id = combined_df_merged$run_id))

# ## PICK LOCI WITH HIGHER He for multivariate analyses::: MAKES WORSE PLOTS
# 
# # 1) calculate heterozygosity of the population (He); pop = province, region
# #import everything into lists
# rds_files <- list.files(pattern = "\\MOIRE-RESULTS.RDS$", full.names = TRUE)
# rds_files <- rds_files[!rds_files %in% "./all__samples_no_further_filtering_MOIRE-RESULTS.RDS"] #check name for later, may not even be needed.
# 
# 
# # Load each RDS file into the list with the file name as the list name
# He_results_list <- list()
# 
# for (file in rds_files) {
#   print(file)
#   file_name <- tools::file_path_sans_ext(basename(file))
#   He_results_list[[file_name]] <- readRDS(file)
# }
# 
# # Loop through each element in He_results_list
# processed_He_results <- data.frame()
# 
# for (i in seq_along(He_results_list)) {
#   
#   # Summarize He
#   He_results <- moire::summarize_he(He_results_list[[i]])
#   He_results$population <- names(He_results_list[i])
#   
#   # Add the processed He results to the list
#   processed_He_results <- rbind(processed_He_results, He_results)
# }
# 
# #formatting categories
# processed_He_results$population <- gsub("_MOIRE-RESULTS", "", processed_He_results$population)
# 
# library(stringr)
# processed_He_results <- processed_He_results %>%
#   mutate(geo = ifelse(str_detect(population, "North|South|Centre"), "region", "province"),
#          year = substr(population, nchar(population) - 3, nchar(population)))
# 
# processed_He_results$population <- gsub("TEST_|_202.*", "", processed_He_results$population)
# processed_He_results$year <- as.numeric(processed_He_results$year)
# 
# ## pick and sort loci by He
# loci_He <- processed_He_results %>% 
#   group_by(locus) %>%
#   summarize(mean_He = mean(post_stat_mean))%>%
#   arrange(desc(mean_He))
# 
# #keep loci He > 0.7 
# loci_to_keep_He <- loci_He[loci_He$mean_He >= 0.0,]$locus
# raref_input <- raref_input[raref_input$locus %in% loci_to_keep_He, ]


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
 

provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
regions <- c("North", "Centre", "South")

#PCA
pca_result <- prcomp(rearranged, scale. = F)

# Extract the principal component scores
pcs <- as.data.frame(pca_result$x)

# Combine principal component scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# Calculate percentage variance explained by each principal component
pc_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create the plot
ggplot(pcs_with_labels, aes(x = PC1, y = PC2, color = factor(pcs_with_labels$province), shape = factor(pcs_with_labels$region))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "",
       x = paste0("PC1: ", round(pc_variance[1], 2), "%\n"),
       y = paste0("PC2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal()


# Perform PCA prsence/absence
pca_result <- prcomp(rearranged_pres_abs, scale. = F)

# Extract the principal component scores
pcs <- as.data.frame(pca_result$x)

# Combine principal component scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# Calculate percentage variance explained by each principal component
pc_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Create the plot
ggplot(pcs_with_labels, aes(x = PC1, y = PC2, color = factor(pcs_with_labels$province), shape = factor(pcs_with_labels$region))) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "",
       x = paste0("PC1: ", round(pc_variance[1], 2), "%\n"),
       y = paste0("PC2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal()


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
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
# pc_variance <- pcoa_result$values / sum(pcoa_result$values) * 100

ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = region)) +
  geom_point(size = 4, alpha = 0.7) +
  # labs(title = "",
  #      x = paste0("Axis 1: ", round(pc_variance[1], 2), "%\n"),
  #      y = paste0("Axis 2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal()

#ggsave("pcoa_sample_freqs.png", pc, width = 12, height = 10, bg = "white")

#PCoA presence/absence

# Compute Bray-Curtis dissimilarity matrix
bray_curtis_dist <- vegdist(rearranged_pres_abs, method = "bray")

# Perform PCoA
pcoa_result <- pcoa(bray_curtis_dist)

# Extract the principal coordinate scores
pcs <- as.data.frame(pcoa_result$vectors)

# Combine principal coordinate scores with region labels
pcs_with_labels <- cbind(pcs, province = pca_labels$province, region = pca_labels$region)

pcs_with_labels$province <- factor(pcs_with_labels$province, levels = provinces)
pcs_with_labels$region <- factor(pcs_with_labels$region, levels = regions)

# # Plot PCoA
# pc_variance <- pcoa_result$values / sum(pcoa_result$values) * 100

ggplot(pcs_with_labels, aes(x = Axis.1, y = Axis.2, color = province, shape = region)) +
  geom_point(size = 4, alpha = 0.7) +
  # labs(title = "",
  #      x = paste0("Axis 1: ", round(pc_variance[1], 2), "%\n"),
  #      y = paste0("Axis 2: ", round(pc_variance[2], 2), "%")) +
  theme_minimal()


#TSNE
set.seed(69)
perplexity <- floor((nrow(rearranged_filtered) - 1) / 3) #highest possible, if needed
tsne_result_freqs <- Rtsne(as.matrix(rearranged_filtered), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = F, max_iter = 1e4, num_threads = 0, perplexity = 100)
set.seed(69)
tsne_result_pres_abs <- Rtsne(as.matrix(rearranged_pres_abs), dims = 2, verbose = TRUE, check_duplicates = FALSE, pca_center = F, max_iter = 1e4, num_threads = 0, perplexity = 100)

# Convert t-SNE results to data frame
tsne_data_freqs <- as.data.frame(tsne_result_freqs$Y)
tsne_data_freqs <- cbind(tsne_data_freqs, province = pca_labels$province, region = pca_labels$region)
tsne_data_pres_abs <- as.data.frame(tsne_result_pres_abs$Y)
tsne_data_pres_abs <- cbind(tsne_data_pres_abs, province = pca_labels$province, region = pca_labels$region)

# Order factors
tsne_data_freqs$province <- factor(tsne_data_freqs$province, levels = provinces)
tsne_data_freqs$region <- factor(tsne_data_freqs$region, levels = regions)
tsne_data_pres_abs$province <- factor(tsne_data_pres_abs$province, levels = provinces)
tsne_data_pres_abs$region <- factor(tsne_data_pres_abs$region, levels = regions)

# Plot t-SNE of freqs
ggplot(tsne_data_freqs, aes(V1, V2, color = province, shape = region)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "t-SNE of Genetic Content (allele frequency)",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

#ggsave("tsne_sample_freqs.png", ts, width = 12, height = 10, bg = "white")

# Plot t-SNE of presence/absence
ggplot(tsne_data_pres_abs, aes(V1, V2, color = province, shape = region)) +
  geom_point(size = 4, alpha = 0.7) +
  labs(title = "t-SNE of Genetic Content (presence/absence of alleles)",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()


#######################################################
# 7.- IBD: Proportion of related pairwise infections using IBD between provinces and regions
#######################################################

combined_df_merged <- readRDS("FINAL_MOIRE_RESULTS/combined_df_merged_2022_only.RDS")
combined_df_merged$province <- gsub(" ", "_", combined_df_merged$province) # for cabo delgado
#rename alleles
combined_df_merged$allele <- paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)


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
provinces <- c("Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
nsite     <- table(meta$province)[provinces]
ord       <- order(factor(meta$province, levels = provinces))
dsmp <- dsmp[ord]
coi  <- coi[ ord]

#calculate ibd
dres0_2022 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                     alpha = 0.05, nr = 1e3)  

#saveRDS(dres0_2022, "dres0_2022.RDS")
# dres0_2022 <- readRDS("dres0_2022.RDS")

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


##############################
# ALLELE FREQ TSNE
#############################

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
ggplot(tsne_data_freqs, aes(V1, V2, color = metadata_province, shape = metadata_region)) + # shape = factor(metadata_province$site)
  geom_point(size = 10, alpha = 0.7) +
  labs(title = "t-SNE of Allele Frequencies",
       x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()



## TO DO
#1) STATISTICS, check nanna's paper (ongoing)
#2) what to do when pop He is super low and Hw is 0.5 or above... 1-fws ends up HUGE
#3) Fst
#4) subsampling


#######################################3
# 11.- pairwise FST  (https://biology.stackexchange.com/questions/40756/calculating-pairwise-fst-from-allele-frequencies)
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

sample_size_regions <- combined_df_merged %>%
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
# 2.- calculate Fst [(Ht - Hs) / Ht] resampling loci to get confidence intervals (95%; 2.5 and 97.5 percentiles)
# 3.- plot mean genome-wide Fst with mean confidence intervals

# Ideas: use only loci with highest He to make the CI lower?


#1) FOR REGIONS
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


library(boot)

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

# Create a list to store bootstrap results
bootstrap_results <- list()

# Loop through each pairwise comparison
for (i in 1:nrow(combinations)) {
  pop1 <- as.character(combinations$pop1[i])
  pop2 <- as.character(combinations$pop2[i])
  
  # Subset data for the current pairwise comparison
  pairwise_data <- fst_results_df[fst_results_df$pop1 == pop1 & fst_results_df$pop2 == pop2, ]
  
  # Define a unique name for the pairwise comparison
  comparison_name <- paste(pop1, pop2, sep = "_vs_")
  
  # Perform bootstrap resampling
  bootstrap_results[[comparison_name]] <- boot(pairwise_data, calculate_FST, R = 10000)
}


# Calculate the mean of bootstrap replicates for each pairwise comparison
mean_FST <- sapply(bootstrap_results, function(result) mean(result$t0))

# Create a function to extract confidence intervals for the mean
get_mean_confidence_intervals <- function(bootstrap_result) {
  ci <- boot.ci(bootstrap_result, type = "perc")
  lower_limit <- ci$percent[4]
  upper_limit <- ci$percent[5]
  return(c(lower_limit, upper_limit))
}

# Extract confidence intervals for the mean of FST
mean_confidence_intervals <- lapply(bootstrap_results, get_mean_confidence_intervals)

# Create a data frame to store the results
mean_FST_df <- data.frame(
  pairwise_comparisons = names(mean_FST),
  mean_FST = ifelse(is.na(mean_FST), 0, mean_FST),
  lower_limit = sapply(mean_confidence_intervals, `[`, 1),
  upper_limit = sapply(mean_confidence_intervals, `[`, 2)
)

# Remove "_2022" from population names and reorder according to the specified order
#mean_FST_df$pairwise_comparisons <- gsub("_2022", "", mean_FST_df$pairwise_comparisons)

# Exclude rows with 0 in the last three columns
mean_FST_df_filtered <- subset(mean_FST_df, lower_limit != 0 | upper_limit != 0 | mean_FST != 0)

#exclude repeated comparisons: CAREFUL!! IF MEANS FOR DIFFERENT COMPARISONS ARE THE SAME IT WILL ELIMINATE ONE. not the best method.
mean_FST_df_filtered<- mean_FST_df_filtered %>%
  distinct(mean_FST, .keep_all = TRUE)

# Reorder pairwise_comparisons based on mean_FST
mean_FST_df_filtered$pairwise_comparisons <- factor(mean_FST_df_filtered$pairwise_comparisons, 
                                                    levels = mean_FST_df_filtered$pairwise_comparisons[order(mean_FST_df_filtered$mean_FST)])

# Create the boxplot
ggplot(mean_FST_df_filtered, aes(x = pairwise_comparisons, y = mean_FST)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = 0.2) +
  labs(x = "Pairwise Comparisons", y = "Average genome-wide Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()+
  coord_flip()  # Flip the coordinates to make the plot horizontal

#ggsave("fst_CI_provinces.png", fstci, width = 12, height = 10, bg = "white")

# Create heatmap for regions comparisons
create_heatmap <- function(data, title) {
  
  data$mean_FST<- round(data$mean_FST, 3)
  
  ggplot(data, aes(x = pop2, y = pop1, fill = mean_FST, label = round(mean_FST, 3))) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(data$mean_FST), max(data$mean_FST))) +  # Adjust scale limits
    theme_minimal()+
    labs(x = "", y = "")
}

# Split pairwise_comparisons into pop1 and pop2
library(tidyr)
mean_FST_df <- separate(mean_FST_df, pairwise_comparisons, into = c("pop1", "pop2"), sep = "_vs_")

# Reorder pop1 and pop2 columns based on regions order
regions <- c("North", "Centre", "South") #ordered from north to south
regions <- rev(regions)

mean_FST_df$pop1 <- factor(mean_FST_df$pop1, levels = regions)
mean_FST_df$pop2 <- factor(mean_FST_df$pop2, levels = regions)

heatmap_2022_regions <- create_heatmap(mean_FST_df)
print(heatmap_2022_regions)


#2) FOR PROVINCES
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


library(boot)

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

# Create a list to store bootstrap results
bootstrap_results <- list()

# Loop through each pairwise comparison
for (i in 1:nrow(combinations)) {
  pop1 <- as.character(combinations$pop1[i])
  pop2 <- as.character(combinations$pop2[i])
  
  # Subset data for the current pairwise comparison
  pairwise_data <- fst_results_df[fst_results_df$pop1 == pop1 & fst_results_df$pop2 == pop2, ]
  
  # Define a unique name for the pairwise comparison
  comparison_name <- paste(pop1, pop2, sep = "_vs_")
  
  # Perform bootstrap resampling
  bootstrap_results[[comparison_name]] <- boot(pairwise_data, calculate_FST, R = 10000)
}


# Calculate the mean of bootstrap replicates for each pairwise comparison
mean_FST <- sapply(bootstrap_results, function(result) mean(result$t0))

# Create a function to extract confidence intervals for the mean
get_mean_confidence_intervals <- function(bootstrap_result) {
  ci <- boot.ci(bootstrap_result, type = "perc")
  lower_limit <- ci$percent[4]
  upper_limit <- ci$percent[5]
  return(c(lower_limit, upper_limit))
}

# Extract confidence intervals for the mean of FST
mean_confidence_intervals <- lapply(bootstrap_results, get_mean_confidence_intervals)

# Create a data frame to store the results
mean_FST_df <- data.frame(
  pairwise_comparisons = names(mean_FST),
  mean_FST = ifelse(is.na(mean_FST), 0, mean_FST),
  lower_limit = sapply(mean_confidence_intervals, `[`, 1),
  upper_limit = sapply(mean_confidence_intervals, `[`, 2)
)

# Remove "_2022" from population names and reorder according to the specified order
#mean_FST_df$pairwise_comparisons <- gsub("_2022", "", mean_FST_df$pairwise_comparisons)

# Exclude rows with 0 in the last three columns
mean_FST_df_filtered <- subset(mean_FST_df, lower_limit != 0 | upper_limit != 0 | mean_FST != 0)

#exclude repeated comparisons: CAREFUL!! IF MEANS FOR DIFFERENT COMPARISONS ARE THE SAME IT WILL ELIMINATE ONE. not the best method.
mean_FST_df_filtered<- mean_FST_df_filtered %>%
  distinct(mean_FST, .keep_all = TRUE)

# Reorder pairwise_comparisons based on mean_FST
mean_FST_df_filtered$pairwise_comparisons <- factor(mean_FST_df_filtered$pairwise_comparisons, 
                                                    levels = mean_FST_df_filtered$pairwise_comparisons[order(mean_FST_df_filtered$mean_FST)])

# Create the boxplot
ggplot(mean_FST_df_filtered, aes(x = pairwise_comparisons, y = mean_FST)) +
  geom_boxplot() +
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = 0.2) +
  labs(x = "Pairwise Comparisons", y = "Average genome-wide Fst") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_minimal()+
  coord_flip()  # Flip the coordinates to make the plot horizontal

# Create heatmap for provinces comparisons
create_heatmap <- function(data, title) {
  
  data$mean_FST<- round(data$mean_FST, 3)
  
  ggplot(data, aes(x = pop2, y = pop1, fill = mean_FST, label = round(mean_FST, 3))) +
    geom_tile() +
    geom_text(color = "black") +
    scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(data$mean_FST), max(data$mean_FST))) +  # Adjust scale limits
    theme_minimal()+
    labs(x = "", y = "")
}

# Split pairwise_comparisons into pop1 and pop2
library(tidyr)
mean_FST_df <- separate(mean_FST_df, pairwise_comparisons, into = c("pop1", "pop2"), sep = "_vs_")

# Reorder pop1 and pop2 columns based on provinces order
provinces <- c("Niassa", "Cabo_Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
provinces <- rev(provinces)

mean_FST_df$pop1 <- factor(mean_FST_df$pop1, levels = provinces)
mean_FST_df$pop2 <- factor(mean_FST_df$pop2, levels = provinces)

heatmap_2022_provinces <- create_heatmap(mean_FST_df)
print(heatmap_2022_provinces)

#ggsave("heatmap_fst_provinces.png", heatmap_2022_regions, width = 12, height = 10, bg = "white")






#SLIGHTLY DIFFERENT CALCULATION, NO CONFIDENCE INTERVALS. KEEPING IT JUST IN CASE...

# #genome-wide average for Hs1 and Hs2 for each pop
# genome_wide_fst_results <- fst_results_df %>%
#   group_by(pop1, pop2) %>%
#   summarize(genome_wide_avg_Hs1 = mean(Hs1),
#             genome_wide_avg_Hs2 = mean(Hs2),
#             n1 = unique(n1),
#             n2 = unique(n2))
# 
# #calculate Hs for populations PER LOCUS
# genome_wide_fst_results$Hs <- (genome_wide_fst_results$genome_wide_avg_Hs1 + genome_wide_fst_results$genome_wide_avg_Hs2) / 2
# genome_wide_fst_results$Hs <- round(as.numeric(genome_wide_fst_results$Hs),3)
# 
# #Calculate Ht (uses sample sizes for each pop)
# genome_wide_fst_results$Ht <- ((genome_wide_fst_results$n1 * genome_wide_fst_results$genome_wide_avg_Hs1)+ (genome_wide_fst_results$n2 * genome_wide_fst_results$genome_wide_avg_Hs2))/(genome_wide_fst_results$n1 + genome_wide_fst_results$n2)
# genome_wide_fst_results$Ht <- round(as.numeric(genome_wide_fst_results$Ht),3)
# 
# #calcualte Fst
# genome_wide_fst_results$Fst <- (genome_wide_fst_results$Ht -genome_wide_fst_results$Hs)/genome_wide_fst_results$Ht #the same as 1- (genome_wide_fst_results$Hs/genome_wide_fst_results$Ht)
# genome_wide_fst_results<- genome_wide_fst_results[!is.na(genome_wide_fst_results$Fst),] ## IDK WHY NAs ARE INTRODUCED DURING COMPARISONS....
# 
# # Function to create heatmap
# create_heatmap <- function(data, title) {
#   ggplot(data, aes(x = pop2, y = pop1, fill = Fst, label = round(Fst, 4))) +
#     geom_tile() +
#     geom_text(color = "black") +
#     scale_fill_gradient(low = "lightblue1", high = "orange", limits = c(min(data$Fst), max(data$Fst))) +  # Adjust scale limits
#     theme_minimal()
# }
# 
# # Create heatmap for provinces comparisons
# 
# provinces <- c("Niassa", "Cabo Delgado", "Nampula", "Zambezia", "Tete", "Manica_Dry", "Manica_Rainy", "Sofala", "Inhambane", "Maputo_Dry", "Maputo_Rainy") #ordered from north to south
# provinces <- rev(provinces)
# 
# # Remove "_2022" from population names and reorder according to the specified order
# genome_wide_fst_results$pop1 <- gsub("_2022", "", genome_wide_fst_results$pop1)
# genome_wide_fst_results$pop2 <- gsub("_2022", "", genome_wide_fst_results$pop2)
# 
# # Reorder pop1 and pop2 columns based on provinces order
# genome_wide_fst_results$pop1 <- factor(genome_wide_fst_results$pop1, levels = provinces)
# genome_wide_fst_results$pop2 <- factor(genome_wide_fst_results$pop2, levels = provinces)
# 
# heatmap_2022_provinces <- create_heatmap(genome_wide_fst_results)
# print(heatmap_2022_provinces)




# linear regression and correlation coeff of He vs eMOI
