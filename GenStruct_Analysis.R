
library(haven)
library(progress)
library(fs)
library(moire)
library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)
library(RColorBrewer)


# Objective 3: Describe genetic structure of P. falciparum population in Mozambique in 2021 and 2022 and investigate the origin of parasites with variants of concern

### INDEX

# a) Determine genetic diversity of Plasmodium falciparum in Mozambique
# • Intra-host diversity
#   ◦ Calculate MOI overall and means per province and region for each year (Naïve and eMOI, which incorporates intra-host relatedness between clones)
#   ◦ Calculate % of polyclonal infections per province and region for each year
#   ◦ Fixation index Fws??
# • Population genetic diversity: 
#   ◦ Calculate He per locus per province and region
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
  
  colnames(df)[1] <- "sample_id"
  
    df <- df %>% ### TEHERE IS A BUG WITH THE MASK THAT GENERATES MORE ALLELES THAN THERE REALLY ARE. THIS SNIPPET OF CODE COLLAPSES REPETITIONSAND SUMS THE READS AND FREQS. CRITICAL!!!
      group_by(sample_id, locus, pseudo_cigar) %>%
      summarize(reads = sum(reads),
                norm.reads.locus = sum(norm.reads.locus)) %>%
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
# 3.- GENOMIC + DB MERGING for He calculation and beyond
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
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region", "run_id")], by="NIDA2", all.y = T)

# delete rows that have NA in "year", "province" and "region" columns: those are samples NOT in the db, thus, not to be incorporated into the analysis.
combined_df_merged <- combined_df_merged %>%
  filter(!is.na(year) & !is.na(province) & !is.na(region))

#sanity check
if( sum(!(combined_df_merged$NIDA2 %in% db$NIDA2)) == 0){
  print("All nidas in combined_merged_df are also in the metadata db ✅")
}

#save allele_data_list
saveRDS(combined_df_merged, "combined_df_merged.RDS")

######################################################################
#----------------------------ANALYZE DATA----------------------------#
######################################################################

#######################################################
# 4.- calculate MOI and eMOI for each run
#######################################################
# calculating from each run even though not all samples are gonna be used for the analysis. more info == better moi calculation so np

# allele_data_list <- readRDS("allele_data_list.RDS") 
# 
# # Iterate over each element in allele_data_list
# for (i in seq_along(allele_data_list)) {
#   run <- names(allele_data_list)[i]
#   df <- allele_data_list[[i]]
#   print(run)
#   
#   # set MOIRE parameters
#   dat_filter <- moire::load_long_form_data(df)
#   burnin <- 1e4
#   num_samples <- 1e4
#   pt_chains <- seq(1, .5, length.out = 20)
#   
#   # run moire
#   mcmc_results <- moire::run_mcmc(
#     dat_filter, is_missing = dat_filter$is_missing,
#     verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
#     pt_chains = pt_chains, pt_num_threads = length(pt_chains),
#     thin = 10
#   )
#   
#   # checkpoint
#   saveRDS(mcmc_results, paste0(run, "_MOIRE-RESULTS.RDS"))
# }

combined_df_merged <- readRDS("combined_df_merged.RDS") 

colnames(combined_df_merged)[1]<- c("sample_id")

# set MOIRE parameters
dat_filter <- moire::load_long_form_data(combined_df_merged)
burnin <- 2.5e3
num_samples <- 2.5e3
pt_chains <- seq(1, .5, length.out = 20)

# run moire
mcmc_results <- moire::run_mcmc(
  dat_filter, is_missing = dat_filter$is_missing,
  verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
  pt_chains = pt_chains, pt_num_threads = length(pt_chains),
  thin = 10)

saveRDS(mcmc_results, "all_1329_samples_MOIRE-RESULTS.RDS")

#######################################################
# 5.- Present MOI/eMOI results overall and means per province and region for each year
#######################################################

#import  moire results:
# rds_files <- list.files(pattern = "\\MOIRE-RESULTS.RDS$", full.names = TRUE)

# moire_results_list <- list()
# 
# # Load each RDS file into the list with the file name as the list name
# for (file in rds_files) {
#   print(file)
#   file_name <- tools::file_path_sans_ext(basename(file))
#   moire_results_list[[file_name]] <- readRDS(file)
# }
# 
# 
# # Initialize an empty list to store the processed coi_results
# processed_coi_results <- data.frame()
# 
# # Loop through each element in moire_results_list
# for (i in seq_along(moire_results_list)) {
#   # Summarize effective coi and naive coi
#   eff_coi <- moire::summarize_effective_coi(moire_results_list[[i]])
#   naive_coi <- moire::summarize_coi(moire_results_list[[i]])
#   
#   # Merge the summaries by sample_id
#   coi_results <- merge(eff_coi, naive_coi, by = "sample_id")   #[c("sample_id", "post_effective_coi_mean", "post_effective_coi_med", "naive_coi")]
#   
#   # Add the processed coi_results to the list
#   processed_coi_results <- rbind(processed_coi_results, coi_results)
# }


mcmc_results <- readRDS("TEST_all_1329_samples_MOIRE-RESULTS.RDS") 

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
# RAREFACTION CURVES

# raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
#                                    year = combined_df_merged$year, 
#                                    province = combined_df_merged$province,
#                                    region = combined_df_merged$region,
#                                    locus = combined_df_merged$locus,
#                                    n.alleles = combined_df_merged$n.alleles,
#                                    allele = paste0(combined_df_merged$locus, "_", combined_df_merged$pseudo_cigar)))
# 
# raref_input <- raref_input %>% distinct()
# 
# #subsetting #INIT LOOP
# sub <- raref_input[raref_input$year == 2022 & raref_input$province =="Maputo",] #iterate this 
# 
# # Cast the dataframe to wide format
# raref_df <- dcast(sub, NIDA2 ~ locus, value.var = "n.alleles")
# raref_df <- raref_df[, -1]
# raref_df <- apply(raref_df, 2, function(x) as.numeric(as.character(x)))
# 
# 
# # CALCULATE CURVE
# accum_curve <-specaccum(raref_df, 'random', permutations = 1000, method = "rarefaction")
# plot(accum_curve, xlab = "Samples")

###########################################
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
plot(accum_curves_2021[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2021 (per Province)", xlim = c(0,85), ylim = c(0,2700))
for (i in 2:length(accum_curves_2021)) {
  lines(accum_curves_2021[[i]], col = colors[i], lw = 1.5)
}
legend(x = 65, y = 550, legend = names(accum_curves_2021), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022 (per Province)", xlim = c(0,200), ylim = c(0,3000))
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
plot(accum_curves_2021[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2021 (per Region)", xlim = c(0,150), ylim = c(0,3000))
for (i in 2:length(accum_curves_2021)) {
  lines(accum_curves_2021[[i]], col = colors[i], lw = 1.5)
}
legend(x = 120, y = 550, legend = names(accum_curves_2021), fill = colors, x.intersp = 0.7, y.intersp = 0.5)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022 (per Region)", xlim = c(0,500), ylim = c(0,4000))
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
pca_labels<- db %>% distinct(NIDA2, year, province, region)

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


# ### MAKE THE PCA A FUNCTION ON WHICH YOU CAN CHOOSE INPUT: (PRESENCE/ABSENCE OR FREQS) AND THE VARIABLES TO TEST (JOINED OR SIMPLE)
# 
# # Perform PCA on rearranged
# pca_result <- prcomp(rearranged_filtered, scale. = FALSE)
# pc_scores <- as.data.frame(pca_result$x)
# pca_data <- cbind(pc_scores, pca_labels)  # Make sure pca_labels is defined
# 
# # Extract PCA results
# pcs <- as.data.frame(pca_result$x)  # principal components
# variance <- pca_result$sdev^2  # variance of each principal component
# prop_variance <- variance / sum(variance)  # proportion of variance explained by each component
# 
# # Generate contrasting colors from RColorBrewer palette
# num_colors <- length(unique(factor(paste0(pca_data$region, "_", pca_data$year))))
# set.seed(690)
# Set1_colors <- brewer.pal(9, "Set1")
# Set2_colors <- brewer.pal(8, "Set2")
# mixed_colors <- c(Set1_colors, Set2_colors)
# random_colors <- sample(mixed_colors, num_colors)
# 
# # Plot PCA with ggplot including sample labels
# ggplot(pca_data, aes(PC1, PC2, color = factor(paste0(region, "_", year)), label = rownames(pca_data))) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(title = "PCA of Genetic Content",
#        x = paste0("PC1 (", round(prop_variance[1] * 100, 2), "%)"),
#        y = paste0("PC2 (", round(prop_variance[2] * 100, 2), "%)")) +
#   scale_color_manual(values = random_colors) +
#   theme_minimal()
 
#TSNE
library(Rtsne)

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

### TO DO:
# 1.- check issue with min(coi) = 0 (something about having alleles with 0 all across samples maybe? would also explain the need to remove 0 columns on the tsne/pca)
# 2.- make image good!

library(dcifer)
pardef <- par(no.readonly = TRUE)

combined_df_merged <- readRDS("combined_df_merged.RDS")

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
coi[coi == 0] <- 1 #COPE, I DON'T KNOW WHY I'M GETTING COI OF ZERO SOMETIMES (probably some formatting issue with the input df... CHECK IN DEPTH. THIS ALLOWS ibDat FUNCTION TO RUN WITHOUT ERRORS, BUT I DON'T KNOW WHAT'S GOING ON

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#calculate ibd
dres0_2021 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                 alpha = 0.05, nr = 1e3)  

par(mar = c(3, 3, 1, 1))
alpha <- 0.05                          # significance level                    
dmat <- dres0_2021[, , "estimate"]
# create symmetric matrix
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))]  
# determine significant, reverse columns for upper triangle
isig <- which(dres0_2021[, , "p_value"] <= alpha, arr.ind = TRUE)[, 2:1] 
plotRel(dmat, isig = isig, draw_diag = TRUE, lwd_diag = 0.5, idlab = TRUE, 
        col_id = c(1:10)[factor(meta$province)]) 


## 2022 samples ##

#format data
dsmp <- formatDat(combined_df_merged_2022, svar = "NIDA2", lvar = "locus", avar = "pseudo_cigar")
str(dsmp, list.len = 2)

# format metadata
meta <- unique(combined_df_merged_2021[c("NIDA2", "region", "province")])
meta <- meta[match(names(dsmp), meta$NIDA2), ]  # order samples as in dsmp

#estimate naive coi
lrank <- 2
coi   <- getCOI(dsmp, lrank = lrank)
min(coi)
coi[coi == 0] <- 1 #COPE, I DON'T KNOW WHY I'M GETTING COI OF ZERO SOMETIMES (probably some formatting issue with the input df... CHECK IN DEPTH. THIS ALLOWS ibDat FUNCTION TO RUN WITHOUT ERRORS, BUT I DON'T KNOW WHAT'S GOING ON

#estimate allele freqs
afreq <- calcAfreq(dsmp, coi, tol = 1e-5) 
str(afreq, list.len = 2)

#calculate ibd
dres0_2022 <- ibdDat(dsmp, coi, afreq,  pval = TRUE, confint = TRUE, rnull = 0, 
                alpha = 0.05, nr = 1e3)  

par(mar = c(3, 3, 1, 1))
alpha <- 0.05                          # significance level                    
dmat <- dres0_2022[, , "estimate"]
# create symmetric matrix
dmat[upper.tri(dmat)] <- t(dmat)[upper.tri(t(dmat))]  
# determine significant, reverse columns for upper triangle
isig <- which(dres0_2022[, , "p_value"] <= alpha, arr.ind = TRUE)[, 2:1] 
plotRel(dmat, isig = isig, draw_diag = TRUE, lwd_diag = 0.5, idlab = TRUE, 
        col_id = c(1:10)[factor(meta$province)]) 

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

# calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
combined_df_merged <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))

# calculate heterozygosity of the population (He): pop = province, region

# calculate fixation index (Fws)

# linear regression and correlation coeff of He vs eMOI
