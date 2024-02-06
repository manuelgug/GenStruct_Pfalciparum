
library(haven)
library(progress)
library(fs)
library(moire)
library(ggplot2)
library(dplyr)

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
    nida_values <- grep(paste(db$NIDA2, collapse = "|"), truncated_values, value = TRUE)
    
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
  
  # Update the modified data frame back in the list
  allele_data_list[[i]] <- df
}

names(allele_data_list) <- runs

saveRDS(allele_data_list, "allele_data_list.RDS")

######################################################################
#----------------------------ANALYZE DATA----------------------------#
######################################################################

#######################################################
# 3.- calculate MOI, eMOI and relatedness for each run
#######################################################

# moire has to be run on all data for MOI/eMOI calculations (not sure... takes MONTHS!! ) AND separately for province and region for each year for He calculations;

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already
if (!("n.alleles" %in% colnames(combined_df))){
  combined_df <- combined_df %>%
    group_by(sample_id, locus) %>%
    mutate(n.alleles = n_distinct(allele))
}

# merge with metadata
colnames(combined_df)[1]<- c("NIDA2")
combined_df_merged <- merge(combined_df, db[c("NIDA2", "year", "province", "region")], by="NIDA2", all.x = T)

# delete rows that have NA in "year", "province" and "region" columns: those are samples NOT in the db, thus, not to be incorporated into the analysis.
combined_df_merged <- combined_df_merged %>%
  filter(!is.na(year) & !is.na(province) & !is.na(region))
  
colnames(combined_df_merged)[1]<- c("sample_id")

# # MOIRE ON ALL SAMPLES (PROBABLY WON'T DO THIS, TAKES FOREVER)
# # set MOIRE parameters
# dat_filter <- moire::load_long_form_data(combined_df_merged)
# burnin <- 1e4
# num_samples <- 1e4
# pt_chains <- seq(1, .5, length.out = 20)
#   
# # run moire
# mcmc_results <- moire::run_mcmc(
#   dat_filter, is_missing = dat_filter$is_missing,
#   verbose = TRUE, burnin = burnin, samples_per_chain = num_samples,
#   pt_chains = pt_chains, pt_num_threads = length(pt_chains),
#   thin = 10
# )
#   
# # checkpoint
# saveRDS(mcmc_results, "ALL_SAMPLES_MOIRE-RESULTS.RDS")
# 
# eff_coi <- moire::summarize_effective_coi(mcmc_results)
# naive_coi <- moire::summarize_coi(mcmc_results)
# he <- moire::summarize_he(mcmc_results)


## CALCULATE SAMPLE SIZES FOR EACH PAIR OF VARIABLES


# MOIRE ON EACH PROVINCE DURING 2021 (loop)
#subset 2021 data, loop through each province

# MOIRE ON EACH REGION DURING 2021 (loop)
#subset 2021 data, loop through each region

# MOIRE ON EACH PROVINCE DURING 2022 (loop)
#subset 2022 data, loop through each province

# MOIRE ON EACH REGION DURING 2022 (loop)
#subset 2022 data, loop through each region



#######################################################
# 4.- Calculate MOI/eMOI overall and means per province and region for each year
#######################################################

#import  moire results (ONLY REGION DATA SINCE IT HAS MORE SAMPLES PER RUN THAN PROVINCE, THUS BETTER MOI ESTIMATION!!):
region_rds_files <- list.files(pattern = "\\MOIRE-RESULTS.RDS$", full.names = TRUE)

moire_results_list <- list()

# Load each RDS file into the list with the file name as the list name
for (file in region_rds_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  moire_results_list[[file_name]] <- readRDS(file)
}

# Initialize an empty list to store the processed coi_results
processed_coi_results <- data.frame()

# Loop through each element in moire_results_list
for (i in seq_along(moire_results_list)) {
  # Summarize effective coi and naive coi
  eff_coi <- moire::summarize_effective_coi(moire_results_list[[i]])
  naive_coi <- moire::summarize_coi(moire_results_list[[i]])
  
  # Merge the summaries by sample_id
  coi_results <- merge(eff_coi, naive_coi, by = "sample_id") 
  
  # Add the processed coi_results to the list
  processed_coi_results <- rbind(processed_coi_results, coi_results)
}

# label mono and poly infections. NOTE: "proportion of polyclonal infections (eMOI>1.1)" from Nanna's manuscript
processed_coi_results$polyclonal_from_ecoi_med <- ifelse(processed_coi_results$post_effective_coi_med > 1.1, "polyclonal", "monoclonal")






#######################################################
# 5.- Calculate He and Fws per locus and means per province and region 
#######################################################



# # calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
# combined_df_merged <- combined_df_merged %>%
#   group_by(NIDA2, locus) %>%
#   mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))

# calculate heterozygosity of the population (He): pop = province, region


# calculate fixation index (Fws)

# linear regression and correlation coeff of He vs eMOI




