
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

# moire has to be run on all data for MOI/eMOI calculations AND separately for province and region for each year; so A TOTAL OF 5 TIME

# Iterate over each element in allele_data_list
for (i in seq_along(allele_data_list)) {
  run <- names(allele_data_list)[i]
  df <- allele_data_list[[i]]
  
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
  saveRDS(mcmc_results, paste0(run, "_MOIRE-RESULTS.RDS"))
}


#######################################################
# 4.- calculate He for each population (per year per region/province)
#######################################################

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
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


## CHECK SAMPLE SIZES FOR EACH PAIR OF VARIABLES: sample size affects He calculation, probably will need rarefactions or something similar
sample_size_provinces <- combined_df_merged %>%
  group_by(year, province) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))

sample_size_regions <- combined_df_merged %>%
  group_by(year, region) %>%
  summarise(unique_NIDA2_count = n_distinct(NIDA2))




###########################################
# RAREFACTION CURVES

raref_input <- as.data.frame(cbind(NIDA2 = combined_df_merged$NIDA2, 
                                   year = combined_df_merged$year, 
                                   province = combined_df_merged$province,
                                   region = combined_df_merged$region,
                                   locus = combined_df_merged$locus,
                                   n.alleles = combined_df_merged$n.alleles))

raref_input <- raref_input %>% distinct()

#subsetting #INIT LOOP
sub <- raref_input[raref_input$year == 2022 & raref_input$province =="Maputo",] #iterate this 

# Cast the dataframe to wide format
raref_df <- dcast(sub, NIDA2 ~ locus, value.var = "n.alleles")
raref_df <- raref_df[, -1]
raref_df <- apply(raref_df, 2, function(x) as.numeric(as.character(x)))


# CALCULATE CURVE
accum_curve <-specaccum(raref_df, 'random', permutations = 1000, method = "rarefaction")
plot(accum_curve, xlab = "Samples")

###########################################
#ACCUMULATION CURVES

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
plot(accum_curves_2021[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2021", xlim = c(0,85), ylim = c(0,2700))
for (i in 2:length(accum_curves_2021)) {
  lines(accum_curves_2021[[i]], col = colors[i], lw = 1.5)
}
legend("bottomright", legend = names(accum_curves_2021), fill = colors)

# Plot the curves for 2022
plot(accum_curves_2022[[1]], col = colors[1], xlab = "Samples", main = "Accumulation Curves for 2022", xlim = c(0,200), ylim = c(0,3000))
for (i in 2:length(accum_curves_2022)) {
  lines(accum_curves_2022[[i]], col = colors[i], lw = 1.5)
}
legend("bottomright", legend = names(accum_curves_2022), fill = colors)

###########################################









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
# 5.- Present MOI/eMOI results overall and means per province and region for each year
#######################################################

#import  moire results:
rds_files <- list.files(pattern = "\\MOIRE-RESULTS.RDS$", full.names = TRUE)

moire_results_list <- list()

# Load each RDS file into the list with the file name as the list name
for (file in rds_files) {
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
  coi_results <- merge(eff_coi, naive_coi, by = "sample_id")   #[c("sample_id", "post_effective_coi_mean", "post_effective_coi_med", "naive_coi")]
  
  # Add the processed coi_results to the list
  processed_coi_results <- rbind(processed_coi_results, coi_results)
}


# label mono and poly infections. NOTE: "proportion of polyclonal infections (eMOI>1.1)" from Nanna's manuscript
processed_coi_results$polyclonal_from_ecoi_med <- ifelse(processed_coi_results$post_effective_coi_med > 1.1, "polyclonal", "monoclonal")

#merge with categorical variables, controls are removed automatically
colnames(processed_coi_results)[1] <- "NIDA2"
processed_coi_results <- merge(processed_coi_results, db[c("NIDA2", "year", "province", "region")], by="NIDA2")


# % polyclonal infections on each province and region per year
polyclonal_percentage_region <- processed_coi_results %>%
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

polyclonal_percentage_province <- processed_coi_results %>%
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
processed_coi_results$year <- factor(processed_coi_results$year)

c <- ggplot(processed_coi_results, aes(x = region, y = post_effective_coi_med, fill = year)) +
  geom_boxplot() +
  labs(x = "Region", y = "Post Effective COI Median") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
c

d <-ggplot(processed_coi_results, aes(x = province, y = post_effective_coi_med, fill = year)) +
  geom_boxplot() +
  labs(x = "Province", y = "Post Effective COI Median") +
  scale_fill_manual(values = c("2021" = "cyan3", "2022" = "orange")) +  # Adjust colors as needed
  theme_minimal()+
  ylim(0, NA)
d


#######################################################
# 6.- Calculate He and Fws per locus and means per province and region 
#######################################################

allele_data_list <- readRDS("allele_data_list.RDS")

# concat all dataframes together
combined_df <- bind_rows(allele_data_list)

# calculate n.alleles for each locus of each sample if not done already during contaminant filtering
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

# calculate heterozygosity of the individual (Hw): 𝐻W = 1 − (n𝑖(1/n𝑖)**2) 
combined_df_merged <- combined_df_merged %>%
  group_by(NIDA2, locus) %>%
  mutate(Hw = 1 - (n.alleles * (1/n.alleles)^2))

# calculate heterozygosity of the population (He): pop = province, region

# calculate fixation index (Fws)

# linear regression and correlation coeff of He vs eMOI
