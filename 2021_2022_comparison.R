

library(ggplot2)
library(dplyr)
library(ggsignif)


regions <- c("North", "Centre", "South")

######################################################################################################
##### COMPARE eCOI #####

coi_2021 <- read.csv("coi_results_2021.csv")
coi_2022 <- read.csv("coi_results_2022.csv")

coi_all <- rbind(coi_2021, coi_2022)

coi_all$region <- factor(coi_all$region, levels = regions)

# Define a function to perform pairwise comparisons for each province
pairwise_comparisons <- function(data) {
  comparisons <- combn(unique(data$year), 2, simplify = FALSE)
  test_results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(post_effective_coi_med ~ year, data = data[data$year %in% pair,])
    p_value <- test$p.value
    data.frame(region = unique(data$region), group1 = pair[1], group2 = pair[2], p_value = p_value)
  })
  do.call(rbind, test_results)
}

# Calculate significance values
significance_data <- coi_all %>%
  group_by(region) %>%
  do(pairwise_comparisons(.)) %>%
  ungroup()

# Prepare the data for geom_signif
significance_data <- significance_data %>%
  mutate(y_position = max(coi_all$post_effective_coi_med) + 0.5)  # Adjust y_position as needed

significance_data

# Plot with significance annotations
ECOI_PLOT <- ggplot(coi_all, aes(x = region, y = post_effective_coi_med, color = as.factor(year)))+
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
  labs(x = "Region", y = "Post Effective COI Median", color = "Year", fill = "Year", title = "Complexity of Infection") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################################################################################################
##### COMPARE %polyclonal #####

perc_polyclonal <-coi_all %>%
  group_by(region, year) %>%
  summarize(perc_polyclonal =  sum(post_effective_coi_med > 1.1) / length(post_effective_coi_med))

perc_polyclonal$region <- factor(perc_polyclonal$region, levels = regions)

# Plot with significance annotations
PERC_POLYCLONAL_PLOT <- ggplot(perc_polyclonal, aes(x = region, y = perc_polyclonal, fill = as.factor(year))) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
  labs(x = "Region", y = "%Polyclonal Samples", fill = "Year", title = "Polyclonal Infections") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######################################################################################################
##### COMPARE He #####

he_2021 <- read.csv("He_results_regions_2021.csv")
he_2021$year  <- 2021 
he_2022 <- read.csv("He_results_regions_2022.csv")
he_2022$year  <- 2022 

he_all <- rbind(he_2021, he_2022)

he_all<- he_all[he_all$geo == "region",]

he_all$population <- factor(he_all$population, levels = regions)


# ##################33
# library(emmeans)
# 
# # Fit separate linear models for each region
# region_models <- list()
# for (population in unique(he_all$population)) {
#   region_data <- he_all[he_all$population == population, ]
#   region_models[[population]] <- lm(post_stat_mean ~ year, data = region_data)
# }
# 
# # Perform pairwise comparisons between years within each region
# pairwise_comparisons <- list()
# for (population in unique(he_all$population)) {
#   region_model <- region_models[[population]]
#   pairwise_comparisons[[population]] <- emmeans(region_model, pairwise ~ year, adjust = "tukey")
# }
# 
# pairwise_comparisons
# ##################33

##################33
library(emmeans)
library(lmerTest)

# Fit separate linear models for each region
region_models <- list()
for (population in unique(he_all$population)) {
  region_data <- he_all[he_all$population == population, ]
  region_models[[population]] <- lme(post_stat_mean ~ year, 
                                     random = ~ 1 | locus, #locus is random effect
                                     data = region_data)
}


# Loop through each element in region_models
model_results <- list()

for (population in names(region_models)) {

  model <- region_models[[population]]
  
  summ <- summary(model)
  an <- anova(model)
  ci_mod <- intervals(model, which = "fixed")
  
  model_results[[paste0(population, "_summary")]] <- summ
  model_results[[paste0(population, "_anova")]] <- an
  model_results[[paste0(population, "_CI")]] <- ci_mod
}


# Perform pairwise comparisons between years within each region
pairwise_comparisons <- list()
for (population in unique(he_all$population)) {
  region_model <- region_models[[population]]
  pairwise_comparisons[[population]] <- emmeans(region_model, pairwise ~ year, adjust = "tukey")
}

pairwise_comparisons
##################33

# 
# # Define a function to perform pairwise comparisons for each province
# pairwise_comparisons <- function(data) {
#   comparisons <- combn(unique(data$year), 2, simplify = FALSE)
#   test_results <- lapply(comparisons, function(pair) {
#     test <- wilcox.test(post_stat_med ~ year, data = data[data$year %in% pair,])
#     p_value <- test$p.value
#     data.frame(region = unique(data$population), group1 = pair[1], group2 = pair[2], p_value = p_value)
#   })
#   do.call(rbind, test_results)
# }
# 
# # Calculate significance values
# significance_data <- he_all %>%
#   group_by(population) %>%
#   do(pairwise_comparisons(.)) %>%
#   ungroup()
# 
# # Prepare the data for geom_signif
# significance_data <- significance_data %>%
#   mutate(y_position = max(he_all$post_stat_med) + 0.5)  # Adjust y_position as needed
# 
# significance_data


# Plot with significance annotations
HE_PLOT <- ggplot(he_all, aes(x = population, y = post_stat_med, color = as.factor(year)))+
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
  labs(x = "Region", y = "Post He Median", color = "Year", fill = "Year", title = "Population Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


######################################################################################################
### COMPARE CONNECTIVITY ###

conn_2021 <- read.csv("connectivity_2021.csv")
conn_2022 <- read.csv("connectivity_2022.csv")

conn_2021$year <- 2021
conn_2022$year <- 2022

conn_all <- rbind(conn_2021, conn_2022)

# Define a function to perform pairwise comparisons for each province
pairwise_comparisons <- function(data) {
  comparisons <- combn(unique(data$year), 2, simplify = FALSE)
  test_results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(estimate ~ year, data = data[data$year %in% pair,])
    p_value <- test$p.value
    data.frame(region = unique(data$conn_regions), group1 = pair[1], group2 = pair[2], p_value = p_value)
  })
  do.call(rbind, test_results)
}

# Calculate significance values
significance_data <- conn_all %>%
  group_by(conn_regions) %>%
  do(pairwise_comparisons(.)) %>%
  ungroup()

# Prepare the data for geom_signif
significance_data <- significance_data %>%
  mutate(y_position = max(conn_all$estimate) + 0.5)  # Adjust y_position as needed

significance_data


# Plot with significance annotations
CONN_PLOT <- ggplot(conn_all, aes(x = conn_regions, y = estimate, color = as.factor(year)))+
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
  labs(x = "Region", y = "IBD", color = "Year", fill = "Year", title = "Connectivity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



######################################################################################################
### COMPARE 1-FWS ###

fws_2021 <- read.csv("He_and_Fws_2021.csv")
fws_2022 <- read.csv("He_and_Fws_2022.csv")

fws_all <- rbind(fws_2021, fws_2022)
fws_all <- fws_all[!is.na(fws_all$fws_region),]

fws_all$region <- factor(fws_all$region, levels = regions)

# remove 5% least variant amplicons to avoid issues when calculating Fws
amps_q95<- round(length(unique(fws_all$locus))*0.95)

he_amps<- fws_all %>%
  group_by(locus) %>%
  summarize(He_region = mean(He_region))%>%
  arrange(desc(He_region))

he_amps_topHe <- he_amps[1:amps_q95,]

fws_all<- fws_all[fws_all$locus %in% he_amps_topHe$locus,]

fws_all <- fws_all[,c("NIDA2", "locus", "year", "region", "fws_region")] 

fws_all <- distinct(fws_all)

fws_all <- fws_all %>%
  group_by(NIDA2, year, region) %>%
  summarize(fws_region = mean(fws_region))


# Define a function to perform pairwise comparisons for each province
pairwise_comparisons <- function(data) {
  comparisons <- combn(unique(data$year), 2, simplify = FALSE)
  test_results <- lapply(comparisons, function(pair) {
    test <- wilcox.test(fws_region ~ year, data = data[data$year %in% pair,])
    p_value <- test$p.value
    data.frame(region = unique(data$region), group1 = pair[1], group2 = pair[2], p_value = p_value)
  })
  do.call(rbind, test_results)
}

# Calculate significance values
significance_data <- fws_all %>%
  group_by(region) %>%
  do(pairwise_comparisons(.)) %>%
  ungroup()

# Prepare the data for geom_signif
significance_data <- significance_data %>%
  mutate(y_position = max(fws_all$fws_region) + 0.5)  # Adjust y_position as needed

significance_data


# Plot with significance annotations
FWS_PLOT <- ggplot(fws_all, aes(x = region, y = fws_region, color = as.factor(year)))+
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub("\\..*", "", x)) +  # Remove year from x-axis labels
  labs(x = "Region", y = "Genome-wide 1-Fws", color = "Year", fill = "Year", title = "In-sample Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###################################################################################################
# FINAL FIGURES

library(gridExtra)

combined_plot <- grid.arrange(
  PERC_POLYCLONAL_PLOT, ECOI_PLOT, FWS_PLOT, CONN_PLOT, HE_PLOT,
  ncol = 2, nrow = 3
)

combined_plot


ggsave("comparison_plots.png", combined_plot, width = 16, height = 12)
