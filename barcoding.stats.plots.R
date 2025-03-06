###### PACKAGES #####
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(broom)
library(tidyverse)
library(vegan)
library(readxl)
library(rstatix)
library(ggpubr)
###### Data Import + Carpentry #####
barcode.data <- read.delim("~/Data/barcode.data.txt", comment.char="#")
colnames(barcode.data) <- c("Island","Population","Habitat",
                            "1.annae.evermanni.like",
                            "3.lutea","6.lobata",
                            "7.murrayensis",
                            "8.australiensis",
                            "9.australiensis.compressa")
barcode.long <- barcode.data %>%
  pivot_longer(cols = starts_with(c("1", "3", "6", "7", "8", "9")),
               names_to = "Clade", values_to = "Count")
island.order <- c("Maug", "Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam")
habitat.order <- c("RiverDelta", "Backreef", "Forereef")
clade.colors <-c ("#F8766D","#00FF00", "#FFA500", "#00A9FF", "#FF00FF", "#9D00FF")
barcode.long<- barcode.data %>%
  pivot_longer(cols = starts_with(c("1.", "3.", "6.", "7.", "8.", "9.")),
               names_to = "Clade", values_to = "Count")

###### Plot absolute and relative abundance by island #####
plot.absolute.island <- ggplot(barcode.long, aes(x = factor(Island, levels = island.order), y = Count, fill = Clade)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Absolute Abundance by Island", y = "Absolute Abundance (n)", x = NULL)+
  scale_fill_manual(values = clade.colors)
# Calculate relative abundance by grouping and normalizing counts by island
barcode.relative.island <- barcode.long %>%
  group_by(Island) %>%
  mutate(Relative.Abundance = Count / sum(Count) * 100)
# Plot relative abundance (percentage)
plot.relative.island <- ggplot(barcode.relative.island, aes(x = factor(Island, levels = island.order), y = Relative.Abundance, fill = Clade)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Relative Abundance by Island", y = "Relative Abundance (%)", x = NULL)+
  scale_fill_manual(values = clade.colors)
combined.plot.island <- plot.absolute.island + plot.relative.island + 
  plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = "bottom")
combined.plot.island
###### Save Plot #####
ggsave("~/Results/Barcoding/distribution.by.island.pdf", 
       plot = combined.plot.island, 
       device = "pdf", 
       width = 12, 
       height = 10)

###### Plot absolute and relative abundance by habitat (southern three islands) #####
barcode.long.south <- filter(barcode.long, Island %in% c("Guam", "Rota", "Saipan"))
plot.absolute.habitat <- ggplot(barcode.long.south, aes(x = factor(Habitat, levels = habitat.order), y = Count, fill = Clade)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Absolute Abundance by Habitat", y = "Absolute Abundance (n)", x = NULL)+
  scale_fill_manual(values = clade.colors)
# Calculate relative abundance by grouping and normalizing counts by island
barcode.relative.habitat.south <- barcode.long.south %>%
  group_by(Habitat) %>%
  mutate(Relative.Abundance = Count / sum(Count) * 100)
# Plot relative abundance (percentage)
plot.relative.habitat <- ggplot(barcode.relative.habitat.south, aes(x = factor(Habitat, levels = habitat.order), y = Relative.Abundance, fill = Clade)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Relative Abundance by Habitat", y = "Relative Abundance (%)", x = NULL)+
  scale_fill_manual(values = clade.colors)
combined.plot.habitat.south <- plot.absolute.habitat + plot.relative.habitat + 
  plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = "bottom")
combined.plot.habitat.south
###### Save Plot #####
ggsave("~/Results/Barcoding/distribution.by.habitat.south.pdf", 
       plot = combined.plot.habitat.south, 
       device = "pdf", 
       width = 12, 
       height = 10)
###### Is there an association between Habitat type and Clade Distribution (Southern 3 islands)? Fisher Exact Test #####
library(dplyr)
library(tidyverse)
data <- read_excel("~/Data/BarcodingSummary.xlsx", sheet = "Conclusions")
# Filter the data
data_filtered <- data %>%
  filter(grepl("^[0-9]+$", Clade), Island %in% c("Guam", "Saipan", "Rota"))
data_filtered$Clade <- as.integer(data_filtered$Clade)
summary.table <- data_filtered %>%
  group_by(Habitat, Clade) %>%
  summarize(Count = n()) %>%
  pivot_wider(names_from = Habitat, values_from = Count, values_fill = 0)
sum.stats<- data.frame(summary.table$Backreef, summary.table$Forereef, summary.table$`River Delta`)
rownames(sum.stats)<-c("I", "II", "III", "VI", "VII", "VIII", "IX", "X")
colnames(sum.stats) <- c("Backreef", "Forereef", "River Delta")
filtered_data <- subset(sum.stats, !(row.names(sum.stats) %in% c("II", "X")))
mosaicplot(filtered_data, main = "Mosaic Plot", color=TRUE)
test <- fisher.test(filtered_data, simulate.p.value = TRUE, B = 1e6)
test$p.value
# Interpret the results
alpha <- 0.05
if (test$p.value < alpha) {
  print("Reject the null hypothesis (the distribution of clades across habitats likely isn't random).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chancen).")
}

###### What is each clade's relationship to island: Habitat type? Nested Anova #####
clades <- c("1.annae.evermanni.like", "3.lutea",
            "6.lobata", "7.murrayensis",
            "8.australiensis", "9.australiensis")
nested.results.list <- list()

# Specify the southernmost islands
southern_islands <- c("Guam", "Rota", "Saipan")

for (clade in clades) {
  # Subset the data to only include the current clade and the southernmost islands
  subset_data <- subset(barcode.long, Clade == clade & Island %in% southern_islands)
  
  # Check if Island and Habitat have more than one level
  if (length(unique(subset_data$Island)) > 1 && length(unique(subset_data$Habitat)) > 1) {
    # Perform ANOVA if both factors have more than one level
    aov.result <- aov(Count ~ Island/Habitat, data = subset_data)
    anova.summary <- tidy(aov.result)
    
    # Filter for significant results
    significant.results <- anova.summary %>%
      filter(p.value < 0.05)
    
    nested.results.list[[clade]] <- significant.results
  } else {
    # Skip the clade if there is only one level for either factor
    warning(paste("Skipping clade", clade, ": Island or Habitat has only one level."))
  }
}

nested.results.list
###### Is there an association between Island and Clade Distribution? Fisher Exact Test #####
data <- read_excel("~/Data/BarcodingSummary.xlsx", sheet = "Conclusions")
# Filter the data
data_filtered <- data %>%
  filter(grepl("^[0-9]+$", Clade))
data_filtered$Clade <- as.integer(data_filtered$Clade)
summary.table <- data_filtered %>%
  group_by(Island, Clade) %>%
  summarize(Count = n()) %>%
  pivot_wider(names_from = Island, values_from = Count, values_fill = 0)
sum.stats<- data.frame(summary.table$Maug, summary.table$Pagan, summary.table$Sarigan, summary.table$Saipan, summary.table$Tinian, summary.table$Rota, summary.table$Guam)
rownames(sum.stats)<-c("I", "II", "III", "VI", "VII", "VIII", "IX", "X")
colnames(sum.stats) <- c("Maug", "Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam")
filtered_data <- subset(sum.stats, !(row.names(sum.stats) %in% c("II", "X")))
mosaicplot(filtered_data, main = "Mosaic Plot", color=TRUE)
test <- fisher.test(filtered_data, simulate.p.value = TRUE, B = 1e6)
test$p.value
# Interpret the results
alpha <- 0.05
if (test$p.value < alpha) {
  print("Reject the null hypothesis (the distribution of clades across islands likely isn't random).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chancen).")
}
###### Is there a significant difference between River Delta and Forereef (Southern 3 islands)? PERMANOVA #####
data <- read_excel("~/Data/BarcodingSummary.xlsx", sheet = "bypopulation.RDvsFR")
set.seed(1)
# Step 1: Reshape to wide format for distance calculation
long_data <- data %>%
  pivot_longer(
    cols = I:IX,  # specify the range of species columns
    names_to = "Species",  # new column for species names
    values_to = "Count"  # new column for counts
  )
wide_data <- long_data %>%
  pivot_wider(
    names_from = Species,  # Species become column names
    values_from = Count,    # Values are counts
    values_fill = 0         # Fill NAs with 0
  )
# Select only numeric columns (excluding Population and Habitat)
distance_matrix <- vegdist(wide_data %>% select(-Population, -Habitat), method = "bray")


# Step 3: Prepare habitat labels
habitat_labels <- wide_data$Habitat  # Extract the habitat labels
# Step 4: Run PERMANOVA using adonis2
permanova_result <- adonis2(distance_matrix ~ habitat_labels)
# View the results
print(permanova_result)
# Interpret the results
alpha <- 0.05
if (permanova_result$`Pr(>F)`[1] < alpha) {
  print("Reject the null hypothesis (the distribution of clades is signifigantly diffferent).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chance).")
}
permanova_result$`Pr(>F)`[1]
###### Is there a significant difference between River Delta and Backreef (Southern 3 islands)? PERMANOVA #####
data <- read_excel("~/Data/BarcodingSummary.xlsx", sheet = "bypopulation.RDvsBR")
set.seed(2)
# Step 1: Reshape to wide format for distance calculation
long_data <- data %>%
  pivot_longer(
    cols = I:IX,  # specify the range of species columns
    names_to = "Species",  # new column for species names
    values_to = "Count"  # new column for counts
  )
wide_data <- long_data %>%
  pivot_wider(
    names_from = Species,  # Species become column names
    values_from = Count,    # Values are counts
    values_fill = 0         # Fill NAs with 0
  )
# Select only numeric columns (excluding Population and Habitat)
distance_matrix <- vegdist(wide_data %>% select(-Population, -Habitat), method = "bray")
# Step 3: Prepare habitat labels
habitat_labels <- wide_data$Habitat  # Extract the habitat labels
# Step 4: Run PERMANOVA using adonis2
permanova_result <- adonis2(distance_matrix ~ habitat_labels)
# View the results
print(permanova_result)
# Interpret the results
alpha <- 0.05
if (permanova_result$`Pr(>F)`[1] < alpha) {
  print("Reject the null hypothesis (the distribution of clades is signifigantly diffferent).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chance).")
}
permanova_result$`Pr(>F)`[1]


###### Is there a significant difference between Backreef and Forereef (Southern 3 islands)? PERMANOVA #####
set.seed(1)
data <- read_excel("~/Data/BarcodingSummary.xlsx", sheet = "bypopulation.BRvsFR")
# Step 1: Reshape to wide format for distance calculation
long_data <- data %>%
  pivot_longer(
    cols = I:IX,  # specify the range of species columns
    names_to = "Species",  # new column for species names
    values_to = "Count"  # new column for counts
  )
wide_data <- long_data %>%
  pivot_wider(
    names_from = Species,  # Species become column names
    values_from = Count,    # Values are counts
    values_fill = 0         # Fill NAs with 0
  )
# Select only numeric columns (excluding Population and Habitat)
distance_matrix <- vegdist(wide_data %>% select(-Population, -Habitat), method = "bray")
# Step 3: Prepare habitat labels
habitat_labels <- wide_data$Habitat  # Extract the habitat labels
# Step 4: Run PERMANOVA using adonis2
permanova_result <- adonis2(distance_matrix ~ habitat_labels)
# View the results
print(permanova_result)
# Interpret the results
alpha <- 0.05
if (permanova_result$`Pr(>F)`[1] < alpha) {
  print("Reject the null hypothesis (the distribution of clades is signifigantly diffferent).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chance).")
}
permanova_result$`Pr(>F)`[1]

###### Is there a significant different among sites within a habitat? - Fisher's exact test #####
clade_counts <- matrix(
  c(0, 0,  # Clade I
    1, 2,  # Clade III
    0, 0,  # Clade VI
    3, 5,  # Clade VII
    5, 7,  # Clade VIII
    4, 5), # Clade IX
  nrow = 6,  # Number of clades
  byrow = TRUE,  # Fill matrix by rows
  dimnames = list(
    Clade = c("Clade I", "Clade III", "Clade VI", "Clade VII", "Clade VIII", "Clade IX"),  # Row names (clades)
    Site = c("Jimmy Dees", "Pago Bay")  # Column names (sites)
  )
)
test <- fisher.test(clade_counts, simulate.p.value = TRUE, B = 1e6)
test$p.value
# Interpret the results
alpha <- 0.05
if (test$p.value < alpha) {
  print("Reject the null hypothesis (the distribution of clades across islands likely isn't random).")
} else {
  print("Fail to reject the null hypothesis (the observed distribution could be due to random chance).")
}
print(test)

###### Is there a significant different among sites within a habitat across islands? - Kruskal-Wallis Test #####
species_counts <- data.frame( # example data structure
  JimmyDees = c(0, 2, 0, 5, 7, 5),
  MarineLab = c(0, 1, 0, 3, 5, 4),
  MerizoPier = c(4, 1, 0, 3, 3, 2),
  RotaBR = c(1, 3, 0, 5, 6, 0),
  SaipanBR = c(0, 1, 10, 1, 9, 0)
)
species_counts_long <- stack(species_counts)
kruskal_test <- kruskal.test(values ~ ind, data = species_counts_long)
kruskal_test

###### Is there a significant different among sites across river deltas? - Kruskal-Wallis Test #####
species_counts <- data.frame( # example data structure
  FouhaBay = c(0, 2, 0, 27, 0, 0),
  Inarahan = c(0, 2, 0, 13, 11, 0),
  Talofofo = c(0, 1, 0, 18, 3, 0)
)
species_counts_long <- stack(species_counts)
kruskal_test <- kruskal.test(values ~ ind, data = species_counts_long)
kruskal_test


###### Is there a significant different among sites across fore reefs? - Kruskal-Wallis Test #####
species_counts <- data.frame( # example data structure
  FouhaBay = c(1, 4, 0, 2, 5, 8),
  Inarahan = c(1,9,2,3,8,0),
  Talofofo = c(0,3,4,0,15,1),
  Ritidian = c(1,4,1,1,19,0),
  Rota = c(1,14,0,3,11,0),
  Saipan = c(0,11,0,2,9,0)
)
species_counts_long <- stack(species_counts)
kruskal_test <- kruskal.test(values ~ ind, data = species_counts_long)
kruskal_test
###### Does Diversity Correlate to Island Size? Linear Regression Model ######
# Data curation + Shannon Calculation
barcode_summarized <- barcode.long %>%
  group_by(Island, Clade) %>%
  summarize(TotalCount = sum(Count, na.rm = TRUE)) %>%
  ungroup()
barcode_wide <- barcode_summarized %>%
  pivot_wider(names_from = Clade, values_from = TotalCount, values_fill = 0)
shannon_diversity <- barcode_wide %>%
  select(-Island) %>% # Exclude Island column for diversity calculation
  apply(1, function(x) diversity(x, index = "shannon"))
barcode_wide <- barcode_wide %>%
  mutate(ShannonDiversity = shannon_diversity)
island_areas <- data.frame(
  Island = c("Maug", "Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam"),
  Area = c(2.1, 47.2, 4.9, 122, 102, 85.5, 541.3) # replace with actual area values
)
barcode_wide <- barcode_wide %>%
  left_join(island_areas, by = "Island")

# Plot it 
ggplot(barcode_wide, aes(x = Area, y = ShannonDiversity)) +
  geom_point(stat = "identity") +
  geom_smooth(method = "lm", color="dodgerblue")+
  labs(x = "Land Area (km )", y = "Shannon Diversity Index") +
  theme_minimal()+
  theme(axis.title.x = element_text(face = "bold", size =18),
        axis.title.y = element_text(face = "bold", size =18))

###### Does Diversity Correlate to distance from the equator? Linear Regression Model ######
# Data curation + Shannon Calculation

distance_equator <- data.frame(
  Island = c("Maug", "Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam"),
  Equator = c(2202.99717, 1991.86955, 1838.70896, 1648.869442, 1678.263884, 1559.219442, 1476.963884) # replace with actual area values
)
barcode_wide <- barcode_wide %>%
  left_join(distance_equator, by = "Island")
# Plot it 
ggplot(barcode_wide, aes(x = Equator, y = ShannonDiversity)) +
  geom_point(stat = "identity") +
  geom_smooth(method = "lm", color="dodgerblue")+
  labs(x = "Distance from Equator (km)", y = "Shannon Diversity Index") +
  theme_minimal()+
  theme(axis.title.x = element_text(face = "bold", size =18),
        axis.title.y = element_text(face = "bold", size =18))

###### Does Relative Abundance of Clade III Correlate to Lattitude? Linear Regression Model ######
relative_abundance <- data.frame(
  Island = c("Maug", "Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam"),
  RelativeAbundance = c(73.91, 47.62, 72.00, 9.09, 27.91, 38.64, 13.55)
)
barcode_wide <- barcode_wide %>%
  left_join(relative_abundance, by = "Island")

# Plot it 
ggplot(barcode_wide, aes(x = Equator, y = RelativeAbundance)) +
  geom_point(stat = "identity") +
  geom_smooth(method = "lm", color="dodgerblue")+
  labs(x = "Distance from Equator (km)", y = "Relative Abundance (clade III)") +
  theme_minimal()+
  theme(axis.title.x = element_text(face = "bold", size =18),
        axis.title.y = element_text(face = "bold", size =18))
model <- lm(RelativeAbundance ~ Equator, data = barcode_wide)
model_summary <- summary(model)

###### Does Diversity Correlate to Availible Reef Area? Linear Regression Model ######
# Data curation + Shannon Calculation
barcode_summarized <- barcode.long %>%
  group_by(Island, Clade) %>%
  summarize(TotalCount = sum(Count, na.rm = TRUE)) %>%
  ungroup()
barcode_wide <- barcode_summarized %>%
  pivot_wider(names_from = Clade, values_from = TotalCount, values_fill = 0)
shannon_diversity <- barcode_wide %>%
  select(-Island) %>% # Exclude Island column for diversity calculation
  apply(1, function(x) diversity(x, index = "shannon"))
barcode_wide <- barcode_wide %>%
  mutate(ShannonDiversity = shannon_diversity)
# Reef area data (Maug has no data, so it is excluded here)
island_reef_area <- data.frame(
  Island = c("Pagan", "Sarigan", "Saipan", "Tinian", "Rota", "Guam"),
  ReefArea = c(1.18, 0.03, 54.19, 17.08, 21.58, 88.87)
)

# Join reef-area data
barcode_wide <- barcode_wide %>%
  left_join(island_reef_area, by = "Island")
# Plot Reef Area vs Shannon Diversity
ggplot(barcode_wide, aes(x = ReefArea, y = ShannonDiversity)) +
  geom_point() +
  geom_smooth(method = "lm", color = "dodgerblue") +
  labs(
    x = "Reef Area (km )",
    y = "Shannon Diversity Index"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18)
  )
model <- lm(ShannonDiversity ~ ReefArea, data = barcode_wide)
model_summary <- summary(model)
model_summary
