## code to prepare `Georgiou2025` dataset goes here

library(phyloseq)
library(tidyverse)
library(vegan)
library(ape)
library(ggpubr)
library(ggplot2)
library(parafac4microbiome)
library(CMTFtoolbox)

# Load data
cytokines_data = read.csv("./data-raw/input_deduplicated_RvdP.csv", sep=" ", header=FALSE) %>% as_tibble()
colnames(cytokines_data) = c("VEGF", "CRP", "GM-CSF", "IL1alpha", "IL1beta", "IL4", "IL6", "IL8", "IL10", "IL12p70", "IL17A", "IFNgamma", "MIP1alpha", "OPG", "TNFalpha", "RANKL")

cytokines_meta_data = read.csv("./data-raw/input_deduplicated_metadata_RvdP.csv", sep=" ", header=FALSE) %>% as_tibble()
colnames(cytokines_meta_data) = c("SubjectID", "Visit", "Gender", "Age", "Pain_noPain", "case_control")

otherMeta = read.csv("./data-raw/Root_meta_data_parafac.txt", sep="\t") %>% as_tibble() %>% select(-Gender)

# Put into cube
cytokineData = parafac4microbiome::reshapeData(cytokines_data,
                                subjectMetadata = cytokines_meta_data$SubjectID,
                                featureMetadata = colnames(cytokines_data),
                                timepointMetadata = cytokines_meta_data$Visit)

# Prepare export
cytokineData$mode1 = cytokineData$mode1 %>%
  mutate(SubjectID = subjectMetadata) %>%
  left_join(cytokines_meta_data %>% select(-Visit) %>% unique()) %>%
  select(SubjectID, Gender, Age, case_control) %>%
  left_join(otherMeta %>% select(SubjectID, PainS_NopainA) %>% unique())

cytokineData$mode2 = cytokineData$mode2 %>% mutate(name = V1) %>% select(name)

cytokineData$mode3 = cytokineData$mode3 %>% mutate(Visit = timepointMetadata, Week = c(-6, -3, 0, 1, 6, 13))

# Load microbiome data - these are case only by default
microbiome_raw = read.csv("./data-raw/20240429_microbiome_counts.csv", sep=" ", header=FALSE)
taxonomy = read.csv("./data-raw/20240429_taxonomy.csv", sep=" ", header=FALSE)
colnames(taxonomy) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "zOTU")
subjectMeta2 = read.csv("./data-raw/20240429_microbiome_sampleMeta.csv", sep=" ", header=FALSE)

# # Remove extra samples
# remove = c("A11-18", "A11-3", "A11-8 36", "A11-10 17", "A11-15 17") # last is new: "A11-8 46"
# mask = !(subjectMeta2[,3] %in% remove)
#
# microbiome_raw = microbiome_raw[mask,]
# subjectMeta2 = subjectMeta2[mask,]

# Prepare export
Tooth_microbiome = list()

Tooth_microbiome$mode1 = subjectMeta2 %>%
  as_tibble() %>%
  mutate(SubjectID = V2, SampleID = V3) %>%
  select(SubjectID, SampleID) %>%
  left_join(cytokineData$mode1)

Tooth_microbiome$mode2 = taxonomy %>% as_tibble()
Tooth_microbiome$data = microbiome_raw

# Prepare export
Georgiou2025 = list()
Georgiou2025$Inflammatory_mediators = cytokineData
Georgiou2025$Tooth_microbiome = Tooth_microbiome

usethis::use_data(Georgiou2025, overwrite = TRUE)
