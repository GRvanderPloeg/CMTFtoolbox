library(tidyverse)
library(ggplot2)
library(stringr)
library(parafac4microbiome)
library(CMTFtoolbox)

# Faecal microbiome
df = read.csv("./data-raw/faecesCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/newTaxonomy_faeces.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./data-raw/faeces_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

# Read infant anthopometrics
infant_anthropometrics = read.csv("./data-raw/infant_anthropometrics.csv") %>% as_tibble()

# Make subject metadata
sampleInfo = sampleInfo %>% left_join(infant_anthropometrics %>% select(RCID, whz.6m))
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis, whz.6m) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# Filter based on sparsity
threshold = 0.75

sparsity = colSums(df==0) / nrow(df)

featureSelection = sparsity <= threshold
taxonomy_filtered = taxonomy[featureSelection,]
df_filtered = df[,featureSelection]

# CLR
df_clr = t(apply(df_filtered+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_clr)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_clr, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_clr),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
}

# As in Jakobsen2025, remove subject 282 (row 90)
X = X[-90,,]
subjectMeta = subjectMeta[-90,,]

# Center and scale
X_cnt = parafac4microbiome::multiwayCenter(X, mode=1)
X_cnt_scl = parafac4microbiome::multiwayScale(X_cnt, mode=2)

faeces_df = X_cnt_scl
faeces_subjectMeta = subjectMeta
faeces_taxonomy = taxonomy_filtered
faeces_timepoints = timepoints

# Milk microbiome
df = read.csv("./data-raw/milkCounts.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/newTaxonomy_milk.csv", header=FALSE, sep=" ") %>% as_tibble()
sampleInfo = read.csv("./data-raw/milk_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("Sample", "RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Filter taxa to taxa with at least 1 non-zero value
featureMask = colSums(df) > 0
df = df[,featureMask]
taxonomy = taxonomy[featureMask,]

# Filter based on sparsity
threshold = 0.85

sparsity = colSums(df==0) / nrow(df)

featureSelection = sparsity <= threshold
taxonomy_filtered = taxonomy[featureSelection,]
df_filtered = df[,featureSelection]

# CLR
df_clr = t(apply(df_filtered+1, 1, function(x){log(x / compositions::geometricmean(x))})) %>% as_tibble()

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_clr)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_clr, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_clr),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
}

# Center and scale
X_cnt = parafac4microbiome::multiwayCenter(X, mode=1)
X_cnt_scl = parafac4microbiome::multiwayScale(X_cnt, mode=2)

milk_df = X_cnt_scl
milk_subjectMeta = subjectMeta
milk_taxonomy = taxonomy_filtered
milk_timepoints = timepoints

# Milk metabolomics
df = read.csv("./data-raw/milkMetabNumeric.csv", header=FALSE, sep=" ") %>% as_tibble()
taxonomy = read.csv("./data-raw/milk_metab_CAS_numbers.csv", header=TRUE, sep=",") %>% as_tibble()
sampleInfo = read.csv("./data-raw/milkMetab_sampleMeta.csv", header=FALSE, sep=" ") %>% as_tibble()
colnames(sampleInfo) = c("RCID", "BMI", "BMI.group", "Days", "Gestational.age", "C.section", "AB.infant", "AB.mother", "Secretor", "Lewis", "subject")

# Make subject metadata
subjectMeta = sampleInfo %>% select(subject, BMI, BMI.group, C.section, Secretor, Lewis) %>% unique() %>% mutate(subject=as.character(subject)) %>% arrange(subject)

# Remove duplicates
drop = c(61,81,146)
subjectMeta = subjectMeta %>% mutate(index=1:nrow(.)) %>% filter(!index %in% drop) %>% select(-index) %>% arrange(subject)

# Log transform
df_log = log(df)

# Make into cube
I = length(unique(sampleInfo$subject))
J = ncol(df_log)
K = length(unique(sampleInfo$Days))
X = array(0L, c(I,J,K))
timepoints = sampleInfo %>% arrange(Days) %>% select(Days) %>% unique() %>% pull()

for(k in 1:K){
  Day = timepoints[k]
  X[,,k] = cbind(df_log, sampleInfo) %>% as_tibble() %>% mutate(subject=as.character(subject)) %>% filter(Days == Day) %>% select(c(colnames(df_log),subject)) %>% right_join(subjectMeta) %>% arrange(subject) %>% select(-colnames(subjectMeta)) %>% as.matrix()
}

# Center and scale
X_cnt = parafac4microbiome::multiwayCenter(X, mode=1)
X_cnt_scl = parafac4microbiome::multiwayScale(X_cnt, mode=2)

milkMetab_df = X_cnt_scl
milkMetab_subjectMeta = subjectMeta
milkMetab_taxonomy = taxonomy
milkMetab_timepoints = timepoints

sharedSubjects = intersect(intersect(faeces_subjectMeta$subject, milk_subjectMeta$subject), milkMetab_subjectMeta$subject)

faeces_homogenized = faeces_df[faeces_subjectMeta$subject %in% sharedSubjects,,]
milk_homogenized = milk_df[milk_subjectMeta$subject %in% sharedSubjects,,]
milkMetab_homogenized = milkMetab_df[milkMetab_subjectMeta$subject %in% sharedSubjects,,]

homogenized_subjectMeta = faeces_subjectMeta %>% filter(subject %in% sharedSubjects) %>% arrange(subject)

datasets = list(faeces_homogenized, milk_homogenized, milkMetab_homogenized)
modes = list(c(1,2,3),c(1,4,5),c(1,6,7))
Z = setupCMTFdata(datasets, modes)


milkMetab_featureMeta = milkMetab_taxonomy %>% mutate(X = make.names(X), CAS.Registry = make.names(CAS.Registry))
milkMetab_featureMeta[70,1] = "tau.Methylhistidine" # Fix non-ascii tau character

Jakobsen2025 = list("Z"=Z,
                    "homogenizedSubjectMetadata"=homogenized_subjectMeta,
                    "faecal_microbiome_taxonomy"=faeces_taxonomy,
                    "milk_microbiome_taxonomy"=milk_taxonomy,
                    "milk_metabolites"=milkMetab_featureMeta)

usethis::use_data(Jakobsen2025, overwrite = TRUE)
