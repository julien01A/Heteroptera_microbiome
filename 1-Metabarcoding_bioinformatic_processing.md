All metabarcoding bioinformatic analyses were performed using a step-by-step workflow with the FROGS pipeline implemented on Galaxy.

To begin, all demultiplexed raw-read files, available at `PASTE THE NCBI LINK`, were compiled in a single `ARCHIVE.tar` file.

The following dependencies and versions were used:
```
#### Galaxy / Job Dependencies ####
# frogs: version 4.10.0
# blast: version 2.10
# cutadapt: version 2.10
# emboss: version 6.6.0
# flash: version 1.2.11
# swarm: version 3.0.0
# vsearch: version 2.17.0
```

The workflow starts with a pre-processing, including read merging, denoising, and dereplication:
```
#### Galaxy / FROGS_1 Pre-process ####
# Tool Parameters
# Sequencer: illumina
# Input type: archive
# TAR archive file: ARCHIVE.tar
# Are reads already merged?: paired
# Reads 1 size: 305
# Reads 2 size: 305
# Mismatch rate: 0.1
# Merge software: flash
# Expected amplicon size: 430
# Would you like to keep unmerged reads?: No, unmerged reads will be excluded.
# Minimum amplicon size: 390
# Maximum amplicon size: 480
# Do the sequences have PCR primers?: true
# 5' primer: TACGGGAGGCAGCAG
# 3' primer: GGATTAGATACCCTGG
```

#########################################
#### Galaxy / FROGS_2 Clustering swarm [Single-linkage clustering on sequences] ####
# Tool Parameters
# Sequences file: FROGS_1 Pre-process: dereplicated.fasta
# Count file: FROGS_1 Pre-process: count.tsv
# FROGS guidelines version: 3.2
# 	Aggregation distance clustering: 1
# 	Refine clustering: Yes, refine clustering with --fastidious swarm option
#########################################

#########################################
#### Galaxy / FROGS_3 Clustering swarm [Single-linkage clustering on sequences] ####
# Tool Parameters
# Sequences file (format: FASTA): FROGS_2 Clustering swarm: seed_sequences.fasta
# Abundance type: biom
# Abundance file (format: BIOM): FROGS_2 Clustering swarm: clustering_abundance.biom
#########################################

#########################################
#### Galaxy / FROGS_Cluster_Stat [Process some metrics on clusters] ####
# Tool Parameters
# Abundance file: FROGS_2 Clustering swarm: clustering_abundance.biom
#########################################

#########################################
#### Galaxy / FROGS_4 Cluster filters [Filters clusters on several criteria] ####
# Tool Parameters
# Sequence file: FROGS_3 Remove Chimera: non_chimera.fasta
# Abundance file: FROGS3_ Remove chimera: non_chimera_abundance.biom
# Minimum prevalence method: all
# 	Minimum prevalence: Not available.
# Minimum cluster abundancy as proportion or count. We recommend to use a proportion of 0.00005.: count
# 	Minimum number of sequences to keep cluster: 2
# N biggest clusters: Not available.
# Search for contaminant clusters.: no
#########################################

#########################################
#### Galaxy / FROGS_5 Taxonomic affiliation [Taxonomic affiliation of each ASV's seed by RDPtools and BLAST] ####
# Tool Parameters
# reference database: silva_138_16S
# Also perform RDP assignation?: Yes
# Taxonomic ranks: Domain Phylum Class Order Family Genus Species
# Sequence file: FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Abundance file: FROGS_4 Cluster filters: clusterFilters_abundance.biom
#########################################

#########################################
#### Galaxy / FROGS BIOM to TSV [Converts a BIOM file in a TSV file 1] ####
# Tool Parameters
# Abundance file: FROGS_5 Taxonomic affiliation: affiliation_abundance.biom
# Sequences file (optional): FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Extract multi-alignments: Yes
#########################################

# Dowload and inspect the .tsv file. Rename the sample names with the final names. Then, reupload the new file called: abundance_cleaned.tsv

#########################################
#### Galaxy / FROGS TSV to BIOM [Converts a TSV file in a BIOM file 1] ####
# Tool Parameters
# Abundance TSV File: 5-abundance_cleand.tsv
# Multi_affiliation TSV File: Yes
#########################################

#########################################
#### Galaxy / FROGS Affiliation postprocess [Aggregates ASVs based on alignment metrics] ####
# Sequence file: FROGS Affiliation postprocess: sequences.fasta
# Abundance file: FROGS Affiliation postprocess: affiliation_abundance.biom
# Minimum prevalence method: all
# 	Minimum prevalence: Not available.
# Minimum cluster abundancy as proportion or count. We recommend to use a proportion of 0.00005.: proportion
# 	Minimum proportion of sequences abundancy to keep cluster: 5e-05
# N biggest clusters: Not available.
# Search for contaminant clusters.: no
#########################################

#########################################
#### Galaxy / FROGS Tree [Reconstruction of phylogenetic tree] ####
# Tool Parameters
# Sequence file: FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Biom file: FROGS_4 Cluster filters: clusterFilters_abundance.biom
#########################################

# Creation of a variable.tsv file containing the species and family name of each sample.

#########################################
#### Galaxy / FROGSSTAT Phyloseq Import Data [from 3 files: biomfile, samplefile, treefile] ####
# Tool Parameters
# Abundance biom file with taxonomical metadata (format: BIOM): FROGS_4 Cluster filters: clusterFilters_abundance.biom
# Metadata associated to samples (format: TSV): variable.tsv
# Taxonomic tree file (format: Newick): FROGS Tree: tree.nwk
# Names of taxonomic levels: Kingdom Phylum Class Order Family Genus Species
# Do you want to normalise your data ? No, keep abundance as it is.
#########################################

# One final file is obtained: 8-Stinkbugs_microbiome.rdata

# Additionnaly, we create a secundary abundance file grouping all the bacteria-identified species forming a same genus together to have a Genus-identification for each sample. This secundary abundance file was created in three steps. 
# Step 1. Create a file identifying each genus for each ASV. This file was called: 5-abundance_split_2.xlsx, and was obtained from abundance_cleaned.tsv using the following R script:
#########################################
#### R ####
library(tidyr)
library(dplyr)
tax <- read.delim("~/5-abundance_cleaned.txt", check.names = FALSE)
tax$Taxonomy_clean <- gsub("\\(.*?\\)", "", tax$rdp_tax_and_bootstrap)  # Remove the confidence scores
tax$Taxonomy_clean <- gsub(";+", ";", tax$Taxonomy_clean)  # Clean up repeated semicolons if any
tax$Taxonomy_clean <- gsub("^;|;$", "", tax$Taxonomy_clean) # Remove any leading/trailing semicolons
tax <- tax %>% select(-rdp_tax_and_bootstrap) # Remove the rdp_tax_and_bootstrap column
tax_split <- tax %>%  separate(Taxonomy_clean, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") # Split the Taxonomy_clean column to obtain ranks
write.table(tax_split, "~/5-abundance_split.xlsx", col=NA, sep="\t",dec=".")
#########################################

# Step 2. Grouping all the ASV of the same genus together to have a Genus-identification for each sample. To do it, we used the file 5-abundance_split_2.xlsx and used the following R script :
#########################################
#### R ####
# library(dplyr)
# library(readxl)
# library(openxlsx)
# input_file <- "~/5-abundance_split.xlsx"
output_file <- "~/5-abundance_aggregate.xlsx"
tax_table <- read_excel(input_file)

# Detect the numeric columns to sum:
tax_cols <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
numeric_cols <- setdiff(names(tax_table), tax_cols)
tax_table[numeric_cols] <- lapply(tax_table[numeric_cols], function(x) as.numeric(x))
aggregated <- tax_table %>%
  group_by(Genus) %>%
  summarise(
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Species = first(Species),
    across(all_of(numeric_cols), sum, na.rm = TRUE)
  )
write.xlsx(aggregated, output_file) #export the new .xslx file
#########################################

# Step3. Manually inspect the 5-abundance_aggregate.xlsx file obtained. Remove the column 1... if present and do the necessary changed to feat with the most appropriate file for statistic analyses (eg. sample in lines, Bacterial genera in columns).
###############################################################################################################################################################################################################################################################################################
