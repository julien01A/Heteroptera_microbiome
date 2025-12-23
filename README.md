# **Stinkbugs microbiome**

## 0. Quality controle of raw-read files

All the raw-read sequencing files were first checked for their quality using three tools: FastQC (v.0.12.1) (<https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>), Nanoplot (v.1.43.0) (<https://github.com/wdecoster/NanoPlot>, <https://doi.org/10.1093/bioinformatics/btad311>) and Nanocomp (v.1.23.1) (<https://github.com/wdecoster/nanocomp>).

For FastQC, we used:
```
### bash ####
for file in *.fastq; do
    fastqc "$file"
done
```

Then, for Nanoplot:
```
mkdir -p Quality_reads
for reads_file in *.fastq.gz; do
    reads_ech=$(basename "$reads_file" .fastq.gz)
    NanoPlot -t 2 --fastq "$reads_file" -o "./Quality_reads/Quality_${reads_ech}" --N50
done
```

Finally, for NanoComp:
```
mkdir -p Quality_NanoComp_all_samples
fastq_files=$(ls *.fastq.gz)
NanoComp -t 2 --fastq $fastq_files -o Quality_NanoComp_all_samples
```

## 1. Bioinformatic processing and taxonomic assignment of 16S metabarcoding raw-read files

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

The workflow starts with a pre-processing called `FROGS_1 Pre-process`, including read merging, denoising, and dereplication:
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

The workflow continues with a clustering step called `FROGS_2 Clustering swarm`, which performs single-linkage clustering on sequences:
```
#### Galaxy / FROGS_2 Clustering swarm ####
# Tool Parameters
# Sequences file: FROGS_1 Pre-process: dereplicated.fasta
# Count file: FROGS_1 Pre-process: count.tsv
# FROGS guidelines version: 3.2
# 	Aggregation distance clustering: 1
# 	Refine clustering: Yes, refine clustering with --fastidious swarm option
```

The workflow then proceeds with a clustering step called `FROGS_3 Remove chimera`, which removes PCR chimera in each sample:
```
#### FROGS_3 Remove chimera ####
# Tool Parameters
# Sequences file (format: FASTA): FROGS_2 Clustering swarm: seed_sequences.fasta
# Abundance type: biom
# Abundance file (format: BIOM): FROGS_2 Clustering swarm: clustering_abundance.biom
```

The workflow then includes a statistical analysis step called `FROGS_Cluster_Stat`, which processes metrics on clusters:
```
#### Galaxy / FROGS_Cluster_Stat [Process some metrics on clusters] ####
# Tool Parameters
# Abundance file: FROGS_2 Clustering swarm: clustering_abundance.biom
```

Next, the workflow applies a filtering step called `FROGS_4 Cluster filters`, where clusters are filtered based on several criteria:
```
#### Galaxy / FROGS_4 Cluster filters ####
# Tool Parameters
# Sequence file: FROGS_3 Remove Chimera: non_chimera.fasta
# Abundance file: FROGS3_ Remove chimera: non_chimera_abundance.biom
# Minimum prevalence method: all
# 	Minimum prevalence: Not available.
# Minimum cluster abundancy as proportion or count. We recommend to use a proportion of 0.00005.: count
# 	Minimum number of sequences to keep cluster: 2
# N biggest clusters: Not available.
# Search for contaminant clusters.: no
```

The workflow then performs a taxonomic assignment step called `FROGS_5 Taxonomic affiliation`, in which each ASV seed is taxonomically affiliated using RDPtools and BLAST:
```
#### Galaxy / FROGS_5 Taxonomic affiliation ####
# Tool Parameters
# reference database: silva_138_16S
# Also perform RDP assignation?: Yes
# Taxonomic ranks: Domain Phylum Class Order Family Genus Species
# Sequence file: FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Abundance file: FROGS_4 Cluster filters: clusterFilters_abundance.biom
```

The workflow includes a conversion step called `FROGS BIOM to TSV`, which converts the BIOM file into a TSV file:
```
#### Galaxy / FROGS BIOM to TSV [Converts a BIOM file in a TSV file 1] ####
# Tool Parameters
# Abundance file: FROGS_5 Taxonomic affiliation: affiliation_abundance.biom
# Sequences file (optional): FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Extract multi-alignments: Yes
```

Dowload and inspect the `.tsv` file. Rename the sample names with the final names. Then, reupload the modified file called: `abundance_cleaned.tsv`

Then, upload the cleaned TSV file `abundance_cleaned.tsv` and convert it back into a BIOM file using FROGS TSV to BIOM:
```
#### Galaxy / FROGS TSV to BIOM ####
# Tool Parameters
# Abundance TSV File: 5-abundance_cleand.tsv
# Multi_affiliation TSV File: Yes
```

The workflow then includes a post-processing step called `FROGS Affiliation postprocess`, which aggregates ASVs based on alignment metrics:
```
#### Galaxy / FROGS Affiliation postprocess ####
Sequence file: FROGS_4 Cluster filters: clusterFilters_sequences.fasta
Abundance file: FROGS_5 Taxonomic affiliation: affiliation_abundance.biom
Is this an amplicon hyper variable in length? No
Minimum identity for aggregation: 99.0
Minimum coverage for aggregation: 96.0
```

The workflow applies a secundary filtering step also called `FROGS_4 Cluster filters`:
```
#### Galaxy / FROGS_4 Cluster filters ####
# Sequence file: FROGS Affiliation postprocess: sequences.fasta
# Abundance file: FROGS Affiliation postprocess: affiliation_abundance.biom
# Minimum prevalence method: all
# 	Minimum prevalence: Not available.
# Minimum cluster abundancy as proportion or count. We recommend to use a proportion of 0.00005.: proportion
# 	Minimum proportion of sequences abundancy to keep cluster: 5e-05
# N biggest clusters: Not available.
# Search for contaminant clusters.: no
```

The workflow then performs a phylogenetic analysis step called `FROGS Tree`, which reconstructs the phylogenetic tree based on 16S sequences of ASVs:
```
#### Galaxy / FROGS Tree ####
# Tool Parameters
# Sequence file: FROGS_4 Cluster filters: clusterFilters_sequences.fasta
# Biom file: FROGS_4 Cluster filters: clusterFilters_abundance.biom
```

Then, we created a `variable.tsv` file containing the `species` and `family` names for each sample and we uploaded it on Galaxy.

Finally, to create one final microbiome data file `8-Stinkbugs_microbiome.rdata` available for numeric analysis, we used 3 files: `biomfile`, `samplefile`, `treefile` using a final step called `FROGSSTAT Phyloseq Import Data` as follow:
```
#### Galaxy / FROGSSTAT Phyloseq Import Data ####
# Tool Parameters
# Abundance biom file with taxonomical metadata (format: BIOM): FROGS_4 Cluster filters: clusterFilters_abundance.biom
# Metadata associated to samples (format: TSV): variable.tsv
# Taxonomic tree file (format: Newick): FROGS Tree: tree.nwk
# Names of taxonomic levels: Kingdom Phylum Class Order Family Genus Species
# Do you want to normalise your data ? No, keep abundance as it is.
```

Additionally, we created a secondary microbiome data file by grouping all bacteria-identified ASV species belonging to the same genus, in order to obtain a bacterial genus-level identification for each sample. This secondary microbiome data file was generated in three steps.

**Step 1.** Create a file identifying each genus for each ASV. This file was called: `5-abundance_split_2.xlsx`, and was obtained from `abundance_cleaned.tsv` using the following R script:
```
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
```

**Step 2.** Group all the ASV of the same genus together to have a Genus-identification for each sample. To do it, we used the file `5-abundance_split_2.xlsx` and used the following R script :
```
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
```

**Step3.** Manually inspect the `5-abundance_aggregate.xlsx` file obtained. Remove the `column 1...` if present and do the necessary changed to feat with the most appropriate file for statistic analyses (eg. sample in lines, Bacterial genera in columns), etc. 

Our two final files `8-Stinkbugs_microbiome.rdata` and `5-Stinkbugs-abundance.xlsx` were available in this GitHub page.

## 2. Microbiome graphical representation and analyses

The R scripts used to analyse and represent the microbiome of stinkbugs are available above.

#### Heatmap of the 15 most abundant bacterial genera 

Here, we aimed to represent the most abundant bacterial genera in terms of the relative number of reads per species. To define the relative number of reads for each sample, we first created a new column in the `5-Stinkbugs-abundance.xlsx` file called `Total_nb_reads`, corresponding to the sum of all bacterial reads for each sample. We then used this column to calculate the relative number of reads per sample for each bacterial genus. The script is as follows:

```
#### R ####
# load the libraries
library(tidyverse)
library(pheatmap)
# load the data file
data <- read.delim("5-Stinkbugs-abundance.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Define columns of metadata 
meta_cols <- c("sample_full_name", "sample", "species", "family","Total_nb_reads")
bact_cols <- setdiff(colnames(data), meta_cols)
# Relative proportion of reads by samples
data_rel <- data %>%
  mutate(across(all_of(bact_cols),~ .x / Total_nb_reads))
# Find the top14 bacterial genera
top14_genera <- data_rel %>%
  summarise(across(all_of(bact_cols), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "genus", values_to = "abundance") %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 14) %>%
  pull(genus)
# Create a column "other" for all the other genera
data_rel <- data_rel %>%
  mutate(Other = rowSums(across(setdiff(bact_cols, top14_genera))))
# Select the Top14 and other columns
microbiome_15 <- data_rel %>%
  select(species, all_of(top14_genera), Other)
#Agregate data by stinkbugs species
microbiome_species <- microbiome_15 %>%
  group_by(species) %>%
  summarise(across(all_of(c(top14_genera, "Other")), mean, na.rm = TRUE))
#Define a species order
species_order <- c("Picromerus_bidens","Graphosoma_italicum","Graphosoma_semipunctatum","Eurydema_oleracea","Eurydema_ornata","Nezara_viridula","Acrosternum_hegeeri","Palomena_prasina","Carpocoris_sp","Dolycoris_baccarum","Rhaphigaster_nebulosa","Halyomorpha_halys","Pentatomid_sp","Pyrrhocoris_apterus","Scantius_aegyptius","Syromastus_rhombeus","Lygaeus_equestris")
# Define the order of columns 
col_order <- microbiome_species %>%
  select(-species) %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "genus", values_to = "abundance") %>%
  arrange(desc(abundance)) %>%
  pull(genus)
col_order <- c(setdiff(col_order, "Other"), "Other")
# Creation of the heatmap
heatmap_matrix <- microbiome_species %>%
  filter(species %in% species_order) %>%
  mutate(species = factor(species, levels = species_order)) %>%
  arrange(species) %>%
  select(species, all_of(col_order)) %>%
  column_to_rownames("species") %>%
  as.matrix()
#plot
p <- pheatmap(heatmap_matrix, scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("#666666", "#FFCC00", "red"))(100), fontsize_row = 10, fontsize_col = 10, border_color = NA, angle_col = "90")
p
```

The raw figure generated in R is available on this GitHub page under the name `Stinkbugs_heatmap_raw_fig.png`. The figure was then manually post-processed for graphical adaptations.
