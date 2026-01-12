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
### bash ####
mkdir -p Quality_reads
for reads_file in *.fastq.gz; do
    reads_ech=$(basename "$reads_file" .fastq.gz)
    NanoPlot -t 2 --fastq "$reads_file" -o "./Quality_reads/Quality_${reads_ech}" --N50
done
```

Finally, for NanoComp:
```
### bash ####
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

### **Heatmap of the 15 most abundant bacterial genera**

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

### **Composition plot of microbiome by samples**

Here, we aimed to represent the most abundant bacterial genera in terms of relative abundance (%) per samples. The script is as follows:

```
#### R ####
# load the libraries
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
# load the data file
data <- read.delim("5-Stinkbugs-abundance.txt", header = TRUE, sep = "\t", check.names = FALSE)
meta_cols <- c("sample_full_name", "sample", "samplebis", "species", "family", "Total_nb_reads")
bact_cols <- setdiff(colnames(data), meta_cols)
# Relative proportion of reads by samples
data_rel <- data %>%
  mutate(across(all_of(bact_cols), ~ .x / Total_nb_reads))
# Find the top20 bacterial genera
top20_genera <- data_rel %>%
  summarise(across(all_of(bact_cols), sum, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "genus", values_to = "abundance") %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 20) %>%
  pull(genus)
# Create a column "other" for all the other genera
data_rel <- data_rel %>%
  mutate(
    Other = rowSums(across(setdiff(bact_cols, top20_genera)))
  )
# Select the Top20 and other columns
composition_df <- data_rel %>%
  select(sample, species, all_of(top20_genera), Other)

# Reshape to long format (required for ggplot)
composition_long <- composition_df %>%
  pivot_longer(
    cols = c(all_of(top20_genera), Other),
    names_to = "genus",
    values_to = "relative_abundance"
  ) %>%
  mutate(relative_abundance = relative_abundance * 100)  # %

# Define species order
species_order <- c("Picromerus_bidens","Graphosoma_italicum","Graphosoma_semipunctatum","Eurydema_oleracea","Eurydema_ornata","Nezara_viridula","Acrosternum_hegeeri",
  "Palomena_prasina","Carpocoris_sp","Dolycoris_baccarum","Rhaphigaster_nebulosa","Halyomorpha_halys","Pentatomid_sp","Pyrrhocoris_apterus","Scantius_aegyptius",
  "Syromastus_rhombeus","Lygaeus_equestris")
sample_order <- composition_long %>%
  distinct(sample, species) %>%
  mutate(species = factor(species, levels = rev(species_order))) %>%
  arrange(species, sample) %>%
  pull(sample)
composition_long <- composition_long %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    genus  = factor(genus, levels = c(setdiff(top20_genera, "Other"), "Other"))  )

#Set colors
n_genus <- length(levels(composition_long$genus))
palette_genus <- c("#b2df8a", "#d95f02", "#7570b3", "#e7298a","#66a61e", "#e6ab02", "#a6761d", "aquamarine","#1f78b4", "#1b9e77", "#fb9a99", "#fdbf6f","#cab2d6", "azure4", "#b15928","#8dd3c7", "#ffffb3", "yellow2","azure2", "#80b1d3", "chocolate4")[1:n_genus]
names(palette_genus) <- levels(composition_long$genus)
palette_genus["Other"] <- "black"

#composition plot
ggplot(composition_long, aes(x = relative_abundance, y = sample, fill = genus)) +
  geom_col(
    width = 0.8,
    position = position_stack(reverse = TRUE)
  ) +
  scale_fill_manual(values = palette_genus) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Relative abundance (%)",
    y = "Sample",
    fill = "Bacterial genus",
    title = "Microbiome composition of stinkbug samples"
  )
```
The raw figure generated in R is available on this GitHub page under the name `Stinkbugs_compositionplot_raw_fig.png`.

### **Chao1 index by samples**

Here, we estimated the Chao1 index for each sample as a proxy of alpha diversity. The Chao1 index is an estimate of the total bacterial genus richness in a sample, accounting for rare or unobserved genera to better approximate the actual number of bacterial genera present. The script is as follows:

```
#### R ####
# load the libraries
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)
# load the data file
data <- read.delim("5-Stinkbugs-abundance.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Define columns of metadata 
meta_cols <- c("sample_full_name", "sample", "samplebis", "species","family", "Total_nb_reads")
bact_cols <- setdiff(colnames(data), meta_cols)
# Calculate Chao1 for each sample
abundance_matrix <- data %>%
  select(sample, all_of(bact_cols)) %>%
  column_to_rownames("sample") %>%
  as.matrix()
chao1 <- estimateR(abundance_matrix)
# Put Chao1 indexes in a dataframe and add it to the metadata
alpha_div <- data.frame(
  sample = colnames(chao1),
  Chao1  = as.numeric(chao1["S.chao1", ]))
alpha_div <- alpha_div %>%
  left_join(
    data %>%
      select(sample, species, family) %>%
      distinct(),by = "sample")
# order by samples
sample_order <- c("42_Picromerus_bidens","33_Graphosoma_italicum","32_Graphosoma_italicum","31_Graphosoma_italicum","30_Graphosoma_italicum",
  "29_Graphosoma_italicum","28_Graphosoma_italicum","27_Graphosoma_italicum","26_Graphosoma_italicum","21_Graphosoma_semipunctatum",
  "20_Graphosoma_semipunctatum","8_Eurydema_oleracea","7_Eurydema_ornata","6_Eurydema_ornata","5_Eurydema_ornata","25_Nezara_viridula",
  "24_Nezara_viridula","23_Nezara_viridula","4_Acrosternum_hegeeri","22_Palomena_prasina","3_Carpocoris_sp","2_Carpocoris_sp","1_Carpocoris_sp",
  "9_Dolycoris_baccarum","43_Dolycoris_baccarum","11_Dolycoris_baccarum","10_Dolycoris_baccarum","45_Rhaphigaster_nebulosa","36_Halyomorpha_halys",
  "35_Halyomorpha_halys","34_Halyomorpha_halys","44_Pentatomid_sp","19_Pyrrhocoris_apterus","18_Pyrrhocoris_apterus","17_Pyrrhocoris_apterus",
  "16_Pyrrhocoris_apterus","14_Scantius_aegyptius","13_Scantius_aegyptius","12_Scantius_aegyptius","15_Syromastus_rhombeus","41_Lygaeus_equestris",
  "40_Lygaeus_equestris","39_Lygaeus_equestris","38_Lygaeus_equestris","37_Lygaeus_equestris")
sample_order <- composition_long %>%
  distinct(sample, species) %>%
  mutate(species = factor(species, levels = rev(species_order))) %>%
  arrange(species, sample) %>%
  pull(sample)
alpha_div <- alpha_div %>%
  mutate(sample = factor(sample, levels = sample_order)) %>%
  arrange(sample)

# connect each dots in the plot
segments_df <- alpha_div %>%
  mutate(
    sample_next = lead(sample),
    Chao1_next = lead(Chao1)
  ) %>%
  filter(!is.na(sample_next))
#plot
ggplot(alpha_div, aes(y = sample, x = Chao1)) +
  geom_segment(data = segments_df, aes(x = Chao1, xend = Chao1_next, y = sample, yend = sample_next), color = "black") +
  geom_point(size = 3, color = "black") +
  scale_x_continuous(breaks = seq(0, 150, by = 10),labels = c("0", "", "", "", "", "50", "", "", "", "", "100", "", "", "", "", "150"), limits = c(0, 150)) +
  geom_segment(aes(x = 0, xend = 150, y = -0.5, yend = -0.5), color = "black", size = 1) +
  theme_minimal(base_size = 12) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_line(color = "black"), panel.grid.major.y = element_blank(),panel.grid.minor = element_blank()) +
  labs(y = "Sample", x = "Chao1 index")
```

The raw figure generated in R is available on this GitHub page under the name `Stinkbugs_chao1_raw_fig.png`.

We then used a Wilcoxon rank sum test with continuity correction to statistically compare Chao 1 indices depending of the proportion of reads of the two most dominant bacteria Pantoea/Erwinia and Serratia. Here is the script used for Pantoea-Erwinia:

```
#### R ####
# Wilcoxon test unilateral
# calcul the proportion of reads matching with Pantoea-Erwinia and create a new column for samples =>30% of reads vs <30%
pantoea_prop <- data %>%
  select(sample, `Pantoea-Erwinia`, Total_nb_reads) %>%
  mutate(
    prop_pantoea = `Pantoea-Erwinia` / Total_nb_reads,
    group_pantoea = ifelse(prop_pantoea >= 0.3, ">=30%", "<30%")
  ) %>%
  select(sample, prop_pantoea, group_pantoea)
alpha_div2 <- alpha_div %>%
  left_join(pantoea_prop, by = "sample") %>%
  mutate(group_pantoea = factor(group_pantoea, levels = c("<30%", ">=30%")))
table(alpha_div2$group_pantoea)
# test
wilcox.test(Chao1 ~ group_pantoea, data = alpha_div2, alternative = "greater")
# plot
boxplot(Chao1 ~ group_pantoea, data = alpha_div2)
```

### **PCoA**

Here, we calculated the Bray-Curtis and Jaccard distances as proxies for beta diversity. The Bray–Curtis distance measures dissimilarity between two samples based on the relative abundances of bacterial genera, ranging from 0 (identical) to 1 (completely different), whereas the Jaccard distance is based on presence–absence data. We then performed a PCoA on this distance matrices, which allows us to reduce the complex multivariate relationships into a few axes that capture the main patterns of variation between samples, making it easier to visualize similarities and differences in microbiome composition. The script is as follows:

```
#### R ####
# load the libraries
library(tidyverse)
library(vegan)
library(ggplot2)
# load the data file
data <- read.delim("5-Stinkbugs-abundance.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Define columns of metadata 
meta_cols <- c("sample_full_name", "sample", "samplebis", "species", "family", "Total_nb_reads")
bact_cols <- setdiff(colnames(data), meta_cols)
# Create the abundance matrix by normalizing
data$sample <- make.unique(as.character(data$sample))
data_rel <- data %>%
  mutate(across(all_of(bact_cols), ~ .x / Total_nb_reads)) %>%
  column_to_rownames("sample")
abundance_matrix <- as.matrix(data_rel[, bact_cols])
# Calculate the Bray-Curtis index
bray_dist <- vegdist(abundance_matrix, method = "bray") # Here, change "bray" by "jaccard" or other distance indexes
# PCoA
pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
# Extract the coordinates and add metadata
scores_pcoa <- as.data.frame(pcoa_res$points)
colnames(scores_pcoa) <- c("PCoA1", "PCoA2")
scores_pcoa$sample <- rownames(scores_pcoa)
scores_pcoa <- scores_pcoa %>%
  left_join(data %>% select(sample, species), by = "sample")
# Add the metadata column "family"
scores_pcoa <- scores_pcoa %>%
  left_join(data %>% select(sample, family), by = "sample") %>%
  mutate(family = factor(family))

# Calculate the proportion of axe1 and axe2
eig <- pcoa_res$eig
prop_var <- eig / sum(eig)
prop_var_PC1 <- round(prop_var[1] * 100, 1)
prop_var_PC2 <- round(prop_var[2] * 100, 1)

# Define color palette for species and shapes for family
palette_species <- c("Picromerus_bidens" = "#1f78b4","Graphosoma_italicum" = "#33a02c","Graphosoma_semipunctatum" = "#e31a1c","Eurydema_oleracea" = "#ff7f00",
 "Eurydema_ornata" = "#6a3d9a", "Nezara_viridula" = "#b15928","Acrosternum_hegeeri" = "#a6cee3","Palomena_prasina" = "#b2df8a","Carpocoris_sp" = "yellow2",
 "Dolycoris_baccarum" = "#fdbf6f","Rhaphigaster_nebulosa" = "#cab2d6","Halyomorpha_halys" = "black","Pentatomid_sp" = "#8dd3c7","Pyrrhocoris_apterus" = "#80b1d3",
 "Scantius_aegyptius" = "pink2","Syromastus_rhombeus" = "#d9d9d9","Lygaeus_equestris" = "#bc80bd")
shapes_family <- c("Pentatomidea" = 15,"Lygaeidae" = 18,"Coreidae" = 17,"Pyrrhocoridae"= 19)

# PCoA plot
ggplot(scores_pcoa, aes(x = PCoA1, y = PCoA2, color = species, shape = family)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey", size = 0.2) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey", size = 0.2) +
  geom_point(size = 7, alpha = 0.9) +
  scale_color_manual(values = palette_species) +
  scale_shape_manual(values = shapes_family) +
  stat_ellipse(aes(group = family), type = "t", linetype = 1, color = "black", alpha = 0.5, size = 0.4) +
  theme_minimal(base_size = 14) +
  labs(
    x = paste0("PCoA1 (", prop_var_PC1, "%)"),
    y = paste0("PCoA2 (", prop_var_PC2, "%)"),
    color = "Species",
    shape = "Family") +
  theme(panel.grid = element_blank(),legend.position = "right")
```

The raw figures generated in R are available on this GitHub page under the names `Stinkbugs_PCoA_raw_fig_BRAY.png` and `Stinkbugs_PCoA_raw_fig_JACCARD.png`.

### **Phylosymbiosis**

To assess whether a phylosymbiosis pattern structures the microbiome of stinkbugs, we first constructed a distance matrix for the stinkug species. To do this, we collected COI sequences from 15 species (see the supplementary table online), aligned them using Clustal Omega (v.1.2.2) (<https://doi.org/10.1002/pro.3290>) implemented in the Unipro UGENE software (v.52.0) (<https://ugene.net/>, <https://doi.org/10.1093/bioinformatics/bts091>), and obtained a final alignment of 211 amino acids by removing gaps '-'. The COI sequences have been deposited in this GitHub repository under the name `coi_aa.fasta`. Using UGENE, we calculated a multiple sequence alignment distance matrix (in %) based on Hamming dissimilarity. This matrix is available as `matrix_aa_distance.txt`.

We then selected 15 representative specimens to compare their microbiota with the phylogenetic alignment, i.e., the multiple sequence alignment distance matrix. One specimen per species was chosen; for species with multiple specimens, we selected the individual with the highest Chao1 index. The microbiota matrix used was `abundance_phylosymbiosis.txt`, also deposited in the GitHub repository. Bray–Curtis and Jaccard distances were used to quantify microbiota dissimilarities. Finally, a Mantel test (spearman's correlation, 9,999 permutations) was performed to statistically estimate the correlation between the phylogenetic distance matrix and the microbiota dissimilarity matrix. The script used is as follows:

```
#### R ####
# load the libraries
library(vegan)
library(ggplot2)
# Load the matrix of phylogenetic distance
host_dist <- read.table("matrix_aa_distance.txt", header = TRUE, row.names = 1, sep = "\t", dec = ",", check.names = FALSE)
host_dist <- as.dist(host_dist)
# Load the representative microbiome matrix
microbiome <- read.table("abundance_phylosymbiosis.txt", header = TRUE, sep = "\t", check.names = FALSE)
rownames(microbiome) <- microbiome$name
microbiome <- microbiome[, colnames(microbiome) != "" & !is.na(colnames(microbiome))]
microbiome_mat <- microbiome %>%
  select(-(name:family))
#Bray
micro_dist <- vegdist(microbiome_mat, method = "bray") # change here with "jaccard" distance
#Mantel
common_taxa <- intersect(labels(host_dist), labels(micro_dist))
host_dist <- as.dist(as.matrix(host_dist)[common_taxa, common_taxa])
micro_dist <- as.dist(as.matrix(micro_dist)[common_taxa, common_taxa])
mantel_res <- mantel(host_dist,  micro_dist, method = "spearman", permutations = 9999)
mantel_res
# plot
mantel_df <- data.frame(
  Host = as.vector(host_dist),
  Microbiome = as.vector(micro_dist)
)
p_mantel <- ggplot(mantel_df, aes(x = Host, y = Microbiome)) + geom_point(size = 5, alpha = 0.6) + geom_smooth(method = "lm", se = TRUE, color="black",fill="grey70") + labs(x = "Phylogenetic distance of stinkbug species", y = "Microbiome dissimilarity (Bray)") + theme_classic(base_size = 12)
p_mantel
```

The raw figures generated in R are available on this GitHub page under the names `mantel_raw_fig_BRAY.png` and `mantel_raw_fig_JACCARD.png`.

### **Network-based analysis of bacterial co-occurrence patterns**

The objective of this analysis is to characterize and compare the structure of bacterial co-occurrence networks associated with two insect host species, *Graphosoma italicum* and *Lygaeus equestris*. By constructing correlation-based networks from microbiota abundance data, this approach aims to: (i) identify patterns of association (positive and negative) between bacterial taxa, (ii) highlight highly connected taxa (network hubs) that may play a key ecological role within each host-associated microbiome, and (iii) compare the topology and key bacterial taxa between the two host species. 

To reduce noise and spurious correlations driven by rare taxa, bacterial genera are filtered based on prevalence with only taxa present in at least 33% of the samples within a given host species retained. Because microbiota data are compositional by nature, a centered log-ratio (CLR) transformation is applied after adding a pseudo-count to handle zeros. This transformation expresses each taxon relative to the geometric mean abundance of all taxa within a sample, thereby reducing spurious correlations and allowing meaningful correlation-based network inference. Then, for each host species, a taxon × taxon Spearman correlation matrix is computed using CLR-transformed data. To focus on biologically meaningful associations, only strong correlations are retained (|ρ|≥0.8). The filtered correlation matrices are converted into undirected weighted graphs, where nodes represent bacterial taxa and edges represent strong positive or negative associations. Isolated nodes (taxa with no remaining associations) are removed, ensuring that the final networks represent interacting bacterial communities. For each network, several node-level properties are computed and used for visualization (eg. degree, node size, etc.). Full co-occurrence networks are visualized for each host species, with edge color indicating the sign of the association (positive vs. negative) and node size reflecting connectivity. The raw figures generated in R are available on this GitHub page under the names `Network_full_with_labels_prevalence.png` and `Network_full_simplified.png`. Subnetworks composed of the top 20 most connected taxa (highest degree) were then extracted for each host species. The raw figure generated was saved under the names `Network_top20_with_labels.png`. Finally, bacterial taxa are ranked according to two complementary metrics: (i) Degree (= identifying highly connected taxa, called structural hubs), (ii) Strength (= sum of absolute edge weights, identifying taxa involved in strong associations). For each host species, the top-ranked taxa are extracted and reported, providing a concise list of candidates for further ecological or functional investigation.

The script used is as follows:

```
#### R ####
# load the libraries
library(tidyverse)
library(compositions)
library(igraph)
# load the data file
data <- read.delim("5-Stinkbugs-abundance.txt", header = TRUE, sep = "\t", check.names = FALSE)
# Select the two species G. italicum and L. equestris
data_filt <- data %>%
  filter(species %in% c("Graphosoma_italicum", "Lygaeus_equestris"))
# Define columns of metadata and microbiota
meta_cols <- c("sample_full_name", "sample", "samplebis", "species", "family", "Total_nb_reads")
bact_cols <- setdiff(colnames(data_filt), meta_cols)
# create a subdatafile per species
GI <- data_filt %>% filter(species == "Graphosoma_italicum") %>% select(all_of(bact_cols))
LE <- data_filt %>% filter(species == "Lygaeus_equestris") %>% select(all_of(bact_cols))
# transform the sudtafiles in presence/absence matrices and filter the subdatafiles to remove rare genera by only keeping thoses presents in 33% of the samples
filter_taxa <- function(mat, prevalence = 0.33) {
  keep <- colMeans(mat > 0) >= prevalence
  mat[, keep]
}
GI_f <- filter_taxa(GI, prevalence = 0.33)
LE_f <- filter_taxa(LE, prevalence = 0.33)
# transformation CLR (Centered Log-Ratio)
clr_transform <- function(mat) {
  mat <- mat + 1  # pseudo-count
  clr(mat)
}
GI_clr <- clr_transform(GI_f)
LE_clr <- clr_transform(LE_f)
# Create a correlation matrix taxon × taxon with Spearman
GI_cor <- cor(GI_clr, method = "spearman")
LE_cor <- cor(LE_clr, method = "spearman")
# only keep strong correlation |ρ| ≥ 0.8
threshold <- 0.8
GI_cor[abs(GI_cor) < threshold] <- 0
LE_cor[abs(LE_cor) < threshold] <- 0
#####################################################
# Preparation to construct the plots
#####################################################
build_graph <- function(cor_mat) {
  graph_from_adjacency_matrix(
    cor_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
}
g_GI <- build_graph(GI_cor)
g_LE <- build_graph(LE_cor)
# remove isolated nodes without connexions
g_GI <- delete_vertices(g_GI, degree(g_GI) < 1)
g_LE <- delete_vertices(g_LE, degree(g_LE) < 1)
# Estimate the number of degree (connexions) of each node (bacterial taxa)
V(g_GI)$degree <- degree(g_GI) 
V(g_LE)$degree <- degree(g_LE)
# Scale normalization of the node's lenght according to the nb of degree
scale_size <- function(x, min = 6, max = 20) { # 
  (x - min(x)) / (max(x) - min(x)) * (max - min) + min
} 
V(g_GI)$size <- scale_size(V(g_GI)$degree)
V(g_LE)$size <- scale_size(V(g_LE)$degree)
### Facultative : Color the node according to bacterial prevalences
# Estimate the prevalence of each bacteria in the sample
GI_prev <- colSums(GI_f > 0) 
LE_prev <- colSums(LE_f > 0)
# Attribute each prevalence to each node
V(g_GI)$prevalence <- GI_prev[V(g_GI)$name] 
V(g_LE)$prevalence <- LE_prev[V(g_LE)$name]
# Set the color of the node proportional to the prevalence
pal <- colorRampPalette(c("white", "steelblue")) 
add_alpha <- function(col, alpha = 0.7) {
  adjustcolor(col, alpha.f = alpha)
}
V(g_GI)$color <- add_alpha(
  pal(100)[cut(V(g_GI)$prevalence, breaks = 100)],
  alpha = 0.7
)
V(g_LE)$color <- add_alpha(
  pal(100)[cut(V(g_LE)$prevalence, breaks = 100)],
  alpha = 0.7
)
# Create the fig's legend. Example of G.italicum
deg_GI <- sort(unique(V(g_GI)$degree))
prev_GI <- sort(unique(V(g_GI)$prevalence))
deg_legend_GI <- round(c(min(deg_GI),
                         median(deg_GI),
                         max(deg_GI)))
prev_legend_GI <- round(c(min(prev_GI),
                          median(prev_GI),
                          max(prev_GI)))
size_legend_GI <- scale_size(deg_legend_GI,
                             min = min(V(g_GI)$size),
                             max = max(V(g_GI)$size))
pal <- colorRampPalette(c("white", "steelblue"))
col_legend_GI <- add_alpha(
  pal(100)[cut(prev_legend_GI, breaks = 100)],
  alpha = 0.7
)
#####################################################
# Plots #
#####################################################
# Full network G.italicum & L.equestris with node's color depending on prevalence
plot_network <- function(g, title) {
  E(g)$color <- ifelse(E(g)$weight > 0, "firebrick", "steelblue")
  E(g)$width <- abs(E(g)$weight) * 4
  V(g)$label <- ifelse(V(g)$degree >= quantile(V(g)$degree, 0.8), V(g)$name, NA)
  plot(g, vertex.size = V(g)$size, vertex.color = V(g)$color, vertex.frame.color = "black",
       vertex.label = V(g)$name, vertex.label.color = "black", vertex.label.family = "Arial",
       vertex.label.cex = 0.7, main = title)
  legend("topleft", 
         legend = c("Positive association","Negative association",
                    paste("Degree =", deg_legend_GI), paste("Prevalence =", prev_legend_GI)),
         col = c("firebrick", "steelblue", rep("black", 3), rep("black", 3)),
         lwd = c(2, 2, rep(NA, 6)),
         pch = c(NA, NA, rep(21, 3), rep(21, 3)),
         pt.cex = c(NA, NA, size_legend_GI / max(V(g_GI)$size) * 3, rep(2, 3)),
         pt.bg = c(NA, NA, rep("grey80", 3), col_legend_GI),
         bty = "n",
         cex = 0.8)
}
par(mfrow = c(1,2))
plot_network(g_GI, "Graphosoma italicum")
plot_network(g_LE, "Lygaeus equestris")

# Simplified full network G.italicum & L.equestris without labels and without node's color depending on prevalence
plot_network_suppl <- function(g, title) {
  E(g)$color <- ifelse(E(g)$weight > 0, "firebrick", "steelblue")
  E(g)$width <- abs(E(g)$weight) * 4
  V(g)$size <- scale_size(V(g)$degree)
  V(g)$color <- "grey20"
  V(g)$label <- NA
  plot(g, vertex.size = V(g)$size, vertex.color = V(g)$color, vertex.frame.color = "black", 
       vertex.label = V(g)$label, vertex.label.color = "black", vertex.label.family = "Arial",
       vertex.label.cex = 0.7, main = title)
  legend("topleft",
         legend = c("Positive association", "Negative association", paste("Degree =", deg_legend_GI)),
         col = c("firebrick", "steelblue", rep("black", 3), rep("black", 3)),
         lwd = c(2, 2, rep(NA, 6)),
         pch = c(NA, NA, rep(21, 3), rep(21, 3)),
         pt.cex = c(NA, NA, size_legend_GI / max(V(g_GI)$size) * 3, rep(2, 3)),
         pt.bg = c(NA, NA, rep("grey80", 3), col_legend_GI),
         bty = "n",
         cex = 0.8)
}
par(mfrow = c(1,2))
plot_network_suppl(g_GI, "Graphosoma italicum")
plot_network_suppl(g_LE, "Lygaeus equestris")

# Simplified network of only the top 20 most connected bacterial taxa
# extraction of the top 20 most connected taxa
keep_GI <- names(sort(degree(g_GI), decreasing = TRUE))[1:20]
g_GI_sub <- induced_subgraph(g_GI, vids = keep_GI)
keep_LE <- names(sort(degree(g_LE), decreasing = TRUE))[1:20]
g_LE_sub <- induced_subgraph(g_LE, vids = keep_LE)
# plot
plot_subnetwork <- function(g_sub, title) {
  E(g_sub)$color <- ifelse(E(g_sub)$weight > 0, "firebrick", "steelblue")
  E(g_sub)$width <- abs(E(g_sub)$weight) * 4
  V(g_sub)$size <- scale_size(V(g_sub)$degree, min = 15, max = 40)
  pal <- colorRampPalette(c("grey10", "grey10"))
  V(g_sub)$color <- add_alpha(pal(100)[cut(V(g_sub)$prevalence, breaks = 100)], alpha = 0.8)
  V(g_sub)$label <- V(g_sub)$name
  plot(g_sub,
       vertex.size = V(g_sub)$size,
       vertex.color = V(g_sub)$color,
       vertex.frame.color = "black",
       vertex.label.color = "black",
       vertex.label.family = "Arial",
       vertex.label.cex = 0,
       main = title,
       font.main = 1)
  legend("topleft",
         legend = c("Positive association", "Negative association"),
         col = c("firebrick", "steelblue"),
         lwd = 2,
         bty = "n")
}
par(mfrow = c(1,2))
plot_subnetwork(g_GI_sub, "Graphosoma italicum - Top 20 hubs")
plot_subnetwork(g_LE_sub, "Lygaeus equestris - Top 20 hubs")
#####################################################
```

We then compared the degrees and strengths of the top 20 most connected bacteria between *G. italicum* and *L. equestris* using Wilcoxon tests:
```
#### R ####
# Identify the ranking of the 20 most connected bacterial taxa
top_taxa_by_degree <- function(g, n = 20) {
  data.frame(
    Taxon = V(g)$name,
    Degree = degree(g),
    Strength = strength(g, weights = abs(E(g)$weight)),
    Prevalence = V(g)$prevalence
  ) %>%
    arrange(desc(Degree)) %>%
    head(n)
}
top_taxa_by_strength <- function(g, n = 20) {
  data.frame(
    Taxon = V(g)$name,
    Degree = degree(g),
    Strength = strength(g, weights = abs(E(g)$weight)),
    Prevalence = V(g)$prevalence
  ) %>%
    arrange(desc(Strength)) %>%
    head(n)
}
top20_GI <- top_taxa_by_degree(g_GI, n = 20) %>%
  mutate(Group = "GI")
top20_LE <- top_taxa_by_degree(g_LE, n = 20) %>%
  mutate(Group = "LE")
top20_all <- bind_rows(top20_GI, top20_LE)
table(top20_all$Group)
#test
wilcox.test(Degree ~ Group, data = top20_all, alternative = "less")
wilcox.test(Strength ~ Group, data = top20_all, alternative = "less")
# plot
boxplot(Degree ~ Group, data = top20_all)
boxplot(Strength ~ Group, data = top20_all)
```
