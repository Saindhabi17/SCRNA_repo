# SCRNA_repo
The Packages I installed with the libraries are listed here: 
```R
# PACKAGES: 

install.packages('BiocManager')

install.packages("Seurat")
library("Seurat")
install.packages("tidyverse")
library("tidyverse")
install.packages("patchwork")
library("patchwork")
install.packages("cowplot")
library("cowplot")
install.packages("HGNChelper")
library("HGNChelper")

install.packages("dplyr")
library(dplyr)
library(tidyr)
library(ggplot2)


BiocManager::install('multtest')
BiocManager::install("ensembldb")
BiocManager::install("GenomicFeatures")
install.packages("GenomicFeatures")
BiocManager::install("MatrixGenerics")
install.packages('metap')
library("metap")
library("AnnotationHub")
library("ensembldb")
library("GenomicFeatures")
library("MatrixGenerics")
```

# Loading Single Cell RNA Seq Count Data
Before starting with vigorous code chunks, I have created directories - 
1. for my own computer
2. for the workstation

They are listed below.

```R
# Directory
filepath_homecomp = "/Users/HP/Downloads/filtered"
filepath_workstation = "/Users/andrew/Downloads/filtered"
```
I have started working with the filtered files that were prepared before using [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) from the FASTQ sequencing data generated from ScRNAseq. 

## Loading Data and Creating a Merged Seurat Object
Here, I have 11 sample files, on which I have executed my work. For this, I have first loaded the data files and then made a seurat object. 
```R
# creating list of samples
samples <- list.files("/Users/andrew/Downloads/filtered/")

# read in 10X data and creating Seurat Object
for (file in samples){
  print(paste0(file))
  seurat_data <- Read10X(data.dir = paste0("/Users/andrew/Downloads/filtered/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# exploring the metadata (for e.g.- SRR12603783) 
View(SRR12603783@meta.data)

# merging all object into single seurat object
merged_seurat <- merge(x = SRR12603780, 
                       y = c(SRR12603781,
                             SRR12603782,
                             SRR12603783,
                             SRR12603784,
                             SRR12603785,
                             SRR12603786,
                             SRR12603787,
                             SRR12603788,
                             SRR12603789,
                             SRR12603790),
                       add.cell.id = samples)
```
# Quality Control
At first the meta data of the merged seurat object was explored.
```R
# exploring merged meta data
View(merged_seurat@meta.data)
```
There are 3 columns in the merged meta data. They are- 
1. orig.ident: The first column contains the sample identity as known. By default it shows the value provided for the project argument when loading in the data.
2. nCount_RNA: This column represents the number of UMIs per cell. UMI (unique molecular identifiers) is used to determine whether a read is a biological or technical duplicate (PCR duplicate). There can be 2 types of duplicates - Biological Duplicates - Reads with different UMIs mapping to the same transcript derived from different molecules, and Technical Duplicates - Reads with the same UMI originated from the same molecule. For the biological duplicates each read should be counted where for the technical duplicates reads should be counted as a single one.
3. nFeature_RNA: This column represents the number of genes detected per cell.
## Recommended Features to Add to the Metadata 
1. Novelty Score: It is the number of genes detected per UMI. More genes detected per UMI, more complex the data will be.
2. Mitochondrial Ratio: This metric will give us a percentage of cell reads originating from the mitochondrial genes (coming from dying cells).
```R
# adding number of genes per UMI (Novelty score) for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# computing & adding percentage mitochondrial ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

# creating metadata dataframe
merged_metadata <- merged_seurat@meta.data

# adding cell IDs to metadata
merged_metadata$cells <- rownames(merged_metadata)

# adding sample type to metadata. The original file could be download from SRA explorer.
SampleType <- c("BLCA", "BLCA", "Normal", "BLCA", "BLCA", "BLCA", "BLCA", "BLCA", "BLCA", "Normal", "Normal")

# sample type with grade (Not tested)
# SampleType <- c("BLCA_LG", "BLCA_LG", "Normal", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "BLCA_HG", "Normal", "Normal")

names(SampleType) <- c("SRR12603789", "SRR12603790", "SRR12603788", "SRR12603787", "SRR12603786", "SRR12603785", "SRR12603784", "SRR12603783", "SRR12603782", "SRR12603781", "SRR12603780")

merged_metadata$sampleType <- stringr::str_replace_all(merged_metadata$orig.ident, SampleType)

# renaming columns
merged_metadata <- merged_metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA,
                sample = sampleType)

# adding merged metadata back to Seurat object
merged_seurat@meta.data <- merged_metadata
```
After this, I have created .RData object that can be loaded at any time. 
```R
save(merged_seurat, file="/Users/andrew/Downloads/filtered/merged_filtered_seurat.RData")
```
## Visualizing the Plots 
In this part of QC, there will be various plots that can help to understand how this QC is going on. 

# cell counts per sample
```R
# visualizing the number of cell counts per sample
merged_metadata %>%
  ggplot(aes(x=seq_folder, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells before QC")
```


# UMIs per sample

# Visualizing the number UMIs/transcripts per cell
merged_metadata %>% 
  ggplot(aes(x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 1000) +
  labs(fill = "Sample")



# Genes detected per cell

# Visualizing the distribution of genes detected per cell via histogram
merged_metadata %>% 
  ggplot(aes(x=nGene, fill= sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 500) +
  labs(fill = "Sample")



# Novelty score

# Visualizing the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
merged_metadata %>%
  ggplot(aes(x=log10GenesPerUMI, fill=sample)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~seq_folder) +
  xlab("Novelty Score") +
  geom_vline(xintercept = 0.8)



# Mitochondrial gene expression detected per cell 

# Visualizing the distribution of mitochondrial gene expression detected per cell
merged_metadata %>%
  ggplot(aes(x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10() + 
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)) + 
  theme_classic() +
  facet_wrap(~seq_folder) +
  geom_vline(xintercept = 0.2)



# Joint filtering: nUMI, nGene and mitoRatio

# Visualizing the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
merged_metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~seq_folder)



# FILTERING 


# Cell-level filtering

# -nUMI > 1000
# -nGene > 500 & < 6000
# -log10GenesPerUMI > 0.8
# -mitoRatio < 0.10

# filtering out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(merged_seurat, 
                          subset= nUMI >= 1000 &
                            nGene >= 500 &
                            nGene <= 6000 & 
                            log10GenesPerUMI > 0.80 & 
                            mitoRatio < 0.10) 

# exploring filtered seurat meta data
View(filtered_seurat@meta.data)



# Gene-level filtering - Keeping only genes which are expressed in 100 or more cells (usually this is 10)

# Extracting counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 100 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 100

# Only keeping those genes expressed in more than 100 cells
filtered_counts <- counts[keep_genes, ]
# Reassigning to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# exploring filtered seurat metadata after gene-level filtering
View(filtered_seurat@meta.data)

# Creating .RData object to load at any time
save(filtered_seurat, file="seurat_filtered.RData")



# Re-assess QC Metrics


# saving filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

# to see drop in filtering cells:
met_before <- data.frame(unclass(table(merged_metadata$seq_folder)))
met_before$QCgroup <- "before"
met_before$cell<- rownames(met_before)
names(met_before)[1] <- 'count'

met_after <- data.frame(unclass(table(metadata_clean$seq_folder)))
met_after$QCgroup <- "after"
met_after$cell<- rownames(met_after)
names(met_after)[1] <- 'count'

# count
cell_count <- data.frame(rbind(met_before, met_after))

# visualization :
cell_count %>% ggplot(aes(x=cell, y=count, fill=QCgroup)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_manual(values = c("#CC6666", "#9999CC")) +
  xlab("samples") +
  ggtitle("nCells count before and after QC")



# Visualizing the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs

metadata_clean %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~seq_folder)




# NORMALIZATION, REGRESSING OUT UNWANTED VARIATION, CLUSTERING

# Exploring sources of unwanted variation

# Evaluating effects of cell cycle and mitochondrial expression - 
# First scoring the cells for cell cycle genes
# Then determining whether cell cycle is a major source of variation in our data set using PCA.

# SCTransform functions can also be used for normalization that will replace other functions single handedly 

# Normalizing the counts
# This normalization method is solely for the purpose of exploring the sources of variation in our data.
seurat_phase <- NormalizeData(filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

# Loading cell cycle markers
load("C:/Users/andrew/AppData/Local/Packages/microsoft.windowscommunicationsapps_8wekyb3d8bbwe/LocalState/Files/S0/95/Attachments/cycle[131].rda")

# Scoring cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# Viewing cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data) 
table(seurat_phase$Phase)

# G1 49590
# G2M 11433
# S 25096

# Identifying the most variable genes and scaling them
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst", 
                                     nfeatures = 2000, 
                                     verbose = TRUE)

# Identifying the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_phase), 10)

# plotting variable features with and without labels

plot1 <- VariableFeaturePlot(seurat_phase)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.new()
plot1 + plot2
library(ggplot2)
ggsave('plot1.png', plot1)
ggsave('plot2.png', plot2)
print(plot2)
print(plot1)


# Checking quartile values for mitoRatio, we will use this variable later to mitigate unwanted source of variation in dataset
summary(seurat_phase@meta.data$mitoRatio)

# Turning mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                                     breaks=c(-Inf, 0.015, 0.025, 0.045, Inf), 
                                     labels=c("Low","Medium","Medium high", "High"))



# Scaling the counts
# This step is essential for PCA , clustering and heatmap generation
seurat_phase <- ScaleData(seurat_phase)
saveRDS(seurat_phase, "seurat_phase.rds")

# Performing PCA
seurat_phase <- RunPCA(seurat_phase)

# Plotting the PCA colored by cell cycle phase
no_split <- DimPlot(seurat_phase,
                    reduction = "pca",
                    group.by= "Phase")
no_split 

with_split <- DimPlot(seurat_phase,
                      reduction = "pca",
                      group.by= "Phase",
                      split.by= "Phase")

with_split

#dev.new(width=25, height=5)

no_split + with_split



# Plotting the PCA colored by mitochondrial expression
no_split <- DimPlot(seurat_phase,
                    reduction = "pca",
                    group.by= "mitoFr")

with_split <- DimPlot(seurat_phase,
                      reduction = "pca",
                      group.by= "mitoFr",
                      split.by= "mitoFr")

no_split + with_split



# INTEGRATION

# adjusting the limit for allowable object sizes within R
options(future.globals.maxSize = 4000 * 1024^2)

# Splitting seurat object by group
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

# then normalizing by SCTansform
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))
}

# to see what the component of the object are. 

split_seurat  

#$Normal
#An object of class Seurat
#46094 features across 20167 samples within 2 assays
#Active assay: SCT (22984 features, 3000 variable features)
# 1 other assay present: RNA
# 1 dimensional reduction calculated: pca

#$BLCA
#An object of class Seurat
#46220 features across 65952 samples within 2 assays
#Active assay: SCT (23110 features, 3000 variable features)
# 1 other assay present: RNA
# 1 dimensional reduction calculated: pca


# Selecting the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 


# Preparing the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

# Finding best buddies (using canonical correlation analysis or CCA) - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrating across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# Checking assays in the object:
split_seurat$Normal@assays

#$RNA
#Assay data with 23110 features for 20167 cells
#Top 10 variable features:
#  IGLC2, IGKC, IGHA1, IGHG1, JCHAIN, CCL21, CCL19, IGHG3, IGLC1, IGHG2 

#$SCT
#SCTAssay data with 22984 features for 20167 cells, and 1 SCTModel(s) 
#Top 10 variable features:
#  TPSB2, S100A9, PLA2G2A, S100A8, TPSAB1, LYZ, CCL4, PTGDS, SPINK1, CFD 

# Running PCA
seurat_integrated <- RunPCA(object = seurat_integrated, verbose = TRUE)

#PC_ 1 
#Positive:  CD52, PTPRC, CCL5, CD3D, KRT19, RPS19, TRAC, HCST, RPS27, S100P 
#           SRGN, SAMSN1, FXYD3, TRBC2, CD2, RGS1, CD7, CXCR4, RPS29, KRT13 
#           ARHGAP15, PSCA, CLDN4, NKG7, CSTB, CORO1A, CST7, GZMA, CD69, RPL41 
#Negative:  IGFBP7, MGP, SPARC, SPARCL1, VIM, CALD1, IFITM3, A2M, COL1A2, COL4A1 
#           TAGLN, COL4A2, COL6A2, BGN, COL3A1, NNMT, TCF4, CCL2, COL1A1, MYL9 
#           DCN, MT2A, CCN2, CCN1, ADAMTS9, SELENOM, IGFBP4, LUM, GNG11, TIMP1 

#PC_ 2 
#Positive:  COL1A2, COL3A1, COL1A1, TAGLN, DCN, LUM, BGN, C1R, COL6A2, C1S 
#           MFAP4, SOD3, RARRES2, TPM2, MYL9, PRKG1, CRYAB, ACTA2, COL6A3, LGALS1 
#           SERPINF1, COL6A1, CALD1, PCOLCE, TIMP1, AEBP1, C11orf96, MEG3, GPC6, CRISPLD2 
#Negative:  PLVAP, PCAT19, CALCRL, AQP1, MCTP1, VWF, RAMP2, PECAM1, LDB2, RAMP3 
#           FLT1, ZNF385D, TCF4, HSPG2, SPARCL1, ADGRL4, CD74, ACKR1, CLDN5, EMCN 
#           SLCO2A1, SELE, DOCK4, CCL14, ECSCR, GNG11, ERG, RNASE1, ADAMTS9, INSR 

#PC_ 3 
#Positive:  SPINK1, ADIRF, CSTB, CCT2, CCND1, S100P, UCA1, FXYD3, YEATS4, KRT19 
#           RAB3IP, IFI27, S100A6, KRT7, PSCA, GAPDH, IRS2, FRS2, FCHSD2, MIF 
#           SNCG, MYO16, CNOT2, KRT18, CLDN4, CD24, AC025159.1, S100A14, S100A2, GDF15 
#Negative:  CCL5, B2M, CD52, PTPRC, SRGN, IL32, HSPA1A, CD3D, NKG7, GZMA 
#           TRAC, HCST, ARHGAP15, CXCR4, RGS1, RGS2, SAMSN1, CCL4, CORO1A, CD7 
#           FYN, CD2, CST7, CD69, DNAJB1, HLA-DPB1, STAT4, HLA-DRA, CD74, PTPN22 

#PC_ 4 
#Positive:  HLA-DRA, CD74, TYROBP, HLA-DRB1, HLA-DPB1, FCER1G, AIF1, HLA-DPA1, FTL, HLA-DQA1 
#           IFI30, LYZ, HLA-DQB1, C1QA, C1QB, LST1, MS4A6A, C1QC, CD14, APOE 
#           TMEM176B, S100A9, HLA-DMA, FCGR2A, SPI1, CST3, CD68, PSAP, FTH1, HLA-DMB 
#Negative:  IL32, CCL5, CD3D, CRIP1, TRAC, FYN, CD2, IGFBP7, GZMA, CD7 
#           COL4A1, NKG7, TRBC2, CALD1, COL4A2, CD3E, SKAP1, TRBC1, MCAM, MYL9 
#           CD247, CAMK4, RGS5, COL18A1, NDUFA4L2, CYTOR, PPP1R16B, ITGA1, ACTA2, PTPRC 

#PC_ 5 
#Positive:  LUM, MMP2, PTGDS, DCN, RARRES2, FBLN1, LSAMP, SERPINF1, COL8A1, PDPN 
#           CTSK, C1S, CLMP, VCAN, APOD, SFRP2, TSHZ2, POSTN, PDGFRA, RARRES1 
#           FAP, EFEMP1, CFD, CXCL1, NBL1, ABI3BP, BICC1, C1R, MFAP4, CTHRC1 
#Negative:  RGS5, ACTA2, NDUFA4L2, MYL9, PPP1R14A, TAGLN, FRZB, CRIP1, CALD1, MYH11 
#           IGFBP7, GJA4, PRKG1, COL18A1, MCAM, TPPP3, MUSTN1, COX4I2, COL4A1, COL4A2 
#           PTP4A3, MYLK, CDH6, MFGE8, SOD3, TYROBP, HEYL, HIGD1B, WFDC1, HLA-DRA 


# Plotting PCA
png(filename = "PCA_integrated.png", width = 16, height = 8.135, units = "in", res = 300)
PCAPlot(seurat_integrated,
        split.by = "sample")
dev.off() 

seurat_integrated
#An object of class Seurat 
#49220 features across 86119 samples within 3 assays 
#Active assay: integrated (3000 features, 3000 variable features)
#2 other assays present: RNA, SCT
#1 dimensional reduction calculated: pca



# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca",
                             verbose = TRUE)

# Plot UMAP 
png(filename = "UMAP_integrated.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated, split.by = "sample")
dev.off()



# CLUSTERING CELLS BASED ON TOP PCS (METAGENES)


# Identify significant PCs

# Explore heatmap of PCs
png(filename = "heatmap_integrated_2.png", width = 16, height = 8.135, units = "in", res = 300)
DimHeatmap(seurat_integrated, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)
dev.off()

# Printing out the most variable genes driving PCs
print(x = seurat_integrated[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# PC_ 1 
# Positive:  CD52, PTPRC, CCL5, CD3D, KRT19 
# Negative:  IGFBP7, MGP, SPARC, SPARCL1, VIM 

# PC_ 2 
# Positive:  COL1A2, COL3A1, COL1A1, TAGLN, DCN 
# Negative:  PLVAP, PCAT19, CALCRL, AQP1, MCTP1 

# PC_ 3 
# Positive:  SPINK1, ADIRF, CSTB, CCT2, CCND1 
# Negative:  CCL5, B2M, CD52, PTPRC, SRGN 

# PC_ 4 
# Positive:  HLA-DRA, CD74, TYROBP, HLA-DRB1, HLA-DPB1 
# Negative:  IL32, CCL5, CD3D, CRIP1, TRAC 

# PC_ 5 
# Positive:  LUM, MMP2, PTGDS, DCN, RARRES2 
# Negative:  RGS5, ACTA2, NDUFA4L2, MYL9, PPP1R14A 

# PC_ 6 
# Positive:  FABP5, RPS19, CRH, LY6D, KRT20 
# Negative:  CCT2, ADIRF, SPINK1, HSPA1A, CCND1 

# PC_ 7 
# Positive:  SPINK1, CRH, LCN15, CCT2, PLA2G2A 
# Negative:  KRT13, LYPD3, PLAUR, OLFM4, SFN 

# PC_ 8 
# Positive:  SPARC, COL4A1, INSR, IGFBP3, COL4A2 
# Negative:  ACKR1, FOS, JUN, ZNF385D, SELE 

# PC_ 9 
# Positive:  LCN15, PLA2G2A, FABP4, CRTAC1, LINC01088 
# Negative:  H19, RPS19, CRH, AP005230.1, HES1 

# PC_ 10 
# Positive:  CRH, CCL5, LY6D, RPS19, ACKR1 
# Negative:  HSPA1A, HSPA1B, DNAJB1, LCN15, HSP90AA1 


# To determine how many Pcs should be considered for clustering:
# Plotting the elbow plot
png(filename = "elbow.png", width = 16, height = 8.135, units = "in", res = 300)
ElbowPlot(object = seurat_integrated, 
          ndims = 40)
dev.off()


# to make it more quantitative :
# Determining percent of variation associated with each PC
pct <- seurat_integrated[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

pct

#[1] 7.9596270 6.4804792 5.5083025 4.4609188 4.0031754 3.5967530 3.3963246 3.0842674
#[9] 2.8526141 2.6327352 2.4536622 2.4075917 2.2740071 2.2180154 2.0701053 2.0233150
#[17] 1.9361711 1.7899755 1.7465302 1.6432391 1.6341522 1.5688985 1.5288258 1.4906526
#[25] 1.4657534 1.3954436 1.3468486 1.2975141 1.2807993 1.2304310 1.2286440 1.2168580
#[33] 1.2001166 1.1807094 1.1542407 1.1514479 1.1161946 1.0873437 1.0653210 1.0633526
#[41] 1.0481415 1.0289363 1.0172385 1.0054600 0.9848788 0.9677587 0.9398681 0.9329160
#[49] 0.9243895 0.9090550

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
cumu

#[1]   7.959627  14.440106  19.948409  24.409328  28.412503  32.009256  35.405581
#[8]  38.489848  41.342462  43.975197  46.428860  48.836451  51.110458  53.328474
#[15]  55.398579  57.421894  59.358065  61.148041  62.894571  64.537810  66.171962
#[22]  67.740861  69.269686  70.760339  72.226092  73.621536  74.968385  76.265899
#[29]  77.546698  78.777129  80.005773  81.222631  82.422748  83.603457  84.757698
#[36]  85.909146  87.025340  88.112684  89.178005  90.241358  91.289499  92.318435
#[43]  93.335674  94.341134  95.326013  96.293771  97.233639  98.166555  99.090945
#[50] 100.000000

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# 40

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
# 20

# Minimum of the two calculation is the optimal number of PC to pick.
pcs <- min(co1, co2)
pcs
# 20


# CLUSTER THE CELLS - VISUALIZATION

# to check what is active assay
DefaultAssay(object = seurat_integrated)

# Determining the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:18)

#Find clusters
# Determining the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.4))


# Exploring resolutions
head(seurat_integrated@meta.data)

# Assigning identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.2"

# Plotting the UMAP
png(filename = "umap_cluster_with_label.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 6)
dev.off()


#saving the file for further use
save(seurat_integrated, file="seurat_integrated.RData")



# CLUSTERING QC

# Segregation of clusters by sample

# Extracting identity and sample information from seurat object to determine the number of cells per cluster per sample

library(dplyr)
library(tidyr)

# n_cells <- FetchData(seurat_integrated, 
#                      vars = c("ident", "orig.ident")) %>%
#         dplyr::count(ident, orig.ident) %>%
#         tidyr::spread(ident, n)

n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "orig.ident"))
n_cells <- dplyr::count(n_cells, ident, orig.ident)
n_cells <- tidyr::spread(n_cells, ident, n)



#Ading sample data from paper; we expect to see samples from same group have more or less similar number of cells in each cluster. 
#So normal samples should show similar patterns: SRR12603780, SRR12603781, and SRR12603788.


sampleData<- data.frame(tibble::tribble(
  ~sample_id, ~gender, ~age, ~Grade, ~Invasiveness, ~Surgery_Type, ~Tumor_size_cm,
  "SRR12603790",     "M",  67L,  "low", "Noninvasive",       "TURBT",          "1.9",
  "SRR12603789",     "M",  70L,  "low", "Noninvasive",       "TURBT",          "2.5",
  "SRR12603787",     "M",  63L, "high", "Noninvasive",  "Cystectomy",          "3.5",
  "SRR12603786",     "F",  59L, "high", "Noninvasive",  "Cystectomy",          "4.7",
  "SRR12603785",     "M",  57L, "high",    "Invasive",  "Cystectomy",          "5.1",
  "SRR12603784",     "M",  75L, "high",    "Invasive",  "Cystectomy",          "4.3",
  "SRR12603783",     "M",  77L, "high",    "Invasive",  "Cystectomy",          "4.5",
  "SRR12603782",     "F",  72L, "high",    "Invasive",  "Cystectomy",          "4.1",
  "SRR12603781",     "M",  67L, "normal",   "normal",       "TURBT",            "-",
  "SRR12603788",     "M",  75L, "normal",   "normal",  "Cystectomy",            "-",
  "SRR12603780",     "M",  63L, "normal",   "normal",  "Cystectomy",            "-"
))


tmp_df<- seurat_integrated@meta.data

merged_df <- merge(tmp_df, sampleData, 
                   by.x = "orig.ident", 
                   by.y = "sample_id", 
                   all.x = TRUE)

seurat_integrated@meta.data<-merged_df

# Viewing table
head(n_cells)
# saving objects (to mark where and when we stored the file)

saveRDS(seurat_integrated, "seurat_integrated.RDS")

# UMAP of cells in each cluster by sample
# This would allow us to see condition specefic clusters
png(filename = "umap_cluster_sample.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
dev.off()



# Segregation of clusters by cell cycle phase (unwanted source of variation) 

# Exploring whether clusters segregate by cell cycle phase

png(filename = "umap_cluster_cell_cycle.png", width = 16, height = 8.135, units = "in", res = 300)
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()
dev.off()



# Segregation of clusters by various sources of uninteresting variation

# Determining metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio")
png(filename = "umap_unwanted_source_clustering.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
dev.off()



# Exploration of the PCs driving the different clusters

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:18),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(seurat_integrated, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(seurat_integrated, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
library(cowplot)
library(tidyverse)
library(HGNChelper)

png(filename = "umap_on_pcs.png", width = 16, height = 8.135, units = "in", res = 300)
map(paste0("PC_", 1:18), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)
dev.off()



# Examine PCA results 
print(seurat_integrated[["pca"]], dims = 1:5, nfeatures = 5)

# PC_ 1 
# Positive:  CD52, PTPRC, CCL5, CD3D, KRT19 
# Negative:  IGFBP7, MGP, SPARC, SPARCL1, VIM 

# PC_ 2 
# Positive:  COL1A2, COL3A1, COL1A1, TAGLN, DCN 
# Negative:  PLVAP, PCAT19, CALCRL, AQP1, MCTP1 

# PC_ 3 
# Positive:  SPINK1, ADIRF, CSTB, CCT2, CCND1 
# Negative:  CCL5, B2M, CD52, PTPRC, SRGN 

# PC_ 4 
# Positive:  HLA-DRA, CD74, TYROBP, HLA-DRB1, HLA-DPB1 
# Negative:  IL32, CCL5, CD3D, CRIP1, TRAC 

# PC_ 5 
# Positive:  LUM, MMP2, PTGDS, DCN, RARRES2 
# Negative:  RGS5, ACTA2, NDUFA4L2, MYL9, PPP1R14A



# Normalizing RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

# Fibroblast
png(filename = "umap_fibroblast.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("IGFBP7", "MGP"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

# Endothelial Cells
png(filename = "umap_endothelial.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("PLVAP", "CALCRL"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off()

# T Cells
png(filename = "umap_t_cells.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c("CCL5", "CD52", "IL32"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)
dev.off() 

##  



# MARKER IDENTIFICATION 


#______________________________ NOT TO BE RUN________________________________
# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = seurat_integrated, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)


# Explecity set the defult object to normalized values
DefaultAssay(seurat_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var= "Invasiveness",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.60)

view(seurat_cluster@meta.data)

# Adding more annotations to the result

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  tryCatch({
    FindConservedMarkers(seurat_integrated,
                         ident.1 = cluster,
                         grouping.var = "Invasiveness",
                         only.pos = TRUE,
                         logfc.threshold = 0.60) %>%
      rownames_to_column(var = "gene") %>%
      left_join(y = unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")) %>%
      cbind(cluster_id = cluster, .)
  },
  error = function(e) {
    message(paste0("Error: ", e$message))
    return(NULL)
  }
  )
}
# this function can be an argument for 'map_dfr' function :
# Iterate function across desired clusters
   

view(conserved_markers)

# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (Normal_avg_log2FC + BLCA_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

head(top10)


# A tibble: 6 × 16
# Groups:   cluster_id [1]
#cluster_id gene  Normal_p_val Normal_avg_log2FC Normal_pct.1 Normal_pct.2 Normal_p_val_adj BLCA_p_val BLCA_avg_log2FC BLCA_pct.1 BLCA_pct.2 BLCA_p_val_adj max_pval minimump_p_val
#<int> <chr>        <dbl>             <dbl>        <dbl>        <dbl>            <dbl>      <dbl>           <dbl>      <dbl>      <dbl>          <dbl>    <dbl>          <dbl>
#1          0 CD2              0              1.08        0.608        0.032                0          0           1.11       0.656      0.057              0        0              0
#2          0 PTPRC            0              1.37        0.917        0.184                0          0           1.31       0.832      0.128              0        0              0
#3          0 CXCR4            0              1.24        0.852        0.233                0          0           1.08       0.645      0.125              0        0              0
#4          0 STAT4            0              1.11        0.716        0.134                0          0           0.963      0.538      0.061              0        0              0
#5          0 FYN              0              1.21        0.931        0.403                0          0           1.11       0.763      0.188              0        0              0
#6          0 SYTL3            0              1.27        0.788        0.133                0          0           0.931      0.549      0.101              0        0              0
# ℹ 2 more variables: description <chr>, avg_fc <dbl>



data.table::fwrite(top10, "blca_top10_conserved_markers.csv")

top10_mod <- data.frame(unclass(table(top10$gene, top10$cluster_id)))

data.table::fwrite(top10_mod, "top10_mod_conserved_markers.csv")

# markers with highest frequency
M <- c("KRT7", "KRT19", "FCER1G", "AIF1", "AQP3", "CCL5", "CD24", "CD3D", "CD52", "CLDN4", "COL1A1", "COL1A2", "CRTAC1", "CXCL8", "DCN", "FABP4", "FABP5", "FXYD3", "GZMA", "HLA-DRA", "IGHA1", "IGHG1", "IGHG3", "IGHG4", "IGKC", "IGLC1", "IGLC2", "IGLC3", "JCHAIN")

View("blca_top10_conserved_markers.csv")

# Plot interesting marker gene expression - basal cells
png(filename = "umap_high_freq_basal_cells.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = c("KRT7", "KRT19", "AQP3", "CD24", "FXYD3", "CXCL8"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - basal cells
png(filename = "violin_high_freq_basal_cells.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("KRT7", "KRT19", "AQP3", "CD24", "FXYD3", "CXCL8"))
dev.off() 


# Plot interesting marker gene expression - plasma cells
png(filename = "umap_high_freq_Plasma_cells.png", width = 26, height = 15.135, units = "in", res = 600)
FeaturePlot(object = seurat_integrated, 
            features = c("IGHA1","IGHG1","IGHG3","IGHG4","IGKC","IGLC1","IGLC3","IGLC3","JCHAIN"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - plasma cells
png(filename = "violin_high_freq_Plasma_cells.png", width = 26, height = 10.135, units = "in", res = 600)
VlnPlot(object = seurat_integrated, 
        features = c("IGHA1","IGHG1","IGHG3","IGHG4","IGKC","IGLC1","IGLC3","IGLC3","JCHAIN"))
dev.off()


# Plot interesting marker gene expression - dendritic cells
png(filename = "umap_high_freq_dc.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = c("FCER1G","AIF1","FABP4"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - dendritic cells
png(filename = "violin_high_freq_dc.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("FCER1G","AIF1","FABP4"))
dev.off()


# Plot interesting marker gene expression - luminal epithelial cells
png(filename = "umap_high_freq_luminal_epithelial.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = c("CLDN4"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - luminal epithelial cells
png(filename = "violin_high_freq_luminal_epithelial.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("CLDN4"))
dev.off()


# Plot interesting marker gene expression - epithelial cells
png(filename = "umap_high_freq_epithelial.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = c("CRTAC1","FXYD3"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - epithelial cells
png(filename = "violin_high_freq_epithelial.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("CRTAC1","FXYD3"))
dev.off()



# Plot interesting marker gene expression - endothelial cells
png(filename = "umap_high_freq_endothelial.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = c("FABP4", "FABP5"),
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()


# Vln plot - cluster 0 - endothelial cells
png(filename = "violin_high_freq_endothelial.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = c("FABP4", "FABP5"))
dev.off()



d <- data.table::fread("blca_top10_conserved_markers.csv")

png(filename = "umap_cluster4_markers.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = d$gene[d$cluster_id == "4"],
            sort.cell = NULL,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()

# Vln plot - cluster 0
png(filename = "violin_cluster4_markers.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = d$gene[d$cluster_id == "4"])
dev.off() 

png(filename = "umap_cluster1_markers.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = d$gene[d$cluster_id == "1"],
            sort.cell = NULL,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()

# Vln plot - cluster 0
png(filename = "violin_cluster1_markers.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = d$gene[d$cluster_id == "1"])
dev.off()

png(filename = "umap_cluster12_markers.png", width = 16, height = 8.135, units = "in", res = 300)
FeaturePlot(object = seurat_integrated, 
            features = d$gene[d$cluster_id == "12"],
            sort.cell = NULL,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)
dev.off()

# Vln plot - cluster 0
png(filename = "violin_cluster12_markers.png", width = 16, height = 8.135, units = "in", res = 300)
VlnPlot(object = seurat_integrated, 
        features = d$gene[d$cluster_id == "12"])
dev.off()



```


## Issues: 
My first problem is to get run the seurat_integrated object that was merged with sampleData. The codes I have used for meging are-
```R
tmp_df<- seurat_integrated@meta.data

merged_df <- merge(tmp_df, sampleData, 
                   by.x = "orig.ident", 
                   by.y = "sample_id", 
                   all.x = TRUE)

seurat_integrated@meta.data<-merged_df
```

After that I have saved the new object file.
```R
# saving objects (to mark where and when we stored the file)

saveRDS(seurat_integrated, "seurat_integrated.RDS")
```

My problem started with marker identification, the code I used is-
```R
# Explecity set the defult object to normalized values
DefaultAssay(seurat_integrated) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(seurat_integrated,
                                                   ident.1 = 0,
                                                   grouping.var= "Invasiveness",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.60)
```
And the error I got- 
```R
Error in cells[[i]] <- rownames(x = object.var[object.var[, 1] == levels.split[i],  : 
  attempt to select less than one element in integerOneIndex
  ```
And then I have setup the AnnotationHub. That ran well, but I can't prepare the conserved markers. The codes for AnnotationHub (which was fine) are-
```R
# Adding more annotations to the result

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  tryCatch({
    FindConservedMarkers(seurat_integrated,
                         ident.1 = cluster,
                         grouping.var = "Invasiveness",
                         only.pos = TRUE,
                         logfc.threshold = 0.60) %>%
      rownames_to_column(var = "gene") %>%
      left_join(y = unique(annotations[, c("gene_name", "description")]),
                by = c("gene" = "gene_name")) %>%
      cbind(cluster_id = cluster, .)
  },
  error = function(e) {
    message(paste0("Error: ", e$message))
    return(NULL)
  }
  )
}
```
After that the code got stuck and was showing error. 
```R
conserved_markers <- map_dfr(c(0:15), get_conserved)
view(conserved_markers)
```
And the error is -
```R
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
Error: attempt to select less than one element in integerOneIndex
```
