# Short summary of experimentation:
# PBMCs were isolated from four healthy human blood donors.
# PBMCs were stained with 5OPRU-hMR1-Tet and positively enriched before sorting of MAIT (Hashtag1) and non-MAIT cells (Hashtag3).
# Libraries were generated following the 10X Genomics Next-GEM v2 5' Kit with Feature Barcode Technology and TCR sequencing.

# Setup -------------------------------------------------------------------

devtools::install:github('ncborcherding/scRepertoire') #  DOI: 10.1101/2024.12.31.630854, Vignette: https://www.borch.dev/uploads/screpertoire/
library(scRepertoire)

# library(reticulate)
# conda_create("r-reticulate") ##If first time using reticulate
# use_condaenv(condaenv="r-reticulate", required=TRUE)
# install_keras()
# library(keras)

devtools::install_github("BorchLab/Trex") # DOI: 10.1038/s41590-024-01888-9, Vignette: https://www.borch.dev/uploads/screpertoire/articles/trex
library(Trex)

required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "viridis",
                       "scales", "stats", "RColorBrewer", "ggpubr", "Seurat", "HGNChelper", "tidyverse")

for(package in required_packages){if(!require(package, character.only=TRUE))
{install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}
rm(required_packages)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Import TCRseq data ------------------------------------------------------

# load annotated Seurat file
int.TNK <- readRDS('int.TNK.rds') # all T and NK cells passing quality control
rm(int.TNK)
cd4tmems <- subset(int.TNK, CD4Tphenotype=='memory T') # only memory CD4+ T cells
cd4tmems$Donor_cluster <- paste0(cd4tmems$donor, '_', cd4tmems$cd4tmemsubset)

# load VDJ annotations from cellranger output
VDJ1 <- read.csv('./MAITseq_TCR_seq_data/D1_filtered_contig_annotations.csv')
VDJ2 <- read.csv('./MAITseq_TCR_seq_data/D2_filtered_contig_annotations.csv')
VDJ3 <- read.csv('./MAITseq_TCR_seq_data/D3_filtered_contig_annotations.csv')
VDJ4 <- read.csv('./MAITseq_TCR_seq_data/D4_filtered_contig_annotations.csv')

# subset VDJ information do match cd4tmem cell seurat object
cd4tmem_ids_D1 <- colnames(cd4tmems)[cd4tmems$donor=='D1']
cd4tmem_ids_D1 <- sub(".*_", "", cd4tmem_ids_D1)
cd4tmem_ids_D2 <- colnames(cd4tmems)[cd4tmems$donor=='D2']
cd4tmem_ids_D2 <- sub(".*_", "", cd4tmem_ids_D2)
cd4tmem_ids_D3 <- colnames(cd4tmems)[cd4tmems$donor=='D3']
cd4tmem_ids_D3 <- sub(".*_", "", cd4tmem_ids_D3)
cd4tmem_ids_D4 <- colnames(cd4tmems)[cd4tmems$donor=='D4']
cd4tmem_ids_D4 <- sub(".*_", "", cd4tmem_ids_D4)

VDJ1cd4tmem <- VDJ1 %>% filter(barcode %in% cd4tmem_ids_D1)
VDJ2cd4tmem <- VDJ2 %>% filter(barcode %in% cd4tmem_ids_D2)
VDJ3cd4tmem <- VDJ3 %>% filter(barcode %in% cd4tmem_ids_D3)
VDJ4cd4tmem <- VDJ4 %>% filter(barcode %in% cd4tmem_ids_D4)

contig_list <- list(VDJ1cd4tmem, VDJ2cd4tmem, VDJ3cd4tmem, VDJ4cd4tmem)
TCR <- loadContigs(contig_list)
TCR <- combineTCR(TCR, samples=c('D1', 'D2', 'D3', 'D4'))

# Add cell subset metadata from Seurat to TCR information
metadata_subset <- cd4tmems@meta.data$cd4tmemsubset
names(metadata_subset) <- rownames(cd4tmems@meta.data)
TCR <- lapply(TCR, function(x) {x <- x %>% mutate(cd4tmemsubset=as.character(metadata_subset[barcode]))})

exportClones(TCR, write.file=TRUE, file.name='TCR_exported_clones.csv')
rm(VDJ1cd4tmem, VDJ2cd4tmem, VDJ3cd4tmem, VDJ4cd4tmem, cd4tmem_ids_D1, cd4tmem_ids_D2, cd4tmem_ids_D3, cd4tmem_ids_D4, VDJ1, VDJ2, VDJ3, VDJ4, contig_list, metadata_subset)

# Data deposition
TCR_file <- do.call(rbind.data.frame,TCR)
TCR_file <- TCR_file %>% rename(CellBarcode=barcode, Donor=sample, TCRa=TCR1, TCRb=TCR2)
saveRDS(TCR_file, file='./Seurat_data/Data deposition/TCR_CD4T_cells_Kammann_et_al.rds')

# combine TCRseq data with RNAseq seurat file
CD4TCR <- combineExpression(TCR, cd4tmems, group.by='cd4tmemsubset', filterNA=TRUE)
CD4TCR@meta.data$cd4tmemsubset <- factor(CD4TCR@meta.data$cd4tmemsubset, levels=c("memory CD4 T", 'CCR10+ memory CD4 T', 'TOX+ TIGIT+ CD4 T', 'CD161+ CD56- CD4 T', 'CD161+ CD56+ CD4 T'))

# TCRseq clonalQuant ---------------------------------------

# by gene
clonalQuant(TCR, cloneCall='gene', chain='both', scale=FALSE,  group.by='cd4tmemsubset')+
  scale_x_discrete(limits=subset_order)+
  scale_fill_manual(values=mycolors_cd4memorycelltypes)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))+
  labs(x='', fill='CD4 T cell subset')
ggsave(filename='./TCR/Plots/clonalQuant_gene_both_scale_no_by_cd4tmemsubset.pdf', height=80, width=100, units='mm')





# TCRseq clonalCompare --------------------------------------------------

# by gene
x <- clonalCompare(CD4TCR, cloneCall='gene', chain='TRA', top.clones=5, graph='alluvial', group.by='cd4tmemsubset')
write_excel_csv(x[['data']], file='clonalCompare_gene_TRA_top5clones_cd4tmemsubset.csv')
clones <- x[['data']] %>% filter(Sample=='CD161+ CD56+ CD4 T') %>% droplevels() %>% pull(clones)

unique_colors <- RColorBrewer::brewer.pal(n=length(clones), name="Set3")
clone_colors <- setNames(unique_colors, clones)

clonalCompare(CD4TCR, cloneCall='gene', chain='TRA', top.clones=5, graph='alluvial', group.by='cd4tmemsubset', relabel.clones=F, highlight.clones=clones)+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))+
  labs(x='', y='Proportion of clonal repertoire (%)', fill='Clones')+
  #scale_x_discrete(limits=subset_order)+
  scale_fill_manual(values=clone_colors)+
  scale_y_continuous(limits=c(0,0.11), breaks=c(0,0.02,0.04,0.06,0.08,0.1), labels=c(0,2,4,6,8,10), expand=c(0,0), guide='prism_offset')
ggsave(filename='./TCR/Plots/clonalCompare_gene_TRA_alluvial_by_cd4tmemsubset.pdf', height=55, width=100, units='mm')






# TCRseq clonalProportion -------------------------------------------------

x <- clonalProportion(CD4TCR, cloneCall='gene', chain='TRA', group.by='cd4tmemsubset')
write_excel_csv(x[['data']], file='./TCR/clonalProportion_gene_TRA_top5clones_cd4tmemsubset.csv')
y <- x[['data']] #%>% filter(Sample=='CD161+ CD56+ CD4 T') %>% droplevels() %>% pull(clones)

splitnumber <- c(1,10,20,50,100,1000)
clonalProportion(CD4TCR, cloneCall='gene', chain='TRA', group.by='cd4tmemsubset', clonalSplit=splitnumber)+
  #scale_y_continuous(limits=c(0,1.1), breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c(0,20,40,60,80,100), expand=c(0,0), guide='prism_offset')+
  labs(x='', y='Relative clonal abundance (%)')+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))
ggsave(filename='./TCR/Plots/clonalProportion_gene_TRA_by_cd4tmemsubset.pdf', height=60, width=80, units='mm')

clonalProportion(CD4TCR, cloneCall='gene', chain='both', group.by='cd4tmemsubset', clonalSplit=splitnumber)+
  #scale_y_continuous(limits=c(0,1.1), breaks=c(0,0.2,0.4,0.6,0.8,1), labels=c(0,20,40,60,80,100), expand=c(0,0), guide='prism_offset')+
  labs(x='', y='Relative clonal abundance (%)')+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8), axis.text.y=element_text(size=8), axis.title=element_text(size=8), legend.text=element_text(size=8), legend.title=element_text(size=8))
ggsave(filename='./TCR/Plots/clonalProportion_gene_both_by_cd4tmemsubset.pdf', height=60, width=80, units='mm')


# TCRseq percentGenes -----------------------------------------------------

# percent of genes detected across all samples
percentGenes(TCR, chain='TRA', gene='Vgene', group.by='cd4tmemsubset')+scale_y_discrete(limits=subset_order)+theme(text=element_text(size=8))
ggsave(filename='./TCR/Plots/percentGenes_TRAV_heatmap_by_cd4tmemsubset.pdf', height=50, width=150, units='mm')
pctgenesdata <- percentGenes(TCR, chain='TRA', gene='Vgene', group.by='cd4tmemsubset')
pctgenesdata <- write_excel_csv(as.data.frame(pctgenesdata[["data"]]), file='


                                percentGenes_TRAV_heatmap_by_cd4tmemsubset.csv')


# TCRseq clonalOverlap ----------------------------------------------------

clonalOverlap(CD4TCR, cloneCall='gene', chain='TRA', group.by='cd4tmemsubset', method='raw')+theme(axis.text.x=element_text(angle=45, hjust=1, size=8), axis.text.y=element_text(size=8), text=element_text(size=8))
ggsave(filename='./TCR/Plots/clonalOverlap_TRAV_rawclonenumber_similarity_by_cd4tmemsubset.pdf', height=80, width=100, units='mm')







# TCRseq Trex -------------------------------------------------------------

###
### Trex on TRA
###

CD4Trex <- runTrex(CD4TCR, chains="TRA", method='encoder', encoder.model="VAE", encoder.input="AF", reduction.name='Trex')
DimPlot(CD4Trex, reduction='Trex')

# Find TRA epitopes
cd4trex <- Trex::annotateDB(CD4Trex, chains="TRA", edit.distance=0) #Levenshtein distance of 0 == exact sequence as reported, <2 is optimal
DimPlot(cd4trex, reduction="rna_umap", group.by="TRA_Epitope.species")+theme_void()+theme(legend.text=element_text(size=8))+labs(title='')
ggsave(filename='./TCR/Plots/DimPlot_rnaUMAP_Trex_TRA_annotateDB_cd4tmemsubset.pdf', height=80, width=240, units='mm')

# UMAP for all virus-reactive T cells
unique(cd4trex$TRA_Epitope.species)
cd4trex$TRA_Epitope.species_clean <- gsub("Brain Cancer;", "", cd4trex$TRA_Epitope.species)
cd4trex$TRA_Epitope.species_clean <- gsub("Brain Cancer", NA, cd4trex$TRA_Epitope.species_clean)
cd4trex$TRA_Epitope.species_clean <- gsub(";", " ", cd4trex$TRA_Epitope.species_clean)
unique(cd4trex$TRA_Epitope.species_clean)


# UMAP of virus-reactive clones
virus_keywords <- c("CMV", "EBV", "HIV", "HCV", "Influenza", "SARS-CoV-2", "YFV") # from TRA_Epitope.species
cd4trex@meta.data <- cd4trex@meta.data %>% mutate(TRA_Epitope.species_clean_virus=case_when(is.na(TRA_Epitope.species_clean) ~ NA, # Preserve NAs
                     grepl(paste(virus_keywords, collapse = "|"), TRA_Epitope.species_clean) ~ "virus-reactive", TRUE~"non-virus-reactive"))

DimPlot(cd4trex, reduction="rna_umap", group.by="TRA_Epitope.species_clean_virus", cols=c('#8B0000', 'black'), na.value='lightgrey')+
  theme_void()+theme(legend.text=element_text(size=8))+labs(title='')
ggsave(filename='./TCR/Plots/DimPlot_rnaUMAP_Trex_TRA_annotateDB_virus_reactive.pdf', height=60, width=60, units='mm')

# UMAP of CMV-reactive clones
cd4trex@meta.data <- cd4trex@meta.data %>% mutate(TRA_Epitope.species_clean_CMV=case_when(is.na(TRA_Epitope.species_clean) ~ NA, # Preserve NAs
                     grepl('CMV', TRA_Epitope.species_clean)~'CMV-reactive', TRUE~'non-CMV-reactive'))
unique(cd4trex$TRA_Epitope.species_clean_CMV)

DimPlot(cd4trex, reduction="rna_umap", group.by="TRA_Epitope.species_clean_CMV", cols=c('black', '#8B0000'), na.value='lightgrey')+
  theme_void()+theme(legend.text=element_text(size=8))+labs(title='')
ggsave(filename='./TCR/Plots/DimPlot_rnaUMAP_Trex_TRA_annotateDB_CMV_reactive.pdf', height=60, width=60, units='mm')

# Reactivity summary table and plot
i <- 1
antigens <- c("CMV", "EBV", "InfluenzaA", "SARS-CoV-2", "SelfAg", "YFV", "HCV", "HIV")
antigen_table <- data.frame(Antigen=antigens, Count=0, x.axis='i')
for (antigen in antigens) {antigen_table$Count[antigen_table$Antigen==antigen] <- sum(grepl(antigen, cd4trex@meta.data$TRA_Epitope.species_clean, ignore.case=TRUE), na.rm=FALSE)}
antigen_table <- antigen_table %>%mutate(Antigen=fct_reorder(Antigen, Count))
write_excel_csv(antigen_table, file='Antigen_table_numberofreactiveclones.csv')

ggplot(antigen_table, aes(x=i, fill=Antigen, y=Count))+
  geom_bar(stat='identity', position='stack', color='black', linewidth=0.25)+
  scale_y_continuous(limits=c(0,264), breaks=c(0,80,160,240), expand=c(0,0), guide='prism_offset')+
  scale_fill_manual(values=c(brewer.pal(7, name='Greys'), 'darkred' ))+
  labs(x='', y='Number of reactive TRA clones')+
  theme(axis.text.x=element_blank(), panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Antigen_table_countsofreactiveclones.pdf', units='mm', width=55, height=38)

###
### Trex on TRB
###

CD4Trex_TRB <- runTrex(CD4TCR, chains="TRB", method='encoder', encoder.model="VAE", encoder.input="AF", reduction.name='Trex')
CD4Trex_TRB <- Trex::annotateDB(CD4Trex_TRB, chains="TRB", edit.distance=0)
DimPlot(CD4Trex_TRB, reduction="rna_umap", group.by="TRB_Epitope.species")+theme_void()+theme(legend.text=element_text(size=8))+labs(title='')

# add TRB annotation to TRA seurat object
cd4trex <- AddMetaData(cd4trex, metadata=CD4Trex_TRB@meta.data$TRB_Epitope.target, col.name='TRB_Epitope.target')
cd4trex <- AddMetaData(cd4trex, metadata=CD4Trex_TRB@meta.data$TRB_Epitope.sequence, col.name='TRB_Epitope.sequence')
cd4trex <- AddMetaData(cd4trex, metadata=CD4Trex_TRB@meta.data$TRB_Epitope.species, col.name='TRB_Epitope.species')
cd4trex <- AddMetaData(cd4trex, metadata=CD4Trex_TRB@meta.data$TRB_Database, col.name='TRB_Database')

# TCRseq virus-reactive clones --------------------------------------------

viral_identifiers <- c("CMV", "EBV", "InfluenzaA", "SARS-CoV-2", "HIV", "YFV", "HCV")
cd4trex@meta.data <- cd4trex@meta.data %>% mutate(virus_reactive=ifelse(is.na(TRA_Epitope.species), 'non-reactive',
                     ifelse(str_detect(TRA_Epitope.species, paste(viral_identifiers, collapse="|")), "virus-reactive", "non-reactive")))

subset_colors <- c('memory CD4 T'='black', 'CD161+ CD56- CD4 T'="#f4af2d", 'TOX+ TIGIT+ CD4 T'="#8b0000", 'CD161+ CD56+ CD4 T'= "#2f8ca3", 'CCR10+ memory CD4 T'="#FF0000",  'non-reactive'='grey')

# for each donor separately
virus_counts_per_donor <- cd4trex@meta.data %>% group_by(cd4tmemsubset, virus_reactive, donor) %>% summarise(Count_by_donor=n()) %>% ungroup() %>% # calculate proportion of reactivity for each donor
  group_by(donor, cd4tmemsubset) %>% mutate(proportion_in_subset=Count_by_donor*100 / sum(Count_by_donor)) %>%
  mutate(fill_color=ifelse(virus_reactive=='virus-reactive', yes=as.character(cd4tmemsubset), 'non-reactive')) %>% ungroup()

virus_counts_per_donor$cd4tmemsubset <- factor(virus_counts_per_donor$cd4tmemsubset, levels=c("memory CD4 T", 'CCR10+ memory CD4 T', 'TOX+ TIGIT+ CD4 T', 'CD161+ CD56- CD4 T', 'CD161+ CD56+ CD4 T'))
virus_counts_per_donor$fill_color <- factor(virus_counts_per_donor$fill_color, levels=c("memory CD4 T", 'CCR10+ memory CD4 T', 'TOX+ TIGIT+ CD4 T', 'CD161+ CD56- CD4 T', 'CD161+ CD56+ CD4 T', "non-reactive"))
write_excel_csv(virus_counts_per_donor, file='./TCR/virus_counts_per_donor_proportion_of_virus_reactive.csv')

# Plot proportions
ggplot(virus_counts_per_donor, aes(x=donor, y=proportion_in_subset, fill=fill_color))+
  geom_bar(stat='identity', position="stack", color='black', show.legend=FALSE)+
  scale_fill_manual(values=subset_colors)+
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100), expand=c(0,0), guide='prism_offset')+
  labs(title='', x="", y="Proportion of virus-reactive TCR clones (%)", fill="")+
  facet_grid(.~cd4tmemsubset)+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), legend.title=element_text(size=8), strip.background=element_blank(), strip.text=element_text(size=8),
        panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Barplot_proportion_per_donor_in_subset_virus_reactive_nonreactive_clones.pdf', height=80, width=120, units='mm')

virus_counts_per_donor %>% dplyr::filter(virus_reactive=='virus-reactive') %>%
  ggplot(aes(x=donor, y=proportion_in_subset, fill=fill_color))+
  geom_bar(stat='identity', position="stack", color='black', show.legend=FALSE)+
  scale_fill_manual(values=subset_colors)+
  scale_y_continuous(limits=c(0,17.6), breaks=c(0,4,8,12,16), expand=c(0,0), guide='prism_offset')+
  labs(title='', x="", y="Proportion of virus-reactive\nTCR clones (%)", fill="")+
  facet_grid(.~cd4tmemsubset)+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), legend.title=element_text(size=8), strip.background=element_blank(), strip.text=element_text(size=8),
        panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Barplot_proportion_in_subset_virus_reactive_clones.pdf', height=80, width=40, units='mm')

# control CD161+ CD56+ CD4 T vs all other memory CD4 T cell subsets
virus_counts_per_donor <- cd4trex@meta.data %>% mutate(cd161poscd56pos_vs_all=ifelse(cd4tmemsubset=='CD161+ CD56+ CD4 T', yes='CD161+ CD56+ CD4 T', no='other memory CD4 T')) %>%
  group_by(cd161poscd56pos_vs_all, virus_reactive, donor) %>% summarise(Count_by_donor=n()) %>% ungroup() %>% # calculate proportion of reactivity for each donor
  group_by(donor, cd161poscd56pos_vs_all) %>% mutate(proportion_in_subset=Count_by_donor*100 / sum(Count_by_donor)) %>%
  mutate(fill_color=ifelse(virus_reactive=='virus-reactive', yes=as.character(cd161poscd56pos_vs_all), 'non-reactive')) %>% ungroup()

ttest_cd161poscd56pos_vs_all <- virus_counts_per_donor %>% filter(virus_reactive=='virus-reactive') %>% t_test(proportion_in_subset~cd161poscd56pos_vs_all, p.adjust.method='none') %>% adjust_pvalue(method='BH') %>% add_significance(p.col='p', cutpoints=c(0,0.001,0.01,0.05,1), symbols=c('***','**','*','')) %>%
  add_xy_position(x='cd161poscd56pos_vs_all', fun='mean', scales='free') %>% filter(p<0.05)
write_excel_csv(ttest_cd161poscd56pos_vs_all, file='./TCR/ttest_cd161poscd56pos_vs_all_virus_reactive_per_donor_by_subset.csv')

virus_counts_per_donor %>% dplyr::filter(virus_reactive=='virus-reactive') %>%
  ggplot(aes(x=cd161poscd56pos_vs_all, y=proportion_in_subset))+
  geom_boxplot(aes(fill=fill_color), color='black', show.legend=FALSE)+
  geom_point(fill='white', shape=21, color='black')+
  scale_fill_manual(values=subset_colors)+
  scale_y_continuous(limits=c(0,17.6), breaks=c(0,4,8,12,16), expand=c(0,0), guide='prism_offset')+
  labs(title='', x="", y="Proportion of virus-reactive\nTCR clones (%)", fill="")+
  add_pvalue(ttest_cd161poscd56pos_vs_all, tip.length=0, bracket.size=0.35, size=3, y.position=c(16), label='p.signif')+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), legend.title=element_text(size=8), strip.background=element_blank(), strip.text=element_text(size=8),
        panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Boxplot_proportion_in_CD161+CD56+CD4T_vs_all_virus_reactive_clones_stats.pdf', height=68, width=25, units='mm')

# TCRseq CMV-reactive clones --------------------------------------------

cd4trex@meta.data <- cd4trex@meta.data %>% mutate(CMV_reactive=ifelse(is.na(TRA_Epitope.species), 'non-CMV-reactive',
                     ifelse(str_detect(TRA_Epitope.species, "CMV"), "CMV-reactive", "non-CMV-reactive")))

# for each donor separately
CMV_counts_per_donor <- cd4trex@meta.data %>% group_by(cd4tmemsubset, CMV_reactive, donor) %>% summarise(Count_by_donor=n()) %>% ungroup() %>% # calculate proportion of reactivity for each donor
  group_by(donor, cd4tmemsubset) %>% mutate(proportion_in_subset=Count_by_donor*100 / sum(Count_by_donor)) %>%
  mutate(fill_color=ifelse(CMV_reactive=='CMV-reactive', yes=as.character(cd4tmemsubset), 'non-CMV-reactive')) %>% ungroup()

CMV_counts_per_donor$cd4tmemsubset <- factor(CMV_counts_per_donor$cd4tmemsubset, levels=c("memory CD4 T", 'CCR10+ memory CD4 T', 'TOX+ TIGIT+ CD4 T', 'CD161+ CD56- CD4 T', 'CD161+ CD56+ CD4 T'))
CMV_counts_per_donor$fill_color <- factor(CMV_counts_per_donor$fill_color, levels=c("memory CD4 T", 'CCR10+ memory CD4 T', 'TOX+ TIGIT+ CD4 T', 'CD161+ CD56- CD4 T', 'CD161+ CD56+ CD4 T', "non-CMV-reactive"))
write_excel_csv(CMV_counts_per_donor, file='./TCR/CMV_counts_per_donor_proportion_of_CMV_reactive.csv')

ggplot(CMV_counts_per_donor, aes(x=donor, y=proportion_in_subset, fill=fill_color))+
  geom_bar(stat='identity', position="stack", color='black', show.legend=FALSE)+
  scale_fill_manual(values=subset_colors)+
  scale_y_continuous(limits=c(0,100), breaks=c(0,20,40,60,80,100), expand=c(0,0), guide='prism_offset')+
  labs(title='', x="", y="Proportion of CMV-reactive TCR clones (%)", fill="")+
  facet_grid(.~cd4tmemsubset)+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), legend.title=element_text(size=8), strip.background=element_blank(), strip.text=element_text(size=8),
        panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Barplot_proportion_per_donor_in_subset_CMV_reactive_nonreactive_clones.pdf', height=80, width=120, units='mm')

# CD161+ CD56+ CD4 T vs all other memory CD4 T cell subsets
CMV_counts_per_donor <- cd4trex@meta.data %>% mutate(cd161poscd56pos_vs_all=ifelse(cd4tmemsubset=='CD161+ CD56+ CD4 T', yes='CD161+ CD56+ CD4 T', no='other memory CD4 T')) %>%
                        group_by(cd161poscd56pos_vs_all, CMV_reactive, donor) %>% summarise(Count_by_donor=n()) %>% ungroup() %>% # calculate proportion of reactivity for each donor
                        group_by(donor, cd161poscd56pos_vs_all) %>% mutate(proportion_in_subset=Count_by_donor*100 / sum(Count_by_donor)) %>%
                        mutate(fill_color=ifelse(CMV_reactive=='CMV-reactive', yes=as.character(cd161poscd56pos_vs_all), 'non-reactive')) %>% ungroup()

ttest_cd161poscd56pos_vs_all <- CMV_counts_per_donor %>% filter(CMV_reactive=='CMV-reactive') %>% t_test(proportion_in_subset~cd161poscd56pos_vs_all, p.adjust.method='none') %>% adjust_pvalue(method='BH') %>% add_significance(p.col='p', cutpoints=c(0,0.001,0.01,0.05,1), symbols=c('***','**','*','')) %>%
                                add_xy_position(x='cd161poscd56pos_vs_all', fun='mean', scales='free') %>% filter(p<0.05)
write_excel_csv(ttest_cd161poscd56pos_vs_all , file='./TCR/ttest_cd161poscd56pos_vs_all_CMV_reactive_per_donor_by_subset.csv')

CMV_counts_per_donor %>% dplyr::filter(CMV_reactive=='CMV-reactive') %>%
  ggplot(aes(x=cd161poscd56pos_vs_all, y=proportion_in_subset))+
  geom_boxplot(aes(fill=fill_color), color='black', show.legend=FALSE)+
  geom_point(fill='white', shape=21, color='black')+
  scale_fill_manual(values=subset_colors)+
  scale_y_continuous(limits=c(0,17.6), breaks=c(0,4,8,12,16), expand=c(0,0), guide='prism_offset')+
  labs(title='', x="", y="Proportion of CMV-reactive\nTCR clones (%)", fill="")+
  add_pvalue(ttest_cd161poscd56pos_vs_all, tip.length=0, bracket.size=0.35, size=3, y.position=c(16), label='p.signif')+
  theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8), legend.text=element_text(size=8),
        axis.title=element_text(size=8), legend.title=element_text(size=8), strip.background=element_blank(), strip.text=element_text(size=8),
        panel.grid.major.y= element_line(color='#d9d9d9', linetype='dotted', linewidth=0.25))
ggsave(filename='./TCR/Plots/Boxplot_proportion_in_CD161+CD56+CD4T_vs_all_CMV_reactive_clones_stats.pdf', height=68, width=25, units='mm')


# Data export ------------------------------------------------------------------

cd4trex_metadata <- cd4trex@meta.data
cd4trex_metadata[is.na(cd4trex_metadata)] <- 'n.a.' # annotate non-matched epitope (n.a.) different to not tested (NA)
int.TNK <- AddMetaData(int.TNK, metadata=cd4trex_metadata)
saveRDS(int.TNK, file='int.TNK')
