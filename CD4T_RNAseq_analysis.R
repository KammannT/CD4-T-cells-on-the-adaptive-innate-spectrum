# Short summary of experimentation:
# PBMCs were isolated from four healthy human blood donors.
# PBMCs were stained with 5OPRU-hMR1-Tet and positively enriched before sorting of MAIT (Hashtag1) and non-MAIT cells (Hashtag3).
# Libraries were generated following the 10X Genomics Next-GEM v2 5' Kit with Feature Barcode Technology and TCR sequencing.

# Setup -------------------------------------------------------------------

devtools::install:github('plger/scDblFinder') # DOI: 10.12688/f1000research.73600.2
devtools::install:github('powellgenomicslab/Nebulosa') # DOI: 10.1093/bioinformatics/btab003
devtools::install_github("BorchLab/Trex") # DOI: 10.1038/s41590-024-01888-9, Vignette: https://www.borch.dev/uploads/screpertoire/articles/trex
libraries <- c('scDblFinder', 'Nebulosa', 'Trex')

required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix",
                       "scales", "stats", "RColorBrewer", "ggpubr",
                       "glmGamPoi", "UCell", "org.Hs.eg.db", "Seurat", "HGNChelper",
                       "ReactomePA", "msigdbr", "DOSE", "fgsea", "clusterProfiler",  "MAST", "tidyverse")

for(package in required_packages){if(!require(package, character.only=TRUE))
{install.packages(package, dependencies=TRUE)
  library(package, character.only=TRUE) }}

lapply(libraries, library, character.only=TRUE)
rm(required_packages, libraries)

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

mycolors_cd4memorycelltypes <- c('memory CD4 T'='#bdbdbd',
                                 'CD161+ CD56- CD4 T'='#f4af2d',
                                 'TOX+ TIGIT+ CD4 T'='#8b0000',
                                 'CD161+ CD56+ CD4 T'='#2f8ca3',
                                 'CCR10+ memory CD4 T'='#FF0000')


# TNK: PCA and UMAP post filtering ---------------------------------------

int.TNK <- readRDS('./Seurat_data/Data deposition/int.TNK.rds')
DefaultAssay(int.TNK) <- 'RNAint'
int.TNK <- ScaleData(int.TNK)
int.TNK <- RunPCA(int.TNK, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
int.TNK <- ProjectDim(int.TNK, reduction='rna_pca')
int.TNK <- RunUMAP(int.TNK, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:10)
int.TNK <- FindNeighbors(int.TNK, dims=1:10, reduction='rna_pca')
int.TNK <- FindClusters (int.TNK, resolution=1)
Idents(int.TNK) <- int.TNK@meta.data[['RNAint_snn_res.1']]

DimPlot(int.TNK, pt.size=0.2, reduction='rna_umap', label=TRUE)+NoLegend()
ggsave(filename='UMAP/TNK/DimPlot_rnaUMAP_TNK.png', path=plotdir, width=100, height=100, units='mm')

DimPlot(int.TNK, pt.size=0.2, reduction='rna_umap', group.by='sctype', shuffle=TRUE, label=TRUE)+NoLegend()
ggsave(filename='UMAP/TNK/DimPlot_rnaUMAP_TNK_by_sctype.png', path=plotdir, width=100, height=100, units='mm')

# DEG analysis
DefaultAssay(int.TNK) <- 'RNA'
int.TNK <- ScaleData(int.TNK, assay='RNA')
int.TNK_rna_markers <- FindAllMarkers(int.TNK, assay='RNA', test.use='wilcox', only.pos=TRUE, min.pct=0.2, logfc.threshold=0.2)

# Identify CD4+ T cells
plot_density(int.TNK, features='adt_CD45RA', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_adt_CD45RA.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(int.TNK, features='adt_CD45RO', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_adt_CD45RO.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(int.TNK, features='rna_CD8A', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_rna_CD8A.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(int.TNK, features='adt_CD8', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_adt_CD8.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(int.TNK, features='rna_CD4', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_rna_CD4.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(int.TNK, features='adt_CD4', reduction='rna_umap', combine=TRUE, joint=FALSE)
ggsave(filename='/Features/TNK/Density_adt_CD4.pdf', path=plotdir, width=80, height=80, units='mm')

colors_clusters <- c('CD4T cells'='#608080', 'NK cells'='#472D7BFF',
                     'regulatory CD4 T cells'='#3E356BFF','naive CD4 T cells'='#8BBEDA', 'memory CD4 T cells'='#3B5698',
                     'naive CD8 T cells'='#F6BB97FF', 'memory CD8 T cells'='#F4875EFF', 'cytotoxic CD8 T cells'='#CB1B4FFF',
                     'gd T cells'='#96DDB5FF')

DimPlot(int.TNK, pt.size=0.2, reduction='rna_umap', cols=colors_clusters, label=FALSE, shuffle=TRUE)
ggsave(filename='UMAP/TNK/DimPlot_rnaUMAP_TNK_celltype_labelled.pdf', path=plotdir, width=100, height=100, units='mm')

# CD4+ T cell sub-clustering ------------------------------------------------------

# Isolate CD4+ T cells from int.TNK and test for contaminating cells
cd4cluster <- subset(int.TNK, celltype=='regulatory CD4 T cells' | celltype=='naive CD4 T cells' | celltype=='memory CD4 T cells')

FeatureScatter(cd4cluster, feature1='hto_Hashtag1', feature2='rna_TRAV1-2')+geom_vline(xintercept=0.5)+geom_hline(yintercept=0.5)+NoLegend()
FeatureScatter(cd4cluster, feature1='rna_CD4', feature2='rna_CD8A')+geom_vline(xintercept=0.5)+geom_hline(yintercept=0.5)+NoLegend()
FeatureScatter(cd4cluster, feature1='rna_CD4', feature2='rna_CD8B')+geom_vline(xintercept=0.5)+geom_hline(yintercept=0.5)+NoLegend()
FeatureScatter(cd4cluster, feature1='adt_CD4', feature2='adt_CD8')+geom_vline(xintercept=0.5)+geom_hline(yintercept=1)+NoLegend()

# filter contaminating cells from CD4+ T cells
cd4s <- subset(cd4cluster, hto_Hashtag1<0.5 & `rna_TRAV1-2`<0.5 & rna_CD8A<0.5 & rna_CD8B<0.5 & adt_CD8<1 & adt_CD4>0.5)

# Re-Integrate CD4+ T cells
cd4s.list <- SplitObject(cd4s, split.by='donor')
DefaultAssay(cd4s) <- 'RNA'
cd4s.list     <- lapply(cd4s.list, FUN=function(x) {x <- NormalizeData(object=x, assay='RNA') } )
cd4s.list     <- lapply(cd4s.list, FUN=function(x) {x <- FindVariableFeatures(object=x, assay='RNA', nfeatures=5000, selection.method='vst') } )
cd4s.list     <- lapply(cd4s.list, FUN=function(x) {x <- ScaleData(object=x, assay='RNA') } )
cd4s.features <- SelectIntegrationFeatures(object.list=cd4s.list, nfeatures=5000)
cd4s.anchors  <- FindIntegrationAnchors(object.list=cd4s.list, anchor.features=cd4s.features)
cd4s          <- IntegrateData(anchorset=cd4s.anchors, normalization.method='LogNormalize', new.assay.name='RNAint')

DefaultAssay(cd4s) <- 'RNAint'
cd4s <- Trex::quietTCRgenes(cd4s) # silence TCR chain genes from VariableFeatures present in RNAint assay
cd4s <- ScaleData(cd4s)
cd4s <- RunPCA(cd4s, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
cd4s <- ProjectDim(cd4s, reduction='rna_pca')
cd4s <- RunUMAP(cd4s, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:10)
cd4s <- FindNeighbors(cd4s, dims=1:10, reduction='rna_pca', )
cd4s <- FindClusters (cd4s, resolution=0.8, algorithm=4) # 4=Leiden clustering

DimPlot(cd4s, pt.size=0.5, reduction='rna_umap', label=TRUE)+theme_void()+NoLegend()
ggsave(filename='/UMAP/Density_CD4s/Dimplot_cd4s_celltype.pdf', path=plotdir, width=80, height=80, units='mm')

DimPlot(cd4s, pt.size=0.5, reduction='rna_umap', group.by='RNAint_snn_res.0.8', label=TRUE)+theme_void()+NoLegend()
ggsave(filename='/UMAP/Density_CD4s/Dimplot_cd4s_cluster.pdf', path=plotdir, width=80, height=80, units='mm')

# Identify memory CD4+ T cells
DefaultAssay(cd4s) <- 'ADT'
cd4s <- ScaleData(cd4s)
adt    <- c('CD4_TotalSeqC', 'CD8_TotalSeqC', 'CD16_TotalSeqC', 'CD27_TotalSeqC', 'CD28_TotalSeqC', 'CD45RA_TotalSeqC', 'CD45RO_TotalSeqC', 'CD56_TotalSeqC', 'CD62L_TotalSeqC', 'CD161_TotalSeqC', 'CCR7_TotalSeqC', 'Va72_TotalSeqC')
adt_features <- sub('_TotalSeqC', '', adt)

for (i in adt_features){
  plot_density(cd4s, features=paste0('adt_',i), reduction='rna_umap', combine=TRUE, joint=FALSE)+NoLegend()
  ggsave(filename=paste0('/UMAP/Density_CD4s/Features_CD4s_rnaUMAP_adt_',i,'.pdf'), path=plotdir, width=80, height=80, units='mm')
}

plot_density(cd4s, features='rna_FOXP3', reduction='rna_umap')# InnateCD4 are high
ggsave(filename='/UMAP/Density_CD4s/FeaturesDensity_cd4s_legend.pdf', path=plotdir, width=80, height=80, units='mm')

# DEG analysis
DefaultAssay(cd4s) <- 'RNA'
cd4s <- ScaleData(cd4s, assay='RNA')
cd4s_rna_markers <- FindAllMarkers(cd4s,  assay='RNA', test.use='wilcox', only.pos=TRUE, min.pct=0.1, logfc.threshold=0.2)

# annotate cell clusters
Idents(cd4s) <- cd4s@meta.data[["RNAint_snn_res.0.8"]]
DimPlot(cd4s, pt.size=0.5, reduction='rna_umap', label=TRUE)
cd4s <- RenameIdents(cd4s, '1'='memory T', '2'='naive T', '3'='naive T', '4'='memory T', '5'='naive T', '6'='naive T', '7'='memory T', '8'='memory T', '9'='Treg', '10'='Treg')
cd4s@meta.data$CD4Tphenotype <- Idents(cd4s)

# memory CD4+ T cell sub-clustering ---------------------------------------

cd4tmems <- subset(cd4s, CD4Tphenotype=='memory T')

cd4tmems.list <- SplitObject(cd4tmems, split.by='donor')
DefaultAssay(cd4tmems) <- 'RNA'
cd4tmems.list     <- lapply(cd4tmems.list, FUN=function(x) {x <- NormalizeData(object=x, assay='RNA') } )
cd4tmems.list     <- lapply(cd4tmems.list, FUN=function(x) {x <- FindVariableFeatures(object=x, assay='RNA', nfeatures=5000, selection.method='vst') } )
cd4tmems.list     <- lapply(cd4tmems.list, FUN=function(x) {x <- ScaleData(object=x, assay='RNA') } )
cd4tmems.features <- SelectIntegrationFeatures(object.list=cd4tmems.list, nfeatures=5000)
cd4tmems.anchors  <- FindIntegrationAnchors(object.list=cd4tmems.list, anchor.features=cd4tmems.features)
cd4tmems          <- IntegrateData(anchorset=cd4tmems.anchors, normalization.method='LogNormalize', new.assay.name='RNAint')
rm(cd4tmems.anchors, cd4tmems.list)

DefaultAssay(cd4tmems) <- 'RNAint'
cd4tmems <- Trex::quietTCRgenes(cd4tmems) # silence TCR chain genes from VariableFeatures in RNAint assay
cd4tmems <- ScaleData(cd4tmems)
cd4tmems <- RunPCA(cd4tmems, npcs=50, reduction.name='rna_pca', reduction.key='rnaPC_')
cd4tmems <- ProjectDim(cd4tmems, reduction='rna_pca')
ElbowPlot(cd4tmems, ndims=50, reduction='rna_pca')
ggsave(filename='PCA/ElbowPlot_RNAint_cd4tmems_PCA.png', path=plotdir, width=80, height=80, units='mm')
cd4tmems <- RunUMAP(cd4tmems, reduction='rna_pca', reduction.name='rna_umap', reduction.key='rnaUMAP_', dims=1:20)
cd4tmems <- FindNeighbors(cd4tmems, dims=1:20, reduction='rna_pca')
cd4tmems <- FindClusters (cd4tmems, resolution=0.5, algorithm=4)

DimPlot(cd4tmems, pt.size=0.5, reduction='rna_umap', label=TRUE)

# Features of memory CD4+ T cells
for (i in adt_features){
  plot_density(cd4tmems, features=paste0('adt_',i), reduction='rna_umap', combine=TRUE, joint=FALSE)+NoLegend()
  ggsave(filename=paste0('/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_adt_',i,'.pdf'), path=plotdir, width=80, height=80, units='mm')
}

plot_density(cd4tmems, features=c('adt_CD56', 'adt_CD161'), joint=TRUE, reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_adt_CD56_adt_CD161.pdf', path=plotdir, width=260, height=80, units='mm')

plot_density(cd4tmems, features='rna_TOX', reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_rna_TOX.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(cd4tmems, features='rna_TIGIT', reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_rna_TIGIT.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(cd4tmems, features='rna_LYAR', reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_rna_LYAR.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(cd4tmems, features='rna_TBX21', reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_rna_TBX21.pdf', path=plotdir, width=80, height=80, units='mm')

plot_density(cd4tmems, features='rna_EOMES', reduction='rna_umap')+NoLegend()
ggsave(filename='/UMAP/Density_cd4tmems/FeaturesDensity_cd4tmems_rnaUMAP_rna_EOMES.pdf', path=plotdir, width=80, height=80, units='mm')

# DEG analysis in memory cd4tmem ------------------------------------------------------------

DefaultAssay(cd4tmems) <- 'RNA'
cd4tmems <- ScaleData(cd4tmems, assay='RNA')
cd4tmems_rna_markers_mast <- FindAllMarkers(cd4tmems, assay='RNA', test.use='MAST', latent.vars='donor', only.pos=FALSE, min.pct=0.1, logfc.threshold=0, return.thresh=1)
cd4tmems_rna_markers_mast <- cd4tmems_rna_markers_mast %>% adjust_pvalue(method='BH', p.col='p_val', output.col='adjust_p') %>% arrange(cluster, desc(avg_log2FC))
write_excel_csv(cd4tmems_rna_markers_mast, file='./Tables/cd4tmems_rna_markers_mast.csv')

top10_cd4tmems_rna_markers <- cd4tmems_rna_markers_mast %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% mutate(rank=base::rank(dplyr::desc(avg_log2FC), ties.method='random')) %>% arrange(cluster, rank, avg_log2FC)
cd4tmems %>% subset(downsample=200) %>% DoHeatmap(features=top10_cd4tmems_rna_markers$gene, assay='RNA', slot='scale.data', size=3, group.bar=TRUE) & NoLegend()
ggsave(filename='/Heatmaps/Top10_RNA_Heatmap_MAST_cd4tmems.png', path=plotdir, width=80, height=220, units='mm')
ggsave(filename='/Heatmaps/Top10_RNA_Heatmap_MAST_cd4tmems.pdf', path=plotdir, width=80, height=220, units='mm')

cd4tmems %>% subset(downsample=200) %>% DoHeatmap(features=cd4tmems_adt_markers$gene, assay='ADT', slot='scale.data', size=4, group.bar=TRUE) & NoLegend()
ggsave(filename='/Heatmaps/ADT_Heatmap_cd4tmems.png', path=plotdir, width=80, height=220, units='mm')

### compute cluster 4 against all other

DefaultAssay(cd4tmems) <- 'RNA'
genelabels <- cd4tmems_rna_markers_mast %>% filter(cluster=='4', adjust_p<0.001) %>% filter(avg_log2FC>0.5 | avg_log2FC<(-0.5)) %>% pull(gene)
x <- cd4tmems_rna_markers_mast %>% filter(cluster=='4')

ggplot(data=x, aes(x=avg_log2FC, y=-log10(adjust_p)))+
  geom_hline(yintercept=-log10(0.001), linetype="dotted", color='black')+
  geom_vline(xintercept=c(-0.5, 0.5), linetype="dotted", color='black')+
  ggrepel::geom_label_repel(data=subset(x, gene %in% genelabels), aes(label=gene), box.padding=0.2, label.padding=0.15, size=2, show.legend=FALSE)+
  geom_point(fill=ifelse(x$avg_log2FC>0.5, yes='#2f8ca3', no='#bdbdbd'), shape=21, size=1, show.legend=FALSE)+
  scale_x_continuous(limits=c(-6.6, 6.6), breaks=c(-6,-4,-2,0,2,4,6), guide='prism_offset', expand=c(0,0) )+
  scale_y_continuous(limits=c(0,172), breaks=c(0,40,80,120,160), guide='prism_offset', expand=c(0.02, 0))+
  theme(panel.grid.major.y=element_blank())+
  labs(x='Log2FC (CD161+ CD56+ CD4 T / memory CD4 T)', y='-log10(adjusted p-value)', fill='', title='')
ggsave(filename='./Plots/Volcano/VolcanoPlot_cd4tmems_CD161+CD56+CD4T_labelled_full.png', width=250, height=200, units='mm')
ggsave(filename='./Plots/Volcano/VolcanoPlot_cd4tmems_CD161+CD56+CD4T_labelled_full.pdf', width=250, height=200, units='mm')

# annotation of memory CD4+ T cells
cd4tmems <- RenameIdents(cd4tmems, '1'='memory CD4 T', '2'='memory CD4 T', '3'='TOX+ TIGIT+ CD4 T', '4'='CD161+ CD56+ CD4 T', '5'='CD161+ CD56- CD4 T', '6'='CCR10+ memory CD4 T')
cd4tmems@meta.data$cd4tmemsubset <- Idents(cd4tmems)

DimPlot(cd4tmems, pt.size=0.5, reduction='rna_umap', label=TRUE, group.by='cd4tmemsubset', cols=c('#bdbdbd', '#8b0000', '#2f8ca3', '#f4af2d','#FF0000'))
ggsave(filename='UMAP/cd4tmems/UMAP_rnaUMAP_cd4tmems_cd4tmemsubset_withlegend.pdf', path=plotdir, width=80, height=80, units='mm')

# GO Overrepresentation analysis (ORA) -----------------------------------------------------------

# DEG between CD161+CD56+ CD4 T cells and all other memory CD4 T cells
cd4tmems$cd4tmem_memory <- cd4tmems$cd4tmemsubset
Idents(cd4tmems) <- cd4tmems$cd4tmem_memory
cd4tmems <- RenameIdents(cd4tmems, 'memory CD4 T'='memory CD4 T',  'TOX+ TIGIT+ CD4 T'='memory CD4 T', 'CD161+ CD56+ CD4 T'='CD161+ CD56+ CD4 T', 'CD161+ CD56- CD4 T'='memory CD4 T', 'CCR10+ memory CD4 T'='memory CD4 T') # collapse factor
cd4tia <- subset(cd4tmems, cd4tmem_memory %in% c('CD161+ CD56+ CD4 T', 'memory CD4 T'))

# DEG between CD161+ CD56+ CD4+ T cells and all other memory CD4+ T cells
DefaultAssay(cd4tia) <- 'RNA'
markers_cd4tia <- FindAllMarkers(cd4tia, assay='RNA', test.use='MAST', latent.vars='donor', only.pos=FALSE, min.pct=0.1, logfc.threshold=0, return.thresh=1)
markers                <- markers_cd4tia %>% filter(avg_log2FC<(-0.5) | avg_log2FC>0.5)
markers_list           <- markers
markers_list$SYMBOL    <- markers_list$gene
markers_list           <- bitr(markers_list$SYMBOL, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db', drop=T)
markers_list$gene      <- markers_list$SYMBOL
markers                <- left_join(markers, markers_list)

# create universe for ORA GO analysis == all genes expressed by cd4tmem cells
cd4tia_ave     <- AverageExpression(cd4tia, return.seurat=TRUE, group.by='cd4tmem_memory')
cd4tia_ave_rna <- GetAssayData(cd4tia_ave, assay='RNA', slot='data')
cd4tia_ave_rna <- as.data.frame(cd4tia_ave_rna)
cd4tia_ave_rna$gene <- rownames(cd4tia_ave_rna)

universe_background <- cd4tia_ave_rna$gene
ensembl             <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
symbol_to_entrez    <- getBM(values=universe_background, attributes=c('external_gene_name', 'ensembl_gene_id', 'entrezgene_id'), filters='external_gene_name', mart=ensembl)
universe_background <- drop_na(symbol_to_entrez)
colnames(universe_background) <- c('SYMBOL', 'ENSEMBL', 'ENTREZID')

top100  <- markers %>% dplyr::filter(adjust_p<0.05 & avg_log2FC>0.5 & cluster=='CD161+ CD56+ CD4 T') %>% top_n(n=100, wt=avg_log2FC) %>% arrange(desc(avg_log2FC))
top100  <- split(top100$ENTREZID, top100$cluster)

CD161CD56cd4tmem_go        <- enrichGO(ont='BP', gene=top100[['CD161+ CD56+ CD4 T']], universe=as.character(universe_background$ENTREZID), OrgDb=org.Hs.eg.db, pAdjustMethod='BH', pvalueCutoff=0.05, readable=TRUE)
CD161CD56cd4tmem_go_result <- CD161CD56cd4tmem_go@result %>% arrange(desc(pvalue)) %>% filter(Count>=5) %>% filter(pvalue<0.05)
write_excel_csv(CD161CD56cd4tmem_go_result, file='./Tables/GO/CD161CD56cd4tmem_go_result.csv')

selected_ora <- c('leukocyte migration', 'positive regulation of cytokine production', 'lymphocyte proliferation',
                  'response to virus', 'cell killing', 'cell-matrix adhesion',
                  'regulation of metal ion transport', 'regulation of transmembrane transport')

CD161CD56cd4tmem_go_result_selected <- CD161CD56cd4tmem_go_result %>% filter(Description %in% selected_ora) %>% mutate(CountinSet=Count*100/91)

ggplot(CD161CD56cd4tmem_go_result_selected, aes(y=reorder(Description, CountinSet), x=CountinSet))+
  geom_point(aes(fill=pvalue), shape=21, color='black', show.legend=TRUE, size=3)+
  scale_fill_viridis()+
  scale_x_continuous(limits=c(4,12.2), breaks=c(4,8,12), expand=c(0,0.5))+
  scale_y_discrete(expand=c(0.02,0.05), guide='prism_offset')+
  labs(x='gene count in gene set (%)', y='', fill='p-value')+
  theme(panel.grid.major.y=element_line(linewidth=0.3, linetype='dotted'),
        axis.text=element_text(size=8, face='plain'),
        axis.title=element_text(size=8, face='plain'))
ggsave(filename='./Plots/ORA/ORA_GO_CD161+CD56+cd4tmemcells.pdf', width=120, height=40, units='mm')

# GO Gene set enrichment analysis (GSEA) ---------------------------------------

# Hallmark gene sets for GSEA
hallmark_human               <- msigdbr(species='Homo sapiens', category='H')
Hallmark_human_list_ENTREZID <- hallmark_human %>% dplyr::select(gs_name, entrez_gene) %>% group_by(gs_name) %>% summarize(all.genes=list(unique(entrez_gene))) %>% deframe()
Hallmark_human_list_SYMBOL   <- hallmark_human %>% dplyr::select(gs_name, gene_symbol) %>% group_by(gs_name) %>% summarize(all.genes=list(unique(gene_symbol))) %>% deframe()

CD161CD56cd4tmem_vs_memorycd4tmem <- cd4tia_ave_rna %>% mutate(SYMBOL=gene) %>% dplyr::select(-gene) %>%
                                     dplyr::mutate(delta_avg_exp=`CD161+ CD56+ CD4 T`-`memory CD4 T`) %>% arrange(desc(delta_avg_exp))

CD161CD56cd4tmem_vs_memorycd4tmem_delta  <- CD161CD56cd4tmem_vs_memorycd4tmem %>% dplyr::select(SYMBOL, delta_avg_exp) %>% arrange(desc(delta_avg_exp)) #always order from highest NES descending!
any(duplicated(CD161CD56cd4tmem_vs_memorycd4tmem_delta$SYMBOL)) # no duplicates
CD161CD56cd4tmem_vs_memorycd4tmem_delta_FC <- CD161CD56cd4tmem_vs_memorycd4tmem_delta$delta_avg_exp
names(CD161CD56cd4tmem_vs_memorycd4tmem_delta_FC) <- CD161CD56cd4tmem_vs_memorycd4tmem_delta$SYMBOL # necessary for fgsea

# Hallmark GSEA
gsea        <- clusterProfiler::GSEA(geneList=CD161CD56cd4tmem_vs_memorycd4tmem_delta_FC, eps=0, pvalueCutoff=0.5, pAdjustMethod='BH', TERM2GENE=dplyr::select(hallmark_human, gs_name, gene_symbol), seed=42)
gsea_result <- data.frame(gsea@result) %>% arrange(desc(NES), p.adjust)
write_csv(gsea_result, file='./Tables/GSEA/gsea.Hallmark_CD161+CD56+CD4+T_vs_memorycd4tmem.csv')

gsea_result %>% filter(p.adjust<0.05) %>% mutate(Description=gsub('HALLMARK_', "", Description), Description=gsub('_', " ", Description)) %>%
  ggplot(aes(x=reorder(Description, NES), y=NES))+
  geom_col(color='black', fill='darkgrey', show.legend=FALSE)+
  coord_flip()+
  labs(x='Gene Set', y='NES')+
  theme(axis.text=element_text(size=10, face='plain'), axis.title=element_text(size=10, face='plain'))
ggsave(file='./Plots/GSEA/Hallmark_GSEA_CD161+CD56+CD4Tmems.png', width=120, height=40, units='mm')
rm(gsea, gsea_result, hallmark_human, Hallmark_human_list_ENTREZID, Hallmark_human_list_SYMBOL)

# GO GSEA
# BP
gse <- gseGO(geneList=CD161CD56cd4tmem_vs_memorycd4tmem_delta_FC, ont='BP',
             OrgDb='org.Hs.eg.db', keyType='SYMBOL', minGSSize=10, maxGSSize=500,
             pvalueCutoff=0.05, pAdjustMethod='BH', seed=42)

gse_result_bp <- gse@result
write_delim(gse_result_bp, file='./Tables/GSEA/GSEA_GO_BP_result.csv')

gse_result_bp_selected <- gse_result_bp %>% filter(p.adjust<0.05) %>% filter(Description %in% c('disruption of cell in another organism', 'cell killing', 'positive chemotaxis', 'regulation of metal ion transport'))
ggplot(gse_result_bp_selected, aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=p.adjust), color='black', show.legend=TRUE)+
  coord_flip()+
  scale_fill_viridis()+
  scale_y_continuous(limits=c(0, abs(max(gse_result$NES)+0.2)))+
  labs(x='Gene Set', y='NES')
ggsave(file='./Plots/GSEA/gseGO_BP_CD161+CD56+cd4tmem_vs_memorycd4tmem.pdf', width=120, height=50, units='mm')

# MF
gse <- gseGO(geneList=CD161CD56cd4tmem_vs_memorycd4tmem_delta_FC, ont='MF',
             OrgDb='org.Hs.eg.db', keyType='SYMBOL', minGSSize=10, maxGSSize=500,
             pvalueCutoff=0.05, pAdjustMethod='BH', seed=42)
gse_result_MF <- gse@result
write_delim(gse_result_bp, file='./Tables/GSEA/GSEA_GO_MF_result.csv')
gse_result_MF_selected <- gse_result_MF %>% filter(p.adjust<0.05) %>% filter(Description %in% c('chemokine activity', 'cytokine activity', 'endopeptidase activity', 'neurotransmitter receptor regulator activity'))
ggplot(gse_result_MF_selected, aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=p.adjust), color='black', show.legend=TRUE)+
  coord_flip()+
  scale_fill_viridis()+
  scale_y_continuous(limits=c(0, abs(max(gse_result$NES)+0.2)))+
  labs(x='Gene Set', y='NES')
ggsave(file='./Plots/GSEA/gseGO_MF_CD161+CD56+cd4tmem_vs_memorycd4tmem.pdf', width=120, height=50, units='mm')

# combined visualization
gse_result_merged <- rbind(gse_result_MF_selected, gse_result_bp_selected)
ggplot(gse_result_merged, aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=pvalue), color='black', show.legend=TRUE)+
  coord_flip()+
  scale_fill_viridis()+
  scale_y_continuous(limits=c(0, abs(max(gse_result$NES)+0.2)))+
  labs(x='Gene Set', y='NES')
ggsave(file='../5.8 InnateCD4seq/Plots/GSEA/gseGO_BP_and_MF_CD161+CD56+cd4tmem_vs_memorycd4tmem.pdf', width=120, height=50, units='mm')

# Reactome pathway gene set enrichment analysis -------------------------

universe_background <- cd4tia_ave_rna$gene %>% as.data.frame()
universe_background$SYMBOL <- universe_background$.
universe_background$. <- NULL
markers_list        <- clusterProfiler::bitr(universe_background$SYMBOL, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db', drop=T)
markers_list$gene   <- markers_list$SYMBOL

CD161CD56cd4tmems_vs_memorycd4tmem <- cd4tia_ave_rna %>% mutate(SYMBOL=gene) %>% dplyr::select(-gene) %>%
                                      dplyr::mutate(delta_avg_exp=`CD161+ CD56+ CD4 T`-`memory CD4 T`) %>% arrange(desc(delta_avg_exp))

CD161CD56cd4tmems_vs_memorycd4tmem_2 <- left_join(CD161CD56cd4tmems_vs_memorycd4tmem, markers_list)

CD161CD56cd4tmems_vs_memorycd4tmem_delta <- CD161CD56cd4tmems_vs_memorycd4tmem_2 %>% dplyr::select(ENTREZID, delta_avg_exp) %>% arrange(desc(delta_avg_exp)) %>% as.data.frame() #always order from highest NES descending!
CD161CD56cd4tmems_vs_memorycd4tmem_delta    <- CD161CD56cd4tmems_vs_memorycd4tmem_delta %>% distinct(ENTREZID, .keep_all=TRUE) %>% drop_na()
CD161CD56cd4tmems_vs_memorycd4tmem_delta_FC <- CD161CD56cd4tmems_vs_memorycd4tmem_delta$delta_avg_exp
names(CD161CD56cd4tmems_vs_memorycd4tmem_delta_FC) <- CD161CD56cd4tmems_vs_memorycd4tmem_delta$ENTREZID

gse <- gsePathway(geneList=CD161CD56cd4tmems_vs_memorycd4tmem_delta_FC,
                  organism='human', minGSSize=5, maxGSSize=800,
                  pvalueCutoff=0.05, pAdjustMethod='BH', eps=0, seed=42)
gse_result <- gse@result %>% filter(setSize>15)
write_delim(gse_result, file='./GSEA/gsePathway_results_20241216.csv')

gse_to_keep <- c('Interleukin-10 signaling', 'Antimicrobial peptides', 'ER-Phagosome pathway', 'Degradation of the extracellular matrix', 'Signaling by Interleukins', 'Selenoamino acid metabolism')

x <- gse_result %>% filter(Description %in% gse_to_keep)

ggplot(data=x, aes(x=reorder(Description, NES), y=NES))+
  geom_col(aes(fill=pvalue), color='black', show.legend=TRUE)+
  scale_fill_viridis()+
  coord_flip()+
  scale_y_continuous(limits=c(-2.5,2.5))+
  labs(x='Pathway', y='NES', fill='p-value')+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=10, face='plain'))
ggsave(file='./Plots/GSEA/gsePathway_Reactome_CD161+CD56+CD4Tvs_memorycd4tmems.pdf', width=130, height=45, units='mm')


# Data export ---------------------------------------------------------

cd4tmems_metadata <- cd4tmems@meta.data
cd4tmems_metadata[is.na(cd4tmems_metadata)] <- 'n.a.'
int.TNK <- AddMetaData(int.TNK, metadata=cd4tmems_metadata)
saveRDS(int.TNK, file='./Data deposition/int.TNK.rds') # save everything in Seurat file (ca. 10 GB)

# split Seurat object
RNA_matrix <- int.TNK@assays$RNA@counts
ADT_matrix <- int.TNK@assays$ADT@counts
HTO_matrix <- GetAssayData(int.TNK, assay='HTO', slot='counts')
metadata   <- int.TNK@meta.data

# write matrix to .rds files
saveRDS(RNA_matrix, file='./Data deposition/RNA_Kammann_et_al.rds')
saveRDS(ADT_matrix, file='./Data deposition/ADT_Kammann_et_al.rds')
saveRDS(HTO_matrix, file='./Data deposition/HTO_Kammann_et_al.rds')
saveRDS(metadata,   file='./Data deposition/Metadata_Kammann_et_al.rds')

# split seurat object for each donor
seurat_list <- SplitObject(int.TNK, split.by='donor')

seurat_D1 <- seurat_list[['D1']]
seurat_D2 <- seurat_list[['D2']]
seurat_D3 <- seurat_list[['D3']]
seurat_D4 <- seurat_list[['D4']]

rm(seurat_list, int.TNK)

# data deposition
seurat_D1_RNA <- GetAssayData(seurat_D1, assay='RNA', slot='counts')
seurat_D1_ADT <- GetAssayData(seurat_D1, assay='ADT', slot='counts')
seurat_D1_HTO <- GetAssayData(seurat_D1, assay='HTO', slot='counts')

seurat_D2_RNA <- GetAssayData(seurat_D2, assay='RNA', slot='counts')
seurat_D2_ADT <- GetAssayData(seurat_D2, assay='ADT', slot='counts')
seurat_D2_HTO <- GetAssayData(seurat_D2, assay='HTO', slot='counts')

seurat_D3_RNA <- GetAssayData(seurat_D3, assay='RNA', slot='counts')
seurat_D3_ADT <- GetAssayData(seurat_D3, assay='ADT', slot='counts')
seurat_D3_HTO <- GetAssayData(seurat_D3, assay='HTO', slot='counts')

seurat_D4_RNA <- GetAssayData(seurat_D4, assay='RNA', slot='counts')
seurat_D4_ADT <- GetAssayData(seurat_D4, assay='ADT', slot='counts')
seurat_D4_HTO <- GetAssayData(seurat_D4, assay='HTO', slot='counts')

# save as rds files
saveRDS(seurat_D1_RNA, file='./Data deposition/Donor1_RNA_Kammann_et_al.rds')
saveRDS(seurat_D2_RNA, file='./Data deposition/Donor2_RNA_Kammann_et_al.rds')
saveRDS(seurat_D3_RNA, file='./Data deposition/Donor3_RNA_Kammann_et_al.rds')
saveRDS(seurat_D4_RNA, file='./Data deposition/Donor4_RNA_Kammann_et_al.rds')

saveRDS(seurat_D1_ADT, file='./Data deposition/Donor1_ADT_Kammann_et_al.rds')
saveRDS(seurat_D2_ADT, file='./Data deposition/Donor2_ADT_Kammann_et_al.rds')
saveRDS(seurat_D3_ADT, file='./Data deposition/Donor3_ADT_Kammann_et_al.rds')
saveRDS(seurat_D4_ADT, file='./Data deposition/Donor4_ADT_Kammann_et_al.rds')

saveRDS(seurat_D1_HTO, file='./Data deposition/Donor1_HTO_Kammann_et_al.rds')
saveRDS(seurat_D2_HTO, file='./Data deposition/Donor2_HTO_Kammann_et_al.rds')
saveRDS(seurat_D3_HTO, file='./Data deposition/Donor3_HTO_Kammann_et_al.rds')
saveRDS(seurat_D4_HTO, file='./Data deposition/Donor4_HTO_Kammann_et_al.rds')

rm(seurat_D1_RNA, seurat_D2_RNA, seurat_D3_RNA, seurat_D4_RNA,
   seurat_D1_ADT, seurat_D2_ADT, seurat_D3_ADT, seurat_D4_ADT,
   seurat_D1_HTO, seurat_D2_HTO, seurat_D3_HTO, seurat_D4_HTO)

