# Unsupervised analysis (UMAP, CytoTree) ---------------------------------------------------

# Setup
if (!require('devtools', quietly=TRUE)) install.packages('devtools')
if (!require('flowCore', quietly=TRUE)) devtools::install:github("RGLab/flowCore") # Hahne et al. (2009), DOI: 10.1186/1471-2105-10-106
if (!require('flowUtils', quietly=TRUE)) devtools::install_github("jspidlen/flowUtils") #Spidlen et al.(2021)
if (!require('CytoTree', quietly=TRUE)) devtools::install_github("JhuangLab/CytoTree") # Dai et al. (2020), DOI: 10.1186/s12859-021-04054-2

required_packages <- c("readxl", "reshape2", "see", "purrr", "ggprism", "rstatix", "scales", "stats", "pheatmap", "RColorBrewer", "Hmisc", "corrplot", "vtable", "ggpubr", "tidyverse")
for(package in required_packages){
  if(!require(package, character.only=TRUE)) {install.packages(package, dependencies=TRUE)
    library(package, character.only=TRUE) }}
umap_libraries <- c("flowUtils", "flowCore", "CytoTree")
lapply(umap_libraries, library, character.only=TRUE)
rm(required_packages)
mycolors <- c('Blood'='#E88984', 'Spleen'='#A9261F', 'Liver'='#813c5e', 'Ileum'='#8c510a', 'Caecum'='#bf812d', 'Colon'='#dfc27d', 'mLN'='#f6e8c3', 'Lung'='#c7eae5', 'LungLN'='#92c5de', 'Skin'='#01662c')

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

mycolors  <- c('Blood'='#E88984',
               'Spleen'='#A9261F',
               'Liver'='#813c5e',
               'Ileum'='#8c510a',
               'Caecum'='#bf812d',
               'Colon'='#dfc27d',
               'mLN'='#f6e8c3',
               'Lung'='#c7eae5',
               'LungLN'='#92c5de',
               'Skin'='#01662c')

# IHOPE analysis
fcs_path  <- '/Users/tobias.kammann/Projects/5. InnateCD4/5.2 InnateCD4_flowjo/5.2 downsampled concat'
fcs_file  <- list.files(fcs_path, pattern='.fcs', full=TRUE)
fcs_data <- runExprsExtract(fcs_file, comp=FALSE, transformMethod='none')
fcs_protein_data <- fcs_data[, 7:33] # get protein values
fcs_meta_data <- fcs_data[,c(35,36,37,39,40)] # cut flowjo's numerical distribution into discrete metadata
fcs_meta_data[,"Donor_pseudocode<NA>"]  <- as.numeric(as.character(cut(fcs_meta_data[,"Donor_pseudocode<NA>"],  breaks=c(5500,
                                                                                                                         6500, #1
                                                                                                                         7500, #2
                                                                                                                         9500, #3
                                                                                                                         10500,#4
                                                                                                                         11500,#5
                                                                                                                         12500,#6
                                                                                                                         13500,#7
                                                                                                                         14500,#8
                                                                                                                         20500,#11
                                                                                                                         22500,#12
                                                                                                                         24500,#13
                                                                                                                         25500), labels=c('6','7','9','10','11','12','13','14','15','18','20','22','24','25'), right=FALSE)))
fcs_meta_data[,"Tissue_pseudocode<NA>"] <- as.numeric(as.character(cut(fcs_meta_data[,"Tissue_pseudocode<NA>"], breaks=c(500,1500,#Blood
                                                                                                                         2500,#Spleen
                                                                                                                         3500,#Liver
                                                                                                                         4500,#Ileum
                                                                                                                         5500,#Caecum
                                                                                                                         6500,#Colon
                                                                                                                         7500,#mLN
                                                                                                                         8500,#Lung
                                                                                                                         9500,#LungLN
                                                                                                                         10500), labels=c(1,2,3,4,5,6,7,8,9,10))))
fcs_meta_data <- as.data.frame(fcs_meta_data)
colnames(fcs_meta_data) <- c('Donor', 'Donor_pseudocode', 'Stimulation', 'Tissue', 'Tissue_pseudocode')
fcs_protein_data_columns <- c("FJComp-APC-A<NA>"='Va7.2', "FJComp-APC-H7-A<NA>"='IL-22', "FJComp-Alexa Fluor 700-A<NA>"='CD69', "FJComp-BB515-A<NA>"='PD-1', "FJComp-BB630-A<NA>"='Gnly', "FJComp-BB700-A<NA>"="CD4", "FJComp-BB755-A<NA>" ="Perforin", "FJComp-BB790-A<NA>"="GzmB",
                              "FJComp-BUV395-A<NA>"="CD103", "FJComp-BUV496-A<NA>"="CD39", "FJComp-BUV563-A<NA>"="HLA-DR", "FJComp-BUV615-A<NA>"="CD27", "FJComp-BUV661-A<NA>"="CD127",  "FJComp-BUV737-A<NA>"="CD56", "FJComp-BUV805-A<NA>"="CD45", "FJComp-BV421-A<NA>"="CXCR5",
                              "FJComp-BV510-A<NA>"="Dump",  "FJComp-BV570-A<NA>"="CD8",  "FJComp-BV605-A<NA>"="IL17A", "FJComp-BV650-A<NA>"="CD3", "FJComp-BV711-A<NA>"="TNF", "FJComp-BV750-A<NA>"="CD62L", "FJComp-BV786-A<NA>"="IFNg", "FJComp-PE-A<NA>"="5OPRU-Tet",
                              "FJComp-PE-CF594-A<NA>"="IL-10", "FJComp-PE-Cy5-A<NA>"="CD161", "FJComp-PE-Cy7-A<NA>"="CXCR3")
colnames(fcs_protein_data)[match(names(fcs_protein_data_columns), colnames(fcs_protein_data))] <- fcs_protein_data_columns

# create CYT object (S4 object)
i4 <- createCYT(raw.data=fcs_protein_data, normalization.method='log')
i4@meta.data$Donor_pseudocode  <- fcs_meta_data$Donor_pseudocode
i4@meta.data$Tissue_pseudocode <- fcs_meta_data$Tissue_pseudocode
i4@meta.data$Tissue       <- case_match(i4@meta.data$Tissue_pseudocode, 1~'Blood', 2~'Spleen', 3~'Liver', 4~'Ileum', 5~'Caecum', 6~'Colon', 7~'mLN', 8~'Lung', 9~'LungLN', 10~'Skin')
i4@meta.data$Tissue       <- factor(i4@meta.data$Tissue, levels=c('Blood', 'Spleen', 'Liver', 'Ileum', 'Caecum', 'Colon', 'mLN', 'Lung', 'LungLN', 'Skin'))
i4@meta.data$Donor        <- factor(i4@meta.data$Donor_pseudocode, levels=c('IHOPE06', 'IHOPE07', 'IHOPE09', 'IHOPE10', 'IHOPE11', 'IHOPE12', 'IHOPE13', 'IHOPE14', 'IHOPE15', 'IHOPE18', 'IHOPE20', 'IHOPE22', 'IHOPE24', 'IHOPE25'))

tissues <-  levels(i4@meta.data$Tissue)
donors  <- levels(i4@meta.data$Donor)
rm(fcs_data, fcs_meta_data, fcs_protein_data)

clustering_markers <- c('CD161','CD27', 'CD39', 'CD56', 'CD62L' ,'CD69', 'CD103', 'CD127', 'CXCR3', 'CXCR5', 'PD-1', 'HLA-DR', 'GzmB', 'Perforin', 'Gnly', 'IL-10', 'IL17A', 'IL-22', 'TNF', 'IFNg')
i4 <- changeMarker(i4, markers=clustering_markers)

i4cyto <- lapply(tissues, function(i){
  cells  <- i4@meta.data[i4@meta.data$Tissue==i, ]
  donors <- unique(cells$Donor_pseudocode)
  sampled_cells_df <- data.frame()
  # Initial sampling per tissue
  for (j in donors) {
    donor_cells     <- cells[cells$Donor_pseudocode==j, ] # cells per donor
    cells_to_sample <- min(nrow(donor_cells), 3000)
    if (cells_to_sample>0) {
      set.seed(42)
      sampled_cells_df <- rbind(sampled_cells_df, donor_cells[sample(nrow(donor_cells), cells_to_sample), ])
    }
  }
  # Downsampling to 10000 cells per tissue if too many
  if (nrow(sampled_cells_df)>10000) {
    set.seed(42)
    sampled_cells_df <- sampled_cells_df[sample(nrow(sampled_cells_df), 10000), ]
  }
  sampled_cells_df
})

i4cyto <- do.call(rbind, i4cyto)
table(i4cyto$Tissue, i4cyto$Donor_pseudocode)

sampled_cells <- rownames(i4cyto) # vector of cells to keep

i4@meta.data$Blood  <- case_match(i4@meta.data$Tissue, 'Blood'~'Blood', tissues[-1]~'other')
i4@meta.data$Spleen <- case_match(i4@meta.data$Tissue, 'Spleen'~'Spleen', tissues[-2]~'other')
i4@meta.data$Liver  <- case_match(i4@meta.data$Tissue, 'Liver'~'Liver', tissues[-3]~'other')
i4@meta.data$Ileum  <- case_match(i4@meta.data$Tissue, 'Ileum'~'Ileum', tissues[-4]~'other')
i4@meta.data$Caecum <- case_match(i4@meta.data$Tissue, 'Caecum'~'Caecum', tissues[-5]~'other')
i4@meta.data$Colon  <- case_match(i4@meta.data$Tissue, 'Colon'~'Colon', tissues[-6]~'other')
i4@meta.data$mLN    <- case_match(i4@meta.data$Tissue, 'mLN'~'mLN', tissues[-7]~'other')
i4@meta.data$Lung   <- case_match(i4@meta.data$Tissue, 'Lung'~'Lung', tissues[-8]~'other')
i4@meta.data$LungLN <- case_match(i4@meta.data$Tissue, 'LungLN'~'LungLN', tissues[-9]~'other')
i4@meta.data$Skin   <- case_match(i4@meta.data$Tissue, 'Skin'~'Skin', tissues[-10]~'other')

i4@meta.data$Blood  <- factor(i4@meta.data$Blood)
i4@meta.data$Spleen <- factor(i4@meta.data$Spleen)
i4@meta.data$Liver  <- factor(i4@meta.data$Liver)
i4@meta.data$Ileum  <- factor(i4@meta.data$Ileum)
i4@meta.data$Caecum <- factor(i4@meta.data$Caecum)
i4@meta.data$Colon  <- factor(i4@meta.data$Colon)
i4@meta.data$mLN    <- factor(i4@meta.data$mLN)
i4@meta.data$Lung   <- factor(i4@meta.data$Lung)
i4@meta.data$LungLN <- factor(i4@meta.data$LungLN)
i4@meta.data$Skin   <- factor(i4@meta.data$Skin)

i4@meta.data$allTissues <- i4@meta.data$Tissue
i4@meta.data$allTissues <- factor(case_match(i4@meta.data$allTissues, tissues~'all'))

i4sub <- subsetCYT(i4, cells=sampled_cells)
saveRDS(i4sub, './i4_all_tissues_sampled.rds')
saveRDS(i4, './i4_all_tissues.rds')

i4_kmeans <- i4sub
i4_kmeans <- runCluster(i4_kmeans, cluster.method='kmeans', verbose=TRUE, iter.max=100) # computes kmeans
i4_kmeans <- processingCluster(i4_kmeans, k=20, perplexity=5, downsampling.size=1, verbose=TRUE) # no downsampling
i4_kmeans <- runUMAP(i4_kmeans, n_neighbors=20, verbose=TRUE)

plot2D(i4_kmeans, item.use=c('UMAP_1', 'UMAP_2'), alpha=1, category='categorical', color.by='Donor_pseudocode', size=1)+theme_transparent()
ggsave('./kmeans/Plots/UMAP_all_tissues_donor_combined_kmeans.pdf', width=100, height=80, units='mm')

plot2D(i4_kmeans, item.use=c('UMAP_1', 'UMAP_2'), color.by='CD161', alpha=0.9, category='numeric', size=1)+theme_transparent()+labs(color='')+
       scale_color_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))+theme(legend.position='none')
ggsave('./kmeans/Plots/UMAP_all_tissues_CD161.pdf', width=100, height=80, units='mm')

plot2D(i4_kmeans, item.use=c('UMAP_1', 'UMAP_2'), color.by='CD56', alpha=0.9, category='numeric', size=1)+theme_transparent()+labs(color='')+
       scale_color_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))+theme(legend.position='none')
ggsave('./kmeans/Plots/UMAP_all_tissues_CD56.pdf', width=80, height=80, units='mm')

i4_kmeans <- buildTree(i4_kmeans, dim.type='umap', dim.use=1:2)
saveRDS(i4_kmeans, file='./i4_kmeans.rds')

# Tree for each marker
for (i in clustering_markers){
  plotTree(i4_kmeans, color.by=i, show.node.name=FALSE, cex.size=0.8, as.tree=FALSE)+
    scale_colour_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))+
    theme(legend.position='none', title=element_blank())
  ggsave(paste0('./kmeans/Plots/Treeplot_on_umap_',i,'.pdf'), width=60, height=60, units='mm')
}

plotTree(i4_kmeans, color.by='CD56', show.node.name=FALSE, cex.size=0.8, as.tree=FALSE)+
  scale_colour_gradientn(colors=c("#4575b4", "#fafac8", "#b2182b"))+
  theme(title=element_blank())
ggsave(paste0('./kmeans/Plots/Treeplot_on_umap_legend.pdf'), width=60, height=60, units='mm')

# Tree for Tissue
plotTree(i4_kmeans, color.by='branch.id', show.node.name=FALSE, cex.size=1)+theme(title=element_blank())+
  scale_color_gradientn(colors=c('#ffb703','lightgrey','#7678ed','#2a9d8f' ))
ggsave(paste0('./kmeans/Plots/Treeplot_allTissues_UMAP.pdf'), width=80, height=80, units='mm')

plotTree(i4_kmeans, color.by='branch.id', show.node.name=TRUE, cex.size=1)+theme(title=element_blank())+
  scale_color_gradientn(colors=c('#ffb703','lightgrey','#7678ed','#2a9d8f' ))
ggsave(paste0('./kmeans/Plots/Treeplot_allTissues_UMAP_with_clusterID.pdf'), width=180, height=180, units='mm')

# Heatmap of tree branches
diff.info <- runDiff(i4_kmeans)
plot1 <- plotBranchHeatmap(i4_kmeans, colorRampPalette(c("#4575b4", "#fafac8", "#b2182b"))(100), scale='row', cellwidth=8, cellheight=8, fontsize=8, linewidth=0.1)
ggsave(plot1, file='./kmeans/Plots/Heatmap_allTissues_treebranches_umap.pdf', width=120, height=100, units='mm')

plot2 <- plotClusterHeatmap(i4_kmeans, colorRampPalette(c("#4575b4", "#fafac8", "#b2182b"))(100),scale='row')
ggsave(plot2, file='./kmeans/Plots/Heatmap_allTissues_treeclusters_umap.pdf', width=120, height=100, units='mm')
rm(plot, plot2)

#plotPieTree
oldstage <- fetchPlotMeta(i4_kmeans)
i4_kmeans@meta.data$stage <- i4_kmeans@meta.data$Tissue
plotPieTree(i4_kmeans, size.by.cell.number=TRUE, cex.size=4)+theme_transparent()+labs(fill='Tissue')+
  scale_fill_manual(values=c('#E88984', '#A9261F', '#813c5e', '#8c510a', '#bf812d', '#dfc27d', '#f6e8c3', '#c7eae5', '#92c5de', '#01662c'))
ggsave('./kmeans/Plots/Treeplot_tissue_contribution.pdf', width=130, height=80, units='mm')

