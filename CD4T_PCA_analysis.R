### Setup
required_packages <- c("readxl", "reshape2", "purrr", "ggprism", "rstatix", "scales", "stats", "RColorBrewer", "tidyverse", "FactoMineR", "factoextra")
for(package in required_packages){
  if(!require(package, character.only=TRUE)) {install.packages(package, dependencies=TRUE)
    library(package, character.only=TRUE) }}

mycolors <- c('Blood'='#E88984', 'Spleen'='#A9261F', 'Liver'='#813c5e', 'Ileum'='#8c510a', 'Caecum'='#bf812d', 'Colon'='#dfc27d', 'mLN'='#f6e8c3', 'Lung'='#c7eae5', 'LungLN'='#92c5de', 'Skin'='#01662c')
color_pattern <- c('CD161+CD56+'='#2f8ca3', 'CD161+CD56-'='#f4af2d', 'CD161-CD56+'='white','CD161-CD56-'='#bdbdbd')

script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

## pre-processing raw data

# threshold    <- 20 # remove all subsets from a tissue if lower than 20 cells
# cd4          <- read_excel('../5.2 InnateCD4_flowjo/InnateCD4subsets_flowjoexport_20240418.xlsx', sheet=1) # load unprocessed flowjo output
# cd4$Tissue      <- factor(cd4$Tissue, levels=c('Blood', 'Spleen', 'Liver', 'Ileum', 'Caecum', 'Colon', 'mLN', 'Lung', 'LungLN', 'Skin'))
# cd4$Stimulation <- case_match(cd4$Stimulation, 'none'~'unstimulated', 'PMA_I'~'PMA+Ionomycin')
# cd4$Stimulation <- factor(cd4$Stimulation, levels=c('unstimulated', 'PMA+Ionomycin'))
# cd4$Count               <- cd4$`Count_CD161+CD56+`
# cd4$`Count_CD161+CD56+` <- NULL

## isolate unstimulated samples, remove donor- and tissue-specific data if insufficient amount of cells identified
# cd4_base <- cd4 %>% dplyr::filter(InnateSubset!='CD161-CD56+') %>% dplyr::filter(Stimulation=='unstimulated') %>% group_by(Donor, Tissue) %>% dplyr::filter(!any(Count<threshold)) %>% select(-`Sample:`, -Stimulation) %>% ungroup %>% droplevels()

## isolate tissues and markers of interest, format for PCA functions
# pca_cd4_base <- cd4_base %>% dplyr::filter(Tissue %in% c('Liver', 'Ileum', 'Caecum', 'Colon')) %>%
#                 mutate(ID=paste0(Donor,Tissue,'_',InnateSubset)) %>%
#                 select(-Count, -CD161_MFI, -CD56_MFI, -CXCR5_MFI, -(25:28)) %>%
#                 column_to_rownames(var='ID') %>% select(-Donor) %>% dplyr::mutate_at(.vars=c(3:7), .funs=as.numeric)
# write_excel_csv(pca_cd4_base,file="./Tables/PCA/pca_cd4_base_data.csv")

## load processed data from supplementary file
pca_cd4_base        <- read_excel('./Supplemental_material_file_data.xlsx', sheet="Fig2a_pca_data") # read flowjo output
pca_cd4_base$Tissue <- factor(pca_cd4_base$Tissue, levels=c('Blood', 'Spleen', 'Liver', 'Ileum', 'Caecum', 'Colon', 'mLN', 'Lung', 'LungLN', 'Skin'))

# PCA at resting conditions for hepatic and intestinal CD4 T cells identified by CD161 and CD56 co-expression pattern
pca_cd4_base <- FactoMineR::PCA(pca_cd4_base, graph=FALSE, scale.unit=TRUE, quali.sup=c('Tissue', 'InnateSubset'))
write_excel_csv(as.data.frame(pca_cd4_base[['var']][['contrib']]), file=paste0('./Tables_PCA_Liver_and_Gut_contribution.csv'))
write_excel_csv(as.data.frame(pca_cd4_base[['eig']]), file=paste0('./Tables_PCA_Liver_and_Gut_Variance.csv'))

fviz_screeplot(pca_cd4_base, choice='variance', addlabels=TRUE, ylim=c(0, 50), ggtheme=theme_tk, barcolor='black', barfill='darkgrey')+labs(title='', x='PC', y='Explained variance (%)')+scale_y_continuous(expand=c(0,0))

### Figure 2 PCA
# PCA plot of CD4+ T cells identified by CD161 and CD56 co-expression pattern, colored by tissue origin (Fig 2a)
x <- fviz_pca_ind(pca_cd4_base, geom.ind='point', geom.var=c('arrow', 'text'), title='', col.var='contrib', gradient.cols=c("#d9d9d9", "#969696", "#000000"), select.var=list(contrib=10))
y <- x + geom_point(aes(fill=pca_cd4_base$call$X$Tissue), color='black', shape=21, size=2, show.legend=FALSE)+
  scale_fill_manual(values=mycolors)+
  scale_x_continuous(limits=c(-8,8), breaks=c(-8,0,8), expand=c(0,0))+
  scale_y_continuous(limits=c(-4,8), breaks=c(-4,0,8), expand=c(0,0))+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_base_liver_intestine_by_tissue.pdf', height=50, width=50, unit='mm')

# Variables plots (Fig 2b)
x <- fviz_pca_biplot(pca_cd4_base, geom.ind='', geom.var=c('arrow', 'text'), title='', col.var='black')
y <- x +
  scale_x_continuous(limits=c(-8,8), breaks=c(-8,0,8), expand=c(0,0))+
  scale_y_continuous(limits=c(-4,8), breaks=c(-4,0,8), expand=c(0,0))+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_base_liver_intestine_contribution.pdf', height=50, width=50, unit='mm')

# PCA plot, colored by co-expression pattern (Fig 2c)
x <- fviz_pca_ind(pca_cd4_base, geom.ind='point', geom.var=c('arrow', 'text'), title='', col.var='contrib', gradient.cols=c("#d9d9d9", "#969696", "#000000"), select.var=list(contrib=10))
y <- x + geom_point(aes(fill=pca_cd4_base$call$X$InnateSubset), color='black', shape=21, size=2, show.legend=FALSE)+
  labs(fill='Celltype')+
  scale_fill_manual(values=color_pattern)+
  scale_x_continuous(limits=c(-8,8), breaks=c(-8,0,8), expand=c(0,0))+
  scale_y_continuous(limits=c(-4,8), breaks=c(-4,0,8), expand=c(0,0))+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_base_liver_intestine_by_coexpressionpattern.pdf', height=50, width=50, unit='mm')

### Figure 3 PCA
## isolate PMA+Ionomycin-stimulated samples, remove donor- and tissue-specific data if insufficient amount of cells identified
# cd4_stim     <- cd4 %>% dplyr::filter(InnateSubset!='CD161-CD56+') %>% dplyr::filter(Stimulation=='PMA+Ionomycin') %>% group_by(Donor, Tissue)%>% dplyr::filter(!any(Count<threshold)) %>% select(-`Sample:`, -Stimulation) %>% ungroup %>% droplevels()
# pca_cd4_stim <- cd4_stim %>% filter(Tissue %in% c('Ileum', 'Caecum', 'Colon')) %>%
#                 mutate(ID=paste0(Donor,Tissue,'_',InnateSubset)) %>%
#                 select(-CXCR5_MFI, -CD161_MFI, -CD56_MFI, -CXCR3_MFI,
#                 -CD27, -CD39, -CD62L, -CD69, -CD103, -CD127_MFI, -CXCR5) %>%
#                 column_to_rownames(var='ID') %>% select(-Donor) %>% dplyr::mutate_at(.vars=2, .funs=as.factor)
# write_excel_csv(pca_cd4_stim, file="./pca_cd4_stim_raw_data.csv")

## load processed data from supplementary file
pca_cd4_stim        <- read_excel('./Supplemental_material_file_data.xlsx', sheet="Fig3_pca_data") # read flowjo output
pca_cd4_stim$Tissue <- factor(pca_cd4_stim$Tissue, levels=c('Blood', 'Spleen', 'Liver', 'Ileum', 'Caecum', 'Colon', 'mLN', 'Lung', 'LungLN', 'Skin'))
pca_cd4_stim <- pca_cd4_stim %>% select(-Donor, -ICC, `HLA-DR`, -`PD-1`, -`IFNg_IL-22`, `IFNg_TNF`, `IL17_IL-22`)

pca_cd4_stim <- FactoMineR::PCA(pca_cd4_stim, graph=FALSE, scale.unit=TRUE, quali.sup=c('Tissue', 'InnateSubset'))
write_excel_csv(as.data.frame(pca_cd4_stim[['var']][['contrib']]), file=paste0('./PCA_stim_gut_contribution.csv'))
write_excel_csv(as.data.frame(pca_cd4_stim[['eig']]), file=paste0('./PCA_stim_liver_variance.csv'))

fviz_screeplot(pca_cd4_stim, choice='variance', addlabels=TRUE, ylim=c(0, 50), ggtheme=theme_tk, barcolor='black', barfill='darkgrey')+labs(title='', x='PC', y='Explained variance (%)')+scale_y_continuous(expand=c(0,0))

## Figure 3 PCA
# PCA plot of PMA+Ionomycin-stimulated CD4+ T cells identified by CD161 and CD56 co-expression pattern, colored by tissue origin (Fig 3a)
x <- fviz_pca_ind(pca_cd4_stim, geom.ind='point', geom.var=c('arrow', 'text'), title='', col.var='black')
y <- x + geom_point(aes(fill=pca_cd4_stim$call$X$Tissue), shape=21, size=2, show.legend=FALSE)+
  scale_fill_manual(values=mycolors)+
  scale_x_continuous(limits=c(-3,6.5), breaks=c(-3,0,6), expand=c(0,0), guide='prism_offset')+
  scale_y_continuous(limits=c(-6,6), breaks=c(-6,0,6), expand=c(0,0), guide='prism_offset')+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_stim_intestines_by_tissue.pdf', height=50, width=50, unit='mm')

x <- fviz_pca_biplot(pca_cd4_stim, geom.ind='', geom.var=c('arrow', 'text'), title='', col.var='black')
y <- x +
  scale_x_continuous(limits=c(-3,6.5), breaks=c(-3,0,6), expand=c(0,0), guide='prism_offset')+
  scale_y_continuous(limits=c(-6,6), breaks=c(-6,0,6), expand=c(0,0), guide='prism_offset')+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_stim_intestines_contribution.pdf', height=50, width=50, unit='mm')

# PCA plot, colored by co-expression pattern of CD161 and CD56
x <- fviz_pca_ind(pca_cd4_stim, geom.ind='point', geom.var=c('arrow', 'text'), title='', col.var='black')
y <- x + geom_point(aes(fill=pca_cd4_stim$call$X$InnateSubset), shape=21, size=2, show.legend=FALSE)+
  scale_fill_manual(values=color_pattern)+
  scale_x_continuous(limits=c(-3,6.5), breaks=c(-3,0,6), expand=c(0,0), guide='prism_offset')+
  scale_y_continuous(limits=c(-6,6), breaks=c(-6,0,6), expand=c(0,0), guide='prism_offset')+
  theme_tk+
  theme(panel.grid.major.y=element_blank())
print(y)
ggsave(filename='./Plot_PCA_stim_intestines_by_coexpressionpattern.pdf', height=50, width=50, unit='mm')
