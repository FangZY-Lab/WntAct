################################################################################molecular characterisation
rm(list=ls())
gc()
library(GSVA)
library(GSEABase)
library(reshape)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/cluster.Rdata")
load("C:/Users/赵定康/Desktop/input/crc_data_all.Rdata")
GeneSets=getGmt("NM_genesets.gmt")
geneset=gsva(expr=as.matrix(crc_data_all), GeneSets, kcdf="Gaussian",method = "ssgsea",min.sz=1,max.sz=10000)
geneset=as.data.frame(t(geneset))
geneset=as.data.frame(scale(geneset))
geneset=geneset[,-c(4,5)]
name=as.data.frame(colnames(geneset))
############################Infiltration ESTIMATE
geneset_zdk=geneset[,c("IMMUNE_ESTIMATE","STROMAL_ESTIMATE")]
colnames(geneset_zdk)=c("Immune infiltration","Stromal infiltration")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"), 
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Infiltration ESTIMATE")
p 
############################Signatures
geneset_zdk=geneset[,c("SERRATED_UP","EMT_CORE_GENES","MYC_TARGETS_ZELLER","TGFB_CORE_GENES","CSC_BATLLE",
                       "MATRIX_REMODEL_REACTOME","WOUND_RESPONSE_GO_BP","CRYPT_BASE","CRYPT_TOP","EPITH_LOBODA","MESENCH_LOBODA")]
colnames(geneset_zdk)=c("Serrated lesion","EMT activation","MYC targets","TGF-β activation","Cancer stem cell",
                        "Matrix remodeling","Wound response","Crypt base","Crypt top","Epithelial","Mesenchymal")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Signatures") 
p 
############################Pathways
geneset_zdk=geneset[,c("FGF_ACT_CP","NOTCH_REACTOME","NFKB_BIOCARTA",
                       "MAPK_KEGG","PI3K_ACT_REACTOME","SRC_ACT_BIOCARTA","JAK_STAT_KEGG","CASPASE_BIOCARTA",
                       "PROTEASOME_KEGG","CELL_CYCLE_BIOCARTA","SHH_KEGG",
                       "INTEGRIN_BETA3_CP","VEGF_VEGFR_REACTOME","TRANSLATION_RIBOS_REACTOME")]
colnames(geneset_zdk)=c("FGF","NOTCH","NF-κB","MAPK","PI3K","SRC","JAK-STAT","Caspases","Proteosome",
                        "Cell cycle","SHH","Integrin-β3","VEGF,VEGFR","Translation ribosome")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Pathways")   
p 
############################Immune
geneset_zdk=geneset[,c("IMMUNE_RESP_GO_BP","IMMUNE_MDSC_ALBELDA",
                       "IMMUNE_NKC_BREAST","IMMUNE_TH1_GALON","IMMUNE_TH17_GOUNARI",
                       "IMMUNE_THF_BREAST","IMMUNE_TREG_GALON","PD1_REACTOME",
                       "COMPLEMENT_COAG_KEGG")]
colnames(geneset_zdk)=c("Immune response","MDSC",
                        "NK cell infiltration","T(H)1 infiltration","T(H)17 activation",
                        "TF(H) infiltration","T(reg) activation","PD1 activation",
                        "Complement activation")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"), 
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Immune")  
p
############################Metabolism
geneset_zdk=geneset[,c("ARACHNOID_METAB_KEGG","FATTY_ACID_METAB_KEGG","FRUTOSE_MANNOSE_METAB_KEGG",
                       "GALACTOSE_METAB_KEGG","GLUTATHIONE_KEGG","GLYCEROPHOSPHOLIPID_METAB_KEGG",
                       "LINOLEIC_METAB_KEGG","LYSOPHOSPHOLIPID_PID","NUCLEOTIDE_METAB_REACTOME",
                       "PENTOSE_GLUC_METAB_KEGG","STARCH_SUCROSE_METAB_KEGG","TYROSINE_METAB_KEGG",
                       "AMINO_SUGAR_NUCLEO_METAB_KEGG","ALANINE_ASPARTATE_GLUTAMATE_KEGG")]
colnames(geneset_zdk)=c("Arachnoid"," Fatty acid","Frutose Mannose",
                        "Galactose","Glutamine","Glycerophospholipid",
                        "Linoleic","Lysophospholipid","Nucleotide",
                        "Pentose","Starch sucrose","Tyrosine",
                        "Sugar,aa,nucleotide","Alainine,aspartate,glutamate")
geneset_zdk_cluster=merge(cluster,geneset_zdk,by.x="Tag",by.y="row.names")
geneset_zdk_cluster=melt(geneset_zdk_cluster)
colnames(geneset_zdk_cluster)=c("Tag","Group","geneset","value")
geneset_zdk_cluster$Group=factor(geneset_zdk_cluster$Group,levels=c("WSS1","WSS2","WSS3","WSS4"))
library(ggplot2)
library(ggsignif)
library(ggpubr)
p = ggplot(geneset_zdk_cluster, aes(geneset, value, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "ssGSEA score") +  
  theme_minimal() +
  theme(  
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),  
    axis.line = element_line(color = 'black', size = 0.8),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"), 
    axis.title.y = element_text(angle = 90, size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black")
  ) +  
  scale_fill_manual(values = c("#d71345", "#45b97c", "#145b7d", "#ffd400")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..)) +  
  ggtitle("Metabolism") 
p