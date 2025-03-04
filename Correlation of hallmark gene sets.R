################################################################################TCGA hallmark ssGSEA####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
library(IOBR)
library(GSEABase)
hallmarker=getGmt("hallmark.gmt")
hallmarker_gmt=lapply(hallmarker, function(x){ x@geneIds })
names(hallmarker_gmt)=names(hallmarker)
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
load("C:/Users/赵定康/Desktop/input/pancancer_exp.Rdata")
fs=c("OV","ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC",
     "KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","LAML",
     "PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM",
     "UCEC","UCS","UVM")
tme_all_hallmarker=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=pancancer_exp[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_hallmarker=calculate_sig_score(eset = pancancer_exp_fs,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
  tme_all_hallmarker=rbind(tme_all_hallmarker,im_hallmarker)
}
save(tme_all_hallmarker,file="tme_all_hallmarker.Rdata")
################################################################################TARGET hallmark ssGSEA####
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
library(IOBR)
library(GSEABase)
hallmarker=getGmt("hallmark.gmt")
hallmarker_gmt=lapply(hallmarker, function(x){ x@geneIds })
names(hallmarker_gmt)=names(hallmarker)
load("C:/Users/赵定康/Desktop/input/TARGET_G.Rdata")
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
load("C:/Users/赵定康/Desktop/input/TARGET.Rdata")
fs=c("AML","ALL_P1","ALL_P2","ALL_P3","CCSK","NBL","OS","RT","WT")
TARGET_all_hallmarker=data.frame()
for(i in 1:length(fs)){
  group_fs=TARGET_G[TARGET_G$Project==fs[i],]
  group_fs=rownames(group_fs)
  pancancer_exp_fs=TARGET[,group_fs]
  pancancer_exp_fs=pancancer_exp_fs[rowSums(pancancer_exp_fs)>0,]
  im_hallmarker=calculate_sig_score(eset = pancancer_exp_fs,signature = hallmarker_gmt,method="ssgsea",mini_gene_count=1)
  TARGET_all_hallmarker=rbind(TARGET_all_hallmarker,im_hallmarker)
}
save(TARGET_all_hallmarker,file="TARGET_all_hallmarker.Rdata")
################################################################################hallmark相关性分析####
rm(list=ls())
gc()
load("C:/Users/赵定康/Desktop/input/tme_all_hallmarker.Rdata")
load("C:/Users/赵定康/Desktop/input/TARGET_all_hallmarker.Rdata")
tme_all=as.data.frame(tme_all_hallmarker)
tme_all=rbind(tme_all,as.data.frame(TARGET_all_hallmarker))
rownames(tme_all)=tme_all[,1]
tme_all=tme_all[,-1]
load("C:/Users/赵定康/Desktop/input/pancancer_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/TARGET_WNT_score.Rdata")
pancancer_wnt=rbind(pancancer_WNT_score,TARGET_WNT_score)
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("C:/Users/赵定康/Desktop/input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=pancancer_group[pancancer_group$Group=="Tumor",]
TARGET_G=TARGET_G[TARGET_G$Group=="Tumor",]
pancancer_group=rbind(pancancer_group,TARGET_G)
fs=unique(pancancer_group$Project)
cancer_cor=data.frame()
for(i in 1:length(fs)){
  group_fs=pancancer_group[pancancer_group$Project==fs[i],]
  tme_all_fs=tme_all[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt[rownames(group_fs),]
  pancancer_wnt_fs=pancancer_wnt_fs[,-c(2,3),drop=F]
  cycle=colnames(tme_all_fs)
  cancer_tme_all=data.frame()
  for(j in 1:length(cycle)){
    cor_test=cor.test(pancancer_wnt_fs$activity_score,tme_all_fs[,j],method="pearson")
    cancer_tme=data.frame(cancer=fs[i],tme=cycle[j],cor=cor_test$estimate,p.value=cor_test$p.value)
    cancer_tme_all=rbind(cancer_tme_all,cancer_tme)
  }
  cancer_cor=rbind(cancer_cor,cancer_tme_all)
}
cancer_cor$cor[is.na(cancer_cor$cor)]=0
cancer_cor$p.value[is.na(cancer_cor$p.value)]=1
cancer_cor$pstar=ifelse(cancer_cor$p.value > 0.05, "", 
                        ifelse(cancer_cor$p.value <= 0.05 & cancer_cor$p.value > 0.01, "*", 
                               ifelse(cancer_cor$p.value <= 0.01 & cancer_cor$p.value > 0.001, "**", "***")))
library(ggplot2)
library(dplyr)
cancer_cor$tme=factor(cancer_cor$tme, levels = colnames(tme_all_fs))
ggplot(cancer_cor, aes(cancer, tme)) +   
  geom_tile(aes(fill = cor)) +  
  scale_fill_gradient2(high = "#f58f98", mid = "white", low = "#90d7ec") +  
  geom_text(aes(label = pstar), color = "black", size = 6) +
  theme_minimal() +
  theme(  
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16)
  ) + 
  labs(fill = paste0(" * p < 0.05", "\n\n", "** p < 0.01", "\n\n", "*** p < 0.001", "\n\n", "Correlation")) +  
  theme(  
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )  