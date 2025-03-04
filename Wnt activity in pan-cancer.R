################################################################################TCGA
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/pancancer_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
WNT_Score_group=merge(pancancer_WNT_score,pancancer_group,by="row.names",all=T)
p = ggplot(WNT_Score_group, aes(Project, activity_score, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
  theme_bw() +  
  theme(  
    legend.position = "right",  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16),  
    axis.text.y = element_text(colour = "black", size = 16),
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#90d7ec", "#f58f98")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
################################################################################TARGET
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/TARGET_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
WNT_Score_group=merge(TARGET_WNT_score,TARGET_G,by="row.names",all=T)
p = ggplot(WNT_Score_group, aes(Project, activity_score, fill = Group, color = Group)) +  
  geom_boxplot(outlier.shape = 21, outlier.fill = 'black', outlier.size = 0.25, color = "black") +  
  labs(x = " ", y = "Wnt/β-catenin pathway activity score") +
  theme_bw() +  
  theme(  
    legend.position = "right",  
    panel.grid.major.x = element_line(linewidth = 1),  
    panel.grid.major.y = element_line(linewidth = 1),  
    panel.grid.minor.x = element_line(linewidth = 0.8),  
    panel.grid.minor.y = element_line(linewidth = 0.8),  
    axis.line.x = element_line(colour = 'black', linewidth = 1.5),  
    axis.line.y = element_line(colour = 'black', linewidth = 1),  
    axis.text.x = element_text(angle = 80, hjust = 1, vjust = 1, colour = "black", size = 16), 
    axis.text.y = element_text(colour = "black", size = 16), 
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(angle = 90, colour = "black", size = 16)
  ) +  
  scale_fill_manual(values = c("#90d7ec", "#f58f98")) +  
  stat_compare_means(aes(group = Group, label = ..p.signif..), method = "wilcox.test") 
p
################################################################################TCGA-TARGET
rm(list=ls())
gc()
library(stringr)
library(ggplot2)
library(ggpubr)
library(limma)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/pancancer_WNT_score.Rdata")
load("C:/Users/赵定康/Desktop/input/TARGET_WNT_score.Rdata")
pancancer_WNT_score=rbind(pancancer_WNT_score,TARGET_WNT_score)
load("C:/Users/赵定康/Desktop/input/pancancer_group.Rdata")
pancancer_group$Project=paste0("TCGA_",pancancer_group$Project)
load("C:/Users/赵定康/Desktop/input/TARGET_G.Rdata")
TARGET_G$Project=paste0("TARGET_",TARGET_G$Project)
pancancer_group=rbind(pancancer_group,TARGET_G)

WNT_Score_group=merge(pancancer_WNT_score,pancancer_group,by="row.names",all=T)
WNT_Score_group_normal=WNT_Score_group[WNT_Score_group$Group=="Normal",]
WNT_Score_group_normal=WNT_Score_group_normal[,c(1,2,5)]
WNT_Score_group_normal$Row.names=ifelse(substr(WNT_Score_group_normal$Row.names,1,4)=="TCGA",
                                        substr(WNT_Score_group_normal$Row.names, 1, 12),
                                        substr(WNT_Score_group_normal$Row.names,1,(nchar(WNT_Score_group_normal$Row.names)-4)))
WNT_Score_group_normal=as.matrix(WNT_Score_group_normal)
rownames(WNT_Score_group_normal)=WNT_Score_group_normal[,1]
WNT_Score_group_normal=WNT_Score_group_normal[,-1]
WNT_Score_group_normal=as.data.frame(avereps(WNT_Score_group_normal))
WNT_Score_group_tumor=WNT_Score_group[WNT_Score_group$Group=="Tumor",]
WNT_Score_group_tumor=WNT_Score_group_tumor[,c(1,2,5)]
WNT_Score_group_tumor$Row.names=ifelse(substr(WNT_Score_group_tumor$Row.names,1,4)=="TCGA",
                                       substr(WNT_Score_group_tumor$Row.names, 1, 12),
                                       substr(WNT_Score_group_tumor$Row.names,1,(nchar(WNT_Score_group_tumor$Row.names)-4)))
WNT_Score_group_tumor=as.matrix(WNT_Score_group_tumor)
rownames(WNT_Score_group_tumor)=WNT_Score_group_tumor[,1]
WNT_Score_group_tumor=WNT_Score_group_tumor[,-1]
WNT_Score_group_tumor=as.data.frame(avereps(WNT_Score_group_tumor))
WNT_Score_group_tumor$ID=rownames(WNT_Score_group_tumor)
WNT_Score_group_normal$ID=rownames(WNT_Score_group_normal)
rownames(WNT_Score_group_tumor)=NULL
rownames(WNT_Score_group_normal)=NULL
WNT_Score_group_tumor$type="Tumor"
WNT_Score_group_normal$type="Normal"
pairs_data=rbind(WNT_Score_group_tumor,WNT_Score_group_normal)
library(dplyr)
pairs_data <- pairs_data %>%
  group_by(ID) %>%
  filter(n() > 1)
pairs_data$activity_score=as.numeric(pairs_data$activity_score)
pairs_data=pairs_data[!pairs_data$Project%in%c("TCGA_SARC","TCGA_SKCM","TCGA_THYM","TARGET_ALL_P2","TARGET_ALL_P3"),]
mytheme <- theme(plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),  
                 axis.line = element_line(color = "black", linewidth = 0.4),  
                 panel.grid.minor = element_blank(),  
                 panel.grid.major = element_line(size = 0.1, color = "#e5e5e5"),  
                 axis.text.y = element_text(color = "black", size = 16),  
                 axis.text.x = element_text(angle = 60, hjust = 1, color = "black", size = 16),  
                 axis.title.y = element_text(color = "black", size = 16),  
                 axis.title.x = element_text(color = "black", size = 16),  
                 legend.position = "none",  
                 panel.spacing = unit(0, "lines"),  
                 plot.title = element_text(color = "black", size = 16))  
fill.color = c(Normal = "#90d7ec", Tumor = "#f58f98")  
p = ggplot(pairs_data, aes(x = type, y = activity_score, fill = type)) +  
  geom_line(aes(group = ID), position = position_dodge(0.2), color = "grey90") +  
  geom_point(aes(group = ID, size = activity_score), pch = 21,  
             size = 1.8, position = position_jitter(w = 0.1)) +  
  stat_summary(fun.data = 'mean_se', geom = "errorbar", color = "red",   
               width = 0.25, size = 0.8, position = position_dodge(0.9)) +  
  facet_wrap(. ~ Project, scales = "free_y", ncol = 10) +  
  scale_fill_manual(values = fill.color) +  
  labs(x = NULL, y = "Wnt/β-catenin pathway activity score") +  
  stat_compare_means(aes(group = type),  
                     paired = TRUE,  
                     method = "wilcox.test",  
                     label = "p.signif",  
                     size = 3) +   
  theme_bw() + mytheme +  
  ggtitle("TCGA+TARGET") 
p