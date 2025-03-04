################################################################################PCA training data quality testing
rm(list=ls())
gc()
list=c("GSE232944_a","GSE232944_b","GSE199835_a","GSE199835_b","GSE199835_c",
       "GSE120106","GSE234403","GSE39904_a","GSE39904_b","GSE128281",
       "GSE62060","GSE33143","GSE18560_a","GSE18560_b","GSE44097_a","GSE44097_b",
       "GSE73906_a","GSE73906_b","GSE73906_c","GSE73906_d","GSE73906_e",
       "GSE73906_f","GSE123807","GSE89124")
for (q in 1:length(list)) {
  load(paste0("C:/Users/赵定康/Desktop/input/",list[q],".Rdata"))
  load(paste0("C:/Users/赵定康/Desktop/input/",list[q],"_G.Rdata"))
}
###############################PCA
library(ggplot2)
library(plyr)
for (i in 1:length(list)) {
  exp=get(list[i])
  pca_data=as.data.frame(exp)
  group=get(paste0(list[i],"_G"))
  Group=as.data.frame(group)
  pca_data=scale(pca_data)
  pca_data=t(scale(t(pca_data),center=TRUE,scale=F))
  pca=prcomp(t(pca_data),center=FALSE,scale.=FALSE)
  PCA_sum=summary(pca)
  PC12=pca$x[,1:2]
  pc=PCA_sum$importance[2,]*100
  PC12=as.data.frame(PC12)
  PC12$samples=row.names(PC12)
  colnames(Group)=c("samples","group")
  df=merge(PC12,Group,by="samples")
  color=c("#f58f98","#90d7ec")
  p1=ggplot(data = df, aes(x = PC1, y = PC2, color = group, shape = group)) +  
    ggtitle(paste0(list[i])) +  
    theme_bw() +  
    geom_point(size = 3) +  
    theme(panel.grid = element_blank()) +  
    geom_vline(xintercept = 0, linetype = 'dotdash', size = 1, color = 'grey') +  
    geom_hline(yintercept = 0, linetype = 'dotdash', size = 1, color = 'grey') +  
    labs(x = paste0("PC1 (", pc[1], "%)"), y = paste0("PC2 (", pc[2], "%)")) +  
    scale_color_manual(values = color) +  
    theme(axis.title.x = element_text(size = 16,color = "black"),  
          axis.title.y = element_text(size = 16, angle = 90,color = "black"),  
          axis.text.y = element_text(size = 16,color = "black"),  
          axis.text.x = element_text(size = 16,color = "black"),  
          plot.title = element_text(size = 20,color = "black"),  
          legend.title = element_text(size = 14,color = "black"),  
          legend.text = element_text(size = 12,color = "black")) 
  library(plyr)
  cluster_border <- ddply(df, 'group', function(df) df[chull(df[[2]], df[[3]]), ])
  p=p1 + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) +
    scale_fill_manual(values = c("#f58f98","#90d7ec"))
  p
  setwd("C:/Users/赵定康/Desktop/图片20250217/sfig1")
  ggsave(paste0(list[i],"_PCA",".pdf"), p, width = 5, height = 3)
}