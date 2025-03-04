################################################################################Mining of significant genes
rm(list=ls())
gc()
library(dplyr)
setwd("C:/Users/赵定康/Desktop/input")
load("C:/Users/赵定康/Desktop/input/up.p.matrix.Rdata")
load("C:/Users/赵定康/Desktop/input/dn.p.matrix.Rdata")
up.p.matrix=up.p.matrix[1:150,length(colnames(up.p.matrix)),drop=F]
up.p.matrix$gene=rownames(up.p.matrix)
dn.p.matrix=dn.p.matrix[1:150,length(colnames(dn.p.matrix)),drop=F]
dn.p.matrix$gene=rownames(dn.p.matrix)
library(ggplot2)
ggplot(up.p.matrix, aes(x = reorder(gene, -p_score), y = p_score, fill = p_score)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) + 
  labs(
    title = "",
    subtitle = "",
    x = "Activation genes (top 150)",  
    y = expression(-log[10](italic(p)[UP]) * "-combined")
  ) +
  theme_minimal(base_size = 16) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"), 
    axis.text = element_text(size = 16, color = "black"),  
    axis.title = element_text(size = 30, color = "black"),  
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 0, hjust = 1, color = "black")
  ) +  
  scale_fill_gradient(low = "#fffffb", high = "#f58f98") +
  coord_flip()
ggplot(dn.p.matrix, aes(x = reorder(gene, -p_score), y = p_score, fill = p_score)) +  
  geom_bar(stat = "identity", color = "black", width = 0.7) + 
  labs(
    title = "",
    subtitle = "",
    x = "Inhibition genes (top 150)",  
    y = expression(-log[10](italic(p)[DN]) * "-combined")
  ) +
  theme_minimal(base_size = 16) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, color = "black"), 
    axis.text = element_text(size = 16, color = "black"),  
    axis.title = element_text(size = 30, color = "black"),  
    panel.grid.major = element_line(size = 0.5, linetype = 'dotted', colour = "grey"),  
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(face = "italic", size = 16, angle = 0, hjust = 1, color = "black")
  ) +  
  scale_fill_gradient(low = "#fffffb", high = "#90d7ec") +
  coord_flip()