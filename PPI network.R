################################################################################PPI network
rm(list=ls())
gc()
setwd("C:/Users/赵定康/Desktop/input")
library(clusterProfiler)
activation_inhibition_geneset=read.gmt("wpags_wpigs.gmt")
colnames(activation_inhibition_geneset)=c("group","gene")
activation_inhibition_geneset=activation_inhibition_geneset[,c("gene","group")]
activation_inhibition_geneset$group = substr(activation_inhibition_geneset$group, 1, 5)  
activation_inhibition_geneset = unique(activation_inhibition_geneset)
library(tidyverse)
library(STRINGdb)
string_db=STRINGdb$new( version="12.0",
                        species=9606,
                        score_threshold=600,
                        input_directory="C:/Users/赵定康/Desktop/input")
data_frame=data.frame(gene = activation_inhibition_geneset$gene)
dat_map=string_db$map(my_data_frame=data_frame, 
                      my_data_frame_id_col_names="gene",
                      removeUnmappedRows=TRUE)
hits=dat_map$STRING_id 
dat_link=string_db$get_interactions(hits)
dat_link$from=dat_map[match(dat_link$from,dat_map$STRING_id),'gene']
dat_link$to=dat_map[match(dat_link$to,dat_map$STRING_id),'gene'] 
colnames(dat_link)=c('node1','node2','combined_score')
dat_link=dat_link %>% distinct(node1, node2, .keep_all = T)
dat_link=dat_link[dat_link$node1%in%activation_inhibition_geneset$gene,]
dat_link=dat_link[dat_link$node2%in%activation_inhibition_geneset$gene,]
write.csv(dat_link,'string_link.csv',row.names = F,quote = F)
write.csv(activation_inhibition_geneset,'group_PPI.csv',row.names = F,quote = F)