################################################################################fGSEA validation data
rm(list=ls())
gc()
library(GSEABase)
library(limma)
library(dplyr)
setwd("C:/Users/赵定康/Desktop/input")
vector=c("GSE102208","GSE114059_a","GSE114059_b",
         "GSE114059_c","GSE156082","GSE158578_a",
         "GSE158578_b","GSE188465","GSE188466",
         "GSE21815","GSE17537","GSE72970",
         "GSE188467","GSE29435","GSE29621",
         "GSE38832","GSE39904","GSE57760",
         "GSE64337_a","GSE64337_b","GSE94714",
         "GSE69687","GSE237130","UCSC.CRC_APC",
         "UCSC.CRC_TP53","UCSC.LIHC_CTNNB1","UCSC.LIHC_AXIN1",
         "UCSC.CRC_TumorNormal","GSE222023","GSE215947",
         "GSE178260","GSE135679","GSE60564",
         "GSE44076","GSE44861","GSE21510",
         "GSE68468","GSE37178","GSE18105",
         "GSE56896","GSE17385","GSE63072_a",
         "GSE63072_b","GSE61687","UCSC.CRC_CMS2",
         "GSE33113","GSE39582","GSE2109",
         "GSE13067","GSE13294","GSE14333",
         "GSE17536","GSE20916","GSE23878",
         "GSE35896","GSE37892","KFSYSCC",
         "syn26720761","syn2623706")
for (q in 1:length(vector)) {
  load(paste0(vector[q],".Rdata"))
  load(paste0(vector[q],"_G.Rdata"))
}
group_HL=read.delim("WNT_HL_V.txt",header=TRUE,sep='\t',check.names=F)
for(i in 1:length(vector)){
  group=get(paste0(vector[i],"_G"))
  group=as.data.frame(group)
  Group_HL=group_HL
  rownames(Group_HL)=Group_HL$Accession
  Group_HL=Group_HL[,-1]
  group_HL_selected=Group_HL[paste0(vector[i]),,drop=F]
  group$group=ifelse(group$group==group_HL_selected$Group1,group_HL_selected$Group1_Status,group_HL_selected$Group2_Status)
  assign(paste0(vector[i],"_G"),group)
}
limma.one.sided = function(fit, lower = FALSE){
  se.coef=sqrt(fit$s2.post) * fit$stdev.unscaled
  df.total=fit$df.prior + fit$df.residual
  pt(fit$t, df = df.total, lower.tail = lower)
}
Get_limma=function(gene=gene,group=group,alternative_limma=alternative){
  design=model.matrix(~0+factor(group$group))
  colnames(design)=levels(factor(group$group))
  rownames(design)=colnames(gene)
  compare=makeContrasts(H-L, levels=design)
  fit=lmFit(gene, design)
  fit=contrasts.fit(fit, compare)
  fit=eBayes(fit)
  if(alternative_limma=="greater"){
    result=limma.one.sided(fit, lower=FALSE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="less"){
    result=limma.one.sided(fit, lower=TRUE)
    result=as.data.frame(result)
    colnames(result)="P.Value"
    result$gene_id=rownames(result)
    topTable=topTable(fit, coef=1, number=Inf)
    topTable$gene_id=rownames(topTable)
    topTable=topTable[,c("gene_id","logFC")]
    result=merge(result,topTable,by="gene_id")
  }else if(alternative_limma=="two.sided"){
    result=topTable(fit, coef=1, number=Inf)
    result=na.omit(result)
    result$gene_id=rownames(result)
  }else{
    print("Please enter a correct method.(greater/less/two.sided)")
  }
  return(result)
}
for (i in 1:length(vector)){
  result=Get_limma(gene=as.data.frame(get(vector[i])),
                   group=as.data.frame(get(paste0(vector[i],"_G"))),
                   alternative_limma="greater")
  assign(paste0(vector[i],"_result"),result)
}
geneSets=getGmt("wpags_wpigs.gmt")
library(fgsea)
for (i in 1:length(vector)){
  alldiff=get(paste0(vector[i],"_result"))
  alldiff=alldiff[order(alldiff$logFC,decreasing = T),]
  id=alldiff$logFC
  names(id)=alldiff$gene_id
  canonical_pathways=do.call(rbind, lapply(geneSets, function(gs) {  
    data.frame(term = gs@setName, gene = gs@geneIds, stringsAsFactors = FALSE)  
  }))
  canonical_pathways=canonical_pathways[canonical_pathways$gene != "", ]
  canonical_pathways$term=as.character(canonical_pathways$term)
  canonical_pathways=canonical_pathways %>% split(x = .$gene, f = .$term)
  fgseaRes=fgsea(pathways=canonical_pathways,
                 stats=id,
                 minSize=1,
                 maxSize=10000,
                 nperm=1000)
  fgsea_result=data.frame(fgseaRes)
  fgsea_result=fgsea_result[,c("pathway","NES")]
  colnames(fgsea_result)[2]=paste0(vector[i])
  if(i==1){
    fgsea_all=fgsea_result
  }else{
    fgsea_all=merge(fgsea_all,fgsea_result,by="pathway",all=T)
  }
}
rownames(fgsea_all)=fgsea_all[,1]
fgsea_all=fgsea_all[,-1]
fgsea_all=as.data.frame(t(fgsea_all))
library(ComplexHeatmap)
library(circlize)
ha = rowAnnotation(Group = factor(c(rep("WPAGSs", 36), rep("WPIGSs", 27))),  
                   col = list(Group = c(WPAGSs = "#f15b6c", WPIGSs = "#009ad6")),  
                   show_annotation_name = FALSE)  
fgsea_all_transposed = t(as.matrix(fgsea_all))  
Heatmap(fgsea_all_transposed,
        name = "NES",   
        col = colorRamp2(c(min(fgsea_all), 0, max(fgsea_all)),   
                         c("#90d7ec", "white", "#f58f98")),  
        left_annotation = ha,
        show_row_names = F,
        show_column_names = T, 
        cluster_rows = TRUE,   
        cluster_columns = TRUE,  
        row_names_gp = gpar(fontsize = 10),  
        column_names_gp = gpar(fontsize = 10),  
        column_names_rot = 60, 
        heatmap_legend_param = list(title = "NES",   
                                    title_gp = gpar(fontsize = 12),   
                                    labels_gp = gpar(fontsize = 12),  
                                    annotation_name_gp = gpar(fontsize = 16)),  
        rect_gp = gpar(col = "#77787b", lwd = 1) 
) 