setwd("D:/资料/大二第三学期/生物信息/文章复现")
library(ggplot2)
library(reshape2)

gene_anno = read.csv("gene_annotation_filtered.csv")
gene_anno = gene_anno[,-which(names(gene_anno) == "X")]

gene = c("upregulated","downregulated")
early = c(nrow(gene_anno[gene_anno['class']=='early_up',]),nrow(gene_anno[gene_anno['class']=='early_down',]))
constant = c(nrow(gene_anno[gene_anno['class']=='constant_up',]),nrow(gene_anno[gene_anno['class']=='constant_down',]))
late = c(nrow(gene_anno[gene_anno['class']=='late_up',]),nrow(gene_anno[gene_anno['class']=='late_down',]))

df = data.frame(
  gene = gene,
  Early = early,
  Constant = constant,
  Late = late
)

df2 = melt(df,id.vars='gene')
color_map = c("Early" = "#F4DC25","Constant" = "#BDD444","Late" = "#00A655")
df2$color = color_map[df2$variable]
df2$variable = factor(df2$variable,levels = c("Late","Constant","Early"))
df2$gene = factor(df2$gene,levels = c("upregulated","downregulated"))
df2 = df2[order(df2$variable),]

pdf("plot/gene_annno_number_comparison.pdf",height = 7, width = 6)
ggplot(data=df2,aes(gene,value,fill=variable,group=gene))+
    geom_bar(stat="identity",position="stack", color=NA) +
    geom_text(aes(label = value), position = position_stack(vjust = 0.85), size = 7) +
    theme_classic() +
    scale_fill_manual(values = c("Late" = "#00A655","Constant" = "#BDD444","Early" = "#F4DC25")) +
    ylab("# of genes") +
    xlab("") +
    theme(axis.text.x = element_text(size = 20,color = "black",angle = 45, hjust = 1),  
        axis.title.x = element_text(size = 25,color = "black"),
        axis.text.y = element_text(size = 20,color = "black"),  
        axis.title.y = element_text(size = 25,color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.title = element_text(size = 30, hjust = 0.5)) 
dev.off()
    
    