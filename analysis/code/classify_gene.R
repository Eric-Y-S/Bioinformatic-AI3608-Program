setwd("D:/资料/大二第三学期/生物信息/文章复现")
library(DESeq2)
library(dplyr)

# read data
count = read.csv("samples_NONE.csv",row.names = 1)
rpm = read.csv("samples_RPM.csv",row.names = 1)
# set condition vectors
h0_crl = c("X0h_no_dox_1","X0h_no_dox_11","X0h_no_dox_12","X0h_no_dox_2")
h12_dox = c("X12h_is_dox_5","X12h_is_dox_6")
h24_dox = c("X24h_is_dox_10","X24h_is_dox_13","X24h_is_dox_14","X24h_is_dox_17",
            "X24h_is_dox_18","X24h_is_dox_19","X24h_is_dox_20","X24h_is_dox_9")
h12_crl = c("X12h_no_dox_3","X12h_no_dox_4")
h24_crl = c("X24h_no_dox_15","X24h_no_dox_16","X24h_no_dox_7","X24h_no_dox_8")

###########################################################################################
# FILTER  - Only genes that had an expression of >1 RPM ...
#         ...for both duplicates of one or more conditions were kept for further analysis.
###########################################################################################
if(TRUE){
  gene_after_filter = c()
  gene_after_filter = unique(c(gene_after_filter,rpm[rpm[,h0_crl] > 1,'gene']))
  print(length(gene_after_filter))
  gene_after_filter = unique(c(gene_after_filter,rpm[rpm[,h12_crl] > 1,'gene']))
  print(length(gene_after_filter))
  gene_after_filter = unique(c(gene_after_filter,rpm[rpm[,h24_crl] > 1,'gene']))
  print(length(gene_after_filter))
  gene_after_filter = unique(c(gene_after_filter,rpm[rpm[,h12_dox] > 1,'gene']))
  print(length(gene_after_filter))
  gene_after_filter = unique(c(gene_after_filter,rpm[rpm[,h24_dox] > 1,'gene']))
  print(length(gene_after_filter))
}

###########################################################################################
# IDENTIFY DEG  - cutoff of Padj = 0.05 and a minimum of 1.5-fold change between both conditions
#               - The dox and no dox conditions were compared separately at 12 and 24 h;
###########################################################################################
count_filtered = count[count$gene %in% gene_after_filter,]
gene = count_filtered$gene
id = rownames(count_filtered)
id2gene = c()
for (i in 1:length(gene)){
  id2gene[id[i]] = gene[i]
}
# 12h comparison
if(TRUE){
  # make analysis data
  count_12h = count_filtered[,c(h12_crl,h12_dox)]
  col_12h = data.frame(row.names = colnames(count_12h),condition = c("CTR","CTR","TREATMENT","TREATMENT"),type=c('CTR_R1','CTR_R2','TREATMENT_R1','TREATMENT_R2'))
  # differentiation analysis
  dds = DESeqDataSetFromMatrix(countData = count_12h, colData = col_12h, design = ~ condition)
  dds = dds[ rowSums(counts(dds)) > 1, ]
  dds$condition = relevel(dds$condition, ref="CTR")
  dds = DESeq(dds)
  res = results(dds)
  resOrdered = res[order(res$padj),]
  #transform into dataframe
  resOrdered = data.frame(resOrdered)
  # add gene names
  resOrdered$gene = id2gene[rownames(resOrdered)]
  # filter
  resOrdered = resOrdered[abs(resOrdered['log2FoldChange']) >= 0.5849625 & resOrdered['padj'] <= 0.05,]
  # write in csv
  write.csv(resOrdered,file="12h_DEG.csv",row.names = FALSE)

}
# 24h comparison
if(TRUE){
  count_24h = count_filtered[,c(h24_crl,h24_dox)]
  col_24h = data.frame(row.names = colnames(count_24h),condition = c(rep("CTR",4),rep("TREATMENT",8)),
                       type=c('CTR_R1','CTR_R2','CTR_R3','CTR_R4',
                              'TREATMENT_R1','TREATMENT_R2','TREATMENT_R3','TREATMENT_R4',
                              'TREATMENT_R5','TREATMENT_R6','TREATMENT_R7','TREATMENT_R8'))
  # differentiation analysis
  dds = DESeqDataSetFromMatrix(countData = count_24h, colData = col_24h, design = ~ condition)
  dds = dds[ rowSums(counts(dds)) > 1, ]
  dds$condition = relevel(dds$condition, ref="CTR")
  dds = DESeq(dds)
  res = results(dds)
  resOrdered = res[order(res$padj),]
  #transform into dataframe
  resOrdered = data.frame(resOrdered)
  # add gene names
  resOrdered$gene = id2gene[rownames(resOrdered)]
  # filter
  resOrdered = resOrdered[abs(resOrdered['log2FoldChange']) >= 0.5849625 & resOrdered['padj'] <= 0.05,]
  # write in csv
  write.csv(resOrdered,file="24h_DEG.csv",row.names = FALSE)
}

###########################################################################################
# MIXED - these two lists of differentially expressed genes were then ...
#       ... pooled to define their kinetics of expression
###########################################################################################
h12_deg = read.csv("12h_DEG.csv")
h24_deg = read.csv("24h_DEG.csv")
pooled = unique(c(h12_deg$gene,h24_deg$gene))
print(length(pooled))

###########################################################################################
# CLASSIFY - Genes were first classified as up- or downregulated using a threshold of ...
#          ... 1.5-fold change ( mean dox mean no dox > 1.5 or < 1 / 1.5, respectively).
###########################################################################################
rpm_deg = rpm[rpm$gene %in% pooled,]
dox_col = c(h12_dox,h24_dox)
crl_col = c(h12_crl,h24_crl)
rpm_deg$mean_dox = rowMeans(rpm_deg[,dox_col])
rpm_deg$mean_crl = rowMeans(rpm_deg[,crl_col])
rpm_deg$ratio = rpm_deg$mean_dox / rpm_deg$mean_crl
rpm_deg$sig = "none"
rpm_deg[rpm_deg['ratio'] > 1.5,]$sig = "up"
rpm_deg[rpm_deg['ratio'] < 2.0/3,]$sig = "down"
up_gene_list = c(rpm_deg[rpm_deg['sig'] == "up",]$gene)
down_gene_list = c(rpm_deg[rpm_deg['sig'] == "down",]$gene)


###########################################################################################
# FURTHER CLASSIFY - classify upregulated genes into early, constant and late activated,
###########################################################################################
## for upregulated genes
# calculate slope
rpm_up = rpm[rpm$gene %in% up_gene_list,]
rpm_up$a_dox = rowMeans(rpm_up[,h12_dox]) - rowMeans(rpm_up[,h0_crl])
rpm_up$a_crl = rowMeans(rpm_up[,h12_crl]) - rowMeans(rpm_up[,h0_crl])
rpm_up$b_dox = rowMeans(rpm_up[,h24_dox]) - rowMeans(rpm_up[,h12_dox])
rpm_up$b_crl = rowMeans(rpm_up[,h24_crl]) - rowMeans(rpm_up[,h12_crl])
rpm_up$a = rpm_up$a_dox - rpm_up$a_crl
rpm_up$b = rpm_up$b_dox - rpm_up$b_crl
# give sig
rpm_up$sig = "constant_up"
rpm_up[rpm_up$a / rpm_up$b > 3,]$sig = "early_up"
rpm_up[rpm_up$b / rpm_up$a > 3,]$sig = "late_up"
rpm_up[rpm_up$a > 0 & rpm_up$b < 0,]$sig = "early_up"
rpm_up[rpm_up$a < 0 & rpm_up$b > 0,]$sig = "late_up"
# get gene
late_up_gene_list = c(rpm_up[rpm_up["sig"] == "late_up",]$gene)
early_up_gene_list = c(rpm_up[rpm_up["sig"] == "early_up",]$gene)
constant_up_gene_list = c(rpm_up[rpm_up["sig"] == "constant_up",]$gene)
length(early_up_gene_list)
length(late_up_gene_list)
length(constant_up_gene_list)
## for downregulated genes
# calculate slope
rpm_down = rpm[rpm$gene %in% down_gene_list,]
rpm_down$a_dox = rowMeans(rpm_down[,h12_dox]) - rowMeans(rpm_down[,h0_crl])
rpm_down$a_crl = rowMeans(rpm_down[,h12_crl]) - rowMeans(rpm_down[,h0_crl])
rpm_down$b_dox = rowMeans(rpm_down[,h24_dox]) - rowMeans(rpm_down[,h12_dox])
rpm_down$b_crl = rowMeans(rpm_down[,h24_crl]) - rowMeans(rpm_down[,h12_crl])
rpm_down$a = rpm_down$a_dox - rpm_down$a_crl
rpm_down$b = rpm_down$b_dox - rpm_down$b_crl
# give sig
rpm_down$sig = "constant_down"
rpm_down[rpm_down$a / rpm_down$b > 3,]$sig = "early_down"
rpm_down[rpm_down$b / rpm_down$a > 3,]$sig = "late_down"
rpm_down[rpm_down$a > 0 & rpm_down$b < 0,]$sig = "early_down"
rpm_down[rpm_down$a < 0 & rpm_down$b > 0,]$sig = "late_down"
# get gene
late_down_gene_list = c(rpm_down[rpm_down["sig"] == "late_down",]$gene)
early_down_gene_list = c(rpm_down[rpm_down["sig"] == "early_down",]$gene)
constant_down_gene_list = c(rpm_down[rpm_down["sig"] == "constant_down",]$gene)
length(early_down_gene_list)
length(late_down_gene_list)
length(constant_down_gene_list)


###########################################################################################
# WRITE INTO CSV
###########################################################################################
if(TRUE){
  write_gene = c(early_up_gene_list)
  write_class = c(rep("early_up",length(early_up_gene_list)))
  write_gene = c(write_gene,early_down_gene_list)
  write_class = c(write_class,rep("early_down",length(early_down_gene_list)))
  write_gene = c(write_gene,late_up_gene_list)
  write_class = c(write_class,rep("late_up",length(late_up_gene_list)))
  write_gene = c(write_gene,late_down_gene_list)
  write_class = c(write_class,rep("late_down",length(late_down_gene_list)))
  write_gene = c(write_gene,constant_up_gene_list)
  write_class = c(write_class,rep("constant_up",length(constant_up_gene_list)))
  write_gene = c(write_gene,constant_down_gene_list)
  write_class = c(write_class,rep("constant_down",length(constant_down_gene_list)))
  df = data.frame(gene = write_gene, class = write_class)
  table(df$class)
  write.csv(df,file="gene_annotation.csv")
  df_filtered = df[!grepl("Rik", df$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("LOC", df_filtered$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("AI", df_filtered$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("Gm", df_filtered$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("Mir", df_filtered$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("AU", df_filtered$gene, ignore.case = FALSE), ]
  df_filtered = df_filtered[!grepl("BC", df_filtered$gene, ignore.case = FALSE), ]
  table(df_filtered$class)
  write.csv(df_filtered,file="gene_annotation_filtered.csv")
}

