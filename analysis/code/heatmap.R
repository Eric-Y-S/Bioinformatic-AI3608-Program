setwd("D:/资料/大二第三学期/生物信息/文章复现")
library(ggplot2)

rpm = read.csv("samples_RPM.csv",sep = ",",header = 1)
rpm = rpm[, -which(names(rpm) == "X")]
gene_anno = read.csv("gene_annotation_filtered.csv")
gene_anno = gene_anno[, -which(names(gene_anno) == "X")]

h0_crl = c("X0h_no_dox_1","X0h_no_dox_11","X0h_no_dox_12","X0h_no_dox_2")
h12_dox = c("X12h_is_dox_5","X12h_is_dox_6")
h24_dox = c("X24h_is_dox_10","X24h_is_dox_13","X24h_is_dox_14","X24h_is_dox_17",
            "X24h_is_dox_18","X24h_is_dox_19","X24h_is_dox_20","X24h_is_dox_9")
h12_crl = c("X12h_no_dox_3","X12h_no_dox_4")
h24_crl = c("X24h_no_dox_15","X24h_no_dox_16","X24h_no_dox_7","X24h_no_dox_8")

up_gene = gene_anno[grepl("up",gene_anno$class),]$gene
rpm_filter = rpm[rpm$gene %in% up_gene,]
rpm_filter$mean_h12_dox = rowMeans(rpm_filter[,h12_dox])
rpm_filter$mean_h24_dox = rowMeans(rpm_filter[,h24_dox])
rpm_filter$ratio = rpm_filter$mean_h24_dox / rpm_filter$mean_h12_dox
rpm_filter_sort <- rpm_filter[order(rpm_filter$ratio), ]

# standization
rpm_filter_sort[h12_dox] = scale(rpm_filter_sort[h12_dox], center = TRUE, scale = TRUE) 
rpm_filter_sort[h12_crl] = scale(rpm_filter_sort[h12_crl], center = TRUE, scale = TRUE) 
rpm_filter_sort[h24_dox] = scale(rpm_filter_sort[h24_dox], center = TRUE, scale = TRUE) 
rpm_filter_sort[h24_crl] = scale(rpm_filter_sort[h24_crl], center = TRUE, scale = TRUE) 
rpm_filter_sort$h12_dox = rowMeans(rpm_filter_sort[h12_dox])
rpm_filter_sort$h12_crl = rowMeans(rpm_filter_sort[h12_crl])
rpm_filter_sort$h24_dox = rowMeans(rpm_filter_sort[h24_dox])
rpm_filter_sort$h24_crl = rowMeans(rpm_filter_sort[h24_crl])

early_up_list = gene_anno[gene_anno['class']=="early_up",]$gene
constant_up_list = gene_anno[gene_anno['class']=="constant_up",]$gene
late_up_list = gene_anno[gene_anno['class']=="late_up",]$gene
# define x and y
x = c('h12_crl','h24_crl','h12_dox','h24_dox')
x = factor(x,levels = x)
y = factor(rpm_filter_sort$gene,levels = rpm_filter_sort$gene)
# transform data
var1 = c(); var2 = c(); value = c(); class = c()
for (i in x){
  for (j in y){
    var1 = c(var1,i)
    var2 = c(var2,j)
    value = c(value,rpm_filter_sort[rpm_filter_sort['gene'] == j,i])
    class = c(class,gene_anno[gene_anno['gene']==j,'class'])
  }
}
df = data.frame(x = var1,y = var2,value = value,class = class)
df$x = factor(df$x, levels = x)
df$y = factor(df$y, levels = y)
color_map = c('early_up' = "darkgreen", 'constant_up' = "orange", "late_up" = "darkred")
df$co = color_map[df$class]

# plot heatmap
scale_c = c(-1.5,1.5)
ggplot(data = df, aes(x = x, y = y, fill= value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#0E46F0", high = "#902135", mid = "white",
                       midpoint = mean(scale_c), limit = scale_c, space = "Lab",
                       name="z-score") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())
        
ggplot(data = df, aes(x = x, y = y, fill= value)) + 
  geom_tile(aes(x=1, y=y, fill=co), width=0.1, height=1) +
  scale_fill_manual(values = c('darkgreen'="darkgreen", 'orange'="orange","darkred"="darkred")) +
  geom_text(data = subset(df, y == 'Lhx1'),
            aes(x = x, y = y, label = 'y'), color = "black",
            vjust = -0.5, size = 3) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
  ) 


##############################################
scale_c = c(-1.5,1.5)
ggplot(data = df[df$y %in% early_up_list,], aes(x = x, y = y, fill= value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#0E46F0", high = "#902135", mid = "white", 
                       midpoint = mean(scale_c), limit = scale_c, space = "Lab", 
                       name="xxx") +
  #geom_text(aes(x = 5, y = seq_along('Lhx1'), label = 'Lhx1'), 
  #          hjust = 0, size = 3) +
  geom_text(data = subset(df, y == 'Lhx1'),
            aes(x = x, y = y, label = 'y'), color = "black",
            vjust = -0.5, size = 3) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
  ) 




  
ggplot(data = df[df$y %in% constant_up_list,], aes(x = x, y = y, fill= value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#0E46F0", high = "#902135", mid = "white", 
                       midpoint = mean(scale_c), limit = scale_c, space = "Lab", 
                       name="xxx") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        ) 

