setwd("D:/资料/大二第三学期/生物信息/文章复现") # remember to change 
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

early_up_gene = gene_anno[gene_anno['class'] == "early_up",]$gene
constant_up_gene = gene_anno[gene_anno['class'] == "constant_up",]$gene
late_up_gene = gene_anno[gene_anno['class'] == "late_up",]$gene
early_down_gene = gene_anno[gene_anno['class'] == "early_down",]$gene
constant_down_gene = gene_anno[gene_anno['class'] == "constant_down",]$gene
late_down_gene = gene_anno[gene_anno['class'] == "late_down",]$gene

# get mean value
rpm$h0 = rowMeans(rpm[,h0_crl])
rpm$h12_crl = rowMeans(rpm[,h12_crl])
rpm$h24_crl = rowMeans(rpm[,h24_crl])
rpm$h12_dox = rowMeans(rpm[,h12_dox])
rpm$h24_dox = rowMeans(rpm[,h24_dox])

selected_col = c('h0','h12_crl','h24_crl','h12_dox','h24_dox')

# early up
pdf("plot/early_up_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% early_up_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Early") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

# constant up
pdf("plot/constant_up_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% constant_up_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Constant") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

# late up
pdf("plot/late_up_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% late_up_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Late") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

# early down
pdf("plot/early_down_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% early_down_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Early") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

# constant down
pdf("plot/constant_down_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% constant_down_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Constant") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

# late down
pdf("plot/late_down_patttern.pdf", height = 8, width = 8)
if(TRUE){
  rpm_filter = rpm[rpm$gene %in% late_down_gene,c('gene',selected_col)]
  for (i in 1:nrow(rpm_filter)){
    maxn = max(rpm_filter[i,selected_col])
    minn = min(rpm_filter[i,selected_col])
    rpm_filter[i,selected_col] = (rpm_filter[i,selected_col] - minn) / (maxn - minn)
  }
  
  mean_0h = mean(rpm_filter[,'h0'])
  mean_12h_crl = mean(rpm_filter[,'h12_crl'])
  mean_24h_crl = mean(rpm_filter[,'h24_crl'])
  mean_12h_dox = mean(rpm_filter[,'h12_dox'])
  mean_24h_dox = mean(rpm_filter[,'h24_dox'])
  
  mean_df = data.frame(
    mean_0h = c(mean_0h),
    mean_12h_crl = c(mean_12h_crl),
    mean_24h_crl = c(mean_24h_crl),
    mean_12h_dox = c(mean_12h_dox),
    mean_24h_dox = c(mean_24h_dox)
  )
  
  ggplot(rpm_filter) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_crl`, yend=`h24_crl`,color = 'darkblue'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=1, xend=2, y=`h0`, yend=`h12_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(aes(x=2, xend=3, y=`h12_dox`, yend=`h24_dox`,color = 'darkred'), linewidth=1, show.legend=F,alpha=0.2) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_crl, yend=mean_24h_crl,color = 'blue'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=1, xend=2, y=mean_0h, yend=mean_12h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    geom_segment(data = mean_df,aes(x=2, xend=3, y=mean_12h_dox, yend=mean_24h_dox,color = 'red'), linewidth=3, show.legend=F,alpha=1) +
    scale_color_manual(values = c('darkblue'="#4D3DA1", 'darkred'="#9D3124","blue"="#41329B","red"="#DC2F1E")) +
    xlim(1,3) +
    ylim(0,1) + 
    scale_x_continuous(breaks = c(1, 2, 3), labels = c("0", "12", "24")) +
    scale_y_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1.0")) +
    xlab("Time(h)") +  
    ylab("Normalized expression") +
    labs(title = "Late") +
    theme_classic() + 
    theme(axis.text.x = element_text(size = 20,color = "black"),  
          axis.title.x = element_text(size = 25,color = "black"),
          axis.text.y = element_text(size = 20,color = "black"),  
          axis.title.y = element_text(size = 25,color = "black"),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 30, hjust = 0.5)) 
}
dev.off()

