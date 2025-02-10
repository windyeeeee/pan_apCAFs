library(Seurat)
library(SeuratObject)
library(BPCells)
library(SeuratDisk)
library(Azimuth)
library(harmony)
library(magrittr)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(colourpicker)
library(colorspace)
library(ggsci)
library(ggthemes)
library(pals)
library(grDevices)
library(RColorBrewer)
library(scales)
library(paletteer)
library(ggtree)
library(treeio)
library(cowplot)


Figure1A{
  # Load Seurat object
  seu_pan <- readRDS("./data/seu_pan.rds")
  
  # Define color palettes for different cell types
  # B Cell
  bc <- c("#00897B")
  
  # Epithelial Cell
  epc <- rev(c("#FFFDE7", "#FFF9C4", "#FFF59D", "#FFF176", "#FFEE58", "#FFE082", "#FFD54F", 
               "#FFCA28", "#FFC107", "#FFB300", "#FFCC80", "#FFB74D", "#FFA726", "#FF9800", 
               "#FB8C00", "#F57C00", "#EF6C00"))
  
  # Myeloid Cell
  myc <- rev(c("#E0F7FA", "#B2EBF2", "#80DEEA", "#4DD0E1", "#26C6DA", "#00BCD4"))
  
  # Endothelial Cell
  enc <- rev(c("#BBDEFB", "#90CAF9", "#64B5F6", "#42A5F5"))
  
  # Fibroblasts
  fic <- rev(c("#FFEBEE", "#FFCDD2", "#EF9A9A", "#E57373", "#EF5350", "#F44336"))
  
  # Mast Cell
  mac <- c("#873C3C")
  
  # Plasma Cell
  plc <- rev(c("#C8E6C9", "#A5D6A7", "#66BB6A", "#43A047", "#388E3C", "#2E7D32", "#1B5E20"))
  
  # T & NK Cell
  tnc <- rev(c("#F9FBE7", "#F0F4C3", "#E6EE9C", "#DCE775", "#D4E157", "#CDDC39", "#C0CA33", "#AFB42B"))
  
  # Combine all color palettes
  color_umap <- c(bc, epc, myc, enc, fic, mac, plc, tnc)
  
  # Prepare UMAP data for plotting
  mat <- data.frame(
    pan_seu@reductions$umap@cell.embeddings, 
    group = pan_seu$seurat_clusters
  )
  
  # Plot UMAP and save as PDF
  pdf("./fig/Figure1A.pdf", width = 8, height = 8)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = group)) +
    geom_point(size = 1e-4, alpha = 0.3) +
    scale_color_manual(values = color_umap) +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.line = element_blank(),      
      axis.text = element_blank(),     
      axis.ticks = element_blank(),      
      axis.title = element_blank()      
    )
  dev.off()
}


Figure1B{
  # Load metadata
  meta <- read.csv("./data/pan_meta.csv", row.names = 1)
  
  # Calculate total cell numbers per disease in Tumor group
  cell_num <- meta %>%
    filter(Group == "Tumor") %>%
    group_by(Disease) %>%
    summarise(count = n()) %>%
    na.omit()
  
  # Calculate cell type counts per disease in Tumor group
  cell_type <- meta %>%
    filter(Group == "Tumor") %>%
    group_by(Type, Disease) %>%
    summarise(count = n()) %>%
    na.omit()
  
  # Merge cell number and cell type data
  df <- merge(cell_num, cell_type, by = "Disease")
  df$Proportion <- df$count.y / df$count.x
  
  # Reorder cell types
  df$Type <- fct_relevel(df$Type, c("Mast cell", "B cell", "Plasma cell", "Endothelial cell", 
                                    "Myeloid cell", "Fibroblasts cell", "T cell", "Epithelial cell"))
  
  # Select relevant columns
  df <- df[, c("Disease", "Type", "Proportion")]
  
  # Reshape data for clustering
  library(reshape2)
  dt <- dcast(df, Disease ~ Type)
  rownames(dt) <- dt[, 1]
  dt <- dt[, -1]
  
  # Perform hierarchical clustering
  tree <- hclust(vegan::vegdist(dt, method = 'bray'), method = 'average')
  
  # Plot phylogenetic tree
  p1 <- ggtree(tree) + 
    geom_tiplab() + 
    xlim(NA, 2)
  
  # Define colors for cell types
  Choose_col <- c("#FFCA28", "#C0CA33", "#EF5350", "#00BCD4", "#90CAF9", "#2E7D38", "#00897B", "#873C3C")
  
  # Plot stacked bar plot
  g1 <- df %>%
    ggplot(aes(y = Proportion, x = Disease, fill = Type, colour = Type)) +
    geom_bar(stat = "identity", width = 0.5, alpha = 0.65) +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_blank()
    ) +
    scale_x_discrete(limits = c("BLCA", "SCC", "PMP", "PDAC", "CESC", "BCC", "BRCA", "COLO", 
                                "STAD", "ESCC", "OV", "ccRCC", "PRAD", "HCC", "NSCLC")) +
    scale_colour_manual(values = rev(Choose_col)) +
    scale_fill_manual(values = rev(Choose_col)) +
    coord_flip() +
    ylab("Frequency") +
    xlab("")
  
  # Combine plots
  g <- ggdraw() +
    draw_plot(p1, 0, 0.06, 0.6, 0.94) +
    draw_plot(g1, 0.2, 0, 0.8, 1)
  
  ggsave("./fig/Figure1B.pdf", g, height = 5, width = 8)
}

Figure1C{

  seu_fb <- readRDS("./data/d_Fb7.rds")
  
  cols <- c("#DCE775", "#E91E63", "#FFE082", "#673AB7", "#4DB6AC", "#EF9A9A", 
            "#A5D6A7", "#FF9800", "#CE93D8", "#AED581", "#3F51B5", "#2196F3", 
            "#FF5722")
  
  g <- DimPlot(
    object = seu_fb, 
    reduction = 'umap', 
    pt.size = 1,          
    raster = FALSE,       
    label = TRUE ) + 
    scale_color_manual(values = cols) +  
    theme_void()               
  
  ggsave("./fig/Figure1C.pdf", g, height = 5, width = 5)
}


Figure1D{
  genes<-c( "COL1A2","PDGFRA","PI16","CXCL8","POSTN",  
            "LRRC15",  "CD74","CD37","CD24","RGS5")             
  
  P<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10")
  
  for (i in 1:length(genes)){
    mat <- data.frame(seu_fb@reductions$umap@cell.embeddings, 
                      group = seu_fb@assays$RNA@data[genes[i],])
    x<-ggplot(mat, aes(x = umap_1, y = umap_2, color = group)) +
      geom_point(size = 1e-3,alpha = 0.1) +
      scale_color_gradient(low = "#e6e4df",high = "#f74343")+
      theme_void()
    assign(P[i],x)
  }
  
  combined_plot <- (p1 + p2 + p3 + p4 + p5) / 
                   (p6 + p7 + p8 + p9 + p10)
  
  ggsave("./fig/Figure1D.pdf", combined_plot, height = 4, width = 12)
}


Figure1E{
  
  fb_dot_data <- read.csv("./data/fb_dot_data.csv",row.names = 1)

  
  genes<-c( "RGS5","HIGD1B","STEAP4","CD36","FABP4",      #C01  RSG5+pericyte     
            "CFD","ADH1B","PI16","WISP2","PLA2G2A",       #C02  PI16+ssCAF
            "FAP","LRRC15", "COL10A1","COL11A1","MMP11",  #C03  LRRC15+myCAF
            "MYH11","PLN","RERGL","BCAM","SORBS2",        #C04  RERGL+SMC  
            "CD74","HLA-DRA","HLA-DPA1","SRGN","CD37",    #C05  CD37+apCAF
            "EZR","ATP1B1","SPP1","CD24","SERPINA1",      #C06 CD24+apCAF
            "DPT","IGF1","ICAM1","RGS2","ARC",            #C07  ICAM1+ssCAF
            "POSTN","PLAT","F3","CCL11","HSD17B2",        #C08  PLAT+myCAF
            "CXCL8","TNFAIP6","CXCL2","HGF","MRPL44",     #C09  HGF+iCAF  
            "IL6","LIF","CSF3","SLC2A1","IL11",           #C10  LIF+iCAF
            "FTX","MEIS1", "PARD3B","FBXL7", "DLG2",      #C11 DLG2+CAF
            "MYLK","CNN1","ACTG2","DES","RAMP1",          #C12 DES+SMC
            "TOP2A","CENPF","PTTG1","UBE2C","MKI67")      #C13 MKI67+CAF

  fb_dot_data$features.plot <- factor(fb_dot_data$features.plot,levels = genes)
  
  
  g <- ggplot(data = fb_dot_data)+ 
    geom_point(aes(x = id,y = features.plot, color = avg.exp.scaled, size = pct.exp))+
    scale_size_area(max_size = 6)+
    scale_color_gradientn(colors = c( "#4a7497",'#cba081', '#EC407A')) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    ylab("") + xlab("") + coord_flip()+  
    theme(panel.background = element_rect(fill = "transparent", color = "black",size = 1), 
                                               legend.key = element_rect(fill = "transparent", color = "transparent"),
                                               text = element_text(color = "black", size = 12,vjust=0.5),
                                               plot.title = element_text(hjust = .5), 
                                               axis.text = element_text(color = "black",size = 12))
  
  ggsave("./fig/Figure1E.pdf", g, height = 4, width = 16)
  
}
  
 
Figure1F{

  signatures<-read_csv("./data/fb_score.csv")
  Genesets<-as.vector(signatures)
  All_genes<-unlist(Genesets)
  
  exprMatrix <- seu_fb@assays$RNA@data
  exprMatrix <- as(exprMatrix, "dgCMatrix")
  meta <- as.data.frame(seu_fb@meta.data)
  
  ##AUC_score
  library(AUCell)
  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
  cells_AUC <- AUCell_calcAUC(Genesets, cells_rankings)
  aucs <- getAUC(cells_AUC)
  
  aucs_t <- data.frame(t(aucs))  
  aucs_t$celltype <- meta[rownames(aucs_t), "seurat_clusters"]
  aucs_t <- aucs_t %>% group_by(celltype) %>% dplyr::summarise_each(funs = mean)
  #write.csv(aucs_t ,"./data/fb_AUCscore.csv")
  
  ####plot
  aucs_t$celltype <- factor(aucs_t$celltype)
  col<-c("#FFCCBC" ,"#FFCC80" ,"#FFE082" ,"#FFF9C4", "#E6EE9C", "#DCEDC8" ,"#C8E6C9", "#B2DFDB", "#80DEEA" ,"#B3E5FC")
  colnames(aucs_t)
  
  p1 <- aucs_t %>% ggplot() + 
    geom_bar(aes(x = celltype, y = Interleukin.Signaling,fill = T), stat = "identity", width = 0.8) + 
    geom_hline(aes(yintercept = median(Interleukin.Signaling)), linetype = "dashed") + 
    ylab(colnames(aucs_t)[2]) +
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[1])  
  p2 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Chemokine.Signaling,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Chemokine.Signaling)),linetype="dashed") + 
    ylab(colnames(aucs_t)[3])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[2])     
  p3 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Extracellular.matrix,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Extracellular.matrix)),linetype="dashed") + 
    ylab(colnames(aucs_t)[4])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[3])    
  p4 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= MHC.II.Antigen.Presentation,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(MHC.II.Antigen.Presentation)),linetype="dashed") + 
    ylab(colnames(aucs_t)[5])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[4])     
  p5 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Angiogenesis,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Angiogenesis)),linetype="dashed") + 
    ylab(colnames(aucs_t)[6])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[5])     
  p6 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Hypoxia,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Hypoxia)),linetype="dashed") + 
    ylab(colnames(aucs_t)[7])+
    xlab("") +
    theme_classic()  +
    scale_fill_manual(values = col[6])    
  p7 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Cell.Proliferation,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Cell.Proliferation)),linetype="dashed") + 
    ylab(colnames(aucs_t)[8])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[7])     
  p8 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= Hedgehog.Signaling,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Hedgehog.Signaling)),linetype="dashed") + 
    ylab(colnames(aucs_t)[9])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[8])     
  p9 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= NF.kappaB.Signaling,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(NF.kappaB.Signaling)),linetype="dashed") + 
    ylab(colnames(aucs_t)[10])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[9])     
  p10 <- aucs_t  %>% ggplot() + geom_bar( aes(x =celltype ,y= TGF.beta.Signaling,fill= T),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(TGF.beta.Signaling)),linetype="dashed") + 
    ylab(colnames(aucs_t)[11])+
    xlab("") +
    theme_classic() +
    scale_fill_manual(values = col[10])     
  

  p <- p1/p2/p3/p4/p5/p6/p7/p8/p9/p10

  ggsave("./fig/Figure1F.pdf",p, height = 24, width =6)
  
}


Figure1G{
  library(pheatmap)
  library(psych)
  exp=AverageExpression(seu_fb)
  exp2<-as.matrix(exp$RNA)
  corr2<-psych::corr.test(x = exp2,y = exp2,method="spearman")
 
  g <- pheatmap::pheatmap(corr2$r)
  ggsave("./fig/Figure1G.pdf",g, height = 6, width =6)
  dev.off()
}


Figure1H{
library(reshape2)
df <- table(seu_fb@meta.data$seurat_clusters,seu_fb@meta.data$Disease) %>% melt()
colnames(df) <- c("Cluster","Group","Number")
write.csv(df,"./data/fb_proportion_data")
df <- na.omit(df)
g <- ggplot(data = df, aes(x = Group, y = Number, fill = Cluster)) +
     geom_bar(stat = "identity", width = 0.8, position = "fill") +
     scale_fill_manual(values = cols) +
     labs(x = "", y = "Ratio") +
     theme_bw() +
     theme(
        panel.grid = element_blank(),  # 移除网格线
        panel.background = element_rect(fill = "transparent", color = "black", size = 1), 
        legend.key = element_rect(fill = "transparent", color = "transparent"),  
        text = element_text(color = "black", size = 12, vjust = 0.5),  
        plot.title = element_text(hjust = 0.5),  
        axis.text = element_text(color = "black", size = 12),  
        axis.text.y = element_text(size = 12, colour = "black"),  
        axis.text.x = element_text(size = 12, colour = "black", hjust = 1, vjust = 1, angle = 45) 
      )
  
ggsave("./fig/Figure1H.pdf",g, height = 6, width =6)
  
}


Figure1I{
  meta<-seu_fb@meta.data
  meta <- read.csv("./data/seu_fb_meta.csv",row.names = 1)
  
  
  dt = meta  %>% group_by(Group,seurat_clusters) %>% dplyr::summarise( count = n())


  test <- reshape2::dcast(dt,Group~seurat_clusters)
  rownames(test)<-test$Group
  test<-test[,-1]
  test <- as.matrix(test)
  
  a=margin.table(test, 1)
  b=margin.table(test, 2)
  c=margin.table(test)
  test = test/(outer(a,b,"*")/c)
  
  tmp <- reshape2::melt(test)
  tmp <- tmp[order(tmp$Var1,tmp$value,tmp$Var2),]
  
  tmp$Var1 <- factor(tmp$Var1,levels=as.character(unique(tmp$Var1)))
  tmp$state <- ifelse(tmp$value>1,"Enrichment","Depletion")
  
  col = c("#F4511E","#1976D2")
  names(col) = c("Enrichment","Depletion")
  
  g <- ggplot(data = tmp)+ 
    geom_point(aes(x = Var1,y = Var2, color = state, size = value))+
    scale_size_area()+
    scale_color_manual(values = col) + 
    cowplot::theme_cowplot() +
    xlab("") +
    ylab("") +
    theme(
      panel.grid = element_blank(),  # 移除网格线
      panel.background = element_rect(fill = "transparent", color = "black", size = 1), 
      legend.key = element_rect(fill = "transparent", color = "transparent"),  
      text = element_text(color = "black", size = 12, vjust = 0.5),  
      plot.title = element_text(hjust = 0.5),  
      axis.text = element_text(color = "black", size = 12),  
      axis.text.y = element_text(size = 12, colour = "black"),  
      axis.text.x = element_text(size = 12, colour = "black", hjust = 1, vjust = 1, angle = 45) 
    )
  
  ggsave("./fig/Figure1I.pdf",g, height = 6, width =3)
  
}

Figure1J{
  library(ggpubr)
  Idents(seu_fb)<-"Group"
  deg_fb <- FindMarkers(seu_fb,ident.1 = "Tumor",ident.2 = "Normal", 
                        test.use = "wilcox",logfc.threshold = 0)
  write.csv(deg_fb,"deg_fb.csv")
  P.Value_t = 1e-10
  
  deg_fb$group = ifelse(deg_fb$p_val_adj < P.Value_t & deg_fb$avg_log2FC <= -1,"down",
                       ifelse(deg_fb$p_val_adj < P.Value_t & deg_fb$avg_log2FC >= 1,"up","stable"))
  deg_fb$name=rownames(deg_fb)
  deg_fb$name[deg_fb$group == "stable"]<-NA
  

  deg_fb$"-log10(p_val_adj)"<- -log10(deg_fb$p_val_adj+1e-300)
  
  p<-ggscatter(deg_fb, x = "avg_log2FC", y = "-log10(p_val_adj)", color = "group",size = 1,
               label = "name", repel = T,palette = c("#39489f", "gray", "#b81f25"))
  g<-p+
    geom_hline(aes(yintercept=10),linetype="dashed")+
    geom_vline(aes(xintercept=-1),linetype="dashed")+
    geom_vline(aes(xintercept=1),linetype="dashed")+
    theme(legend.position = 'none')  + 
    theme(panel.background = element_rect(fill = "transparent", color = "black",size = 1), 
        legend.key = element_rect(fill = "transparent", color = "transparent"),
        text = element_text(color = "black", size = 6,vjust=0.5),
        plot.title = element_text(hjust = .5), 
        axis.text = element_text(color = "black",size = 6))
  
  ggsave("./fig/Figure1J.pdf",g, height = 6, width =6)
}
