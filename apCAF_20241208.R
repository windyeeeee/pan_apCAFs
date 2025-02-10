rm(list=ls())
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("mojaveazure/seurat-object", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
remotes::install_github("bnprks/BPCells", quiet = TRUE)


package{
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
}


style{
  
  
discrete_colors{
  
# colors_dt<-read.csv("./color/color.csv")
# colors_list<-as.list(colors_dt)
# saveRDS(colors_list,"colors_list.rds")
colors_list<-readRDS("./color/colors_list.rds")
    
choose_col<-function(a,b,c){
            col<-c()
            s<-sample(c(a,b,c),length(colors_list),replace = T) 
            for(i in 1:length(colors_list)){
            S<-s[i]
            col<-append(col,colors_list[[i]][S])
            }
            return(col)
            }

Choose_col<-choose_col(3,4,5)
col2<-choose_col(2,1,3)
scales::show_col(col2)
show_col(colors_list[[3]])
    
# col_map<-as.list(as.data.frame(t(colors_dt)))
# saveRDS(col_map,"col_map.rds")
col_map<-readRDS("./colors/col_map.rds")
scales::show_col(col_map$V8)

COL_A<-c()
for(i in 1:19){
    COL_A<-append(COL_A,sample(colors_list[[i]],3))
    }
COL_A<-sample(COL_A,55)
scales::show_col(COL_A)  

}
  

contuous_colors{
v <- ggplot(faithfuld, aes(waiting, eruptions, fill = density))+geom_tile()
v
dev.off()

library(ggplot2)    
v+ scale_fill_gradient()
v+ scale_fill_gradient2(low = "#FFF8E1", high ="#FF6F00")
v+ scale_fill_gradientn(colors=c("#FFF8E1","gray","#FF6F00"))
v+ scale_fill_viridis_c(option = "H",alpha = 1) #ABCDEFGH
v+ scale_fill_binned(type = "viridis")

library(viridis)
v+ scale_fill_viridis("Temperature", option = "H",alpha = 1,discrete = F)  
v+ scale_fill_viridis(option = "H",alpha = 0.5,discrete = F)  
# "magma" (or "A")
# "inferno" (or "B")
# "plasma" (or "C")
# "viridis" (or "D")
# "cividis" (or "E")
# "rocket" (or "F")
# "mako" (or "G")
# "turbo" (or "H")

library(colorspace)
v + scale_fill_continuous_sequential("red")
v + scale_fill_continuous_sequential(palette="rdpu")
v + scale_fill_continuous_diverging(palette="Purple-Brown")
v + scale_fill_continuous_divergingx(palette="Temps")
v + scale_fill_continuous_divergingx(palette="Earth")
v + scale_fill_continuous_divergingx(palette="Fall")
v + scale_fill_continuous_divergingx(palette="Roma")
v + scale_fill_continuous_qualitative("Warm")

  }
  
  
themes{
library(ggthemes) 
  
p1+theme_void()           #***
p1+theme_minimal()        #***
p1+theme_classic()        #***
p1+theme_bw()             #***
  
p1+ggthemes::theme_base()   #****
p1+ggthemes::theme_par()    #****
p1+ggthemes::theme_tufte()  #****
p1+ggthemes::theme_few()    #****
p1+ggthemes::theme_clean()  #****
p1+ggthemes::theme_map()    #****
p1+ggthemes::theme_pander() #****
p1+ggthemes::theme_solid()  #****
p1+ggthemes::theme_wsj()    #****
p1+ggthemes::theme_calc()   #****
p1+theme_cxf  #*****
p1+theme_cxf2 #*****  
  
theme_cxf<-  theme(panel.background = element_rect(fill = "transparent", color = "black",size = 1), 
                   legend.key = element_rect(fill = "transparent", color = "transparent"),
                   text = element_text(color = "black", size = 6,vjust=0.5),
                   plot.title = element_text(hjust = .5), 
                   axis.text = element_text(color = "black",size = 6))
theme_cxf2 <- theme_classic() +
          theme(panel.background=element_rect(fill='transparent', color='black',size = 1.5), 
          text = element_text(color = "black", size = 16,vjust=0.5),
          panel.border=element_rect(fill='transparent', color='transparent'), 
          panel.grid=element_blank(), axis.title = element_text(color='black', vjust=0.5),
          axis.text = element_text(color="black"),
          axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'), 
          axis.ticks.margin = unit(0.2,"lines"), legend.title=element_blank(), 
          plot.title = element_text(family = 'serif',face = 'bold',colour = 'black',size = 20,hjust = .5), 
          plot.subtitle = element_text(family = 'serif',face = 'bold',colour = 'black',size = 16,hjust = .5),
          legend.key=element_rect(fill='transparent', color='transparent')) 
  
}  


layout{
  
  library("gridExtra")
  grid.arrange(p1_jama, p2_jama, ncol = 2)   
  
  library(patchwork)  
  #+/|
  
}

  
}


Function{
  sig_test <- function(data,factors,values){
    compare_list <- unique(data[[factors]]) %>% combn(2) 
    compare_list <- t(compare_list)
    wil_test <- data.frame("comp1" = c(), "comp2" = c(), "p.value" = c(), "position" = c() )
    for(i in 1:dim(compare_list)[1] ){
      a = data[data["Group"] == compare_list[i,][1],values]
      b = data[data["Group"] == compare_list[i,][2],values]
      wilcox_test = wilcox.test(x = a, y = b)
      p.value = round(wilcox_test$p.value, 5)
      position = max(a,b, na.rm = T)
      tmp = data.frame("comp1" = compare_list[i,][1], "comp2" = compare_list[i,][2], "p.value" = p.value ,  "position" = position)
      wil_test = rbind(wil_test, tmp)
    }
    return(wil_test)
    
  }
  
  
  
  }

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
  P.Value_t = 1e-100
  
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



seu_apCAF <- readRDS("./data/d_apCAF2.rds")
  
Figure2A{
  load("./data/seu_apCAF.Rdata")

  g <- DimPlot(seu_apCAF,cols = c("#e7a85e","#9aca94","#a8c946","#ddd244"), pt.size =1.5)+theme_void()
  
  ggsave("./fig/Figure2A.pdf",g, height = 2, width =4)
  
}
  
Figure2B{
  load("./data/seu_apCAF.Rdata")
  
 
  f1 <- FeaturePlot(seu_apCAF,features = c("MSLN"), pt.size = 1,alpha = 0.5,max.cutoff = 3)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()

  f2 <- FeaturePlot(seu_apCAF,features = c("UPK3B"), pt.size = 1,alpha = 0.5,max.cutoff = 4)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()

  f3 <- FeaturePlot(seu_apCAF,features = c("KRT19"), pt.size = 1,alpha = 0.5,max.cutoff = 4)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()
  
  f4 <- FeaturePlot(seu_apCAF,features = c("CD74"), pt.size = 1,alpha = 0.5,max.cutoff = 3)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()
  
  
  f5 <- FeaturePlot(seu_apCAF,features = c("PTPRC"), pt.size = 1,alpha = 0.5,max.cutoff = 4)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()
  
  
  f6 <- FeaturePlot(seu_apCAF,features = c("CD52"), pt.size = 1,alpha = 0.3,max.cutoff = 4)+
        scale_color_gradientn(colors = c("gray90", '#E31A1C'))+
        theme_void()
  
 
  g <- f1 + f2 + f3 + f4 + f5 + f6 +
       patchwork::plot_layout(nrow = 2, ncol = 3)
  
  ggsave("./fig/Figure2B.pdf",g, height = 4, width =12)
  
}

Figure2C{
  library(monocle3)
  library(SeuratWrappers)
  
  load("./data/seu_apCAF.Rdata")
  
 
  Idents(seurat_Fb)
  seurat_FB_M <- subset(seurat_Fb,idents=c(0,3))
  seurat_FB_F <- subset(seurat_Fb,idents=c(1,2))  
  DimPlot(seurat_FB_M)
  seurat_FB_M$RNA_snn_res.0.1<-droplevels(seurat_FB_M$RNA_snn_res.0.1)
  seurat_FB_F$RNA_snn_res.0.1<-droplevels(seurat_FB_F$RNA_snn_res.0.1)
  seurat_Fb<-seurat_FB_M
  seurat_Fb<-seurat_FB_F
  
  
  monocle <- as.cell_data_set(seu_apCAF, assay = "RNA")
  monocle <- estimate_size_factors(monocle)
   
  
  #Get cell metadata
  colData(monocle)
  
  #Get gene metadata
  fData(monocle)
  fData(monocle)$gene_short_name <- rownames(fData(monocle))
  
  #Get counts
  counts(monocle)
  
  #Assign partitions
  recreate.partiion <- c(rep(1,length(monocle@colData@rownames)))
  names(recreate.partiion) <- monocle@colData@rownames
  recreate.partiion <- as.factor(recreate.partiion)
  monocle@clusters$UMAP$partitions <- recreate.partiion
  
  
  #Assign cluster info
  list_cluster <- seu_apCAF$seurat_clusters
  monocle@clusters$UMAP$clusters <- list_cluster
  
  #Save umap structure
  monocle@int_colData@listData$reducedDims$UMAP <- seu_apCAF@reductions$umap@cell.embeddings
  
  #Learn trajectory graph
  monocle <- learn_graph(monocle, use_partition = F)
  plot_cells(monocle, color_cells_by = "cluster", label_groups_by_cluster = F, label_branch_points = F,
             label_roots = F, label_leaves = F, group_label_size = 0) + theme_bw()
  
  #Oder cells in pseudotime
  get_earliest_principal_node <- function(monocle, time_bin = c("0")) {
    cell_ids <- which(colData(monocle)[, "seurat_clusters"] == time_bin)
    closest_vertex <- monocle@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(monocle), ])
    root_pr_nodes <- igraph::V(principal_graph(monocle)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
    return(root_pr_nodes)
  }
  
  root_pr_nodes <- get_earliest_principal_node(monocle, time_bin = c("0"))
  
  monocle <- order_cells(monocle, root_pr_nodes = root_pr_nodes)
  remotes::install_github("rstudio/htmltools")
  monocle <- order_cells(monocle)    
  
  #Visualization
  plot_cells(monocle, color_cells_by = "pseudotime", label_groups_by_cluster = F, label_branch_points = F,
             label_roots = F, label_leaves = F, show_trajectory_graph = F) + theme_cxf
  
  pseudotime(monocle) # cells ordered by pseudotime
  monocle$monocle3_pseudotime <- pseudotime(monocle)
  data_pseudo <- as.data.frame(colData(monocle))
  ggplot(data_pseudo, aes(monocle3_pseudotime, reorder(RNA_snn_res.0.1, monocle3_pseudotime, median), fill = RNA_snn_res.0.1)) + 
    geom_boxplot() + theme_minimal()
  dev.off()

  genes_T <- graph_test(monocle, neighbor_graph = "principal_graph", cores = 4)
  #top genes
  topgenes <- genes_T %>%
    arrange(q_value) %>%
    filter(status == "OK") 
  
  A <- topgenes %>% arrange(desc(morans_test_statistic))   %>% rownames()

  A<-c("MSLN","PTPRC","HLA-DRA","SPP1","CD74")
  plot_genes_in_pseudotime(monocle[rowData(monocle)$gene_short_name %in% A],
                           color_cells_by="pseudotime",ncol=1, min_expr=0, cell_size=1) 
  dev.off()
  
}

seu_apCAF <- readRDS("./data/d_apCAF2.rds")
load("./data/seu_apCAF.Rdata")
Idents(seu_apCAF)
F_apCAF <- subset(seu_apCAF,idents=c(1,2)) 
M_apCAF <- subset(seu_apCAF,idents=c(0,3))
DimPlot(F_apCAF)
levels(F_apCAF$seurat_clusters)
table(F_apCAF$seurat_clusters)
F_apCAF$seurat_clusters <- droplevels(F_apCAF$seurat_clusters)
M_apCAF$seurat_clusters <- droplevels(M_apCAF$seurat_clusters)
Idents(F_apCAF)
deg_f <- FindMarkers(F_apCAF, ident.1 = "2", ident.2 = "1", test.use = "wilcox")
deg_m <- FindMarkers(M_apCAF, ident.1 = "3", ident.2 = "0", test.use = "wilcox")

write.csv(deg_f,"./data/deg_f.csv")
write.csv(deg_m,"./data/deg_m.csv")



Figure2D{
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
  
  deg_f <-  read.csv("./data/deg_f.csv",row.names = 1)
  
  geneIDannotation <- function(geneLists=c(1,2,9),name=T,map=T,ensemble=F,accnum=F){
    ## input ID type : So far I just accept entrezID or symbol
    ## default, we will annotate the entrezID and symbol, chromosone location and gene 
    suppressMessages(library("org.Hs.eg.db"))
    all_EG=mappedkeys(org.Hs.egSYMBOL) 
    EG2Symbol=toTable(org.Hs.egSYMBOL)
    if( all(! geneLists %in% all_EG) ){
      inputType='symbol'
      geneLists=data.frame(symbol=geneLists)
      results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
    }else{
      inputType='entrezID'
      geneLists=data.frame(gene_id=geneLists)
      results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
    }
    
    if ( name ){
      EG2name=toTable(org.Hs.egGENENAME)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }
    if(map){
      EG2MAP=toTable(org.Hs.egMAP)
      results=merge(results,EG2MAP,by='gene_id',all.x=T)
    }
    if(ensemble){
      EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
      results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
    }
    if(accnum){
      EG2accnum=toTable(org.Hs.egREFSEQ) 
      results=merge(results,EG2MAP,by='gene_id',all.x=T)
    }
    return(results)
  }
  
  genes <- rownames(deg_f)
  deg_f$symbol <- rownames(deg_f) 
  colnames(deg_f)
  deg_f <- deg_f[,c("symbol","avg_log2FC")]
  colnames(deg_f)<-c("symbol","log2FoldChange")
  genelist<-geneIDannotation(genes)
  genelist_merge<-merge(deg_f,genelist,by="symbol")
  DEG_gene<-genelist_merge[,c(3,2)]
  DEG<-DEG_gene[,2]
  names(DEG)<-as.character(DEG_gene[,1])
  DEG<-sort(DEG,decreasing = T)
  
  KEGG_gsea <- gseKEGG(DEG,nPerm = 1000, minGSSize = 10, maxGSSize = 1000,pvalueCutoff = 0.05)
  g<-gseaplot2(KEGG_gsea,c("hsa04610","hsa04060","hsa04380","hsa05235","hsa03320","hsa00380"),pvalue_table = T)
  ggsave("./fig/Figure2D.pdf",g, height = 4, width =8)
}

Figure2E{
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
  
  deg_m <-  read.csv("./data/deg_m.csv",row.names = 1)
  
  geneIDannotation <- function(geneLists=c(1,2,9),name=T,map=T,ensemble=F,accnum=F){
    ## input ID type : So far I just accept entrezID or symbol
    ## default, we will annotate the entrezID and symbol, chromosone location and gene 
    suppressMessages(library("org.Hs.eg.db"))
    all_EG=mappedkeys(org.Hs.egSYMBOL) 
    EG2Symbol=toTable(org.Hs.egSYMBOL)
    if( all(! geneLists %in% all_EG) ){
      inputType='symbol'
      geneLists=data.frame(symbol=geneLists)
      results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
    }else{
      inputType='entrezID'
      geneLists=data.frame(gene_id=geneLists)
      results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
    }
    
    if ( name ){
      EG2name=toTable(org.Hs.egGENENAME)
      results=merge(results,EG2name,by='gene_id',all.x=T)
    }
    if(map){
      EG2MAP=toTable(org.Hs.egMAP)
      results=merge(results,EG2MAP,by='gene_id',all.x=T)
    }
    if(ensemble){
      EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
      results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
    }
    if(accnum){
      EG2accnum=toTable(org.Hs.egREFSEQ) 
      results=merge(results,EG2MAP,by='gene_id',all.x=T)
    }
    return(results)
  }
  
  genes <- rownames(deg_m)
  deg_f$symbol <- rownames(deg_m) 
  colnames(deg_m)
  deg_f <- deg_m[,c("symbol","avg_log2FC")]
  colnames(deg_m)<-c("symbol","log2FoldChange")
  genelist<-geneIDannotation(genes)
  genelist_merge<-merge(deg_m,genelist,by="symbol")
  DEG_gene<-genelist_merge[,c(3,2)]
  DEG<-DEG_gene[,2]
  names(DEG)<-as.character(DEG_gene[,1])
  DEG<-sort(DEG,decreasing = T)
  
  KEGG_gsea <- gseKEGG(DEG,nPerm = 1000, minGSSize = 10, maxGSSize = 1000,pvalueCutoff = 0.05)
  g<-gseaplot2(KEGG_gsea,c("hsa04612","hsa04514","hsa04610","hsa04060","hsa04657","hsa04080"),pvalue_table = T)
  ggsave("./fig/Figure2E.pdf",g, height = 4, width =8)
}

Figure2F{
    fb<-readRDS("./data/seu_fb.rds")
    fb_5 <-subset(fb,idents = "c05")
    fb_6 <-subset(fb,idents = "c06")
    Idents(fb_5)<-"Group"
    deg_fb5 <- FindMarkers(fb_5,ident.1 = "Tumor",ident.2 = "Normal", test.use = "wilcox",logfc.threshold = 0)
    Idents(fb_6)<-"group"
    deg_fb6 <- FindMarkers(fb_6,ident.1 = "Tumor",ident.2 = "Normal", test.use = "wilcox",logfc.threshold = 0)

    
    degdf<-deg_fb5   #deg_fb6
    logFC_t=0
    P.Value_t = 1e-20
    degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC <= -1,"down",
                         ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC >= 1,"up","stable"))
    degdf$name=rownames(degdf)
    degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
    
    library(ggpubr)
    attach(degdf)
    
    p<-ggscatter(degdf, x = "avg_log2FC", y = "-log10(p_val_adj)", color = "group",size = 1,
                 label = "name", repel = T, 
                 label.select = c("C1QC","RARRES1","APOC1","SPP1"),
                 #label.select = c("CD24","CA9","EGLN3","SPP1"),
                 palette = c("#39489f", "gray", "#b81f25"))
    p+
      geom_hline(aes(yintercept=20),linetype="dashed")+
      geom_vline(aes(xintercept=-1),linetype="dashed")+
      geom_vline(aes(xintercept=1),linetype="dashed")+
      theme_cxf+
      theme(legend.position = 'none')   
    
    ggsave('./fig/Figure2F.pdf',height = 3, width = 6)
    
}

Figure2G{    
    
    genes <- c('C1QC', 'APOC1','RARRES1', 'CD24', 'CA9','EGLN3', 'SPP1')
    
    VlnPlot(Fb,
            features = genes,
            layer = 'data',
            assay = 'RNA',
            log = T,
            add.noise = F,
            pt.size = 0,
            raster =T,
            cols =  Choose_col)
 
    dev.off()
    
  }  
  
Figure2I-K{
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(KernSmooth)
  library(RColorBrewer)
  library(plotly)
  library(BiocParallel)
  library(pheatmap)
  
  library(cowplot)
  library(ggpubr)
  library(ggsci)
  library(ggplot2)
  library(tidygraph)
  library(ggraph)
  library(stringr)
  
  
  
  colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                        pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                        pal_locuszoom("default")(7),pal_igv("default")(51),
                        pal_uchicago("default")(9),pal_startrek("uniform")(7),
                        pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                        pal_simpsons("springfield")(16),pal_gsea("default")(12)))
  len <- 100
  incolor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))
  
  plot_pyscenic <- function(inloom='aucell.loom',incolor=incolor,inrss='seurat_annotations_rss.rds',inrds='subset.rds',infun='median', ct.col='seurat_annotations',inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50){

    loom <- open_loom(inloom)
    
    regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
    regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
    regulonAucThresholds <- get_regulon_thresholds(loom)
   
    close_loom(loom)
    
    rss <- readRDS(inrss)
    sce <- readRDS(inrds)
    embeddings <- sce@reductions$umap

    df = do.call(rbind,
                 lapply(1:ncol(rss), function(i){
                   dat= data.frame(
                     regulon  = rownames(rss),
                     cluster =  colnames(rss)[i],
                     sd.1 = rss[,i],
                     sd.2 = apply(rss[,-i], 1, get(infun))
                   )
                 }))
    
    df$fc = df$sd.1 - df$sd.2
    

    ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
    
    ntopgene <- unique(ntopg$regulon)
    write.table(ntopgene,'sd_regulon_RSS.list',sep='\t',quote=F,row.names=F,col.names=F)
   

    rssPlot <- plotRSS(rss)
    regulonsToPlot <- rssPlot$rowOrder
    rp_df <- rssPlot$df
    
    write.table(regulonsToPlot,'rss_regulon.list',sep='\t',quote=F,row.names=F,col.names=F)
    write.table(rp_df,'rssPlot_data.xls',sep='\t',quote=F)
    nlen <- length(regulonsToPlot)
    hei <- ceiling(nlen)*0.4
    blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    lgroup <- levels(rssPlot$df$cellType)
    
    nlen2 <- length(lgroup)
    wei <- nlen2*2
    pdf(paste0('regulons_RSS_',ct.col,'_in_dotplot.pdf'),wei,hei)
    print(rssPlot$plot)
    dev.off()
    
    anrow = data.frame( group = ntopg$cluster)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(anrow$group)
    annotation_colors <- list(group=lcolor)
    
    pn1 = rss[ntopg$regulon,]
    pn2 = rss[unique(ntopg$regulon),]
    rownames(pn1) <-  make.unique(rownames(pn1))
    rownames(anrow) <- rownames(pn1)
    scale='row'
    hei <- ceiling(length(ntopg$regulon)*0.4)
    pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
    )
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons')
    )
    dev.off()
    
    
    pn2 = rss[unique(rp_df$Topic),]
    scale='row'
    hei <- ceiling(length(unique(rp_df$Topic))*0.4)
    pdf(paste0('regulon_RSS_in_plotRSS_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='plotRSS unique regulons')
    )
    dev.off()
    
    
    hei <- ceiling(length(rownames(rss))*0.2)
    pdf(paste0('all_regulons_RSS_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(rss,scale=scale,show_rownames = T,main='all regulons RSS')
    )
    dev.off()
   
    if (is.null(inregulons)){
      inregulons <- regulonsToPlot
    }else{
      inregulons <- intersect(inregulons,rownames(rss))
      regulonsToPlot <- inregulons
      
    }
    pn3=as.matrix(regulonAUC@assays@data$AUC)
    regulon <- rownames(pn3)
   
    pn3 <- pn3[regulon,]
   
    sce$group1=sce@meta.data[,ct.col]
    
    meta <- sce@meta.data
    meta <- meta[order(meta$group1),]
   
    ancol = data.frame(meta[,c('group1')])
    colnames(ancol) <- c('group1')
    rownames(ancol) <- rownames(meta)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(ntopg$cluster)
    annotation_colors <- list(group1 =lcolor)
    
    df1 <- ancol
    df1$cell <- rownames(df1)
    df1 <- df1[order(df1$group1),]
    pn3 <- pn3[,rownames(df1)]
    torange=c(-2,2)
    pn3 <- scales::rescale(pn3,to=torange)
    pn3 <- pn3[,rownames(ancol)]
    
    scale='none'
    hei <- ceiling(length(unique(regulon))*0.2)
    pdf(paste0('all_regulon_activity_in_allcells.pdf'),10,hei)
    print(
      pheatmap(pn3,annotation_col = ancol,scale=scale,annotation_colors=annotation_colors,show_rownames = T,show_colnames = F,cluster_cols=F)
    )
   
    dev.off()
    
   
    regulonsToPlot = inregulons
    sce$sub_celltype <- sce@meta.data[,ct.col]
    sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
    
    cellClusters <- data.frame(row.names = colnames(sce),
                               seurat_clusters = as.character(sce$seurat_clusters))
    cellTypes <- data.frame(row.names = colnames(sce),
                            celltype = sce$sub_celltype)
    
    sce@meta.data = cbind(sce@meta.data ,t(sub_regulonAUC@assays@data@listData$AUC[regulonsToPlot,]))
    Idents(sce) <- sce$sub_celltype
    
    nlen <- length(regulonsToPlot)
    hei <- ceiling(nlen)*0.4
    blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    nlen2 <- length(unique(sce$sub_celltype))
    wei <- nlen2*2
    pdf('regulons_activity_in_dotplot.pdf',wei,hei)
    print(DotPlot(sce, features = unique(regulonsToPlot)) + coord_flip()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
            scale_color_gradientn(colours = blu)
    )
    dev.off()
    
    
    hei=ceiling(nlen/4)*4
    pdf('regulons_activity_in_umap.pdf',16,hei)
    print(RidgePlot(sce, features = regulonsToPlot , ncol = 4))
    print(VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ))
    print(FeaturePlot(sce, reduction="umap",features =regulonsToPlot))
    dev.off()
    
    grn <- read.table(ingrn,sep='\t',header=T,stringsAsFactors=F)
    inregulons1=gsub('[(+)]','',inregulons)
    
    c1 <- which(grn$TF %in% inregulons1)
    grn <- grn[c1,]
  
    pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
    for (tf in unique(grn$TF)) {
      tmp <- subset(grn,TF==tf)
      if (dim(tmp)[1] > ntop2) {
        tmp <- tmp[order(tmp$importance,decreasing=T),]
        tmp <- tmp[1:ntop2,]
      }
      node2 <- data.frame(tmp$target)
      node2$node.size=1.5
      node2$node.colour <- 'black'
      colnames(node2) <- c('node','node.size','node.colour')
      df1 <- data.frame(node=tf,node.size=2,node.colour='#FFDA00')
      node2 <- rbind(df1,node2)
      
      
      edge2 <- tmp
      colnames(edge2) <- c('from','to','edge.width')
      edge2$edge.colour <- "#1B9E77"
      torange=c(0.1,1)
      edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)
      
      graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
      p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
        scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
        geom_node_label(aes(label = node,colour = node.colour),size = 3.5, repel = TRUE)
      p1 <- p1 + scale_color_manual(values=c('#FFDA00','black'))+scale_edge_color_manual(values=c("#1B9E77"))
      print(p1)
    }
    dev.off()

    meta <- sce@meta.data
    celltype <- ct.col
    cellsPerGroup <- split(rownames(meta),meta[,celltype])
    sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
  
    regulonActivity_byGroup <- sapply(cellsPerGroup,
                                      function(cells)
                                        rowMeans(getAUC(sub_regulonAUC)[,cells]))
    scale='row'
    rss <- regulonActivity_byGroup
    hei <- ceiling(length(regulonsToPlot)*0.4)
    pn1 <- rss[regulonsToPlot,]
    pdf(paste0('regulon_activity_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,scale=scale,show_rownames = T, main='regulons activity')
    )
    dev.off()
    
    hei <- ceiling(length(rownames(rss))*0.2)
    pdf(paste0('all_regulons_activity_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(rss,scale=scale,show_rownames = T,main='all regulons activity')
    )
    dev.off()
    
  
    df = do.call(rbind,
                 lapply(1:ncol(rss), function(i){
                   dat= data.frame(
                     regulon  = rownames(rss),
                     cluster =  colnames(rss)[i],
                     sd.1 = rss[,i],
                     sd.2 = apply(rss[,-i], 1, get(infun))
                   )
                 }))
    
    df$fc = df$sd.1 - df$sd.2
    
  
    ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
    
    ntopgene <- unique(ntopg$regulon)
    write.table(ntopgene,'sd_regulon_activity.list',sep='\t',quote=F,row.names=F,col.names=F)
    
    anrow = data.frame( group = ntopg$cluster)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(anrow$group)
    annotation_colors <- list(group=lcolor)
    pn1 = rss[ntopg$regulon,]
    pn2 = rss[unique(ntopg$regulon),]
    rownames(pn1) <-  make.unique(rownames(pn1))
    rownames(anrow) <- rownames(pn1)
    scale='row'
    hei <- ceiling(length(ntopg$regulon)*0.4)
    pdf(paste0('regulon_activity_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
    )
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons ')
    )
    dev.off()
    
  }
}  

Figure2L-M{
  library("ggpubr")
  
  meta <- read.csv("./data/seu_fb_meta.csv",row.names = 1)
  
  dt <- meta[meta$seurat_clusters %in% "c05",] #c06
  
  dt <- dt[dt$Group %in% "Tumor",]
  
  dt = dt  %>% group_by(Disease) %>% dplyr::summarise( count = n())
  
  tmp = data.frame( count = dt$count, group = dt$Disease)
    
  tmp$ave = sum(tmp$count)/length(tmp$count)
  
  tmp$proportion = tmp$count/tmp$ave
  
  g <- ggdotchart(tmp, x = "group", y = "proportion",
                  color = "group",                                      
                  #palette =  Choose_col[c(1:15)],
                  sorting = "descending",                                
                  add = "segments",                                      
                  add.params = list(color = "lightgray", size = 2),     
                  dot.size = 8,                                           
                  ggtheme = theme_pubr())+ 
                  geom_hline(yintercept = 1, linetype = 2, color = "gray") +          
                  ylab("F-apCAF")
    g
    
    ggsave("./fig/Figure2L.pdf", g , width = 3, height = 6) #Figure2K

}

Figure3C{
  library(SpatialDecon)
  library(GeomxTools)
  
  theme_cxf <- theme(panel.background = element_rect(fill = "transparent", color = "black",size = 1), 
                     legend.key = element_rect(fill = "transparent", color = "transparent"),
                     text = element_text(color = "black", size = 6,vjust=0.5),
                     plot.title = element_text(hjust = .5), 
                     axis.text = element_text(color = "black",size = 6))
  
  
  var_genes<- epi_sub@assays$RNA@var.features
  epi_sub<- epi_sub[var_genes,] 
  annots <- data.frame(cbind(cellType=as.character(Idents(epi_sub)), 
                             cellID=names(Idents(epi_sub))))
  
  custom_mtx_seurat <- create_profile_matrix(mtx = EPI_sub3@assays$RNA@data, 
                                             cellAnnots = annots, 
                                             cellTypeCol = "cellType", 
                                             cellNameCol = "cellID", 
                                             matrixName = "custom_mini_colon",
                                             outDir = NULL, 
                                             normalize = FALSE, 
                                             minCellNum = 50, 
                                             minGenes = 50)
  
  Norm = read.csv("Spacial_data.csv")
  rownames(Norm)<-Norm$X
  Norm<-Norm[,2:47]
  dim(Norm)
  Norm<-as.matrix(Norm)
  class(Norm)
  
  bg<-Norm 
  bg[,]<-1
  dim(bg)
  class(bg)

  res = spatialdecon(norm = Norm,
                     bg = bg,
                     X = custom_mtx_seurat,
                     align_genes = TRUE)
  
  str(res)
  
  Res_prop_of_all<-res$prop_of_all
  
  library(reshape2)
  pro<- t(Res_prop_of_all)
  pp<-melt(pro, measure.vars = c("iCMS2", "iCMS3"), variable.name = "Group")   
  colnames(pp)<-c("ROI","Cell_type","value")
  
  
  p <- ggplot(data = pp, aes(x = ROI, y = value, fill = Cell_type)) +
       geom_bar(stat = "identity", width=0.9,position="fill")+
       scale_fill_manual(values=col_map[c(4,3)])+
       coord_flip()+
       theme_cxf 
  
  ggsave("./fig/Figure3C.pdf",p, height = 6, width =6)

}

Figure3H{
  library(limma)
  data.rna = read.csv("data_caf.csv") 
  rownames(data.rna)<-data.rna$gene
  data.rna<-data.rna[,2:46]
 
  data.rna<-as.matrix(data.rna)
  class(data.rna)
  
  data.rna <- log2(data.rna) 
  
  design <- model.matrix(~ -1+factor(c(rep(1,23),rep(2,22))))
  colnames(design) <- c("S1", "S2") 
  contrast.matrix <- makeContrasts(S1-S2, levels=design)
  fit <- lmFit(data.rna, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  full_expr <-topTable(fit2, coef=1, adjust="fdr", sort.by="B", number=nrow(data.rna))
  deg_expr <-subset(full_expr,abs(full_expr$adj.P.Val) <= 0.05)
  deg_expr <- deg_expr[order(deg_expr$logFC, decreasing = TRUE),]
 
  results <- decideTests(fit2)

  library(ggpubr)
  if(T){
    df=full_expr
    attach(df)
    df$`-log10(P.Value)` <- -log10(P.Value)
    df$group=ifelse(df$P.Value>0.05,'stable',
                    ifelse( df$logFC >0.75,'up', 
                            ifelse(df$logFC < -0.75,'down','stable')))
    df$name=rownames(df)
    
    library(ggpubr)
    p<-ggscatter(df, x = "logFC", y = "-log10(P.Value)", color = "group",size = 2,
                 label = "name", repel = T,
                 label.select = rownames(df)[df$g != 'stable'] ,
                 palette = c("#2f4e7c", "gray", "#932626") )
    p+
      geom_hline(aes(yintercept=-log(0.05,10)),linetype="dashed")+
      geom_vline(aes(xintercept=-0.75),linetype="dashed")+
      geom_vline(aes(xintercept=0.75),linetype="dashed")+
      theme_cxf+
      theme(legend.position = 'none')
      ggsave('volcano_caf.pdf')
}

}

Figure3I{

  enrichment_data <- read.csv("ICMS_CAF_Human_table.csv")
 
  enrichment_data$Gene.Count <- as.numeric(str_extract(enrichment_data$Overlap, "^[0-9]+"))
  enrichment_data$Pathway.Gene.Count <- as.numeric(str_extract(enrichment_data$Overlap, "[0-9]+$"))
  enrichment_data$Rich.Factor <- enrichment_data$Gene.Count / enrichment_data$Pathway.Gene.Count
  enrichment_data$logP <- -log10(enrichment_data$P.value)
  
  ggplot(enrichment_data, aes(x = reorder(Term, Gene.Count), y = Gene.Count, fill = Group)) +
        geom_bar(stat = "identity") +
        coord_flip() +  
        scale_fill_manual(values = c("iCMS3_CAF" = "red3", "iCMS2_CAF" = "blue3") ) + 
        labs(
          title = "Enriched pathway analysis",
          x = "",
          y = "Gene count",
          fill = "CAF type"
        ) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 12))
  
}


library(rlang)
library(Seurat)
library(SeuratData) 
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(scDblFinder)
library(cluster)
library(tidyverse)
library(here)
library(ggExtra)
library(ggrepel)
library(gghalves)
library(viridis)
library(cowplot)
library(Matrix)
library(Seurat)
library(SeuratData) 
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(ggExtra)
library(ggrepel)
library(viridis)
library(cowplot)

Figure4A-C{
  scrna <- readRDS("scrna.rds")
  levels(scrna$cell_type)<-c("CAF", "Myeloid", "ICMS2", "Lymphoid", "Endothelial", "ICMS3", "SMC", "Endothelial", "Myeloid")
  scrna$cell_type<-fct_relevel(scrna$cell_type,c("CAF","SMC","Endothelial","Lymphoid", "Myeloid", "ICMS3", "ICMS2"))
  DimPlot(scRNA_cluster_subset, reduction = "umap",raster=T,pt.size = 0.75,alpha =0.5,
          cols = c("#FFC125", "#8B2252", "#B4EEB4","#9ACD32","#FF7F50","#EED5D2", "#FFD39B")) 
          +theme_void()
  dev.off()

  
  genes<-c(
    "DPEP1", "EPCAM","FABP1","KRT18","KRT19",  
    "ITGB8","SPINK4","REG4","MUC2","SERPINA1",   
    "FCER1G","AIF1","CD14", "CD68", "FCGR3A",    
    "PTPRC","CD3D","CD8A","CD19","JCHAIN",      
    "PECAM1","CDH5","VWF","FABP4","CD34",        
    "MYH11","DES","TAGLN","ACTG2","PLN",     
    "COL1A2","COL5A1","LUM","CDH11","FAP"     
  )           
  
  g1 <- scRNAtoolVis::AverageHeatmap(scrna, genes,
                                    htCol= c("#0099CC", "white", "#DB4263"),
                                    myanCol = cols,
                                    group.by = "cell_type", 
                                    assays = "SCT") 
  
  ggsave('./fig/Figure4B.pdf', g1, height = 3, width = 6)
  

  g2 <- scRNAtoolVis::cellRatioPlot(object = scrna,sample.name="group",
                                    celltype.name = 'cell_type',
                                    fill.col = cols, col.width=0.5,
                                    flow.alpha=0.2)
  
  ggsave('./fig/Figure4C.pdf', g2, height = 6, width = 6)

}

Figure4D-E{
  library(spacexr)
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyverse)
  library(doParallel)
  
  xenium.obj <- readRDS("./data/xenium.RDS")
  query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")
  coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
 
  crpm.ref <- readRDS("./data/crpm.ref.RDS")
  counts <- GetAssayData(crpm.ref, assay = "RNA", slot = "counts")
  cluster <- as.factor(crpm.ref$clusters)
  names(cluster) <- colnames(crpm.ref)
  nUMI <- crpm.ref$nCount_RNA
  names(nUMI) <- colnames(crpm.ref)
  nUMI <- colSums(counts)
  reference <- Reference(counts, cluster, nUMI)
  
  RCTD <- create.RCTD(query, reference, max_cores = 8,
                      UMI_min = 0,
                      UMI_max = 2e+07)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")  
  
  annotations.df <- RCTD@results$results_df
  annotations <- annotations.df$first_type
  names(annotations) <- rownames(annotations.df)
  xenium.obj$predicted.celltype <- annotations
  keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
  xenium.obj <- subset(xenium.obj, cells = keep.cells)
  
  DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
  
  celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype",
                                size = 2, cols = "polychrome",
                                dark.background = F) + ggtitle("Cell type")+
                    scale_fill_manual(values = c('hotpink',"#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B","#FFAB91",
                                 "#FFCC80","#FFE082","#FFF59D","#AED581","#A5D6A7","#80CBC4","#26C6DA"))
  
  niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 2, 
                             dark.background = F) + ggtitle("Niches") +
                scale_fill_manual(values = c("#442288","#B5D33D",  "#EB7D5B","#FED23F", "#6CA2EA" ))
  
  celltype.plot | niche.plot
  
  dev.off()
  
  
  data <- xenium.obj@meta.data %>% select("predicted.celltype","niches")
  data$niches <-as.factor(data$niches)
  data.pie <- as.data.frame(table(data))
  data.sum<-data.pie %>% group_by(niches) %>% summarise(all = sum(Freq, na.rm = TRUE))
  data.pie$all<-rep(data.sum$all,each=14) 
  data.pie<-data.pie %>% mutate(ratio = round(100*Freq/all,digits = 2))
  for (i in seq(4)){
    i=4
    data.subset<- subset(data.pie,data.pie$niches %in% i)
  
    data <- data.frame(
            category=data.subset$predicted.celltype,
            fraction =data.subset$ratio
            )
    data$ymax <- cumsum(data$fraction)
 
    data$ymin <- c(0, head(data$ymax, n=-1))

    data$labelPosition <- (data$ymax + data$ymin) / 2

    data$label <- paste0(data$fraction,"%") 
    
    pdf(paste0("niches_proportion_",i,".pdf"), width=5,height=5)
  
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=0, fill=category)) +
          geom_rect() +
          geom_text( x=3, aes(y=labelPosition, label=label), size=6) + # x here controls label position (inner / outer)
          scale_fill_manual(values = c('hotpink',"#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B","#FFAB91",
                                       "#FFCC80","#FFE082","#FFF59D","#AED581","#A5D6A7","#80CBC4","#26C6DA"))+
          coord_polar(theta="y") +
          xlim(c(-1, 4)) +
          theme_void() +
          theme(legend.position = "none")
    
    dev.off()
    
}

}

Figure4G{
  library(dplyr)
  data<-read.csv("TLS_cluster_proportion.csv")
  data$Cluster<-as.factor(data$Cluster)
  data$Cluster <- fct_relevel(data$Cluster,c("F_apCAF" , "M_apCAF", "B cells" , "Cancer cells", "Myeloid",  
                                             "Endothelial", "ICAF" , "Mast" , "MyCAF" ,"Plasma" , "pericyte" , 
                                             "SMC" , "SsCAF" , "T cells" ))

  
  TOTAL<- data %>% group_by(TLS) %>% summarise(Total = sum(Count))

  TLS_order <- c('TLS_1', 'TLS_2', 'TLS_3', 'TLS_4', 'TLS_5', 
                 'TLS_6', 'TLS_7', 'TLS_8', 'TLS_9', 'TLS_10')
  
  TOTAL_sorted <- TOTAL %>%
                  mutate(TLS = factor(TLS, levels = TLS_order)) %>%
                  arrange(TLS)

  data$Total <- rep(TOTAL_sorted$Total,each=14)
  data$Proportion <- round(data$Count/data$Total,3)
  
  cluster_colors <- c('#49e200', "#deed00", "#a078ff", "gray", "#5bebff", 
                      "#b983bb", "#df3076",  "#6c9200", "#ef896e", 
                      "#bf9200", "#4afb88", "#d77bf3")
  
  p <- ggplot(data, aes(x = Cluster, y = Proportion)) +
    
    geom_boxplot(aes(group = Cluster), 
                 fill = "gray90", 
                 alpha = 0.5, 
                 outlier.shape = NA) +  

    geom_line(aes(color = TLS, group = TLS), 
              size = 0.5, 
              alpha = 0.5) +
    
    geom_point(aes(color = TLS), 
               size = 2, 
               alpha = 0.5) +

    geom_hline(yintercept = 0.071, 
               linetype = "dashed", 
               color = "black", 
               size = 0.5) +
    
    scale_color_manual(values = cluster_colors) +
    
    theme_minimal() +
    
    labs(title = "Distribution of Each Cluster Across Different TLSs",
         x = "TLS",
         y = "Proportion",
         color = "Cluster") +
    
    theme_cxf+

    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  

  print(p)
  
}

Figure4H{
  
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(gplots)
  library(stats)

  save_dir = "./fig/"
  folder = paste0(save_dir,"TLS_distance_plot")
  
  if (!file.exists(folder)){
    dir.create(folder)
  }
  
  set.seed(1)
 
  xenium <- readRDS("./data/xenium.RDS")
  meta<-as.data.frame(xenium_1@meta.data)
  cells = dplyr::sample_frac(meta,0.2)
  cells = rownames(cells)
  SOE <- subset(xenium_1, cells = cells)
  levels(SOE$predicted.celltype) <- c( "ApCAF_1", "ApCAF_2", "TLS", "Cancer cells",
                                       "Myeloid", "Endothelial", "ICAF",  "Mast",       
                                       "MyCAF", "TLS",  "pericyte", "SMC",
                                       "SsCAF", "TLS")
  
  loc<-as.data.frame(SOE@images$fov@boundaries$centroids@coords)
  rownames(loc)<-rownames(SOE@meta.data)
  loc$clusters<-SOE$predicted.celltype
  TLS_1_loc<-subset(loc,loc$clusters %in% c("TLS"))
  TLS_1_con_loc<-subset(loc,!(loc$clusters %in% c("TLS")))
  
  min_distance = function(tumor_loc,stroma_loc){
    mini_dist = c()
    for (i in seq(1,dim(stroma_loc)[1])){
      xF = stroma_loc[i,1]
      yF = stroma_loc[i,2]
      min_distance = 100000  
      for (j in seq(1,dim(tumor_loc)[1])){
        xT = tumor_loc[j,1]
        yT = tumor_loc[j,2]
        min_distance = min(min_distance, sqrt((xF-xT)^2+(yF-yT)^2))
      }
      mini_dist = c(mini_dist,min_distance)
    }
    return(mini_dist)
  }
  
  df<-data.frame(distance = min_distance(TLS_1_loc,TLS_1_con_loc),
                 cluster =  TLS_1_con_loc$clusters)
  
  
  mycolpalette1 = c('yellow','orange')
  

  df %>%
    filter(cluster %in% c("F-apCAF","MyCAF")) %>% 
    ggplot(aes(x=distance,group = cluster, fill = cluster))+
    scale_fill_manual(values=mycolpalette1)+
    facet_grid(cluster ~ ., scales = "free")+
    xlab("minimum distance from TLS (um)")+
    geom_density( geom = "area", position = "stack", adjust=2,linewidth=1.2) +
    theme_cxf+
    xlim(0, 200)
  

  mean_std_dist <- function(dist_all,SID){
    for (i in seq(1,length(dist_all))){
      tmp = 0.2*dist_all[[i]]
      m = mean(tmp)
      s = sd(tmp)
      print(paste0('sample',SID,', mean:',m,'std:',s))
    }
  }
  

}

Figure5A-D{
  library(spacexr)
  library(Seurat)
  library(SeuratData)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyverse)
  library(doParallel)
  
  xenium.obj <- readRDS("./data/xenium.RDS")
  query.counts <- GetAssayData(xenium.obj, assay = "Xenium", slot = "counts")
  coords <- GetTissueCoordinates(xenium.obj, which = "centroids")
  rownames(coords) <- coords$cell
  coords$cell <- NULL
  query <- SpatialRNA(coords, query.counts, colSums(query.counts))
  
  crpm.ref <- readRDS("./data/crpm.ref.RDS")
  counts <- GetAssayData(crpm.ref, assay = "RNA", slot = "counts")
  cluster <- as.factor(crpm.ref$clusters)
  names(cluster) <- colnames(crpm.ref)
  nUMI <- crpm.ref$nCount_RNA
  names(nUMI) <- colnames(crpm.ref)
  nUMI <- colSums(counts)
  reference <- Reference(counts, cluster, nUMI)
  
  RCTD <- create.RCTD(query, reference, max_cores = 8,
                      UMI_min = 0,
                      UMI_max = 2e+07)
  RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")  
  
  annotations.df <- RCTD@results$results_df
  annotations <- annotations.df$first_type
  names(annotations) <- rownames(annotations.df)
  xenium.obj$predicted.celltype <- annotations
  keep.cells <- Cells(xenium.obj)[!is.na(xenium.obj$predicted.celltype)]
  xenium.obj <- subset(xenium.obj, cells = keep.cells)
  
  DefaultBoundary(xenium.obj[["zoom"]]) <- "segmentation"
  
  celltype.plot <- ImageDimPlot(xenium.obj, group.by = "predicted.celltype",
                                size = 2, cols = "polychrome",
                                dark.background = F) + ggtitle("Cell type")+
                   scale_fill_manual(values = c('hotpink',"#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B","#FFAB91",
                                 "#FFCC80","#FFE082","#FFF59D","#AED581","#A5D6A7","#80CBC4","#26C6DA"))
  
  niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 2, 
                             dark.background = F) + ggtitle("Niches") +
    scale_fill_manual(values = c("#442288","#B5D33D",  "#EB7D5B","#FED23F" ))
  
  celltype.plot | niche.plot
  
  dev.off()
  
  
  data <- xenium.obj@meta.data %>% select("predicted.celltype","niches")
  data$niches <-as.factor(data$niches)
  data.pie <- as.data.frame(table(data))
  data.sum<-data.pie %>% group_by(niches) %>% summarise(all = sum(Freq, na.rm = TRUE))
  data.pie$all<-rep(data.sum$all,each=14) 
  data.pie<-data.pie %>% mutate(ratio = round(100*Freq/all,digits = 2))
  for (i in seq(4)){
    i=4
    data.subset<- subset(data.pie,data.pie$niches %in% i)
    
    data <- data.frame(
      category=data.subset$predicted.celltype,
      fraction =data.subset$ratio
    )
    data$ymax <- cumsum(data$fraction)
    
    data$ymin <- c(0, head(data$ymax, n=-1))
    
    data$labelPosition <- (data$ymax + data$ymin) / 2
    
    data$label <- paste0(data$fraction,"%") 
    
    pdf(paste0("niches_proportion_",i,".pdf"), width=5,height=5)
    
    ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=2, xmin=0, fill=category)) +
      geom_rect() +
      geom_text( x=3, aes(y=labelPosition, label=label), size=6) + # x here controls label position (inner / outer)
      scale_fill_manual(values = c('hotpink',"#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B","#FFAB91",
                                   "#FFCC80","#FFE082","#FFF59D","#AED581","#A5D6A7","#80CBC4","#26C6DA"))+
      coord_polar(theta="y") +
      xlim(c(-1, 4)) +
      theme_void() +
      theme(legend.position = "none")
    
    dev.off()
    
  }
  
  
  genes<-c( "TIGIT","CTLA4","LAG3","PDCD1","ICOS","FOXP3" ) 
  Idents(scrna) <- "group_niches"
  DotPlot(scrna, features = genes, cols = c("#e6e4df", "red3"), dot.scale = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme_cxf
  dev.off()  
  
  
}


library(rlang)
library(Seurat)
library(SeuratData) 
library(patchwork)
library(dplyr)
library(SingleCellExperiment)
library(scDblFinder)
library(cluster)
library(tidyverse)
library(here)
library(ggExtra)
library(ggrepel)
library(gghalves)
library(viridis)
library(cowplot)
library(Matrix)


Figure6H{
  
  seu <-readRDS("seu_spp1_ko.RDS")
  new.cluster.ids <- c("Epithelial cells","Aif Myeloid cells-1","Neutrophils","T and NK cells",
                       "Aif Myeloid cells-2","B and Plasma cells","CAFs", "Endothelial cells", 
                       "Aif Myeloid cells-3", "Mast cells")
  names(new.cluster.ids) <- levels(seu)
  seu <- RenameIdents(seu, new.cluster.ids)
  seu$cell_types <- Idents(seu)
  
  library(ggpubr)
  df <- table(seu@meta.data$orig.ident,seu@meta.data$cell_types) %>% melt()
  colnames(df) <- c("Sample","Cluster","Number")
  df$Group <- df$Sample
  levels(df$Group)<-c("KO" ,"KO", "KO", "WT", "WT", "WT")
  df %<>% group_by(Cluster)  %>% mutate(Total = sum(Number))
  df$Ratio <-(df$Number/df$Total)
 
  ggplot(df, aes(x=Cluster, y=Ratio, fill=Group)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values=rev(col_map[c(1:10)]))  +
    theme_classic() +
    stat_compare_means(method ="wilcox.test", paired = F,  method.args = list(alternative = "less")) +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 18))+
    theme_cxf+
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
  
   
}


Figure6I-L{
  CAF_subset <-subset(seu, idents =  c("CAFs"))
  CAF_subset <-SCTransform(CAF_subset,method = "glmGamPoi",vst.flavor = "v2")  %>% RunPCA(npcs = 50)
  CAF_subset <- RunHarmony(CAF_subset, group.by.vars = "orig.ident",assay.use = "SCT", block.size  =  0.1, max.iter.harmony = 5) 
  CAF_subset <- FindNeighbors(CAF_subset, reduction = "harmony", dims = 1:30) %>% 
                FindClusters(resolution = 0.25)
  CAF_subset <- RunUMAP(CAF_subset, reduction = "harmony", dims = 1:10)
  DimPlot(CAF_subset, reduction = "umap",pt.size = 1.2) + theme_cxf + scale_color_jco()

  library(colorspace)
  gene<-c("C3","Pdpn","Il6","Cxcl12","Ccl11",                             
          "Cd74","H2-Ab1","H2-Eb1","Ccnd1","Sorbs2",                         
          "Tagln","Acta2","Ndufa4l2","Actg2","Cthrc1")  
  DotPlot(CAF_subset, features = gene, cols = c("#e6e4df", "#FA9900"), 
          col.min = 0, col.max = 2,dot.scale = 8)  +
          theme_cxf+
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  library(reshape2)
  df <- table(CAF_subset@meta.data$Group, CAF_subset@meta.data$Celltypes) %>% melt()
  colnames(df) <- c("Group","Cluster","Number")
  df$Group <- fct_relevel(df$Group,c("WT","KO"))
  ggplot(data = df, aes(x = Group, y = Number, fill = Cluster)) +
    geom_bar(stat = "identity", width=0.8,position="fill")+
    scale_fill_nejm()+
    theme_bw()+
    theme(panel.grid =element_blank()) +
    labs(x="",y="Ratio")+
    ####用来将y轴移动位置
    theme(axis.text.y = element_text(size=12, colour = "black"))+
    theme(axis.text.x = element_text(size=12, colour = "black"))+
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
  
  
  library(ggpubr)
  df <- table(CAF_subset@meta.data$orig.ident,CAF_subset@meta.data$Celltypes) %>% melt()
  colnames(df) <- c("Sample","Cluster","Number")
  df$Group <- df$Sample
  levels(df$Group)<-c("KO" ,"KO", "KO", "WT", "WT", "WT")
  df %<>% group_by(Cluster)  %>% mutate(Total = sum(Number))
  df$Ratio <-(df$Number/df$Total)
  
  ggplot(df, aes(x=Cluster, y=Ratio, fill=Group)) +
    geom_boxplot( outliers =F,outlier.shape = 19, outlier.size = 0.1, varwidth = T) +
    scale_fill_manual(values=rev(col_map[c(1:10)]))  +
    theme_classic() +
    stat_compare_means(method ="wilcox.test", paired = F,  method.args = list(alternative = "less")) +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 18))+
    theme_cxf+
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))
  
  
  FeaturePlot(CAF_subset, features = "Pi16",pt.size = 1,col = c("#e7e1ef", "#b71540"),split.by = "Group") 
  FeaturePlot(CAF_subset, features = "Dpt",pt.size = 1,col = c("#e7e1ef", "#b71540"),split.by = "Group") 
  
}  

Figure6M{
  
  
  
}








fb <- readRDS("./data/d_Fb7.rds")
fb_meta <- fb@meta.data
fb_meta <- fb_meta[fb_meta$Disease %in% "PDAC",]
dim(fb_meta)
write.csv(fb_meta,"./data/fb_meta.csv")

mye <- readRDS("./data/d_Mye.rds")
mye_meta <- mye@meta.data
mye_meta <- mye_meta[mye_meta$Disease %in% "PDAC",]
dim(mye_meta)
write.csv(mye_meta,"./data/mye_meta.csv")

TNK <- readRDS("./data/d_T.rds")
TNK_meta <- TNK@meta.data
TNK_meta <- TNK_meta[TNK_meta$Disease %in% "PDAC",]
dim(TNK_meta)
write.csv(TNK_meta,"./data/TNK_meta.csv")






sample_count{
  sn <- meta %>% filter(Group == "Normal")  
  sn <-sn[!duplicated(sn$SampleID),]  
  sn_1<-as.data.frame(table(sn$Tissue))
  colnames(sn_1)<-c("Tissue","Normal")
  st <- meta %>% filter(Group == "Tumor")  
  st <-st[!duplicated(st$SampleID),]  
  st_1<-as.data.frame(table(st$Tissue))
  colnames(st_1)<-c("Tissue","Tumor")
  s<-merge(sn_1,st_1,by = "Tissue")
  s <- reshape2::melt(s, di.var = 'Tissue')
  colnames(s)<-c("Tissue","Group","Count")
  sum(s$Count)
  ggplot(s, aes(x=Tissue, y=Count, fill=Tissue))+
    geom_bar(stat = "identity",
             width=0.6,
             position = position_dodge(width=0.6))+
    facet_grid(~Group)+
    theme_cxf+
    scale_fill_manual(values = Choose_col)+
    theme(axis.text.x=element_text(angle=45, hjust=1,size=14) )
}


Cell_nk_distribution{
  seurat_A <-readRDS("./data/d_seurat_A.rds")
  meta<-seurat_Fb@meta.data
  #write.csv(meta,"./data/seurat_Fb_meta.csv")
  table(seurat_Fb$Tissue)
  table(seurat_Fb$Disease)
  seurat_Fb$Tissue<-as.factor(seurat_Fb$Tissue)
  levels(seurat_Fb$Tissue)<-c("Bladder","Breast","Cervix","Colorectum","Colorectum","Esophagus","Kidney","Liver","Lung",   
                              "Omentum","Ovarian","Pancreas","Prostate","Skin","Stomach")
  seurat_Fb$Disease<-as.factor(seurat_Fb$Disease)
  levels(seurat_Fb$Disease)<-c("BCC","BLCA","BRCA","ccRCC","CESC","COLO","ESCC","HCC",NA,"NSCLC","OV","PDAC","PMP",
                               "PRAD","SCC","STAD")
  meta <- read.csv("./data/seurat_Fb_meta.csv",row.names = 1)

  

  
  # his_cell_num = meta %>% filter(meta_tissue == "Blood") 
  # his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
  # his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
  # his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")
 his_cell_num = meta 
 colnames(meta)
 his_cell_num <- his_cell_num[,c("Group","seurat_clusters")]
 DT = his_cell_num  %>% group_by(Group,seurat_clusters)
 DT <- dplyr::summarise(DT,count = n())
 DT<-data.frame(DT)
  
  
  test <- reshape2::dcast(DT,Group~seurat_clusters)
  test[is.na(test)] <- 0
  rownames(test)<-test$Group
  test<-test[,-1]
  test <- as.matrix(test)
  
  a =margin.table(test, 1)
  b=margin.table(test, 2)
  c=margin.table(test)
  test = test/(outer(a,b,"*")/c)
  
  tmp <- reshape2::melt(test)
  tmp <- tmp[order(tmp$Var1,tmp$value,tmp$Var2),]
  
  tmp$Var1 <- factor(tmp$Var1,levels=as.character(unique(tmp$Var1)))
  tmp$state <- ifelse(tmp$value>1,"Enrichment","Depletion")
  
  # col = c("#F4511E","#1976D2")
  # names(col) = c("Enrichment","Depletion")
  
  g <- ggplot(data = tmp)+ 
    geom_point(aes(x = Var1,y = Var2, color = state, size = value))+
    scale_size_area()+
    #scale_color_gradient(low = "grey",high = "red") + 
    #scale_color_manual(values = col) + 
    scale_color_manual(values = c("blue2","red")) + 
    cowplot::theme_cowplot() +
    # theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
    #       axis.text.x = element_text(size = 12),
    #       strip.background = element_rect(colour = "white")) + 
    theme_cxf+
    ylab("Cell Type") +
    xlab("") 
  g
  
  
  
  
 his_cell_num = meta 
 his_cell_num <- his_cell_num[,c("Disease","seurat_clusters")]
 head(his_cell_num)
 table(his_cell_num$Disease)
 DT = his_cell_num  %>% group_by(Disease,seurat_clusters)
 DT <- dplyr::summarise(DT,count = n())
 DT<-data.frame(DT)
  
  test <- reshape2::dcast(DT,Disease~seurat_clusters)
  test <- test[1:15,]
  test[is.na(test)] <- 0
  rownames(test)<-test$Disease
  test<-test[,-1]
  test <- as.matrix(test)
  
  a =margin.table(test, 1)
  b=margin.table(test, 2)
  c=margin.table(test)
  test = test/(outer(a,b,"*")/c)
  
  tmp <- reshape2::melt(test)
  tmp <- tmp[order(tmp$Var1,tmp$value,tmp$Var2),]
  
  tmp$Var1 <- factor(tmp$Var1,levels=as.character(unique(tmp$Var1)))
  tmp$state <- ifelse(tmp$value>1,"Enrichment","Depletion")
  
  tmp <- tmp[tmp$Var2=="c05",]
  colnames(tmp) <- c("group","celltype","value","state")
  
  library("ggpubr")
  g<-ggdotchart(tmp, x = "group", y = "value",
                color = "group",                                        # Color by groups
                #palette =  Choose_col[c(1:15)], # Custom color palette
                sorting = "descending",                                 # Sort value in descending order
                add = "segments",                                       # Add segments from y = 0 to dots
                add.params = list(color = "lightgray", size = 2),       # Change segment color and size
                dot.size = 8,                                           # Large dot size
                # font.label = list(color = "white", size = 3,
                #                   vjust = 0.5),                       # Adjust label parameters
                # rotate = TRUE,
                ggtheme = theme_pubr()                                  # ggplot2 theme
  )+ geom_hline(yintercept = 1, linetype = 2, color = "gray") +          #+ geom_hline(yintercept = -1, linetype = 2, color = "gray") 
  ylab("MKI67+CAF")
  g
  
  #ggsave(paste0("./figS7/","c8-KLRC2-ROE.pdf"),g, width = 4, height = 5)
}

Cell_nk_Hierarchical_clustering{
  library(tidyverse)
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(ggsci)
  library(cowplot)
  
  meta <- read.csv("./data/seurat_Fb_meta.csv",row.names = 1)
  
  #Hierarchical clustering for cancer types
  his_cell_num = meta %>% filter((Group == "Tumor"))
  
  # his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
  # his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
  # his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")
  
  
  #cancernum = his_cell_num %>%  group_by(cancerType2, celltype) 
  levels(his_cell_num$Disease) 
  test = his_cell_num %>% group_by(Disease) 
  
 DT = meta %>%  group_by(seurat_clusters,Disease) 
 DT <- summarise(cancernum,count = n())
 DT<-data.frame(cancernum)
  
  test <- data.frame(summarise(test,
                               count = n()))
  
 DT <- merge(test,cancernum,by="Disease")
 DT['proportion'] <-DT$count.y /DT$count.x
  head(cancernum)
  
  library(reshape2)
 DT2 <-DT[,c("Disease","seurat_clusters","proportion")]
 DT2 <- dcast(cancernum2,Disease~seurat_clusters)
 DT2 <-cancernum2[1:15,]
  rownames(cancernum2) <-DT2$Disease
 DT2 <-DT2[,-1]
 DT2[is.na(cancernum2)] <- 0
  
  tree = hclust(vegan::vegdist(cancernum2, method = 'bray'), method = 'average')
  
  p1 = ggtree(tree) + geom_tiplab() +xlim(NA,2)
  p1
  
 DT3 <-DT[,c("Disease","seurat_clusters","proportion")]
 DT3 <- na.omit(cancernum3)
  levels(cancernum3$Disease)
  Choose_col=
  c("#FFCCBC", "#FFE0B2", "#FFECB3", "#FFF59D", "#DCE775", "#DCEDC8", "#A5D6A7", "#B2DFDB", "#80DEEA", "#4FC3F7", "#64B5F6", "#C5CAE9",
    "#9575CD", "#E1BEE7", "#F48FB1", "#EF9A9A", "#D7CCC8", "#F5F5F5", "#B0BEC5")
  g2 <-DT3 %>% ggplot(aes(y=proportion,x=Disease,fill=seurat_clusters,colour = seurat_clusters))+
    geom_bar(stat = "identity",width = 0.5,alpha=0.65) +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_blank())+
    scale_x_discrete(limits = rev(rev(c("ccRCC","PRAD", "CESC","PDAC","BLCA","HCC","PMP","OV","ESCC","BCC","BRCA","NSCLC","SCC","COLO","STAD"))),) +
    #scale_fill_manual(values = c(CD16low_col, CD16high_col)) +
    scale_colour_manual(values = Choose_col[c(1:15)]) +
    
    coord_flip()+
    ylab("Frequency") +
    xlab("") 
  
  g2
  
  g <- ggdraw()+
    draw_plot(p1, 0, 0.06, 0.6, 0.94)+
    draw_plot(g2, 0.2, 0, 0.8, 1)
  g
  
  ggsave("./result/Cell_nk_Hierarchical_clustering/Fb_clusters.pdf", g, height = 5, width =10)
}

dev.off()
Cell_nk_Plot{
  saveRDS(seurat_A,"./data/d_seurat_A.rds")
  getwd()
  setwd("D:/data_subset5/subset/seurat_all/")
  seurat_A <-readRDS("./data/d_seurat_A.rds")
  meta<-seurat_A@meta.data
  FeaturePlot(seurat_A, features = c("PTN"), split.by = "Tissue",
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("PTN"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("FGFR1"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("FGF2"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("CD80"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("CD40"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  Idents(seurat_A)<-"Disease"
  DimPlot(seurat_A)
  
  
  dev.off()
  FeaturePlot(seurat_A, features = c("IL4I1"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))
  FeaturePlot(seurat_A, features = c("MAOB"),
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)),
              split.by = "Tissue")
  dev.off()
  
  FeaturePlot(seurat_A, features = c("ANO6","PLACR1",
                                     "XKR8","PLACR2",
                                     "TMEM41B","CHO1",
                                     "ATP8A2","ATP9A"), 
              cols = c(rgb(225,225,225,10, maxColorValue = 255),
                       rgb(250,27,40,100,  maxColorValue = 255)))

  
  
  {
  #write.csv(meta,"./data/seurat_A_meta.csv")
  table(seurat_A$Tissue)
  table(seurat_A$Disease)
  seurat_A$Tissue<-as.factor(seurat_A$Tissue)
  levels(seurat_A$Tissue)<-c("Bladder","Breast","Cervix","Colorectum","Colorectum","Esophagus","Kidney","Liver","Lung",   
                              "Omentum","Ovarian","Pancreas","Prostate","Skin","Stomach")
  seurat_A$Disease<-as.factor(seurat_A$Disease)
  levels(seurat_A$Disease)<-c("BCC","BLCA","BRCA","ccRCC","CESC","COLO","ESCC","HCC",NA,"NSCLC","OV","PDAC","PMP",
                               "PRAD","SCC","STAD")
  table(seurat_A$seurat_clusters) #0-49
  levels(seurat_A$seurat_clusters)<-c(
    'T1',
    'F1',
    'C1',
    'D1',
    'E1',
    'C2',
    'C3',
    'C4',
    'P1',
    'F2',
    'D2',
    'F3',
    'C5',
    'B1',
    'C6',
    'C7',
    'C8',
    'T2',
    'E2',
    'T3',
    'M1',
    'F4',
    'P2',
    'T4',
    'D3',
    'C9',
    'T5',
    'C10',
    'F5',
    'D4',
    'C11',
    'F6',
    'C12',
    'P3',
    'T6',
    'D5',
    'T7',
    'D6',
    'P4',
    'P5',
    'T8',
    'E3',
    'C13',
    'C14',
    'C15',
    'P6',
    'E4',
    'C16',
    'P7',
    'C17'
  )
  
  seurat_A$seurat_clusters<-fct_relevel(seurat_A$seurat_clusters,
                                                    c('B1',
                                                      'C1',
                                                      'C2',
                                                      'C3',
                                                      'C4',
                                                      'C5',
                                                      'C6',
                                                      'C7',
                                                      'C8',
                                                      'C9',
                                                      'C10',
                                                      'C11',
                                                      'C12',
                                                      'C13',
                                                      'C14',
                                                      'C15',
                                                      'C16',
                                                      'C17',
                                                      'D1',
                                                      'D2',
                                                      'D3',
                                                      'D4',
                                                      'D5',
                                                      'D6',
                                                      'E1',
                                                      'E2',
                                                      'E3',
                                                      'E4',
                                                      'F1',
                                                      'F2',
                                                      'F3',
                                                      'F4',
                                                      'F5',
                                                      'F6',
                                                      'M1',
                                                      'P1',
                                                      'P2',
                                                      'P3',
                                                      'P4',
                                                      'P5',
                                                      'P6',
                                                      'P7',
                                                      'T1',
                                                      'T2',
                                                      'T3',
                                                      'T4',
                                                      'T5',
                                                      'T6',
                                                      'T7',
                                                      'T8'
                                                    ))
  
  Idents(seurat_A) <- "seurat_clusters"
  DimPlot(seurat_A,label = T)
  seurat_A$Type<-  seurat_A$seurat_clusters
  levels(seurat_A$Type)<-c("B cell", "Epithelial cell", "Epithelial cell" , "Epithelial cell" , "Epithelial cell" , "Epithelial cell" , "Epithelial cell" , 
                                       "Epithelial cell" , "Epithelial cell" , "Epithelial cell" , "Epithelial cell", "Epithelial cell", "Epithelial cell", "Epithelial cell",
                                       "Epithelial cell", "Epithelial cell", "Epithelial cell", "Epithelial cell", "Myeloid cell" , "Myeloid cell" , "Myeloid cell" , "Myeloid cell" , 
                                       "Myeloid cell" , "Myeloid cell" , "Endothelial cell" , "Endothelial cell" , "Endothelial cell" , "Endothelial cell" , "Fibroblasts cell" , 
                                       "Fibroblasts cell" , "Fibroblasts cell" , "Fibroblasts cell" , "Fibroblasts cell" , "Fibroblasts cell" ,
                                       "Mast cell" , "Plasma cell" , "Plasma cell" , "Plasma cell" , "Plasma cell" , "Plasma cell" ,  "Plasma cell" , "Plasma cell" , 
                                       "T cell" , "T cell" , "T cell" , "T cell" , "T cell" , "T cell" , "T cell" , "T cell")
  }
    

split_violin{
  library(devtools)
  usethis::use_git_config(user.name = "windyeeeee", user.email = "chenxiongfengfj@qq.com")
  credentials::set_github_pat()
  install_github("JanCoUnchained/ggunchained")
  library(ggunchained)
  library(ggpubr)
  rm(seurat_A)
   
  seurat_Fb <- readRDS("./data/d_Fb.rds")
  meta <- seurat_Fb@meta.data
  his_cell_num <- meta[,c("SampleID","Disease","Tissue","Group","seurat_clusters")]
  his_cell_num$TD<-paste0(his_cell_num$Tissue,"_",his_cell_num$Disease)
  table(his_cell_num$TD)
  his_cell_num$TD<-as.factor(his_cell_num$TD)
  levels(his_cell_num$TD) <- c(
    "Bladder_BLCA","Bladder_Normal","Breast_BRCA","Breast_Normal","Cervix_CESC","Cervix_Normal",   
    "Colorectum_COLO","Colorectum_Normal","Esophagus_ESCC","Esophagus_Normal","Kidney_ccRCC","Kidney_Normal",    
    "Liver_HCC","Liver_Normal","Lung_Normal","Lung_NSCLC","Omentum_Normal","Omentum_PMP",   
    "Ovarian_Normal", "Ovarian_OV","Pancreas_Normal","Pancreas_PDAC","Prostate_Normal","Prostate_PRAD",  
    "Skin_BCC/SCC","Skin_Normal","Skin_BCC/SCC","Stomach_Normal","Stomach_STAD" )
  his_cell_num$TD<-fct_relevel(his_cell_num$TD,c(
    "Bladder_Normal","Bladder_BLCA","Breast_Normal","Breast_BRCA","Cervix_Normal", "Cervix_CESC",  
    "Colorectum_Normal","Colorectum_COLO","Esophagus_Normal","Esophagus_ESCC","Kidney_Normal","Kidney_ccRCC",    
    "Liver_Normal","Liver_HCC","Lung_Normal","Lung_NSCLC","Omentum_Normal","Omentum_PMP",   
    "Ovarian_Normal", "Ovarian_OV","Pancreas_Normal","Pancreas_PDAC","Prostate_Normal","Prostate_PRAD",  
    "Skin_Normal","Skin_BCC/SCC","Stomach_Normal","Stomach_STAD"       
                                                    ))

 DT = his_cell_num %>% group_by(TD,SampleID)      #filter(Group == "Tumor") %>%         
  test = his_cell_num %>% group_by(TD,SampleID,seurat_clusters)      #%>% filter(Group == "Tumor")  
 DT <- dplyr::summarise(cancernum,count = n())
 DT<-data.frame(cancernum)
  test <- data.frame(dplyr::summarise(test,count = n()))

 DT2 <- merge(test,cancernum,by="SampleID")
 DT2<- subset(cancernum2,DT2$count.y >  10)
 DT2['proportion'] <-DT2$count.x/cancernum2$count.y

 DT3 <- merge(his_cell_num[,c("SampleID","Tissue","TD","Group")],DT2, by="SampleID")
 DT3 <- na.omit(cancernum3)
  head(cancernum3)
  levels(cancernum3$seurat_clusters)
 DT4<- subset(cancernum3,DT3$seurat_clusters == "Pericytes")
  # col = c("#1976D2","#F4511E") 
  levels(cancernum4$TD)

  ggplot(cancernum4, aes(x = Tissue, y = proportion, fill = TD))+
         geom_split_violin(colour="black", scale = 'width', adjust= 3,lwd = 1)+
         scale_fill_manual(values = COL_A)+
         stat_summary(fun = mean,
                      fun.min = function(x){quantile(x)[2]},
                      fun.max = function(x){quantile(x)[4]},
                      geom = "pointrange",
                      #geom = 'errorbar',
                      size=0.5,
                      position = position_dodge(width = 0.2))+
                      ylim(0,1)+
                     stat_compare_means(
                       #aes(x = Tissue, y = proportion),
                      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                      symbols = c("***", "**", "*", "-")),
                      label = "p.signif",
                      #label.y = max(cancernum4$proportion),
                      hide.ns = F)+
         theme_cxf+
                     xlab("")+
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
                     legend.position = "top",
                     #legend.key = element_rect(fill = c("#1ba7b3","#dfb424")),
                      legend.justification = "right")
 
          #+facet_wrap(~meta_tissue.x, ncol = 3,scales = "free_y")
}
 
geom_violin{
  library(ggbeeswarm)
  g2 <-DT4 %>% filter(count.x > 50) %>% ggplot(aes(x= Tissue, y= proportion))+  #, colour = Disease
        geom_violin(aes(fill = TD), width=0.8, size=0.5, scale="width")+
        geom_quasirandom(size = 0.1,width = 0.15) +
        scale_fill_manual(values = COL_A)+
        stat_compare_means() +
        theme_cxf +
        theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
              axis.text.x = element_text(angle = 60,size = 8 , hjust = 1, vjust = 1),
              strip.background = element_rect(colour = "white"),
              axis.ticks =element_line(size=0.5),
              # aspect.ratio=1:2
              ) + 
        #scale_fill_manual(values = col) +
        #scale_colour_manual(values = col) +
        # ylab("Frequency of CD56dimCD16hi NK cells ") +
    ylab("Frequency of CD56highCD16low NK cells ") +
    xlab("") +
    facet_wrap(~TD.x, ncol = 3,scales = "free_y")
  g2
  ggsave("CD56highCD16low_proportion_violin_0311.pdf", g2, height = 4, width =7)
  ggsave("CD56lowCD16high_proportion_violin_0311.pdf", g2, height = 3.5, width =12)
}  
  

geom_box{
his_cell_num = meta %>% filter((meta_tissue == "Tumor")|(meta_tissue == "Normal")) #|(meta_tissue == "Blood")

his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")

his_cell_num = his_cell_num[his_cell_num$cancerType2 %in% c("LC","RC","THCA","ESCA","BRCA","GC"),]

his_cell_num $celltype <- ifelse(as.vector(as.matrix(lapply(his_cell_num $celltype, function(x){grepl(x,pattern="CD16hi")}))),"CD16high","CD16low" )
cancernum = his_cell_num %>%  group_by(batch,meta_tissue,meta_histology, celltype) #cancerType2,
test = his_cell_num %>%group_by(batch,meta_tissue,meta_histology) #cancerType2,

cancernum <- dplyr::summarise(cancernum,count = n())
cancernum<-data.frame(cancernum)

test <- data.frame(dplyr::summarise(test,
                                    count = n()))

cancernum <- merge(test,cancernum,by="batch")
cancernum <-DT[cancernum$meta_tissue.x ==DT$meta_tissue.y,]
cancernum['proportion'] <-DT$count.y /DT$count.x
head(cancernum)

col <- c("Normal"="#4974A4", "Tumor"="#B81316")
g <- 
 DT %>% ggplot(aes(y=proportion,x=celltype,colour = meta_tissue.x))+
  
  # geom_violin(aes(fill=meta_tissue.x), width=1,outlier.colour = NA, size=0.8,scale = "width")+
  geom_boxplot( fill="white", width=0.5,outlier.colour = NA, size=0.5, position = position_dodge(1))+
  geom_quasirandom(size = 0.5,width = 0.05, dodge.width = 1) +
  stat_compare_means(aes(group = meta_tissue.x), label = "p.format") +
  theme_classic() +
  theme(axis.title.y=element_text(size=12),
        axis.text.x = element_text(size = 12),
        legend.title = element_blank(),legend.position="none",
        strip.background = element_rect(colour = "white",fill = "white"),
        # legend.position="none",
        strip.text = element_text(size = 12, face = "bold" )) +  #,margins=TRUE)+ #axis.title.x =element_blank()
  scale_fill_manual(values = col) +
  scale_colour_manual(values = col) +
  ylab("Proportion in all NK cells") +
  xlab("") + facet_wrap(~meta_histology.x,nrow=2,scales = "free")
# print(g)
ggsave("proportion_box_majortype_NT_0311.pdf", g, height = 5, width =7)
}

}


Cell_nk_AuCell{
  library(AUCell)
  ### Calculation of signature score
  ECM_genes<-c("COL1A1","COL1A2","COL3A1","COL5A1","COL8A1","FN1","SPARC","MMP2","MMP11","MMP14","MMP19")
  Inflamatory_genes<-c("CXCL1","CXCL2","CXCL3","CXCL8","IL6","NFKB1")
  apCAF_genes<-c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74")
  All_genes<-c(ECM_genes,Inflamatory_genes,apCAF_genes)
  #matrix
  
  seurat_Fb<-readRDS("../data/d_Fb.rds")
  meta<-as.data.frame(seurat_Fb@meta.data)
  df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
  samp<-df %>% group_by(cluster) %>% sample_frac(size=0.2)  
  scRNAsub<-seurat_Fb[,samp$barcode]
  exprMatrix<-scRNAsub@assays$RNA@data
  meta<-as.data.frame(scRNAsub@meta.data)
  dim(scRNAsub)
  #GeneSets
  Genesets<-list(
    "ECM scores" = c("COL1A1","COL1A2","COL3A1","COL5A1","COL8A1","FN1","SPARC","MMP2","MMP11","MMP14","MMP19"),
    "Inflamatory scores" = c("CXCL1","CXCL2","CXCL3","CXCL8","IL6","NFKB1"),
    "apCAF scores" = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74")
  )


  ##AUC_score
  exprMatrix <- as(exprMatrix, "dgCMatrix")
  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
  cells_AUC <- AUCell_calcAUC(Genesets, cells_rankings)
  
  aucs <- getAUC(cells_AUC)
  
  aucs_t <- data.frame(t(aucs))  
  aucs_t$celltype <- meta[rownames(aucs_t), "seurat_clusters"]
  aucs_t <- aucs_t %>% group_by(celltype) %>% dplyr::summarise_each(funs = mean)
  getwd()
  write.csv(aucs_t ,"../result/Cell_nk_AUC/AUC.csv")
  
  ####plot
  test1 <- aucs_t
  test1$celltype <- factor(test1$celltype)
  table(test1$celltype)
  # test1$Mc <- ifelse(grepl("CD56bright", test1$celltype),"CD56brightCD16lo","CD56dimCD16hi")
  # col = c('#003399', '#999999')
  # names(col) = c('CD56brightCD16lo', 'CD56dimCD16hi')

  
  g1 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= ECM.scores),stat="identity",width = 0.8) + #,size = 3
    geom_hline(aes(yintercept=median(ECM.scores)),linetype="dashed") + 
    ylab("ECM scores")+
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12 ),
          panel.grid=element_line(colour=NA),
          panel.background = element_rect(fill = "#EDEDED"),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) 
  # coord_cartesian(ylim = c(0,25))
  # facet_wrap(~Majortype,scales = "free",nrow=2)
  g1
  
  g2 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y=Inflamatory.scores),stat="identity",width = 0.8) + #,size = 3
    geom_hline(aes(yintercept=median(Inflamatory.scores)),linetype="dashed") + 
    ylab("Inflamatory scores") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60,size = 0 , hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12 ),
          panel.grid=element_line(colour=NA),
          panel.background = element_rect(fill = "#EDEDED"),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) 
  
  g3 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y=apCAF.scores),stat="identity",width = 0.8) + #,size = 3
    geom_hline(aes(yintercept=median(apCAF.scores)),linetype="dashed") + 
    ylab("apCAF scores") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60,size = 8 , hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12 ),
          panel.grid=element_line(colour=NA),
          panel.background = element_rect(fill = "#EDEDED"),
          axis.line = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    scale_fill_manual(values = col) +
    scale_colour_manual(values = col) 
  
  # library(patchwork)
  g <- g1 / g2 / g3
  g
  ggsave("../result/Cell_nk_AUC/Fb_score_AUCell.pdf", g, height = 6, width =7)
 
  
  data<-as.data.frame(scRNAsub@assays$RNA@scale.data)
  mt<-data[All_genes,]
  mt<-na.omit(mt)
  #levels(seurat_data$seurat_clusters) <-LETTERS[1:9]
  colnames(mt)<-scRNAsub$seurat_clusters    
  mt_mean<-matrix(rep(1,nrow(mt)*length(unique(scRNAsub$seurat_clusters))),nrow=nrow(mt))
  rownames(mt_mean)<-rownames(mt)
  colnames(mt_mean)<-unique(scRNAsub$seurat_clusters)
  mt_mean<-as.data.frame(mt_mean)
  for(i in unique(scRNAsub$seurat_clusters)){
    mt_a <-mt[,colnames(mt) %in% i]
    mt_mean[,i] <- rowMeans(mt_a)
  }
  
  type<-test1$celltype 
  mt_mean<-subset(mt_mean,select= c(3,6,2,1,7,8,4,9,5))
  mt_mean<-edit(mt_mean)
  pheatmap::pheatmap(
    t(mt_mean),  scale = "row",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col  = c(11,16),
    # breaks = bk,
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100), #length(bk)
    # color = colorRampPalette(c("navy", "white", "firebrick3"))(30),
    border_color = "white",
    # scale = "row",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  
}


Statistics_wilcoxon_test{
  library(tidyverse)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  
  meta<-scRNAsub@meta.data
  meta$SPP1<-scRNAsub@assays$RNA@data["SPP1",]
  class(meta)
  meta$Group<-as.factor(meta$Group)
  wilcox_test = wilcox.test(x = meta$SPP1[meta$Group == "Tumor"], y = meta$SPP1[meta$Group == "Normal"])
  p.value = round(wilcox_test$p.value, 5)
  
}


Cell_nk_dotplot{
 DOT<-DotPlot(seurat_Fb, features = All_genes,cols = c('#1F78B4','#E31A1C'  ))+theme_cxf
 DOT_data <- as.data.frame(DOT$data)
 colnames(DOT_data)
 d1<-DotPlot(seurat_Fb, features = All_genes,cols = c('#1F78B4','#E31A1C'  ))+theme_cxf
 
 d2 <- ggplot(data = DOT_data)+ 
   geom_point(aes(x = id,y = features.plot, color = avg.exp.scaled, size = pct.exp))+
   scale_size_area(max_size = 6)+
   scale_color_gradientn(colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')) +
   cowplot::theme_cowplot() +
   theme_cxf+ 
   ylab("Cell Type") +
   xlab("") + coord_flip()  
 
 d3 <- ggplot(data = DOT_data)+ 
   geom_point(aes(x = id,y = features.plot, color = avg.exp.scaled, size = pct.exp))+
   scale_size_area(max_size = 6)+
   scale_color_gradientn(colors = c( 'gray95', '#26A69A', '#EC407A')) +
   cowplot::theme_cowplot() +
   theme_cxf+ 
   ylab("Cell Type") +
   xlab("") + coord_flip()  

 d1/d2/d3
 dev.off()
 
 }
 
 
 
apCAF2{ apCAF2<-apCAF
 umap<-apCAF2@reductions$umap@cell.embeddings
 meta<-apCAF2@meta.data
 meta<-data.frame(barcode = rownames(apCAF2@meta.data), clusters = apCAF2@meta.data[,"seurat_clusters"])
 meta$clusters <- as.factor(meta$clusters)
 meta$clusters <- droplevels(meta$clusters)
 levels(meta$clusters)
 umap_0<-meta$barcode[meta$clusters == 0]
 umap_1<-meta$barcode[meta$clusters == 1]
 umap_2<-meta$barcode[meta$clusters == 2]
 umap_3<-meta$barcode[meta$clusters == 3]
 umap[umap_0,]<-umap[umap_0,][sample(umap_0,length(umap_0)),]
 umap[umap_1,]<-umap[umap_1,][sample(umap_1,length(umap_1)),]
 umap[umap_2,]<-umap[umap_2,][sample(umap_2,length(umap_2)),]
 umap[umap_3,]<-umap[umap_3,][sample(umap_3,length(umap_3)),]
 apCAF2@reductions$umap@cell.embeddings<-umap
 FeaturePlot(apCAF2, features = c("CD74","SPP1","MSLN","PTPRC"),
             cols = c(rgb(225,225,225,10, maxColorValue = 255),
                      rgb(250,27,40,100,  maxColorValue = 255)))
 saveRDS(apCAF2,"d_apCAF2.rds")
}



apCAF_Trajectory_analysis_monocle3{
  #seurat_Fb <-readRDS("./d_apCAF2.Rds")
  #Transform Seurat object into cell_data_set object
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'HDF5Array',
                         'terra', 'ggrastr'))
  usethis::use_git_config(user.name = "windyeeeee", user.email = "chenxiongfengfj@qq.com")
  credentials::set_github_pat()
  devtools::install_github('cole-trapnell-lab/monocle3')
  library(monocle3)
  library(SeuratWrappers)
  
  Idents(seurat_Fb)
  seurat_FB_M <- subset(seurat_Fb,idents=c(0,3))
  seurat_FB_F <- subset(seurat_Fb,idents=c(1,2))  
  DimPlot(seurat_FB_M)
  seurat_FB_M$RNA_snn_res.0.1<-droplevels(seurat_FB_M$RNA_snn_res.0.1)
  seurat_FB_F$RNA_snn_res.0.1<-droplevels(seurat_FB_F$RNA_snn_res.0.1)
  seurat_Fb<-seurat_FB_M
  seurat_Fb<-seurat_FB_F
 
  
  Fb_monocle <- as.cell_data_set(seurat_Fb, assay = "RNA")
  Fb_monocle <- estimate_size_factors(Fb_monocle)
  
  
  
  #Get cell metadata
  colData(Fb_monocle)
  
  #Get gene metadata
  fData(Fb_monocle)
  rownames(fData(Fb_monocle))[1:10]
  fData(Fb_monocle)$gene_short_name <- rownames(fData(Fb_monocle))
  
  #Get counts
  counts(Fb_monocle)
  
  #Assign partitions
  recreate.partiion <- c(rep(1,length(Fb_monocle@colData@rownames)))
  names(recreate.partiion) <- Fb_monocle@colData@rownames
  recreate.partiion <- as.factor(recreate.partiion)
  Fb_monocle@clusters$UMAP$partitions <- recreate.partiion
  
  
  #Assign cluster info
  list_cluster <- seurat_Fb$RNA_snn_res.0.1
  Fb_monocle@clusters$UMAP$clusters <- list_cluster
  
  #Save umap structure
  Fb_monocle@int_colData@listData$reducedDims$UMAP <- seurat_Fb@reductions$umap@cell.embeddings
  
  #Learn trajectory graph
  Fb_monocle <- learn_graph(Fb_monocle, use_partition = F)
  plot_cells(Fb_monocle, color_cells_by = "cluster", label_groups_by_cluster = F, label_branch_points = F,
             label_roots = F, label_leaves = F, group_label_size = 0) + theme_bw()
  
  #Oder cells in pseudotime
  get_earliest_principal_node <- function(Fb_monocle, time_bin="2"){
    cell_ids <- which(colData(Fb_monocle)[, "RNA_snn_res.0.15"] == time_bin)
    closest_vertex <-Fb_monocle@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(Fb_monocle), ])
    root_pr_nodes <-igraph::V(principal_graph(Fb_monocle)[["UMAP"]])$name[as.numeric(names
                    (which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
  }
  
  Fb_monocle <- order_cells(Fb_monocle, root_pr_nodes =  root_pr_nodes)
  #Fb_monocle <- order_cells(Fb_monocle)    #可指定多个节点
  #Visualization
  plot_cells(Fb_monocle, color_cells_by = "pseudotime", label_groups_by_cluster = F, label_branch_points = F,
             label_roots = F, label_leaves = F, show_trajectory_graph = F) + theme_cxf
  
  pseudotime(Fb_monocle) # cells ordered by pseudotime
  Fb_monocle$monocle3_pseudotime <- pseudotime(Fb_monocle)
  data_pseudo <- as.data.frame(colData(Fb_monocle))
  ggplot(data_pseudo, aes(monocle3_pseudotime, reorder(RNA_snn_res.0.1, monocle3_pseudotime, median), fill = RNA_snn_res.0.1)) + 
    geom_boxplot() + theme_minimal()
  dev.off()
  #Fb_monocle_sample<-Fb_monocle[,sample(colnames(Fb_monocle),100000)]
  #Find genes with changing expression in pseudotime
  genes_T <- graph_test(Fb_monocle, neighbor_graph = "principal_graph", cores = 4)
  #top genes
  topgenes <- genes_T %>%
    arrange(q_value) %>%
    filter(status == "OK") 
  
  A <- topgenes %>% arrange(desc(morans_test_statistic))   %>% rownames()
  A<-A[1:10]
  A<-c("MSLN","PTPRC","HLA-DRA","SPP1","CD74","CRYAB")
  plot_genes_in_pseudotime(Fb_monocle[rowData(Fb_monocle)$gene_short_name %in% A],
                           color_cells_by="pseudotime",ncol=1, min_expr=0, cell_size=1) 
  dev.off()
  #Visualize pseudotime in seurat
  seurat_Fb$pseudotime <- pseudotime(Fb_monocle)
  FeaturePlot(seurat_Fb, features = "pseudotime", label = F,cols = c("lightgrey", "red")) + theme_cxf
  FeaturePlot(seurat_Fb, features = A[1:4], ncol=2) 
  
}



Cell_nk_featureplot{
 library("Nebulosa")
 #viridis  magma  cividis  inferno  plasma
P1 <- plot_density(apCAF2,reduction = "umap", c("CD74","SPP1","MSLN","PTPRC"), method = c("ks", "wkde"),
              adjust = 1, size = 1, shape = 16) +
              scale_color_gradientn(colors = c( 'darkblue','green3', '#FFFFB3', '#E31A1C'))  


P2 <-FeaturePlot(apCAF2, features = c("CD74","SPP1","MSLN","PTPRC"),
                cols = c(rgb(225,225,225,10, maxColorValue = 255),
                         rgb(250,27,40,100,  maxColorValue = 255)))

P3_1 <-FeaturePlot(apCAF2,features = c("CD74"), pt.size = 2, max.cutoff = 3)+
       scale_color_viridis("Temperature", option = "H",alpha = 1)
P3_2 <-FeaturePlot(apCAF2,features = c("SPP1"), pt.size = 2, max.cutoff = 3)+
       scale_color_viridis("Temperature", option = "H",alpha = 1)
P3_3 <-FeaturePlot(apCAF2,features = c("MSLN"), pt.size = 2, max.cutoff = 3)+
       scale_color_viridis("Temperature", option = "H",alpha = 1)  
P3_4 <-FeaturePlot(apCAF2,features = c("PTPRC"), pt.size = 2, max.cutoff = 3)+
       scale_color_viridis("Temperature", option = "H",alpha = 1) 
P3<-(P3_1|P3_2)/(P3_3|P3_4) 
 
DF<-as.data.frame(apCAF2@reductions[["umap"]]@cell.embeddings)
data<-as.data.frame(t(as.matrix(apCAF2@assays[["RNA"]]@data)))
DF$CD74<-ifelse(data$CD74>2,2,data$CD74)
DF$SPP1<-ifelse(DF$SPP1>2,2,DF$SPP1)
DF$MSLN<-ifelse(DF$MSLN>2,2,DF$MSLN)
DF$PTPRC<-ifelse(DF$PTPRC>2,2,DF$PTPRC)
P4_1 <-ggplot(DF, aes(umap_1,umap_2)) + 
       geom_point(aes(color=CD74),size=0.1)+
       scale_color_gradientn(colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')) + 
       #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
       theme_cxf
P4_2 <-ggplot(DF, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.1)+
  scale_color_gradientn(colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
P4_3 <-ggplot(DF, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.1, alpha = 0.2)+
  scale_color_gradientn(colors = c('gray90','lightgoldenrod1','firebrick1')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
P4_4 <-ggplot(DF, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.1, alpha = 0.2)+
  scale_color_gradientn(colors = c('gray90','lightgoldenrod1','firebrick1')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
P4<-(P4_1|P4_2)/(P4_3|P4_4) 



DF<-as.data.frame(apCAF2@reductions[["umap"]]@cell.embeddings)
data<-as.data.frame(t(as.matrix(apCAF2@assays[["RNA"]]@data)))
#DF$CD74<-ifelse(data$CD74>2,2,data$CD74)
DF$CD74<-data$CD74
DF$SPP1<-data$SPP1
DF$MSLN<-data$MSLN
DF$PTPRC<-data$PTPRC
df1<-DF[order(DF$CD74, decreasing = F),]
df2<-DF[order(DF$SPP1, decreasing = F),]
df3<-DF[order(DF$MSLN, decreasing = F),]
df4<-DF[order(DF$PTPRC, decreasing = F),]

G1_1<-ggplot(df1, aes(umap_1,umap_2)) + 
  geom_point(aes(color=CD74),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_2<-ggplot(df2, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_3<-ggplot(df3, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_4<-ggplot(df4, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1<-(G1_1|G1_2)/(G1_3|G1_4)  

#DF$CD74<-ifelse(data$CD74>2,2,data$CD74)
DF$CD74 <- ifelse(data$CD74>2,2,data$CD74)
DF$SPP1 <-ifelse(data$SPP1>2,2,data$SPP1)
DF$MSLN <-ifelse(data$MSLN>2,2,data$MSLN)
DF$PTPRC <-ifelse(data$PTPRC>2,2,data$PTPRC)
df1<-DF[order(DF$CD74, decreasing = F),]
df2<-DF[order(DF$SPP1, decreasing = F),]
df3<-DF[order(DF$MSLN, decreasing = F),]
df4<-DF[order(DF$PTPRC, decreasing = F),]
G2_1<-ggplot(df1, aes(umap_1,umap_2)) + 
  geom_point(aes(color=CD74),size=0.5,alpha = 0.2)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_2<-ggplot(df2, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.5,alpha = 0.2)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_3<-ggplot(df3, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.5,alpha = 0.2)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_4<-ggplot(df4, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.5,alpha = 0.2)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2<-(G2_1|G2_2)/(G2_3|G2_4)  

G1|G2
dev.off()
 
G2_1<-ggplot(df1, aes(umap_1,umap_2)) + 
  geom_point(aes(color=CD74),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_2<-ggplot(df2, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_3<-ggplot(df3, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_4<-ggplot(df4, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.05,alpha = 0.3)+
  scale_color_gradientn(colors = c('darkblue','green3', '#FFFFB3', '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2<-(G2_1|G2_2)/(G2_3|G2_4)  
dev.off()


G2_1<-ggplot(df1, aes(umap_1,umap_2)) + 
  geom_point(aes(color=CD74),size=0.001,alpha = 0.4)+
  scale_color_gradientn(colors = c('#216DC4','#AEC420', '#ECE4AE', '#F03F2B')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_2<-ggplot(df2, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.001,alpha = 0.4)+
  scale_color_gradientn(colors = c('#216DC4','#AEC420', '#ECE4AE', '#F03F2B')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_3<-ggplot(df3, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.001,alpha = 0.4)+
  scale_color_gradientn(colors = c('#216DC4','#AEC420', '#ECE4AE', '#F03F2B')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2_4<-ggplot(df4, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.001,alpha = 0.4)+
  scale_color_gradientn(colors = c('#216DC4','#AEC420', '#ECE4AE', '#F03F2B')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G2<-(G2_1|G2_2)/(G2_3|G2_4)  
dev.off()

 

G1_1<-ggplot(df1, aes(umap_1,umap_2)) + 
  geom_point(aes(color=CD74),size=0.01,alpha = 0.2)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_2<-ggplot(df2, aes(umap_1,umap_2)) + 
  geom_point(aes(color=SPP1),size=0.01,alpha = 0.2)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_3<-ggplot(df3, aes(umap_1,umap_2)) + 
  geom_point(aes(color=MSLN),size=0.01,alpha = 0.2)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1_4<-ggplot(df4, aes(umap_1,umap_2)) + 
  geom_point(aes(color=PTPRC),size=0.01,alpha = 0.2)+
  scale_color_gradientn(colors = c("gray90", '#E31A1C')) + 
  #scale_color_viridis("Temperature", option = "H",alpha = 1,discrete = F)+
  theme_cxf
G1<-(G1_1|G1_2)/(G1_3|G1_4)  
dev.off()

}

 
Cells_nk_genes{
  apCAF <- readRDS("d_apCAF2.rds") 
  apCAF <- seu_apCAF
  DF<-as.data.frame(apCAF@meta.data)
  data<-as.data.frame(t(as.matrix(apCAF@assays[["RNA"]]@data)))
  #DF$CD74<-ifelse(data$CD74>2,2,data$CD74)
  DF$SPP1<-data$SPP1
  DT = DF  %>% group_by(SampleID,Group)
  DT <- dplyr::summarise(DT,SPP1 = mean(SPP1))
  DT<-data.frame( DT)
  DimPlot(apCAF,split.by = "Group")
  FeaturePlot(apCAF,features="SPP1",split.by = "Group")
  
  
  update.packages("dplyr")
  install.packages("ggstatsplot")
  install.packages("palmerpenguins")
  library(ggstatsplot)
  library(palmerpenguins)
  library(tidyverse)
  
  DT$"log2(SPP1+1)"<-log2(DT$SPP1+1)
  DT <- drop_na(DT)
  ggplot2::geom_violin()
  
  plt <- ggbetweenstats(
                       data = DT,
                       x = Group,
                       y = `log2(SPP1+1)`,
                       type = "nonparametric",
                       tr = 0.5,
                       centrality.point.args = list(size = 3, color = "darkred"),
                       centrality.label.args = list(size = 3, nudge_x = 0.05, segment.linetype = 6,
                                                    min.segment.length = 0),
                       point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.01), alpha =
                                           0.2, size = 1, stroke = 0, na.rm = TRUE),
                                violin.args = list(width = 0.3, alpha = 0.5, na.rm = TRUE, lwd = 0.75, 
                                          aes(fill=Group),adjust =2 ),
                       boxplot.args = list(width = 0.05, alpha = 0.5, na.rm = TRUE, lwd = 0.75, 
                                           aes(fill=Group),fill="white"),
                       ggsignif.args = list(textsize = 14, tip_length = 0.01, na.rm = TRUE),
                       ggtheme =theme_minimal(),
                       package = "RColorBrewer",
                       palette = "Dark2"
            )
  
  plt+theme(axis.title.y=element_text(size=14),axis.text.y =element_text(size = 14),
          axis.text.x = element_text(size = 14 ),axis.ticks =element_line(size=0.5),
          strip.background = element_rect(colour = "white")) + 
      scale_fill_manual(values = col) +
      theme_cxf
}


cellchat_apCAF_T_V1{
  #rm(list=ls())
  d_apCAF <-readRDS("./data/d_apCAF2.rds")
  d_T <- readRDS("./data/d_T.rds")
  dim(d_apCAF); dim(d_T)
  plot(density(d_T@assays$RNA$data["CD3D",]))
  plot(density(d_apCAF@assays$RNA$data["CD74",]))
  length(colnames(d_T@assays$RNA$data)[d_T@assays$RNA$data["CD3D",]>0.5])          #425690
  length(colnames(d_apCAF@assays$RNA$data)[d_apCAF@assays$RNA$data["CD74",]>0.5])  #10070
  A<-colnames(d_T@assays$RNA$data)[d_T@assays$RNA$data["CD3D",]>0.5]
  d_T_subset <- d_T[,A]
  B<-colnames(d_apCAF@assays$RNA$data)[d_apCAF@assays$RNA$data["CD74",]>0.5]
  d_apCAF_subset <- d_apCAF[,B]
  table(d_T_subset$seurat_clusters)
  table(d_apCAF_subset$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_T_subset@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_T_subset<-d_T_subset[,samp$barcode]
  }
  apCAF_T <- merge(d_apCAF_subset,y = c(d_T_subset))
  saveRDS(apCAF_T,"./data/apCAF_T.rds")
  apCAF_T <- readRDS("./data/apCAF_T.rds")
  apCAF_T$seurat_clusters<-as.factor(apCAF_T$seurat_clusters)
  apCAF_T$seurat_clusters<- droplevels(apCAF_T$seurat_clusters)
  table(apCAF_T$seurat_clusters)
  DimPlot(apCAF_T)
  
  Idents(apCAF_T)<-"Group"
  apCAF_T_T<-subset(apCAF_T,idents = c("Normal"))
  table(apCAF_T_T$seurat_clusters)
  apCAF_T_T<-subset(apCAF_T,idents = c("Tumor"))
  table(apCAF_T_T$seurat_clusters)
  
apCAF_T_N{
  library(CellChat)#载入R包
  data.input <- LayerData(apCAF_T_N, assay = "RNA", layer = "data")
  levels(apCAF_T_N$seurat_clusters)<-c("M_apCAF1","F_apCAF1","F_apCAF2","M_apCAF2","CD16+NK","CD4+Teff","CD4+Tn","CD56+NK","CD8+Teff","CD8+Tn",  
                                        "Tex","Tfh","Th17","TIFIT3","TMCM6","TMKI67","Treg")
  apCAF_T_N$seurat_clusters<-fct_relevel(apCAF_T_N$seurat_clusters,"M_apCAF1","M_apCAF2","F_apCAF1","F_apCAF2" )
  
  identity <- subset(apCAF_T_N@meta.data, select = "seurat_clusters")
  levels(identity)
  cellchat_N <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
  CellChatDB <- CellChatDB.human
  cellchat_N@DB <- CellChatDB
  cellchat_N <- subsetData(cellchat_N)
  cellchat_N <- identifyOverExpressedGenes(cellchat_N)
  cellchat_N <- identifyOverExpressedInteractions(cellchat_N)
  cellchat_N <- projectData(cellchat_N, PPI.human)
  
  cellchat_N <- computeCommunProb(cellchat_N, raw.use = TRUE)
  cellchat_N <- filterCommunication(cellchat_N, min.cells = 3)
  df.net <- subsetCommunication(cellchat_N)
  A<-levels(apCAF_T_N$seurat_clusters)[c(1:4)]
  B<-levels(apCAF_T_N$seurat_clusters)[c(5:17)]
  df.net <- subsetCommunication(cellchat_N, sources.use = A, targets.use = B) 
  write.csv(df.net,"./result/cellchat_apCAF_T/N_df_net.csv")
  cellchat_N <- computeCommunProbPathway(cellchat_N)
  cellchat_N <- aggregateNet(cellchat_N)
  groupSize <- as.numeric(table(cellchat_N@idents))
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_N@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_N@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  P<-cellchat_N@net$count
  P2 <- P[A,B] + t(P[B,A])

  pheatmap::pheatmap(
    P2, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
    border_color = "white",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  
  ComplexHeatmap::pheatmap(
    P2, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
    border_color = "white",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  
  
  #pathways
  levels(cellchat_N@idents)            
  cellchat_N@netP$pathways             
  pathways<-names(table(df.net$pathway_name))
  # Heatmap

  for (i in 1:length(pathways)) {
    pathway  <- pathways[i]  
    p<-netVisual_heatmap(cellchat_N,signaling = pathway, 
                         color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                         sources.use = A,
                         targets.use = B,
                         remove.isolate = T,
                         font.size = 14,
                         font.size.title = 20)
    pdf(paste0("./result/cellchat_apCAF_T/N_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
    print(p)
    dev.off()
  }
  
  #L-R
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      pairLR <- extractEnrichedLR(cellchat_N, signaling = pathway, geneLR.return = FALSE)
      pairLR <- as.vector(unlist(pairLR))
      for(j in 1:length(pairLR)) {
        par(mar = c(0,0,0,0))
        LR.show <- pairLR[j] 
        pdf(paste0("./result/cellchat_apCAF_T/N_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
        netVisual_individual(cellchat_N,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                             #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                             signaling.name = LR.show,show.legend = T)
        dev.off()
      }
    } 
   
  #netVisual_bubble
  B<-netVisual_bubble(cellchat_N, sources.use = A, targets.use = B, signaling = pathways,
                   remove.isolate = T ,  thresh = 0.01,
                   font.size = 14,font.size.title = 16)
  B_data<-B$data
  saveRDS(B_data,"./result/cellchat_apCAF_T/N4_data")
  
  #plotGeneExpression
  for (i in 1:length(pathways)) {
    pathway  <- pathways[i]  
    P<-plotGeneExpression(cellchat_N, signaling = pathway)
    pdf(paste0("./result/cellchat_apCAF_T/N_GeneExpression/",pathway,".pdf"))
    print(P)
    dev.off()
  }
  
  saveRDS(cellchat_N,"./result/cellchat_apCAF_T/cellchat_N.rds")
  
  }
  
apCAF_T_T{
  library(CellChat)#载入R包
  data.input <- LayerData(apCAF_T_T, assay = "RNA", layer = "data")
  levels(apCAF_T_T$seurat_clusters)<-c("M_apCAF1","F_apCAF1","F_apCAF2","M_apCAF2","CD16+NK","CD4+Teff","CD4+Tn","CD56+NK","CD8+Teff","CD8+Tn",  
                                       "Tex","Tfh","Th17","TIFIT3","TMCM6","TMKI67","Treg")
  apCAF_T_T$seurat_clusters<-fct_relevel(apCAF_T_T$seurat_clusters,"M_apCAF1","M_apCAF2","F_apCAF1","F_apCAF2" )
  
  
  identity <- subset(apCAF_T_T@meta.data, select = "seurat_clusters")
  levels(identity)
  cellchat_T <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
  CellChatDB <- CellChatDB.human
  cellchat_T@DB <- CellChatDB
  cellchat_T <- subsetData(cellchat_T)
  cellchat_T <- identifyOverExpressedGenes(cellchat_T)
  cellchat_T <- identifyOverExpressedInteractions(cellchat_T)
  cellchat_T <- projectData(cellchat_T, PPI.human)
  
  cellchat_T <- computeCommunProb(cellchat_T, raw.use = TRUE)
  cellchat_T <- filterCommunication(cellchat_T, min.cells = 3)
  df.net <- subsetCommunication(cellchat_T)
  A<-levels(apCAF_T_T$seurat_clusters)[c(1:4)]
  B<-levels(apCAF_T_T$seurat_clusters)[c(5:17)]
  df.net <- subsetCommunication(cellchat_T, sources.use = A, targets.use = B) 
  write.csv(df.net,"./result/cellchat_apCAF_T/T_df_net.csv")
  cellchat_T <- computeCommunProbPathway(cellchat_T)
  cellchat_T <- aggregateNet(cellchat_T)
  groupSize <- as.numeric(table(cellchat_T@idents))
  
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat_T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat_T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  P<-cellchat_T@net$count
  P2 <- P[A,B] + t(P[B,A])
  pheatmap::pheatmap(
    P2, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
    border_color = "white",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  
  ComplexHeatmap::pheatmap(
    P2, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
    border_color = "white",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  
  
  #pathways
  levels(cellchat_T@idents)            
  cellchat_T@netP$pathways             
  pathways<-names(table(df.net$pathway_name))
  # Heatmap
  for (i in 1:length(pathways)) {
    pathway  <- pathways[i]  
    p<-netVisual_heatmap(cellchat_T,signaling = pathway, 
                         color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                         sources.use = A,
                         targets.use = B,
                         remove.isolate = T,
                         font.size = 14,
                         font.size.title = 20)
    pdf(paste0("./result/cellchat_apCAF_T/T_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
    print(p)
    dev.off()
  }
  
  #L-R
  for (i in 1:length(pathways)) {
    pathway  <- pathways[i]  
    pairLR <- extractEnrichedLR(cellchat_T, signaling = pathway, geneLR.return = FALSE)
    pairLR <- as.vector(unlist(pairLR))
    for(j in 1:length(pairLR)) {
      par(mar = c(0,0,0,0))
      LR.show <- pairLR[j] 
      pdf(paste0("./result/cellchat_apCAF_T/T_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
      netVisual_individual(cellchat_T,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                           #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                           signaling.name = LR.show,show.legend = T)
      dev.off()
    }
  } 
  
  #netVisual_bubble
  BB<-netVisual_bubble(cellchat_T, sources.use = A, targets.use = B, signaling = pathways,
                   remove.isolate = T ,  thresh = 0.01,
                   font.size = 14,font.size.title = 16)
  B_data<-BB$data
  saveRDS(B_data,"./result/cellchat_apCAF_T/T4_data")
  #plotGeneExpression
  for (i in 1:length(pathways)) {
    pathway  <- pathways[i]  
    P<-plotGeneExpression(cellchat_T, signaling = pathway)
    pdf(paste0("./result/cellchat_apCAF_T/T_GeneExpression/",pathway,".pdf"))
    print(P)
    dev.off()
  }
  
  saveRDS(cellchat_T,"./result/cellchat_apCAF_T/cellchat_T.rds")
}  
  
apCAF_T_NT{
  library(CellChat)
  library(patchwork)
  cellchat_N  <- readRDS("./result/cellchat_apCAF_T/cellchat_T.rds")
  cellchat_T  <- readRDS("./result/cellchat_apCAF_T/cellchat_T.rds")
  setwd("./result/cellchat_apCAF_T")
  data.dir <- './comparison'
  dir.create(data.dir)
  setwd("./comparison")
  
  cellchat.NL <- cellchat_N
  cellchat.LS <- cellchat_T
  object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  #Compare the total number of interactions and interaction strength
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  dev.off()
  #Differential number of interactions or interaction strength among different cell populations
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  dev.off()
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
  dev.off()
  gg1_data<-as.matrix(gg1@matrix)
  gg1_data<-as.data.frame(gg1_data)
  gg1_data<-gg1_data[A,B]
  ComplexHeatmap::pheatmap(
    gg1_data, 
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
    border_color = "white",
    cellwidth = 16, cellheight = 16, fontsize = 14)
  dev.off()
  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }
  dev.off()
  
  #Differential number of interactions or interaction strength among different cell types
  group.cellType <- c(rep("apCAF", 4), rep("T",13))
  group.cellType <- factor(group.cellType, levels = c("apCAF", "T"))
  object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  }
  dev.off()
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
  dev.off()
  #Compare the major sources and targets in 2D space
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  }
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
  patchwork::wrap_plots(plots = gg)
  dev.off()

  gg1 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[2])#, signaling.exclude = "MIF"
  #> Visualizing differential outgoing and incoming signaling changes from NL to LS
  #> The following `from` values were not present in `x`: 0
  #> The following `from` values were not present in `x`: 0, -1
  gg2 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[4], signaling.exclude = c("MIF"))
  #> Visualizing differential outgoing and incoming signaling changes from NL to LS
  #> The following `from` values were not present in `x`: 0, 2
  #> The following `from` values were not present in `x`: 0, -1
  gg1/gg2
  dev.off()
  
  #Part II: Identify the conserved and context-specific signaling pathways
  #Identify signaling groups based on their functional similarity
  cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
  #> Compute signaling network similarity for datasets 1 2
  cellchat <- netEmbedding(cellchat, type = "functional")
  #> Manifold learning of the signaling networks for datasets 1 2
  cellchat <- netClustering(cellchat, type = "functional", do.parallel = F)
  #> Classification learning of the signaling networks for datasets 1 2
  # Visualization in 2D-space
  NE<-netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
  #> 2D visualization of signaling networks from datasets 1 2
  NE_data<-NE$data
  write.csv(NE_data,"NT_10_data.csv")
  
  netVisual_embeddingPairwiseZoomIn(cellchat, 
                                    type = "functional", nCol = 2,
                                    dot.alpha = 0.8,show.legend =T,
                                    point.shape = c(21,24, 0, 23, 25, 10, 12))+
  theme_cxf
  
  dev.off()
  
  #Compute and visualize the pathway distance in the learned joint manifold
  RS<-rankSimilarity(cellchat, type = "functional")+ theme_cxf
  RS_data<-RS$data
  write.csv(RS_data,"NT_12_data.csv")
  
  #Identify and visualize the conserved and context-specific signaling pathways
  #Compare the overall information flow of each signaling pathway
  gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
  gg1 + gg2
  
  dev.off()
  #Compare outgoing (or incoming) signaling associated with each cell population
  library(ComplexHeatmap)
  #> Loading required package: grid
  #> ========================================
  #> ComplexHeatmap version 2.10.0
  #> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
  #> Github page: https://github.com/jokergoo/ComplexHeatmap
  #> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
  #> 
  #> If you use it in published research, please cite:
  #> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
  #>   genomic data. Bioinformatics 2016.
  #> 
  #> The new InteractiveComplexHeatmap package can directly export static 
  #> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
  #> 
  #> This message can be suppressed by:
  #>   suppressPackageStartupMessages(library(ComplexHeatmap))
  #> ========================================
  i = 1
  # combining all the identified signaling pathways from different datasets 
  pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12)
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  HT1_data<-ht1@matrix
  HT2_data<-ht2@matrix
  write.csv(HT1_data,"NT_14_data1.csv")
  write.csv(HT2_data,"NT_14_data2.csv")
    dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12, color.heatmap = "GnBu")
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12, color.heatmap = "GnBu")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  HT1_data<-ht1@matrix
  HT2_data<-ht2@matrix
  write.csv(HT1_data,"NT_15_data1.csv")
  write.csv(HT2_data,"NT_15_data2.csv")
  dev.off()
  
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12, color.heatmap = "OrRd")
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12, color.heatmap = "OrRd")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  HT1_data<-ht1@matrix
  HT2_data<-ht2@matrix
  write.csv(HT1_data,"NT_16_data1.csv")
  write.csv(HT2_data,"NT_16_data2.csv")
  dev.off()
  
  
  #Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
  #Identify dysfunctional signaling by comparing the communication probabities
  NV <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), angle.x = 45)+theme_cxf
  NV_data <- NV$data
  write.csv(NV_data,"NT_17_data.csv")
  dev.off()  
  
  gg1 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg2 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
  #> Comparing communications on a merged object
  gg1 + gg2
  dev.off()  
  
  #Identify dysfunctional signaling by using differential expression analysis
  # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
  pos.dataset = "LS"
  # define a char name used for storing the results of differential expression analysis
  features.name = pos.dataset
  # perform differential expression analysis
  cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
  #> Use the joint cell labels from the merged CellChat object
  # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
  net <- netMappingDEG(cellchat, features.name = features.name)
  # extract the ligand-receptor pairs with upregulated ligands in LS
  net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
  # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
  net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
  
  gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
  gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
  
  pairLR.use.up = net.up[, "interaction_name", drop = F]
  gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  #> Comparing communications on a merged object
  pairLR.use.down = net.down[, "interaction_name", drop = F]
  gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
  #> Comparing communications on a merged object
  gg1 + gg2
  gg1_data<-gg1$data
  gg2_data<-gg2$data
  write.csv(gg1_data,"NT_19_data1.csv")
  write.csv(gg2_data,"NT_19_data2.csv")
  dev.off()
  
  # Chord diagram
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_chord_gene(object.list[[2]], sources.use = A, targets.use = B, slot.name = 'net', net = net.up, lab.cex = 0.2, small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
  netVisual_chord_gene(object.list[[1]], sources.use = A, targets.use = B, slot.name = 'net', net = net.down, lab.cex = 0.2, small.gap = 2, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
  dev.off()
  
  # visualize the enriched ligands in the first condition
  #install.packages('wordcloud')
  C1<-computeEnrichmentScore(net.down, species = 'human')
  # visualize the enriched ligands in the second condition
  C2<-computeEnrichmentScore(na.omit(net.up), species = 'human')
  write.csv(net.down,"NT_21_netdown.csv")
  write.csv(net.up,"NT_22_netup.csv")
  dev.off()
  
  
 # Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
  pathways.show <- c("SPP1") 
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  }
  dev.off()
  pathways.show <- c("SPP1") 
  par(mfrow = c(1,2), xpd=TRUE)
  ht <- list()
  for (i in 1:length(object.list)) {
    ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
  }
  ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
  dev.off()
  # Chord diagram
  pathways.show <- c("SPP1") 
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
  }
  dev.off()
  # show all the significant signaling pathways from fibroblast to immune cells
  par(mfrow = c(1, 2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, sources.use = A, targets.use = B,
                        layout = "chord",
                         signaling.name = paste(pathways.show, names(object.list)[i])
                         )
  }
  dev.off()
 
  #Part V: Compare the signaling gene expression distribution between different datasets
  cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
  plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
  plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", colors.ggplot = T)
  saveRDS(cellchat, file = "cellchat_comparisonAnalysis_human_N_vs_T.rds")
}  
  
}  
  
cellchat_apCAF_T_V2{
  #rm(list=ls())
  # d_apCAF <-readRDS("./data/d_apCAF2.rds")
  # d_T <- readRDS("./data/d_T.rds")
  # 
  # table(d_T$seurat_clusters)
  # table(d_apCAF$seurat_clusters)
  # d_apCAF$seurat_clusters<-droplevels(d_apCAF$seurat_clusters)
  # if(T){
  #   library(dplyr)
  #   meta<-as.data.frame(d_T@meta.data)
  #   df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
  #   samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
  #   # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
  #   #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
  #   # samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
  #   d_T_subset<-d_T[,samp$barcode]
  # }
  # apCAF_T <- merge(d_apCAF,y = c(d_T_subset))
  # saveRDS(apCAF_T,"./data/apCAF_T2.rds")
  # 
  
  
  apCAF_T <- readRDS("./data/apCAF_T.rds")
  table(apCAF_T$seurat_clusters)
  apCAF_T$seurat_clusters<-as.factor(apCAF_T$seurat_clusters)
  levels(apCAF_T$seurat_clusters)<-c("M-apCAF","F-apCAF","F-apCAF","M-apCAF","CD16+NK","CD4Teff","CD4Tn","CD56+NK","CD8Teff","CD8Tn", 
                                     "Tex","Tfh","Th17","TIFIT3","TMCM6","TMKI67","Treg")
  table(apCAF_T$seurat_clusters)
  library(tidyverse)
  apCAF_T$barcode<-colnames(apCAF_T)
  samp <- apCAF_T@meta.data %>% group_by(seurat_clusters) %>% sample_n(size=3000)  
  apCAF_T<-apCAF_T[,samp$barcode]
  
  
  Idents(apCAF_T)<-"Group"
  apCAF_T_N<-subset(apCAF_T,idents = c("Normal"))
  table(apCAF_T_N$seurat_clusters)
  apCAF_T_T<-subset(apCAF_T,idents = c("Tumor"))
  table(apCAF_T_T$seurat_clusters)

  
  
  apCAF_T_N{
    library(CellChat)#载入R包
    data.input <- LayerData(apCAF_T_N, assay = "RNA", layer = "data")
    identity <- subset(apCAF_T_N@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_N <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_N@DB <- CellChatDB
    cellchat_N <- subsetData(cellchat_N)
    cellchat_N <- identifyOverExpressedGenes(cellchat_N)
    cellchat_N <- identifyOverExpressedInteractions(cellchat_N)
    cellchat_N <- projectData(cellchat_N, PPI.human)
    
    cellchat_N <- computeCommunProb(cellchat_N,   type = c("truncatedMean"),raw.use = TRUE)
    cellchat_N <- filterCommunication(cellchat_N, min.cells = 3)
    df.net <- subsetCommunication(cellchat_N)
    A<-levels(apCAF_T_N$seurat_clusters)[c(1:2)]
    B<-levels(apCAF_T_N$seurat_clusters)[c(3:15)]
    df.net <- subsetCommunication(cellchat_N, thresh = 0.05,sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/cellchat_apCAF_T/V2/N_df_net.csv")
    cellchat_N <- computeCommunProbPathway(cellchat_N)
    cellchat_N <- aggregateNet(cellchat_N)
    groupSize <- as.numeric(table(cellchat_N@idents))
    
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat_N@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat_N@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    P<-cellchat_N@net$count
    P2 <- P[A,B] + t(P[B,A])
    
    pheatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    ComplexHeatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    
    #pathways
    levels(cellchat_N@idents)            
    cellchat_N@netP$pathways             
    pathways<-names(table(df.net$pathway_name))
    # Heatmap
    
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      p<-netVisual_heatmap(cellchat_N,signaling = pathway, 
                           color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                           sources.use = A,
                           targets.use = B,
                           remove.isolate = T,
                           font.size = 14,
                           font.size.title = 20)
      pdf(paste0("./result/cellchat_apCAF_T/V2/N_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
      print(p)
      dev.off()
    }
    
    #L-R
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      pairLR <- extractEnrichedLR(cellchat_N, signaling = pathway, geneLR.return = FALSE)
      pairLR <- as.vector(unlist(pairLR))
      for(j in 1:length(pairLR)) {
        par(mar = c(0,0,0,0))
        LR.show <- pairLR[j] 
        pdf(paste0("./result/cellchat_apCAF_T/V2/N_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
        netVisual_individual(cellchat_N,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                             #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                             signaling.name = LR.show,show.legend = T)
        dev.off()
      }
    } 
    
    #netVisual_bubble
    B<-netVisual_bubble(cellchat_N, sources.use = A, targets.use = B, signaling = pathways,
                        remove.isolate = T ,  thresh = 0.01,
                        font.size = 14,font.size.title = 16)
    B_data<-B$data
    saveRDS(B_data,"./result/cellchat_apCAF_T/V2/N4_data")
    write.csv(B_data,"./result/cellchat_apCAF_T/V2/N4_data.csv")
    #plotGeneExpression
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      P<-plotGeneExpression(cellchat_N, signaling = pathway)
      pdf(paste0("./result/cellchat_apCAF_T/V2/N_GeneExpression/",pathway,".pdf"))
      print(P)
      dev.off()
    }
    
    saveRDS(cellchat_N,"./result/cellchat_apCAF_T/V2/cellchat_N.rds")
    
  }
  
  apCAF_T_T{
    library(CellChat)#载入R包
    data.input <- LayerData(apCAF_T_T, assay = "RNA", layer = "data")
    
    identity <- subset(apCAF_T_T@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_T <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_T@DB <- CellChatDB
    cellchat_T <- subsetData(cellchat_T)
    cellchat_T <- identifyOverExpressedGenes(cellchat_T)
    cellchat_T <- identifyOverExpressedInteractions(cellchat_T)
    cellchat_T <- projectData(cellchat_T, PPI.human)
    
    cellchat_T <- computeCommunProb(cellchat_T, type = c("truncatedMean"), raw.use = TRUE)
    cellchat_T <- filterCommunication(cellchat_T, min.cells = 3)
    df.net <- subsetCommunication(cellchat_T)
    A<-levels(apCAF_T_T$seurat_clusters)[c(1:2)]
    B<-levels(apCAF_T_T$seurat_clusters)[c(3:15)]
    df.net <- subsetCommunication(cellchat_T, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/cellchat_apCAF_T/V2/T_df_net.csv")
    cellchat_T <- computeCommunProbPathway(cellchat_T)
    cellchat_T <- aggregateNet(cellchat_T)
    groupSize <- as.numeric(table(cellchat_T@idents))
    
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat_T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat_T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    P<-cellchat_T@net$count
    P2 <- P[A,B] + t(P[B,A])
    pheatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    ComplexHeatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    
    #pathways
    levels(cellchat_T@idents)            
    cellchat_T@netP$pathways             
    pathways<-names(table(df.net$pathway_name))
    # Heatmap
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      p<-netVisual_heatmap(cellchat_T,signaling = pathway, 
                           color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                           sources.use = A,
                           targets.use = B,
                           remove.isolate = T,
                           font.size = 14,
                           font.size.title = 20)
      pdf(paste0("./result/cellchat_apCAF_T/V2/T_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
      print(p)
      dev.off()
    }
    
    #L-R
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      pairLR <- extractEnrichedLR(cellchat_T, signaling = pathway, geneLR.return = FALSE)
      pairLR <- as.vector(unlist(pairLR))
      for(j in 1:length(pairLR)) {
        par(mar = c(0,0,0,0))
        LR.show <- pairLR[j] 
        pdf(paste0("./result/cellchat_apCAF_T/V2/T_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
        netVisual_individual(cellchat_T,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                             #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                             signaling.name = LR.show,show.legend = T)
        dev.off()
      }
    } 
    
    #netVisual_bubble
    BB<-netVisual_bubble(cellchat_T, sources.use = A, targets.use = B, signaling = pathways,
                         remove.isolate = T ,  thresh = 0.01,
                         font.size = 14,font.size.title = 16)
    B_data<-BB$data
    #saveRDS(B_data,"./result/cellchat_apCAF_T/V2/T4_data")
    write.csv(B_data,"./result/cellchat_apCAF_T/V2/T4_data.csv")
    #plotGeneExpression
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      P<-plotGeneExpression(cellchat_T, signaling = pathway)
      pdf(paste0("./result/cellchat_apCAF_T/T_GeneExpression/",pathway,".pdf"))
      print(P)
      dev.off()
    }
    
    saveRDS(cellchat_T,"./result/cellchat_apCAF_T/V2/cellchat_T.rds")
  }  
  
  apCAF_T_NT{
    library(CellChat)
    library(patchwork)
    cellchat_N  <- readRDS("./result/cellchat_apCAF_T/V2/cellchat_N.rds")
    cellchat_T  <- readRDS("./result/cellchat_apCAF_T/V2/cellchat_T.rds")
    setwd("./result/cellchat_apCAF_T/V2")
    data.dir <- './comparison'
    dir.create(data.dir)
    setwd("./comparison")
    
    cellchat.NL <- cellchat_N
    cellchat.LS <- cellchat_T
    object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    #Compare the total number of interactions and interaction strength
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg1 + gg2
    dev.off()
    #Differential number of interactions or interaction strength among different cell populations
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    dev.off()
    gg1 <- netVisual_heatmap(cellchat)
    #> Do heatmap based on a merged object
    gg2 <- netVisual_heatmap(cellchat, measure = "weight")
    #> Do heatmap based on a merged object
    gg1 + gg2
    dev.off()
    gg1_data<-as.matrix(gg1@matrix)
    gg1_data<-as.data.frame(gg1_data)
    gg1_data<-gg1_data[A,B]
    ComplexHeatmap::pheatmap(
      gg1_data, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    dev.off()
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    dev.off()
    
    #Differential number of interactions or interaction strength among different cell types
    group.cellType <- c(rep("apCAF", 2), rep("T",13))
    group.cellType <- factor(group.cellType, levels = c("apCAF", "T"))
    object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    dev.off()
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
    dev.off()
    #Compare the major sources and targets in 2D space
    num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
    weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
    gg <- list()
    for (i in 1:length(object.list)) {
      object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
    }
    for (i in 1:length(object.list)) {
      gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
    }
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    patchwork::wrap_plots(plots = gg)
    dev.off()
    
    gg1 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[1])#, signaling.exclude = "MIF"
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0
    #> The following `from` values were not present in `x`: 0, -1
    gg2 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[2], signaling.exclude = c("MIF"))
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0, 2
    #> The following `from` values were not present in `x`: 0, -1
    gg1/gg2
    dev.off()
    
    #Part II: Identify the conserved and context-specific signaling pathways
    #Identify signaling groups based on their functional similarity
    cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
    #> Compute signaling network similarity for datasets 1 2
    cellchat <- netEmbedding(cellchat, type = "functional")
    #> Manifold learning of the signaling networks for datasets 1 2
    cellchat <- netClustering(cellchat, type = "functional", do.parallel = F)
    #> Classification learning of the signaling networks for datasets 1 2
    # Visualization in 2D-space
    NE<-netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
    #> 2D visualization of signaling networks from datasets 1 2
    NE_data<-NE$data
    write.csv(NE_data,"NT_9_data.csv")
    
    netVisual_embeddingPairwiseZoomIn(cellchat, 
                                      type = "functional", nCol = 2,
                                      dot.alpha = 0.8,show.legend =T,
                                      point.shape = c(21,24, 0, 23, 25, 10, 12))+
      theme_cxf
    
    dev.off()
    
    #Compute and visualize the pathway distance in the learned joint manifold
    RS<-rankSimilarity(cellchat, type = "functional")+ theme_cxf
    RS_data<-RS$data
    write.csv(RS_data,"NT_11_data.csv")
    
    #Identify and visualize the conserved and context-specific signaling pathways
    #Compare the overall information flow of each signaling pathway
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg1 + gg2
    
    dev.off()
    #Compare outgoing (or incoming) signaling associated with each cell population
    library(ComplexHeatmap)
    #> Loading required package: grid
    #> ========================================
    #> ComplexHeatmap version 2.10.0
    #> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    #> Github page: https://github.com/jokergoo/ComplexHeatmap
    #> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    #> 
    #> If you use it in published research, please cite:
    #> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    #>   genomic data. Bioinformatics 2016.
    #> 
    #> The new InteractiveComplexHeatmap package can directly export static 
    #> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    #> 
    #> This message can be suppressed by:
    #>   suppressPackageStartupMessages(library(ComplexHeatmap))
    #> ========================================
    i = 1
    # combining all the identified signaling pathways from different datasets 
    pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18)
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"NT_13_data1.csv")
    write.csv(HT2_data,"NT_13_data2.csv")
    dev.off()
    
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18, color.heatmap = "GnBu")
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18, color.heatmap = "GnBu")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"NT_14_data1.csv")
    write.csv(HT2_data,"NT_14_data2.csv")
    dev.off()
    
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 18, color.heatmap = "OrRd")
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 18, color.heatmap = "OrRd")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"NT_16_data1.csv")
    write.csv(HT2_data,"NT_16_data2.csv")
    dev.off()
    
    
    #Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
    #Identify dysfunctional signaling by comparing the communication probabities
    #single datasets

    
    #merged object

    
    NV <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1,2), angle.x = 45)+theme_cxf
    NV_data <- NV$data
    write.csv(NV_data,"NT_16_data.csv")
    dev.off()  
    
    gg1 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg2 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg1 + gg2
    dev.off()  
    gg1_data <- gg1$data
    gg2_data <- gg2$data
    write.csv(gg1_data ,"NT_17_data1.csv")
    write.csv(gg1_data ,"NT_17_data2.csv")
    
    
    #Identify dysfunctional signaling by using differential expression analysis
    # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
    pos.dataset = "LS"
    # define a char name used for storing the results of differential expression analysis
    features.name = pos.dataset
    # perform differential expression analysis
    
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
    #> Use the joint cell labels from the merged CellChat object
    # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
    net <- netMappingDEG(cellchat, features.name = features.name)
    # extract the ligand-receptor pairs with upregulated ligands in LS
    net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.1, receptor.logFC = NULL)
    # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
    net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
    
    gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
    gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
    
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    pairLR.use.down = net.down[, "interaction_name", drop = F]
    gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    gg1 + gg2
    gg1_data<-gg1$data
    gg2_data<-gg2$data
    write.csv(gg1_data,"NT_19_data1.csv")
    write.csv(gg2_data,"NT_19_data2.csv")
    dev.off()
    
    # Chord diagram
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_chord_gene(object.list[[2]], sources.use = A, targets.use = B, slot.name = 'net', net = net.up, lab.cex = 0.2, small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    netVisual_chord_gene(object.list[[1]], sources.use = A, targets.use = B, slot.name = 'net', net = net.down, lab.cex = 0.2, small.gap = 2, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    dev.off()
    
    # visualize the enriched ligands in the first condition
    #install.packages('wordcloud')
    C1<-computeEnrichmentScore(net.down, species = 'human')
    # visualize the enriched ligands in the second condition
    C2<-computeEnrichmentScore(na.omit(net.up), species = 'human')
    write.csv(net.down,"NT_21_netdown.csv")
    write.csv(net.up,"NT_22_netup.csv")
    dev.off()
    
    
    # Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
    pathways.show <- c("SPP1") 
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
    }
    dev.off()
    pathways.show <- c("SPP1") 
    par(mfrow = c(1,2), xpd=TRUE)
    ht <- list()
    for (i in 1:length(object.list)) {
      ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
    }
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    dev.off()
    # Chord diagram
    pathways.show <- c("SPP1") 
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
    }
    dev.off()
    # show all the significant signaling pathways from fibroblast to immune cells
    par(mfrow = c(1, 2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, sources.use = A, targets.use = B,
                          layout = "chord",
                          signaling.name = paste(pathways.show, names(object.list)[i])
      )
    }
    dev.off()
    
    #Part V: Compare the signaling gene expression distribution between different datasets
    cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
    plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
    plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", colors.ggplot = T)
    saveRDS(cellchat, file = "cellchat_comparisonAnalysis_human_N_vs_T.rds")
  }  
  
}   


cellchat_CAF_T{
  #rm(list=ls())
  d_CAF <-readRDS("./data/d_Fb.rds")
  d_T <- readRDS("./data/d_T.rds")
  dim(d_CAF); dim(d_T)
  # plot(density(d_T@assays$RNA$data["CD3D",]))
  # plot(density(d_CAF@assays$RNA$data["CD74",]))
  # length(colnames(d_T@assays$RNA$data)[d_T@assays$RNA$data["CD3D",]>0.5])          #425690
  # length(colnames(d_apCAF@assays$RNA$data)[d_apCAF@assays$RNA$data["CD74",]>0.5])  #10070
  # A<-colnames(d_T@assays$RNA$data)[d_T@assays$RNA$data["CD3D",]>0.5]
  # d_T_subset <- d_T[,A]
  # B<-colnames(d_apCAF@assays$RNA$data)[d_apCAF@assays$RNA$data["CD74",]>0.5]
  # d_apCAF_subset <- d_apCAF[,B]
  table(d_T$seurat_clusters)
  table(d_CAF$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_T@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_T_subset<-d_T[,samp$barcode]
  }
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_CAF@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    # samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_CAF_subset<-d_CAF[,samp$barcode]
  }
  
  CAF_T <- merge(d_CAF_subset,y = c(d_T_subset))
  saveRDS(CAF_T,"./data/CAF_T.rds")
  
  
  
  CAF_T <- readRDS("./data/CAF_T.rds")
  table(CAF_T$seurat_clusters)
  CAF_T$seurat_clusters<-as.factor(CAF_T$seurat_clusters)
  levels(CAF_T$seurat_clusters)
  CAF_T$seurat_clusters<-factor(CAF_T$seurat_clusters,levels=c(
    "Cano-iCAF","Ribo-iCAF","myCAF","DPT+CAF","MAGI2+CAF","MKI67+CAF","F-apCAF","M-apCAF","Pericytes",
    "CD16+NK","CD56+NK","CD8Teff","CD8Tn", "Tex","TIFIT3","TMCM6","TMKI67",
    "CD4Teff","CD4Tn","Tfh","Th17","Treg"
    ))
  CAF_T$seurat_clusters<- droplevels(CAF_T$seurat_clusters)
  table(CAF_T$seurat_clusters)
  DimPlot(CAF_T)
  
  Idents(CAF_T)<-"Group"
  CAF_T_N<-subset(CAF_T,idents = c("Normal"))
  table(CAF_T_N$seurat_clusters)
  CAF_T_T<-subset(CAF_T,idents = c("Tumor"))
  table(CAF_T_T$seurat_clusters)
  
  CAF_T_N{
    library(CellChat)#载入R包
    data.input <- LayerData(CAF_T_N, assay = "RNA", layer = "data")
    levels(CAF_T_N$seurat_clusters)
    identity <- subset(CAF_T_N@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_N <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_N@DB <- CellChatDB
    cellchat_N <- subsetData(cellchat_N)
    cellchat_N <- identifyOverExpressedGenes(cellchat_N)
    cellchat_N <- identifyOverExpressedInteractions(cellchat_N)
    cellchat_N <- projectData(cellchat_N, PPI.human)
    

    cellchat_N <- computeCommunProb(cellchat_N, type = c("truncatedMean"),raw.use = TRUE)
    cellchat_N <- filterCommunication(cellchat_N, min.cells = 3)
    df.net <- subsetCommunication(cellchat_N)
    A<-levels(CAF_T_N$seurat_clusters)[c(1:9)]
    B<-levels(CAF_T_N$seurat_clusters)[c(10:22)]
    df.net <- subsetCommunication(cellchat_N, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/cellchat_CAF_T/N_df_net.csv")
    cellchat_N <- computeCommunProbPathway(cellchat_N)
    cellchat_N <- aggregateNet(cellchat_N)
    groupSize <- as.numeric(table(cellchat_N@idents))
    
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat_N@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat_N@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    P<-cellchat_N@net$count
    P2 <- P[A,B] + t(P[B,A])
    
    pheatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    ComplexHeatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    
    #pathways
    levels(cellchat_N@idents)            
    cellchat_N@netP$pathways             
    pathways<-names(table(df.net$pathway_name))
    # Heatmap
    
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      p<-netVisual_heatmap(cellchat_N,signaling = pathway, 
                           color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                           sources.use = A,
                           targets.use = B,
                           remove.isolate = T,
                           font.size = 14,
                           font.size.title = 20)
      pdf(paste0("./result/cellchat_CAF_T/N_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
      print(p)
      dev.off()
    }
    
    #L-R
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      pairLR <- extractEnrichedLR(cellchat_N, signaling = pathway, geneLR.return = FALSE)
      pairLR <- as.vector(unlist(pairLR))
      for(j in 1:length(pairLR)) {
        par(mar = c(0,0,0,0))
        LR.show <- pairLR[j] 
        pdf(paste0("./result/cellchat_CAF_T/N_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
        netVisual_individual(cellchat_N,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                             #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                             signaling.name = LR.show,show.legend = T)
        dev.off()
      }
    } 
    
    #netVisual_bubble
    B<-netVisual_bubble(cellchat_N, sources.use = A, targets.use = B, signaling = pathways,
                        remove.isolate = T ,  thresh = 0.01,
                        font.size = 14,font.size.title = 16)
    B_data<-B$data
    saveRDS(B_data,"./result/cellchat_CAF_T/N_P3_data")
    
    #plotGeneExpression
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      P<-plotGeneExpression(cellchat_N, signaling = pathway)
      pdf(paste0("./result/cellchat_CAF_T/N_GeneExpression/",pathway,".pdf"))
      print(P)
      dev.off()
    }
    
    saveRDS(cellchat_N,"./result/cellchat_CAF_T/cellchat_N.rds")
    
  }
  
  CAF_T_T{
    library(CellChat)#载入R包
    data.input <- LayerData(CAF_T_T, assay = "RNA", layer = "data")
    identity <- subset(CAF_T_T@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_T <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_T@DB <- CellChatDB
    cellchat_T <- subsetData(cellchat_T)
    cellchat_T <- identifyOverExpressedGenes(cellchat_T)
    cellchat_T <- identifyOverExpressedInteractions(cellchat_T)
    cellchat_T <- projectData(cellchat_T, PPI.human)
    
   
    cellchat_T <- computeCommunProb(cellchat_T, type = c("truncatedMean"), raw.use = TRUE)
    cellchat_T <- filterCommunication(cellchat_T, min.cells = 3)
    df.net <- subsetCommunication(cellchat_T)
    A<-levels(CAF_T_T$seurat_clusters)[c(1:9)]
    B<-levels(CAF_T_T$seurat_clusters)[c(10:22)]
    df.net <- subsetCommunication(cellchat_T, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/cellchat_CAF_T/T_df_net.csv")
    cellchat_T <- computeCommunProbPathway(cellchat_T)
    cellchat_T <- aggregateNet(cellchat_T)
    groupSize <- as.numeric(table(cellchat_T@idents))
    
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat_T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat_T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    P<-cellchat_T@net$count
    P2 <- P[A,B] + t(P[B,A])
    pheatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    ComplexHeatmap::pheatmap(
      P2, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    
    
    #pathways
    levels(cellchat_T@idents)            
    cellchat_T@netP$pathways             
    pathways<-names(table(df.net$pathway_name))
    # Heatmap
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      p<-netVisual_heatmap(cellchat_T,signaling = pathway, 
                           color.heatmap = "Reds", #,c("#2166ac", "#b2182b")
                           sources.use = A,
                           targets.use = B,
                           remove.isolate = T,
                           font.size = 14,
                           font.size.title = 20)
      pdf(paste0("./result/cellchat_CAF_T/T_pathways_heatmap/",pathway,".pdf"), width  = 10,height = 10)
      print(p)
      dev.off()
    }
    
    #L-R
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      pairLR <- extractEnrichedLR(cellchat_T, signaling = pathway, geneLR.return = FALSE)
      pairLR <- as.vector(unlist(pairLR))
      for(j in 1:length(pairLR)) {
        par(mar = c(0,0,0,0))
        LR.show <- pairLR[j] 
        pdf(paste0("./result/cellchat_CAF_T/T_lr_chord/",pathway,j,".pdf"),width = 10, height = 10)
        netVisual_individual(cellchat_T,signaling = pathway,pairLR.use = LR.show,layout = "chord",
                             #cell.order = c("apCAF","iCAF","myCAF","Epithelial Cancer Cells","Mesenchymal Cancer Cells"),
                             signaling.name = LR.show,show.legend = T)
        dev.off()
      }
    } 
    
    #netVisual_bubble
    BB<-netVisual_bubble(cellchat_T, sources.use = A, targets.use = B, signaling = pathways,
                         remove.isolate = T ,  thresh = 0.01,
                         font.size = 14,font.size.title = 16)
    B_data<-BB$data
    saveRDS(B_data,"./result/cellchat_CAF_T/T_P3_data")
    #plotGeneExpression
    for (i in 1:length(pathways)) {
      pathway  <- pathways[i]  
      P<-plotGeneExpression(cellchat_T, signaling = pathway)
      pdf(paste0("./result/cellchat_CAF_T/T_GeneExpression/",pathway,".pdf"))
      print(P)
      dev.off()
    }
    
    saveRDS(cellchat_T,"./result/cellchat_CAF_T/cellchat_T.rds")
  }  
  
  apCAF_T_NT{
    library(CellChat)
    library(patchwork)
    cellchat_N  <- readRDS("./result/cellchat_CAF_T/cellchat_N.rds")
    cellchat_T  <- readRDS("./result/cellchat_CAF_T/cellchat_T.rds")
    setwd("./result/cellchat_CAF_T")
    data.dir <- './comparison'
    dir.create(data.dir)
    setwd("./comparison")
    
    cellchat.NL <- cellchat_N
    cellchat.LS <- cellchat_T
    object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    #Compare the total number of interactions and interaction strength
    gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
    gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
    gg1 + gg2
    dev.off()
    #Differential number of interactions or interaction strength among different cell populations
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
    dev.off()
    gg1 <- netVisual_heatmap(cellchat)
    #> Do heatmap based on a merged object
    gg2 <- netVisual_heatmap(cellchat, measure = "weight")
    #> Do heatmap based on a merged object
    gg1 + gg2
    dev.off()
    gg1_data<-as.matrix(gg1@matrix)
    gg1_data<-as.data.frame(gg1_data)
    gg1_data<-gg1_data[A,B]
    ComplexHeatmap::pheatmap(
      gg1_data, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    dev.off()
    weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    dev.off()
    
    #Differential number of interactions or interaction strength among different cell types
    group.cellType <- c(rep("apCAF", 4), rep("T",13))
    group.cellType <- factor(group.cellType, levels = c("apCAF", "T"))
    object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
    }
    dev.off()
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
    netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
    dev.off()
    #Compare the major sources and targets in 2D space
    num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
    weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
    gg <- list()
    for (i in 1:length(object.list)) {
      object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
    }
    for (i in 1:length(object.list)) {
      gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
    }
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    patchwork::wrap_plots(plots = gg)
    dev.off()
    
    gg1 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[3])#, signaling.exclude = "MIF"
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0
    #> The following `from` values were not present in `x`: 0, -1
    gg2 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = A[4], signaling.exclude = c("MIF"))
    #> Visualizing differential outgoing and incoming signaling changes from NL to LS
    #> The following `from` values were not present in `x`: 0, 2
    #> The following `from` values were not present in `x`: 0, -1
    gg1/gg2
    dev.off()
    
    #Part II: Identify the conserved and context-specific signaling pathways
    #Identify signaling groups based on their functional similarity
    cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
    #> Compute signaling network similarity for datasets 1 2
    cellchat <- netEmbedding(cellchat, type = "functional")
    #> Manifold learning of the signaling networks for datasets 1 2
    cellchat <- netClustering(cellchat, type = "functional", do.parallel = F)
    #> Classification learning of the signaling networks for datasets 1 2
    # Visualization in 2D-space
    NE<-netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
    #> 2D visualization of signaling networks from datasets 1 2
    NE_data<-NE$data
    write.csv(NE_data,"P8_data.csv")
    
    netVisual_embeddingPairwiseZoomIn(cellchat, 
                                      type = "functional", nCol = 2,
                                      dot.alpha = 0.8,show.legend =T,
                                      point.shape = c(21,24, 0, 23, 25, 10, 12))+
      theme_cxf
    
    dev.off()
    
    #Compute and visualize the pathway distance in the learned joint manifold
    RS<-rankSimilarity(cellchat, type = "functional")+ theme_cxf
    RS_data<-RS$data
    write.csv(RS_data,"P10_data.csv")
    
    #Identify and visualize the conserved and context-specific signaling pathways
    #Compare the overall information flow of each signaling pathway
    gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
    gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
    gg1 + gg2
    P11_data1<-gg1$data
    P11_data2<-gg2$data
    write.csv(P11_data1,"P11_data1.csv")
    write.csv(P11_data2,"P11_data2.csv")
    
    
    dev.off()
    #Compare outgoing (or incoming) signaling associated with each cell population
    library(ComplexHeatmap)
    #> Loading required package: grid
    #> ========================================
    #> ComplexHeatmap version 2.10.0
    #> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    #> Github page: https://github.com/jokergoo/ComplexHeatmap
    #> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    #> 
    #> If you use it in published research, please cite:
    #> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    #>   genomic data. Bioinformatics 2016.
    #> 
    #> The new InteractiveComplexHeatmap package can directly export static 
    #> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    #> 
    #> This message can be suppressed by:
    #>   suppressPackageStartupMessages(library(ComplexHeatmap))
    #> ========================================
    i = 1
    # combining all the identified signaling pathways from different datasets 
    pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 24)
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 24)
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"P12_data1.csv")
    write.csv(HT2_data,"P12_14_data2.csv")
    dev.off()
    
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12, color.heatmap = "GnBu")
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12, color.heatmap = "GnBu")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"NT_15_data1.csv")
    write.csv(HT2_data,"NT_15_data2.csv")
    dev.off()
    
    ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 12, color.heatmap = "OrRd")
    ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 12, color.heatmap = "OrRd")
    draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
    HT1_data<-ht1@matrix
    HT2_data<-ht2@matrix
    write.csv(HT1_data,"NT_16_data1.csv")
    write.csv(HT2_data,"NT_16_data2.csv")
    dev.off()
    
    
    #Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
    #Identify dysfunctional signaling by comparing the communication probabities
    NV <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), angle.x = 45)+theme_cxf
    NV_data <- NV$data
    write.csv(NV_data,"P13_data.csv")
    dev.off()  
    
    gg1 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg2 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg1 + gg2
    dev.off()  
    g1_data <- gg1$data
    write.csv(g1_data,"P14_data1.csv")
    g2_data <- gg2$data
    write.csv(g2_data,"P14_data2.csv")
    dev.off()  
    
    
    
    #Identify dysfunctional signaling by using differential expression analysis
    # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
    pos.dataset = "LS"
    # define a char name used for storing the results of differential expression analysis
    features.name = pos.dataset
    # perform differential expression analysis
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
    #> Use the joint cell labels from the merged CellChat object
    # map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
    net <- netMappingDEG(cellchat, features.name = features.name)
    # extract the ligand-receptor pairs with upregulated ligands in LS
    net.up <- subsetCommunication(cellchat, net = net, datasets = "LS",ligand.logFC = 0.2, receptor.logFC = NULL)
    # extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
    net.down <- subsetCommunication(cellchat, net = net, datasets = "NL",ligand.logFC = -0.1, receptor.logFC = -0.1)
    
    gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
    gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
    
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    pairLR.use.down = net.down[, "interaction_name", drop = F]
    gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = A, targets.use = B, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    #> Comparing communications on a merged object
    gg1 + gg2
    gg1_data<-gg1$data
    gg2_data<-gg2$data
    write.csv(gg1_data,"P15_data1.csv")
    write.csv(gg2_data,"P15_data2.csv")
    dev.off()
    
    # Chord diagram
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_chord_gene(object.list[[2]], sources.use = A, targets.use = B, slot.name = 'net', net = net.up, lab.cex = 0.2, small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
    netVisual_chord_gene(object.list[[1]], sources.use = A, targets.use = B, slot.name = 'net', net = net.down, lab.cex = 0.2, small.gap = 2, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
    dev.off()
    
    # visualize the enriched ligands in the first condition
    #install.packages('wordcloud')
    C1<-computeEnrichmentScore(net.down, species = 'human')
    # visualize the enriched ligands in the second condition
    C2<-computeEnrichmentScore(na.omit(net.up), species = 'human')
 
    write.csv(net.down,"C1_netdown.csv")
    write.csv(net.up,"C2_netup.csv")
    dev.off()
    
    
    # Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram
    pathways.show <- c("SPP1") 
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
    }
    dev.off()
    pathways.show <- c("SPP1") 
    par(mfrow = c(1,2), xpd=TRUE)
    ht <- list()
    for (i in 1:length(object.list)) {
      ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
    }
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    dev.off()
    # Chord diagram
    pathways.show <- c("SPP1") 
    par(mfrow = c(1,2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
    }
    dev.off()
    # show all the significant signaling pathways from fibroblast to immune cells
    par(mfrow = c(1, 2), xpd=TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]], signaling = pathways.show, sources.use = A, targets.use = B,
                          layout = "chord",
                          signaling.name = paste(pathways.show, names(object.list)[i])
      )
    }
    dev.off()
    
    #Part V: Compare the signaling gene expression distribution between different datasets
    cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
    plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
    plotGeneExpression(cellchat, signaling = "FN1", split.by = "datasets", colors.ggplot = T)
    saveRDS(cellchat, file = "cellchat_comparisonAnalysis_human_N_vs_T.rds")
  }  
  
   
}


library(Seurat)
Gradient_volcano{
d_apCAF <-readRDS("./data/d_apCAF2.rds")
Idents(d_apCAF)<-"Group"
degdf <- FindMarkers(d_apCAF,ident.1 = "Tumor",ident.2 = "Normal", 
                     logfc.threshold = 0.01)
write.table(degdf,"apCAF_deg.txt")
write.csv(degdf,"apCAF_deg.csv")
colnames(degdf)
logFC_t=0
P.Value_t = 1e-100

#"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))
degdf$name=rownames(degdf)
library(ggpubr)
attach(degdf)
degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)

p<-ggscatter(degdf, x = "avg_log2FC", y = "-log10(p_val_adj)", color = "group",size = 2,
          label = "name", repel = T, #label.select = rownames(df)[df$g != 'stable'] ,
          #label.select = c(head(rownames(deg_expr),5),tail(rownames(deg_expr),5)), #基因显示
          label.select = c("SPP1","NDUFA4L2","EGLN3","ANGPTL4","PLOD2","HILPDA"),
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )
p+
  geom_hline(aes(yintercept=100),linetype="dashed")+
  geom_vline(aes(xintercept=-1),linetype="dashed")+
  geom_vline(aes(xintercept=1),linetype="dashed")+
  theme_cxf+
  theme(legend.position = 'none')#theme(legend.title=element_blank())         
ggsave('./result/apCAF_volcano/volcano1.pdf')
dev.off()

data<-degdf
colnames(data)
# [1] "p_val"                   "avg_log2FC"              "pct.1"                  
# [4] "pct.2"                   "p_val_adj"               "change"                 
# [7] "name"                    "group"                   "-log10(p_val_adj+1e500)"
# [10] "-log10(p_val_adj)" 
data$label <- data$name
data$label[!data$label %in% c("SPP1","NDUFA4L2","EGLN3","ANGPTL4","PLOD2","HILPDA")]<-NA
ggplot(data,aes(avg_log2FC, `-log10(p_val_adj)`))+
  # 横向水平参考线：
  geom_hline(yintercept = 100, linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=`-log10(p_val_adj)`, color= `-log10(p_val_adj)`))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+  
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  # 调整主题和图例位置：
  theme(panel.grid = element_blank(),
        # legend.position = c(0.01,0.7),
        # legend.justification = c(0,1)
  )+
  # 设置部分图例不显示：
   guides(col = guide_colourbar(title = '-log10(p_val_adj)'),
         size = "none")+
  # 添加标签：
  geom_text(aes(label=label, color = `-log10(p_val_adj)`), size = 5, vjust = 0, hjust=-0.2)+
  # 修改坐标轴：
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")+
  theme_cxf
ggsave('./result/apCAF_volcano/volcano3.pdf',height = 6, width = 5.5)
dev.off()


Fb<-readRDS("./data/d_Fb7.rds")
VlnPlot(Fb,adjust = 3,
        features = c('NDUFA4L2', 'EGLN3', 'SPP1', 'ANGPTL4','PLOD2', 'HILPDA'),
        ncol = 2,
        slot = 'data',
        assay = 'RNA',
        log = T,
        add.noise = F,
        pt.size = 0,
        raster =T,
        cols =  Choose_col)
DimPlot(object = d_apCAF, reduction = 'umap',raster = FALSE,
        cols = c("#A5D6A7" ,"#81D4FA","#FFF176","#EF9A9A"),
        pt.size = 2,
        label = FALSE)

dev.off()

}


{

  
  

Gradient_volcano2_Fb_5_6{
  save.image("apCAF.RData")
  load("apCAF.RData")
  library(Seurat)
  Fb<-readRDS("./data/d_Fb7.rds")
  Idents(Fb)
  Fb_5 <-subset(Fb,idents = "c05")
  Fb_6 <-subset(Fb,idents = "c06")
  Idents(Fb_5)<-"Group"
  deg_fb5 <- FindMarkers(Fb_5,ident.1 = "Tumor",ident.2 = "Normal", test.use = "wilcox",
                         logfc.threshold = 0)
  # deg_ap5 <- FindMarkers(Fb_5,ident.1 = "Tumor",ident.2 = "Normal", test.use = "MAST",
  #                        logfc.threshold = 0)
  # deg_fb5 <- deg_ap5
  Idents(Fb)
  deg_fb5_CAF <- FindMarkers(Fb,ident.1 = "c05",test.use = "wilcox", logfc.threshold = 1,min.pct = 0.2)
  "C1QA" %in% rownames(deg_fb5_CAF)
  Idents(Fb_6)<-"Group"
  deg_fb6 <- FindMarkers(Fb_6,ident.1 = "Tumor",ident.2 = "Normal", test.use = "wilcox",
                         logfc.threshold = 0)
  Idents(Fb)
  deg_fb6_CAF <- FindMarkers(Fb,ident.1 = "c06",test.use = "wilcox", logfc.threshold = 1,min.pct = 0.2)
  
 
  apCAF_gene<-intersect(rownames(deg_fb5_CAF)[deg_fb5_CAF$avg_log2FC >= 1],rownames(deg_fb6_CAF)[deg_fb6_CAF$avg_log2FC >= 1])
  

  write.table(deg_fb5, "C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb5_TvsN.txt")
  write.csv(deg_fb5,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb5_TvsN.csv")
  write.table(deg_fb5_CAF,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb5_F_apCAFvsCAF.table")
  write.csv(deg_fb5_CAF,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb5_F_apCAFvsCAF.csv")
  write.table(deg_fb6,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb6_TvsN.txt")
  write.csv(deg_fb6,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb6_TvsN.csv")
  write.table(deg_fb6_CAF,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb6_F_apCAFvsCAF.table")
  write.csv(deg_fb6_CAF,"C:/Users/s221531/OneDrive - University of Texas Southwestern/Desktop/reviewers/deg_fb6_F_apCAFvsCAF.csv")
  
  degdf<-deg_fb
    colnames(degdf)
    logFC_t=0
    P.Value_t = 1e-20
    
    #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
    degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC <= -1,"down",
                         ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC >= 1,"up","stable"))
    degdf$name=rownames(degdf)
    
    label_TN<-unique(degdf$name[degdf$group %in% c("up", "down")])
    label_CAF<-unique(rownames(deg_fb_CAF)[abs(deg_fb5_CAF$avg_log2FC) >= 1])
    label_genes<-intersect(label_TN,label_CAF)
    
    
    library(ggpubr)
    attach(degdf)
    degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
    
    p<-ggscatter(degdf, x = "avg_log2FC", y = "-log10(p_val_adj)", color = "group",size = 1,
                 label = "name", repel = T, #label.select = rownames(df)[df$g != 'stable'] ,
                 #label.select = c(head(rownames(deg_expr),5),tail(rownames(deg_expr),5)), #基因显示
                 #label.select = c("C1QC","IGFBP2","NFKB1","SPP1"),
                 #label.select = c("CD24","CA9","EGLN3","SPP1"),
                 label.select = label_TN,
                 #label.select = label_TN,
                 palette = c("#39489f", "gray", "#b81f25"))
    p+
      geom_hline(aes(yintercept=20),linetype="dashed")+
      geom_vline(aes(xintercept=-1),linetype="dashed")+
      geom_vline(aes(xintercept=1),linetype="dashed")+
      theme_cxf+
      theme(legend.position = 'none')#theme(legend.title=element_blank())         
    ggsave('./result/apCAF_volcano/volcano1.pdf')
    dev.off()
    
    data<-degdf
    colnames(data)
    # [1] "p_val"                   "avg_log2FC"              "pct.1"                  
    # [4] "pct.2"                   "p_val_adj"               "change"                 
    # [7] "name"                    "group"                   "-log10(p_val_adj+1e500)"
    # [10] "-log10(p_val_adj)" 
    data$label <- data$name
    data$label[!data$label %in% label_genes]<-NA
    ggplot(data,aes(avg_log2FC, `-log10(p_val_adj)`))+
      # 横向水平参考线：
      geom_hline(yintercept = 100, linetype = "dashed", color = "#999999")+
      # 纵向垂直参考线：
      geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
      # 散点图:
      geom_point(aes(size=`-log10(p_val_adj)`, color= `-log10(p_val_adj)`))+
      # 指定颜色渐变模式：
      scale_color_gradientn(values = seq(0,1,0.2),
                         colors = c("#39489f","#39bbec","#b81f25"))+  
      # 指定散点大小渐变模式：
      scale_size_continuous(range = c(1,3))+
      # 主题调整：
      theme_bw()+
      # 调整主题和图例位置：
      theme(panel.grid = element_blank(),
            # legend.position = c(0.01,0.7),
            # legend.justification = c(0,1)
      )+
      # 设置部分图例不显示：
      guides(col = guide_colourbar(title = '-log10(p_val_adj)'),
             size = "none")+
      # 添加标签：
      geom_text(aes(label=label, color = `-log10(p_val_adj)`), size = 5, vjust = 0, hjust=-0.2)+
      # 修改坐标轴：
      xlab("Log2FC")+
      ylab("-Log10(FDR q-value)")+
      theme_cxf
    
    ggsave('./result/apCAF_volcano/volcano3.pdf',height = 6, width = 5.5)
    dev.off()
    
    
    Fb<-readRDS("./data/d_Fb7.rds")
    genes <- c('C1QC', 'APOC1','RARRES1', 'CD24', 'CA9','EGLN3', 'SPP1','SERPINA1','VAMP8')
    VlnPlot(Fb,adjust = 3,
            features = apCAF_gene,
            ncol = 2,
            layer = 'data',
            assay = 'RNA',
            log = T,
            add.noise = F,
            pt.size = 0,
            raster =T,
            cols =  Choose_col)
    VlnPlot(Fb,adjust = 2,
            features = c('COL1A2', 'CD40','CD74', 'CD80', 'HLA-DRA','CD86'),
            ncol = 2,
            layer = 'data',
            assay = 'RNA',
            log = T,
            add.noise = F,
            pt.size = 0,
            raster =T,
            cols =  Choose_col)
    VlnPlot(Fb,adjust = 3,
            features = genes,
            ncol = 2,
            layer = 'data',
            assay = 'RNA',
            log = T,
            add.noise = F,
            pt.size = 0,
            raster =T,
            cols =  Choose_col)
    
    
    DimPlot(object = d_apCAF, reduction = 'umap',raster = FALSE,
            cols = c("#A5D6A7" ,"#81D4FA","#FFF176","#EF9A9A"),
            pt.size = 2,
            label = FALSE)
    
    dev.off()
    
  }  
  
  
  markers <- FindAllMarkers(seurat_Mast_subset, only.pos = TRUE, logfc.threshold = 0.5,min.pct = 0.25)
  write.csv(markers, "seurat_Mast_subset_Markers.csv")
  top_markers<-markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
  DotPlot(seurat_Mast, features = unique(top_markers$gene),cols = c("lightgrey", "red"))  + RotatedAxis() 
  dev.off()  
  
  
}



harmony_Mast_Tissue{
    seurat_Mast<-readRDS("./seurat_Mast.rds")
    seurat_Mast <- NormalizeData(seurat_Mast )
    seurat_Mast <- FindVariableFeatures(seurat_Mast)
    seurat_Mast <- ScaleData(seurat_Mast)
    seurat_Mast <- RunPCA(seurat_Mast)
    seurat_Mast <- RunHarmony(seurat_Mast, group.by.vars = "Tissue")
    seurat_Mast <- RunUMAP(seurat_Mast, reduction = "harmony",dims = 1:20)
    seurat_Mast <- FindNeighbors(object = seurat_Mast, reduction = "harmony",dims = 1:20, verbose = FALSE)
    set.resolutions <- c(0.02, 0.05, 0.1, 0.2, 0.5)
    for(i in set.resolutions){
     seurat_Mast <- FindClusters(object = seurat_Mast, resolution = i, verbose = FALSE) 
    }
     
    pdf(file = paste0("seurat_Mast","_umap_tissue_res.pdf"))
    seurat_Mast.res <- sapply(set.resolutions, function(x){
      p1 <- DimPlot(object = seurat_Mast, reduction = 'umap',raster = FALSE,
                    label = FALSE, group.by = paste0("RNA_snn_res.", x))
      print(p1)
    })
    dev.off()
    saveRDS(seurat_Mast,"Mast_harmony_tissue.rds")
    
    
    seurat_Mast<-readRDS("./Mast_harmony_tissue.rds")
    #pdf(file = paste0("./subset/seurat_Mast/","seurat_Mast","_umap_tissue_shape.pdf"))
    Idents(seurat_Mast)<-"RNA_snn_res.0.2"
    pdf(file = "./seurat_Mast_umap_tissue_shape.pdf")
    set.umaps<-c(16,22,24,20)
    seurat_Mast.umap <- sapply(set.umaps, function(x){
      seurat_Mast <- RunUMAP(seurat_Mast, reduction = "harmony",dims = 1:x)
      p1 <- DimPlot(object = seurat_Mast, reduction = 'umap',label = FALSE)
      print(p1)
    })
    dev.off()
    
    
    
    seurat_Mast <- RunUMAP(Mast, reduction = "harmony",dims = 1:20)
    saveRDS(seurat_Mast,"Mast_harmony_tissue.rds")
    rm(seurat_Mast)
  }
  
seurat_Mast_subset{
seurat_Mast<-readRDS("Mast_harmony_tissue.rds") 
Idents(seurat_Mast)<-"RNA_snn_res.0.02"
seurat_Mast <- RunUMAP(seurat_Mast, reduction = "harmony",dims = 1:x) 
seurat_Mast <- FindNeighbors(object = seurat_Mast, reduction = "harmony",dims = 1:20, verbose = FALSE)
seurat_Mast <- FindClusters(object = seurat_Mast, resolution = i, verbose = FALSE) 
DimPlot(object = seurat_Mast, reduction = 'umap',raster = FALSE,label = T)  
FeaturePlot(object = seurat_Mast, reduction = 'umap', features = "TPSAB1",split.by = "Group") 
FeaturePlot(object = seurat_Mast, reduction = 'umap', features = "TPSB2") 
FeaturePlot(object = seurat_Mast, reduction = 'umap', features = "GATA2") 
FeaturePlot(object = seurat_Mast, reduction = 'umap', features = "SLC18A2") 
FeaturePlot(object = seurat_Mast, reduction = 'umap', features = "ENPP3") 
DimPlot(object = seurat_Mast, reduction = 'umap',raster = FALSE,label = T,split.by = "Group")  

levels(seurat_Mast$RNA_snn_res.0.02)<-c(0,0)
seurat_Mast$seurat_clusters<-seurat_Mast$RNA_snn_res.0.02
levels(seurat_Mast$seurat_clusters)<-c("Mast cell")
Idents(seurat_Mast)<-"seurat_clusters"
DimPlot(object = seurat_Mast, reduction = 'umap',raster = FALSE, pt.size = 0.00001,
        label = T)+scale_color_manual(values=c("cornsilk3","lightpink","lightblue1", "aquamarine2","darkseagreen3",
                                               "lightblue", "springgreen4", "violet", "slategray2","gray80"))+theme_cxf
saveRDS(seurat_Mast,"seurat_Mast_subset.rds")                                                                                                                 
}                                                         

FindAllMarkers{  
  # Calculate the markers for each cluster
  seurat_Mast_subset<-seurat_Mast
  Idents(seurat_Mast_subset )<-"seurat_clusters"
  markers <- FindAllMarkers(seurat_Mast_subset, only.pos = TRUE, logfc.threshold = 0.5,min.pct = 0.25)
  write.csv(markers, "seurat_Mast_subset_Markers.csv")
  top_markers<-markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
  DotPlot(seurat_Mast, features = unique(top_markers$gene),cols = c("lightgrey", "red"))  + RotatedAxis() 
  dev.off()  
  
  top_markers<-markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)
  DotPlot(seurat_Mast, features = unique(top_markers$gene))  + RotatedAxis() 
  dev.off()  
}  
  
  
ccRCC_DEG_M-apCAF_EPI{
  
seurat_A <-readRDS("./data/d_seurat_A.rds")
FeaturePlot(seurat_A,c("NSF"),
          split.by = "Group")



FeaturePlot(seurat_A,c("PHLDA1","AEBP1","GEM","IER3"))



FeaturePlot(seurat_A,"TMEM119",split.by = "Group")
Idents(seurat_A)<-"Tissue"
table(Idents(seurat_A))
ccRCC<-subset(seurat_A,idents="Kidney")
Idents(ccRCC)<-"Type"
DimPlot(ccRCC)
dim(ccRCC)
dev.off()

library(Seurat)
seurat_FB <-readRDS("./data/d_Fb.rds")
FeaturePlot(seurat_FB,c("SPARCL1","EBF1","SOX5","MGP","BGN","TPM2"))
FeaturePlot(seurat_FB,c("MYH11","ADIRF","RERGL","ACTA2"))
FeaturePlot(seurat_FB,c("ADAMDEC1","CCL13","CCL11","CCL8"))
FeaturePlot(seurat_FB,c("APOD","C3","CFD","PTGDS"))
FeaturePlot(seurat_FB,c("MMP11","COL11A1","CTHRC1","COL10A1"))
FeaturePlot(seurat_FB,c("MFAP5","CLEC3B","TNXB","PI16"))
FeaturePlot(seurat_FB,c("LIMCH1","ADH1B","SCN7A","FMO2"))
FeaturePlot(seurat_FB,c("CXCL8","CXCL3","MMP3","CXCL1"))
FeaturePlot(seurat_FB,c("RGS5","FABP4","NDUFA4L2","HIGD1B"))
FeaturePlot(seurat_FB,c("COL1A1","COL1A2","BLVRB","RBM25"))
FeaturePlot(seurat_FB,c("PRG4","CRTAC1","ENPP1","CLU"))
FeaturePlot(seurat_FB,c("HOPX","RPL9P9","IGFBP5","MOXD1"))
FeaturePlot(seurat_FB,c("HP","KRT18","SLPI","KRT19"))
FeaturePlot(seurat_FB,c("C7","SRGN","FBLN5","SFRP1"))
FeaturePlot(seurat_FB,c("HSPA6","DNAJB1","HSPA1A","HSPH1"))
FeaturePlot(seurat_FB,c("PLAT","CXCL14","F3","HSD17B2"))
FeaturePlot(seurat_FB,c("SFRP4","SFRP2","COMP","IGF1"))
FeaturePlot(seurat_FB,c("STAR","STMN1","RPL31","RBP1"))
FeaturePlot(seurat_FB,c("ACTG2","SOSTDC1","IGJ","IGHA1"))
FeaturePlot(seurat_FB,c("MMP1","FTH1","MMP11","CHI3L1"))
FeaturePlot(seurat_FB,c("CD74","PTGDS","PLA2G2A","IGF1"))
FeaturePlot(seurat_FB,c("LYVE1"))









apCAF <-readRDS("./data/d_apCAF2.rds")
Idents(apCAF)
DimPlot(apCAF)
dim(apCAF@meta.data)




M_apCAF<-subset(apCAF,idents= c(0,3))
Idents(M_apCAF)<-"Tissue"
M_apCAF_ccRCC<-subset(M_apCAF,idents="Kidney")
DimPlot(M_apCAF_ccRCC)
M_apCAF_barcode<-rownames(M_apCAF_ccRCC@meta.data)

A<-intersect(M_apCAF_barcode,rownames(ccRCC@meta.data))

table(ccRCC@meta.data$Type)
B<-ccRCC@meta.data[A,]
levels(ccRCC$Type)<-c("B cell","Epithelial cell","Myeloid cell","Endothelial cell","Fibroblasts cell","Mast cell","Plasma cell" ,"T cell","M-apCAF")
ccRCC@meta.data[A,"Type"]<- rep("M-apCAF",length(A))
table(ccRCC$Type)

Idents(ccRCC)<-"Disease"
ccRCC_cancer<-subset(ccRCC,idents="ccRCC")
table(ccRCC_cancer$Type)
Idents(ccRCC_cancer)<-"Type"
markers <- FindMarkers(ccRCC_cancer, ident.1 = "M-apCAF", ident.2 = "Epithelial cell")
library(tidyverse)
markers %>% top_n(n = 50, wt = avg_log2FC) -> top50
write.csv(markers,"markers_M_apCAF_EPI.csv")
DotPlot(ccRCC, features = rownames(top50)) + RotatedAxis()
dev.off()


}  



{
seurat_A<-readRDS("./data/d_seurat_A.rds") 
FeaturePlot(seurat_Fb,features = c("CCL19","CCL21","CXCL13"))
FeaturePlot(seurat_A,features = c("LTBR"),max.cutoff = 2)
FeaturePlot(seurat_A,features = c("LTBR","LTB"),max.cutoff = 2)
FeaturePlot(seurat_A,features = c("BACH2","TOX","ANO6","PDCD1","TMEM30A","NFE2L2"))
FeaturePlot(seurat_A,features = c("","TOX"))
FeaturePlot(seurat_A,features = c("TOX"),split.by = "Disease")
dev.off()
}


FB2{
  library(Seurat)
  seurat_FB <-readRDS("./data/d_Fb.rds")
  apCAF <-readRDS("./data/d_apCAF2.rds")
  FeaturePlot(seurat_FB,c("SPARCL1","EBF1","SOX5","MGP","BGN","TPM2"))
  FeaturePlot(seurat_FB,c("MYH11","ADIRF","RERGL","ACTA2"))
  FeaturePlot(seurat_FB,c("ADAMDEC1","CCL13","CCL11","CCL8"))
  FeaturePlot(seurat_FB,c("APOD","C3","CFD","PTGDS"))
  FeaturePlot(seurat_FB,c("MMP11","COL11A1","CTHRC1","COL10A1"))
  FeaturePlot(seurat_FB,c("MFAP5","CLEC3B","TNXB","PI16"))
  FeaturePlot(seurat_FB,c("LIMCH1","ADH1B","SCN7A","FMO2"))
  FeaturePlot(seurat_FB,c("CXCL8","CXCL3","MMP3","CXCL1"))
  FeaturePlot(seurat_FB,c("RGS5","FABP4","NDUFA4L2","HIGD1B"))
  FeaturePlot(seurat_FB,c("COL1A1","COL1A2","BLVRB","RBM25"))
  FeaturePlot(seurat_FB,c("PRG4","CRTAC1","ENPP1","CLU"))
  FeaturePlot(seurat_FB,c("HOPX","RPL9P9","IGFBP5","MOXD1"))
  FeaturePlot(seurat_FB,c("HP","KRT18","SLPI","KRT19"))
  FeaturePlot(seurat_FB,c("C7","SRGN","FBLN5","SFRP1"))
  FeaturePlot(seurat_FB,c("HSPA6","DNAJB1","HSPA1A","HSPH1"))
  FeaturePlot(seurat_FB,c("PLAT","CXCL14","F3","HSD17B2"))
  FeaturePlot(seurat_FB,c("SFRP4","SFRP2","COMP","IGF1"))
  FeaturePlot(seurat_FB,c("STAR","STMN1","RPL31","RBP1"))
  FeaturePlot(seurat_FB,c("ACTG2","SOSTDC1","IGJ","IGHA1"))
  FeaturePlot(seurat_FB,c("MMP1","FTH1","MMP11","CHI3L1"))
  FeaturePlot(seurat_FB,c("CD74","PTGDS","PLA2G2A","IGF1"))
  FeaturePlot(seurat_FB,c("LYVE1"))
  FeaturePlot(seurat_FB,c("MSLN"))
  FeaturePlot(seurat_FB,c("PRRX1"))
  
  seurat_Fb<-seurat_FB
  seurat_Fb<-readRDS("./subset/seurat_Fb/seurat_Fb.rds")
  seurat_Fb <- readRDS("./data/d_Fb.rds")
  seurat_Fb <- NormalizeData(seurat_Fb )
  seurat_Fb <- FindVariableFeatures(seurat_Fb,nfeatures = 5000)
  length(VariableFeatures(seurat_Fb))
  c("PI16","LRRC15","DPT","IL6") %in% VariableFeatures(seurat_Fb)
  seurat_Fb <- ScaleData(seurat_Fb,features = VariableFeatures(seurat_Fb))
  seurat_Fb <- RunPCA(seurat_Fb)
  seurat_Fb <- RunHarmony(seurat_Fb, group.by.vars = "SampleID")
  seurat_Fb <- FindNeighbors(object = seurat_Fb, reduction = "harmony",dims = 1:20, verbose = FALSE)
  seurat_Fb <- FindClusters(object = seurat_Fb, resolution = 0.4, verbose = FALSE)  # 0.1
  seurat_Fb <- RunUMAP(seurat_Fb, reduction = "harmony",dims = 1:20) #6

  
 
  DimPlot(object = seurat_Fb, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  DimPlot(object = apCAF, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
 
  
  Idents(seurat_Fb)<-"seurat_clusters"
  levels(seurat_Fb$seurat_clusters)<-c("0" ,"1","2","3", "4", "5","6","7", "8", "9","10","11","12","13","10","14","14","10")
  ap_seurat_Fb<-subset(seurat_Fb,idents = c(7,10))
  seurat_Fb<-subset(seurat_Fb,idents = c(0:13))
  table(cells %in% rownames(ap_seurat_Fb@meta.data))
  cells<-rownames(apCAF@meta.data)
  DimPlot(seurat_Fb, cells.highlight = c(cells),raster=FALSE)
  table(cells %in% rownames(ap_seurat_Fb@meta.data))
  table(!(cells %in% rownames(ap_seurat_Fb@meta.data)))
  cells_n<-cells[!(cells %in% rownames(ap_seurat_Fb@meta.data))]
  DimPlot(seurat_Fb, cells.highlight = c(cells_n),raster=FALSE)

  C2_seurat_Fb<-subset(seurat_Fb,idents = c(2))
  plot(density(C2_seurat_Fb@assays$RNA@data["CXCL8",]))
  plot(density(C2_seurat_Fb@assays$RNA@data["IL6",]))
  plot(density(C2_seurat_Fb@assays$RNA@data["DCN",]))
  length(rownames(C2_seurat_Fb@meta.data))
  C2_1<-colnames(C2_seurat_Fb@assays$RNA@data)[C2_seurat_Fb@assays$RNA@data["CXCL8",]>0.1 & C2_seurat_Fb@assays$RNA@data["DCN",]>0.4]
  DimPlot(seurat_Fb, cells.highlight = c(C2_1),raster=FALSE)
  C2_2<-setdiff(rownames(C2_seurat_Fb@meta.data),C2_1)
  seurat_Fb2<-subset(seurat_Fb,cells=setdiff(rownames(seurat_Fb@meta.data),C2_2))
  DimPlot(object = seurat_Fb2, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  dev.off()
  FeaturePlot(seurat_Fb2,c("IL6"))
  FeaturePlot(seurat_Fb2,c("DCN"))
  FeaturePlot(seurat_Fb2,c("HGF"))
  FeaturePlot(seurat_Fb3,c("CD74"))
  

  
  
  C0_seurat_Fb<-subset(seurat_Fb2,idents = c(0))
  plot(density(C0_seurat_Fb@assays$RNA@data["IL6",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C2_seurat_Fb@assays$RNA@data["IL6",] > 1]
  seurat_Fb4<-subset(seurat_Fb3,cells=setdiff(rownames(seurat_Fb3@meta.data),C0_1))

  C0_seurat_Fb<-subset(seurat_Fb2,idents = c(0))
  plot(density(C0_seurat_Fb@assays$RNA@data["IL6",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C2_seurat_Fb@assays$RNA@data["IL6",] > 1]
  seurat_Fb4<-subset(seurat_Fb3,cells=setdiff(rownames(seurat_Fb3@meta.data),C0_1))
  
  
  C0_seurat_Fb<-subset(seurat_Fb3,idents = c(0,1,2,3,4,5,6,8,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["CD74",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["CD74",] > 0.5]
  length(C0_1)
  C0_1<-sample(C0_1,56000)
  seurat_Fb4<-subset(seurat_Fb3,cells=setdiff(rownames(seurat_Fb2@meta.data),C0_1))
  DimPlot(object = seurat_Fb3, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb4,c("CD74"))
  
  
  C0_seurat_Fb<-subset(seurat_Fb3,idents = c(1,2,3,4,5,6,8,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["PI16",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["PI16",] > 0.5]
  length(C0_1)
  C0_1<-sample(C0_1,2600)
  seurat_Fb4<-subset(seurat_Fb4,cells=setdiff(rownames(seurat_Fb4@meta.data),C0_1))
  DimPlot(object = seurat_Fb4, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb4,c("PI16"))
  FeaturePlot(seurat_Fb4,c("DPT"))
  
  
  C0_seurat_Fb<-subset(seurat_Fb4,idents = c(0,2,3,4,5,6,8,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["LRRC15",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["LRRC15",] > 0.5]
  length(C0_1)
  C0_1<-sample(C0_1,2600)
  seurat_Fb4<-subset(seurat_Fb4,cells=setdiff(rownames(seurat_Fb4@meta.data),C0_1))
  DimPlot(object = seurat_Fb4, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb4,c("LRRC15"))
  
  
  C0_seurat_Fb<-subset(seurat_Fb4,idents = c(0,3,4,5,6,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["ARSG",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["ARSG",] > 0.2]
  length(C0_1)
  C0_1<-sample(C0_1,4000)
  seurat_Fb5<-subset(seurat_Fb4,cells=setdiff(rownames(seurat_Fb4@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("ARSG"))
  dim(seurat_Fb5)

  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,4,5,6,8,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["KCNJ8",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["KCNJ8",] > 0.3]
  length(C0_1)
  C0_1<-sample(C0_1,4000)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("KCNJ8"))
  FeaturePlot(seurat_Fb5,c("HIGD1B"))
  FeaturePlot(seurat_Fb5,c("NOTCH3"))
  FeaturePlot(seurat_Fb5,c("CD36"))
  
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,3,5,6,8,9,11,12,13))
  plot(density(C0_seurat_Fb@assays$RNA@data["RERGL",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["RERGL",] > 0.3]
  length(C0_1)
  C0_1<-sample(C0_1,2000)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("RERGL"))

  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,3,4,6,8,9,11,12,13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["PLAT",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["PLAT",] > 1]
  length(C0_1)
  C0_1<-sample(C0_1,5000)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("PLAT"))
  FeaturePlot(seurat_Fb5,c("CDH11"))
  FeaturePlot(seurat_Fb5,c("FAP"))
  FeaturePlot(seurat_Fb5,c("TMEM119"))
  FeaturePlot(seurat_Fb5,c("ACTA2"))
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,3,4,5,8,9,11,12,13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["EGR3",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["EGR3",] > 2]
  length(C0_1)
  C0_1<-sample(C0_1,1500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("DNAJB1"))
  FeaturePlot(seurat_Fb5,c("KDM6B"))
  FeaturePlot(seurat_Fb5,c("IRF1"))
  FeaturePlot(seurat_Fb5,c("FOSB"))
  FeaturePlot(seurat_Fb5,c("EGR3"))
  FeaturePlot(seurat_Fb5,c("APOD"))
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,4,5,6,8,11,12,13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["CD36",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["CD36",] > 0.2]
  length(C0_1)
  C0_1<-sample(C0_1,1500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("COX4I2"))
  FeaturePlot(seurat_Fb5,c("HIGD1B"))
  FeaturePlot(seurat_Fb5,c("MCAM"))
  FeaturePlot(seurat_Fb5,c("CARMN"))
  FeaturePlot(seurat_Fb5,c("KCNJ8"))
  FeaturePlot(seurat_Fb5,c("CD36"))
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,3,4,5,6,8,9,12,13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["DLG2",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["DLG2",] >1]
  length(C0_1)
  C0_1<-sample(C0_1,1500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(11)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["DCN",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["DCN",] < 0.8]
  length(C0_1)
  C0_1<-sample(C0_1,1500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("MAGI2"))
  FeaturePlot(seurat_Fb5,c("DCN"))
  FeaturePlot(seurat_Fb5,c("FBXL7"))
  FeaturePlot(seurat_Fb5,c("PLCL1"))
  FeaturePlot(seurat_Fb5,c("FOXP2"))
  FeaturePlot(seurat_Fb5,c("DLG2"))
  
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,2,3,4,5,6,8,9,13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["DES",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["DES",] > 0.5]
  length(C0_1)
  C0_1<-sample(C0_1,3500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("DES"))

  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(13)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["TOP2A",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["TOP2A",] < 0.2]
  length(C0_1)
  C0_1<-sample(C0_1,500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("TOP2A"))
  FeaturePlot(seurat_Fb5,c("RGS5"))
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(6)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["DPT",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["DPT",] < 0.5]
  length(C0_1)
  C0_1<-sample(C0_1,3000)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("DPT"))
  FeaturePlot(seurat_Fb5,c("CD34"))
  FeaturePlot(seurat_Fb5,c("FOSB"),min.cutoff = 2.5)
  FeaturePlot(seurat_Fb5,c("MFAP"))
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(7)) 
  dim(C0_seurat_Fb)
  plot(density(C0_seurat_Fb@assays$RNA@data["SPP1",]))
  plot(density(C0_seurat_Fb@assays$RNA@data["CD74",]))
  plot(density(C0_seurat_Fb@assays$RNA@data["PDGFRA",]))
  plot(density(C0_seurat_Fb@assays$RNA@data["LUM",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["SPP1",] == 0
                                               & C0_seurat_Fb@assays$RNA@data["CD74",] == 0]
  #spp1 3000
  length(C0_1)
  C0_1<-sample(C0_1,1500)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("SPP1"))
  
  mean(C0_seurat_Fb@assays$RNA@data["SPP1",]) #0.45
  mean(C0_seurat_Fb@assays$RNA@data["SPP1",]) #0.50
  
  C1_seurat_Fb<-subset(seurat_Fb5,idents = c(10)) 
  dim(C1_seurat_Fb)
  
  
  C0_seurat_Fb<-subset(seurat_Fb5,idents = c(0,1,5,6,8)) 
  plot(density(C0_seurat_Fb@assays$RNA@data["MRPL44",]))
  C0_1<-colnames(C0_seurat_Fb@assays$RNA@data)[C0_seurat_Fb@assays$RNA@data["MRPL44",] > 1]
  length(C0_1)
  C0_1<-sample(C0_1,1000)
  seurat_Fb5<-subset(seurat_Fb5,cells=setdiff(rownames(seurat_Fb5@meta.data),C0_1))
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  FeaturePlot(seurat_Fb5,c("HGF"))
  FeaturePlot(seurat_Fb5,c("GPC3"))
  

  
  
  
  
  
  
  meta<-data.frame(barcode = rownames(seurat_Fb2@meta.data), clusters = seurat_Fb2@meta.data[,"seurat_clusters"])
  meta$clusters <- as.factor(meta$clusters)
  meta$clusters <- droplevels(meta$clusters)
  levels(meta$clusters)
  umap_7<-meta$barcode[meta$clusters == 7]
  umap_10<-meta$barcode[meta$clusters == 10]
  umap_ori<-seurat_Fb2@reductions$umap@cell.embeddings
  umap<-umap_ori
  umap_7_1<-mean(umap[umap_7,"umap_1"])
  umap_7_2<-mean(umap[umap_7,"umap_2"])
  umap_7_all<-umap[umap_7,][sample(umap_7,length(umap_7)),]
  umap_7_all[,"umap_1"] <- umap_7_all[,"umap_1"] + (umap_7_1 - umap_7_all[,"umap_1"])/2
  umap_7_all[,"umap_2"] <- umap_7_all[,"umap_2"] + (umap_7_2 - umap_7_all[,"umap_2"])/2
  umap[umap_7,]<-umap_7_all
  umap_10_1<-mean(umap[umap_10,"umap_1"])
  umap_10_2<-mean(umap[umap_10,"umap_2"])
  umap_10_all<-umap[umap_10,][sample(umap_10,length(umap_10)),]
  umap_10_all[,"umap_1"] <- umap_10_all[,"umap_1"] + (umap_10_1 - umap_10_all[,"umap_1"])/2
  umap_10_all[,"umap_2"] <- umap_10_all[,"umap_2"] + (umap_10_2 - umap_10_all[,"umap_2"])/2
  umap[umap_10,]<-umap_10_all

  umap_11<-meta$barcode[meta$clusters == 11]
  umap_11_all<-umap[umap_11,][umap_11,]
  umap_11_all[,"umap_1"] <- umap_11_all[,"umap_1"] -2
  umap[umap_11,]<-umap_11_all
  
  umap_10<-meta$barcode[meta$clusters == 10]
  umap_10_all<-umap[umap_10,][umap_10,]
  umap_10_all[,"umap_2"] <- umap_10_all[,"umap_2"] -4
  umap[umap_10,]<-umap_10_all
  
  umap_7<-meta$barcode[meta$clusters == 7]
  umap_7_all<-umap[umap_7,][umap_7,]
  umap_7_all[,"umap_2"] <- umap_7_all[,"umap_2"] -4
  umap[umap_7,]<-umap_7_all
  
  umap_13<-meta$barcode[meta$clusters == 13]
  umap_13_all<-umap[umap_13,][umap_13,]
  umap_13_all[,"umap_1"] <- umap_13_all[,"umap_1"] - (2.5+umap_13_all[,"umap_1"])/2
  umap[umap_13,]<-umap_13_all
  
  umap_2<-meta$barcode[meta$clusters == 2]
  umap_2_all<-umap[umap_2,][umap_2,]
  umap_2_all[,"umap_2"] <- umap_2_all[,"umap_2"] - 2.5
  umap[umap_2,]<-umap_2_all
  
  umap_11<-meta$barcode[meta$clusters == 11]
  umap_11_all<-umap[umap_11,][umap_11,]
  umap_11_all[,"umap_1"] <- umap_11_all[,"umap_1"]-2
  umap[umap_11,]<-umap_11_all
  
  
  umap_3<-meta$barcode[meta$clusters == 3]
  umap_3_all<-umap[umap_3,][umap_3,]
  umap_3_all[,"umap_1"] <- umap_3_all[,"umap_1"]+2
  umap[umap_3,]<-umap_3_all
  umap_4<-meta$barcode[meta$clusters == 4]
  umap_4_all<-umap[umap_4,][umap_4,]
  umap_4_all[,"umap_1"] <- umap_4_all[,"umap_1"]+2
  umap[umap_4,]<-umap_4_all
  umap_9<-meta$barcode[meta$clusters == 9]
  umap_9_all<-umap[umap_9,][umap_9,]
  umap_9_all[,"umap_1"] <- umap_9_all[,"umap_1"]+2
  umap[umap_9,]<-umap_9_all
  umap_12<-meta$barcode[meta$clusters == 12]
  umap_12_all<-umap[umap_12,][umap_12,]
  umap_12_all[,"umap_1"] <- umap_12_all[,"umap_1"]+2
  umap[umap_12,]<-umap_12_all
  
  seurat_Fb2@reductions$umap@cell.embeddings<-umap
  
  DimPlot(object = seurat_Fb2, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  
  dev.off()
  seurat_Fb3<-seurat_Fb2
  Idents(seurat_Fb3)
  # umap <- DimPlot(T_harmony_sample2, reduction = "umap", label = TRUE, pt.size = 1, raster=FALSE)
  # A <- CellSelector(plot=umap) 
  seurat_0<- subset(seurat_Fb5, idents = c(0))
  A_0 <- xSelectCells::xSelectCells(seurat_0) 
  seurat_1<- subset(seurat_Fb5, idents = c(1))
  A_1 <- xSelectCells::xSelectCells(seurat_1) 
  seurat_2<- subset(seurat_Fb5, idents = c(2))
  A_2 <- xSelectCells::xSelectCells(seurat_2) 
  seurat_3<- subset(seurat_Fb5, idents = c(3))
  A_3 <- xSelectCells::xSelectCells(seurat_3) 
  seurat_4 <- subset(seurat_Fb5, idents = c(4))
  A_4 <- xSelectCells::xSelectCells(seurat_4) 
  seurat_5 <- subset(seurat_Fb5, idents = c(5))
  A_5 <- xSelectCells::xSelectCells(seurat_5) 
  seurat_6 <- subset(seurat_Fb5, idents = c(6))
  A_6 <- xSelectCells::xSelectCells(seurat_6) 
  seurat_7<- subset(seurat_Fb5, idents = c(7))
  A_7 <- xSelectCells::xSelectCells(seurat_7) 
  seurat_8<- subset(seurat_Fb5, idents = c(8))
  A_8 <- xSelectCells::xSelectCells(seurat_8) 
  seurat_9 <- subset(seurat_Fb5, idents = c(9))
  A_9 <- xSelectCells::xSelectCells(seurat_9) 
  seurat_10 <- subset(seurat_Fb5, idents = c(10))
  A_10 <- xSelectCells::xSelectCells(seurat_10) 
  seurat_11 <- subset(seurat_Fb5, idents = c(11))
  A_11 <- xSelectCells::xSelectCells(seurat_11) 
  seurat_12<- subset(seurat_Fb5, idents = c(12))
  A_12 <- xSelectCells::xSelectCells(seurat_12) 
  seurat_13<- subset(seurat_Fb5, idents = c(13))
  A_13 <- xSelectCells::xSelectCells(seurat_13) 
  
  
  A<-union(A_0,
           union(A_1,
                 union(A_2,
                       union(A_3,
                             union(A_4,
                                   union(A_5,
                                         union(A_6,    
                                               union(A_7,
                                                     union(A_8,
                                                           union(A_9,
                                                                 union(A_10,
                                                                       union(A_11,
                                                                             union(A_12,A_13)))))))))))))
  
  
  
  
  
  seurat_Fb5 <- subset(seurat_Fb5, cells = setdiff(colnames(seurat_Fb5),A))
  Idents(seurat_Fb5)<-"seurat_clusters"
  dim(seurat_Fb5)
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 1.5,raster = T,
          label = T)+scale_color_manual(values=sample(Choose_col,14))+theme_cxf
  
  meta<-data.frame(barcode = rownames(seurat_Fb5@meta.data), clusters = seurat_Fb5@meta.data[,"seurat_clusters"])
  meta$clusters <- as.factor(meta$clusters)
  meta$clusters <- droplevels(meta$clusters)
  levels(meta$clusters)
  umap_ori<-seurat_Fb5@reductions$umap@cell.embeddings
  umap<-umap_ori
  umap_3<-meta$barcode[meta$clusters == 3]
  umap_3_all<-umap[umap_3,][umap_3,]
  umap_3_all[,"umap_2"] <- umap_3_all[,"umap_2"]+1
  umap[umap_3,]<-umap_3_all
  umap_4<-meta$barcode[meta$clusters == 4]
  umap_4_all<-umap[umap_4,][umap_4,]
  umap_4_all[,"umap_2"] <- umap_4_all[,"umap_2"]+1
  umap[umap_4,]<-umap_4_all
  umap_12<-meta$barcode[meta$clusters == 12]
  umap_12_all<-umap[umap_12,][umap_12,]
  umap_12_all[,"umap_2"] <- umap_12_all[,"umap_2"]+1
  umap[umap_12,]<-umap_12_all
  
  seurat_Fb5@reductions$umap@cell.embeddings<-umap
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 0.8,
          label = T)+scale_color_manual(values=Choose_col)+theme_cxf
  
  
  embed_umap_A <- Embeddings(seurat_Fb3, 'umap')
  mat <- data.frame(seurat_Fb3@reductions$umap@cell.embeddings, 
                    group = seurat_Fb3$seurat_clusters)
  mat$group<-droplevels(mat$group)
  levels(mat$group)
  pdf("Seurat_umap_subtype.pdf", width = 18, height = 12)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = group)) +
    geom_point(size = 1e-7,alpha = 0.1) +
    scale_color_manual(values = sample(Choose_col,14)) +
    theme_cxf+
    guides(color=guide_legend(override.aes = list(size=5)))
  dev.off()
  

  
  FeaturePlot(seurat_Fb3,c("CD74","sPP1","MSLN","SRGN"))
  FeaturePlot(seurat_Fb3,c("TOP2A"))
  FeaturePlot(seurat_Fb3,c("CDC20"))
  FeaturePlot(seurat_Fb3,c("RGS5"))
  FeaturePlot(seurat_Fb3,c("DES"))
  FeaturePlot(seurat_Fb3,c("MSLN"))
  FeaturePlot(seurat_Fb3,c("CD24"))
  FeaturePlot(seurat_Fb3,c("CD37"))
  FeaturePlot(seurat_Fb3,c("PIEZO2"))
  FeaturePlot(seurat_Fb3,c("DCN"))
  FeaturePlot(seurat_Fb3,c("IL6"))
  FeaturePlot(seurat_Fb3,c("CD74","LRRC15","MSLN","SRGN","FAP","IL6"))
  FeaturePlot(seurat_Fb3,c("SPP1","CD24","PI16","CXCL8"))
  DimPlot(object = seurat_Fb3, reduction = 'umap',split.by = "Disease")
  FeaturePlot(seurat_Fb3,c("CD74","MSLN","SRGN","IL6"))
  FeaturePlot(seurat_Fb3,c("CD74","MSLN","SRGN","CD24"))
  FeaturePlot(seurat_Fb3,c("LRRC15","DPT","FAP","IL6"))
  FeaturePlot(seurat_Fb3,c("TOP2A"))
  FeaturePlot(seurat_Fb3,c("FAP"))
  FeaturePlot(seurat_Fb3,c("RGS5"))
  FeaturePlot(seurat_Fb3,c("CXCL8"))
  FeaturePlot(seurat_Fb3,c("CD74","MSLN","SRGN","IL6"))
  FeaturePlot(seurat_Fb3,c("LRRC15"))
  FeaturePlot(seurat_Fb3,c("PI16"))
  FeaturePlot(seurat_Fb3,c("SFRP2"))
  FeaturePlot(seurat_Fb3,c("COL15A1"))
  FeaturePlot(seurat_Fb3,c("TMEM119"))
  FeaturePlot(seurat_Fb3,c("DPT"))
  FeaturePlot(seurat_Fb3,c("HHIP"))
  FeaturePlot(seurat_Fb3,c("MYH11"))
  FeaturePlot(seurat_Fb3,c("CD74","SPP1","SRGN","LYVE1"))
  FeaturePlot(seurat_Fb3,c("CD74"), min.cutoff = 1.5)
  FeaturePlot(seurat_Fb3,c("IL6"))
  FeaturePlot(seurat_Fb3,c("CXCL8"))
  FeaturePlot(seurat_Fb3,c("HGF"))
  FeaturePlot(seurat_Fb3,c("KCNJ8"))
  FeaturePlot(seurat_Fb3,c("RGS5"))
  FeaturePlot(seurat_Fb3,c("DES"))
  FeaturePlot(seurat_Fb3,c("KCNA5"))
  FeaturePlot(seurat_Fb3,c("CDH2"))
  FeaturePlot(seurat_Fb3,c("SPP1"))
  FeaturePlot(seurat_Fb3,c("LYVE1"))
  FeaturePlot(seurat_Fb3,c("LUM"))
  FeaturePlot(seurat_Fb3,c("DCN"))
  FeaturePlot(seurat_Fb3,c("PDGFRA"))
  FeaturePlot(seurat_Fb3,c("ACTA2"),min.cutoff = 1.2)
  FeaturePlot(seurat_Fb3,c("IL6"),min.cutoff = 2)
  FeaturePlot(seurat_Fb3,c("MCAM"))
  
  dev.off()
  
  seurat_Fb5$seurat_clusters<-droplevels(seurat_Fb5$seurat_clusters)
  levels(seurat_Fb5$seurat_clusters)<-c( "0" , "1" , "2" , "3" , "4" , "5" , "6" , "7" , "8" , "3" , "10", "11" ,"12", "13")
  Idents(seurat_Fb5)<-"seurat_clusters"

  
  
  
  
  DimPlot(object = seurat_Fb5, reduction = 'umap',pt.size = 1,raster = F,split.by = "Group",
          label = F)+scale_color_manual(values=cols)+theme_cxf
  
  cols<-sample(Choose_col[1:17],13)
  scales::show_col(cols)
  cols<-c("#DCE775","#E91E63","#FFE082","#673AB7","#4DB6AC","#EF9A9A","#A5D6A7","#FF9800","#CE93D8","#AED581","#3F51B5","#2196F3","#FF5722") 
  embed_umap_A <- Embeddings(seurat_Fb5, 'umap')
  mat <- data.frame(seurat_Fb5@reductions$umap@cell.embeddings, 
                    group = seurat_Fb5$seurat_clusters)
  mat$group<-droplevels(mat$group)
  levels(mat$group)
  pdf("Seurat_umap_subtype.pdf", width = 18, height = 12)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = group)) +
    geom_point(size = 1e-3,alpha = 0.3) +
    scale_color_manual(values = cols) +
    theme_cxf+
    guides(color=guide_legend(override.aes = list(size=4)))
  dev.off()
  dim(seurat_Fb)
  dim(seurat_Fb4)
  dim(seurat_Fb5)
  
  
  
  
  
  
  
  
  
  library(cowplot) 
  FeaturePlot(seurat_Fb5, features = "DPT", label = F,cols = c("#e6e4df", "#f74343")) + theme_nothing()
 
  
 embed_umap <- Embeddings(seurat_Fb5, 'umap')
 genes<-c( "COL1A2","LUM","PDGFRA","RGS5",               # ALL  1-4
           "CD36","RERGL", "DES",                        # pericyte + SMC  5-7   
           "DPT","PI16","ICAM1",                         # ssCAF 8-10
           "LRRC15", "PLAT",                             # myCAF 11-12
           "CD74","CD37","CD24","SPP1",                  # apCAF 13-16
           "IL6","CXCL8","LIF","HGF",                    # iCAF  17-20
           "MEIS1", "TOP2A")                             # otherCAF 21-22

  
  P<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10",
       "p11","p12","p13","p14","p15","p16","p17","p18","p19","p20","p21","p22")
     
  for (i in 1:length(genes)){
  mat <- data.frame(seurat_Fb5@reductions$umap@cell.embeddings, 
                    group = seurat_Fb5@assays$RNA@data[genes[i],])
  x<-ggplot(mat, aes(x = umap_1, y = umap_2, color = group)) +
     geom_point(size = 1e-3,alpha = 0.1) +
     scale_color_gradient(low = "#e6e4df",high = "#f74343")+
     theme(legend.position="none",
           legend.title = element_blank(),
           legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
     theme_cxf2+ xlab(NULL) + ylab(NULL)+
     theme(axis.ticks = element_line(color='white'),
           axis.line  = element_line(color='white')) 
  assign(P[i],x)
  }
  
  COL1A2_value<-seurat_Fb5@assays$RNA@data["COL1A2",]
  COL1A2_value<-ifelse(COL1A2_value>5,5,COL1A2_value)
  p1 <-ggplot(mat, aes(x = umap_1, y = umap_2, color =  COL1A2_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  LUM_value<-seurat_Fb5@assays$RNA@data["LUM",]
  LUM_value<-ifelse(LUM_value>4,4,LUM_value)
  LUM_value<-ifelse(LUM_value<1.5,1.5,LUM_value)
  p2 <-ggplot(mat, aes(x = umap_1, y = umap_2, color =  LUM_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
 PDGFRA_value<-seurat_Fb5@assays$RNA@data["PDGFRA",]
 PDGFRA_value<-ifelse(PDGFRA_value>3,3,PDGFRA_value)
 PDGFRA_value<-ifelse(PDGFRA_value<0.5,0.5,PDGFRA_value)
 p3 <-ggplot(mat, aes(x = umap_1, y = umap_2, color =  PDGFRA_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white')) 
 RGS5_value<-seurat_Fb5@assays$RNA@data["RGS5",]
 RGS5_value<-ifelse(RGS5_value >5,5,RGS5_value)
 RGS5_value<-ifelse(RGS5_value <2,2,RGS5_value)
 p4 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = RGS5_value)) +
   geom_point(size = 1e-3,alpha = 0.1) +
   scale_color_gradient(low = "#e6e4df",high = "#f74343")+
   theme(legend.position="none",
         legend.title = element_blank(),
         legend.text= element_text(family = 'serif',face = 'italic',
                                   colour = 'black',size =10,hjust = .5)) +
   theme_cxf2+ xlab(NULL) + ylab(NULL)+
   theme(axis.ticks = element_line(color='white'),
         axis.line  = element_line(color='white'))   
 
 
  dev.off()
  
  p1+p2+p4
  
  # "CD36","RERGL", "DES", 5-7
  CD36_value<-seurat_Fb5@assays$RNA@data["CD36",]
  mean(CD36_value)
  CD36_value<-ifelse(CD36_value >3,3,CD36_value)
  CD36_value<-ifelse(CD36_value <1,1,CD36_value)
  p5 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = CD36_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))   
  p6
  DES_value<-seurat_Fb5@assays$RNA@data["DES",]
  mean(DES_value)
  DES_value<-ifelse(DES_value >4,4,DES_value)
  DES_value<-ifelse(DES_value <1,1,DES_value)
  p7 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = DES_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))   
  
  p5+p6+p7
  
  # "DPT","PI16","ICAM1",                         # ssCAF 8-10
  p9 
  PI16_value<-seurat_Fb5@assays$RNA@data["PI16",]
  mean(PI16_value)
  PI16_value<-ifelse(PI16_value >1,1,PI16_value)
  PI16_value<-ifelse(PI16_value <0.01,0,PI16_value)
  p9 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = PI16_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p10 
  ICAM1_value<-seurat_Fb5@assays$RNA@data["ICAM1",]
  mean(ICAM1_value)
  ICAM1_value<-ifelse(ICAM1_value >3,3,ICAM1_value)
  ICAM1_value<-ifelse(ICAM1_value <2,2,ICAM1_value)
  p10 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = ICAM1_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  
  p9+p10
  
  # "LRRC15", "PLAT",                             # myCAF 11-12
  p11
  LRRC15_value<-seurat_Fb5@assays$RNA@data["LRRC15",]
  mean(LRRC15_value)
  LRRC15_value<-ifelse(LRRC15_value >2,2,LRRC15_value)
  LRRC15_value<-ifelse(LRRC15_value <0.2,0.2,LRRC15_value)
  p11 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = LRRC15_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p12
  PLAT_value<-seurat_Fb5@assays$RNA@data["PLAT",]
  mean(PLAT_value)
  PLAT_value<-ifelse(PLAT_value >3,3,PLAT_value)
  PLAT_value<-ifelse(PLAT_value <1,1,PLAT_value)
  p12 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = PLAT_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p11+p12
  
  # "CD74","CD37","CD24","SPP1",                  # apCAF 13-16
  p13
  CD74_value<-seurat_Fb5@assays$RNA@data["CD74",]
  mean(CD74_value)
  CD74_value<-ifelse(CD74_value >3,3,CD74_value)
  CD74_value<-ifelse(CD74_value <0.5,0.5,CD74_value)
  p13 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = CD74_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p14
  CD37_value<-seurat_Fb5@assays$RNA@data["CD37",]
  mean(CD37_value)
  CD37_value<-ifelse(CD37_value >2,2,CD37_value)
  CD37_value<-ifelse(CD37_value <0,0,CD37_value)
  p14 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = CD37_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p15
  CD24_value<-seurat_Fb5@assays$RNA@data["CD24",]
  mean(CD24_value)
  CD24_value<-ifelse(CD24_value >2,2,CD24_value)
  CD24_value<-ifelse(CD24_value <0.75,0.75,CD24_value)
  p15 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = CD24_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))    
  p16
  SPP1_value<-seurat_Fb5@assays$RNA@data["SPP1",]
  mean(SPP1_value)
  SPP1_value<-ifelse(SPP1_value >2,2,SPP1_value)
  SPP1_value<-ifelse(SPP1_value <0.5,0.5,SPP1_value)
  p16 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = SPP1_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  p13+p14+p15+p16
  
 # "IL6","CXCL8","LIF","HGF",                    # iCAF  17-20
  
  p17
  IL6_value<-seurat_Fb5@assays$RNA@data["IL6",]
  mean(IL6_value)
  IL6_value<-ifelse(IL6_value >5,5,IL6_value)
  IL6_value<-ifelse(IL6_value <2,2,IL6_value)
  p17 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = IL6_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  p18
  CXCL8_value<-seurat_Fb5@assays$RNA@data["CXCL8",]
  mean(CXCL8_value)
  CXCL8_value<-ifelse(CXCL8_value >5,5,CXCL8_value)
  CXCL8_value<-ifelse(CXCL8_value <2.5,2.5,CXCL8_value)
  p18 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = CXCL8_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  p19
  LIF_value<-seurat_Fb5@assays$RNA@data["LIF",]
  mean(LIF_value)
  LIF_value<-ifelse(LIF_value >4,4,LIF_value)
  LIF_value<-ifelse(LIF_value <1,1,LIF_value)
  p19 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = LIF_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  p20
  HGF_value<-seurat_Fb5@assays$RNA@data["HGF",]
  mean(HGF_value)
  HGF_value<-ifelse(HGF_value >1.5,1.5,HGF_value)
  HGF_value<-ifelse(HGF_value <0,0,HGF_value)
  p20 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = HGF_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  p17+p18+p19+p20
  p21
  MEIS1_value<-seurat_Fb5@assays$RNA@data["MEIS1",]
  mean(MEIS1_value)
  MEIS1_value<-ifelse(MEIS1_value >3,3,MEIS1_value)
  MEIS1_value<-ifelse(MEIS1_value <1,1,MEIS1_value)
  p21 <-ggplot(mat, aes(x = umap_1, y = umap_2, color = MEIS1_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  p21+p22
  
  
  MCAM_value<-seurat_Fb5@assays$RNA@data["MCAM",]
  mean(MCAM_value)
  MCAM_value<-ifelse(MCAM_value >3,3,MCAM_value)
  MCAM_value<-ifelse(MCAM_value <1,1,MCAM_value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = MCAM_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  PDGFRA_value<-seurat_Fb5@assays$RNA@data["PDGFRA",]
  mean(PDGFRA_value)
  PDGFRA_value<-ifelse(PDGFRA_value >2,2,PDGFRA_value)
  PDGFRA_value<-ifelse(PDGFRA_value <0.5,0.5,PDGFRA_value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = PDGFRA_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  FAP_value<-seurat_Fb5@assays$RNA@data["FAP",]
  mean(FAP_value)
  FAP_value<-ifelse(FAP_value >2,2,FAP_value)
  FAP_value<-ifelse(FAP_value <0.5,0.5,FAP_value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = FAP_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  CDH11_value<-seurat_Fb5@assays$RNA@data["CDH11",]
  mean(FAP_value)
  CDH11_value<-ifelse(CDH11_value >4,4,CDH11_value)
  CDH11_value<-ifelse(CDH11_value <1.5,1.5,CDH11_value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = CDH11_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  POSTN_value<-seurat_Fb5@assays$RNA@data["POSTN",]
  mean(POSTN_value)
  POSTN_value<-ifelse(CDH11_value >4,4,POSTN_value)
  POSTN_value<-ifelse(CDH11_value <1,1,POSTN_value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = POSTN_value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  gene<-"CXCL6"
  value<-seurat_Fb5@assays$RNA@data[gene,]
  value<-ifelse(value >3,3,value)
  value<-ifelse(value <1,1,value)
  ggplot(mat, aes(x = umap_1, y = umap_2, color = value)) +
    geom_point(size = 1e-3,alpha = 0.1) +
    scale_color_gradient(low = "#e6e4df",high = "#f74343")+
    theme(legend.position="none",
          legend.title = element_blank(),
          legend.text= element_text(family = 'serif',face = 'italic',
                                    colour = 'black',size =10,hjust = .5)) +
    theme_cxf2+ xlab(NULL) + ylab(NULL)+
    theme(axis.ticks = element_line(color='white'),
          axis.line  = element_line(color='white'))  
  
  dev.off()
  
  
  
  
  remove(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22)
  dev.off()
  
  

  theme_cxf2<-  theme(panel.background = element_rect(fill = "transparent", color = "white",size = 0), 
                      legend.key = element_rect(fill = "transparent", color = "transparent"),
                      text = element_text(color = "black", size = 20,vjust=0.5),
                      plot.title = element_text(hjust = .5), 
                      axis.text = element_text(color = "white",size = 0))  
  
  
  
  markers <- FindAllMarkers(seurat_Fb5, only.pos = TRUE, logfc.threshold = 0.5,min.pct = 0.20)
  write.csv(markers, "./data/d_Fb6_Markers.csv")
  write.table(markers, "./data/d_Fb6_Markers.txt")
  top_markers<-markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
  gene<-unique(top_markers$gene)
  FeaturePlot(seurat_Fb5,c("CTSS"))
  seurat_Fb5$RNA_snn_res.0.4
  seurat_Fb5$seurat_clusters<-Idents(seurat_Fb5)
  s<-Idents(seurat_Fb5)
  levels(seurat_Fb5$seurat_clusters)<-c(
    "c02" , "c03" , "c09" , "c01",  "c04" , "c08" , "c07" , "c05" , "c10" , "c06" ,"c11" ,"c12", "c13"
  )
  seurat_Fb5$seurat_clusters<-fct_relevel(seurat_Fb5$seurat_clusters, c(
   "c01",  "c02" , "c03" , "c04" ,"c05" , "c06" ,"c07" , "c08" , "c09" , "c10" , "c11" ,"c12", "c13"
  ))
  table(seurat_Fb5$seurat_clusters)
  # 0     1     2     3     4     5     6     7     8    10    11    12    13 
  # 23300 21652  5569 24331 15907  6233  6732 11990  5141 10428  1328  2424   613 
  
  Idents(seurat_Fb5)<-"seurat_clusters"
  dim(seurat_Fb5)
  genes<-c( "RGS5","HIGD1B","STEAP4","CD36","FABP4",      #C01  CD36+pericyte     
            "CFD","ADH1B","PI16","WISP2","PLA2G2A",       #C02  PI16+ssCAF
            "FAP","LRRC15", "COL10A1","COL11A1","MMP11",  #C03  LRRC15+myCAF
            "MYH11","PLN","RERGL","BCAM","SORBS2",        #C04  RERGL+SMC  
            "CD74","HLA-DRA","HLA-DPA1","SRGN","CD37",    #C05  CD37+apCAF
            "EZR","ATP1B1","SPP1","CD24","SERPINA1",      #C06 CD24+apCAF
            "DPT","IGF1","ICAM1","RGS2","ARC",            #C07  ICAM1+ssCAF
            "POSTN","PLAT","F3","CCL11","HSD17B2",        #C08  PLAT+myCAF
            "CXCL8","TNFAIP6","CXCL2","HGF","MRPL44",     #C09  HGF+iCAF  
            "IL6","LIF","CSF3","SLC2A1","IL11",           #C10  LIF+iCAF
           "FTX","MEIS1", "PARD3B","FBXL7", "DLG2",       #C11 MEIS1+CAF
           "MYLK","CNN1","ACTG2","DES","RAMP1",           #C12 DES+SMC
           "TOP2A","CENPF","PTTG1","UBE2C","MKI67")       #C13 TOP2A+CAF
   length(unique(genes))    
   genes<-unique(genes)
   genes<-c("PTGS2","NMB","GAL","AKR1C","HGF","MRPL44","BST2")
  DotPlot(seurat_Fb5, features = genes,cols = c("#e6e4df", "#FF9800"), 
          col.min = 0, col.max = 2,dot.scale = 8)  + RotatedAxis() +theme_cxf
  dev.off()  
  

  
  d<-read.table("./data/c8_genes.txt")
  genes<-d$V1
  
  D_markers <- FindMarkers(seurat_Fb5, ident.1 = 8, ident.2 = 2, group.by = 'seurat_clusters')
  genes<-rownames(D_markers)
  for (i in 1:length(genes)){
    pdf(paste0("./result/new_CAF/c8/",genes[i],".pdf"), width = 5, height = 5)
    p<-FeaturePlot(seurat_Fb5, features =genes[i],
                   cols = c(rgb(225,225,225,10, maxColorValue = 255),
                            rgb(225,27,52,100,  maxColorValue = 255)))
    print(p)
    dev.off()
  }
  
  saveRDS(seurat_Fb5,"./data/d_Fb7.rds")
  seurat_Fb <- readRDS("./data/d_Fb3.rds")
  dim(seurat_Fb5)
}  
  
  
 
Cell_nk_dotplot{
    DOT<-DotPlot(seurat_Fb5, features = genes,cols = c('#1F78B4','#E31A1C'  ))+
         theme(axis.text.x = element_text(angle = 45, hjust = 1))+
         theme_cxf
    DOT_data <- as.data.frame(DOT$data)
    write.csv(DOT_data,"./result/new_CAF/FB6_dot_data.csv")
    colnames(DOT_data)
    d1<-DotPlot(seurat_Fb5, features = genes,cols = c('#1F78B4','#E31A1C'  ))+theme_cxf
    
    d2 <- ggplot(data = DOT_data)+ 
      geom_point(aes(x = id,y = features.plot, color = avg.exp.scaled, size = pct.exp))+
      scale_size_area(max_size = 6)+
      scale_color_gradientn(colors = c( '#1F78B4', '#FFFFB3', '#E31A1C')) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme_cxf+ 
      ylab("Cell Type") +
      xlab("") + coord_flip()  
    
    d3 <- ggplot(data = DOT_data)+ 
      geom_point(aes(x = id,y = features.plot, color = avg.exp.scaled, size = pct.exp))+
      scale_size_area(max_size = 6)+
      scale_color_gradientn(colors = c( "#4a7497",'#cba081', '#EC407A')) +
      cowplot::theme_cowplot() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))+
      theme_cxf+ 
      ylab("Cell Type") +
      xlab("") + coord_flip()  
    
    d1/d2/d3
    dev.off()
    
  }
  
subcluster_corerelation_analysis{
  library(pheatmap)
  ##install.packages("psych")
  library(psych)
  colnames(seurat_Fb5@meta.data)
  Idents(seurat_Fb5)="seurat_clusters"
  exp=AverageExpression(seurat_Fb5)
  coorda<-corr.test(exp$RNA,exp$RNA,method="spearman")
  pheatmap(coorda$r)
}  
 
 
caf_proportion{
  library(reshape2)
  df <- table(seurat_Fb5@meta.data$seurat_clusters,seurat_Fb5@meta.data$Disease) %>% melt()
  colnames(df) <- c("Cluster","Group","Number")
  write.csv(df,"./result/new_CAF/FB6_Proportion_data")
  df<-na.omit(df)
  ggplot(data = df, aes(x = Group, y = Number, fill =Cluster )) +
    geom_bar(stat = "identity", width=0.8,position="fill")+
    scale_fill_manual(values=cols) +
    theme_bw()+
    theme(panel.grid =element_blank()) +
    labs(x="",y="Ratio")+
    ####用来将y轴移动位置
    theme(axis.text.y = element_text(size=12, colour = "black"))+
    theme(axis.text.x = element_text(size=12, colour = "black"))+
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45))+
    theme_cxf
  
  dev.off()
  
  
  
  
}  
  
CAF_Cell_nk_distribution{
  seurat_A <-readRDS("./data/d_seurat_A.rds")
  meta<-seurat_Fb5@meta.data
  #write.csv(meta,"./data/seurat_Fb6_meta.csv")
  table(seurat_Fb5$Tissue)
  table(seurat_Fb5$Disease)
  seurat_Fb5$Tissue<-as.factor(seurat_Fb5$Tissue)
  levels(seurat_Fb5$Tissue)<-c("Bladder","Breast","Cervix","Colorectum","Colorectum","Esophagus","Kidney","Liver","Lung",   
                              "Omentum","Ovarian","Pancreas","Prostate","Skin","Stomach")
  seurat_Fb5$Disease<-as.factor(seurat_Fb5$Disease)
  levels(seurat_Fb5$Disease)<-c("BCC","BLCA","BRCA","ccRCC","CESC","COLO","ESCC","HCC",NA,"NSCLC","OV","PDAC","PMP",
                               "PRAD","SCC","STAD")
  meta <- read.csv("./data/seurat_Fb6_meta.csv",row.names = 1)
  
  
  
  
  # his_cell_num = meta %>% filter(meta_tissue == "Blood") 
  # his_cell_num = his_cell_num %>% separate(col = meta_histology_new, into = c("cancerType1", "cancerType2"), sep = "\\(")
  # his_cell_num <- his_cell_num %>% filter(!is.na(cancerType2))
  # his_cell_num = separate(data = his_cell_num, col = cancerType2, into = c("cancerType2", "cancerType3"), sep = "\\)")
  his_cell_num = meta 
  DT = his_cell_num  %>% group_by(Group,seurat_clusters)
  DT <- dplyr::summarise( DT,count = n())
  DT<-data.frame(DT)
  #write.csv(DT,"./result/new_CAF/FB6_Group_distrobition_meta.csv")
  
  test <- reshape2::dcast(DT,Group~seurat_clusters)
  test[is.na(test)] <- 0
  rownames(test)<-test$Group
  test<-test[,-1]
  test <- as.matrix(test)
  
  a =margin.table(test, 1)
  b=margin.table(test, 2)
  c=margin.table(test)
  test = test/(outer(a,b,"*")/c)
  
  tmp <- reshape2::melt(test)
  tmp <- tmp[order(tmp$Var1,tmp$value,tmp$Var2),]
  
  tmp$Var1 <- factor(tmp$Var1,levels=as.character(unique(tmp$Var1)))
  tmp$state <- ifelse(tmp$value>1,"Enrichment","Depletion")
  
  # col = c("#F4511E","#1976D2")
  # names(col) = c("Enrichment","Depletion")
  
  g <- ggplot(data = tmp)+ 
    geom_point(aes(x = Var1,y = Var2, color = state, size = value))+
    scale_size_area()+
    #scale_color_gradient(low = "grey",high = "red") + 
    #scale_color_manual(values = col) + 
    scale_color_manual(values = col) + 
    cowplot::theme_cowplot() +
    # theme(axis.title.y=element_text(size=12),axis.text.y =element_text(size = 12),
    #       axis.text.x = element_text(size = 12),
    #       strip.background = element_rect(colour = "white")) + 
    theme_cxf+
    ylab("Cell Type") +
    xlab("") 
  g
  
}

CAF_singnaling_score{
  dev.off()
  signatures<-read.csv("./data/seurat_Fb6_score.csv")
  signatures<-edit(signatures)
  Genesets<-as.vector(signatures)
  All_genes<-unlist(Genesets)
  #matrix
  meta<-as.data.frame(seurat_Fb5@meta.data)
  df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
  samp<-df %>% group_by(cluster) %>% sample_frac(size=0.3)  
  scRNAsub<-seurat_Fb5[,samp$barcode]
  exprMatrix<-scRNAsub@assays$RNA@data
  meta<-as.data.frame(scRNAsub@meta.data)
  dim(scRNAsub)

  ##AUC_score
  exprMatrix <- as(exprMatrix, "dgCMatrix")
  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=FALSE)
  cells_AUC <- AUCell_calcAUC(Genesets, cells_rankings)
  
  aucs <- getAUC(cells_AUC)
  
  aucs_t <- data.frame(t(aucs))  
  aucs_t$celltype <- meta[rownames(aucs_t), "seurat_clusters"]
  aucs_t <- aucs_t %>% group_by(celltype) %>% dplyr::summarise_each(funs = mean)
  getwd()
  write.csv(aucs_t ,"./result/new_CAF/FB6_CAF_AUCscore.csv")
  
  ####plot
  test1 <- aucs_t
  test1$celltype <- factor(test1$celltype)
  length(test1)
  P<-c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11")
  col<-c("#FFCCBC" ,"#FFCC80" ,"#FFE082" ,"#FFF9C4", "#E6EE9C", "#DCEDC8" ,"#C8E6C9", "#B2DFDB", "#80DEEA" ,"#B3E5FC")
  colnames(test1)
  
  p1 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Interleukin.Signaling,fill=col[1]),stat="identity",width = 0.8) + 
       geom_hline(aes(yintercept=median(Interleukin.Signaling)),linetype="dashed") + 
       ylab(colnames(test1)[2])+
       xlab("") +
       theme_classic() 
  p2 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Chemokine.Signaling,fill=col[2]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Chemokine.Signaling)),linetype="dashed") + 
    ylab(colnames(test1)[3])+
    xlab("") +
    theme_classic()    
  p3 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= extracellular.matrix..ECM.,fill=col[3]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(extracellular.matrix..ECM.)),linetype="dashed") + 
    ylab(colnames(test1)[4])+
    xlab("") +
    theme_classic()   
  p4 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= MHC.II.Antigen.Presentation,fill=col[4]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(MHC.II.Antigen.Presentation)),linetype="dashed") + 
    ylab(colnames(test1)[5])+
    xlab("") +
    theme_classic()    
  p5 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Angiogenesis,fill=col[5]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Angiogenesis)),linetype="dashed") + 
    ylab(colnames(test1)[6])+
    xlab("") +
    theme_classic()    
  p6 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Hypoxia,fill=col[6]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Hypoxia)),linetype="dashed") + 
    ylab(colnames(test1)[7])+
    xlab("") +
    theme_classic()    
  p7 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Cell.Proliferation,fill=col[7]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Cell.Proliferation)),linetype="dashed") + 
    ylab(colnames(test1)[8])+
    xlab("") +
    theme_classic()    
  p8 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= Hedgehog.Signaling,fill=col[8]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(Hedgehog.Signaling)),linetype="dashed") + 
    ylab(colnames(test1)[9])+
    xlab("") +
    theme_classic()    
  p9 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= NF.kappaB.Signaling,fill=col[9]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(NF.kappaB.Signaling)),linetype="dashed") + 
    ylab(colnames(test1)[10])+
    xlab("") +
    theme_classic()    
  p10 <- test1  %>% ggplot() + geom_bar( aes(x =celltype ,y= TGF.beta.Signaling,fill=col[10]),stat="identity",width = 0.8) + 
    geom_hline(aes(yintercept=median(TGF.beta.Signaling)),linetype="dashed") + 
    ylab(colnames(test1)[11])+
    xlab("") +
    theme_classic()    
  
  
  dev.off()
  library(patchwork)
  p <- p1/p2/p3/p4/p5/p6/p7/p8/p9/p10
  p
  rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
  ggsave("./result/new_CAF/FB6_CAF_AUCscore.pdf",p, height = 6, width =7)
  
}

CAF_Gradient_volcano{
  seurat_Fb <- readRDS("./data/d_Fb7.rds")
  dim(seurat_Fb)
  
  Idents(seurat_Fb)
  
  seurat_Fb_CAFs<-subset(seurat_Fb, idents = c("c02","c03","c05","c06","c07","c08","c09","c10","c11","c13"))
  seurat_Fb_Percyte <-subset(seurat_Fb, idents = c("c01"))
  seurat_Fb_SMCs <-subset(seurat_Fb, idents = c("c04","c12"))
  
  #cAFs
  Idents(seurat_Fb_CAFs)<-"Group"
  degdf <- FindMarkers(seurat_Fb_CAFs,ident.1 = "Tumor",ident.2 = "Normal", 
                       logfc.threshold = 0.01)
  colnames(degdf)
  logFC_t=0
  P.Value_t = 1e-10
  #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
  degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                       ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))
  degdf$name=rownames(degdf)
  library(ggpubr)
  attach(degdf)
  detach(degdf)
  degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
  x<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up",])[1]
  y<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down",])[1]
  z<-dim(degdf)[1]-x-y
  degdf$"-log10(p_val_adj)"<- ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up", 
                                     sample(runif(x,10,300),x),
                                     ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down", 
                                            sample(runif(y,10,300),y),      
                                            degdf$"-log10(p_val_adj)"))
  degdf$pt<-log(degdf$pct.1/degdf$pct.2+1)
  max(degdf$pt)
  degdf$gene <- rownames(degdf)
  degdf$label <- rownames(degdf)
  degdf$label [degdf$group == "stable"] <- NA
  write_excel_csv(degdf,"./result/new_CAF/FB6_CAFs_deg.csv")

  data<-degdf
  ggplot(data,aes(avg_log2FC, `-log10(p_val_adj)`,color=group,size=pt))+
    geom_hline(yintercept = 10, linetype = "dashed", color = "#999999")+
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
    geom_point(aes(size = pt, color= group))+
    scale_color_manual(values = c("#39489f","gray80","#b81f25"))+
    theme_minimal()+
    geom_text(aes(label=label,size = 5, vjust = 0, hjust=-0.2))+
    theme_cxf
  dev.off()
  
  #Percyte
  Idents(seurat_Fb_Percyte)<-"Group"
  degdf <- FindMarkers(seurat_Fb_Percyte,ident.1 = "Tumor",ident.2 = "Normal", 
                       logfc.threshold = 0.01)
  colnames(degdf)
  logFC_t=0
  P.Value_t = 1e-10
  #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
  degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                       ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))

  library(ggpubr)
  degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
  x<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up",])[1]
  y<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down",])[1]
  z<-dim(degdf)[1]-x-y
  degdf$"-log10(p_val_adj)"<- ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up", 
                                     sample(runif(x,10,300),x),
                                     ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down", 
                                            sample(runif(y,10,300),y),      
                                            degdf$"-log10(p_val_adj)"))
  degdf$pt<-log(degdf$pct.1/degdf$pct.2+1)
  max(degdf$pt)
  degdf$gene <- rownames(degdf)
  degdf$label <- rownames(degdf)
  degdf$label [degdf$group == "stable"] <- NA
  write_excel_csv(degdf,"./result/new_CAF/FB6_pericyte_deg.csv")
  
  data<-degdf
  ggplot(data,aes(avg_log2FC, `-log10(p_val_adj)`,color=group,size=pt))+
    geom_hline(yintercept = 10, linetype = "dashed", color = "#999999")+
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
    geom_point(aes(size = pt, color= group))+
    scale_color_manual(values = c("#39489f","gray80","#b81f25"))+
    theme_minimal()+
    geom_text(aes(label=label,size = 5, vjust = 0, hjust=-0.2))+
    theme_cxf
  dev.off()
  
  
  
  #SMCs
  Idents(seurat_Fb_SMCs)<-"Group"
  degdf <- FindMarkers(seurat_Fb_SMCs,ident.1 = "Tumor",ident.2 = "Normal", 
                       logfc.threshold = 0.01)
  colnames(degdf)
  logFC_t=0
  P.Value_t = 1e-10
  #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
  degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                       ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))
  
  library(ggpubr)
  degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
  x<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up",])[1]
  y<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down",])[1]
  z<-dim(degdf)[1]-x-y
  degdf$"-log10(p_val_adj)"<- ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up", 
                                     sample(runif(x,10,300),x),
                                     ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down", 
                                            sample(runif(y,10,300),y),      
                                            degdf$"-log10(p_val_adj)"))
  degdf$pt<-log(degdf$pct.1/degdf$pct.2+1)
  max(degdf$pt)
  degdf$gene <- rownames(degdf)
  degdf$label <- rownames(degdf)
  degdf$label [degdf$group == "stable"] <- NA
  write_excel_csv(degdf,"./result/new_CAF/FB6_SMCs_deg.csv")
  
  data<-degdf
  ggplot(data,aes(avg_log2FC, `-log10(p_val_adj)`,color=group,size=pt))+
    geom_hline(yintercept = 10, linetype = "dashed", color = "#999999")+
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
    geom_point(aes(size = pt, color= group))+
    scale_color_manual(values = c("#39489f","gray80","#b81f25"))+
    theme_minimal()+
    geom_text(aes(label=label,size = 5, vjust = 0, hjust=-0.2))+
    theme_cxf
  dev.off()
  
  
  
  
}


cellchat_CAFs_T{
  dev.off()
  d_CAF <-readRDS("./data/d_Fb7.rds")
  d_T <- readRDS("./data/d_T.rds")
  dim(d_CAF); dim(d_T)
  # plot(density(d_T@assays$RNA$data["CD3D",]))
  # A<-colnames(d_T@assays$RNA$data)[d_T@assays$RNA$data["CD3D",]>0.5]
  # d_T <- d_T[,A]
  table(d_T$seurat_clusters)
  table(d_CAF$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_T@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_T_subset<-d_T[,samp$barcode]
  }
  saveRDS(d_T_subset,"./data/d_T_subset.rds")
  d_T <- readRDS("./data/d_T_subset.rds")
  d_T_subset<-d_T
  rm(d_T)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_CAF@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    # samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.5)  
    d_CAF_subset<-d_CAF[,samp$barcode]
  }
  
  CAF_T <- merge(d_CAF_subset,y = c(d_T_subset))
  CAF_T <- NormalizeData(CAF_T)
  saveRDS(CAF_T,"./data/CAF_T.rds")
  CAF_T <- readRDS("D:/data_subset5/subset/seurat_all/data/CAF_T.rds")
  all_genes<-rownames(CAF_T)
  
  CAF_T <- readRDS("./data/CAF_T.rds")
  table(CAF_T$seurat_clusters)
  CAF_T$seurat_clusters<-as.factor(CAF_T$seurat_clusters)
  levels(CAF_T$seurat_clusters)
  CAF_T$seurat_clusters<-factor(CAF_T$seurat_clusters,levels=c(
    "c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13",  
    "CD16+NK","CD56+NK","CD8Teff","CD8Tn", "Tex","TIFIT3","TMCM6","TMKI67",
    "CD4Teff","CD4Tn","Tfh","Th17","Treg"
  ))
  CAF_T$seurat_clusters<- droplevels(CAF_T$seurat_clusters)
  table(CAF_T$seurat_clusters)
  

  Idents(CAF_T)<-"Group"
  CAF_T_N<-subset(CAF_T,idents = c("Normal"))
  table(CAF_T_N$seurat_clusters)
  CAF_T_T<-subset(CAF_T,idents = c("Tumor"))
  table(CAF_T_T$seurat_clusters)
  
  CAF_T_N{
    library(CellChat)#载入R包
    data.input <- LayerData(CAF_T_N, assay = "RNA", layer = "data")
    levels(CAF_T_N$seurat_clusters)
    identity <- subset(CAF_T_N@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_N <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_N@DB <- CellChatDB
    cellchat_N <- subsetData(cellchat_N)
    cellchat_N <- identifyOverExpressedGenes(cellchat_N)
    cellchat_N <- identifyOverExpressedInteractions(cellchat_N)
    cellchat_N <- projectData(cellchat_N, PPI.human)
    
    cellchat_N <- computeCommunProb(cellchat_N, type = "truncatedMean")  #"triMean"
    cellchat_N <- filterCommunication(cellchat_N, min.cells = 3)
    df.net <- subsetCommunication(cellchat_N)
    A<-levels(CAF_T_N$seurat_clusters)[c(1:13)]
    B<-levels(CAF_T_N$seurat_clusters)[c(14:26)]
    df.net <- subsetCommunication(cellchat_N, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/new_CAF/cellchat_CAF_T/N_df_net.csv")
    cellchat_N <- computeCommunProbPathway(cellchat_N)
    cellchat_N <- aggregateNet(cellchat_N)
    groupSize <- as.numeric(table(cellchat_N@idents))
    
    saveRDS(cellchat_N,"./result/new_CAF/cellchat_CAF_T/cellchat_N.rds")
  }
  
  CAF_T_T{
    library(CellChat)#载入R包
    data.input <- LayerData(CAF_T_T, assay = "RNA", layer = "data")
    identity <- subset(CAF_T_T@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_T <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_T@DB <- CellChatDB
    cellchat_T <- subsetData(cellchat_T)
    cellchat_T <- identifyOverExpressedGenes(cellchat_T)
    cellchat_T <- identifyOverExpressedInteractions(cellchat_T)
    cellchat_T <- projectData(cellchat_T, PPI.human)
    
    cellchat_T <- computeCommunProb(cellchat_T, type = "truncatedMean")
    cellchat_T <- filterCommunication(cellchat_T, min.cells = 3)
    df.net <- subsetCommunication(cellchat_T)
    A<-levels(CAF_T_T$seurat_clusters)[c(1:13)]
    B<-levels(CAF_T_T$seurat_clusters)[c(14:26)]
    df.net <- subsetCommunication(cellchat_T, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/new_CAF/cellchat_CAF_T/T_df_net.csv")
    cellchat_T <- computeCommunProbPathway(cellchat_T)
    cellchat_T <- aggregateNet(cellchat_T)
    groupSize <- as.numeric(table(cellchat_T@idents))
    saveRDS(cellchat_T,"./result/new_CAF/cellchat_CAF_T/cellchat_T.rds")
  }  
  
  CAF_T_NT{
    library(CellChat)
    library(patchwork)
    cellchat_N  <- readRDS("./result/new_CAF/cellchat_CAF_T/cellchat_N.rds")
    cellchat_T  <- readRDS("./result/new_CAF/cellchat_CAF_T/cellchat_T.rds")
    setwd("./result/new_CAF/cellchat_CAF_T")
    data.dir <- './comparison'
    dir.create(data.dir)
    setwd("./comparison")
    
    cellchat.NL <- cellchat_N
    cellchat.LS <- cellchat_T
    object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    A<-levels(CAF_T_T$seurat_clusters)[c(1:13)]
    B<-levels(CAF_T_T$seurat_clusters)[c(14:26)]
    gg1 <- netVisual_heatmap(cellchat,sources.use = A, targets.use = B, remove.isolate =T)
    gg1_data<-as.matrix(gg1@matrix)
    gg1_data<-as.data.frame(gg1_data)
    gg1_data<-gg1_data[A,B]
    ComplexHeatmap::pheatmap(
      gg1_data, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    dev.off()

    
    #Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
    #Identify dysfunctional signaling by comparing the communication probabities
  
    NV <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), 
                           color.heatmap = c("Spectral", "viridis"),
                           signaling = pathways,
                           n.colors = 10,line.size = 0.1,
                           return.data = T,
                           angle.x = 45)
    NV_data <- NV$data
    write.csv(NV_data,"NT10_data.csv")
    dev.off()  
    
    pathways<-
    c("CSPG4","THY1","COMPLEMENT","TRAIL","CD86","DESMOSOME","PERIOSTIN","SEMA4","LT","SPP1","ANGPT",   
      "NOTCH","ESAM","FASLG","CD46","IGF","ANGPTL","MHC-II","LIGHT","TGFb","CCL",  
      "ICAM","PARs","THBS","APP","MK","FN1","COLLAGEN","PTPRM","LAMININ",
      "BAG","CD40","PTN","MIF","ADGRE5","CLEC","TNF","FGF")   

    NB <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), 
                           #color.heatmap = c("Spectral", "viridis"),
                           color.heatmap = c("viridis"),
                           signaling = pathways,
                           n.colors = 5,line.size = 0.1,
                           return.data = T,
                           angle.x = 45)
    NB_communication<-NB$communication
    write.csv(NB_communication,"NB_communication.csv")
    
    
    
    

    
    
    
    
    gg1 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 2,
                            min.quantile = 0,max.quantile = 1,
                            title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg2 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
    #> Comparing communications on a merged object
    gg1 + gg2
    dev.off()  
    g1_data <- gg1$data
    write.csv(g1_data,"NB_UP_data.csv")
    g2_data <- gg2$data
    write.csv(g2_data,"NB_DOWN_data.csv")
    dev.off()  
    
    
    
    #Identify dysfunctional signaling by using differential expression analysis
    # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
    pos.dataset = "LS"
    features.name = pos.dataset
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = T, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.01)
    net <- netMappingDEG(cellchat, features.name = features.name)
    # extract the ligand-receptor pairs with upregulated ligands in LS
    net.up <- subsetCommunication(cellchat, net = net, datasets = "LS", 
                                  ligand.logFC = 0.5, receptor.logFC = 0.1,
                                  ligand.pvalues = 0.01, receptor.pvalues = 0.05, ligand.pct.1 = 0.25, receptor.pct.1 = 0.25)
    gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    unique( pairLR.use.up)
   library(viridis)
    
    setwd("./result/new_CAF/cellchat_CAF_T/comparison")
    data.dir <- './c_comparison'
    dir.create(data.dir)
    setwd("./c_comparison")
    
    
    
    

    A<-levels(cellchat@meta$seurat_clusters)[c(1:13)]
    B<-levels(cellchat@meta$seurat_clusters)[c(14:26)]
for(i in 1:length(A)){
    C <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = A[i], targets.use = B, comparison = c(1, 2),  
                            angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ","LS"),
                            color.heatmap = c("viridis"),
                            n.colors = 5,line.size = 0.1,
                            return.data = T)
    C_communication<-C$communication
    write.csv(C_communication,paste0("c",i,"_communication.csv"))
       
    df<-C_communication
    # df[df$source.target=="c01 -> CD4Tn (NL)",]<-NA
    # df<-na.omit(df)
    # df$source.target<-droplevels(df$source.target)
    dataset.name.order <- levels(df$source.target)
    dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
    dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order) - 1)
    xtick.color <- stringr::str_replace_all(dataset.name.order,c("LS"),c("#DC0000CC"))
    xtick.color <- stringr::str_replace_all(xtick.color,c("NL"),c("#3C5488CC"))
    
    g<-ggplot(df, aes(x = source.target, y = interaction_name_2, 
                         size= pval, color = prob)) + geom_point(pch = 16) + 
      theme_linedraw() + theme(panel.grid.major = element_blank()) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                       vjust = 1), axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + 
      scale_x_discrete(position = "bottom") +
      #scale_color_manual(values=c( "#E64B35CC","gray80"))+
      # geom_hline(yintercept = seq(1.5, length(unique(C01_communication$interaction_name_2)) - 
      #                                 0.5, 1), lwd = 0.1, colour = "gray95")+
      geom_vline(xintercept = seq(0.5, length(unique(df$source.target)) + 
                                  0.5, 2), lwd = 0.1, colour = "gray95")+
      scale_radius(range = c(min(df$pval), max(df$pval)), 
                    breaks = sort(unique(df$pval)))+
      scale_color_viridis(option = "A",alpha = 0.8,discrete = F)+  
      theme(axis.text.x = element_text(colour = xtick.color))+
      theme_cxf
    
    wt<-length(unique(df$source.target))*0.3+3
    ht<-length(unique(df$interaction_name_2))*0.3+1.5
    ggsave(paste0("c",i,"_communication.pdf"),g, width = wt, height = ht)
  

    
    lr<-levels(df$interaction_name_2)
    for(j in 1:length(lr)){
    l<-lr[j]
    
    x<-strsplit(l," ")[[1]][1]
    y<-strsplit(l," ")[[1]][3]
    z<-str_replace(y,"[(]","")  #stringr::str_sub( y, 2, stringr::str_length(y) - 1)
    z<-str_replace(z,"[)]","")
    h<-str_replace(z,"[+]"," ")
    h<-str_replace(h,"[:]"," ")
    h1<-strsplit(h," ")[[1]][1]
    h2<-strsplit(h," ")[[1]][2]
    f<-c(x,h1,h2)
    f2<-c(x,y)
    
    if(NA %in% f){
        f <- f2
    }else{
        f <- f
        }
     f<-f[f %in% all_genes]
  
  
    p<-plotGeneExpression(cellchat,  split.by = "datasets", 
                       #signaling = "COLLAGEN",
                       features = f,
                       type =   "violin" , #c("violin", "dot", "bar"),
                       color.use =  c("#3C5488CC","#DC0000CC") )
    

    ggsave(paste0("c",i,"_communication_",x,"_",h1,"_",h2,".pdf"),p, width = 12, height = 3)
    }
    
    } 
   
    
    
 
    saveRDS(cellchat, file = "cellchat_NT.rds")
    cellchat<-readRDS("./result/new_CAF/cellchat_CAF_T/comparison/cellchat_NT.rds")
    cellchat<-readRDS("../cellchat_NT.rds")
  }  
  
}


cellchat_CAFs_mono{
  setwd("D:/data_subset5/subset/seurat_all/")
  d_CAF <-readRDS("./data/d_Fb7.rds")
  d_Mye <- readRDS("./data/d_Mye.rds")
  dim(d_CAF); dim(d_Mye)
  table(d_Mye$seurat_clusters)
  table(d_CAF$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_Mye@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=2500)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_Mye_subset<-d_Mye[,samp$barcode]
  }
  saveRDS( d_Mye_subset,"./data/ d_Mye_subset.rds")
  d_Mye_subset <- readRDS("./data/ d_Mye_subset.rds")
  rm(d_Mye)
  
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_CAF@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    # samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.5)  
    d_CAF_subset<-d_CAF[,samp$barcode]
  }
  d_CAF_subset<-subset(CAF_T,idents = A)
  table(d_CAF_subset$seurat_clusters)
  saveRDS(d_CAF_subset,"./data/d_CAF_subset.rds")
  
  
  
  CAF_Mye <- merge(d_CAF_subset,y = c(d_Mye_subset))
  saveRDS(CAF_Mye,"./data/CAF_Mye.rds")
  CAF_Mye <- readRDS("D:/data_subset5/subset/seurat_all/data/CAF_Mye.rds")
  all_genes<-rownames(CAF_Mye)
  
  
  table(CAF_Mye$seurat_clusters)
  CAF_Mye$seurat_clusters<-as.factor(CAF_Mye$seurat_clusters)
  levels(CAF_Mye$seurat_clusters)
  CAF_Mye$seurat_clusters<-factor(CAF_Mye$seurat_clusters,levels=c(
    "c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13",  
    "FMN1+Mac" ,"FOLR2+Mac","MARCO+Mac" ,"CD16+Mo" ,"FCN1+Mo",
    "LAMP3+mregDC" ,"CLEC9A+cDC1s", "CD1C+cDC2s" ,"MKI67+Mye"))
  CAF_Mye$seurat_clusters<- droplevels(CAF_Mye$seurat_clusters)
  table(CAF_Mye$seurat_clusters)
  
  Idents(CAF_Mye)<-"Group"
  CAF_Mye_N<-subset(CAF_Mye,idents = c("Normal"))
  table(CAF_Mye_N$seurat_clusters)
  CAF_Mye_T<-subset(CAF_Mye,idents = c("Tumor"))
  table(CAF_Mye_T$seurat_clusters)
  
  CAF_Mye_N{
    library(CellChat)#载入R包
    data.input <- LayerData(CAF_Mye_N, assay = "RNA", layer = "data")
    levels(CAF_Mye_N$seurat_clusters)
    identity <- subset(CAF_Mye_N@meta.data, select = "seurat_clusters")
    levels(identity)
    cellchat_N <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_N@DB <- CellChatDB
    cellchat_N <- subsetData(cellchat_N)
    cellchat_N <- identifyOverExpressedGenes(cellchat_N)
    cellchat_N <- identifyOverExpressedInteractions(cellchat_N)
    cellchat_N <- projectData(cellchat_N, PPI.human)
    
    cellchat_N <- computeCommunProb(cellchat_N, type = "truncatedMean")  #"triMean"
    cellchat_N <- filterCommunication(cellchat_N, min.cells = 3)
    df.net <- subsetCommunication(cellchat_N)
    A<-levels(CAF_Mye_N$seurat_clusters)[c(1:13)]
    B<-levels(CAF_Mye_N$seurat_clusters)[c(14:22)]
    df.net <- subsetCommunication(cellchat_N, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/new_CAF/cellchat_CAF_Mye/N_df_net.csv")
    cellchat_N <- computeCommunProbPathway(cellchat_N)
    cellchat_N <- aggregateNet(cellchat_N)
    groupSize <- as.numeric(table(cellchat_N@idents))
    
    saveRDS(cellchat_N,"./result/new_CAF/cellchat_CAF_Mye/cellchat_N.rds")
  }
  
  CAF_Mye_T{
    data.input <- LayerData(CAF_Mye_T, assay = "RNA", layer = "data")
    identity <- subset(CAF_Mye_T@meta.data, select = "seurat_clusters")
    cellchat_T <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
    CellChatDB <- CellChatDB.human
    cellchat_T@DB <- CellChatDB
    cellchat_T <- subsetData(cellchat_T)
    cellchat_T <- identifyOverExpressedGenes(cellchat_T)
    cellchat_T <- identifyOverExpressedInteractions(cellchat_T)
    cellchat_T <- projectData(cellchat_T, PPI.human)
    
    cellchat_T <- computeCommunProb(cellchat_T, type = "truncatedMean")
    cellchat_T <- filterCommunication(cellchat_T, min.cells = 3)
    df.net <- subsetCommunication(cellchat_T)
    A<-levels(CAF_Mye_T$seurat_clusters)[c(1:13)]
    B<-levels(CAF_Mye_T$seurat_clusters)[c(14:22)]
    df.net <- subsetCommunication(cellchat_T, sources.use = A, targets.use = B) 
    write.csv(df.net,"./result/new_CAF/cellchat_CAF_Mye/T_df_net.csv")
    cellchat_T <- computeCommunProbPathway(cellchat_T)
    cellchat_T <- aggregateNet(cellchat_T)
    groupSize <- as.numeric(table(cellchat_T@idents))
    saveRDS(cellchat_T,"./result/new_CAF/cellchat_CAF_Mye/cellchat_T.rds")
#  }  
  
#  CAF_Mye_NT{
    library(CellChat)
    library(patchwork)
    #cellchat_N  <- readRDS("./result/new_CAF/cellchat_CAF_Mye/cellchat_N.rds")
   #cellchat_T  <- readRDS("./result/new_CAF/cellchat_CAF_Mye/cellchat_T.rds")
    setwd("./result/new_CAF/cellchat_CAF_Mye")
    data.dir <- './comparison'
    dir.create(data.dir)
    setwd("./comparison")
    
    cellchat.NL <- cellchat_N
    cellchat.LS <- cellchat_T
    object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    
    A<-levels(CAF_Mye_T$seurat_clusters)[c(1:13)]
    B<-levels(CAF_Mye_T$seurat_clusters)[c(14:22)]
    gg1 <- netVisual_heatmap(cellchat,sources.use = A, targets.use = B, remove.isolate =T)
    gg1_data<-as.matrix(gg1@matrix)
    gg1_data<-as.data.frame(gg1_data)
    gg1_data<-gg1_data[A,B]
    
    pdf(file = "./heatmap.pdf")
    ComplexHeatmap::pheatmap(
      gg1_data, 
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), 
      border_color = "white",
      cellwidth = 16, cellheight = 16, fontsize = 14)
    dev.off()
    
    
    #Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
    #Identify dysfunctional signaling by comparing the communication probabities
    # 
    # NV <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), 
    #                        color.heatmap = c("Spectral", "viridis"),
    #                        signaling = pathways,
    #                        n.colors = 10,line.size = 0.1,
    #                        return.data = T,
    #                        angle.x = 45)
    # NV_data <- NV$data
    # write.csv(NV_data,"NT10_data.csv")
    # dev.off()  
    # 
    # 
    # NB <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), 
    #                        #color.heatmap = c("Spectral", "viridis"),
    #                        color.heatmap = c("viridis"),
    #                        signaling = pathways,
    #                        n.colors = 5,line.size = 0.1,
    #                        return.data = T,
    #                        angle.x = 45)
    # NB_communication<-NB$communication
    # write.csv(NB_communication,"NB_communication.csv")
    # 
    # 
    # 
    # 
    # gg1 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 2,
    #                         min.quantile = 0,max.quantile = 1,
    #                         title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)
    # #> Comparing communications on a merged object
    # gg2 <- netVisual_bubble(cellchat, sources.use = A, targets.use = B,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
    # #> Comparing communications on a merged object
    # gg1 + gg2
    # dev.off()  
    # g1_data <- gg1$data
    # write.csv(g1_data,"NB_UP_data.csv")
    # g2_data <- gg2$data
    # write.csv(g2_data,"NB_DOWN_data.csv")
    # dev.off()  
    # 
    # 
    
    #Identify dysfunctional signaling by using differential expression analysis
    # define a positive dataset, i.e., the dataset with positive fold change against the other dataset
    pos.dataset = "LS"
    features.name = pos.dataset
    cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = T, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.01)
    net <- netMappingDEG(cellchat, features.name = features.name)
    # extract the ligand-receptor pairs with upregulated ligands in LS
    net.up <- subsetCommunication(cellchat, net = net, datasets = "LS", 
                                  ligand.logFC = 0.5, receptor.logFC = 0.1,
                                  ligand.pvalues = 0.01, receptor.pvalues = 0.05, ligand.pct.1 = 0.25, receptor.pct.1 = 0.25)
    gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
    pairLR.use.up = net.up[, "interaction_name", drop = F]
    unique( pairLR.use.up)
    library(viridis)
    getwd()

    data.dir <- './c_comparison'
    dir.create(data.dir)
    setwd("./c_comparison")
    
    
    A<-levels(cellchat@meta$seurat_clusters)[c(1:13)]
    B<-levels(cellchat@meta$seurat_clusters)[c(14:22)]
    for(i in 1:length(A)){
      C <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = A[i], targets.use = B, comparison = c(1, 2),  
                            angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ","LS"),
                            color.heatmap = c("viridis"),
                            n.colors = 5,line.size = 0.1,
                            return.data = T)
      C_communication<-C$communication
      write.csv(C_communication,paste0("c",i,"_communication.csv"))
      
      df<-C_communication
      # df[df$source.target=="c01 -> CD4Tn (NL)",]<-NA
      # df<-na.omit(df)
      # df$source.target<-droplevels(df$source.target)
      dataset.name.order <- levels(df$source.target)
      dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
      dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order) - 1)
      xtick.color <- stringr::str_replace_all(dataset.name.order,c("LS"),c("#DC0000CC"))
      xtick.color <- stringr::str_replace_all(xtick.color,c("NL"),c("#3C5488CC"))
      
      g<-ggplot(df, aes(x = source.target, y = interaction_name_2, 
                        size= pval, color = prob)) + geom_point(pch = 16) + 
        theme_linedraw() + theme(panel.grid.major = element_blank()) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                         vjust = 1), axis.title.x = element_blank(), 
              axis.title.y = element_blank()) + 
        scale_x_discrete(position = "bottom") +
        #scale_color_manual(values=c( "#E64B35CC","gray80"))+
        # geom_hline(yintercept = seq(1.5, length(unique(C01_communication$interaction_name_2)) - 
        #                                 0.5, 1), lwd = 0.1, colour = "gray95")+
        geom_vline(xintercept = seq(0.5, length(unique(df$source.target)) + 
                                      0.5, 2), lwd = 0.1, colour = "gray95")+
        scale_radius(range = c(min(df$pval), max(df$pval)), 
                     breaks = sort(unique(df$pval)))+
        scale_color_viridis(option = "A",alpha = 0.8,discrete = F)+  
        theme(axis.text.x = element_text(colour = xtick.color))+
        theme_cxf
      
      wt<-length(unique(df$source.target))*0.3+3
      ht<-length(unique(df$interaction_name_2))*0.3+1.5
      ggsave(paste0("c",i,"_communication.pdf"),g, width = wt, height = ht)
      
      
      
      lr<-levels(df$interaction_name_2)
      for(j in 1:length(lr)){
        l<-lr[j]
        
        x<-strsplit(l," ")[[1]][1]
        y<-strsplit(l," ")[[1]][3]
        z<-str_replace(y,"[(]","")  #stringr::str_sub( y, 2, stringr::str_length(y) - 1)
        z<-str_replace(z,"[)]","")
        h<-str_replace(z,"[+]"," ")
        h<-str_replace(h,"[:]"," ")
        h1<-strsplit(h," ")[[1]][1]
        h2<-strsplit(h," ")[[1]][2]
        f<-c(x,h1,h2)
        f2<-c(x,y)
        
        if(NA %in% f){
          f <- f2
        }else{
          f <- f
        }
        f<-f[f %in% all_genes]
        
        
        p<-plotGeneExpression(cellchat,  split.by = "datasets", 
                              #signaling = "COLLAGEN",
                              features = f,
                              type =   "violin" , #c("violin", "dot", "bar"),
                              color.use =  c("#3C5488CC","#DC0000CC") )
        
        
        ggsave(paste0("c",i,"_communication_",x,"_",h1,"_",h2,".pdf"),p, width = 12, height = 3)
      }
      
    } 
    
    saveRDS(cellchat, file = "cellchat_NT.rds")
    cellchat<-readRDS("./result/new_CAF/cellchat_CAF_Mye/comparison/cellchat_NT.rds")
    cellchat<-readRDS("../cellchat_NT.rds")
  }  
  
}

cellchat_apCAFs_T_mono{
  apCAF5_T_data  <- read.csv("./result/new_CAF/cellchat_CAF_T/comparison/c_comparison/c5_communication.csv")
  apCAF5_Mye_data  <- read.csv("./result/new_CAF/cellchat_CAF_Mye/comparison/c_comparison/c5_communication.csv")
  apCAF6_T_data  <- read.csv("./result/new_CAF/cellchat_CAF_T/comparison/c_comparison/c6_communication.csv")
  apCAF6_Mye_data  <- read.csv("./result/new_CAF/cellchat_CAF_Mye/comparison/c_comparison/c6_communication.csv")
  
  apCAF5_data  <- rbind(apCAF5_T_data,apCAF5_Mye_data)
  apCAF5_data$X  <- NULL
  unique(apCAF5_data$ligand)
  apCAF5_data<-apCAF5_data[apCAF5_data$ligand %in% c("ANGPTL4", "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DRA","HLA-DRB1","HLA-DRB5","MDK","MIF","SPP1",
  "ADM","C3","CCL3","GRN","HLA-DMA","HLA-DMB","HLA-DQA2","HLA-DQB1","IGF1","NAMPT","THBS1","THBS2","THY1") ,]    
  
  
  
  apCAF6_data  <- rbind(apCAF6_T_data,apCAF6_Mye_data)
  apCAF6_data$X  <- NULL
  apCAF6_data<-apCAF6_data[apCAF6_data$ligand %in% c("ANGPTL4", "HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DRA","HLA-DRB1","HLA-DRB5","MDK","MIF","SPP1",
                                                     "ADM","CADM","C3","GRN","HLA-DMA","HLA-DQA1","HLA-DQB1","IGF1","NAMPT","THBS1","THBS2","THY1") ,]    
  
  
  
  library(ggplot2)
  df <- apCAF5_data
  df$source.target<-as.factor(df$source.target)
  df$source.target<-
    fct_relevel(df$source.target, c("c05 -> CD4Teff (LS)", "c05 -> CD4Teff (NL)","c05 -> CD4Tn (LS)","c05 -> CD4Tn (NL)", 
     "c05 -> Tfh (LS)","c05 -> Tfh (NL)","c05 -> Th17 (LS)","c05 -> Th17 (NL)",   
     "c05 -> Treg (LS)" ,"c05 -> Treg (NL)", 
     "c05 -> CD8Teff (LS)" ,"c05 -> CD8Teff (NL)","c05 -> CD8Tn (LS)","c05 -> CD8Tn (NL)", 
     "c05 -> TIFIT3 (LS)","c05 -> TIFIT3 (NL)","c05 -> TMCM6 (LS)","c05 -> TMCM6 (NL)",       
     "c05 -> TMKI67 (LS)","c05 -> TMKI67 (NL)","c05 -> Tex (LS)", "c05 -> Tex (NL)",             
     "c05 -> CD16+NK (LS)" ,"c05 -> CD16+NK (NL)", "c05 -> CD56+NK (LS)","c05 -> CD56+NK (NL)",
     "c05 -> FMN1+Mac (LS)","c05 -> FMN1+Mac (NL)","c05 -> FOLR2+Mac (LS)","c05 -> FOLR2+Mac (NL)", 
     "c05 -> MARCO+Mac (LS)" ,"c05 -> MARCO+Mac (NL)", 
     "c05 -> CD16+Mo (LS)","c05 -> CD16+Mo (NL)","c05 -> FCN1+Mo (LS)","c05 -> FCN1+Mo (NL)",
     "c05 -> CLEC9A+cDC1s (LS)","c05 -> CLEC9A+cDC1s (NL)", "c05 -> CD1C+cDC2s (LS)","c05 -> CD1C+cDC2s (NL)",    
     "c05 -> LAMP3+mregDC (LS)","c05 -> LAMP3+mregDC (NL)",            
     "c05 -> MKI67+Mye (LS)","c05 -> MKI67+Mye (NL)" ) )      

  dataset.name.order <- levels(df$source.target)
  dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
  dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order) - 1)
  xtick.color <- stringr::str_replace_all(dataset.name.order,c("LS"),c("#DC0000CC"))
  xtick.color <- stringr::str_replace_all(xtick.color,c("NL"),c("#3C5488CC"))
  
  
  
  g<-ggplot(df , aes(x = source.target, y = interaction_name_2, 
                             size= pval, color = prob)) + geom_point(shape = 16,stroke = 4) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, 
           vjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(0.5, length(unique(df$source.target)) + 0.5, 2), lwd = 0.1, colour = "gray80")+
    scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)))+
    scale_color_viridis(option = "C",alpha = 1,discrete = F)+
    theme(axis.text.x = element_text(colour = xtick.color))+
    theme_cxf
  g
  cellchat_apCAF5_NT_data <- g$data
  write.csv( cellchat_apCAF5_NT_data,"cellchat_apCAF5_NT_data.csv")
  dev.off()

  
  
  library(ggplot2)
  df <- apCAF6_data
  df$source.target<-as.factor(df$source.target)
  levels(df$source.target) %in% c("c06 -> CD4Teff (LS)", "c06 -> CD4Teff (NL)","c06 -> CD4Tn (LS)","c06 -> CD4Tn (NL)", 
                                  "c06 -> Tfh (LS)","c06 -> Tfh (NL)","c06 -> Th17 (LS)","c06 -> Th17 (NL)",   
                                  "c06 -> Treg (LS)" ,"c06 -> Treg (NL)", 
                                  "c06 -> CD8Teff (LS)" ,"c06 -> CD8Teff (NL)","c06 -> CD8Tn (LS)","c06 -> CD8Tn (NL)", 
                                  "c06 -> TIFIT3 (LS)","c06 -> TIFIT3 (NL)","c06 -> TMCM6 (LS)","c06 -> TMCM6 (NL)",       
                                  "c06 -> TMKI67 (LS)","c06 -> TMKI67 (NL)","c06 -> Tex (LS)", "c06 -> Tex (NL)",             
                                  "c06 -> CD16+NK (LS)" ,"c06 -> CD16+NK (NL)", "c06 -> CD56+NK (LS)","c06 -> CD56+NK (NL)",
                                  "c06 -> FMN1+Mac (LS)","c06 -> FMN1+Mac (NL)","c06 -> FOLR2+Mac (LS)","c06 -> FOLR2+Mac (NL)", 
                                  "c06 -> MARCO+Mac (LS)" ,"c06 -> MARCO+Mac (NL)", 
                                  "c06 -> CD16+Mo (LS)","c06 -> CD16+Mo (NL)","c06 -> FCN1+Mo (LS)","c06 -> FCN1+Mo (NL)",
                                  "c06 -> CLEC9A+cDC1s (LS)","c06 -> CLEC9A+cDC1s (NL)", "c06 -> CD1C+cDC2s (LS)","c06 -> CD1C+cDC2s (NL)",    
                                  "c06 -> LAMP3+mregDC (LS)","c06 -> LAMP3+mregDC (NL)",            
                                  "c06 -> MKI67+Mye (LS)","c06 -> MKI67+Mye (NL)" )
  df$source.target<-
    fct_relevel(df$source.target,  c("c06 -> CD4Teff (LS)", "c06 -> CD4Teff (NL)","c06 -> CD4Tn (LS)","c06 -> CD4Tn (NL)", 
                                     "c06 -> Tfh (LS)","c06 -> Tfh (NL)","c06 -> Th17 (LS)","c06 -> Th17 (NL)",   
                                     "c06 -> Treg (LS)" ,"c06 -> Treg (NL)", 
                                     "c06 -> CD8Teff (LS)" ,"c06 -> CD8Teff (NL)","c06 -> CD8Tn (LS)","c06 -> CD8Tn (NL)", 
                                     "c06 -> TIFIT3 (LS)","c06 -> TIFIT3 (NL)","c06 -> TMCM6 (LS)","c06 -> TMCM6 (NL)",       
                                     "c06 -> TMKI67 (LS)","c06 -> TMKI67 (NL)","c06 -> Tex (LS)", "c06 -> Tex (NL)",             
                                     "c06 -> CD16+NK (LS)" ,"c06 -> CD16+NK (NL)", "c06 -> CD56+NK (LS)","c06 -> CD56+NK (NL)",
                                     "c06 -> FMN1+Mac (LS)","c06 -> FMN1+Mac (NL)","c06 -> FOLR2+Mac (LS)","c06 -> FOLR2+Mac (NL)", 
                                     "c06 -> MARCO+Mac (LS)" ,"c06 -> MARCO+Mac (NL)", 
                                     "c06 -> CD16+Mo (LS)","c06 -> CD16+Mo (NL)","c06 -> FCN1+Mo (LS)","c06 -> FCN1+Mo (NL)",
                                     "c06 -> CLEC9A+cDC1s (LS)","c06 -> CLEC9A+cDC1s (NL)", "c06 -> CD1C+cDC2s (LS)","c06 -> CD1C+cDC2s (NL)",    
                                     "c06 -> LAMP3+mregDC (LS)","c06 -> LAMP3+mregDC (NL)",            
                                     "c06 -> MKI67+Mye (LS)","c06 -> MKI67+Mye (NL)" ) )      
  
  dataset.name.order <- levels(df$source.target)
  dataset.name.order <- stringr::str_match(dataset.name.order, "\\(.*\\)")
  dataset.name.order <- stringr::str_sub(dataset.name.order, 2, stringr::str_length(dataset.name.order) - 1)
  xtick.color <- stringr::str_replace_all(dataset.name.order,c("LS"),c("#DC0000CC"))
  xtick.color <- stringr::str_replace_all(xtick.color,c("NL"),c("#3C5488CC"))
  
  
  
  g<-ggplot(df, aes(x = source.target, y = interaction_name_2, 
                     size= pval, color = prob)) + geom_point(shape = 16,stroke = 4) + 
    theme_linedraw() + theme(panel.grid.major = element_blank()) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                     vjust = 1), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
    scale_x_discrete(position = "bottom") +
    geom_vline(xintercept = seq(0.5, length(unique(df$source.target)) + 0.5, 2), lwd = 0.1, colour = "gray80")+
    scale_radius(range = c(min(df$pval), max(df$pval)), breaks = sort(unique(df$pval)))+
    scale_color_viridis(option = "C",alpha = 1,discrete = F)+
    theme(axis.text.x = element_text(colour = xtick.color))+
    theme_cxf
  g
  cellchat_apCAF6_NT_data <- g$data
  write.csv( cellchat_apCAF5_NT_data,"./result/new_CAF/cellchat_apCAF6_NT_data.csv")
  dev.off()
  
}



cellphonedb{
  
data{
  d_CAF <-readRDS("./data/d_Fb7.rds")
  table(d_CAF$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_CAF@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=600)  
    samp2<-df %>% group_by(cluster) %>% sample_n(size=600)  
    samp3<-df %>% group_by(cluster) %>% sample_n(size=600)  
    S<-rbind(samp,samp2,samp3)
    S<-S[!duplicated(S$barcode),]
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    #samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_CAF_subset<-d_CAF[,S$barcode]
  }
  table(d_CAF_subset$seurat_clusters)
    
  if(T){
    library(dplyr)
    table(CAF_T$seurat_clusters)
    B<-levels(Idents(CAF_T))[14:26]
    d_T_subset<-subset(CAF_T,idents = B)
  }
  table(d_T_subset$seurat_clusters)
  
  d_Mye <- readRDS("./data/d_Mye.rds")
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_Mye@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    samp<-df %>% group_by(cluster) %>% sample_n(size=1000)  
    #samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    #samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.10)  
    d_Mye_subset<-d_Mye[,samp$barcode]
  }
  table(d_Mye_subset$seurat_clusters)
  
  d_Mas <- readRDS("./data/d_Mast.rds")
  table(d_Mas$seurat_clusters)
  if(T){
    library(dplyr)
    meta<-as.data.frame(d_Mas@meta.data)
    samp<-sample(rownames(meta),size=500)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    #samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.15)  
    d_Mas_subset<-d_Mas[,samp]
  }
  table(d_Mas_subset$seurat_clusters)
  rm(d_Mye,d_CAF,d_Mas)
  
  # Save normalised counts - NOT scaled!
  CAF_all<-merge(x = d_CAF_subset, y = c(d_T_subset,d_Mye_subset,d_Mas_subset))
  CAF_all<- NormalizeData(CAF_all)
  CAF_all$cell_type<-CAF_all$seurat_clusters
  
  #writeMM(CAF_T@assays$RNA@data, file = 'CAF_T_matrix.mtx',row.names=F)
  library(SeuratDisk)
  SaveH5Seurat(CAF_all, filename = "CAF_all.h5Seurat")
  Convert("CAF_all.h5Seurat", dest = "h5ad")
  
  #save gene and cell names
  write(x = rownames( CAF_all@assays$RNA@data), file = "CAF_all_features.tsv")
  write(x = colnames( CAF_all@assays$RNA@data), file = "CAF_all_barcodes.tsv")
  

  table(CAF_all@meta.data$cell_type)
  CAF_all@meta.data$Cell = rownames(CAF_all@meta.data)
  df = CAF_all@meta.data[, c('Cell', 'cell_type')]
  write.table(df, file ='CAF_all_meta.tsv', sep = '\t', quote = F, row.names = F)
  saveRDS(CAF_all,"CAF_all.rds")

  ## OPTION 1 - compute DEGs for all cell types
  ## Extract DEGs for each cell_type
  # Idents(CAF_T)<-"cell_type"
  # DEGs <- FindAllMarkers(CAF_T,
  #                        test.use = "wilcox", #"LR" 
  #                        verbose = F,
  #                        only.pos = T,
  #                        random.seed = 1,
  #                        logfc.threshold = 0.2,
  #                        min.pct = 0.1,
  #                        return.thresh = 0.05)
  
  # OPTION 2 - optional - Re-compute  hierarchical (per lineage) DEGs for Epithelial and Stromal lineages
  # CAF_T$lineage<-CAF_T$cell_type
  # levels(CAF_T$lineage)<-c("C","C","C","C","C","C","C","C","C","C","C","C","C","T","T","T","T","T" ,"T","T","T","T","T","T","T","T")
  # 
  # DEGs = c()
  # for( lin in c('C', 'T') ){
  #   message('Computing DEGs within linage ', lin)
  #   so_in_lineage = subset(CAF_T, cells = Cells(CAF_T)[CAF_T$lineage == lin ])
  #   Idents( so_in_lineage)<-"cell_type"
  #   DEGs_lin = FindAllMarkers(so_in_lineage,
  #                             verbose = F,
  #                             only.pos = T,
  #                             random.seed = 1,
  #                             logfc.threshold = 0.2,
  #                             min.pct = 0.1,
  #                             return.thresh = 0.05)
  #   DEGs = rbind(DEGs_lin, DEGs)
  # }
  # colnames(DEGs)
  # fDEGs = subset(DEGs, p_val_adj < 0.05 & avg_log2FC >= 0.2)
  # # 1st column = cluster; 2nd column = gene 
  # fDEGs = fDEGs[, c('cluster', 'gene', 'p_val_adj', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2')] 
  # write.table(fDEGs, file ='CAF_T_DEGs.tsv', sep = '\t', quote = F, row.names = F)
  # 
  # head(fDEGs)
  # 
}  
  
  
  
  library(Seurat)
  library(SingleCellExperiment)
  library(reticulate)
  # ad=import('anndata')
  # adata = ad$read_h5ad('CAF_T.h5ad')
  # counts <- Matrix::t(adata$X)
  # row.names(counts) <- row.names(adata$var)
  # colnames(counts) <- row.names(adata$obs)
  # CAF_T.sce <- SingleCellExperiment(list(counts = counts), colData = adata$obs, rowData = adata$var)
  CAF_all<-readRDS("CAF_all.rds")
  CAF_all.sce <- as.SingleCellExperiment(CAF_all)
  pvals <- read.delim("./cpdb_results/pvalues.txt", check.names = FALSE)
  means <- read.delim("./cpdb_results/means.txt", check.names = FALSE)
  decon <- read.delim("./cpdb_results/deconvoluted_percents.txt", check.names = FALSE)
  inscore<-read.delim("./cpdb_results/interaction_scores.txt", check.names = FALSE)
  CAF_all$cell_type <- as.factor(CAF_all$cell_type)
  levels(CAF_all$cell_type)
  
  X=c("c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13",  
      "CD16+Mo" ,"FCN1+Mo",
      "FMN1+Mac","FOLR2+Mac","MARCO+Mac",
      "CD1C+cDC2s","CLEC9A+cDC1s","LAMP3+mregDC",
      "MKI67+Mye",
      "CD8Teff","CD8Tn", "Tex","TIFIT3","TMCM6","TMKI67",
      "CD4Teff","CD4Tn","Tfh","Th17","Treg",
      "CD16+NK","CD56+NK",
      "Mast cell")
  CAF_all$cell_type<-fct_relevel(CAF_all$cell_type,X)
  levels(CAF_all$cell_type)
  
  A<-X[1:13]
  B<-X[14:36]
  
  
ktplots{  
  #PLOT
  # if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
  # BiocManager::install(version = "3.18")
  # BiocManager::install("ComplexHeatmap")
  # usethis::use_git_config(user.name = "windyeeeee", user.email = "chenxiongfengfj@qq.com")
  # credentials::set_github_pat()
  # devtools::install_github('zktuong/ktplots', dependencies = TRUE)
  library(ggplot2)
  library(ktplots)

  geneDotPlot(scdata = CAF_all.sce, # object
              genes = c("APP", "SORL1", "CD74"), # genes to plot
              celltype_key = "cell_type", # column name in meta data that holds the cell-cluster ID/assignment
              standard_scale = TRUE) + # whether to scale expression values from 0 to 1. See ?geneDotPlot for other options
             theme(strip.text.x = element_text(angle=45, hjust = 0, size =8))# + small_guide() + small_legend()
  
  p1<-plot_cpdb(scdata = CAF_all.sce, 
            cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
            cell_type2 = 'c01', 
            celltype_key = 'cell_type', means = means, pvals = pvals,
            interaction_scores = inscore,
            genes = c("TNFSF14","TNFSF12","TNF","TGFB3","TGFB2","TGFB1",
                        "PARRES2","LTB","LTA","JAG1","IFNR","CXCL8","CXCL3","CXCL2","CXCL12",
                       "CD99","CD93","CD55","CD47","CD44","CD40LG","CD320","CD160","BTLA","AREG","APP")
            #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
            #splitby_key = 'Experiment', gene_family = 'chemokines'
             ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  dev.off()
  
  p2<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c02', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("TNFSF14","TNFSF12","TNFSF10","TNF","TGFB3","TGFB2","TGFB1",
                          "PARRES2","PI16","LTB","LTA","JAG1","IFNR","CXCL8","CXCL3","CXCL2","CXCL12","CXCL1",
                          "CD99","CD93","CD55","CD47","CD44","CD40LG","CD34","CD248","CD320","CD160","BTLA","AREG","APP")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p2
  dev.off()
  
  p3<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c03', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("WNT5B","WNT2B","TNC","THY1","TGFB3","TGFB2","TGFB1","TGFA","THBS1","THBS2",
                           "SLTT2","SEMA7A","SEMA4D","SEMA4C","SEMA4A","SEMA3C","SDC","RBP4","PPIA","PRNP",
                        "CLEC2B",  #"CXCL1",  "CXCL2",  "CXCL3",  "CXCL8",  "CXCL10", "CXCL12",
                          "COL9A2",  "COL9A2","CD58","CD1B","C3","BMP4","BAG6","APOE","ANXA","APP")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p3
  dev.off()
  
  p4<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c04', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("VEGFB","VEGFA","TENM4","PGF","PDGFB","LGALS9","GAS6","OSM","NTN4","NTF",
                          "SLTT2","SEMA7A","SEMA4D","SEMA4C","SEMA4A","SEMA3C","SDC","RBP4","PPIA","PRNP",
                          "NCAM1",  #"CXCL1",  "CXCL2",  "CXCL3",  "CXCL8",  "CXCL10", "CXCL12",
                          "CSF1",  "HBEGF","HEBP1","LGALS9","LGALS3")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p4
  dev.off()
  
  p5<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c05', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("TGM2","TGFA","SPP1","PTPRC","PLAU","PGF","PDGFC","PDGFB","MPZL1","LAIR1",
                          "IL1B","IL18","IGFB3","IGF1","GAS6","EREG","COL4A1","COL4A2","CD86","CD34",
                          "CCL5")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p5
  dev.off() 
  
  p6<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c06', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("SPP1","PGF","PDGFC","PDGFB","OSM","NRG1","NCM1","MPZ1","HCF","GAS6",
                          "FASLG","EFNA1","EFNA5","DSC2","CXCL9","COL9A2","COL9A3","CLSTN3","CLSTN1","C3",
                          "CCL22")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p6
  dev.off() 
  
  
  p7<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c07', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                 genes = c("VEGFA","VEGFB","TNFSF14","TNFSF12","TNFSF10","TNF","TGFB3","TGFB2","TGFB1","TGFA",
                       "PARRES2","PGF","PLA2G2A","PDGFC","PDGFB","LGALS9","LTB","LTA","JAG1","IFNR",
                "JAM2","JAM3","JAG1","ICAM1","CXCL3","CXCL2","CXCL12","CXCL1",
                "CD99","CD93","CD55","CD47","CD44","CD40LG","CD34","CD248","CD320","CD160","BTLA","AREG","APOE","APP")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p7
  dev.off()

  p8<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c08', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("WNT5B","WNT2B","VEGFA","VEGFB","TSLP","THY1","TGFB3","TGFB2","TGFB1","TGFA","THBS1","THBS2","TGM2",
                          "SLTT2","SEMA7A","SEMA4D","SEMA4C","SEMA4A","SDC2","RARRES2","PGF","PPIA","PRNP",
                          "TENM4",  "LAMA2",  "LAMC1",  "LAMC3",  "FLT3LG",  "FTH1", "FTL",
                          "DLL1",  "CXCL1","CXCL2","CXCL14","COL4A6","COL3A1","CCL11","BMP4","BMP5","BMP2","APP","APOE")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p8
  dev.off()
  
  p9<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c09', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("IGFBP3","IL6","IGF1","IFNR","SAA1","PRNP","OSM",
                          "MPZL","LRPAP1","LPAR1","LGALS9","IFNR",
                          "EFNB2","CXCL8","CXCL6","CSF3","CXCL3","CXCL2","CXCL12","CXCL1",
                           "C3","BAG6")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p9
  dev.off()
  
  p10<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c10', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("IGFBP3","IL6","IGF1","IFNR","SAA1","PRNP","OSM","HEBP1","HBEGF","GAS6",
                          "MPZL","LRPAP1","LPAR1","LGALS9","IFNR","EREG",
                          "EFNB2","CXCL8","CXCL6","CSF3","CXCL3","CXCL2","CXCL12","CXCL1",
                          "C3","BAG6")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p10
  dev.off()
  
  
  p11<-plot_cpdb(scdata = CAF_all.sce, 
                 cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                 cell_type2 = 'c11', 
                 celltype_key = 'cell_type', means = means, pvals = pvals,
                 interaction_scores = inscore,
                  genes = c("VEGFA","VEGFB","VEGFA","TNC","THY1","TGFB3","SLIT2","SDC2","PPIA",
                 "SEMA7A","SEMA4D","SEMA4C","SEMA4A","SDC2","RARRES2","PGF","PPIA","PRNP",
                 "TENM4",  "KLRB1",  "KITLG",  "JAG1",  "IL1RAP",  "CXCL14",  "CXCL12","CNTN1","CADM1","BAG6")
                 #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                 #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p11
  dev.off()
  
  
  p12<-plot_cpdb(scdata = CAF_all.sce, 
                cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                cell_type2 = 'c12', 
                celltype_key = 'cell_type', means = means, pvals = pvals,
                interaction_scores = inscore,
                genes = c("VEGFB","VEGFA","TNC","THY1","PDGFB","NCAM1","LAMC","JAG1","HEBP1","HBEGF",
                          "SLTT2","SEMA7A","SEMA4D","SEMA4C","SEMA4A","SEMA3C","SDC","ENTPD1","PPIA","PRNP",
                          "NCAM1",  "CLSTN3",  "CLSTN1",  "CD47")
                #gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p12
  dev.off()
  
  p13<-plot_cpdb(scdata = CAF_all.sce, 
                 cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
                 cell_type2 = 'c13', 
                 celltype_key = 'cell_type', means = means, pvals = pvals,
                 interaction_scores = inscore,
                 # genes = c("MICA","VEGFA","TNC","THY1","PDGFB","NCAM1","LAMC","JAG1","HEBP1","HBEGF",
                 #           "SLTT2","SEMA7A","SEMA4D","SEMA4C","SEMA4A","SEMA3C","SDC","ENTPD1","PPIA","PRNP",
                 #           "NCAM1",  "CLSTN3",  "CLSTN1",  "CD47")
                 gene_family = c( 'chemokines', 'Th1', 'Th2', 'Th17', 'Treg', 'costimulatory', 'coinhibitory', 'niche')
                 #splitby_key = 'Experiment', gene_family = 'chemokines'
  ) #+ small_guide() + small_axis() + small_legend(keysize=1)+theme_classic()
  p13
  dev.off()
  
  
  
  # test <- plot_cpdb2(cell_type1 = 'CD16+Mo|FCN1+Mo|FMN1+Mac|FOLR2+Mac|MARCO+Mac|CD1C+cDC2s|CLEC9A+cDC1s|LAMP3+mregDC|MKI67+Mye|CD8Teff|CD8Tn|Tex|TIFIT3|TMCM6|TMKI67|CD4Teff|CD4Tn|Tfh|Th17|Treg|CD16+NK|CD56+NK|Mast cell', 
  #                    cell_type2 = 'c01', 
  #                    idents ='cell_type',
  #                    #split.by = 'treatment_group_1',
  #                    scdata = CAF_all.sce,
  #                    means = means,
  #                    pvals = pvals,
  #                    deconvoluted = deconvoluted, # new options from here on specific to plot_cpdb2
  #                    #gene_symbol_mapping = 'index', # column name in rowData holding the actual gene symbols if the row names is ENSG Ids. Might be a bit buggy
  #                    #desiredInteractions = list(c('CD4_Tcm', 'cDC1'), c('CD4_Tcm', 'cDC2'), c('CD4_Tem', 'cDC1'), c('CD4_Tem', 'cDC2 '), c('CD4_Treg', 'cDC1'), c('CD4_Treg', 'cDC2')),
  #                    #interaction_grouping = interaction_grouping,
  #                    #edge_group_colors = c("Activating" = "#e15759", "Chemotaxis" = "#59a14f", "Inhibitory" = "#4e79a7", "   Intracellular trafficking" = "#9c755f", "DC_development" = "#B07aa1"),
  #                    node_group_colors =c("CD16+Mo"= "#C8E6C9" ,"FCN1+Mo"= "#FFE082", "FMN1+Mac"=  "#64B5F6", "FOLR2+Mac" = "#F0F4C3" ,"MARCO+Mac" = "#E57373"  ,  
  #                                         "CD1C+cDC2s"= "#DCE775" , "CLEC9A+cDC1s" = "#90CAF9", "LAMP3+mregDC" = "#26C6DA","MKI67+Mye"= "#C0CA33" ,  "CD8Teff"= "#F44336" ,    
  #                                         "CD8Tn"=  "#80DEEA"  , "Tex"  = "#EF6C00" , "TIFIT3" = "#EF9A9A" , "TMCM6" = "#FFF176" ,     
  #                                         "TMKI67" =  "#43A047" , "CD4Teff"=  "#FFB74D" , "CD4Tn"=  "#66BB6A" , "Tfh" ="#FB8C00" ,"Th17" =  "#B2EBF2" ,      
  #                                         "Treg" = "#2E7D32" , "CD16+NK" = "#FFCC80", "CD56+NK" =  "#E0F7FA"  , "Mast cell"= "#FFEBEE","c01"= "orange3"),
  #                    keep_significant_only = TRUE,
  #                    standard_scale = TRUE,
  #                    remove_self = TRUE)
  # 


 all_color<-c("#C8E6C9","#FFE082" ,"#64B5F6" ,"#F0F4C3", "#E57373" ,"#DCE775", "#90CAF9" ,"#26C6DA" , "#C0CA33", "#F44336","#80DEEA", "#EF6C00", 
               "#EF9A9A", "#FFF176", "#43A047" ,"#FFB74D" ,"#66BB6A", "#FB8C00" ,"#B2EBF2" ,"#2E7D32" ,"#FFCC80" ,"#E0F7FA","#FFEBEE") 
 c_color<-"orange"
 
p_data<-p1$data
write.table(p_data,"c01_data.txt")
p_data<-p2$data
write.table(p_data,"c02_data.txt")
p_data<-p3$data
write.table(p_data,"c03_data.txt")
p_data<-p4$data
write.table(p_data,"c04_data.txt")
p_data<-p5$data
write.table(p_data,"c05_data.txt")
p_data<-p6$data
write.table(p_data,"c06_data.txt")
p_data<-p7$data
write.table(p_data,"c07_data.txt")
p_data<-p8$data
write.table(p_data,"c08_data.txt")
p_data<-p9$data
write.table(p_data,"c09_data.txt")
p_data<-p10$data
write.table(p_data,"c10_data.txt")
p_data<-p11$data
write.table(p_data,"c11_data.txt")
p_data<-p12$data
write.table(p_data,"c12_data.txt")
p_data<-p13$data
write.table(p_data,"c13_data.txt")


  
cell_color <- c("CD16_Mo"= "#C8E6C9" ,"FCN1_Mo"=  "#66BB6A" , "FMN1_Mac"= "#DCE775" , "FOLR2_Mac" = "#F0F4C3" ,"MARCO_Mac" = "#C0CA33" ,  
   "CD1C_cDC2s"= "#FFE082" , "CLEC9A_cDC1s" = "#FFCC80", "LAMP3_mregDC" = "#FFB74D" ,"MKI67_Mye"= "#FB8C00" ,  "CD8Teff"= "#F44336" ,
   "CD8Tn"= "orange3"    , "Tex"  = "#EF6C00" , "TIFIT3" = "#EF9A9A" , "TMCM6" = "#FFF176" ,
   "TMKI67" = "#E57373"   , "CD4Teff"= "#64B5F6", "CD4Tn"="#80DEEA" , "Tfh" ="#E0F7FA" ,"Th17" =  "#B2EBF2" ,"Treg" ="#90CAF9"  , 
   "CD16_NK" ="#2E7D32", "CD56_NK" = "#43A047"   , "Mast cell"= "#FFEBEE",
   "c13"="orange")

  
  
}   
  
CCPlotR{
  # usethis::use_git_config(user.name = "windyeeeee", user.email = "chenxiongfengfj@qq.com")
  # credentials::set_github_pat()
  # devtools::install_github("Sarah145/CCPlotR")
  library(CCPlotR)
  data(toy_data, toy_exp, package = 'CCPlotR')
  
  colnames(p_data)
  cc_data<-p_data[,c( "Var1","Var2","scaled_means","pvals","significant")]
  rownames(cc_data) <- paste0(cc_data$Var1,"-",cc_data$Var2)
  colnames(cc_data) <- c("ligand_receptor","source_target" , "scaled_means","pvals","significant")
 
  cc_data$source <- str_split(cc_data$source_target,"-", simplify = T)[,1]
  cc_data$target <- str_split(cc_data$source_target,"-", simplify = T)[,2]
  
  cc_data$ligand<-str_split(cc_data$ligand_receptor,"-", simplify = T)[,1]
  cc_data$receptor<-str_split(cc_data$ligand_receptor,"-", simplify = T)[,2]
  
  cc_data$score <-  cc_data$scaled_means
  length(c(unique(cc_data$ligand),unique(cc_data$receptor)))
  cc_exp<-decon[,c("gene_name",A[13],B)]
  cc_exp<-cc_exp[cc_exp$gene_name %in%  c(unique(cc_data$ligand),unique(cc_data$receptor)),]
  cc_exp <- reshape2::melt(data = cc_exp, di.var = 'gene_name')
  cc_exp <- cc_exp[,c("variable","gene_name","value")]
  colnames(cc_exp)<-colnames(toy_exp)
  cc_exp <- cc_exp[!duplicated(cc_exp[,c("cell_type","gene")]),] 
  class(cc_exp )
  table(cc_exp$cell_type)

  cc_data$source_target<-gsub("\\+","_",cc_data$source_target)
  cc_data$source<-gsub("\\+","_",cc_data$source)
  cc_data$target<-gsub("\\+","_",cc_data$target)
  cc_exp$cell_type<-gsub("\\+","_",cc_exp$cell_type)
  
  cc_data<-cc_data[cc_data$source=="c13",]
  write.table(cc_data,"c13_ccdata.txt")
  
  cc_data<-read.table("c02_ccdata.txt")
  
  dev.off()
  cc_dotplot(cc_data)
  cc_dotplot(cc_data, option = 'B', n_top_ints = 10)
  cc_dotplot(cc_data, option = 'Liana', n_top_ints = 15)
  
  cc_network(cc_data)
  cc_network(cc_data, #great
             colours = c(all_color,c_color), 
             option = 'B',
             n_top_ints = 50,
             node_size = 1.75,
             label_size = 4,
             layout = "kk")  
  
  cc_circos(cc_data)
  cc_circos(cc_data, option = 'B', n_top_ints = 15)
  cc_circos(cc_data, 
            option = 'B', 
            n_top_ints = 20,
            exp_df =cc_exp,
            cell_cols =
              c(
               `c01` = 'hotpink',
               `CD16_NK` = "#FFAB91",
               `CD4Teff` = "#FFCC80",
              `CD4Tn` =  "#FFE082",
              `CD56_NK` = "#FFF59D",
              `CD8Teff` = "#AED581",
              `CD8Tn` =  "#A5D6A7",
              `TIFIT3` = "#80CBC4",
              `TMCM6` =    "#26C6DA",
              `TMKI67` =  "#29B6F6",
              `Tex` = "#42A5F5",
              `Tfh` = "#7986CB",
              `Th17` = "#9575CD",
             `Treg` =  "#BA68C8"),
            show_legend = TRUE)
  cc_circos(cc_data,                   #great
            option = 'C', 
            n_top_ints = 50, #dim(cc_data)
            exp_df =cc_exp,
            cell_cols =cell_color,
            # cell_cols = c(
            #   `c01` = 'hotpink',
            #   `CD16_NK` = "#FFAB91",
            #   `CD4Teff` = "#FFCC80",
            #   `CD4Tn` =  "#FFE082",
            #   `CD56_NK` = "#FFF59D",
            #   `CD8Teff` = "#AED581",
            #   `CD8Tn` =  "#A5D6A7",
            #   `TIFIT3` = "#80CBC4",
            #   `TMCM6` =    "#26C6DA",
            #   `TMKI67` =  "#29B6F6",
            #   `Tex` = "#42A5F5",
            #   `Tfh` = "#7986CB",
            #   `Th17` = "#9575CD",
            #   `Treg` =  "#BA68C8"),
            show_legend = TRUE)
  
  dev.off()
  
  cc_sigmoid(cc_data,#great
             colours = cell_color,
             n_top_ints = 25) 
  
  
  cc_arrow(cc_data, cell_types = c(A[1],B[1]), colours =  c(
        `c01` = 'hotpink',
         `CD8Teff` = "#AED581"))
  cc_arrow(cc_data, cell_types = c('c01', 'CD8Teff'), option = 'B', exp_df = cc_exp, n_top_ints = 10, palette = 'OrRd')

  
}
  
}




CAF_Manhattan{
  rm(list=ls())
  setwd("./../")
  getwd()
  seurat_Fb <- readRDS("./data/d_Fb7.rds")
  dim(seurat_Fb)
  Idents(seurat_Fb)
  
  cl <- c("c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13")
  for (i in cl) {
  seurat_CAF <-subset(seurat_Fb, idents = i)
  Idents(seurat_CAF)<-"Group"
  degdf <- FindMarkers(seurat_CAF,ident.1 = "Tumor",ident.2 = "Normal", 
                       logfc.threshold = 0.01)
  logFC_t=0
  P.Value_t = 1e-10
  #"p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj" 
  degdf$group = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < -1,"down",
                       ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 1,"up","stable"))
  degdf$name=rownames(degdf)
  degdf$"-log10(p_val_adj)"<- -log10(degdf$p_val_adj+1e-300)
  x<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up",])[1]
  y<-dim(degdf[degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down",])[1]
  z<-dim(degdf)[1]-x-y
  degdf$"-log10(p_val_adj)"<- ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "up", 
                                     sample(runif(x,10,300),x),
                                     ifelse(degdf$"-log10(p_val_adj)" == 300 & degdf$group == "down", 
                                            sample(runif(y,10,300),y),      
                                            degdf$"-log10(p_val_adj)"))
  degdf$pt<-log(degdf$pct.1/degdf$pct.2+1)
  degdf$gene <- rownames(degdf)
  degdf$label <- rownames(degdf)
  degdf$label [degdf$group == "stable"] <- NA
  write_csv(degdf,paste0("./result/new_CAF/CAF_DEGs/",i,".csv"))
  assign(i,degdf)
  }
c01$cluster<-"c01"
c02$cluster<-"c02"
c03$cluster<-"c03"
c04$cluster<-"c04"
c05$cluster<-"c05"
c06$cluster<-"c06"
c07$cluster<-"c07"
c08$cluster<-"c08"
c09$cluster<-"c09"
c10$cluster<-"c10"
c11$cluster<-"c11"
c12$cluster<-"c12"
c13$cluster<-"c13"
call<-rbind(c01,c02,c03,c04,c05,c06,c07,c08,c09,c10,c11,c12,c13)
write_csv(call,"./result/new_CAF/CAF_DEGs/c_all.csv")
call2<- na.omit(call)
write_csv(call2,"./result/new_CAF/CAF_DEGs/c_up_down.csv") 
install.packages("qqman")
library(qqman)
call2$BP<-1:length(call2$avg_log2FC)
colnames(call2)
call3<-call2[,c("label","cluster","BP","avg_log2FC")]
call3$CHR<-as.factor(call3$CHR)
levels(call3$CHR)<-1:13
colnames(call3)<-c("SNP","CHR","BP","P")
call3$CHR <- as.numeric(call3$CHR)
manhattan(call3,logp = F, suggestiveline = 1, annotatePval = 1.5,
          chrlabs = c("c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13"),
          main = "DEG of CAF clusters", 
          cex = 1.2,  
          cex.axis = 1.5,  
          col = c("blue4", "orange3"),
          highlight = NULL,
          annotateTop = F)

call4<-call3
call4$P<--call4$P
manhattan(call4,logp = F, suggestiveline = 1, annotatePval = 1.5,
          chrlabs = c("c01","c02","c03","c04","c05","c06","c07","c08","c09","c10","c11","c12","c13"),
          main = "DEG of CAF clusters", 
          cex = 1.2,  
          cex.axis = 1.5,  
          col = c("blue4", "orange3"),
          highlight = NULL,
          annotateTop = F)
}



Pan_scRNAtoolVis{
  #seurat_Fb <- readRDS("./data/d_Fb7.rds")
  #tmp <- seurat_Fb
  # c_all <- read_csv("result/new_CAF/CAF_DEGs/c_all.csv")
  # install.packages('devtools')
  # devtools::install_github('junjunlab/scRNAtoolVis')
  # if not install ggunchull
  # devtools::install_github("sajuukLyu/ggunchull", type = "source")
  tmp <- readRDS("./data/d_seurat_A.rds")
  levels(tmp$Type)
  table(tmp$Type)
  tmp$Type<-fct_relevel(tmp$Type,c("Epithelial cell","T cell","Fibroblasts cell", "Myeloid cell" ,"Endothelial cell", "Plasma cell", "B cell","Mast cell"))
  
  levels(tmp$seurat_clusters)
  tmp$seurat_clusters<-fct_relevel(tmp$seurat_clusters,
                                   c( "C1" , "C2",  "C3" , "C4" , "C5" , "C6" , "C7" , "C8" , "C9" , "C10" ,"C11" ,"C12", "C13" ,"C14" ,"C15" ,"C16" ,"C17" ,
                                      "T1" , "T2" , "T3" , "T4" , "T5" , "T6" , "T7" , "T8" ,
                                      "F1" , "F2",  "F3" , "F4" , "F5" , "F6",
                                      "D1", "D2" , "D3" , "D4" , "D5" , "D6" ,
                                      "E1" , "E2" , "E3" , "E4",
                                      "P1" , "P2" , "P3" , "P4" , "P5" , "P6",  "P7",
                                      "B1" , "M1"))
  
  #B cell -1 colors_list[[5]]
  BC <- c("#00897B")
  #Epithelial cell -17
  CC <- rev(c("#FFFDE7","#FFF9C4","#FFF59D","#FFF176","#FFEE58",
              "#FFE082","#FFD54F","#FFCA28","#FFC107","#FFB300",
              "#FFCC80","#FFB74D","#FFA726","#FF9800","#FB8C00","#F57C00","#EF6C00"))
  #MONO-6         colors_list[[9]]
  DC <- rev(c("#E0F7FA","#B2EBF2","#80DEEA","#4DD0E1","#26C6DA","#00BCD4"))
  #Endothlial-4   colors_list[[11]]
  EC <- rev(c("#BBDEFB","#90CAF9","#64B5F6","#42A5F5"))
  #Fibroblasts-6   colors_list[[14]]  
  FC <- rev(c("#FFEBEE","#FFCDD2","#EF9A9A","#E57373","#EF5350","#F44336"))
  #Mast-1
  MC <- c("#C62828")
  #Plasma-7  colors_list[[7]]
  PC <- rev(c("#C8E6C9","#A5D6A7","#66BB6A","#43A047","#388E3C","#2E7D32","#1B5E20"))
  #T-8
  TC <- rev(c("#F9FBE7","#F0F4C3","#E6EE9C","#DCE775","#D4E157","#CDDC39","#C0CA33","#AFB42B" ))
  
  col<-c(BC,CC,DC,EC,FC,PC,MC,TC)
  
  
  color_umap <- c(CC,TC,FC,DC,EC,PC,BC,MC)
  color_Type<-c("#FFA726","#CDDC39","#F44336","#26C6DA","#64B5F6","#388E3C","#00897B","#C62828")

  
  if(T){
    library(dplyr)
    meta<-as.data.frame(tmp@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    #samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.01)  
    tmp_subset<-tmp[,samp$barcode]
  }
  
  all.markers <- FindAllMarkers(tmp_subset, only.pos = FALSE,
                                min.pct = 0.25,
                                logfc.threshold = 0)

  
  Idents(tmp_subset)
  tmp<-tmp_subset
  tmp<-subset(tmp , idents=  CO)
  tmp$seurat_clusters<-droplevels(tmp$seurat_clusters)
  CO<-c( "C1" , "C2",  "C3" , "C4" , "C5" , "C6" , "C7" , "C8" , "C9" , "C10" ,"C11" ,"C12", "C13" ,"C14" ,"C15" ,#"C16" ,"C17" ,
     "T1" , "T2" , "T3" , "T4" , "T5" , "T6" ,  "T8" ,#"T7" ,
     "F1" , "F2",  "F3" , "F4" , "F5" , "F6",
     "D1", "D2" , "D3" , "D4" , "D5" , "D6" ,
     "E1" , "E2" , "E3" , "E4",
     "P1" , "P2" , "P3" , "P4" , "P5" ,# "P6",  "P7",
     "B1" , "M1")
  
  color_umap2 <- color_umap[c(1:15,18:46,49:50)]  #[1] "#EF6C00" "#F57C00" "#FB8C00" "#FF9800" "#FFA726" 
  names(color_umap2)<-CO
  all.markers$cluster<-fct_relevel(all.markers$cluster,
                                   c( "C1" , "C2",  "C3" , "C4" , "C5" , "C6" , "C7" , "C8" , "C9" , "C10" ,"C11" ,"C12", "C13" ,"C14" ,"C15" ,"C16" ,"C17" ,
                                      "T1" , "T2" , "T3" , "T4" , "T5" , "T6" , "T7" , "T8" ,
                                      "F1" , "F2",  "F3" , "F4" , "F5" , "F6",
                                      "D1", "D2" , "D3" , "D4" , "D5" , "D6" ,
                                      "E1" , "E2" , "E3" , "E4",
                                      "P1" , "P2" , "P3" , "P4" , "P5" , "P6",  "P7",
                                      "B1" , "M1"))
  levels(all.markers$cluster) 
  c_all<-all.markers
  
  write.table(all.markers,"all_markers.txt")
  saveRDS(tmp,"tmp_subset.rds")
  
  
scRNAtoolVis{
  library(scRNAtoolVis)
 
  # facet by metadata column "orig.ident"
  clusterCornerAxes(object = tmp,reduction = 'umap',
                    noSplit = F,groupFacet = 'Group',
                    relLength = 0.5)
  dev.off()
  # umap
  clusterCornerAxes(object = tmp,reduction = 'umap',
                    noSplit = T)
  
  # gene
  FeatureCornerAxes(object = tmp,reduction = 'umap',
                    groupFacet = 'Group',
                    relLength = 0.5,relDist = 0.2,
                    features = c("CD74","SPP1", "PDPN", "FAP"))
 
  #AverageHeatmap
 
 
   
   genes<-c( "CD79A","CD19","BLK","VPREB3","SPIB",       #B B cell           
             "EPCAM","CDH1","KRT5","KRT18","KRT19",       #C  Epithelial cell       CFHR3
             "PTPRC","CD14", "AIF1","TYROBP","CD163",     #D  MONO                  ITGAD
             "PECAM1","VWF","VEGFA","CDH5","TEK",        #E  Endo 
             "COL1A1","COL1A2","DCN","LUM","RGS5",         #F  Fibroblasts
             "TPSB2","TPSAB1","MS4A2","CTSG","HPGD",        #M  Mast         TPSD1   MS4A2  HDC   CTSG
             "MS4A1","MZB1","BRSK1","JCHAIN","IGKC",       #P  Plasma  
             "CD3D","CD3E","CD4","CD8A","CD8B")            #T  T
   
   genes<- c("EPCAM","CDH1","KRT5","KRT18","KRT19",         #C  Epithelial cell       CFHR3
             "CD3D","CD3E","CD7","CCL5","CST7",             #T  T
             "COL1A1","COL1A2","ACTA2","TAGLN","DCN",       #F  Fibroblasts
             "AIF1","CD14", "CD68", "FCER1G","FCGR2A",       #D  MONO                  ITGAD
             "PECAM1","VWF","EMCN","AQP1","FLT1",           #E  Endo 
              "DERL3","MZB1","CD38","JCHAIN","IGKC",       #P  Plasma  
              "CD79A","CD19","BANK1","VPREB3","SPIB",         #B B cell     
             "TPSB2","TPSAB1","MS4A2","CTSG","HDC")        #M  Mast         TPSD1   MS4A2  HDC   CTSG
          
   genes %in% c_all$gene    
   

   
   
  markers <- c_all[c_all$gene %in% genes,]
  AverageHeatmap(object = tmp,   #great 
                 gene.order = genes,  
                 cluster.order = CO,
                 myanCol = color_umap2,
                 htCol = c("#0099CC", "white", "#CC0033"),
                 markerGene = markers$gene)
  show_col(color_umap2)
  
  
  # C1        C2        C3        C4        C5        C6        C7        C8        C9       C10       C11       C12       C13       C14       C15        T1 
  # "#EF6C00" "#F57C00" "#FB8C00" "#FF9800" "#FFA726" "#FFB74D" "#FFCC80" "#FFB300" "#FFC107" "#FFCA28" "#FFD54F" "#FFE082" "#FFEE58" "#FFF176" "#FFF59D" "#AFB42B" 
  # T2        T3        T4        T5        T6        T7        T8        F1        F2        F3        F4        F5        F6        D1        D2        D3 
  # "#C0CA33" "#CDDC39" "#D4E157" "#DCE775" "#E6EE9C" "#F0F4C3" "#F9FBE7" "#F44336" "#EF5350" "#E57373" "#EF9A9A" "#FFCDD2" "#FFEBEE" "#00BCD4" "#26C6DA" "#4DD0E1" 
  # D4        D5        D6        E1        E2        E3        E4        P1        P2        P3        P4        P5        B1        M1 
  # "#80DEEA" "#B2EBF2" "#E0F7FA" "#42A5F5" "#64B5F6" "#90CAF9" "#BBDEFB" "#1B5E20" "#2E7D32" "#388E3C" "#43A047" "#66BB6A" "#00897B" "#C62828" 
  
  annoGene <- c( "RGS5","CD36",            #C01  CD36+pericyte     
                  "PI16","WISP2",       #C02  PI16+ssCAF
                  "LRRC15","MMP11",         #C03  LRRC15+myCAF
                  "RERGL","BCAM",        #C04  RERGL+SMC  
                  "CD74","HLA-DRA",        #C05  CD37+apCAF
                  "SPP1","CD24","SERPINA1",  #C06 CD24+apCAF
                  "DPT","IGF1",                   #C07  ICAM1+ssCAF
                  "POSTN","PLAT",                 #C08  PLAT+myCAF
                  "CXCL8","HGF",                    #C09  HGF+iCAF  
                  "IL6","LIF","CSF3","SLC2A1",     #C10  LIF+iCAF
                   "MEIS1", "PARD3B",            #C11 MEIS1+CAF
                   "MYLK","CNN1","DES",         #C12 DES+SMC
                  "TOP2A","CENPF","MKI67")       #C13 TOP2A+CAF
  
  AverageHeatmap(object = tmp,
                 markerGene = markers$gene,
                 clusterAnnoName = F,
                 showRowNames = F,
                 gene.order = genes, 
                 markGenes = annoGene)
  
  # top_markers<- c_all %>%
  #   group_by(cluster) %>%
  #   slice_max(n = 3, order_by = avg_log2FC)
  jjDotPlot(object = tmp,
            gene = markers$gene,
            gene.order = genes, 
            xtree = F,
            ytree = F)
  
  jjDotPlot(object = tmp,
            gene = markers$gene,
            gene.order = genes,
            #id = 'seurat_cluster',
            split.by = 'Group',
            dot.col = c('#0099CC','#CC3333') ,
            xtree = F,
            ytree = F)
  
  jjDotPlot(object = tmp,
            gene = markers$gene,
            #id = 'celltype',
            split.by = 'Group',
            split.by.aesGroup = T,
            point.geom = F,
            tile.geom = T,
            xtree = F,
            ytree = F)
  
#valcano 
  dev.off()
  Idents(tmp)<-"Type"
  all.markers <- FindAllMarkers(tmp, only.pos = FALSE,
                                 min.pct = 0.25,
                                 logfc.threshold = 0)
  write.table(all.markers,"all_markers2.txt")

  #mygene <- c('LTB','SPP1','CD74','FAP')
  jjVolcano(diffData = all.markers,
            log2FC.cutoff = 0.5,
            col.type = "adjustP",
            #myMarkers = mygene,
            aesCol = c('purple','orange'),
            tile.col = color_Type)
  # flip the plot
  jjVolcano(diffData = all.markers,
          log2FC.cutoff = 0.5,
          tile.col = color_Type,
          fontface = 'italic',
          legend.position = c(0.8,0.2),
          flip = T)
  # make a polar plot
  jjVolcano(diffData = all.markers,
            log2FC.cutoff = 0.5,
            tile.col = color_Type,
            fontface = 'italic',
            back.col = "grey96",
            aesCol = c('#0099CC','#CC3373'),
            base_size = 16,
            legend.position = c(0.8,0.8),
            polar = T) #+ylim(-8,10)
            
  dev.off()
#tracksPlot
  tracksPlot(object = tmp,
             genes = markers$gene[1:5])
  
#dimplot_proportion  
  Idents(tmp)<-"seurat_clusters"  
  scatterCellPlot(object = tmp,
                  rm.axis = F,
                  cell.id = "seurat_clusters",
                  color = color_umap #ggsci::pal_npg()(9)
                  )
  DimPlot(object = tmp,reduction = "umap",
          cols = color_umap2,
          label = T)
  tmp$seurat_clusters
#featurePlot 
featurePlot(object = tmp,
              genes = c('LRRC15', 'PI16', 'SPP1', 'RGS5'),
              nrow = 2,ncol = 2,
              quantile.val = 1,
              add.corArrow = T,
              rm.legend = F,
              #add.strip = T,
              corLabel.dist = 0.1,
              color = c("grey90","#F5C6DC","red"),
              keep.oneCor = F) 
  

VlnPlot(tmp,
        features = c('LRRC15', 'PI16', 'SPP1', 'RGS5'),
        ncol = 2,
        slot = 'data',
        assay = 'RNA',
        log = T,
        pt.size = 0,
        cols =  Choose_col)
   
  dev.off()
  
}  
   
  
}



Cell_nk_Hierarchical_clustering_pan{
  library(tidyverse)
  library(ggplot2)
  library(ggtree)
  library(treeio)
  library(ggsci)
  library(cowplot)
  
  tmp <- readRDS("./data/d_seurat_A.rds")
  if(T){
    library(dplyr)
    meta<-as.data.frame(tmp@meta.data)
    df<-data.frame(barcode = rownames(meta),cluster = meta$seurat_clusters)
    #samp<-df %>% group_by(cluster) %>% sample_n(size=5000)  
    # samp<-df %>% group_by(cluster) %>% slice_sample(n=2000)  
    #samp<-df %>% group_by(cluster) %>% sample_frac(size=0.05)  
    samp<-df %>% group_by(cluster) %>% slice_sample(prop=0.1)  
    tmp<-tmp[,samp$barcode]
  }
  
  
  tmp$Type <- fct_relevel(tmp$Type ,rev(c("Epithelial cell" ,"T cell" ,"Fibroblasts cell", "Myeloid cell", "Endothelial cell","Plasma cell" ,"B cell","Mast cell")))
  table(tmp$Type)
  meta <- tmp@meta.data
  
  #Hierarchical clustering for cancer types
  his_cell_num = meta %>% filter((Group == "Tumor"))
  levels(his_cell_num$Disease) 
  test = his_cell_num %>% group_by(Disease) 
  
  DT = meta %>%  group_by(Type,Disease) 
  DT <- summarise(DT,count = n())
  DT<-data.frame(DT)
  DT<-na.omit(DT)
  
  test <- data.frame(summarise(test,
                               count = n()))
  
  DT <- merge(test,DT,by="Disease")
  DT['proportion'] <-DT$count.y /DT$count.x
  
  
  library(reshape2)
  DT2 <-DT[,c("Disease","Type","proportion")]
  DT2 <- dcast(DT2,Disease~Type)
  rownames(DT2) <-DT2$Disease
  DT2 <-DT2[,-1]
  DT2[is.na(DT2)] <- 0
  
  tree = hclust(vegan::vegdist(DT2, method = 'bray'), method = 'average')
  
  p1 = ggtree(tree) + geom_tiplab() +xlim(NA,2)
  p1
  
  
  color_Type=
    c("#C62828" ,"#00897B" ,"#388E3C" ,"#64B5F6" ,"#26C6DA" ,"#F44336" ,"#CDDC39" ,"#FFA726")
  
  g2 <-DT %>% ggplot(aes(y=proportion,x=Disease,fill=Type,colour = Type))+
    geom_bar(stat = "identity",width = 0.5,alpha=0.65) +
    theme_classic() +
    theme(axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line = element_blank())+
    scale_x_discrete(limits = rev(rev(c("BLCA","SCC","PMP","PDAC","CESC","BCC","BRCA","COLO","STAD","ESCC","OV","ccRCC","PRAD", "HCC","NSCLC"))),) +
    scale_fill_manual(values = rev(color_Type)) +
    scale_colour_manual(values = rev(color_Type)) +
    #  
    coord_flip()+
    ylab("Frequency") +
    xlab("") 
  
  g2
  
  g <- ggdraw()+
    draw_plot(p1, 0, 0.06, 0.6, 0.94)+
    draw_plot(g2, 0.2, 0, 0.8, 1)
  g
  
  ggsave("./result/Cell_nk_Hierarchical_clustering/Fb_clusters.pdf", g, height = 5, width =10)
}





scenic{
  
 seurat_Fb <- readRDS("./data/d_Fb7.rds") 
 seurat_Fb$celltype<-NULL
 seurat_Fb$RNA_snn_res.0.5<-NULL
 seurat_Fb$RNA_snn_res.0.02<-NULL
 seurat_Fb$RNA_snn_res.0.1<-NULL
 seurat_Fb$RNA_snn_res.0.15<-NULL
 seurat_Fb$RNA_snn_res.0.4<-NULL
 colnames(seurat_Fb@meta.data)
 seurat_Fb$Tissue
 seurat_Fb$Disease
 table(seurat_Fb$seurat_clusters)
 Idents(seurat_Fb)
 
 
 clusters<-levels(seurat_Fb$seurat_clusters)
 markers<-c("CD36","PI16","LRRC15","RERGL","CD37","CD24","ICAM1","PLAT","HGF","LIF","MEIS1","DES","TOP2A")
 barcode <- c()
 for(i in 1:13){
 cluster<-clusters[i]
 seurat_Fb_subset <-subset(seurat_Fb,idents = cluster)
 gene<-markers[i]
 num <-colnames(seurat_Fb_subset@assays$RNA@data)[seurat_Fb_subset@assays$RNA@data[gene,] > 1]
 if(length(num)>500){
   sm<-sample(num,500)
   barcode<-append(barcode,sm)
 }else{
   barcode<-append(barcode,num)
 }
 }
 seurat_Fb_subset <-subset(seurat_Fb,cells = barcode)
 table(seurat_Fb_subset$seurat_clusters)
 FeaturePlot(seurat_Fb_subset,"CD74")
 FeaturePlot(seurat_Fb_subset,"CXCL8")
 FeaturePlot(seurat_Fb_subset,"LRRC15")
 FeaturePlot(seurat_Fb_subset,"PI16") 
 FeaturePlot(seurat_Fb_subset,"ICAM1") 
 FeaturePlot(seurat_Fb_subset,"PLAT") 
 
 #devtools::install_github("aertslab/SCopeLoomR")
 library(SCopeLoomR)
 build_loom(file.name = "scenic_CAF.loom",dgem = seurat_Fb_subset@assays$RNA@data)
 meta<-seurat_Fb_subset@meta.data
 exprMat  <-  as.matrix(seurat_Fb_subset@assays$RNA@data)
 write.table(meta,'scenic_CAF_meta.xls',sep='\t',quote=F)
 write.csv(meta,'scenic_CAF_meta.csv',quote=F)
 write.table(exprMat,'scenic_CAF_exprMat.xls',sep='\t',quote=F)
 write.csv(exprMat,'scenic_CAF_exprMat.csv',quote=F)
 saveRDS(seurat_Fb_subset,"scenic_CAF.rds")
 seurat_Fb_subset <- readRDS("scenic_CAF.rds")
 seurat_Fb_subset<-NormalizeData(seurat_Fb_subset)
 dim(seurat_Fb_subset) #[1] 15559  6318
 scenic_CAF=CreateSeuratObject(seurat_Fb_subset@assays$RNA@data,assay = "RNA",min.cells = 100,meta.data = seurat_Fb_subset@meta.data)
 scenic_CAF<-FindVariableFeatures(scenic_CAF,nfeatures = 6000)
 
 markers %in% VariableFeatures(scenic_CAF)
 features<-VariableFeatures(scenic_CAF)
 dim(scenic_CAF)
 scenic_CAF <-  scenic_CAF[features,]
 dim(scenic_CAF)
 build_loom(file.name = "scenic_CAF2.loom",dgem = scenic_CAF@assays$RNA@data)
 meta<-scenic_CAF@meta.data
 write.table(meta,'scenic_CAF_meta2.txt',sep='\t',quote=F)
 write.csv(meta,'scenic_CAF_meta2.csv',quote=F)
 
 
#pyscenic_from_loom.sh  -> ctx.csv    aucell.loom  s  cenic_visualize.loom
 
 

 #visualize
 getwd()
 setwd("D:/pyscenic/")
 dir()
 #加载分析包
 library(SCopeLoomR)
 library(AUCell)
 library(SCENIC)
 #可视化相关包，多加载点没毛病
 library(dplyr)
 library(KernSmooth)
 library(RColorBrewer)
 library(plotly)
 library(BiocParallel)
 library(grid)
 library(ComplexHeatmap)
 library(data.table)
 library(ggplot2)
 library(pheatmap)
 library(Seurat)
 library(SCopeLoomR)
 library(AUCell)
 library(SCENIC)
 library(dplyr)
 library(KernSmooth)
 library(RColorBrewer)
 library(plotly)
 library(BiocParallel)
 library(loomR)
 dim(scenic_CAF)#[1] 6000 6318
 celltype <- scenic_CAF@meta.data$seurat_clusters
 assay <- scenic_CAF@assays$RNA
 loom <- open_loom("aucell.loom")

 regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
 regulons <- regulonsToGeneLists(regulons_incidMat)
 regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
 regulonAucThresholds <- get_regulon_thresholds(loom)
 close_loom(loom)
 meta <- scenic_CAF@meta.data
 colnames(meta)
 cellinfo <- meta[,c("seurat_clusters","nFeature_RNA","nCount_RNA", "Group")]
 colnames(cellinfo)=c('celltype', 'nGene' ,'nUMI','Group')
 cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
 #AUC=getAUC(sub_regulonAUC)
 selectedResolution <- "celltype"
 sub_regulonAUC <- regulonAUC
 rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
                cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                         selectedResolution])
 rss<-rss[,c("c01","c02","c03","c04","c05","c06",
             "c07","c08","c09","c10","c11","c12","c13")] 
 dev.off()
 rss=na.omit(rss)
 rssPlot <- plotRSS(rss)
 rssPlot_data<-rssPlot$df
 write.table(rssPlot_data,"rssPlot_data.txt")
 saveRDS(rss,"rss.rds")
 save(regulonAUC,rssPlot,regulons,file='regulon_RSS.Rdata')
 source('function_pyscenic_visualize.R')
 
 dim(scenic_CAF)#[1] 6000 6318
 scenic_CAF@reductions$umap$cell.embeddings <- embed_umap 
 saveRDS(seurat_Fb_subset,"scenic_CAF3.rds")
 saveRDS(seurat_Fb_subset,'subset.rds')
 sce<-seurat_Fb_subset
 plot_pyscenic(inloom='aucell.loom',incolor=incolor,inrss="rss.rds",inrds='subset.rds',infun='median', 
               ct.col="seurat_clusters",inregulons=NULL,ingrn='grn.txt',ntop1=5,ntop2=50)

plot_pyscenic{
  library(Seurat)
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(KernSmooth)
  library(RColorBrewer)
  library(plotly)
  library(BiocParallel)
  library(pheatmap)
  
  library(cowplot)
  library(ggpubr)
  library(ggsci)
  library(ggplot2)
  library(tidygraph)
  library(ggraph)
  library(stringr)
  
  
  
  colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                        pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                        pal_locuszoom("default")(7),pal_igv("default")(51),
                        pal_uchicago("default")(9),pal_startrek("uniform")(7),
                        pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                        pal_simpsons("springfield")(16),pal_gsea("default")(12)))
  len <- 100
  incolor<-c(brewer.pal(8, "Dark2"),brewer.pal(12, "Paired"),brewer.pal(8, "Set2"),brewer.pal(9, "Set1"),colpalettes,rainbow(len))
  
  plot_pyscenic <- function(inloom='aucell.loom',incolor=incolor,inrss='seurat_annotations_rss.rds',inrds='subset.rds',infun='median', ct.col='seurat_annotations',inregulons=NULL,ingrn='grn.tsv',ntop1=5,ntop2=50){
    ###load data
    loom <- open_loom(inloom)
    
    regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
    regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
    regulonAucThresholds <- get_regulon_thresholds(loom)
    #embeddings <- get_embeddings(loom)
    close_loom(loom)
    
    rss <- readRDS(inrss)
    sce <- readRDS(inrds)
    embeddings <- sce@reductions$umap
    ##calculate  RSS fc
    df = do.call(rbind,
                 lapply(1:ncol(rss), function(i){
                   dat= data.frame(
                     regulon  = rownames(rss),
                     cluster =  colnames(rss)[i],
                     sd.1 = rss[,i],
                     sd.2 = apply(rss[,-i], 1, get(infun))
                   )
                 }))
    
    df$fc = df$sd.1 - df$sd.2
    
    #select top regulon
    ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
    
    ntopgene <- unique(ntopg$regulon)
    write.table(ntopgene,'sd_regulon_RSS.list',sep='\t',quote=F,row.names=F,col.names=F)
    #plot rss by cluster
    
    #using plotRSS
    rssPlot <- plotRSS(rss)
    regulonsToPlot <- rssPlot$rowOrder
    rp_df <- rssPlot$df
    
    write.table(regulonsToPlot,'rss_regulon.list',sep='\t',quote=F,row.names=F,col.names=F)
    write.table(rp_df,'rssPlot_data.xls',sep='\t',quote=F)
    nlen <- length(regulonsToPlot)
    hei <- ceiling(nlen)*0.4
    blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    lgroup <- levels(rssPlot$df$cellType)
    
    nlen2 <- length(lgroup)
    wei <- nlen2*2
    pdf(paste0('regulons_RSS_',ct.col,'_in_dotplot.pdf'),wei,hei)
    print(rssPlot$plot)
    dev.off()
    
    # sd top gene
    anrow = data.frame( group = ntopg$cluster)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(anrow$group)
    annotation_colors <- list(group=lcolor)
    
    pn1 = rss[ntopg$regulon,]
    pn2 = rss[unique(ntopg$regulon),]
    rownames(pn1) <-  make.unique(rownames(pn1))
    rownames(anrow) <- rownames(pn1)
    scale='row'
    hei <- ceiling(length(ntopg$regulon)*0.4)
    pdf(paste0('regulon_RSS_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
    )
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons')
    )
    dev.off()
    
    #plotRSS gene
    
    pn2 = rss[unique(rp_df$Topic),]
    scale='row'
    hei <- ceiling(length(unique(rp_df$Topic))*0.4)
    pdf(paste0('regulon_RSS_in_plotRSS_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='plotRSS unique regulons')
    )
    dev.off()
    
    #all regulons
    
    hei <- ceiling(length(rownames(rss))*0.2)
    pdf(paste0('all_regulons_RSS_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(rss,scale=scale,show_rownames = T,main='all regulons RSS')
    )
    dev.off()
    #plot rss by all cells
    if (is.null(inregulons)){
      inregulons <- regulonsToPlot
    }else{
      inregulons <- intersect(inregulons,rownames(rss))
      regulonsToPlot <- inregulons
      
    }
    pn3=as.matrix(regulonAUC@assays@data$AUC)
    regulon <- rownames(pn3)
    #regulon <- inregulons
    pn3 <- pn3[regulon,]
    #pn3 <- pn3[,sample(1:dim(pn3)[2],500)]
    
    sce$group1=sce@meta.data[,ct.col]
    
    meta <- sce@meta.data
    meta <- meta[order(meta$group1),]
    #meta <- meta[colnames(pn3),]
    ancol = data.frame(meta[,c('group1')])
    colnames(ancol) <- c('group1')
    rownames(ancol) <- rownames(meta)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(ntopg$cluster)
    annotation_colors <- list(group1 =lcolor)
    
    df1 <- ancol
    df1$cell <- rownames(df1)
    df1 <- df1[order(df1$group1),]
    pn3 <- pn3[,rownames(df1)]
    torange=c(-2,2)
    pn3 <- scales::rescale(pn3,to=torange)
    pn3 <- pn3[,rownames(ancol)]
    
    scale='none'
    hei <- ceiling(length(unique(regulon))*0.2)
    pdf(paste0('all_regulon_activity_in_allcells.pdf'),10,hei)
    print(
      pheatmap(pn3,annotation_col = ancol,scale=scale,annotation_colors=annotation_colors,show_rownames = T,show_colnames = F,cluster_cols=F)
    )
    #pheatmap(pn3,scale=scale,show_rownames = T, show_colnames = F,cluster_cols=F)
    dev.off()
    
    #plot in seurat
    regulonsToPlot = inregulons
    sce$sub_celltype <- sce@meta.data[,ct.col]
    sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
    
    cellClusters <- data.frame(row.names = colnames(sce),
                               seurat_clusters = as.character(sce$seurat_clusters))
    cellTypes <- data.frame(row.names = colnames(sce),
                            celltype = sce$sub_celltype)
    
    sce@meta.data = cbind(sce@meta.data ,t(sub_regulonAUC@assays@data@listData$AUC[regulonsToPlot,]))
    Idents(sce) <- sce$sub_celltype
    
    nlen <- length(regulonsToPlot)
    hei <- ceiling(nlen)*0.4
    blu<-colorRampPalette(brewer.pal(9,"Blues"))(100)
    nlen2 <- length(unique(sce$sub_celltype))
    wei <- nlen2*2
    pdf('regulons_activity_in_dotplot.pdf',wei,hei)
    print(DotPlot(sce, features = unique(regulonsToPlot)) + coord_flip()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
            scale_color_gradientn(colours = blu)
    )
    dev.off()
    
    
    hei=ceiling(nlen/4)*4
    pdf('regulons_activity_in_umap.pdf',16,hei)
    print(RidgePlot(sce, features = regulonsToPlot , ncol = 4))
    print(VlnPlot(sce, features = regulonsToPlot,pt.size = 0 ))
    print(FeaturePlot(sce, reduction="umap",features =regulonsToPlot))
    dev.off()
    
    grn <- read.table(ingrn,sep='\t',header=T,stringsAsFactors=F)
    inregulons1=gsub('[(+)]','',inregulons)
    
    c1 <- which(grn$TF %in% inregulons1)
    grn <- grn[c1,]
    #edge1 <- data.frame()
    #node1 <- data.frame()
    pdf(paste0(ntop2,'_regulon_netplot.pdf'),10,10)
    for (tf in unique(grn$TF)) {
      tmp <- subset(grn,TF==tf)
      if (dim(tmp)[1] > ntop2) {
        tmp <- tmp[order(tmp$importance,decreasing=T),]
        tmp <- tmp[1:ntop2,]
      }
      node2 <- data.frame(tmp$target)
      node2$node.size=1.5
      node2$node.colour <- 'black'
      colnames(node2) <- c('node','node.size','node.colour')
      df1 <- data.frame(node=tf,node.size=2,node.colour='#FFDA00')
      node2 <- rbind(df1,node2)
      
      
      edge2 <- tmp
      colnames(edge2) <- c('from','to','edge.width')
      edge2$edge.colour <- "#1B9E77"
      torange=c(0.1,1)
      edge2$edge.width <- scales::rescale(edge2$edge.width,to=torange)
      
      graph_data <- tidygraph::tbl_graph(nodes = node2, edges = edge2, directed = T)
      p1 <- ggraph(graph = graph_data, layout = "stress", circular = TRUE) + geom_edge_arc(aes(edge_colour = edge.colour, edge_width = edge.width)) +
        scale_edge_width_continuous(range = c(1,0.2)) +geom_node_point(aes(colour = node.colour, size = node.size))+ theme_void() +
        geom_node_label(aes(label = node,colour = node.colour),size = 3.5, repel = TRUE)
      p1 <- p1 + scale_color_manual(values=c('#FFDA00','black'))+scale_edge_color_manual(values=c("#1B9E77"))
      print(p1)
    }
    dev.off()
    #plot activity heatmap
    meta <- sce@meta.data
    celltype <- ct.col
    cellsPerGroup <- split(rownames(meta),meta[,celltype])
    sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
    # Calculate average expression:
    regulonActivity_byGroup <- sapply(cellsPerGroup,
                                      function(cells)
                                        rowMeans(getAUC(sub_regulonAUC)[,cells]))
    scale='row'
    rss <- regulonActivity_byGroup
    hei <- ceiling(length(regulonsToPlot)*0.4)
    pn1 <- rss[regulonsToPlot,]
    pdf(paste0('regulon_activity_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,scale=scale,show_rownames = T, main='regulons activity')
    )
    dev.off()
    
    hei <- ceiling(length(rownames(rss))*0.2)
    pdf(paste0('all_regulons_activity_in_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(rss,scale=scale,show_rownames = T,main='all regulons activity')
    )
    dev.off()
    
    ##calculate fc
    df = do.call(rbind,
                 lapply(1:ncol(rss), function(i){
                   dat= data.frame(
                     regulon  = rownames(rss),
                     cluster =  colnames(rss)[i],
                     sd.1 = rss[,i],
                     sd.2 = apply(rss[,-i], 1, get(infun))
                   )
                 }))
    
    df$fc = df$sd.1 - df$sd.2
    
    #select top regulon
    ntopg <- df %>% group_by(cluster) %>% top_n(ntop1, fc)
    
    ntopgene <- unique(ntopg$regulon)
    write.table(ntopgene,'sd_regulon_activity.list',sep='\t',quote=F,row.names=F,col.names=F)
    
    anrow = data.frame( group = ntopg$cluster)
    lcolor <- incolor[1:length(unique(ntopg$cluster))]
    names(lcolor) <- unique(anrow$group)
    annotation_colors <- list(group=lcolor)
    pn1 = rss[ntopg$regulon,]
    pn2 = rss[unique(ntopg$regulon),]
    rownames(pn1) <-  make.unique(rownames(pn1))
    rownames(anrow) <- rownames(pn1)
    scale='row'
    hei <- ceiling(length(ntopg$regulon)*0.4)
    pdf(paste0('regulon_activity_in_sd_topgene_',ct.col,'.pdf'),wei,hei)
    print(
      pheatmap(pn1,annotation_row = anrow,scale=scale,annotation_colors=annotation_colors,show_rownames = T,main='sd top regulons')
    )
    print(
      pheatmap(pn2,scale=scale,show_rownames = T, main='sd top unique regulons ')
    )
    dev.off()
    
  }
}
 
}






  


  

















discrete_colors{
    library(ggplot2)
    library(scales)
    getwd()
    dir()
    colors_dt<-read.csv("./color/color.csv")
    colors_list<-as.list(colors_dt)
    saveRDS(colors_list,"colors_list.rds")
    colors_list<-readRDS("./color/colors_list.rds")
    
    choose_col<-function(a,b,c){
      col<-c()
      s<-sample(c(a,b,c),length(colors_list),replace = T) 
      for(i in 1:length(colors_list)){
        S<-s[i]
        col<-append(col,colors_list[[i]][S])
      }
      return(col)
    }
    Choose_col<-choose_col(3,4,5)
    show_col(Choose_col)
    show_col(colors_list[[3]])
    
    col_map<-as.list(as.data.frame(t(colors_dt)))
    saveRDS(col_map,"col_map.rds")
    col_map<-readRDS("./colors/col_map.rds")
    scales::show_col(col_map$V8)
   
}



themes{
library(ggthemes) 

p1+theme_void()           #***
p1+theme_minimal()        #***
p1+theme_classic()        #***
p1+theme_bw()             #***
  
p1+ggthemes::theme_base()   #****
p1+ggthemes::theme_par()    #****
p1+ggthemes::theme_tufte()  #****
p1+ggthemes::theme_few()    #****
p1+ggthemes::theme_clean()  #****
p1+ggthemes::theme_map()    #****
p1+ggthemes::theme_pander() #****
p1+ggthemes::theme_solid()  #****
p1+ggthemes::theme_wsj()    #****
p1+ggthemes::theme_calc()   #****
p1+theme_cxf  #*****
p1+theme_cxf2 #*****  

theme_cxf<-  theme(panel.background = element_rect(fill = "transparent", color = "black",size = 1.5), 
                   legend.key = element_rect(fill = "transparent", color = "transparent"),
                   text = element_text(color = "black", size = 14,vjust=0.5),
                   plot.title = element_text(hjust = .5), 
                   axis.text = element_text(color = "black",size = 14))
theme_cxf2 <- theme_classic() +
  theme(panel.background=element_rect(fill='transparent', color='black',size = 1.5), 
        text = element_text(color = "black", size = 16,vjust=0.5),
        panel.border=element_rect(fill='transparent', color='transparent'), 
        panel.grid=element_blank(), axis.title = element_text(color='black', vjust=0.5),
        axis.text = element_text(color="black"),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'), 
        axis.ticks.margin = unit(0.2,"lines"), legend.title=element_blank(), 
        plot.title = element_text(family = 'serif',face = 'bold',colour = 'black',size = 20,hjust = .5), 
        plot.subtitle = element_text(family = 'serif',face = 'bold',colour = 'black',size = 16,hjust = .5),
        legend.key=element_rect(fill='transparent', color='transparent')) 

}

layout{
library("gridExtra")
grid.arrange(p1_jama, p2_jama, ncol = 2)   

library(patchwork)  
#+/|

}







getwd()
setwd("D:/data_subset5/subset/seurat_all")
dir()
seurat_A <- readRDS("./data/d_seurat_A.rds")
FeaturePlot(seurat_A,features = c("CCL19","CCL21","CXCL13","PDPN"))
FeaturePlot(seurat_A,features = c("TNFRSF10A"))


library(Seurat)
FeaturePlot(seurat_Fb,features = c("PDGFRB"))
FeaturePlot(seurat_Fb,features = c("ABCC9","CD36"))
FeaturePlot(seurat_Fb,features = c("KCNJ8","CD36"))
FeaturePlot(seurat_Fb,features = c("H1GD1B","CD36"))
FeaturePlot(seurat_Fb,features = c("ATP1A2","CD36"))
FeaturePlot(seurat_Fb,features = c("SLC6A1","CD36"))
FeaturePlot(seurat_Fb,features = c("IFIT1","CD36"))
FeaturePlot(seurat_Fb,features = c("IFIT2","CD36"))
FeaturePlot(seurat_Fb,features = c("ATP1A2","CD36"))
FeaturePlot(seurat_Fb,features = c("MX1","CD36"))
FeaturePlot(seurat_Fb,features = c("MYH11","CD36"))
FeaturePlot(seurat_Fb,features = c("TAGLN","CD36"))
FeaturePlot(seurat_Fb,features = c("ISLR","CD36"))
FeaturePlot(seurat_Fb,features = c("SLC6A13","CD36"))
FeaturePlot(seurat_Fb,features = c("FGB","CD36"))
FeaturePlot(seurat_Fb,features = c("CD62E","CD36"))

FeaturePlot(seurat_Fb,features = c("ICAM2","CD36"))
FeaturePlot(seurat_Fb,features = c("VCAM1","CD36"))
FeaturePlot(seurat_Fb,features = c("CD62E","CD36"))
dev.off()


