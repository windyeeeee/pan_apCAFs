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
install.packages("devtools")
devtools::install_local("PATH/TO/DIRECTORY/CytoTRACE_0.3.3.tar.gz")

pip install scanoramaCT
pip install numpy

library(CytoTRACE) 
results <- CytoTRACE(spp1_ko_caf_expr, ncores = 8, subsamplesize = 1000)
plotCytoTRACE(results, phenotype = marrow_10x_pheno)
}