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