library(spacexr)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(doParallel)
library(tibble)
library(data.table)
library(ggplot2)
library(ggrastr)
library(ggpubr)


Figure5A-D{
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