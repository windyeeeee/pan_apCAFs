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
