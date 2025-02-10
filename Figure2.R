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
