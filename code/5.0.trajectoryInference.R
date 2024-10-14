library(tidyverse)
library(Seurat)

cds <- c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")
names(cds) <- c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")
cds <- lapply(cds, function(x){
  paste0("~/ST/data/SC/trajectory/",x,"/monocle3.rds") %>% readRDS()
})
pseudo <- lapply(cds, function(x){
  x@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] %>% as.data.frame()
}) %>% bind_rows()

sce.mer <- readRDS("~/ST/data/SC/sce.mer.RDS")

pseudo <- cbind(sce.mer@meta.data[,"cellIdent_sub",drop = F],pseudo)
names(pseudo)[2] <- "pseudo"
pdf("~/ST/result/SC/plot/trajectory/fig4b.pdf",height = 5,width = 8)
for (i in c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")) {
  pseudo_sub <- pseudo[grep(i,pseudo$cellIdent_sub),]
  pseudo_sub <- pseudo_sub[order(pseudo_sub$cellIdent_sub),]
  fills <- pal_d3("category20")(length(unique(pseudo_sub$cellIdent_sub)))
  names(fills) <- unique(pseudo_sub$cellIdent_sub)
  lv <- pseudo_sub %>% group_by(cellIdent_sub) %>% summarise(median=median(pseudo)) %>% arrange(median)
  lv$cellIdent_sub <- as.character(lv$cellIdent_sub)
  pseudo_sub$cellIdent_sub <- factor(pseudo_sub$cellIdent_sub,levels = lv$cellIdent_sub)
  p <- ggplot(pseudo_sub,aes(cellIdent_sub,pseudo,fill=cellIdent_sub)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = fills[lv$cellIdent_sub]) +
    ggtitle(i)
  print(p)
}
dev.off()

