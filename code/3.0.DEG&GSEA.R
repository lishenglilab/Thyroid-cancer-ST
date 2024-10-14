library(tidyverse)
library(Seurat)
library(clusterProfiler)

ST_merge <- readRDS("~/ST/data/ST/ST_merge.RDS")
Idents(ST_merge) <- "cancerType"

markers <- FindAllMarkers(ST_merge)
markers <- readRDS("~/ST/data/ST/deg/markerDiffCancerType.RDS")
hallmark <- read.gmt("/home/reference/GenePathway/HALLMARKS/h.all.v2022.1.Hs.symbols.gmt")
gsea <- list()
for (i in unique(markers$cluster)) {
  genetable <- markers %>% filter(cluster %in% i)
  genelist <- genetable$avg_log2FC
  names(genelist) <- genetable$gene
  genelist <- sort(genelist,decreasing = T)
  gsea[[i]] <- GSEA(genelist,TERM2GENE = hallmark)
}
saveRDS(gsea,"~/ST/data/ST/deg/gseaCancerType.RDS")
