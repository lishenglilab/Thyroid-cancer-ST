setwd("~/Project/ST/data/ST/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

filename <- dir(getwd())
ST_samples <- paste(filename,"/filtered_feature_bc_matrix",sep = "")
names(ST_samples) <- filename
ST_samples <- lapply(ST_samples,Read10X)
for (i in names(ST_samples)) {
  ST_samples[[i]] <- CreateSeuratObject(ST_samples[[i]],project = i,
                                        assay = "Spatial")
  ST_samples[[i]]$slice <- i
  ST_samples[[i]]$region <- "thyroid"
  img <- Seurat::Read10X_Image(image.dir = paste(i,"/spatial",sep = ""))
  Seurat::DefaultAssay(object = img) <- 'Spatial'
  img <- img[colnames(x = ST_samples[[i]])]
  ST_samples[[i]][['image']] <- img
  ST_samples[[i]] <- SCTransform(ST_samples[[i]], assay = "Spatial")
  ST_samples[[i]] <- RunPCA(ST_samples[[i]], assay = "SCT") %>% 
    FindNeighbors(., reduction = "pca", dims = 1:30) %>%
    FindClusters() %>%
    RunUMAP(., reduction = "pca", dims = 1:30)
}

ST_split <- ST_samples

ST_split <- lapply(ST_split, function(x){
  FindNeighbors(x, reduction = "pca", dims = 1:30) %>%
    FindClusters(.,resolution = seq(0.2,0.8,.05))
})

saveRDS(ST_split,"ST_split.RDS")

ST_merge <- Reduce(merge,ST_split)

saveRDS(ST_merge,"ST_merge.RDS")

pdf("fig1b.pdf",height = 4,width = 5.5)
DimPlot(ST_merge,group.by = "orig.ident")+ggtitle("samples")
dev.off()

pdf("fig1c.pdf",height = 4,width = 15)
VlnPlot(ST_merge,features = "nCount_Spatial",pt.size = 0,group.by = "orig.ident")+
  NoLegend()+ggtitle("nUMI")
VlnPlot(ST_merge,features = "nFeature_Spatial",pt.size = 0,group.by = "orig.ident")+
  NoLegend()+ggtitle("nGene")
dev.off()

pdf("fig1d.pdf",height = 4,width = 5)
FeaturePlot(ST_merge,features = "nCount_Spatial",order = T)+ggtitle("nGene")
FeaturePlot(ST_merge,features = "nFeature_Spatial",order = T)+ggtitle("nUMI")
dev.off()

markers <- readRDS("markerDiffCancerType.RDS")
features <- markers %>% group_by(cluster) %>% top_n(10,avg_log2FC) %>% .$gene
Idents(ST_merge) <- "samples"
avg <- AverageExpression(ST_merge,assays = "SCT",features = features)[["SCT"]]
avg <- t(scale(t(avg)))
colnames(avg) <- substr(colnames(avg),1,regexpr("[-]",colnames(avg)) - 1)

colanno <- HeatmapAnnotation(cancerType = colnames(avg),
                             col = list(cancerType = c("N" = "#F8766D",
                                                       "PTC" = "#7CAE00",
                                                       "LPTC" = "#00BFC4",
                                                       "ATC" = "#C77CFF")))
col_fun <- circlize::colorRamp2(breaks = c(min(avg),0,max(avg)),colors = c("blue","white","red"))

splits <- rep(1:3,each = 3)

pdf("fig1e.pdf",height = 12,width = 8)
Heatmap(avg,name = "Z score",
        show_column_names = F,col = col_fun,top_annotation = colanno,
        cluster_rows = F,cluster_columns = F,rect_gp = gpar(col = "white", lwd = 1))
dev.off()

gsea <- readRDS("gseaCancerType.RDS")

gsea_table <- lapply(gsea, function(x)x@result)
for (i in names(gsea_table)) {
  gsea_table[[i]]$group <- i
}

gsea_table <- Reduce(rbind,gsea_table)

gsea_table$ID <- gsea_table$ID %>% gsub("HALLMARK_","",.) %>% gsub("_"," ",.) %>% 
  tolower() %>% Hmisc::capitalize()

pw <- as.data.frame(table(gsea_table$ID))
pw_ord <- ddply(gsea_table,"ID",function(x)x[1,"group",drop = F])
pw_ord$group <- factor(pw_ord$group,levels = c("N","PTC","LPTC","ATC"))
pw <- left_join(pw,pw_ord,by = c("Var1"="ID")) %>% .[order(.$Freq,.$group),]

gsea_table$ID <- factor(gsea_table$ID,levels = as.character(pw$Var1))

gsea_table$group <- factor(gsea_table$group,levels = c("N","PTC","LPTC","ATC"))

pdf("fig1f.pdf",height = 6,width = 9)
ggplot(gsea_table, aes(x = ID, y = group)) + 
  geom_point(aes(size = -(log10(p.adjust)), color = NES)) +
  geom_point(aes(size = -(log10(p.adjust))),shape = 1) +
  theme_bw(base_size = 14) +
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
  ylab(NULL)+ggtitle("hallmark") + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        axis.ticks = element_blank(),panel.grid = element_blank())
dev.off()
