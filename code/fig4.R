setwd("~/Project/ST/data/SC/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

sce.mer <- readRDS("sce.mer.RDS")
ST_merge <- readRDS("ST_merge.RDS")

cds <- c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")
names(cds) <- c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")
cds <- lapply(cds, function(x){
  paste0(x,"/monocle3.rds") %>% readRDS()
})
pseudo <- lapply(cds, function(x){
  x@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] %>% as.data.frame()
}) %>% bind_rows()

pseudo <- cbind(sce.mer@meta.data[,"cellIdent_sub",drop = F],pseudo)
names(pseudo)[2] <- "pseudo"
pdf("fig4b.pdf",height = 5,width = 8)
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

cds <- readRDS("monocle3.rds")
Thyrocyte <- readRDS("Thyrocyte.RDS")
pseudo <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] %>% as.data.frame()
score <- Thyrocyte@meta.data[,c("cellIdent_sub","TDS1","BRAF1","RAS1","nCount_RNA")]
score <- cbind(score,pseudo)
names(score)[6] <- "pseudoTime"
score <- score[order(score$nCount_RNA),]

pdf("fig4c.pdf",height = 5,width = 6)
ggplot(score, aes(x = pseudoTime, y = TDS1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("TDS")
ggplot(score, aes(x = pseudoTime, y = BRAF1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("BRAF")
ggplot(score, aes(x = pseudoTime, y = RAS1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("RAS")
dev.off()

cds <- readRDS("monocle3.rds")
Thyrocyte <- readRDS("Thyrocyte.RDS")
pseudo <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] %>% as.data.frame()
score <- Thyrocyte@meta.data[,c("cellIdent_sub","TDS1","BRAF1","RAS1","nCount_RNA")]
score <- cbind(score,pseudo)
names(score)[6] <- "pseudoTime"
score <- score[order(score$nCount_RNA),]

pdf("fig4de.pdf",height = 5,width = 6)
ggplot(score, aes(x = pseudoTime, y = TDS1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("TDS")
ggplot(score, aes(x = pseudoTime, y = BRAF1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("BRAF")
ggplot(score, aes(x = pseudoTime, y = RAS1)) + 
  geom_point(aes(color = log10(nCount_RNA)))+
  scale_color_gradient(low = "white",high = "blue")+
  geom_smooth(method = "lm",color = "black") +
  ggpubr::stat_cor() + theme_classic() + ylab("RAS")
dev.off()


pseudo <- pseudotime(cds) %>% as.data.frame()
names(pseudo) <- "pseudoTime"
pseudo <- cbind(pseudo,Thyrocyte@meta.data[,"cellIdent_sub",drop = F])
scores <- cbind(pseudo,Thyrocyte@meta.data[,c("TDS1","BRAF1","RAS1")])

scores_ct <- ddply(scores,"cellIdent_sub",function(x){
  data.frame(pseudoTime=median(x[["pseudoTime"]]),
             TDS=median(x[["TDS1"]]),
             BRAF=median(x[["BRAF1"]]),
             RAS=median(x[["RAS1"]]))
})
rownames(scores_ct) <- scores_ct$cellIdent_sub
scores_ct$cellIdent_sub <- NULL
out.dist <- dist(scores_ct,method="euclidean")
out.hclust <- hclust(out.dist, method = "ward.D")
out.id <- cutree(out.hclust,k=3) %>% as.data.frame()

library(factoextra)
pdf("fig4fg.pdf",height=5,width=5)
fviz_nbclust(x = tmp,FUN = hcut,method="silhouette")
fviz_nbclust(x = tmp,FUN = hcut,method="wss")
fviz_cluster(list(data = tmp, cluster = sub))+theme_classic()
fviz_dend(out.hclust,k=3,cex = 0.5, k_colors = c("#00AFBB", "#E7B800", "#FC4E07"),color_labels_by_k = TRUE,rect = TRUE)
dev.off()

markers <- readRDS("~/Project/ST/data/SC/TDS/state_markers.RDS")
features <- markers %>% group_by(cluster) %>% top_n(20,avg_log2FC) %>% .$gene %>% unique()
Thyrocyte <- NormalizeData(Thyrocyte)
exp <- AverageExpression(Thyrocyte,assays = "RNA",features = features,group.by = "state")[["RNA"]]
exp <- t(scale(t(exp)))
pdf("fig4h.pdf",height = 4,width = 10)
DotPlot(Thyrocyte,features = features) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red")
dev.off()

norm_weights <- readRDS("norm_weights_trajectory_sub.RDS")
norm_weights$barcode <- rownames(norm_weights)
weight_sub <- reshape2::melt(norm_weights,id = "barcode")
names(weight_sub)[2:3] <- c("cluster","ratio")
weight_sub <- plyr::ddply(weight_sub,"barcode",function(x){
  x[x$ratio == max(x$ratio),] %>% 
    .[!(duplicated(.$barcode)),]
})
weight_sub <- weight_sub[grep("cluster",weight_sub$cluster),]
sampleInfo <- cbind(barcode = rownames(ST_merge@meta.data),ST_merge@meta.data[,"samples",drop = F])
weight_sub <- left_join(weight_sub,sampleInfo,by = "barcode")
weight_sub$samples <- as.character(weight_sub$samples)
names(weight_sub)[4] <- "sample"
rownames(weight_sub) <- weight_sub$barcode

coords <- lapply(ST_merge@images, function(x)x@coordinates[rownames(x@coordinates) %in% rownames(norm_weights),])

neighbor <- apply(weight_sub, 1, function(x){
  samples <- x[["sample"]]
  rows <- coords[[samples]][x[["barcode"]],"row"]
  cols <- coords[[samples]][x[["barcode"]],"col"]
  surround <- coords[[samples]] %>% 
    filter(((abs(row - rows) + abs(col - cols)) == 2)&(abs(row - rows) != 2))
  return(surround)
})

for (i in names(neighbor)) {
  if (nrow(neighbor[[i]]) < 6) {
    neighbor[[i]] <- NULL
  }
}

for (i in names(neighbor)) {
  neighbor[[i]]$core <- weight_sub[i,"barcode"]
  neighbor[[i]]$neighbor <- rownames(neighbor[[i]])
}

neighbor <- bind_rows(neighbor)

neighbor <- left_join(neighbor,norm_weights,by = c("neighbor" = "barcode"))
core <- plyr::ddply(neighbor,"core",function(x){
  colMeans(x[,8:ncol(x)])
})

cancerInfo <- cbind(core = rownames(ST_merge@meta.data),ST_merge@meta.data[,"CancerType",drop = F])
core <- merge(cancerInfo,core,by = "core")
core <- merge(weight_sub[,c("barcode","cluster")],core,by.x = "barcode",by.y ="core")
names(core)[1] <- "core"
core_mean <- core %>% group_by(CancerType,cluster) %>% summarise_if(is.numeric,~ mean(.))
core_mean <- plyr::ddply(core,"CancerType",function(x){
  plyr::ddply(x,"cluster",function(xx){
    colMeans(xx[,4:ncol(xx)])
  })
})

core_mean <- core_mean[core_mean$CancerType != "N",]
core_mean$cluster <- factor(core_mean$cluster,levels = paste0("cluster",1:3))
core_mean <- core_mean[order(core_mean$CancerType,core_mean$cluster),]
anno <- core_mean[,1:2]
htdata <- as.matrix(core_mean[,3:(ncol(core_mean)-3)])
rownames(htdata) <- anno$cluster

rowanno <- ComplexHeatmap::rowAnnotation(
  df = anno,
  col = list(CancerType = c("ATC" = "#F8766D",
                            "LPTC" = "#C77CFF",
                            "PTC" = "#00BFC4"),
             cluster = c("cluster1" = "#00AFBB",
                         "cluster2" = "#FC4E07",
                         "cluster3" = "#E7B800"))
)

colanno <- ComplexHeatmap::HeatmapAnnotation(
  cellType = colnames(htdata),
  col = list(cellType = c("B" = "#BC3C29FF",
                          "Tcells" = "#0072B5FF",
                          "Fibroblast" = "#7876B1FF",
                          "Endothelial" = "#6F99ADFF",
                          "NK" = "#FFDC91FF",
                          "Dendritic" = "#E18727FF",
                          "Granulocyte" = "#EE4C97FF",
                          "Macrophage" = "#20854EFF"))
)
col_fun <- circlize::colorRamp2(breaks = c(0,max(htdata)),colors = c("white","#01665e"))

splits <- rep(1:3,each = 3)

htdata_immu <- t(scale(t(htdata[,-(grep("^Endothelial|^Fibroblast",colnames(htdata)))])))
htdata_pa <- t(scale(t(htdata[,grep("^Endothelial|^Fibroblast",colnames(htdata))])))

htdata_immu <- htdata_immu[,order(colnames(htdata_immu))] %>% .[,c(1:14,19:26,15:18,27:35,39:46,36:38)]
htdata_pa <- htdata_pa[,order(colnames(htdata_pa))] %>% .[,c(1:10,16:23,11:15)]

col_fun <- circlize::colorRamp2(breaks = c(-4,0,5),colors = c("blue","white","red"))

pdf("fig4ik.pdf",height = 8,width = 12)
ComplexHeatmap::Heatmap(
  htdata_immu,
  name = "z.score",
  left_annotation = rowanno,
  show_row_names = F,
  column_title = NULL,
  column_names_rot = 45,
  cluster_rows = F,
  cluster_columns = F,
  row_split = splits,
  row_title = NULL,
  col = col_fun
)
ComplexHeatmap::Heatmap(
  htdata_pa,
  name = "z.score",
  left_annotation = rowanno,
  show_row_names = F,
  column_title = NULL,
  column_names_rot = 45,
  cluster_rows = F,
  cluster_columns = F,
  row_split = splits,
  row_title = NULL,
  col = col_fun
)
dev.off()
