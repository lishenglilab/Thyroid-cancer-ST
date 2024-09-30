setwd("~/Project/ST/data/SC/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

sce.mer <- readRDS("sce.mer.RDS")
ST_merge <- readRDS("ST_merge.RDS")

sce.mer.split <- SplitObject(sce.mer,"cellIdent")
pdf("fig3a.pdf",width = 6,height = 4)
for (i in names(sce.mer.split)) {
  myLabel <- paste0(names(table(sce.mer.split[[i]]$cellIdent_sub))," (",table(sce.mer.split[[i]]$cellIdent_sub),")")
  p <- DimPlot(sce.mer.split[[i]],group.by = "cellIdent_sub",raster=FALSE) + 
    scale_color_d3(palette = "category20",label = myLabel)
  print(p)
}
dev.off()

markers <- list(
  B = readRDS("/B_marker.RDS"),
  Dendritic = readRDS("Dendritic_marker.RDS"),
  Endothelial = readRDS("Endothelial_marker.RDS"),
  Fibroblast = readRDS("Fibroblast_marker.RDS"),
  Granulocyte = readRDS("Granulocyte_marker.RDS"),
  Macrophage = readRDS("Macrophage_marker.RDS"),
  NK = readRDS("NK_marker.RDS"),
  Tcells = readRDS("Tcells_marker.RDS"),
  Thyrocyte = readRDS("Thyrocyte_marker.RDS")
)

pdf("fig3b.pdf",width = 15,height = 10)
for (i in names(sce.mer.split)) {
  markers <- markers[[i]] %>% group_by(cluster) %>% top_n(5,avg_log2FC) %>% .$gene %>% unique()
  p <- DotPlot(object = sce.mer.split[[i]],features = markers,group.by = "cellIdent_sub")+ 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+
    scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red")
  print(p)
}
dev.off()


cellratio <- sce.mer@meta.data
cellratio$cancerType <- factor(cellratio$cancerType,levels = c("ATC","LPTC","PTC","Normal"))

pdf("fig3c.pdf",height = 8,width = 10)
for (i in unique(cellratio$cellIdent)) {
  p <- ggplot() + 
    geom_bar(data = cellratio[cellratio$cellIdent == i,], aes(x = cancerType, fill = cellIdent_sub),position = "fill") + 
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,size=15,hjust = 1),
          axis.title.x = element_blank(),
          axis.text.y=element_text(size=16), 
          axis.title.y=element_text(size = 20),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size=1), 
          legend.text=element_text(size=16),
          legend.title=element_text(size=18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_manual(values = pal_d3("category20")(20))+
    ylab("Cell Ratio")+ggtitle(i)
  print(p)
}
dev.off()

gsea_table <- readRDS("gsea_table.RDS")

pdf("fig3e.pdf",height = 18,width = 14)
for (j in names(gsea_table)) {
  ord <- unique(gsea_table[[j]]$celltype)
  gsea_table[[j]]$celltype <- factor(gsea_table[[j]]$celltype,levels = ord)
  gsea_table[[j]] <- gsea_table[[j]][order(gsea_table[[j]]$ID,gsea_table[[j]]$celltype),]
  if (i == "go") {
    fp <- gsea_table[[j]] %>% group_by(celltype) %>% top_n(-10,p.adjust)
  } else {
    fp <- gsea_table[[j]]
  }
  pw <- as.data.frame(table(fp$ID))
  pw_ord <- ddply(fp,"ID",function(x)x[1,"celltype",drop = F])
  pw <- left_join(pw,pw_ord,by = c("Var1"="ID")) %>% .[order(.$Freq,.$celltype),]
  fp$ID <- factor(fp$ID,levels = pw$Var1)
  p <- ggplot(fp, aes(x = celltype, y = ID)) + 
    geom_point(aes(size = -(log10(p.adjust)), color = NES)) +
    theme_bw(base_size = 14) +
    scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") +
    ylab(NULL)+ggtitle(i) + 
    theme(axis.text.x = element_text(angle = 45,hjust = 1))
  print(p)
}
dev.off()

tds <- c("DIO1","DIO2","DUOX1","DUOX2","FOXE1","GLIS3","NKX2-1","PAX8",
         "SLC26A4","SLC5A5","SLC5A8","TG","THRA","THRB","TPO","TSHR")
Thyrocyte <- readRDS("Thyrocyte.RDS")
Thyrocyte <- AddModuleScore(Thyrocyte,list(TDS=tds),name = "TDS")
BRS <- list(BRAF=readLines("BRAF.txt"),
            RAS=readLines("RAS.txt"))
Thyrocyte <- AddModuleScore(Thyrocyte,list(BRAF=BRS[["BRAF"]]),name = "BRAF")
Thyrocyte <- AddModuleScore(Thyrocyte,list(RAS=BRS[["RAS"]]),name = "RAS")

score <- Thyrocyte@meta.data[,c("cellIdent_sub","TDS1","BRAF1","RAS1")]

pdf("fig3fgh.pdf",height = 5,width = 8)
for (i in c("TDS1","BRAF1","RAS1")) {
  score_sub <- score[,c("cellIdent_sub",i)]
  score_sub <- score_sub[order(score_sub$cellIdent_sub),]
  fills <- pal_d3("category20")(length(unique(score_sub$cellIdent_sub)))
  names(fills) <- unique(score_sub$cellIdent_sub)
  lv <- ddply(score_sub,"cellIdent_sub",function(x)median(x[,i])) %>% arrange(V1)
  lv$cellIdent_sub <- as.character(lv$cellIdent_sub)
  score_sub$cellIdent_sub <- factor(score_sub$cellIdent_sub,levels = lv$cellIdent_sub)
  p <- ggplot(score_sub,aes(cellIdent_sub,score_sub[,i],fill=cellIdent_sub)) + 
    geom_boxplot(outlier.shape = NA) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = fills[lv$cellIdent_sub]) +
    ylab(substr(i,1,nchar(i) - 1))
  print(p)
}
dev.off()
