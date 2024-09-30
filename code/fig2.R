setwd("~/Project/ST/data/ST/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

ST_merge <- readRDS("ST_merge.RDS")

pdf("fig2a.pdf",height = 5,width = 4.5)
for (i in names(ST_merge@images)) {
  p <- SpatialDimPlot(ST_merge,images = i,crop = F,pt.size.factor = 1.2) + ggtitle(i)
  print(p)
}
dev.off()

markers <- c("IGKC","IGHM","MS4A1","CD79A","CD3D","IL32","CD3G","TRAC","CD3E","IL7R","KLRB1",
             "TG","EPCAM","KRT18","KRT19","LYZ","TYROBP","FCER1G","S100A9",
             "CD14","ACTA2","COL1A2","COL1A1","COL3A1","TAGLN","PECAM1",
             "VWF","RAMP2","CD34","CDH5","STMN1","MKI67","TOP2A","HMGB2")
pdf("fig2b.pdf",height = 5,width = 12)
DotPlot(ST_merge,features = markers)+
  scale_color_gradient2(low = "blue", mid = "#D9D9D9", high = "red") + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf("fig2c.pdf",height = 5,width = 6)
for (i in names(ST_merge@meta.data)[91:97]) {
  ST_merge@meta.data[,i] <- ifelse(is.na(ST_merge@meta.data[,i]),0,ST_merge@meta.data[,i])
  p <- FeaturePlot(ST_merge,i,order = T)+labs(colour = "Ratio")
  print(p)
}
dev.off()

norm_weights <- readRDS("RCTD/norm_weights.RDS")
spotInfo <- ST_merge@meta.data[,c("samples","clusters")]
spotInfo <- spotInfo[rownames(norm_weights),]
spotInfo <- cbind(spotInfo,norm_weights)

clusterInfo <- reshape2::melt(spotInfo[,2:9],value.name = "clusters")
names(clusterInfo) <- c("clusters","cellType","ratio")
clusterInfo <- plyr::ddply(clusterInfo,"clusters",function(x){
  plyr::ddply(x,"cellType",function(xx){
    mean(xx$ratio)
  })
})
names(clusterInfo)[3] <- "ratio"
pdf("fig2d.pdf",height = 5,width = 6)
for (i in unique(clusterInfo$clusters)) {
  subInfo <- clusterInfo[clusterInfo$clusters == i,]
  mylabel <- paste0(subInfo$cellType,"(",round(subInfo$ratio * 100,2),"%)")
  p <- ggplot(subInfo,aes(x="",y = ratio, fill = cellType)) + 
    geom_col(colour = "white") + coord_polar("y") +
    scale_fill_nejm(label = mylabel) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) + ggtitle(i)
  print(p)
}
dev.off()

sample_weight <- ddply(norm_weights,"samples",function(x){
  (colSums(x[,1:7]) / nrow(x))*100
})
sample_weight <- reshape2::melt(sample_weight,value.name = "samples")
names(sample_weight) <- c("samples","cellType","ratio")
pdf("fig2e.pdf",height = 5,width = 7)
ggplot(sample_weight,aes(x=samples,y = ratio, fill = cellType)) + 
  geom_col() + theme_classic()+scale_fill_nejm() + theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

pdf("fig2f.pdf",width=25,height=23)
for (i in names(ST_merge@images)) {
  coords <- ST_merge@images[[i]]@coordinates
  coords$barcode <- rownames(coords)
  coords <- merge(coords,norm_weights,by = "barcode")
  names(coords)[3:4]=c("y","x")
  coords$y <- max(coords$y) - coords$y
  if (max(coords$x) > max(coords$y)) {
    coords$x <- (coords$x / max(coords$x)) * max(coords$y)
  } else if (max(coords$x) < max(coords$y)) {
    coords$y <- (coords$y / max(coords$y)) * max(coords$x)
  }
  if (i == "PTC-5") {
    r <- 0.3
  }else if (i == "ATC-4") {
    r <- 0.25
  }else{
    r <- 0.5
  }
  plt <- vizAllTopics(theta = coords[,7:13], 
                      pos = coords[,3:4], 
                      topicOrder=seq(ncol(coords[,7:13])), 
                      topicCols=rainbow(ncol(coords[,7:13])), 
                      groups = NA, 
                      group_cols = NA, 
                      r = r, 
                      lwd = 0, showLegend = TRUE, plotTitle = "scatterpies")+
    scale_fill_nejm()
  print(plt)
}
dev.off()
