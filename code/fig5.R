setwd("~/Project/ST/data/SC/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

ST_merge <- readRDS("ST_merge.RDS")

norm_weights <- readRDS("RCTD/norm_weights.RDS")

coords <- lapply(ST_merge@images,function(x)x@coordinates)
LE <- ST_merge@meta.data[,"class",drop = F]
coords <- lapply(coords, function(x){
  cbind(x,LE[rownames(x),,drop = F]) %>% 
    .[rownames(.) %in% rownames(norm_weights),] %>% 
    cbind(.,norm_weights[rownames(.),]) %>% {.$col <- ((.$col / sum(.$col)) * sum(.$row));.}
})

coords <- lapply(coords, function(x){
  x$Tumor <- ifelse(x$class == "Tumor",1,0)
  x$Normal <- ifelse(x$class == "Normal",1,0)
  x[x$class %in% c("Tumor","Normal"),7:13] <- 0
  return(x)
})

for (i in names(coords)) {
  if (!("LE" %in% coords[[i]]$class)) {
    coords[[i]] <- NULL
  }else{
    coords[[i]] <- coords[[i]][coords[[i]]$class %in% c("LE"),]
  }
}

pdf("fig5b.pdf",height = 10,width = 16)
for (i in c("ATC","LPTC","PTC")) {
  coordSub <- bind_rows(coords[grep(i,names(coords))])
  total <- colSums(coordSub[,7:13]) %>% as.data.frame()
  names(total) <- "ratio"
  total$ratio <- total$ratio / sum(total$ratio)
  total$cellType <- rownames(total)
  total$cellType <- factor(total$cellType,levels = c("B","Tcells","Thyrocyte","Myeloid","Fibroblast","Endothelial","NK"))
  myLabel <- paste0(total$cellType," (",round(total$ratio * 100,2),"%)")
  p <- ggplot(total,aes(x="",y = ratio, fill = cellType)) + 
    geom_col(colour = "white") + coord_polar("y") +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) + scale_fill_nejm(label = myLabel) + 
    ggtitle(i)
  print(p)
}
dev.off()


cellInfo <- ST_merge@meta.data[rownames(norm_weights),c("CancerType","class")]
norm_weights <- cbind(cellInfo,norm_weights)

cellRatio <- norm_weights %>% group_by(CancerType,class) %>% summarise_if(is.numeric, ~ sum(.))
cellRatio <- cellRatio[cellRatio$CancerType != "N",]

tmp <- list()
for (i in c("Normal","LE","Tumor")) {
  tmp[[i]] <- cellRatio[cellRatio$class == i,]
  rowname <- tmp[[i]]$CancerType
  tmp[[i]] <- tmp[[i]][,c(-1,-2)]
  rownames(tmp[[i]]) <- rowname
  tmp[[i]] <- t(tmp[[i]])
}
cellRatio <- lapply(tmp, as.data.frame)

Roe <- list()
for (i in names(cellRatio)) {
  for (j in rownames(cellRatio[[i]])) {
    cellSub <- cellRatio[[i]][j,]
    other <- cellRatio[[i]][!rownames(cellRatio[[i]]) %in% j,] %>% colSums()
    cellSub <- rbind(cellSub,other)
    rownames(cellSub) <- c(j,"other")
    chisq <- chisq.test(cellSub)
    RoeSub <- chisq$observed / chisq$expected
    RoeSub <- cbind(RoeSub[1,,drop = F],pval=chisq$p.value) %>% as.data.frame()
    Roe[[i]][[j]] <- RoeSub
  }
  Roe[[i]] <- bind_rows(Roe[[i]])
  Roe[[i]]$padj <- p.adjust(Roe[[i]]$pval,method = "BH")
}

cellRatio <- norm_weights %>% group_by(CancerType,class) %>% summarise_if(is.numeric, ~ sum(.))
cellRatio <- cellRatio[cellRatio$CancerType != "N",]

tmp <- list()
for (i in c("ATC","LPTC","PTC")) {
  tmp[[i]] <- cellRatio[cellRatio$CancerType == i,]
  rowname <- tmp[[i]]$class
  tmp[[i]] <- tmp[[i]][,c(-1,-2)]
  rownames(tmp[[i]]) <- rowname
  tmp[[i]] <- t(tmp[[i]])
}
cellRatio <- lapply(tmp, as.data.frame)

Roe <- list()
for (i in names(cellRatio)) {
  for (j in rownames(cellRatio[[i]])) {
    cellSub <- cellRatio[[i]][j,]
    other <- cellRatio[[i]][!rownames(cellRatio[[i]]) %in% j,] %>% colSums()
    cellSub <- rbind(cellSub,other)
    rownames(cellSub) <- c(j,"other")
    chisq <- chisq.test(cellSub)
    RoeSub <- chisq$observed / chisq$expected
    RoeSub <- cbind(RoeSub[1,,drop = F],pval=chisq$p.value) %>% as.data.frame()
    Roe[[i]][[j]] <- RoeSub
  }
  Roe[[i]] <- bind_rows(Roe[[i]])
  Roe[[i]]$padj <- p.adjust(Roe[[i]]$pval,method = "BH")
}

color <- c("#F8766D","#00BFC4")
names(color) <- c("padj<0.05","Not Sig")
pdf("fig5def.pdf")
for (cellType in c("B","Dendritic","Endothelial","Fibroblast","Granulocyte","Macrophage","NK","Tcells","Thyrocyte")) {
  RoeSub <- Roe[["LE"]][grep(cellType,rownames(Roe[["LE"]])),]
  padj <- cbind(Var1 = rownames(RoeSub),padj = RoeSub[,"padj",drop = F])
  RoeSub <- reshape2::melt(as.matrix(RoeSub[,1:3]))
  RoeSub <- left_join(RoeSub,padj,by = "Var1")
  RoeSub$group <- ifelse(RoeSub$padj < 0.05,"padj<0.05","Not Sig")
  RoeSub$group <- factor(RoeSub$group,levels = c("padj<0.05","Not Sig"))
  p <- ggplot(RoeSub,aes(x=Var1,y=value)) + 
    geom_segment(aes(x=Var1,xend=Var1,y=0,yend=value),
                 size=1.5,color="#C9CACA",linetype="solid")+
    geom_point(aes(fill = group),size=5,shape=21) + 
    scale_fill_manual(values = color[unique(RoeSub$group)]) +
    coord_flip() + 
    facet_grid(~Var2) + 
    theme_classic() + 
    theme(axis.title = element_blank())
  print(p)
}
dev.off()

markers <- readRDS("leadingEdge/markers_LE.RDS")

pdf("fig5g.pdf")
for (i in names(markers)) {
  for (j in unique(markers[[i]]$cluster)) {
    markerSub <- filter(markers[[i]],cluster == j)
    markerSub <- markerSub %>% filter(p_val_adj < 0.05) %>% arrange(-avg_log2FC) %>% 
      mutate(rank = 1:nrow(.)) %>% mutate(trend = ifelse(avg_log2FC > 0,"up","down"))
    pointOut <- rbind(head(markerSub,5),tail(markerSub,5))
    p <- ggplot() + 
      geom_point(data = markerSub,mapping = aes(x=rank,y=avg_log2FC,color = avg_log2FC))+
      scale_color_gradient2(low = "navy",high = "firebrick3",mid = "white",midpoint = 0)+
      geom_point(data = pointOut,mapping = aes(x=rank,y=avg_log2FC),size = 3,color = "black")+
      geom_point(data = pointOut,mapping = aes(x=rank,y=avg_log2FC),size = 2,color = "yellow")+
      ggrepel::geom_text_repel(data = pointOut,mapping = aes(x=rank,y=avg_log2FC,label = gene),size = 3)+
      theme_classic() + ggtitle(paste(i,j))
    print(p)
  }
}
dev.off()

upreg <- lapply(markers, function(x){
  x[x$cluster == "LE" & x$avg_log2FC > log2(1.5) & x$p_val_adj < 0.05,"gene"]
})

downreg <- lapply(markers, function(x){
  x[x$cluster == "LE" & x$avg_log2FC < -log2(1.5) & x$p_val_adj < 0.05,"gene"]
})

pdf("fig5h.pdf",height = 6,width = 6)
ggvenn(upreg)+ggtitle("UpRegulate: avg_log2FC > log2(1.5) & p_val_adj < 0.05")
ggvenn(downreg)+ggtitle("DownRegulate: avg_log2FC < -log2(1.5) & p_val_adj < 0.05")
dev.off()

gsea_hallmark <- readRDS("/leadingEdge/gsea_hallmark.RDS")
pdf("fig5i.pdf",height = 20,width = 11)
for (i in names(gsea_hallmark)) {
  plotlist <- list()
  for (j in rownames(gsea_hallmark[[i]][["LE"]]@result)) {
    NES <- gsea_hallmark[[i]][["LE"]]@result %>% .[j,"NES"]
    gsdata <- gsInfo(gsea_hallmark[[i]][["LE"]], j)
    k <- 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                     term)
      gsdata[idx, "ymin"] <- k
      gsdata[idx, "ymax"] <- k + 1
      k <- k + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) + 
      geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax), color = ifelse(NES > 0,"red","blue")) + 
      xlab(NULL) + ylab(NULL) + 
      theme_bw() + 
      theme(legend.position = "none", 
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            plot.background = element_blank())
    plotlist[[j]] <- plot_grid(plotlist = list(gridExtra::tableGrob(NES),p2,gridExtra::tableGrob(j)), 
                               ncol = 3,rel_widths = c(1,2,3)) + 
      theme(plot.background = element_rect(color = "black"))
  }
  plotlist <- plot_grid(plotlist = plotlist,ncol = 1)
  titles <- ggdraw() + draw_label(paste(i,"LE"), fontface='bold')
  plotlist <- plot_grid(titles,plotlist,rel_heights = c(0.05,1),ncol = 1)
  print(plotlist)
}
dev.off()
