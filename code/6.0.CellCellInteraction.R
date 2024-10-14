setwd("~/Project/ST/data/SC/")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)

ST_merge <- readRDS("ST_merge.RDS")

lr_merge <- readRDS("LR/lr_merge.RDS")
names(lr_merge) <- c("Ligand","Receptor")
rownames(lr_merge) <- paste0(lr_merge$Ligand,"_",lr_merge$Receptor)
coords <- lapply(ST_merge@images, function(x)x@coordinates)

Exp_tot <- as.matrix(ST_merge@assays[["SCT"]]@data)

lr_cor <- apply(lr_merge, 1, function(lr){
  Exp <- t(Exp_tot[as.character(lr),])
  cells <- rownames(Exp[Exp[,1] > 0,])
  names(cells) <- cells
  neighbor <- lapply(cells, function(x){
    samples <- as.character(ST_merge@meta.data[x,"samples"])
    rows <- coords[[samples]][x,"row"]
    cols <- coords[[samples]][x,"col"]
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
    neighbor[[i]]$core <- i
    neighbor[[i]]$neighbor <- rownames(neighbor[[i]])
  }
  
  neighbor <- bind_rows(neighbor)
  Exp <- as.data.frame(Exp)
  Exp$neighbor <- rownames(Exp)
  
  neighbor <- left_join(neighbor,Exp[,2:3],by = "neighbor")
  names(neighbor)[8] <- "receptor"
  neighbor <- plyr::ddply(neighbor,"core",function(x){
    mean(x$receptor)
  })
  names(neighbor)[2] <- "receptor"
  neighbor <- left_join(neighbor,Exp[,c(1,3)],by = c("core" = "neighbor"))
  names(neighbor)[3] <- "ligand"
  neighbor$cancerType <- substr(neighbor$core,1,regexpr("[-]",neighbor$core) - 1)
  correlation <- plyr::ddply(neighbor[neighbor$cancerType != "N",],"cancerType",function(x){
    x <- rcorr(as.matrix(x[,c("ligand","receptor")]),type = "spearman")
    return(data.frame(cor = x[["r"]][1,2],pval=x[["P"]][1,2]))
  })
  rownames(correlation) <- correlation$cancerType
  return(data.frame(ATC_cor = correlation["ATC","cor"],ATC_p = correlation["ATC","pval"],
                    LPTC_cor = correlation["LPTC","cor"],LPTC_p = correlation["LPTC","pval"],
                    PTC_cor = correlation["PTC","cor"],PTC_p = correlation["PTC","pval"]))
})

for (i in names(lr_cor)) {
  rownames(lr_cor[[i]]) <- i
}
lr_cor <- bind_rows(lr_cor)

lr_cor <- list(ATC = lr_cor[,1:2],
               LPTC = lr_cor[,3:4],
               PTC = lr_cor[,5:6])
for (i in names(lr_cor)) {
  names(lr_cor[[i]]) <- c("cor","P")
  lr_cor[[i]]$Padj <- p.adjust(lr_cor[[i]]$P,method = "BH")
  lr_cor[[i]]$cancerType <- i
}

lr_cor <- lapply(lr_cor, function(x){
  filter(x,Padj < 0.05 & cor > 0)
})
for (i in names(lr_cor)) {
  lr_cor[[i]]$LRpair <- rownames(lr_cor[[i]])
}
lr_cor <- bind_rows(lr_cor)
lr_cor$Padj <- lr_cor$Padj + min(lr_cor$Padj[lr_cor$Padj != 0])
lr_cor$Ligand <- substr(lr_cor$LRpair,1,regexpr("_",lr_cor$LRpair) - 1)
lr_cor$Receptor <- substring(lr_cor$LRpair,regexpr("_",lr_cor$LRpair) + 1)

pdf("fig6a.pdf",height = 30,width = 10)
ggplot(data = lr_cor,aes(x = "cancerType",y = LRpair)) + 
  geom_point(aes(size = -log10(Padj),color = cor)) + 
  scale_color_gradient(low = "white",high = "red") + 
  geom_point(aes(size = -log10(Padj)),shape = 1) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  facet_grid(~cancerType)
dev.off()

norm_weights <- readRDS("RCTD/norm_weights_sub.RDS")
norm_weights$barcode <- rownames(norm_weights)
weight <- reshape2::melt(norm_weights,id = "barcode")
weight <- plyr::ddply(weight,"barcode",function(x){
  x[x$value == max(x$value),]
})

ST_sub <- subset(ST_merge,cells = weight$barcode)
rownames(weight) <- weight$barcode
weight <- weight[rownames(ST_sub@meta.data),]
ST_sub$cellIdent_sub <- weight$variable
Idents(ST_sub) <- "cellIdent_sub"

ST_split <- SplitObject(ST_sub,split.by = "CancerType")

lr_merge <- readRDS("LR/lr_merge.RDS")
names(lr_merge) <- c("from","to")
rownames(lr_merge) <- paste0(lr_merge$from,"_",lr_merge$to)
Exp <- lapply(ST_split[c("ATC","LPTC","PTC")],function(x){
  AverageExpression(object = x,assays = "SCT",features = unique(unlist(lr_merge)))[["SCT"]]
})

pct <- lapply(ST_split[c("ATC","LPTC","PTC")], function(x){
  exp <- t(as.matrix(x@assays[["SCT"]]@data))[,unique(unlist(lr_merge))]
  meta <- x@meta.data[rownames(exp),"cellIdent_sub",drop = F]
  exp <- cbind(meta,exp)
  pct <- plyr::ddply(exp,"cellIdent_sub",function(xx){
    colSums(xx[,2:71] > 0) / nrow(exp)
  })
  rownames(pct) <- pct$cellIdent_sub
  pct$cellIdent_sub <- NULL
  return(t(pct))
})

cellGroup <- unique(as.data.frame(t(combn(rep(colnames(norm_weights),2),2))))
cellGroup <- cellGroup[order(cellGroup$V1,cellGroup$V2),]
rownames(cellGroup) <- paste0(cellGroup$V1,"-",cellGroup$V2)
names(cellGroup) <- c("sender","receiver")

lr_score <- list()

for (i in names(pct)) {
  group_sub <- cellGroup %>% 
    filter(sender %in% colnames(pct[[i]]) & 
             receiver %in% colnames(pct[[i]]))
  lr_score[[i]] <- apply(group_sub, 1, function(cell){
    apply(lr_merge[lr_merge$from %in% rownames(pct[[i]]) & lr_merge$to %in% rownames(pct[[i]]),], 1, function(lr){
      min(((Exp[[i]][lr[[1]],cell[[1]]])*(pct[[i]][lr[[1]],cell[[1]]])),
          ((Exp[[i]][lr[[2]],cell[[2]]])*(pct[[i]][lr[[2]],cell[[2]]])))
    })
  })
}

commu <- list()
for (i in names(lr_score)) {
  commu[[i]] <- lr_score[[i]][rownames(lr_score[[i]]) %in% lr_cor$LRpair[lr_cor$cancerType == i],] %>% 
    reshape2::melt() %>% plyr::ddply("Var2",function(xx)sum(xx$value))
}

commu_sub <- lapply(commu, function(x){
  x[grep("^Fibroblast|Macrophage",x$Var2),] %>% 
    filter(substr(Var2,1,regexpr("[-]",Var2) - 1) != substring(Var2,regexpr("[-]",Var2) + 1)) %>% 
    top_n(100,V1)
})


commu_mat <- lapply(commu_sub, function(x){
  x <- x %>% mutate(sender = substr(Var2,1,regexpr("[-]",Var2) - 1)) %>% 
    mutate(receiver = substring(Var2,regexpr("[-]",Var2) + 1)) %>% .[,2:4] %>%
    reshape2::dcast(sender~receiver,value.var = "V1")
  rownames(x) <- x$sender
  x$sender <- NULL
  return(x)
})

commu_mat <- lapply(commu_mat, function(x){
  x[,rownames(x)[!(rownames(x) %in% colnames(x))]] <- NA
  x[colnames(x)[!(colnames(x) %in% rownames(x))],] <- NA
  x[is.na(x)] <- 0
  x <- x[order(rownames(x)),]
  x <- x[,order(colnames(x))]
  return(as.matrix(x))
})

markers <- readRDS("deg/markers_vs_normal.RDS")
for (i in names(markers)) {
  markers[[i]]$gene <- rownames(markers[[i]])
  markers[[i]]$cancerType <- i
}

markerL <- bind_rows(markers) %>% filter(gene %in% lr_cor$Ligand)
markerL$p_val_adj <- markerL$p_val_adj + min(markerL$p_val_adj[markerL$p_val_adj != 0])
markerR <- bind_rows(markers) %>% filter(gene %in% lr_cor$Receptor)
markerR$p_val_adj <- markerR$p_val_adj + min(markerR$p_val_adj[markerR$p_val_adj != 0])
markerR$gene <- factor(markerR$gene,levels = unique(lr_cor$Receptor[order(lr_cor$LRpair)]))

pdf("figbbc.pdf",height = 20,width = 10)
ggplot(markerL,aes(x = abs(avg_log2FC),y = gene)) + 
  geom_col(aes(fill = avg_log2FC)) + 
  scale_fill_gradient2(low = "blue", mid = "#D9D9D9", high = "red") + 
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.background = element_blank()) + 
  facet_grid(~cancerType) + ylab("Ligand")
ggplot(markerR,aes(x = abs(avg_log2FC),y = gene)) + 
  geom_col(aes(fill = avg_log2FC)) + 
  scale_fill_gradient2(low = "blue", mid = "#D9D9D9", high = "red") + 
  theme_bw() + theme(panel.grid = element_blank(),
                     panel.background = element_blank()) + 
  facet_grid(~cancerType) + ylab("Receptor")
dev.off()

pdf("fig6hij.pdf",height = 8,width = 8)
for (i in names(commu_mat)) {
  p <- CellChat::netVisual_circle(as.matrix(commu_mat[[i]]), 
                                  vertex.weight = rep(1,nrow(commu_mat[[i]])), 
                                  weight.scale = T, 
                                  label.edge= F, 
                                  title.name = i)
  print(p)
}
dev.off()
