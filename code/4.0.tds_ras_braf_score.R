library(tidyverse)
library(SeuratObject)

tds <- c("DIO1","DIO2","DUOX1","DUOX2","FOXE1","GLIS3","NKX2-1","PAX8",
         "SLC26A4","SLC5A5","SLC5A8","TG","THRA","THRB","TPO","TSHR")
Thyrocyte <- readRDS("~/ST/data/SC/subcluster/Thyrocyte.RDS")
Thyrocyte <- AddModuleScore(Thyrocyte,list(TDS=tds),name = "TDS")
BRS <- list(BRAF=readLines("~/ST/data/SC/BRAF.txt"),
            RAS=readLines("~/ST/data/SC/RAS.txt"))
Thyrocyte <- AddModuleScore(Thyrocyte,list(BRAF=BRS[["BRAF"]]),name = "BRAF")
Thyrocyte <- AddModuleScore(Thyrocyte,list(RAS=BRS[["RAS"]]),name = "RAS")

score <- Thyrocyte@meta.data[,c("cellIdent_sub","TDS1","BRAF1","RAS1")]

pdf("~/ST/result/SC/plot/fig3fgh.pdf",height = 5,width = 8)
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


####correlation of tranjactory and TDS/BRS
cds <- readRDS("~/ST/data/SC/trajectory/Thyrocyte/monocle3.rds")
Thyrocyte <- readRDS("~/ST/data/SC/subcluster/Thyrocyte.RDS")
pseudo <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]] %>% as.data.frame()
score <- Thyrocyte@meta.data[,c("cellIdent_sub","TDS1","BRAF1","RAS1","nCount_RNA")]
score <- cbind(score,pseudo)
names(score)[6] <- "pseudoTime"
score <- score[order(score$nCount_RNA),]

pdf("~/ST/result/SC/plot/fig4cde.pdf",height = 5,width = 6)
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
