library(tidyverse)
library(ggsci)
library(ggpubr)

cellFraction <- readRDS("~/ST/data/cellFraction.RDS")
norm_weights <- readRDS("~/Project/ST/data/ST/RCTD/norm_weights.RDS")

ST_merge <- readRDS("~/Project/ST/data/ST/ST_merge.RDS")
sampleInfo <- ST_merge@meta.data[,c("samples","CancerType")]
sampleInfo$barcode <- rownames(sampleInfo)
dat_mer <- data.frame()

for (j in colnames(norm_weights)) {
  dat <- list()
  for (i in names(cellFraction)) {
    dat[[i]] <- cellFraction[[i]][,j,drop = F]
    colnames(dat[[i]]) <- "cellFraction"
    dat[[i]]$method <- i
    dat[[i]]$barcode <- rownames(dat[[i]])
  }
  dat[["RCTD"]] <- norm_weights[,j,drop = F]
  colnames(dat[["RCTD"]]) <- "cellFraction"
  dat[["RCTD"]]$method <- "RCTD"
  dat[["RCTD"]]$barcode <- rownames(dat[["RCTD"]])
  dat <- bind_rows(dat)
  dat <- left_join(dat,sampleInfo,by = "barcode")
  dat <- dat %>% group_by(method,CancerType,samples) %>% summarise(cellFraction = mean(cellFraction)) %>% mutate(cellType = j)
  dat$method <- factor(dat$method,levels = c("RCTD","cell2location","Tangram","DestVI","CARD","CytoSPACE","Redeconve"))
  dat_mer <- rbind(dat_mer,dat)
}

pdf("~/Project/ST/result/revise/cellFractionBarplot.pdf",height = 5,width = 8)
for (i in c("N","PTC","LPTC","ATC")) {
  p <- ggplot(dat_mer %>% filter(CancerType == i),aes(method,cellFraction,fill=method))+
    geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
    stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    facet_grid(~cellType,scales = 'free') + 
    ggtitle(i) + 
    ggsci::scale_fill_nejm()
  print(p)
}
dev.off()
