library(tidyverse)
library(Seurat)
library(harmony)

dir_name <- dir("/data/ST/data/")
dir_path <- paste("/data/ST/data",dir_name,sep = "/")

names(dir_path) <- dir_name
sce <- lapply(dir_path, Read10X)
for (i in names(sce)) {
  sce[[i]] <- CreateSeuratObject(counts = sce[[i]],
                                 min.cells = 3, 
                                 min.features = 200,
                                 project = i )
  sce[[i]][["percent.mt"]] <- PercentageFeatureSet(object = sce[[i]], 
                                                   pattern = "^MT-")
  sce[[i]] <- RenameCells(sce[[i]],add.cell.id = i)
  sce[[i]] <- subset(
    sce[[i]],
    subset =  nFeature_RNA > 200 &
      nCount_RNA > 1000 & 
      percent.mt < 20
  )
  sce[[i]]<- NormalizeData(sce[[i]]) %>% 
    ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% 
    FindNeighbors() %>% 
    FindClusters()
}

sce.mer <- Reduce(merge,sce)
sce.mer <- merge(sce.mer,c(GSE184362_PTC_merge,GSE148673_ATC_merge))

sce.mer <- NormalizeData(sce[[i]]) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  harmony_pip(.,group.by = c("orig.ident","batch","tissue_source")) %>% 
  RunUMAP(.,reduction = "harmony",dims = 1:30) %>% 
  FindNeighbors(., reduction = "harmony", dims = 1:30) %>% 
  FindClusters() %>% RunTSNE(., reduction = "harmony", dims = 1:30)

saveRDS(sce.mer,"/data/ST/data/sce.mer.RDS")