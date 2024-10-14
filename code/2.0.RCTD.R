library(spacexr)

ST_merge <- readRDS("~/ST/data/ST/ST_merge.RDS")
total_count <- as_matrix(ST_merge@assays[["Spatial"]]@counts)
puck <- list()
for (i in unique(ST_merge$orig.ident)) {
  coords <- ST_merge@images[[i]]@coordinates[,c("col","row")]
  names(coords) <- c("xcoord","ycoord")
  counts <- total_count[,rownames(coords)]
  nUMI <- colSums(counts)
  puck[[i]] <- SpatialRNA(coords, counts, nUMI)
}

sce.mer <- readRDS("~/ST/data/SC/sce.mer.RDS")
counts <- sce.mer@assays[["RNA"]]@counts
cell_types <- Idents(sce.mer)
nUMI <- sce.mer$nCount_RNA
names(nUMI) <- rownames(sce.mer@meta.data)
reference <- Reference(counts, cell_types, nUMI,n_max_cells = 1000000)

myRCTD <- lapply(puck, function(x){
  create.RCTD(x, reference, max_cores = 20) %>% run.RCTD(., doublet_mode = 'full')
})
saveRDS(myRCTD,"~/ST/data/ST/RCTD/myRCTD.RDS")

myRCTD <- readRDS("~/ST/data/ST/RCTD/myRCTD.RDS")
norm_weights <- lapply(myRCTD, function(x){
  normalize_weights(x@results[["weights"]]) %>% 
    as.matrix() %>% as.data.frame()
}) %>% bind_rows()
saveRDS(norm_weights,"~/ST/data/ST/RCTD/norm_weights.RDS")

