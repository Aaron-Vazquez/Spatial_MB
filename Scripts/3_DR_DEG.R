
rm(list=ls())

library(Seurat)
library(ggplot2)
library(dplyr)

MB.integrated <- readRDS("./Sobject_MB/MB_integrated.rds")

MB.integrated <- RunPCA(MB.integrated, verbose = FALSE, seed.use = 12)
MB.integrated <- FindNeighbors(MB.integrated, dims = 1:50)
MB.integrated <- FindClusters(MB.integrated, verbose = FALSE, resolution = 0.8)
MB.integrated <- RunUMAP(MB.integrated, dims = 1:50, seed.use = 12L)

pth.Res <- "./Results"
if (!dir.exists(pth.Res)) {dir.create(pth.Res)}

DimPlot(MB.integrated, reduction = "umap", group.by = c("ident", "Dx","Histo"))
ggsave(paste0(pth.Res,"/umap.png"))

SpatialDimPlot(MB.integrated)
ggsave(paste0(pth.Res,"/umap_sp.png"))

Idents(MB.integrated) <- "Dx"
#DimPlot(MB.integrated)

MB.markers <- FindAllMarkers(MB.integrated, only.pos = F, min.pct = 0.25, 
                             logfc.threshold = 0.25)

write.csv(MB.markers,paste0(pth.Res,"./DGE_MB_GrpMol.csv",quote = FALSE, row.names = FALSE))

Idents(MB.integrated) <- "seurat_clusters"
#DimPlot(MB.integrated)

MB.markers <- FindAllMarkers(MB.integrated, only.pos = F, min.pct = 0.25, 
                             logfc.threshold = 0.25)

write.csv(MB.markers,paste0(pth.Res,"./DGE_MB_Clust.csv",quote = FALSE, row.names = FALSE))