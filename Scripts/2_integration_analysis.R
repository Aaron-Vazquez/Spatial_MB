rm(list = ls())

library(Seurat)
library(ggplot2)

#####
fls <- list.files("./")
tags <- gsub("09062023[:_:]","\\1",fls)
tags <- gsub("[_]","-", tags)
MetDat <- read.csv("./Metadata.csv", row.names = 1)

Data.MB <- lapply(1:length(tags), function(x) {
  print(x)
  Data <- readRDS(paste("./Sobject_MB/",tags[x],".rds", sep = ""))
  Data@meta.data[["orig.ident"]] <- tags[x]
  Idents(object = Data) <- "orig.ident"
#### compute percentage of mitochondrial genes
  Data <- PercentageFeatureSet(Data, "^MT-", col.name = "percent_mito")
  VlnPlot(Data, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 2) + NoLegend()
  Data <- Data[, Data$nFeature_Spatial > 500 & Data$percent_mito < 25 & 
                 Data$percent_hb < 20]
  
  C = Data@assays$Spatial@counts
  C@x = C@x/rep.int(colSums(C), diff(C@p))
  most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
  boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
          col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
  
  Data <- SCTransform(Data, assay = "Spatial", verbose = FALSE, method = "poisson")
  return(Data)
})

names(Data.MB) <- tags

options(future.globals.maxSize = 10000 * 1024^2)  # set allowed size to 2K MiB
MB.features = SelectIntegrationFeatures(Data.MB, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = Data.MB, anchor.features = MB.features,
                              verbose = FALSE)

int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
                                      verbose = FALSE, anchor.features = MB.features)
rm(Data.MB, st.list)
gc()

MB.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
                               verbose = FALSE)

rm(int.anchors, st.list)
gc()
saveRDS(Data.MB, "./Sobject_MB/MB_integrated.rds")

