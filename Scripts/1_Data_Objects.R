
rm(list=ls())
library(Seurat)

#####
fls <- list.files("./")
tags <- gsub("09062023[:_:]","\\1",fls)
tags <- gsub("[_]","-", tags)
MetDat <- read.csv("./Metadata.csv", row.names = 1)
####
outputDIR = "./Sobject_MB"
if (!dir.exists(outputDIR)) {dir.create(outputDIR)}

invisible( lapply(1:length(tags), function(x){
  
Data = Load10X_Spatial(data.dir = paste("./",fls[x],"/outs",sep =""),
                         filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial",
                         slice = tags[x],
                         filter.matrix = TRUE,
                         to.upper = FALSE,
                         image = NULL)
Data <- AddMetaData(
  object = Data,
  metadata = MetDat[tags[x],1],
  col.name = "Patient")

Data <- AddMetaData(
  object = Data,
  metadata = MetDat[tags[x],2],
  col.name = "Sample")

Data <- AddMetaData(
  object = Data,
  metadata = MetDat[tags[x],3],
  col.name = "Dx")

Data <- AddMetaData(
  object = Data,
  metadata = MetDat[tags[x],4],
  col.name = "Histo")

Data <- AddMetaData(
  object = Data,
  metadata = MetDat[tags[x],5],
  col.name = "Run")

saveRDS(Data, paste("./Sobject_MB/",tags[x],".rds", sep = ""))
}))
