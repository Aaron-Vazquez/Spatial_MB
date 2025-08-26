
rm(list=ls())

library(Seurat)
library(patchwork)
library(dplyr)
library(cowplot)
library(rhdf5)
library(Matrix)
library(stringr)

#####
setwd("~/MedulloBlastoma")
pth.Res <- "~/MedulloBlastoma/Results"
MB.integrated <- readRDS("./Data/Sobject_MB/MB_Integrated.rds"))
#####
Idents(MB.integrated) <- "GrpMol"

dir = MB.integrated@active.ident
dir.x = gsub(" ", "", x = dir)

lbs = unique(dir.x)
for (y in lbs) {
  tag = which(dir.x == y)
  dir.x[tag] = unlist(lapply(1:length(tag), function(x) paste(y,x, sep = "_")))  
}
names(dir.x) = names(dir)

dir.create("./Results/ML", recursive = TRUE)
MyData <- MB.integrated@assays$SCT@counts

###### Saving Data for the machine learning analysis
## Count Matrix
writeMM(MyData,"./Results/ML/Data_counts.txt")
## Genes
write(x = rownames(MyData), file = "./Results/ML/genes.tsv")
## Samples
write(x = dir.x, file = "./Results/ML/labels.tsv")


MyData <- MB.integrated@assays$SCT@data

###### Saving Data for the machine learning analysis
## Count Matrix
writeMM(MyData,"./Results/ML/Data_data.txt")