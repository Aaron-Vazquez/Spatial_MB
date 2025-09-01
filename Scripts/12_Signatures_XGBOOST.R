rm(list = ls())

library(Giotto)
library(Seurat)
library(RColorBrewer)
library(stringr)
library(dplyr)
library(ggplot2)
library(sparrow)
library(pheatmap)
library(dplyr)
library(plotly) 
library(stats) 

fls <- list.files("./Sobject_MB/")
fls <- fls[grep("^09|^B",fls)]
tags <- gsub("\\..*","",fls)



SHH_genes = read.csv("./Results/ML/SHAP_values/SHAP_SHH_Important.csv")[,2]
G3_genes = read.csv("./Results/ML/SHAP_values/SHAP_Grp3_Important.csv")[,2]
G4_genes = read.csv("./Results/ML/SHAP_values/SHAP_Grp4_Important.csv")[,2]

#####################
##### Only positive SHAP Values according to dependency plot

SHH_genes <- SHH_genes[c(1,3,5,8)] 
G4_genes <- G4_genes[c(1:4)]
G3_genes <- G3_genes[c(1,5,9,10)]


#######################
python_path = 'C:/ProgramData/miniconda3'
save_dir <- './Results/ML/figures/'

for (i in 1:length(tags)){
  
  Data.MB <- readRDS(paste0("./Sobject_MB/",fls[i]))
  

  SHH_genes_1 = SHH_genes[which(SHH_genes %in% rownames(Data.MB))]
  G3_genes_1 = G3_genes[which(G3_genes %in% rownames(Data.MB))]
  G4_genes_1 = G4_genes[which(G4_genes %in% rownames(Data.MB))]
  
  ref = list(SHH = SHH_genes_1, G3=G3_genes_1, G4=G4_genes_1)
  
  marker_matrix <- makeSignMatrixPAGE(names(ref), ref)
  
  #######################
  raw_expr <- Data.MB@assays[["SCT"]]@counts
  norm_expr <- Data.MB@assays[["SCT"]]@data
  scale_expr <- Data.MB@assays[["SCT"]]@scale.data
  cell_meta <- Data.MB@meta.data
  image_info <- Data.MB@images
  coords <- image_info[[1]]@coordinates
  #####
  giotto_coords <- coords[,c('imagerow', 'imagecol')]
  colnames(giotto_coords) <- c('row_pxl', 'col_pxl')
  giotto_coords[,2] <- -giotto_coords[,2]
  myinst=createGiottoInstructions(save_plot=T, show_plot=T, 
                                  save_dir = save_dir, python_path=python_path)
  
  
  giotto <- createGiottoObject(raw_exprs = raw_expr, 
                               norm_expr = norm_expr,
                               norm_scaled_expr = scale_expr,
                               spatial_locs = giotto_coords, 
                               instructions = myinst, 
                               cell_metadata = cell_meta)
  
  enrich_results <- runPAGEEnrich(giotto, marker_matrix, expression_values='normalized',
                                  p_value=F, n_times=100, min_overlap_genes =2,return_gobject = F)
  
  cell_types = colnames(marker_matrix)
  
  
  enrich_scores <- as.data.frame( enrich_results$matrix )
  rownames(enrich_scores) <- enrich_scores[,1]
  enrich_scores <- enrich_scores[,!(colnames(enrich_scores) %in% c('cell_ID'))]
  colnames(enrich_scores) <- str_c(colnames(enrich_scores), 'enrich scores', sep=' ')
  
  # Replace inf values with max not inf #
  for (coli in 1:ncol(enrich_scores)) {
    isinf <- is.infinite(enrich_scores[,coli])
    enrich_scores[isinf,coli] <- max(enrich_scores[isinf==F,coli])
  }
  

  #enrich_scores <- as.data.frame( giotto$matrix )
  colnames(enrich_scores) = str_remove(colnames(enrich_scores)," enrich scores")

  enrich_scores_scld = as.data.frame(scale_rows(enrich_scores))
pheatmap(as.matrix(enrich_scores), scale = "row", show_rownames = F, filename = paste0(save_dir,tags[i],"_heatmap_shap.png"))
  
  Data.MB@meta.data$SHH <- enrich_scores_scld$SHH[match(colnames(Data.MB),rownames(enrich_scores_scld))]
  Data.MB@meta.data$G4 <- enrich_scores_scld$G4[match(colnames(Data.MB),rownames(enrich_scores_scld))]
  Data.MB@meta.data$G3 <- enrich_scores_scld$G3[match(colnames(Data.MB),rownames(enrich_scores_scld))]

  
  SpatialFeaturePlot(Data.MB,
                     features = c("G4","SHH","G3"), slot = "data", image.alpha = 0.5, keep.scale = "all", alpha = 1, pt.size.factor = 2.2, stroke = 0.25, image.scale = "lowres")+ 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  ggsave(paste0(save_dir,tags[i],"_spatial_shap.png"), width = 7, height = 7,
         dpi=500)

  
  
}








#library(umap)
#iris.umap = umap(enrich_scores, n_components = 2, random_state = 15) 
