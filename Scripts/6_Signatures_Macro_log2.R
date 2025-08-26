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



fls <- list.files("./Sobject_MB")
fls <- fls[grep("^G|^S",fls)]
tags <- gsub("[:.:]rds$","",fls)


python_path = 'C:/ProgramData/miniconda3'
save_dir <- './Results/New_Macro/'
### DEG 
Gns <- read.csv("./Results/signature.csv",header = T)
Grps <- unique(Gns$cluster)

logFc_vals <- c(0.8,1,1.5,2,2.5,3,1)

for (i in 1:7){
  #### Giotto object
  Data.MB <- readRDS(paste0("./Sobject_MB/",fls[i]))
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
  ### Save matrix
  scores_FC = matrix(NA, nrow = length(logFc_vals), ncol = 3, 
                     dimnames = list(paste0("log2FC_",logFc_vals),
                                       c("G3","G4","SHH")))

  for (p in 7:length(logFc_vals)){
    if (p != length(logFc_vals)) {
      Gns %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC)) %>% 
      dplyr::filter(avg_log2FC >= logFc_vals[p]) %>%
  #    slice_head(n = 20) %>%
      ungroup() -> top20
    }else{
      Gns %>%
        group_by(cluster) %>%
        arrange(desc(avg_log2FC)) %>% 
        dplyr::filter(avg_log2FC >= 1) %>%
        slice_head(n = 20) %>%
        ungroup() -> top20
    }
    G4_genes = top20$gene[ which(top20$cluster == Grps[1])]
    SHH_genes = top20$gene[ which(top20$cluster == Grps[2])]
    G3_genes = top20$gene[ which(top20$cluster == Grps[3])]

    
    SHH_genes = SHH_genes[which(SHH_genes %in% rownames(Data.MB))]
    G3_genes = G3_genes[which(G3_genes %in% rownames(Data.MB))]
    G4_genes = G4_genes[which(G4_genes %in% rownames(Data.MB))]
    
    ref = list(SHH = SHH_genes, G3=G3_genes, G4=G4_genes)
    
    marker_matrix <- makeSignMatrixPAGE(names(ref), ref)
    gns_val <- ifelse(p %in% c(5,6,7,8),0,2)
    enrich_results <- runPAGEEnrich(giotto, marker_matrix, 
                                    expression_values='normalized',
                                    p_value=F, n_times=100, 
                                    min_overlap_genes =gns_val,                                    return_gobject = F)

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
  
    colnames(enrich_scores) = str_remove(colnames(enrich_scores)," enrich scores")
    
    
    
    enrich_scores_scld = as.data.frame(scale_rows(enrich_scores))
    
    bin_enrich <- lapply(1:nrow(enrich_scores_scld), function(x) {
      scr = matrix(data = NA, nrow = 1,ncol = 3)
      scr[1,which(enrich_scores_scld[x,] != max(enrich_scores_scld[x,]))] = 0
      scr[1,which(enrich_scores_scld[x,] == max(enrich_scores_scld[x,]))] = 1
      return(scr)
    })
    bin_enrich <-  matrix(unlist(bin_enrich), ncol = 3, byrow = TRUE, 
                          dimnames = list(rownames(enrich_scores_scld),
                                          colnames(enrich_scores_scld)))
    
    pheatmap(as.matrix(enrich_scores), scale = "row", show_rownames = F, 
             filename = paste0(save_dir,tags[i],"_",
                               rownames(scores_FC)[p],"_heatmap_signature.png"))
    
    scores_FC[p,"G3"] = sum(bin_enrich[,"G3"])
    scores_FC[p,"G4"] = sum(bin_enrich[,"G4"])
    scores_FC[p,"SHH"] = sum(bin_enrich[,"SHH"])
    scores_FC[p,]= scores_FC[p,]/sum(scores_FC[p,])
  }
  write.csv(scores_FC,paste0(save_dir,tags[i],"_scores_signature_1.csv"))
}
