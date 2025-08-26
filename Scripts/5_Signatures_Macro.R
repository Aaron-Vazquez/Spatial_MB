rm(list = ls())

combine_plots <- function(gg, heatmap_obj, plotly_obj, 
                          output_file = "combined_2x2_plot.png",
                          width = 800, height = 600, dpi = 300,
                          label_fontsize = 40, label_offset = 20,tag, grp) {
  library(ggplot2)
  library(pheatmap)
  library(plotly)
  library(webshot2)
  library(htmlwidgets)
  library(magick)
  library(grid)
  
  #  webshot2::install_phantomjs()
  
  # Define sizes
  full_width <- width * 2
  half_width <- width
  row_height <- height
  
  ggplot_height <- row_height * 1.3
  # Temp files
  gg_file <- tempfile(fileext = ".png")
  heatmap_file <- tempfile(fileext = ".png")
  html_file <- tempfile(fileext = ".html")
  plotly_file <- tempfile(fileext = ".png")
  
  # Save ggplot (A)
  ggsave(gg_file, plot = gg, width = full_width / dpi, height = ggplot_height / dpi, dpi = dpi)
  
  # Save pheatmap (B)
  png(heatmap_file, width = half_width, height = row_height)
  grid.newpage()
  grid.draw(heatmap_obj$gtable)
  dev.off()
  
  # Save plotly (C)
  saveWidget(plotly_obj, html_file, selfcontained = TRUE)
  webshot(html_file, plotly_file, vwidth = half_width, vheight = row_height)
  
  # Read and resize
  img_gg <- image_read(gg_file) %>%
    image_resize(paste0(full_width, "x", round(ggplot_height), "!"))
  
  img_heat <- image_read(heatmap_file) %>%
    image_resize(paste0(half_width, "x", row_height, "!"))
  
  img_plotly <- image_read(plotly_file) %>%
    image_resize(paste0(half_width, "x", row_height, "!"))
  
  # Annotate each image with A, B, C
  label_img <- function(img, label) {
    image_annotate(img, label, size = label_fontsize, color = "black",
                   location = paste0("+", label_offset, "+", label_offset),
                   gravity = "northwest")
  }
  
  img_gg <- label_img(img_gg, "A")
  img_heat <- label_img(img_heat, "B")
  img_plotly <- label_img(img_plotly, "C")
  
  # Combine into rows
  row1 <- img_gg
  row2 <- image_append(c(img_heat, img_plotly))
  
  # Add vertical padding above all rows for the title
  title_padding <- 100  # height of white space for title
  total_width <- image_info(row1)$width
  
  blank_bar <- image_blank(width = total_width, height = title_padding, color = "white")
  
  # Stack: title padding → ggplot → heatmap+plotly
  final_image <- image_append(c(blank_bar, row1, row2), stack = TRUE)
  
  # Now add left-aligned title in the padding area
  final_image <- image_annotate(
    final_image,
    text = paste0(tag," True label:", grp),  # your title here
    gravity = "northwest",          # top-left corner
    location = "+40+20",            # offset from corner: right and down
    size = 60,
    color = "black",
    weight = 700
  )  # Save
  image_write(final_image, output_file)
  message("✅ Combined 2x2 plot with labels saved to: ", output_file)
}


library(Giotto)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggplotify)
library(plotly) 
library(stats)
library(cowplot)
library(notly)

Gns <- read.csv("./Results/signature.csv",header = T)

Gns %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  dplyr::filter(avg_log2FC >= 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20

Grps <- unique(Gns$cluster)

fls <- list.files("./Sobject_MB")
fls <- fls[grep("^G|^S",fls)]
tags <- gsub("[:.:]rds$","",fls)

SHH_genes = top20$gene[ which(top20$cluster == "SHH")]
G3_genes = top20$gene[ which(top20$cluster == "Group 3")]
G4_genes = top20$gene[ which(top20$cluster == "Group 4")]

#######################
python_path = 'C:/ProgramData/miniconda3'
save_dir <- './Results/'

for (i in 1:length(tags)){
  
  Data.MB <- readRDS(paste0("./Sobject_MB/",fls[i]))
  
  SHH_genes = top20$gene[ which(top20$cluster == "SHH")]
  G3_genes = top20$gene[ which(top20$cluster == "Group 3")]
  G4_genes = top20$gene[ which(top20$cluster == "Group 4")]
  
  SHH_genes = SHH_genes[which(SHH_genes %in% rownames(Data.MB))]
  G3_genes = G3_genes[which(G3_genes %in% rownames(Data.MB))]
  G4_genes = G4_genes[which(G4_genes %in% rownames(Data.MB))]
  
  ref = list(SHH = SHH_genes, G3=G3_genes, G4=G4_genes)
  
  marker_matrix <- makeSignMatrixPAGE(names(ref), ref)
  
  #######################
  
  raw_expr <- Data.MB@assays[["SCT"]]@counts
  norm_expr <- Data.MB@assays[["SCT"]]@data
  scale_expr <- Data.MB@assays[["SCT"]]@scale.data
  cell_meta <- Data.MB@meta.data
  image_info <- Data.MB@images
  coords <- image_info[[1]]@coordinates
#################

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
                                  p_value=F, n_times=100, min_overlap_genes =2, 
                                  return_gobject = F)
  

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
  bin_enrich <-  matrix(unlist(bin_enrich), ncol = 3, byrow = TRUE)
  
 fig2 <-pheatmap(as.matrix(enrich_scores), scale = "row", show_rownames = F, 
                 silent = TRUE,
                 angle_col=0,fontsize = 20,         # overall font size
                 fontsize_row = 20,     # row label font
                 fontsize_col = 20,     # column label font
                 fontsize_number = 10)
  
 
  
  ### scaled
  Data.MB@meta.data$SHH <- enrich_scores_scld$SHH[match(colnames(Data.MB),rownames(enrich_scores_scld))]
  Data.MB@meta.data$G4 <- enrich_scores_scld$G4[match(colnames(Data.MB),rownames(enrich_scores_scld))]
  Data.MB@meta.data$G3 <- enrich_scores_scld$G3[match(colnames(Data.MB),rownames(enrich_scores_scld))]

  
  fig1 <- SpatialFeaturePlot(Data.MB,
                     features = c("G4","SHH","G3"), slot = "data", image.alpha = 0.5, keep.scale = "all", alpha = 1, pt.size.factor = 5, stroke = 0.25, image.scale = "lowres", ncol = 3)+ 
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  

  
  
  X <- enrich_scores
  
  axis = list(showline=FALSE, 
              
              zeroline=FALSE, 
              
              gridcolor='#ffff', 
              
              ticklen=4)
  
  fig3 <-  enrich_scores %>%  
    
    plot_ly()  %>%  
    
    add_trace(  
      
      type = 'splom',  
      
      dimensions = list( 
        
        list(label = 'G3',values=~G3),  
        
        list(label = 'G4',values=~G4),  
        
        list(label ='SHH',values=~SHH)  ))  %>% 
    layout( xaxis = list(tickfont = list(size = 20)),
            yaxis = list(tickfont = list(size = 20)), 
            font = list(size = 20))
  
  combine_plots(fig1, fig2, fig3, 
                output_file = paste0(save_dir,"Supp_",tags[i],".png"),
                tag = tags.blind[i], grp = Grps[i])
}



