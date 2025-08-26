rm(list=ls())

library(EnhancedVolcano)
library(ggplot2)
library(cowplot)

Dat.Gns <-  read.csv("./Results/DGE_Macro.csv") 
Grps = c("SHH","Group 3", "Group 4")
DEG.3 <- Dat.Gns[which(Dat.Gns$cluster == "Group 3"),]
DEG.4 <- Dat.Gns[which(Dat.Gns$cluster == "Group 4"),]
DEG.SHH <- Dat.Gns[which(Dat.Gns$cluster == "SHH"),]

 p1 <- EnhancedVolcano(DEG.3[,c("gene", "avg_log2FC", "p_val_adj")],
                lab = DEG.3[,c("gene")],
                subtitle = "",#bquote(italic(EnhancedVolcano)),
                caption = "", #paste0("total = ", nrow(toptable), " variables"),
                xlab = bquote(~Log[2] ~ "FC"),
                ylab = bquote(~-Log[10] ~ "P"),
                FCcutoff = 1,
                pCutoff = 1e-02,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Group 3 MB')
 p2 <- EnhancedVolcano(DEG.4[,c("gene", "avg_log2FC", "p_val_adj")],
                       lab = DEG.4[,c("gene")],
                       subtitle = "",#bquote(italic(EnhancedVolcano)),
                       caption = "", #paste0("total = ", nrow(toptable), " variables"),
                       xlab = bquote(~Log[2] ~ "FC"),
                       ylab = bquote(~-Log[10] ~ "P"),
                       FCcutoff = 1,
                       pCutoff = 1e-02,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = 'Group 4 MB')
 p3 <- EnhancedVolcano(DEG.SHH[,c("gene", "avg_log2FC", 
                                  "p_val_adj")],
                       lab = DEG.SHH[,c("gene")],
                       subtitle = "",#bquote(italic(EnhancedVolcano)),
                       caption = "", #paste0("total = ", nrow(toptable), " variables"),
                       xlab = bquote(~Log[2] ~ "FC"),
                       ylab = bquote(~-Log[10] ~ "P"),
                       FCcutoff = 1,
                       pCutoff = 1e-02,
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = 'SHH MB')

 

p = plot_grid(p3, p1, p2, labels = c('', '',''), ncol = 3)

ggsave("./Results/Volcano.png", units = "in",
       height = 7, width = 21, plot = p) 


DEG.3 = DEG.3[which(DEG.3$avg_log2FC >= 1 & DEG.3$p_val_adj <= 1E-2),]

DEG.4 = DEG.4[which(DEG.4$avg_log2FC >= 1 & DEG.4$p_val_adj <= 1E-2),]

DEG.SHH = DEG.SHH[which(DEG.SHH$avg_log2FC >= 1 & DEG.SHH$p_val_adj <= 1E-2),]

write.csv(DEG.Gns,"./Results/signature.csv",quote = F, row.names = F)

library(VennDiagram)
library(RColorBrewer)

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(DEG.3$gene, DEG.4$gene, DEG.SHH$gene),
  category.names = c("Group 3" , "Group 4" , "SHH"),
  filename = './Results/venn_diagramm_DGE.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-5, 0, 180),
  cat.dist = c(0.05, 0.0, 0.0),
  cat.fontfamily = "sans",
  rotation = 1
)

#####################
##### Enrichment
#####################


library(clusterProfiler)
library(ggplot2)
library(GSEABase)
library(msigdbr)


enrich <- function (genes, db, smpl) {
  print(smpl)
  genes <- genes[which (genes$p_val_adj <= 0.05 ),]
  
  clts <- unique(genes$cluster)
  clts <- clts[order(clts)]
  
  ### report
    tst = genes
    
    tst = tst[order(tst$avg_log2FC, decreasing = TRUE),]
    gns <- tst[,c(7,2)]
    gns <- as.vector(gns$avg_log2FC)
    names(gns) <- tst$gene
    
    em <- GSEA(gns, TERM2GENE = db, by='fgsea', minGSSize = 10,
               pvalueCutoff = .05,nPerm=10000) 

    
    if (!isEmpty(em@result)){
      #################
      em@result[,'nes'] <- -em@result[,'NES']
      result <- em@result
      result <- result[order(result[,'NES'], decreasing = T),]
    }
  write.csv(results, paste0("./Results/GSEA_",smpl,".csv"), quote = F,row.names = F)
}


cats <- c('C2','C5', 'H')
subcats <- c('CP','GO:BP', "")

ref = lapply (1:3, function(x)
  msigdbr(species = "human", category=cats[x], subcategory=subcats[x]))
db_t2g <- bind_rows(ref)

db= db_t2g %>% dplyr::select(gs_name, gene_symbol)

enrich(genes = DEG.3,db,smpl ="Grp3")
enrich(genes = DEG.4,db,smpl ="Grp4")
enrich(genes = DEG.SHH,db,smpl ="SHH")
