### Load shap values for global analyzed data
###
rm(list=ls())

multimerge <- function (mylist) {
  ## mimics a recursive merge or full outer join
  
  unames <- unique(unlist(lapply(mylist, rownames)))
  
  n <- length(unames)
  
  out <- lapply(mylist, function(df) {
    
    tmp <- matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] <- as.matrix(df)
    rm(df); gc()
    
    return(tmp)
  })
  
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
  
  bigout <- do.call(cbind, out)
  colnames(bigout) <- paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(bigout)
}

library(ggplot2) 

MB.labs <- c("SHH","Group3","Group4")
MB.cols <- c("gray40","firecrick1","cyan4")

for (i in c(1:3)){
  SHAP <- lapply(0:4, function(x){ 
    y = read.csv(paste0("./Results/ML/SHAP_values/SHAP_",MB.labs[i],"_",x,".csv"))
    y$Gene <- gsub("',\\)","", y$Gene)
    y$Gene <- gsub("\\('","", y$Gene)
    y=y[-1]
    write.csv(y,paste0("./Results/New_Macro/Results/ML/SHAP_values/SHAP_Grupo4_cln_",
                     as.character(x),".csv"),row.names = F)
  })
  
  
  df <- multimerge(SHAP)
  df <- df[,-c(1,3,5,7,9)]
  colnames(df) <- c("Rn_1","Rn_2","Rn_3","Rn_4","Rn_5")
  #df = as.numeric(df)
  
  
  df = as.data.frame(df)
  df[is.na(df)] <- 0
  df2 <- data.frame(apply(df, 2, function(x) as.numeric(x)))
  rownames(df2) <- rownames(df) 
  avg = rowMeans(df2)
  std = apply(df2,1,sd)
  SHAP_st = data.frame(gene = rownames(df),
                       avg = rowMeans(df2),
                       std = apply(df2,1,sd))
  
  SHAP_st <- SHAP_st[order(SHAP_st$avg,decreasing=TRUE),]
  #### SHH 0.1
  #### Grp3 0.01
  #### Grp 0.1
  SHAP_flt = SHAP_st[SHAP_st$avg>=0.01,]
  #### SHH gray40
  #### Grp4 firecrick1
  #### Grp3  cyan4
  
  
  
  pp <- ggplot(SHAP_flt, aes(x=factor(gene,levels = rev(SHAP_st$gene)), y=avg) )+
    geom_bar(stat="identity", color="black", fill=MB.cols[i])+
    theme(text = element_text(size=20))+
    labs(y = "Mean SHAP value", x = "") +
    geom_errorbar(aes(ymin=avg, ymax=avg+std), 
                  width=.2,position=position_dodge(.9))+coord_flip()
  
  ggsave(paste0("./Results/ML/figures/SHAP_",MB.labs[i],".png"),
         plot = pp, device = "png",dpi = 300, height = 6, width = 6.5, 
         units = "in",bg = "white")
  
  write.csv(SHAP_flt$gene, 
            file = paste0("./Results/ML/SHAP_values/SHAP_",MB.labs[i],"_Important.csv"))
}
