rm(list=ls())

library(ggplot2)
library(reshape2)

fls <- list.files("./Sobject_MB")
fls <- fls[grep("^G|^S",fls)]
tags <- gsub("[:.:]rds$","",fls)
Grps <- sub("-.*", "", tags)

save_dir <- './Results'


top20 <- matrix(data = NA, nrow =length(tags) , ncol = 3, 
                dimnames = list(tags, c("G3","G4","SHH")))

for (i in 1:length(tags)){
  Dat <- read.csv(paste0(save_dir,tags[i],"_scores_signature_1.csv"))
  nr <- nrow(Dat)
  Dat$X[nr] <- "Top20"
  top20[tags[i],"G3"] <- Dat$G3[nr]
  top20[tags[i],"G4"] <- Dat$G4[nr]
  top20[tags[i],"SHH"] <- Dat$SHH[nr]
  Dat <- melt(Dat,value.name = "Prop", variable.name = "Group")
  
  Dat$Prop = Dat$Prop*100 
  
  ggplot(data=Dat, aes(x=X, y=Prop, fill=Group)) +
    geom_bar(stat="identity", position=position_dodge())+ylim(0,100)+
    labs(x=bquote(~Log[2] ~ "(Fold Change)"), 
         y="% of Spots",
         title = paste(tags[i]," ","True Label: ",Grps[i]))+
    scale_fill_manual(values=c("#E69F00", "#00B2EE", "#EE2C2C"))+
    scale_x_discrete(labels=c("log2FC_0.8"= "0.8","log2FC_1"="1.0",
                              "log2FC_1.5"="1.5","log2FC_2"="2.0",
                              "log2FC_2.5"="2.5","log2FC_3"="3.0",
                              "log2FC_3.5"="3.5","log2FC_4"="4.0",
                              "Top20"= "Top20"))+
    labs(fill = "")+ theme_classic()+
    theme(legend.text=element_text(size=12),
          axis.title=element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          legend.position="bottom")
  
  ggsave(filename = paste0(save_dir,"BarPlot_",tags[i],".png"),
         width = 4, height = 3)
}

library(dplyr)
top20 <- cbind(tags,as.data.frame(top20))
top20$tags <- tags
top20.1 <- melt(top20,value.name = "Prop", variable.name = "Group")
top20.1$Prop <- top20.1$Prop*100

ggplot(data=top20.1, aes(x=tags, y=Prop, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+ylim(0,85)+
  labs(x="", 
       y="% of Spots")+
  scale_fill_manual(values=c("#E69F00", "#00B2EE", "#EE2C2C"))+
  scale_x_discrete(labels=c("log2FC_0.8"= "0.8","log2FC_1"="1.0",
                            "log2FC_1.5"="1.5","log2FC_2"="2.0",
                            "log2FC_2.5"="2.5","log2FC_3"="3.0",
                            "log2FC_3.5"="3.5","log2FC_4"="4.0",
                            "Top20"= "Top20"))+
  labs(fill = "")+ theme_classic()+
  theme(legend.text=element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        legend.position="bottom")+coord_flip()

ggsave(filename = paste0(save_dir,"BarPlot_top20.png"),
       width = 4, height = 3)
