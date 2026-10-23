DATA  <- "data"
OUT   <- file.path("output", "FigS6")
CACHE <- file.path("output", "FigS6")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

DENSITY_LEVELS   <- c("10re", "40re", "70re", "100re")
DENSITY_LABELS   <- c("10%", "40%", "70%", "100%")
TREAT_LEVELS <- c("CK", "C", "N", "pH", "T", "W")
TREAT_LABELS <- c("CK", "C", "N", "pH", "W", "D")

{
  library(reshape2)
  library(plyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)

  group = read.csv(file.path(DATA, "amplicon/prokaryote/group.csv"), header = T)
  TAX  <- read.delim(file.path(DATA, "amplicon/prokaryote/metadata.txt"), row.name = 1, check.names = FALSE)
  TAX <- TAX %>%
    mutate_all(~ ifelse(. != "Other", substr(., 4, nchar(.)), .))

  otulist = list.files(path = file.path(DATA, "amplicon/prokaryote/by_ecosystem"), pattern="*.csv")

  for (i in 1:length(otulist)) {
    name <- substr(otulist[i],1,nchar(otulist[i])-4)
    OTU <-read.csv(file.path(DATA, "amplicon/prokaryote/by_ecosystem", otulist[i]), header = T, row.names = 1)
    mycol= c("#BFBFBF","#CD3333","#DB7093","#CD853F","#6959CD","#FFC0CB","#9ACD32","#FFD700","#FFDAB9")

    OTU_Tax = merge(OTU,TAX,by="row.names")
    colnames(OTU_Tax)[98]
    ClassRSW=aggregate(x=OTU_Tax[,2:97],by=list(OTU_Tax$phylum),FUN='sum')
    Class=ClassRSW[,2:97]
    row.names(Class)=ClassRSW[,1]
    for (j in 1:nrow(ClassRSW)){
      if (ClassRSW[j,1] %in% c("Proteobacteria","Actinobacteriota","Acidobacteriota","Chloroflexi",
             "Planctomycetota","Firmicutes","Bacteroidota","Verrucomicrobiota")){

      }else{
        ClassRSW[j,1]="Other"
      }
    }

    other_rows <- ClassRSW[ClassRSW[,1] == "Other", ]

    mean_values <- colSums(other_rows[, -1])

    mean_df <- other_rows[1,]
    mean_df[colnames(other_rows)[1]] <- "Other"

    for (cn in colnames(other_rows)[-1]) {
      mean_df[cn] <- mean_values[cn]
    }

    non_other_rows <- ClassRSW[ClassRSW[,1] != "Other", ]

    Class1 <- rbind(non_other_rows, mean_df)
    colnames(Class1)[1] = 'Class'
    tax_list = c("Acidobacteriota","Actinobacteriota","Bacteroidota","Chloroflexi","Firmicutes","Planctomycetota","Proteobacteria","Verrucomicrobiota","Other")
    indices <- which(Class1[,1] %in% tax_list)
    mycol <- mycol[indices]

    Class0.01<-melt(Class1,id='Class')
    names(Class0.01)[2]<-'sample'

    Class_15 <- ddply(Class0.01,'sample',transform,percent = value/sum(value)*100)
    data = merge(Class_15,group,by = "sample")

    mean <- aggregate(data$percent, by=list(data$Class, data$sample), FUN=mean)
    colnames(mean) = c("Class", "Sample",  "Mean")
    group2 = group
    mean = merge(mean,group2,by.x = "Sample",by.y = "sample")
    mean$Mean <- mean$Mean/4

    mean$Class <- as.character(mean$Class)
    mean$Class[which(mean$Class=="Other")] <- "ZOther"
    rowno <- as.character(unique(mean$Class))
    rowno <- rowno[order(substring(rowno, 1, 1))]

    .s6_lvl <- rev(rowno)
    mean$Class<-factor(unique(mean$Class),levels=.s6_lvl,labels=sub("^Z","",.s6_lvl))
    mean$group = factor(mean$group, levels = DENSITY_LEVELS, labels = DENSITY_LABELS)
    mean$group1 = factor(mean$group1, levels = TREAT_LEVELS, labels = TREAT_LABELS)

    pALL<-ggplot(mean,aes(x=group1,y=Mean,fill=Class))+
      geom_bar(stat="identity",position="stack",width = 0.8)+
      facet_wrap(~group, scales = 'free_x', ncol = 2) +

      ggtitle("")+
      labs(x="",y="Relative abundance (%)")+
      scale_fill_manual(values = mycol)+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"),
            panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.length=unit(0.1,'cm'),
            legend.text = element_text(size = 8, color="black"),
            legend.key = element_rect( fill = "lemonchiffon"),
            legend.key.size = unit(8, "pt"),
            legend.position = "none",
            strip.text = element_text(size = 6, color="black"),
            axis.title = element_text(size = 8,color="black"),

            axis.text.x = element_text(size = 6, color="black", angle = 45,
                                       hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 6,color="black"),
            axis.title.y = element_text(colour="black", size = 8),
            axis.title.x = element_text(size = 8, color="black"),
      )+
      guides(fill=guide_legend(title=NULL))

    pALL
    outfilename <- paste0(name,"_stacked_bar.pdf")

    ggsave(outfilename, pALL, width = 4.5, height = 8.0, path = OUT, units = "cm")
    if (!exists(".figS6_leg_16S")) {

      p_leg <- pALL +
        theme(legend.key.size = unit(8, "pt")) +
        guides(fill = guide_legend(title = NULL, ncol = 2))
      .save_legend(p_leg, "FigS6_16S_stack_legend", OUT, orient = "right")
      .figS6_leg_16S <- TRUE
    }
  }

  group = read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header = T)
  TAX  <- read.delim(file.path(DATA, "amplicon/fungi/metadata.txt"), row.name = 1, check.names = FALSE)

  TAX <- TAX %>%
    mutate_all(~ ifelse(grepl("^[a-z]__", .), substr(., 4, nchar(.)), .))

  otulist = list.files(path = file.path(DATA, "amplicon/fungi/by_ecosystem"), pattern="*.csv")

  for (i in 1:length(otulist)) {
    name <- substr(otulist[i],1,nchar(otulist[i])-4)
    OTU <-read.csv(file.path(DATA, "amplicon/fungi/by_ecosystem", otulist[i]), header = T, row.names = 1)

    mycol= c("#BFBFBF","#000080","#00BFFF","#EEAEEE","#CD3333","#DB7093","#CD853F","#6959CD","#FFC0CB","#9ACD32","#FFD700","#FFDAB9","#2FFAB9")

    OTU_Tax = merge(OTU,TAX,by="row.names")
    colnames(OTU_Tax)[98]
    ClassRSW=aggregate(x=OTU_Tax[,2:97],by=list(OTU_Tax$class),FUN=sum)
    Class=ClassRSW[,2:97]
    row.names(Class)=ClassRSW[,1]
    ClassRSW[,1]

    Class<-Class/colSums(Class)
    abundance=0.01
    idx = order(rowMeans(Class), decreasing = T)
    Class = Class[idx,]
    idx = rowMeans(Class) > abundance
    idx
    Class1 = Class[idx,]
    Class1['Other',] <- 1-colSums(Class1)
    Class1$Class<-factor(row.names(Class1),levels=rev(rownames(Class1)))

    Class0.01<-melt(Class1,id='Class')
    names(Class0.01)[2]<-'sample'

    Class_15 <- ddply(Class0.01,'sample',transform,percent = value/sum(value)*100)
    data = merge(Class_15,group,by = "sample")

    mean <- aggregate(data$percent, by=list(data$Class, data$sample), FUN=mean)
    colnames(mean) = c("Class", "Sample",  "Mean")
    group2 = group
    mean = merge(mean,group2,by.x = "Sample",by.y = "sample")

    mean[mean == "unclassified"] <- "Other"
    mean[mean == "unidentified"] <- "Other"
    mean$Mean <- mean$Mean/4

    mean$Class <- as.character(mean$Class)
    mean$Class[which(mean$Class=="Other")] <- "ZOther"
    rowno <- as.character(unique(mean$Class))

    rowno <- rowno[!grepl("Other", rowno)]
    rowno <- rowno[order(substring(rowno, 1, 1))]
    .s6_lvl <- rev(c(rowno, "ZOther"))
    mean$Class<-factor(unique(mean$Class),levels=.s6_lvl,labels=sub("^Z","",.s6_lvl))
    mean$group = factor(mean$group, levels = DENSITY_LEVELS, labels = DENSITY_LABELS)
    mean$group1 = factor(mean$group1, levels = TREAT_LEVELS, labels = TREAT_LABELS)

    pALL<-ggplot(mean,aes(x=group1,y=Mean,fill=Class))+
      geom_bar(stat="identity",position="stack",width = 0.8)+
      facet_wrap(~group, scales = 'free_x', ncol = 2) +

      ggtitle("")+
      labs(x="",y="Relative abundance (%)")+
      scale_fill_manual(values = mycol)+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"),
            panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.length=unit(0.1,'cm'),
            legend.text = element_text(size = 8, color="black"),
            legend.key = element_rect( fill = "lemonchiffon"),
            legend.key.size = unit(8, "pt"),
            legend.position = "none",
            strip.text = element_text(size = 6, color="black"),
            axis.title = element_text(size = 8,color="black"),

            axis.text.x = element_text(size = 6, color="black", angle = 45,
                                       hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 6,color="black"),

            axis.title.y = element_text(colour="black", size = 8),
            axis.title.x = element_text(size = 8, color="black"),
      )+
      guides(fill=guide_legend(title=NULL))

    pALL
    outfilename <- paste0(name,"_stacked_bar.pdf")

    ggsave(outfilename, pALL, width = 4.5, height = 8.0, path = OUT, units = "cm")
    if (!exists(".figS6_leg_ITS")) {

      p_leg <- pALL +
        theme(legend.key.size = unit(8, "pt")) +
        guides(fill = guide_legend(title = NULL, ncol = 2))
      .save_legend(p_leg, "FigS6_ITS_stack_legend", OUT, orient = "right")
      .figS6_leg_ITS <- TRUE
    }
  }

  mycol= c("#BFBFBF","#9ACD32","#00BFFF","#9B30FF","#000080","#6959CD","#CD853F","#CD3333","#DB7093","#EEAEEE","#FFC0CB","#FFD700","#FFDAB9")
}
