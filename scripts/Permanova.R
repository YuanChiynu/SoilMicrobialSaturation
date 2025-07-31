{
  setwd("path/to/your/file")
  design <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  
  setwd("path/to/your/file")
  design <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  library(vegan)
  library(dplyr)
  
  resultslis <- c()
  for (otuname in otulist) {
    otu <- read.csv(otuname,row.names = 1)
    group <- design[design$sample %in% colnames(otu),]
    otu<-(otu/unique(colSums(otu)))^0.5
    otu <- as.data.frame(t(otu))
    results1<- adonis2(formula = otu ~ group+group1,data = group, permutations = 999, method = "bray")
    write.csv(results1,file = paste0(otuname,"_permanova.csv"))
    print(otuname)
    print(results1)
    pvalue <- results1[c(1,2),5,drop=F]
    pvalue <- as.data.frame(t(pvalue))
    resultcatcher <- as.data.frame(c(otuname,pvalue))
    colnames(resultcatcher) <- c("OTU","group", "group.group1")
    resultslis <- rbind(resultslis,resultcatcher)
  }
}

{
  setwd("path/to/your/file")
  otu<-read.csv("Verify_rare_16S.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  setwd("path/to/your/file")
  otu<-read.csv("Verify_rare_ITS.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  otu<-(otu/unique(colSums(otu)))^0.5
  otu <- as.data.frame(t(otu))
  colnum <- ncol(otu)
  otu$sample <- rownames(otu)
  group_fin <- group[group$sample %in% otu$sample,]
  group_fin <- group_fin[ group_fin$treat %in% c("10%-D","40%-D","70%-D","100%-D","100%"),]
  otu <- otu[otu$sample %in% group_fin$sample,]
  otu <- otu[,1:colnum]
  results1<- adonis2(formula = otu ~ treat,data = group_fin, permutations = 999, method = "bray")
  results1
}
