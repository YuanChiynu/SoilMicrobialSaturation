{
  library(dplyr)
  setwd("path/to/your/file")
  
  data <- read.delim("otutab_raw_16s.txt", header=TRUE, sep="path/to/your/file",row.names=1)
  OTU <- data
  
  
  n <- 9 
  cols <- names(data) 
  long_cols <- cols[nchar(cols) > n] 
  new_cols <- unique(substr(long_cols, 3, nchar(long_cols)))
  
  
  olddata <- data[-c(grep("F1|F2",colnames(data)))]
  f1data <- data[c(grep("F1",colnames(data)))]
  f2data <- data[c(grep("F2",colnames(data)))]
  
  uniqdata <- olddata[c(grep(new_cols[1],colnames(olddata)))] 
  for (sampcot in 2:length(new_cols)) {
    tempdata <- olddata[c(grep(new_cols[sampcot],colnames(olddata)))]
    uniqdata <- cbind(uniqdata,tempdata)
  }
  removlist <- as.vector(colnames(uniqdata))
  oldsuc <-olddata[ , !colnames(olddata) %in% removlist]
  
  supllementry <- oldsuc[,c(1,2)]
  for (sampcot in 1:length(new_cols)){
    f1temp <- f1data[c(grep(new_cols[sampcot],colnames(f1data)))]
    f2temp <- f2data[c(grep(new_cols[sampcot],colnames(f2data)))]
    oldtemp <- uniqdata[c(grep(new_cols[sampcot],colnames(uniqdata)))]
    tempD <- cbind(f1temp,f2temp,oldtemp)
    if (colSums(tempD)[1]>50000) {
      supllementry <- cbind(supllementry,f1temp)
    }else if (colSums(tempD)[2]>50000) {
      supllementry <- cbind(supllementry,f2temp)
    }else if ((colSums(tempD)[1]+colSums(tempD)[2])>50000) {
      samplename <- substr(colnames(tempD)[1],3,nchar(colnames(tempD)[1]))
      tempD$sum_2 <- rowSums(tempD[,c(1,2)])
      tempD <- tempD[,c(1,ncol(tempD))]
      colnames(tempD)[ncol(tempD)] <- samplename
      supllementry <- cbind(supllementry,tempD[2])
    }else{
      samplename <- substr(colnames(tempD)[1],3,nchar(colnames(tempD)[1]))
      tempD$sum_3 <- rowSums(tempD[,c(1,2,3)])
      tempD <- tempD[,c(1,ncol(tempD))]
      colnames(tempD)[ncol(tempD)] <- samplename
      supllementry <- cbind(supllementry,tempD[2])
    }
  }
  supllementry <- supllementry[,-c(1,2)]
  
  otufilter <- cbind(oldsuc,supllementry)
  write.csv(otufilter,file='OTU_raw_16S_all.csv')
  
  
  otutab<- otufilter
  otutab <- t(otutab)
  otutab <- as.data.frame(otutab)
  
  rownames(otutab)
  colnames(otutab)
  colnum<-ncol(otutab)
  otutab$sum <- rowSums(otutab[, 1:colnum])
  otutab[,colnum+1]
  otutab1<-otutab[order(otutab[,colnum+1]),] 
  head(otutab1[,colnum+1])
  otutab <- otutab1[,-colnum+1]
  
  library(permute)
  library(vegan)
  otu<-otutab
  raremax1 <- min(rowSums(otu))
  newOTU <- rrarefy(otu, raremax1)
  newOTU1 = t(newOTU)
  newOTU1 <- as.data.frame(newOTU1)
  write.csv(newOTU1,file='OTU_rare_16S_all.csv')
  
}

{
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  
  otu <- t(otu)
  otu <- as.data.frame(otu)
  colnum <- ncol(otu)
  colnumP1 <- colnum+1
  otu$name <- rownames(otu)
  
  class(otu[1,(colnumP1)])
  
  group <- read.csv('分组.csv', header=TRUE)
  sample_otu <- merge(otu, group, by.x = 'name',by.y = 'names', all.x = TRUE)
  rownames(sample_otu) <-sample_otu[,1]
  sample_otu <- sample_otu[,-1]
  colnumP2 <- ncol(sample_otu)
  
  records<-c()
  for (ecosys in c("G","F","C")) {
    otuputName <-"16S"
    sub_eco <- subset(sample_otu,ecosystem==ecosys)
    out_filenamelv1<-sapply(otuputName, function(x)paste(x,ecosys,sep = "_")) 
    out_filename<-sapply(out_filenamelv1, function(x)paste(x,".csv",sep = ""))
    otudata <- sub_eco[-c(colnumP1:colnumP2)]
    otudata<-t(otudata)
    setwd("path/to/your/file")
    write.csv(otudata, out_filename, row.names = T)
    printrecords<-c(out_filename,nrow(sub_eco))
    records<- rbind(records,printrecords)
  }
  print(records)
  records<-c()
  for (ecosys in c("G","F","C")) {
    otuputName <-"16S"
    sub_eco <- subset(sample_otu,ecosystem==ecosys)
    out_filenamelv1<-sapply(otuputName, function(x)paste(x,ecosys,sep = "_")) 
    for (saturation in c("100%re","70%re","40%re","10%re")) {
      out_filename<-out_filenamelv1
      sub_eco_satru <- subset(sub_eco,group==saturation)
      out_filenamelv2<-sapply(out_filename, function(x)paste(x,saturation,sep = "_"))
      out_filename<-sapply(out_filenamelv2, function(x)paste(x,".csv",sep = ""))
      otudata <- sub_eco_satru[-c(colnumP1:colnumP2)]
      otudata<-t(otudata)
      setwd("path/to/your/file")
      write.csv(otudata, out_filename, row.names = T)
      printrecords<-c(out_filename,nrow(sub_eco_satru))
      records<- rbind(records,printrecords)
      
    }
  }
  print(records)
  
  records<-c()
  for (ecosys in c("C","S","T")) {
    otuputName <-"16S"
    sub_eco <- subset(sample_otu,ecosystem==ecosys)
    out_filenamelv1<-sapply(otuputName, function(x)paste(x,ecosys,sep = "_")) 
    for (saturation in c("100%re","70%re","40%re","10%re")) {
      out_filename<-out_filenamelv1
      sub_eco_satru <- subset(sub_eco,group==saturation)
      out_filenamelv2<-sapply(out_filename, function(x)paste(x,saturation,sep = "_")) 
      for (treatment in c("C","pH","W","N","CK","T" )) {
        out_filename<-out_filenamelv2
        sub_eco_satru_treat <- subset(sub_eco_satru,group1==treatment)
        out_filename<-sapply(out_filename, function(x)paste(x,treatment,sep = "_"))
        out_filename<-sapply(out_filename, function(x)paste(x,".csv",sep = ""))
        otudata <- sub_eco_satru_treat[-c(40118:40121)]
        otudata<-t(otudata)
        write.csv(otudata, out_filename, row.names = T)
        printrecords<-c(out_filename,nrow(sub_eco_satru_treat))
        records<- rbind(records,printrecords)
      }
    }
  }
  print(records)

  records<-c()
  for (treat in c("C","pH","W","N","CK","T")) {
    otuputName <-"ITS"
    sub_eco <- subset(sample_otu,group1==treat)
    out_filenamelv1<-sapply(otuputName, function(x)paste(x,treat,sep = "_")) 
    out_filename<-sapply(out_filenamelv1, function(x)paste(x,".csv",sep = ""))
    otudata <- sub_eco[-c(colnumP1:colnumP2)]
    otudata<-t(otudata)
    setwd("path/to/your/file")
    write.csv(otudata, out_filename, row.names = T)
    printrecords<-c(out_filename,nrow(sub_eco))
    records<- rbind(records,printrecords)
  }
  
  print(records)
  
}

{
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(agricolae)

  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  for (i in c(1:length(otulist))){
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-4)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    otu <- t(otu)
    sumrow<-rowSums(otu)[1]
    otu <- (otu/sumrow)^0.5
    otu <- as.data.frame(otu)
    
    otudata <- otu
    nmds1 <- metaMDS(otudata, autotransform =FALSE)
    nmds1.stress <- nmds1$stress
    nmds1.point <- data.frame(nmds1$point)
    nmds1.species <- data.frame(nmds1$species)
    
    
    sample_site <- nmds1.point[1:2]
    sample_site$names <- rownames(sample_site)
    names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
    sample_site <- merge(sample_site, group, by.x = 'names',  by.y = 'sample',all.x = TRUE)
    
    sample_site$group <- factor(sample_site$group,levels=c("100%re","70%re","40%re","10%re"))

    my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
    p <- ggplot(sample_site, aes(x=NMDS1, y=NMDS2, group = group)) +
      geom_point(aes(color = group, shape = group1), size = 4, alpha =2) + 
      scale_shape_manual(values = c(15,3,16,17,18,4)) + 
      scale_color_manual(values = c("
      stat_ellipse(data=sample_site,
                   geom = "polygon",level=0.95,
                   linetype = 1,size=0.8,
                   aes(fill=group),
                   alpha=0.1,
                   show.legend = T)+
      scale_fill_manual(values = c("
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) + 
      labs(x = 'NMDS 1', y = 'NMDS 2', title = paste('Stress =', round(nmds1$stress,3))) +
      theme(plot.title = element_text(hjust = 0.5,size = 20,face = "plain",vjust = -7))+  
      theme(axis.text.x = element_text(size = 20,face = "plain"))+
      theme(axis.text.y = element_text(size = 20,face = "plain"))+
      theme(axis.title.x = element_text(size = 20, face = "plain"))+
      theme(axis.title.y = element_text(size = 20, face = "plain"))+
      theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = 1))+
      theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = 1))+
      theme(legend.position = "right")+
      my_theme
    p 
    
    out_filename<-sapply(savename, function(x)paste(x,"NMDS.pdf",sep = "")) 
    ggsave(p, file=out_filename, width = 6, height = 6)
  }
  
  
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  for (i in c(1:length(otulist))){
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-7)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    otu <- t(otu)
    sumrow<-ncol(otu)
    otu <- (otu/sumrow)^0.5
    otu <- as.data.frame(otu)
    
    otudata <- otu
    nmds1 <- metaMDS(otudata, autotransform =FALSE)
    nmds1.stress <- nmds1$stress
    nmds1.point <- data.frame(nmds1$point)
    nmds1.species <- data.frame(nmds1$species)
    
    
    sample_site <- nmds1.point[1:2]
    sample_site$names <- rownames(sample_site)
    names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
    sample_site <- merge(sample_site, group, by = 'names', all.x = TRUE)
    
    sample_site$group1 <- factor(sample_site$group1,levels=c("C","CK","N","pH","T","W"))
    sample_site$group1
    
    p <- ggplot(sample_site, aes(NMDS1, NMDS2, group = group)) +
      geom_point(aes(shape = group1), size = 6, alpha =2) + 
      scale_shape_manual(values = c(15,17,3,18,16,8)) + 
      scale_color_manual(values = c("
      theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) + 
      theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank()) + 
      labs(x = 'NMDS 1', y = 'NMDS 2', title = paste('Stress =', round(nmds1$stress,3))) +
      theme(plot.title = element_text(hjust = 0.5,size = 20,face = "plain",vjust = -7))+  
      theme(axis.text.x = element_text(size = 20,face = "plain"))+
      theme(axis.text.y = element_text(size = 20,face = "plain"))+
      theme(axis.title.x = element_text(size = 20, face = "plain"))+
      theme(axis.title.y = element_text(size = 20, face = "plain"))+
      theme(legend.text=element_text(size=20, face = "plain"))+
      theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
      theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = 1))+
      theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = 1))
    p 
    out_filename<-sapply(savename, function(x)paste(x,".pdf",sep = "")) 
    ggsave(p, file=out_filename, width = 6.5, height = 5)
  }
  
}

{
  library(permute)
  library(vegan)
  library(reshape)
  library(ggpubr)
  library(ggplot2)
  library(agricolae)
  library(data.table)

  
  setwd("path/to/your/file")
  ado1<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  
  setwd("path/to/your/file")
  ado1<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  
  otusum<-unname(colSums(ado1)[1])
  
  ado=(t(ado1)/otusum)^0.5
  prok.dist <- vegdist(ado)
  
  prok.dist <- as.matrix(prok.dist)
  prok.dist <- as.data.frame(prok.dist)
  prok.dist$X1 <- rownames(prok.dist)
  prok.matr<- reshape2::melt(prok.dist,id.vars = "X1", variable.name = "X2", value.name = "score") 
  
  colnames(prok.matr)[2] <-"X2"
  colnames(prok.matr)[3] <-"similarity" 
  prok.matr$similarity <- 1-prok.matr$similarity
  reselect <-prok.matr
  
  group1<-read.csv("相似性柱状图分组1.csv", header=TRUE)
  group2<-read.csv("相似性柱状图分组2.csv", header=TRUE)
  simi_dis <- merge(reselect, group1, by.x = 'X1',by.y = 'names', all.x = TRUE)
  simi_dis <- merge(simi_dis, group2, by.x = 'X2',by.y = 'names_2', all.x = TRUE)

  fin_dis_zhu<-subset(simi_dis,ecosystem_2==ecosystem)
  fin_dis_zhu<-subset(fin_dis_zhu,group_2==group)
  fin_dis_zhu<-subset(fin_dis_zhu,group1_2!=group1)
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="ecosystem")] <-"site" 
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="group")] <-"Group" 
  
  results_summary <- fin_dis_zhu %>%
    group_by(site,Group) %>%  
    summarise(
      mean_similarity = mean(similarity, na.rm = TRUE),  
      sd_similarity   = sd(similarity, na.rm = TRUE)      
    ) %>%
    ungroup()  
  write.csv(fin_dis_zhu,file='fin_dis_zhu.csv')
  
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="ecosystem")] <-"site" 
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="group")] <-"Group" 
  phi<-read.csv(file='fin_dis_zhu.csv',header = T,row.names = 1)

  reresult <- c()
  for (ecosys in unique(phi$site)){
    sub_phi <- subset(phi,site==ecosys)
    mod1 = aov( similarity~Group, data= sub_phi)
    summary(mod1)
    re = LSD.test(mod1,"Group",alpha = 0.05)
    re1 = re$groups
    re1$groupfactor <- paste0(ecosys,rownames(re1))
    print(re1)
    reresult <- rbind(reresult,re1)
  }
  reresult
  compare_means(similarity~Group, data=phi, group.by = "site",method = "anova")
  phi$site = factor( phi$site , levels = c("C","F","G"))
  phi$Group = factor(phi$Group, levels = c("10%re","40%re","70%re","100%re"))
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  p_16s <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",
                  palette = c("
    ylab("Similarity")+
    xlab("Ecosystem")+
    ylim(0,0.7)+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.65, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  p_16s
  setwd("path/to/your/file")
  ggsave('similarity_16s.pdf', p_16s, width = 6, height = 6)
  
  its_p <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",
                  palette = c("
    ylab("Similarity")+
    xlab("Ecosystem")+
    ylim(0,0.8)+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.70, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  its_p
  setwd("path/to/your/file")
  ggsave('similarity_its.pdf', its_p, width = 9, height = 9)
  
  c("
}

{
  library(permute)
  library(vegan)
  library(reshape)
  library(ggpubr)
  library(ggplot2)
  library(agricolae)
  
  setwd("path/to/your/file")
  
  setwd("path/to/your/file")

  phi<-read.csv("fin_dis_zhu.csv", header=TRUE)
  
  reresult <- c()
  for (ecosys in unique(phi$site)){
    sub_phi <- subset(phi,site==ecosys)
    mod1 = aov( differ~Group, data= sub_phi)
    summary(mod1)
    re = LSD.test(mod1,"Group",alpha = 0.05)
    re1 = re$groups
    re1$groupfactor <- paste0(ecosys,rownames(re1))
    print(re1)
    reresult <- rbind(reresult,re1)
  }
  reresult
  compare_means(differ~Group, data=phi, group.by = "site",method = "anova")
  phi$Group = factor(phi$Group, levels = c("10%re","40%re","70%re","100%re"))
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  phi$differ <- round(phi$differ,4)
  
  p_16s <- ggbarplot(phi, x="site", y="differ", add = "mean_se",color = "black",fill = "Group",
                     palette = c("
    ylab("Bray-Curtis")+
    xlab("Ecosystem")+
    ylim(0,0.7)+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.65, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  p_16s
  
  ggsave('Bray-Curtis_16s.pdf', p_16s, width = 9, height = 9)
  
  its_p <- ggbarplot(phi, x="site", y="differ", add = "mean_se",color = "black",fill = "Group",
                     palette = c("
    ylab("Bray-Curtis")+
    xlab("Ecosystem")+
    ylim(0,0.781)+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.69, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  its_p
  
  ggsave('Bray-Curtis_its.pdf', its_p, width = 9, height = 9)
  
  c("
}

{
  library(Hmisc)
  library(psych)
  library(igraph)
 
  
    
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  
  
  
  otulist = list.files(pattern="*.csv")
  
  resultdeprition <- c()
  for (i in c(1:length(otulist))) {
    
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-4)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    otu1 <- otu
    otu1[otu1>0] <- 1
    oturownum <- ncol(otu1)/2
    otu <- otu[which(rowSums(otu1) >= oturownum),]
    
    
    rowsum <- as.data.frame(rowSums(otu))
    sum <- colSums(rowsum)
    otu$RA <- rowsum/sum
    otu <- subset(otu[which(otu$RA>0.0001),], select = -RA)
    
    
    otu <- t(otu)
    otu <- scale(otu)
    
    rcorr_otu <- Hmisc::rcorr(as.matrix(otu), type = 'spearman')
    
    p <- rcorr_otu$P
    p <- p.adjust(p, method = 'BH')
    
    
    r <- rcorr_otu$r
    r[which(r < 0)]
    r[abs(r) < 0.8] <- 0
    
    
    p[p>=0.05] <- -1
    p[p<0.05 & p>=0] <- 1
    p[p==-1] <- 0
    
    
    z <- r * p
    z[which(z < 0)]
    
    
    
    g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected', diag = FALSE)
    g
    
    g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
    names(degree(g)[degree(g) == 0])
    
    E(g)$correlation <- E(g)$weight
    E(g)$weight <- abs(E(g)$weight)
    
    
    taxonomy <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
    taxonomy <- taxonomy[as.character(V(g)$name),]
    V(g)$kingdom <- taxonomy$kingdom
    V(g)$phylum <- taxonomy$phylum
    V(g)$class <- taxonomy$class
    V(g)$order <- taxonomy$order
    V(g)$family <- taxonomy$family
    V(g)$genus <- taxonomy$genus
    V(g)$species <- taxonomy$species
    
    
    
    edge <- data.frame(as_edgelist(g))    
    
    df <- as.data.frame(E(g)$correlation)
    df[df>0] <- 1
    df[df<0] <- -1
    colnames(df) <- c('cor')
    
    edge_list <- data.frame(
      source = edge[[1]],
      target = edge[[2]],
      weight = E(g)$weight,
      correlation = E(g)$correlation,
      cor = df
    )
    
    head(edge_list)
    out_filename<-sapply(savename, function(x)paste(x,"edge_screen0.0001_r0.8p0.05.csv",sep = "")) 
    write.table(edge_list,out_filename, sep = ',', row.names = FALSE, quote = FALSE)
    
    deprition <- otulist[i]
    posedgenum <- nrow(subset(edge_list,cor==1))
    negedgenum <- nrow(edge_list)-posedgenum
    deprition <-cbind(deprition,posedgenum,negedgenum)
    
    
    node <- data.frame(
      id = names(V(g)),
      kingdom =V(g)$kingdom,
      phylum =V(g)$phylum,
      class =V(g)$class,
      order =V(g)$order,
      family =V(g)$family,
      genus =V(g)$genus,
      species =V(g)$species
    )
    out_filename<-sapply(savename, function(x)paste(x,"node_screen0.0001_r0.8p0.05.csv",sep = "")) 
    write.table(node, out_filename, sep = ',', row.names = FALSE, quote = FALSE)
    
    
    nodenum <- nrow(node)
    deprition <-cbind(deprition,nodenum)
    
    
    num.edges = length(E(g)) 
    num.edges
    num.vertices = length(V(g))
    num.vertices
    connectance = edge_density(g,loops=FALSE)
    connectance
    average.degree = mean(igraph::degree(g))
    average.degree
    average.path.length = average.path.length(g) 
    average.path.length
    diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
    diameter
    edge.connectivity = edge_connectivity(g)
    edge.connectivity
    clustering.coefficient = transitivity(g) 
    clustering.coefficient
    no.clusters = no.clusters(g)
    no.clusters
    centralization.betweenness = centralization.betweenness(g)$centralization 
    centralization.betweenness
    centralization.degree = centralization.degree(g)$centralization
    centralization.degree
    
    network_property <-c(num.edges=num.edges,num.vertices=num.vertices,connectance=connectance,
                         average.degree=average.degree,average.path.length=average.path.length,
                         diameter=diameter,edge.connectivity=edge.connectivity,clustering.coefficient=clustering.coefficient,
                         no.clusters=no.clusters,centralization.betweenness=centralization.betweenness,
                         centralization.degree=centralization.degree)
    network_property <- data.frame(network_property)
    out_filename<-sapply(savename, function(x)paste(x,"network_property.csv",sep = "")) 
    write.csv(network_property, out_filename, row.names = T)
    resultdeprition <- rbind(resultdeprition,deprition)
  }
  
  resultdeprition
  write.csv(resultdeprition, "边、点信息.csv", row.names = T)
}

{
  library(ggvenn)
  library(ggplot2)
  library(dplyr)
  library(VennDiagram)
  
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  
  otulist = list.files(pattern="*.csv")
  for (i in c(1:length(otulist))){
    name <- substr(otulist[i],1,nchar(otulist[i])-31)
    data <- read.csv(otulist[i], header=TRUE)
    
    count <- 0
    for (i in 1:nrow(data)) {
      sorted_row <- data[i, ] %>% arrange(source,target) 
      if (sorted_row$target < sorted_row$source) { 
        data[i, c("source", "target")] <- sorted_row[c("target", "source")]
        count <- count+1
      }
    }
    print(paste0(name,"改变次数",count))
    outfile <- paste0(name,"_order_edge.csv")
    write.csv(data, outfile, row.names = T)
  }
  S10 <- read.csv("16S_G_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("16S_G_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("16S_G_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("16S_G_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                      ),
               fill = c("
               print.mode="percent",
               filename = NULL)
  pdf("16S_Grass.pdf")
  grid.draw(vennplot)
  dev.off()
  S10 <- read.csv("16S_F_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("16S_F_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("16S_F_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("16S_F_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                     ),
                    fill = c("
                    print.mode="percent",
                    filename = NULL)
  pdf("16S_Forest.pdf")
  grid.draw(vennplot)
  dev.off()
  S10 <- read.csv("16S_C_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("16S_C_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("16S_C_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("16S_C_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                     ),
                    fill = c("
                    print.mode="percent",
                    filename = NULL)
  pdf("16S_Corpland.pdf")
  grid.draw(vennplot)
  dev.off()
  
  
  S10 <- read.csv("ITS_C_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("ITS_C_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("ITS_C_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("ITS_C_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                      ),
                      fill = c("
                      print.mode="percent",
                      filename = NULL)
  pdf("ITS_Grass.pdf")
  grid.draw(vennplot)
  dev.off()
  S10 <- read.csv("ITS_S_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("ITS_S_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("ITS_S_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("ITS_S_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                      ),
                      fill = c("
                      print.mode="percent",
                      filename = NULL)
  pdf("ITS_Forest.pdf")
  grid.draw(vennplot)
  dev.off()
  S10 <- read.csv("ITS_T_10%re_order_edge.csv", header=TRUE)
  S10 <- paste0(S10$source,S10$target)
  
  S40 <- read.csv("ITS_T_40%re_order_edge.csv", header=TRUE)
  S40 <- paste0(S40$source,S40$target)
  
  S70 <- read.csv("ITS_T_70%re_order_edge.csv", header=TRUE)
  S70 <- paste0(S70$source,S70$target)
  
  S100 <- read.csv("ITS_T_100%re_order_edge.csv", header=TRUE)
  S100 <- paste0(S100$source,S100$target)
  venndata <- list(S10=S10,S40=S40,S70=S70,S100=S100)
  vennplot<- venn.diagram(x=list(S10=S10,
                      S40=S40,
                      S70=S70,
                      S100=S100
                      ),
                      fill = c("
                      print.mode="percent",
                      filename = NULL)
  
  pdf("ITS_Corpland.pdf")
  grid.draw(vennplot)
  dev.off()
}

{
  library(reshape2)
  library(plyr)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)
  
  setwd("path/to/your/file")
  group = read.csv("分组.csv", header = T)
  TAX  <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  TAX <- TAX %>% 
    mutate_all(~ ifelse(. != "Other", substr(., 4, nchar(.)), .))
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")

  for (i in 1:length(otulist)) {
    name <- substr(otulist[i],1,nchar(otulist[i])-4)
    OTU <-read.csv(otulist[i], header = T, row.names = 1)
    mycol= c("
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
    for (k in 2:97) {
      mean_df[colnames(other_rows)[k]] <- mean_values[k]
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
    mean$Class<-factor(unique(mean$Class),levels=rev(rowno))
    mean$group = factor(mean$group,levels = c("10%re","40%re","70%re","100%re"))
    mean$group1 = factor(mean$group1,levels = c("CK","C","N","pH","T","W"))
    
    pALL<-ggplot(mean,aes(x=group1,y=Mean,fill=Class))+
      geom_bar(stat="identity",position="stack",width = 0.8)+ 
      facet_wrap(~group, scales = 'free_x', ncol = 2) +
      ggtitle("")+
      labs(x="",y="Relative abundance(%)")+
      scale_fill_manual(values = mycol)+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"), 
            panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.length=unit(0.05,'cm'),
            legend.text=element_text(size=8, color="black"),
            legend.key = element_rect( fill = "lemonchiffon"),
            legend.key.size = unit(8, "pt"),
            legend.position = "right",
            axis.title=element_text(size=8,color="black"),
            axis.text.x = element_text(size = 8,color="black"),
            axis.text.y = element_text(size = 8,color="black"),
            axis.title.y = element_text(colour="black", size = 8),
            axis.title.x  = element_text(size = 8, color="black"),
      )+
      guides(fill=guide_legend(title=NULL))
    
    
    pALL
    outfilename <- paste0(name,"_堆叠图.pdf")
    ggsave(outfilename, pALL, width = 10, height = 10,units = "cm")
  }
  
  setwd("path/to/your/file")
  group = read.csv("分组.csv", header = T)
  TAX  <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  TAX <- TAX %>% 
    mutate_all(~ ifelse(. != "Other", substr(., 4, nchar(.)), .))
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  
  for (i in 1:length(otulist)) {
    name <- substr(otulist[i],1,nchar(otulist[i])-4)
    OTU <-read.csv(otulist[i], header = T, row.names = 1)
    mycol= c("
    
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
    mean[mean == "c__unidentified"] <- "Other"
    mean$Mean <- mean$Mean/4
    
    dev.new()
    mean$Class <- as.character(mean$Class)
    mean$Class[which(mean$Class=="Other")] <- "ZOther"
    rowno <- as.character(unique(mean$Class))
    rowno <- rowno[order(substring(rowno, 1, 1))][!grepl("Other", rowno)]
    rowno <- c(rowno,"Other")
    mean$Class<-factor(unique(mean$Class),levels=rev(rowno))
    mean$group = factor(mean$group,levels = c("10%re","40%re","70%re","100%re"))
    mean$group1 = factor(mean$group1,levels = c("CK","C","N","pH","T","W"))
    

    pALL<-ggplot(mean,aes(x=group,y=Mean,fill=Class))+
      geom_bar(stat="identity",position="stack",width = 0.8)+ 
      facet_wrap(~group1, scales = 'free_x', ncol = 3) +
      ggtitle("")+
      labs(x="",y="Relative abundance(%)")+
      scale_fill_manual(values = mycol)+
      theme_bw()+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"), 
            panel.grid = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.length=unit(0.2,'cm'),
            legend.text=element_text(size=19, color="black"),
            legend.key = element_rect( fill = "lemonchiffon"),
            legend.key.size = unit(19, "pt"),
            legend.position = "right",
            axis.title=element_text(size=19,color="black"),
            axis.text.x = element_text(size = 14,color="black"),
            axis.text.y = element_text(size = 14,color="black"),
            axis.title.y = element_text(colour="black", size = 16),
            axis.title.x  = element_text(size = 16, color="black"),
      )+
      guides(fill=guide_legend(title=NULL))
    
    
    pALL
    outfilename <- paste0(name,"_堆叠图.pdf")
    ggsave(outfilename, pALL, width = 10, height = 10)
    dev.off()
  }
  
  mycol= c("
}

{
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(d3Network)
  library(grid)
  

  setwd("path/to/your/file")
  setwd("path/to/your/file")
  S10 <- read.csv("ITS_T_10%re.csv", header=TRUE,row.names = 1)
  S10 <- S10/5.8318
  
  S40 <- read.csv("ITS_T_40%re.csv", header=TRUE,row.names = 1)
  S40 <- S40/5.8318
  
  S70 <- read.csv("ITS_T_70%re.csv", header=TRUE,row.names = 1)
  S70 <- S70/5.8318
  
  S100 <- read.csv("ITS_T_100%re.csv", header=TRUE,row.names = 1)
  S100 <- S100/5.8318
  
  S100testresultfin <-c()
  S100names <- rownames(S100)
  S100mean <- rowMeans(S100)
  for (i in 1:nrow(S100)) { 
    otuname <- S100names[i]
    S100RAmean <- S100mean[i]
    result <- t.test(S100[i,],alternative ="less", mu = 1,conf.level = 0.95)
    S100rare_p <- result$p.value 
    S100rare_p <-formatC(S100rare_p, digits = 5, format = "f")
    temp <- cbind(otuname,S100RAmean,S100rare_p)
    S100testresultfin <- rbind(S100testresultfin,temp)
  }
  S100testresultfin <- data.frame(S100testresultfin)
  S100testresultfin$S100type <- ifelse(S100testresultfin$S100RAmean == 0, "NoneOccurence",
                                ifelse(S100testresultfin$S100RAmean < 0.1 & S100testresultfin$S100rare_p > 0.05, "RareNosig",
                                      ifelse(S100testresultfin$S100RAmean < 0.1 & S100testresultfin$S100rare_p <= 0.05, "RareSig", "Abundunt")))
  
  S70testresultfin <-c()
  S70names <- rownames(S70)
  S70mean <- rowMeans(S70)
  for (i in 1:nrow(S70)) { 
    otuname <- S70names[i]
    S70RAmean <- S70mean[i]
    result <- t.test(S70[i,],alternative ="less", mu = 1,conf.level = 0.95)
    S70rare_p <- result$p.value 
    S70rare_p <-formatC(S70rare_p, digits = 5, format = "f")
    temp <- cbind(otuname,S70RAmean,S70rare_p)
    S70testresultfin <- rbind(S70testresultfin,temp)
  }
  S70testresultfin <- data.frame(S70testresultfin)
  S70testresultfin$S70type <- ifelse(S70testresultfin$S70RAmean == 0, "NoneOccurence",
                                   ifelse(S70testresultfin$S70RAmean < 0.1 & S70testresultfin$S70rare_p > 0.05, "RareNosig",
                                          ifelse(S70testresultfin$S70RAmean < 0.1 & S70testresultfin$S70rare_p <= 0.05, "RareSig", "Abundunt")))
  
  S40testresultfin <-c()
  S40names <- rownames(S40)
  S40mean <- rowMeans(S40)
  for (i in 1:nrow(S40)) { 
    otuname <- S40names[i]
    S40RAmean <- S40mean[i]
    result <- t.test(S40[i,],alternative ="less", mu = 1,conf.level = 0.95)
    S40rare_p <- result$p.value 
    S40rare_p <-formatC(S40rare_p, digits = 5, format = "f")
    temp <- cbind(otuname,S40RAmean,S40rare_p)
    S40testresultfin <- rbind(S40testresultfin,temp)
  }
  S40testresultfin <- data.frame(S40testresultfin)
  S40testresultfin$S40type <- ifelse(S40testresultfin$S40RAmean == 0, "NoneOccurence",
                                   ifelse(S40testresultfin$S40RAmean < 0.1 & S40testresultfin$S40rare_p > 0.05, "RareNosig",
                                          ifelse(S40testresultfin$S40RAmean < 0.1 & S40testresultfin$S40rare_p <= 0.05, "RareSig", "Abundunt")))
  S10testresultfin <-c()
  S10names <- rownames(S10)
  S10mean <- rowMeans(S10)
  for (i in 1:nrow(S10)) { 
    otuname <- S10names[i]
    S10RAmean <- S10mean[i]
    result <- t.test(S10[i,],alternative ="less", mu = 1,conf.level = 0.95)
    S10rare_p <- result$p.value 
    S10rare_p <-formatC(S10rare_p, digits = 5, format = "f")
    temp <- cbind(otuname,S10RAmean,S10rare_p)
    S10testresultfin <- rbind(S10testresultfin,temp)
  }
  S10testresultfin <- data.frame(S10testresultfin)
  S10testresultfin$S10type <- ifelse(S10testresultfin$S10RAmean == 0, "NoneOccurence",
                                   ifelse(S10testresultfin$S10RAmean < 0.1 & S10testresultfin$S10rare_p > 0.05, "RareNosig",
                                          ifelse(S10testresultfin$S10RAmean < 0.1 & S10testresultfin$S10rare_p <= 0.05, "RareSig", "Abundunt")))
  sanjiDate <- merge(S10testresultfin, S40testresultfin, by = "row.names", all = TRUE)
  rownames(sanjiDate) <- sanjiDate[,1]
  sanjiDate <- sanjiDate[,c(5,9)]
  sanjiDate <- merge(sanjiDate, S70testresultfin, by = "row.names", all = TRUE)
  rownames(sanjiDate) <- sanjiDate[,1]
  sanjiDate <- sanjiDate[,c(2,3,7)]
  sanjiDate <- merge(sanjiDate, S100testresultfin, by = "row.names", all = TRUE)
  rownames(sanjiDate) <- sanjiDate[,1]
  sanjiDate <- sanjiDate[,c(2,3,4,8)]
  
  
  
  sanjiDate$summary <- paste(sanjiDate[,4],sanjiDate[,3],sanjiDate[,2],sanjiDate[,1])
  sanjiDate <- sanjiDate %>% arrange(summary)
  counttemp <- sanjiDate[1,]
  counttemp$count <- 1
  rownum <- 1
  for (cot in 2:nrow(sanjiDate)) {
    if (sanjiDate[cot,5]==counttemp[rownum,5]) {
      counttemp[rownum,6] <- counttemp[rownum,6]+1
    }else{
      rownum <- rownum+1
      counttemp[rownum,c(1:5)] <- sanjiDate[cot,c(1:5)]
      counttemp[rownum,6] <- 1
    }
  }
  counttemp[counttemp == "NoneOccurence"] <- NA
  counttemp <- na.omit(counttemp)
  counttemp <- subset(counttemp,S100type!="Abundunt")
  counttemp <- counttemp %>% arrange(S100type,S70type,S40type,S10type)
  counttemp <-counttemp[c(1:7),]
  
  

  
  p1 <- ggplot(data = counttemp,
         aes(axis1 = S100type, axis2 = S70type,
             y = count)) +
    scale_x_discrete(limits = c("100%", "70%", "40%", "10%")) +
    xlab("Demographic") +
    geom_alluvium(aes(fill = S70type),curve_type="xspline",knot.pos=0.05)+
    scale_fill_manual(values = c(RareSig = "
    geom_stratum(width=0.1) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
    scale_y_continuous(breaks = NULL)+
    theme(panel.grid.major = element_line(colour=NA),
          legend.position = "None"
    )
  
  p2 <- ggplot(data = counttemp,
         aes(axis2 = S70type, axis3 = S40type,
             y = count)) +
    scale_x_discrete(limits = c("100%", "70%", "40%", "10%")) +
    xlab("Demographic") +
    geom_alluvium(aes(fill = S40type),curve_type="xspline",knot.pos=0.05)+
    scale_fill_manual(values = c(RareSig = "
    geom_stratum(width=0.1) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
    scale_y_continuous(breaks = NULL)+
    theme(panel.grid.major = element_line(colour=NA),
          legend.position = "None"
    )
  
  p3 <- ggplot(data = counttemp,
         aes(axis3 = S40type, axis4 = S10type,
             y = count)) +
    scale_x_discrete(limits = c("100%", "70%", "40%", "10%")) +
    xlab("Demographic") +
    geom_alluvium(aes(fill = S10type),curve_type="xspline",knot.pos=0.05)+
    scale_fill_manual(values = c(RareSig = "
    geom_stratum(width=0.1) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
    scale_y_continuous(breaks = NULL)+
    theme(panel.grid.major = element_line(colour=NA),
          legend.position = "None"
    )
  ggsave(p1, file="ITS_Corpland_100-70_sankey.pdf", width = 6, height = 4.5)
  ggsave(p2, file="ITS_Corpland_70-40_sankey.pdf", width = 6, height = 4.5)
  ggsave(p3, file="ITS_Corpland_40-10_sankey.pdf", width = 6, height = 4.5)
  
  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  for (name in c("ITS_C_","ITS_S_","ITS_T_")) {
    outname <- paste0(name,"10%re.csv")
    S10 <- read.csv(outname, header=TRUE,row.names = 1)
    S10 <- S10/5.8318
    
    outname <- paste0(name,"40%re.csv")
    S40 <- read.csv(outname, header=TRUE,row.names = 1)
    S40 <- S40/5.8318
    
    outname <- paste0(name,"70%re.csv")
    S70 <- read.csv(outname, header=TRUE,row.names = 1)
    S70 <- S70/5.8318
    
    outname <- paste0(name,"100%re.csv")
    S100 <- read.csv(outname, header=TRUE,row.names = 1)
    S100 <- S100/5.8318
    
    S100testresultfin <-c()
    S100names <- rownames(S100)
    S100mean <- rowMeans(S100)
    for (i in 1:nrow(S100)) { 
      otuname <- S100names[i]
      S100RAmean <- S100mean[i]
      result <- t.test(S100[i,],alternative ="less", mu = 1,conf.level = 0.95)
      S100rare_p <- result$p.value 
      S100rare_p <-formatC(S100rare_p, digits = 5, format = "f")
      temp <- cbind(otuname,S100RAmean,S100rare_p)
      S100testresultfin <- rbind(S100testresultfin,temp)
    }
    S100testresultfin <- data.frame(S100testresultfin)
    S100testresultfin$S100type <- ifelse(S100testresultfin$S100RAmean == 0, "NoneOccurence",
                                         ifelse(S100testresultfin$S100RAmean < 0.1 & S100testresultfin$S100rare_p > 0.05, "RareNosig",
                                                ifelse(S100testresultfin$S100RAmean < 0.1 & S100testresultfin$S100rare_p <= 0.05, "RareSig", "Abundunt")))
    
    S70testresultfin <-c()
    S70names <- rownames(S70)
    S70mean <- rowMeans(S70)
    for (i in 1:nrow(S70)) { 
      otuname <- S70names[i]
      S70RAmean <- S70mean[i]
      result <- t.test(S70[i,],alternative ="less", mu = 1,conf.level = 0.95)
      S70rare_p <- result$p.value 
      S70rare_p <-formatC(S70rare_p, digits = 5, format = "f")
      temp <- cbind(otuname,S70RAmean,S70rare_p)
      S70testresultfin <- rbind(S70testresultfin,temp)
    }
    S70testresultfin <- data.frame(S70testresultfin)
    S70testresultfin$S70type <- ifelse(S70testresultfin$S70RAmean == 0, "NoneOccurence",
                                       ifelse(S70testresultfin$S70RAmean < 0.1 & S70testresultfin$S70rare_p > 0.05, "RareNosig",
                                              ifelse(S70testresultfin$S70RAmean < 0.1 & S70testresultfin$S70rare_p <= 0.05, "RareSig", "Abundunt")))
    
    S40testresultfin <-c()
    S40names <- rownames(S40)
    S40mean <- rowMeans(S40)
    for (i in 1:nrow(S40)) { 
      otuname <- S40names[i]
      S40RAmean <- S40mean[i]
      result <- t.test(S40[i,],alternative ="less", mu = 1,conf.level = 0.95)
      S40rare_p <- result$p.value 
      S40rare_p <-formatC(S40rare_p, digits = 5, format = "f")
      temp <- cbind(otuname,S40RAmean,S40rare_p)
      S40testresultfin <- rbind(S40testresultfin,temp)
    }
    S40testresultfin <- data.frame(S40testresultfin)
    S40testresultfin$S40type <- ifelse(S40testresultfin$S40RAmean == 0, "NoneOccurence",
                                       ifelse(S40testresultfin$S40RAmean < 0.1 & S40testresultfin$S40rare_p > 0.05, "RareNosig",
                                              ifelse(S40testresultfin$S40RAmean < 0.1 & S40testresultfin$S40rare_p <= 0.05, "RareSig", "Abundunt")))
    S10testresultfin <-c()
    S10names <- rownames(S10)
    S10mean <- rowMeans(S10)
    for (i in 1:nrow(S10)) { 
      otuname <- S10names[i]
      S10RAmean <- S10mean[i]
      result <- t.test(S10[i,],alternative ="less", mu = 1,conf.level = 0.95)
      S10rare_p <- result$p.value 
      S10rare_p <-formatC(S10rare_p, digits = 5, format = "f")
      temp <- cbind(otuname,S10RAmean,S10rare_p)
      S10testresultfin <- rbind(S10testresultfin,temp)
    }
    S10testresultfin <- data.frame(S10testresultfin)
    S10testresultfin$S10type <- ifelse(S10testresultfin$S10RAmean == 0, "NoneOccurence",
                                       ifelse(S10testresultfin$S10RAmean < 0.1 & S10testresultfin$S10rare_p > 0.05, "RareNosig",
                                              ifelse(S10testresultfin$S10RAmean < 0.1 & S10testresultfin$S10rare_p <= 0.05, "RareSig", "Abundunt")))
    sanjiDate <- merge(S10testresultfin, S40testresultfin, by = "row.names", all = TRUE)
    rownames(sanjiDate) <- sanjiDate[,1]
    sanjiDate <- sanjiDate[,c(5,9)]
    sanjiDate <- merge(sanjiDate, S70testresultfin, by = "row.names", all = TRUE)
    rownames(sanjiDate) <- sanjiDate[,1]
    sanjiDate <- sanjiDate[,c(2,3,7)]
    sanjiDate <- merge(sanjiDate, S100testresultfin, by = "row.names", all = TRUE)
    rownames(sanjiDate) <- sanjiDate[,1]
    sanjiDate <- sanjiDate[,c(2,3,4,8)]
    
    
    
    sanjiDate$summary <- paste(sanjiDate[,4],sanjiDate[,3],sanjiDate[,2],sanjiDate[,1])
    sanjiDate <- sanjiDate %>% arrange(summary)
    counttemp <- sanjiDate[1,]
    counttemp$count <- 1
    rownum <- 1
    for (cot in 2:nrow(sanjiDate)) {
      if (sanjiDate[cot,5]==counttemp[rownum,5]) {
        counttemp[rownum,6] <- counttemp[rownum,6]+1
      }else{
        rownum <- rownum+1
        counttemp[rownum,c(1:5)] <- sanjiDate[cot,c(1:5)]
        counttemp[rownum,6] <- 1
      }
    }
    counttemp[counttemp == "NoneOccurence"] <- NA
    counttemp <- na.omit(counttemp)
    counttemp <- subset(counttemp,S100type!="Abundunt")
    counttemp <- counttemp %>% arrange(S100type,S70type,S40type,S10type)
    counttemp <-counttemp[c(1:7),]
    
    
    outname <- paste0(name,"sankeyresult.csv")
    write.csv(counttemp,file =outname )
  }
  
  
  
  
  ggplot(data = counttemp,
         aes(axis1 = S100type, axis2 = S70type, axis3 = S40type,axis4 = S10type,
             y = count)) +
    scale_x_discrete(limits = c("100%", "70%", "40%", "10%")) +
    xlab("Demographic") +
    geom_alluvium(aes(fill = S10type),curve_type="xspline",knot.pos=0.05)+
    scale_fill_manual(values = c(RareSig = "
    geom_stratum(width=0.1) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
    scale_y_continuous(breaks = NULL)+
    theme(panel.grid.major = element_line(colour=NA),
          legend.position = "None"
    )
  
  
  counttemp$S10type <- paste0("S10",counttemp$S10type)
  counttemp$S40type <- paste0("S40",counttemp$S40type)
  counttemp$S70type <- paste0("S70",counttemp$S70type)
  counttemp$S100type <- paste0("S100",counttemp$S100type)
  sankeydata <- c()
  for (col in 4:2) {
    for (row in 1:nrow(counttemp)) {
      temp <- cbind(counttemp[row,c(col,col-1,6)],colnames(counttemp)[col])
      colnames(temp) <- c("Source","Target","Count","Group")
      sankeydata <- rbind(sankeydata,temp)
    }
  }
  tempsankey<- sankeydata %>%
    group_by(Group,Source, Target) %>%
    summarise(Count = sum(Count))
  write.csv(tempsankey,file = "sankey.csv")

}

{
  setwd("path/to/your/file")
  rare16s<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  
  rare<-t(rare16s)
  richness <- specnumber(rare)
  
  shannon <- diversity(rare, index = "shannon")
  chao1=estimateR(rare)[2,]
  
  invsimpson <- diversity(rare, index = "invsimpson")
  
  alpha = data.frame(richness,shannon,invsimpson,chao1,check.names = T)
  write.csv(alpha,file='16S_alpha.csv')
  
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(agricolae)
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  for (i in c(1:length(otulist))){
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-4)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    rare <- t(otu)
    richness <- specnumber(rare)
    chao1=estimateR(rare)[2,]
    
    alpha = data.frame(richness,chao1,check.names = T)
    alpha$names <- row.names(alpha)
    simi_dis <- merge(alpha, group, by.x = 'names',by.y = 'names', all.x = TRUE)
    out_filename<-sapply(savename, function(x)paste(x,"_alpha.csv",sep = "")) 
    write.csv(simi_dis,file=out_filename)
  }
  
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(agricolae)
  
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist

  
  for (i in c(1:length(otulist))){
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-4)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    rare <- t(otu)
    richness <- specnumber(rare)
    shannon <- diversity(rare, index = "shannon")
    
    alpha = data.frame(richness,shannon,check.names = T)
    alpha$names <- row.names(alpha)
    simi_dis <- merge(alpha, group, by.x = 'names',by.y = 'sample', all.x = TRUE)
    Group=c("10%re","40%re","70%re","100%re")
    
    
    simi_dis$Group <- factor(simi_dis$group,levels=c("10%re","40%re","70%re","100%re"))
    
    reresult <- c()
    for (ecosys in unique(simi_dis$ecosystem)){
      sub_phi <- subset(simi_dis,ecosystem==ecosys)
      mod1 = aov( richness~Group, data= sub_phi)
      summary(mod1)
      re = LSD.test(mod1,"Group",alpha = 0.05)
      re1 = re$groups
      re1$groupfactor <- paste0(ecosys,rownames(re1))
      print(re1)
      reresult <- rbind(reresult,re1)
    }
    reresult
    out_filename<-paste0(savename,"_richness_P.csv")
    write.csv( reresult,file=out_filename)
    p2 <- ggbarplot(simi_dis, x="ecosystem", y="richness", add = "mean_se",color = "black",fill = "Group",
                    palette = c("
      xlab("Ecosystem")+
      ylab("richness")+
      ylim(0,12000)+
      stat_compare_means(aes(group=Group), label = "p.signif", label.y = 10000, size = 10,method="anova")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black")) 
    
    p2
    out_filename<-sapply(savename, function(x)paste(x,"_richness.pdf",sep = "")) 
    ggsave( file=out_filename,p2, width = 6.5, height = 5)
    
    reresult <- c()
    for (ecosys in unique(simi_dis$ecosystem)){
      sub_phi <- subset(simi_dis,ecosystem==ecosys)
      mod1 = aov( shannon~Group, data= sub_phi)
      summary(mod1)
      re = LSD.test(mod1,"Group",alpha = 0.05)
      re1 = re$groups
      re1$groupfactor <- paste0(ecosys,rownames(re1))
      print(re1)
      reresult <- rbind(reresult,re1)
    }
    reresult
    out_filename<-paste0(savename,"_shannon_P.csv")
    write.csv( reresult,file=out_filename)

    p4 <- ggbarplot(simi_dis, x="ecosystem", y="shannon", add = "mean_se",color = "black",fill = "Group",
                    palette = c("
      xlab("Ecosystem")+
      ylab("shannon")+
      ylim(0,10)+
      stat_compare_means(aes(group=Group), label = "p.signif", label.y = 8, size = 10,method="anova")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black")) 
    
    p4 
    out_filename<-sapply(savename, function(x)paste(x,"_shannon.pdf",sep = "")) 
    ggsave(file=out_filename,p4,  width = 6.5, height = 5)
  }
  
  
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', header=TRUE)
  setwd("path/to/your/file")
  otulist = list.files(pattern="*.csv")
  otulist
  
  
  for (i in c(1:length(otulist))){
    name <-otulist[i]
    savename <- substr(name,1,nchar(name)-4)
    otu<-read.csv(otulist[i], header=TRUE, row.names=1)
    rare <- t(otu)
    richness <- specnumber(rare)
    shannon <- diversity(rare, index = "shannon")
    
    alpha = data.frame(richness,shannon,check.names = T)
    alpha$names <- row.names(alpha)
    simi_dis <- merge(alpha, group, by.x = 'names',by.y = 'sample', all.x = TRUE)
    Group=c("10%re","40%re","70%re","100%re")
    
    
    simi_dis$Group <- factor(simi_dis$group,levels=c("10%re","40%re","70%re","100%re"))
    
    reresult <- c()
    for (ecosys in unique(simi_dis$ecosystem)){
      sub_phi <- subset(simi_dis,ecosystem==ecosys)
      mod1 = aov( richness~Group, data= sub_phi)
      summary(mod1)
      re = LSD.test(mod1,"Group",alpha = 0.05)
      re1 = re$groups
      re1$groupfactor <- paste0(ecosys,rownames(re1))
      print(re1)
      reresult <- rbind(reresult,re1)
    }
    reresult
    out_filename<-paste0(savename,"_richness_P.csv")
    write.csv( reresult,file=out_filename)
    p2 <- ggbarplot(simi_dis, x="ecosystem", y="richness", add = "mean_se",color = "black",fill = "Group",
                    palette = c("
      xlab("Ecosystem")+
      ylab("richness")+
      ylim(0,2000)+
      stat_compare_means(aes(group=Group), label = "p.signif", label.y = 1500, size = 10,method="anova")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black")) 
    
    p2
    out_filename<-sapply(savename, function(x)paste(x,"_richness.pdf",sep = "")) 
    ggsave( file=out_filename,p2, width = 6.5, height = 5)
    
    reresult <- c()
    for (ecosys in unique(simi_dis$ecosystem)){
      sub_phi <- subset(simi_dis,ecosystem==ecosys)
      mod1 = aov( shannon~Group, data= sub_phi)
      summary(mod1)
      re = LSD.test(mod1,"Group",alpha = 0.05)
      re1 = re$groups
      re1$groupfactor <- paste0(ecosys,rownames(re1))
      print(re1)
      reresult <- rbind(reresult,re1)
    }
    reresult
    out_filename<-paste0(savename,"_shannon_P.csv")
    write.csv( reresult,file=out_filename)
    
    p4 <- ggbarplot(simi_dis, x="ecosystem", y="shannon", add = "mean_se",color = "black",fill = "Group",
                    palette = c("
      xlab("Ecosystem")+
      ylab("shannon")+
      ylim(0,6)+
      stat_compare_means(aes(group=Group), label = "p.signif", label.y = 5, size = 10,method="anova")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black")) 
    
    p4 
    out_filename<-sapply(savename, function(x)paste(x,"_shannon.pdf",sep = "")) 
    ggsave(file=out_filename,p4,  width = 6.5, height = 5)
  }
}

{
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  setwd("path/to/your/file")
  similarity <- read.csv("fin_dis_zhu.csv", header=TRUE,row.names = 1)
  alphaD <- read.csv("16S_alpha.csv", header=TRUE)
  
  setwd("path/to/your/file")
  similarity <- read.csv("fin_dis_zhu.csv", header=TRUE,row.names = 1)
  alphaD <- read.csv("ITS_alpha.csv", header=TRUE)
  
  similarity <- similarity %>% arrange(X1,X2)
  tempsankey<- similarity %>%
    group_by(X1) %>%
    reframe(similarity = mean(similarity),site=unique(site),Group=unique(Group),group1=unique(group1))
  merge_data <- merge(alphaD,tempsankey,by.x="X",by.y="X1")
  dot_data <- merge_data
  dot_data <- subset(merge_data,site=="T")
  dot_data$Group <- factor(dot_data$Group,levels=c("100%re","70%re","40%re","10%re"))

  
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  p_richness <- ggplot(data=dot_data,mapping=aes( x=richness, y=similarity))+
    geom_point(aes(color =Group ,shape =group1 ), size = 3, alpha =2) + 
    scale_shape_manual(values = c(15,17,18,16,3,8)) + 
    scale_color_manual(values = c("
    geom_smooth(method = 'lm', formula = y ~  x + I(x ^ 2), se = T,colour="
    stat_cor(data=dot_data, method = "pearson",r.digits = 3,p.accuracy=0.01)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size = 20,face = "plain",vjust = -7))+  
    theme(axis.text.x = element_text(size = 15,face = "plain"))+
    theme(axis.text.y = element_text(size = 15,face = "plain"))+
    theme(axis.title.x = element_text(size = 15, face = "plain"))+
    theme(axis.title.y = element_text(size = 15, face = "plain"))+
    theme(legend.position = "None")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
    theme(axis.ticks.x=element_line(color="black",linewidth=1.1,lineend = 1))+
    theme(axis.ticks.y=element_line(color="black",linewidth=1.1,lineend = 1))+
    my_theme
  p_richness
  ggsave(p_richness, file="Corp_its_相似与丰富度.pdf", width = 6, height = 6)
  
  p_chao1 <-ggplot(data=dot_data,mapping=aes(  x=richness, y=similarity,color=Group))+
    geom_point(aes(color =Group ,shape =group1 ), size = 3, alpha =2) + 
    scale_shape_manual(values = c(15,17,18,16,3,8)) + 
    scale_color_manual(values = c("
    geom_smooth(method = 'lm', formula = y ~  x + I(x ^ 2), se = T,linewidth=1.1,alpha=0.5)+
    stat_cor(data=dot_data, method = "pearson",r.digits = 3,p.accuracy=0.01)+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5,size = 20,face = "plain",vjust = -7))+  
    theme(axis.text.x = element_text(size = 15,face = "plain"))+
    theme(axis.text.y = element_text(size = 15,face = "plain"))+
    theme(axis.title.x = element_text(size = 15, face = "plain"))+
    theme(axis.title.y = element_text(size = 15, face = "plain"))+
    theme(legend.position = "None")+
    theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+
    theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = 1))+
    theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = 1))+
    my_theme
  p_chao1
  ggsave(p_chao1, file="16s_相似与丰富度.pdf", width = 6, height = 6)
  
  
  c("
  
  c("

}

{
  library(dplyr)
  library(ggpubr)
  library(ggplot2)
  library(agricolae)
  setwd("path/to/your/file")
  copynum <- read.csv("16S_qPCR.csv", header=TRUE)

  setwd("path/to/your/file")
  copynum <- read.csv("ITS_qPCR.csv", header=TRUE)

  copynum <- na.omit(copynum)
  meancopynum<- copynum %>%
    group_by(sample_id) %>%
    reframe(sq = mean(sq))
  write.csv(meancopynum,file = "16S_copynum.csv")
  
  write.csv(meancopynum,file = "ITS_copynum.csv")
  
  setwd("path/to/your/file")
  fin_copynum <- read.csv("16S_copynum.csv", header=TRUE)
  

  setwd("path/to/your/file")
  fin_copynum <- read.csv("ITS_copynum.csv", header=TRUE)
  
  
  fin_copynum <-   fin_copynum[,c(9,10,11,13)]
  
  
  fin_copynum$copynum <- log10(fin_copynum$copynum)
  
  reresult <- c()
  for (ecosys in unique(fin_copynum$ecosystem)){
    sub_phi <- subset(fin_copynum,ecosystem==ecosys)
    mod1 = aov( copynum~group, data=fin_copynum)
    summary(mod1)
    re = LSD.test(mod1,"group",alpha = 0.05)
    re1 = re$groups
    re1$groupfactor <- paste0(ecosys,rownames(re1))
    print(re1)
    reresult <- rbind(reresult,re1)
  } 
  reresult
  
  fin_copynum$group <- factor(fin_copynum$group,levels=c("100%re","70%re","40%re","10%re"))
  compare_means(copynum~group, data=fin_copynum, group.by = "ecosystem",method = "anova")
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p_16s <- ggbarplot(fin_copynum, x="ecosystem", y="copynum", add = "mean_se",color = "black",fill = "group",
                     palette = c("
    ylab("")+
    xlab("Ecosystem")+
    ylim(0,12)+
    stat_compare_means(aes(group=group), label = "p.signif", label.y = 10.5, size = 8,method = "anova")+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.position = "None")+
    my_theme
  p_16s
  ggsave('qPCR_16s.pdf', p_16s, width = 6, height = 6)
  
  p_its <- ggbarplot(fin_copynum, x="ecosystem", y="copynum", add = "mean_se",color = "black",fill = "group",
                     palette = c("
    ylab("")+
    xlab("Ecosystem")+
    ylim(0,10)+
    stat_compare_means(aes(group=group), label = "p.signif", label.y = 8.5, size = 8,method = "anova")+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.position = "None")+
    my_theme
  p_its
  ggsave('qPCR_its.pdf', p_its, width = 6, height = 6)
}

{
  library(ggpubr)
  library(ggplot2)
  library(agricolae)
  
  setwd("path/to/your/file")
  fin_copynum <- read.csv("验证_16s_qPCR.csv", header=TRUE)
  
  
  setwd("path/to/your/file")
  fin_copynum <- read.csv("验证_its_qPCR.csv", header=TRUE)
  
  
  fin_copynum$copynum <- log10(fin_copynum$copynum)
  fin_copynum$group <- factor(fin_copynum$group,levels=c("0%","10%-D","40%-D","70%-D","100%-D","100%"))
  compare_means(copynum~group, data=fin_copynum, group.by = "Ecosystem",method = "anova")
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p_16s <- ggbarplot(fin_copynum, x="Ecosystem", y="copynum", add = "mean_se",color = "black",fill = "group",
                     palette = c("
    ylab("Copies")+
    ylim(0,12)+
    stat_compare_means(aes(group=group), label = "p.signif", label.y = 10.5, size = 8)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.position = "None")+
    my_theme
  p_16s
  ggsave('验证qPCR_16s.pdf', p_16s, width = 6, height = 6)
  
  p_its <- ggbarplot(fin_copynum, x="Ecosystem", y="copynum", add = "mean_se",color = "black",fill = "group",
                     palette = c("
    ylab("Copies")+
    
    ylim(0,12)+
    stat_compare_means(aes(group=group), label = "p.signif", label.y = 8.5, size = 8)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.position = "None")+
    my_theme
  p_its
  ggsave('验证qPCR_its.pdf', p_its, width = 6, height = 6)
}

{
  rm(list=ls())
  library(tidyverse)
  library(microeco)
  library(magrittr)
  
  
  library(snowfall) 
  library(parallel)

  

  setwd("path/to/your/file")
  setwd("path/to/your/file")
  tax_table <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  sample_table <- read.csv('LEfSE分组.csv', row.names = 1, check.names = FALSE)
  wd <- ("path/to/your/file")
  wd <- ("path/to/your/file")
  
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  
  otulist = list.files(pattern="*.csv")
  
  cpus = 12
  sfInit(parallel = TRUE, cpus)
  sfExport("tax_table","sample_table","otulist","wd")
  sfLibrary(tidyverse)
  sfLibrary(microeco)
  sfLibrary(magrittr)
  
  
  paracacl_lefse <- function(otudivnum){
    feature_table<- read.csv(otulist[otudivnum], row.names = 1)
    name <- substr(otulist[otudivnum],1,nchar(otulist[otudivnum])-7)
    otu_colname <- colnames(feature_table)
    sample_table$rowselect <- rownames(sample_table)
    sample_table <- sample_table[sample_table$rowselect %in% otu_colname, ]
    
    head(feature_table)[1:6,1:6]; head(sample_table)[1:6, ]; head(tax_table)[,1:6]
    
    dataset <- microtable$new(sample_table = sample_table,
                              otu_table = feature_table, 
                              tax_table = tax_table)
    dataset
    lefse <- trans_diff$new(dataset = dataset, 
                            method = "lefse", 
                            group = "group1", 
                            alpha =  0.05, 
                            lefse_subgroup = NULL)
    head(lefse$res_diff)
    pbox <- lefse$plot_diff_bar(use_number = 1:30, 
                                width = 0.8, 
                                group_order = c("CK","C","N","pH","T","W" )) +
      ggsci::scale_color_npg() +
      ggsci::scale_fill_npg()
    
    pbox
    outfilename <- paste0(name,"差异柱状图.pdf")
    ggsave(outfilename, path = wd, pbox, width = 10, height = 8)

    ptree<- lefse$plot_diff_cladogram(use_taxa_num = 200, 
                                      use_feature_num = 50, 
                                      clade_label_level = 5, 
                                      group_order = c( "CK","C","N","pH","T","W" ))
    ptree
    outfilename <- paste0(name,"分类树状图.pdf")
    ggsave(outfilename, path = wd,ptree, width = 10, height = 8)

  }
  
  sfLapply(1:12, paracacl_lefse) 
  
  use_labels <- c("c__Deltaproteobacteria", "c__Actinobacteria", "o__Rhizobiales", "p__Proteobacteria", "p__Bacteroidetes", 
                  "o__Micrococcales", "p__Acidobacteria", "p__Verrucomicrobia", "p__Firmicutes", 
                  "p__Chloroflexi", "c__Acidobacteria", "c__Gammaproteobacteria", "c__Betaproteobacteria", "c__KD4-96",
                  "c__Bacilli", "o__Gemmatimonadales", "f__Gemmatimonadaceae", "o__Bacillales", "o__Rhodobacterales")
  lefse$plot_diff_cladogram(use_taxa_num = 200, 
                            use_feature_num = 50, 
                            select_show_labels = use_labels)
}

{
  library(ggpubr)
  library(dplyr)
    
  setwd("path/to/your/file")

  otulist = list.files(pattern="*.csv")
  
  for (nst_result in otulist[1:3]) {
    data<- read.csv(nst_result, row.names = 1)
    name <- substr(nst_result,1,nchar(nst_result)-4)
    data <- data %>%
      filter(name1 != name2)
    group1 <- read.csv("path/to/your/file",  header=T)
    group2 <- read.csv("path/to/your/file",  header=T)
    simi_dis <- merge(data, group1, by.x = 'name1',by.y = 'names', all.x = TRUE)
    simi_dis <- merge(simi_dis, group2, by.x = 'name2',by.y = 'names_2', all.x = TRUE)
    simi_dis <- simi_dis %>%
      filter(ecosystem == ecosystem_2)
    simi_dis <- simi_dis %>%
      filter(group.y == group_2)
    mytheme = theme_classic() + theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
      theme(axis.title.y= element_text(size=8))+theme(axis.title.x = element_text(size = 8))+
      theme(legend.title=element_text(size=8),legend.text=element_text(size=8))
    my_comparisons <- list(c("10%re", "100%re"), c("40%re", "100%re"),c("70%re", "100%re"))
    
    col=c("
    
    p = ggplot(simi_dis,aes(fill = group.y,x=factor(group.y,level=c("10%re","40%re","70%re","100%re" )), y=NST.ij.ruzicka))+
      geom_violin(position = position_dodge(width = 0.1), scale = 'width')+  
      geom_violin(position = position_dodge(width = 0.1), scale = 'width')+  
      geom_hline(aes(yintercept=0.5), colour="gray", linetype="dashed")+
      geom_boxplot(alpha=1,outlier.size=0, size=0.3, width=0.3,fill="white") +
      scale_fill_manual(values = col)+
      labs(x="Struation", y="Normalized Stochasticity ratio(NST)", color=group)+
      scale_x_discrete(limits=c("10%re","40%re","70%re","100%re" ))
   
    p2=p+mytheme+
      stat_compare_means(comparisons=my_comparisons,  label ="p.signif",method = "wilcox.test")+
      theme(legend.position = "none")
    p2 
    outputname <- paste0(name,"NST_violin.pdf")
    ggsave(outputname, p2, width = 6, height =6,path = "path/to/your/file",units = "cm") 
  }
  
  for (nst_result in otulist[4:6]) {
    data<- read.csv(nst_result, row.names = 1)
    name <- substr(nst_result,1,nchar(nst_result)-4)
    group1 <- read.csv("path/to/your/file",  header=T)
    group2 <- read.csv("path/to/your/file",  header=T)
    simi_dis <- merge(data, group1, by.x = 'name1',by.y = 'names', all.x = TRUE)
    simi_dis <- merge(simi_dis, group2, by.x = 'name2',by.y = 'names_2', all.x = TRUE)
    simi_dis <- simi_dis %>%
      filter(ecosystem == ecosystem_2)
    simi_dis <- simi_dis %>%
      filter(group.y == group_2)
    mytheme = theme_classic() + theme(axis.text.x = element_text(size = 8),axis.text.y = element_text(size = 8))+
      theme(axis.title.y= element_text(size=8))+theme(axis.title.x = element_text(size = 8))+
      theme(legend.title=element_text(size=8),legend.text=element_text(size=8)
            )
    my_comparisons <- list(c("10%re", "100%re"), c("40%re", "100%re"),c("70%re", "100%re"))
    
    col=c("
    
    p = ggplot(simi_dis,aes(fill = group.y,x=factor(group.y,level=c("10%re","40%re","70%re","100%re" )), y=NST.ij.ruzicka))+
      geom_violin(position = position_dodge(width = 0.1), scale = 'width')+  
      geom_violin(position = position_dodge(width = 0.1), scale = 'width')+  
      geom_hline(aes(yintercept=0.5), colour="gray", linetype="dashed")+
      geom_boxplot(alpha=1,outlier.size=0, size=0.3, width=0.3,fill="white") +
      scale_fill_manual(values = col)+
      labs(x="Struation", y="Normalized Stochasticity ratio(NST)", color=group)+
      scale_x_discrete(limits=c("10%re","40%re","70%re","100%re" ))
    
    p2=p+mytheme+
      stat_compare_means(comparisons=my_comparisons,  label ="p.signif",method = "wilcox.test")+
      theme(legend.position = "none")
    p2 
    outputname <- paste0(name,"NST_violin.pdf")
    ggsave(outputname, p2, width = 6, height =6,path = "path/to/your/file",units = "cm") 
  }
}

{
  
  library(vegan)
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  group <- read.csv('分组.csv', row.names = 1, header=T)
  
  setwd("path/to/your/file")
  setwd("path/to/your/file")
  
  otulist = list.files(pattern="*.csv")
  otu <- read.csv(otulist[1], row.names = 1)
  otu_colname <- colnames(otu)
  group$rowselect <- rownames(group)
  group <- group[group$rowselect %in% otu_colname, ]
  otu <- data.frame(t(otu))
  adonis_result_otu <- adonis2(otu~group1, group, permutations = 999, distance = 'bray' )
  adonis_result_otu
  otuput <- data.frame(adonis_result_otu$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
  otuput <- cbind(rownames(otuput), otuput)
  names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
  write.table(otuput, file = 'PERMANOVA.result_all.txt', row.names = FALSE, sep = "path/to/your/file", quote = FALSE, na = '')
}

{
  library(permute)
  library(vegan)
  library(ggpubr)
  library(ggplot2)
  library(agricolae)
  library(reshape2)
  
  setwd("path/to/your/file")
  ado1<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  
  setwd("path/to/your/file")
  ado1<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  
  otusum<-unname(colSums(ado1)[1])
  
  ado=(t(ado1)/otusum)^0.5
  prok.dist <- vegdist(ado)
  
  prok.dist <- as.matrix(prok.dist)
  prok.dist <- as.data.frame(prok.dist)
  prok.dist$X1 <- rownames(prok.dist)
  prok.matr<- melt(prok.dist,id.vars = "X1", variable.name = "X2", value.name = "score") 
  
  colnames(prok.matr)[2] <-"X2"
  colnames(prok.matr)[3] <-"similarity" 
  prok.matr$similarity <- 1-prok.matr$similarity
  reselect <-prok.matr
  
  group1<-read.csv("相似性柱状图分组1.csv", header=TRUE)
  group2<-read.csv("相似性柱状图分组2.csv", header=TRUE)
  simi_dis <- merge(reselect, group1, by.x = 'X1',by.y = 'names', all.x = TRUE)
  simi_dis <- merge(simi_dis, group2, by.x = 'X2',by.y = 'names_2', all.x = TRUE)
  
  fin_dis_zhu<-subset(simi_dis,ecosystem_2==ecosystem)
  fin_dis_zhu<-subset(fin_dis_zhu,group_2==group)
  fin_dis_zhu<-subset(fin_dis_zhu,group1_2==group1)
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="ecosystem")] <-"site" 
  colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="group")] <-"Group" 
  
  write.csv(fin_dis_zhu,file='fin_dis_within-treat.csv')
  phi<-read.csv("fin_dis_within-treat.csv", header=TRUE, row.names=1)
  
  reresult <- c()
  for (ecosys in unique(phi$site)){
    sub_phi <- subset(phi,site==ecosys)
    mod1 = aov( similarity~Group, data= sub_phi)
    summary(mod1)
    re = LSD.test(mod1,"Group",alpha = 0.05)
    re1 = re$groups
    re1$groupfactor <- paste0(ecosys,rownames(re1))
    print(re1)
    reresult <- rbind(reresult,re1)
  }
  reresult
  compare_means(similarity~Group, data=phi, group.by = "site",method = "anova")
  phi$Group = factor(phi$Group, levels = c("10%re","40%re","70%re","100%re"))
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  p_16s <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",
                     palette = c("
    ylab("Similarity")+
    xlab("Ecosystem")+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.65, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  p_16s
  
  ggsave('within-treat-similarity_16s.pdf', p_16s, width = 9, height = 9)
  
  its_p <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",
                     palette = c("
    ylab("Similarity")+
    xlab("Ecosystem")+
    ylim(0,0.8)+
    stat_compare_means(aes(group=Group), label = "p.signif", label.y = 0.70, size = 10)+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  its_p
  
  ggsave('within-treat-similarity_its.pdf', its_p, width = 9, height = 9)
}

{
  library(dplyr)

  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  otu_table <- t(otu)
  otu_table <- as.data.frame(otu_table)
  otu_table$name <- rownames(otu_table)
  sample_otu <- merge(otu_table, group, by.x = 'name',by.y = 'sample', all.x = TRUE)
  rownames(sample_otu) <-sample_otu[,1]
  sample_otu <- sample_otu[,-1]

  filtered_data <- sample_otu[sample_otu$group1 %in% c("CK", "C"), ]
  out_catch <-t(filtered_data)[,1:2]
  out_catch <- as.data.frame(out_catch)
  out_catch <- out_catch[-c(44779:44781),]
  for (ecosystem in c("C", "F", "G")) {
    for (group in c("10%", "40%", "70%", "100%")) {
        inner_catch <- c()
        subset_data <- filtered_data %>%
          filter(ecosystem == ecosystem, group == group)
        subset_data <- subset_data[,-(44779:44780)]
        subset_data_CK <- subset_data %>%
          filter( group1 == "CK")
        subset_data_CK <- subset_data_CK[,-44779]
        subset_data_C <- subset_data %>%
          filter( group1 == "C")
        subset_data_C <- subset_data_C[,-44779]
        inner_catch <- cbind(colMeans(subset_data_CK),colMeans(subset_data_C))
        inner_catch <- as.data.frame(inner_catch)
        inner_catch$differ <- inner_catch[,1]-inner_catch[,2]

        inner_catch$RK[inner_catch$differ > 0] <- "R"
        inner_catch$RK[inner_catch$differ < 0] <- "K"
        inner_catch$RK[inner_catch$differ == 0] <- "N"
        out_catch <- cbind.data.frame(out_catch,inner_catch$RK)
        group_name <- paste(ecosystem, group, sep = "-")
        colnames(out_catch)[ncol(out_catch)] <- group_name
       
    }
  }
  out_catch <- out_catch[,-c(1:2)]
  write.csv(out_catch, file = "rk_otubase.csv", row.names = T)
  out_catch$rk <- apply(out_catch, 1, function(row) paste(row, collapse = ""))
  rk_otu <- out_catch[,13,drop=F]
  taxonomy <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  TAXA_RK <- merge(rk_otu,taxonomy,by="row.names")
  result <- table(TAXA_RK$rk, TAXA_RK$`phylum	`)
  result <- as.data.frame(result)
  write.csv(result, file = "rk_门统计.csv", row.names = T)
  data<-read.csv("rk_门统计.csv", header=TRUE, row.names=1)
  filter_data <- data[data$Var2 %in% c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", 
                                       "Bdellovibrionota", "Chloroflexi", "Myxococcota", "Planctomycetota", 
                                       "Proteobacteria", "Verrucomicrobiota"), ]
  filter_data$Var2 = factor(filter_data$Var2, levels = c("Acidobacteriota", "Actinobacteriota", "Armatimonadota", "Bacteroidota", 
                                                         "Bdellovibrionota", "Chloroflexi", "Myxococcota", "Planctomycetota", 
                                                         "Proteobacteria", "Verrucomicrobiota"))
  library(ggpubr)
  library(ggplot2)
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p_rk <- ggbarplot(filter_data, x="Var2", y="Freq",color = "black",fill = "Var1",
                     palette = c("
    ylab("总数")+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(angle=90,size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  p_rk
  ggsave('rk_16s.pdf', p_rk, width = 9, height =18)
  
  
  otu_table <- t(otu)
  otu_table <- as.data.frame(otu_table)
  otu_table$name <- rownames(otu_table)
  sample_otu <- merge(otu_table, group, by.x = 'name',by.y = 'sample', all.x = TRUE)
  rownames(sample_otu) <-sample_otu[,1]
  sample_otu <- sample_otu[,-1]
  
  filtered_data <- sample_otu[sample_otu$group1 %in% c("CK", "C"), ]
  out_catch <-t(filtered_data)[,1:2]
  out_catch <- as.data.frame(out_catch)
  out_catch <- out_catch[-c(9518:9521),]
  for (ecosystem in c("Grass", "Corland", "forest")) {
    for (group in c("10%", "40%", "70%", "100%")) {
      inner_catch <- c()
      subset_data <- filtered_data %>%
        filter(ecosystem == ecosystem, group == group)
      subset_data <- subset_data[,-(9518:9520)]
      subset_data_CK <- subset_data %>%
        filter( group1 == "CK")
      subset_data_CK <- subset_data_CK[,-9518]
      subset_data_C <- subset_data %>%
        filter( group1 == "C")
      subset_data_C <- subset_data_C[,-9518]
      inner_catch <- cbind(colMeans(subset_data_CK),colMeans(subset_data_C))
      inner_catch <- as.data.frame(inner_catch)
      inner_catch$differ <- inner_catch[,1]-inner_catch[,2]
      
      inner_catch$RK[inner_catch$differ > 0] <- "R"
      inner_catch$RK[inner_catch$differ < 0] <- "K"
      inner_catch$RK[inner_catch$differ == 0] <- "N"
      out_catch <- cbind.data.frame(out_catch,inner_catch$RK)
      group_name <- paste(ecosystem, group, sep = "-")
      colnames(out_catch)[ncol(out_catch)] <- group_name
      
    }
  }
  out_catch <- out_catch[,-c(1:2)]
  write.csv(out_catch, file = "rk_otubase.csv", row.names = T)
  out_catch$rk <- apply(out_catch, 1, function(row) paste(row, collapse = ""))
  rk_otu <- out_catch[,13,drop=F]
  taxonomy <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  TAXA_RK <- merge(rk_otu,taxonomy,by="row.names")
  result <- table(TAXA_RK$rk, TAXA_RK$class)
  result <- as.data.frame(result)
  write.csv(result, file = "rk_纲统计.csv", row.names = T)
  data<-read.csv("rk_纲统计.csv", header=TRUE, row.names=1)

  
  filter_data <- data[data$Var2 %in% c("Agaricomycetes", "Dothideomycetes", "Eurotiomycetes","GS27",  "Leotiomycetes", 
                                       "Orbiliomycetes", "Pezizomycetes", "Sordariomycetes"), ]
  filter_data$Var2 = factor(filter_data$Var2, levels = c("Agaricomycetes", "Dothideomycetes", "Eurotiomycetes","GS27",  "Leotiomycetes", 
                                                         "Orbiliomycetes", "Pezizomycetes", "Sordariomycetes"))
  library(ggpubr)
  library(ggplot2)
  my_theme <- theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  p_rk <- ggbarplot(filter_data, x="Var2", y="Freq",color = "black",fill = "Var1",
                    palette = c("
    ylab("总数")+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(angle=90,size = 20,color="black",vjust = 1.0, hjust = 1.0))+ 
    my_theme
  p_rk
  ggsave('rk_its.pdf', p_rk, width = 9, height =18)
}

{
  library(esc)
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  otu <- as.data.frame(t(otu))
  sumcot <- unique(rowSums(otu))
  colcounter <- ncol(otu)
  otu$names <-rownames(otu)
  input = merge(otu,group,by.x = "names",by.y ="sample" )
  rownames(input) <- input$names
  resultslis <- c()
  for (ecosys in unique(input$ecosystem)) {
    tempheat <- subset(input,input$ecosystem==ecosys)
    for (satruation in unique(tempheat$group)){
      temp_satru <- subset(tempheat,tempheat$group==satruation)
      data <- temp_satru[,c(2:(colcounter+1))]
      data <- data[, colSums(data) != 0]
      gap <- tempheat[,c(1,(colcounter+1):ncol(temp_satru))]
      for (col in 1:ncol(data)) {
        data_otu <- data[,col,drop=F]
        rare_otu <- data_otu
        rare_otu <- rare_otu/sumcot
        all_negative <- all(rare_otu[[1]] < 0.0001)
        print(rare_otu[1])
        print(all_negative)
        if (all_negative) {
          rare <- "Rare"
        }else{
          rare <- "Not"
        }
        data_otu$name <- rownames(data_otu)
        fin_data <-  merge(data_otu,group,by.x = "name",by.y ="sample" )
        for (treat in c("C","N","pH","T","W")) {
          select_treat <- c("CK",treat)
          hedges_fin <- fin_data[fin_data$group1 %in% select_treat,]
          
          sd_hedges_fin_ck <- sd(hedges_fin[hedges_fin$group1 %in% "CK",2])
          mean_hedges_fin_ck <-  mean(hedges_fin[hedges_fin$group1 %in% "CK",2])
          
          sd_hedges_fin_treat <- sd(hedges_fin[hedges_fin$group1 %in% treat,2])
          mean_hedges_fin_treat <-  mean(hedges_fin[hedges_fin$group1 %in% treat,2])
          
          hedges_g_10<- esc_mean_sd(grp1m = sd_hedges_fin_ck, grp1sd =mean_hedges_fin_ck, grp1n = 4,grp2m = sd_hedges_fin_treat, grp2sd =mean_hedges_fin_treat, grp2n = 4, es.type = "g")
          resultcatcher <- c(ecosys,satruation,colnames(hedges_fin)[2],hedges_g_10$es,treat,rare)
          resultslis <- rbind(resultslis,resultcatcher)
          }
      }
    }
  }
  colnames(resultslis) <- c("ecosystem","satruation","OTU_ID","hedges","treat","rare")
  write.csv(resultslis, file = "ITS_hedges_G.csv", row.names = F)
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  otu <- otu/colSums(otu)[1]
  otu <- as.data.frame(t(otu))
  colcounter <- ncol(otu)
  
  otu$sample <- rownames(otu)
  input = merge(otu,group,by.x = "sample",by.y ="sample" )
  for (ecosys in unique(input$ecosystem)) {
    tempheat <- subset(input,input$ecosystem==ecosys)
    rownames(tempheat) <- tempheat$sample
    tempheat <- tempheat[,2:(colcounter+1)]
    tempdata <- tempheat[, colSums(tempheat) != 0]
    tempdata <- as.data.frame(t(tempdata))
    tempdata$meanRA <- rowSums(tempdata)
    tempdata$otu_ID <- rownames(tempdata)
    sumra_data <- tempdata[,c("otu_ID","meanRA")]
    sumra_data <- merge(sumra_data,taxonomy,by.x = "otu_ID",by.y ="FeatureID" )
    sumra_result<- sumra_data %>%
      group_by(phylum) %>%
      summarise(sumra = sum(meanRA))
    sumra_result$sumra <- round(sumra_result$sumra,digits = 5)
    sumra_result <- sumra_result[order(sumra_result$sumra,decreasing = T),]
    write.csv(sumra_result,file =paste0("生态系统",ecosys,"相对多度排名.csv"))
  }
  

  setwd("path/to/your/file")
  data<-read.csv("16S_hedges_G.csv", header=TRUE)
  taxonomy <- read.delim('metadata.txt',check.names = FALSE)

  colnames(taxonomy) <- c("FeatureID", "kingdom", "phylum" , "class",   "order" ,  "family",  "genus" ,  "species")
 
  
  data <- na.omit(data)

  input = merge(data,taxonomy,by.x = "OTU_ID",by.y ="FeatureID" )
  library(dplyr)
  tempsankey<- input %>%
    group_by(ecosystem,satruation, phylum,treat,rare) %>%
    summarise(hedgesg = mean(hedges))
  library(pheatmap)
  library(tidyr)
  tempsankey$type <- paste0(tempsankey$ecosystem,tempsankey$satruation,tempsankey$treat,tempsankey$rare)
  for (ecosys in unique(tempsankey$ecosystem)) {
    tempdata<-  subset(tempsankey,tempsankey$ecosystem==ecosys)
    taxagroup <- read.csv(file =paste0("生态系统",ecosys,"相对多度排名.csv"),header=T ,row.names = 1)
    selecttaxa <- taxagroup$phylum[1:14]
    tempheat <- tempdata[tempdata$phylum %in% selecttaxa,]
    tempheat <- tempheat[,3:7]
    tempheat$phylum <- substr(tempheat$phylum,4,nchar(tempheat$phylum))
    
    otherdata <- tempdata[!tempdata$phylum %in% selecttaxa,]
    otherdata<- otherdata %>%
      group_by(ecosystem,satruation,treat,rare) %>%
      summarise(hedgesg = mean(hedgesg))
    otherdata$phylum <- "Others"
    otherdata$type <- paste0(otherdata$ecosystem,otherdata$satruation,otherdata$treat,otherdata$rare)
    otherdata <- otherdata[,c("phylum","treat","hedgesg","rare","type")]
    
    tempheat <- rbind(tempheat,otherdata)
    tempheat <- tempheat[c("phylum","hedgesg","type")]
    
    heatdata <- spread(tempheat,key="type",value ="hedgesg" )
    heatdata <- as.data.frame(heatdata)

    row_index <- which(heatdata$phylum== "Others",arr.ind =F)
    
    heatdata <- rbind(heatdata[-row_index, ], heatdata[row_index, ])
    
    rownames(heatdata) <- heatdata$phylum
    heatdata <- heatdata[,-1]
    sortcol <- c(paste0(ecosys,"10%reCRare"),paste0(ecosys,"10%reCNot"),paste0(ecosys,"40%reCRare"),paste0(ecosys,"40%reCNot"),
                 paste0(ecosys,"70%reCRare"),paste0(ecosys,"70%reCNot"),paste0(ecosys,"100%reCRare"),paste0(ecosys,"100%reCNot"),
                 paste0(ecosys,"10%reNRare"),paste0(ecosys,"10%reNNot"),paste0(ecosys,"40%reNRare"),paste0(ecosys,"40%reNNot"),
                 paste0(ecosys,"70%reNRare"),paste0(ecosys,"70%reNNot"),paste0(ecosys,"100%reNRare"),paste0(ecosys,"100%reNNot"),
                 paste0(ecosys,"10%repHRare"),paste0(ecosys,"10%repHNot"),paste0(ecosys,"40%repHRare"),paste0(ecosys,"40%repHNot"),
                 paste0(ecosys,"70%repHRare"),paste0(ecosys,"70%repHNot"),paste0(ecosys,"100%repHRare"),paste0(ecosys,"100%repHNot"),
                 paste0(ecosys,"10%reTRare"),paste0(ecosys,"10%reTNot"), paste0(ecosys,"40%reTRare"),paste0(ecosys,"40%reTNot"),
                 paste0(ecosys,"70%reTRare"),paste0(ecosys,"70%reTNot"),paste0(ecosys,"100%reTRare"),paste0(ecosys,"100%reTNot"),
                 paste0(ecosys,"10%reWRare"),paste0(ecosys,"10%reWNot"),paste0(ecosys,"40%reWRare"),paste0(ecosys,"40%reWNot"),
                 paste0(ecosys,"70%reWRare"),paste0(ecosys,"70%reWNot"),paste0(ecosys,"100%reWRare"),paste0(ecosys,"100%reWNot"))
    heatdata <- heatdata[,sortcol]
    pheatmap(
      heatdata,
      scale = "row",
      cellwidth = 10, cellheight = 10,
      color = colorRampPalette(c("
      gaps_col = c(8,16,24,32,40),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "black",
      filename = paste0("稀有_生态系统_",ecosys,"_top12taxa.pdf")
    ) 
  }

  
  setwd("path/to/your/file")
  data<-read.csv("ITS_hedges_G.csv", header=TRUE)
  taxonomy <- read.delim('metadata.txt', check.names = FALSE)
  
  colnames(taxonomy) <- c("FeatureID", "kingdom", "phylum" , "class",   "order" ,  "family",  "genus" ,  "species")
  
  
  data <- na.omit(data)
  
  input = merge(data,taxonomy,by.x = "OTU_ID",by.y ="FeatureID" )
  library(dplyr)
  tempsankey<- input %>%
    group_by(ecosystem,satruation, class,treat,rare) %>%
    summarise(hedgesg = mean(hedges))
  library(pheatmap)
  library(tidyr)
  tempsankey$type <- paste0(tempsankey$ecosystem,tempsankey$satruation,tempsankey$treat,tempsankey$rare)
  for (ecosys in unique(tempsankey$ecosystem)) {
    tempdata <-  subset(tempsankey,tempsankey$ecosystem==ecosys)
    taxagroup <- read.csv(file =paste0("生态系统",ecosys,"相对多度排名.csv"),header=T ,row.names = 1)
    selecttaxa <- taxagroup$class[1:14]
    tempheat <- tempdata[tempdata$class %in% selecttaxa,]
    tempheat <- tempheat[,3:7]
    tempheat$class <- substr(tempheat$class,4,nchar(tempheat$class))
    
    otherdata <- tempdata[!tempdata$class %in% selecttaxa,]
    otherdata<- otherdata %>%
      group_by(ecosystem,satruation,treat,rare) %>%
      summarise(hedgesg = mean(hedgesg))
    otherdata$class <- "Others"
    otherdata$type <- paste0(otherdata$ecosystem,otherdata$satruation,otherdata$treat,otherdata$rare)
    otherdata <- otherdata[,c("class","treat","hedgesg","rare","type")]
    
    tempheat <- rbind(tempheat,otherdata)
    tempheat <- tempheat[c("class","hedgesg","type")]
    
    heatdata <- spread(tempheat,key="type",value ="hedgesg" )
    heatdata <- as.data.frame(heatdata)
    
    row_index <- which(heatdata$class== "Others",arr.ind =F)
    
    heatdata <- rbind(heatdata[-row_index, ], heatdata[row_index, ])
    
    rownames(heatdata) <- heatdata$class
    heatdata <- heatdata[,-1]
    sortcol <- c(paste0(ecosys,"10%reCRare"),paste0(ecosys,"10%reCNot"),paste0(ecosys,"40%reCRare"),paste0(ecosys,"40%reCNot"),
                  paste0(ecosys,"70%reCRare"),paste0(ecosys,"70%reCNot"),paste0(ecosys,"100%reCRare"),paste0(ecosys,"100%reCNot"),
                  paste0(ecosys,"10%reNRare"),paste0(ecosys,"10%reNNot"),paste0(ecosys,"40%reNRare"),paste0(ecosys,"40%reNNot"),
                  paste0(ecosys,"70%reNRare"),paste0(ecosys,"70%reNNot"),paste0(ecosys,"100%reNRare"),paste0(ecosys,"100%reNNot"),
                  paste0(ecosys,"10%repHRare"),paste0(ecosys,"10%repHNot"),paste0(ecosys,"40%repHRare"),paste0(ecosys,"40%repHNot"),
                  paste0(ecosys,"70%repHRare"),paste0(ecosys,"70%repHNot"),paste0(ecosys,"100%repHRare"),paste0(ecosys,"100%repHNot"),
                  paste0(ecosys,"10%reTRare"),paste0(ecosys,"10%reTNot"), paste0(ecosys,"40%reTRare"),paste0(ecosys,"40%reTNot"),
                  paste0(ecosys,"70%reTRare"),paste0(ecosys,"70%reTNot"),paste0(ecosys,"100%reTRare"),paste0(ecosys,"100%reTNot"),
                  paste0(ecosys,"10%reWRare"),paste0(ecosys,"10%reWNot"),paste0(ecosys,"40%reWRare"),paste0(ecosys,"40%reWNot"),
                  paste0(ecosys,"70%reWRare"),paste0(ecosys,"70%reWNot"),paste0(ecosys,"100%reWRare"),paste0(ecosys,"100%reWNot"))
    heatdata <- heatdata[,sortcol]
    pheatmap(
      heatdata,
      scale = "row",
      gaps_col = c(4,8,12,16),
      color = colorRampPalette(c("
      cellwidth = 10, cellheight = 10,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      border_color = "black",
      filename = paste0("生态系统_",ecosys,"_top12taxa.pdf")
    ) 
  }
  
  
  
}

{
  library(Hmisc)
  library(psych)
  library(igraph)

  setwd("path/to/your/file")
  otu16s<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  otu16s <- as.data.frame(t(otu16s))
  otu16s$names <- rownames(otu16s)
  group <- read.csv('分组.csv', header=TRUE)
  data_16s <- merge(otu16s,group,by.x ="names" ,by.y = "sample")
  rownames(data_16s) <- data_16s$mental
  data_16s <- data_16s[,c(2:44779)]
  colnames(data_16s) <- paste0(colnames(data_16s),"16s")
  metadata16s <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  rownames(metadata16s) <- paste0(rownames(metadata16s),"16s")
  
  setwd("path/to/your/file")
  otuits<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  otuits <- as.data.frame(t(otuits))
  otuits$names <- rownames(otuits)
  group <- read.csv('分组.csv', header=TRUE)
  data_its <- merge(otuits,group,by.x ="names" ,by.y = "sample")
  rownames(data_its) <- data_its$mental
  data_its <- data_its[,c(2:9518)]
  colnames(data_its) <- paste0(colnames(data_its),"its")
  metadataits <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
  rownames(metadataits) <- paste0(rownames(metadataits),"its")
  
  data_16s$names <- rownames(data_16s)
  data_its$names <- rownames(data_its)
  multidomin <- merge(data_16s,data_its,by ="names")
  rownames(multidomin) <- multidomin$names
  multidomin <- multidomin[,-1]
  colnames(metadata16s) <- colnames(metadataits)
  taxo <- rbind(metadata16s,metadataits)
  
  
  colcounter <- ncol(multidomin)+1
  multidomin$names <- rownames(multidomin)
  multidomin_cna <- merge(multidomin,group,by.x ="names" ,by.y = "mental")
  rownames(multidomin_cna) <- multidomin_cna$names
  setwd("path/to/your/file")
  resultdeprition <- c()
  for (ecosys in unique(multidomin_cna$ecosystem)) {
    tempdata <- subset(multidomin_cna,ecosystem==ecosys)
    for (satruation in unique(multidomin_cna$group)) {
      otu<-subset(tempdata,group==satruation)
      otu <- otu[,c(2:colcounter)]
      otu <- as.data.frame(t(otu))
      otu1 <- otu
      otu1[otu1>0] <- 1
      oturownum <- ncol(otu1)/2
      otu <- otu[which(rowSums(otu1) >= oturownum),]
      
      
      rowsum <- as.data.frame(rowSums(otu))
      sum <- colSums(rowsum)
      otu$RA <- rowsum/sum
      otu <- subset(otu[which(otu$RA>0.0001),], select = -RA)
      
      
      otu <- t(otu)
      otu <- scale(otu)
      
      rcorr_otu <- Hmisc::rcorr(as.matrix(otu), type = 'spearman')
      
      p <- rcorr_otu$P
      p <- p.adjust(p, method = 'BH')
      
      
      r <- rcorr_otu$r
      r[which(r < 0)]
      r[abs(r) < 0.8] <- 0
      
      
      p[p>=0.05] <- -1
      p[p<0.05 & p>=0] <- 1
      p[p==-1] <- 0
      
      
      z <- r * p
      z[which(z < 0)]
      
      
      
      g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected', diag = FALSE)
      g
      
      g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
      names(degree(g)[degree(g) == 0])
      
      E(g)$correlation <- E(g)$weight
      E(g)$weight <- abs(E(g)$weight)
      
      
      taxonomy <- taxo[as.character(V(g)$name),]
      V(g)$kingdom <- taxonomy$kingdom
      V(g)$phylum <- taxonomy$phylum
      V(g)$class <- taxonomy$class
      V(g)$order <- taxonomy$order
      V(g)$family <- taxonomy$family
      V(g)$genus <- taxonomy$genus
      V(g)$species <- taxonomy$species
      
      
      
      edge <- data.frame(as_edgelist(g))    
      
      df <- as.data.frame(E(g)$correlation)
      df[df>0] <- 1
      df[df<0] <- -1
      colnames(df) <- c('cor')
      
      edge_list <- data.frame(
        source = edge[[1]],
        target = edge[[2]],
        weight = E(g)$weight,
        correlation = E(g)$correlation,
        cor = df
      )
      
      head(edge_list)
      out_filename<-paste0(ecosys,satruation, "_edge_screen0.0001_r0.8p0.05.csv") 
      write.table(edge_list,out_filename, sep = ',', row.names = FALSE, quote = FALSE)
      
      deprition <- paste0(ecosys,satruation)
      posedgenum <- nrow(subset(edge_list,cor==1))
      negedgenum <- nrow(edge_list)-posedgenum
      deprition <-cbind(deprition,posedgenum,negedgenum)
      
      
      node <- data.frame(
        id = names(V(g)),
        kingdom =V(g)$kingdom,
        phylum =V(g)$phylum,
        class =V(g)$class,
        order =V(g)$order,
        family =V(g)$family,
        genus =V(g)$genus,
        species =V(g)$species
      )
      out_filename<-paste0(ecosys,satruation, "_node_screen0.0001_r0.8p0.05.csv") 
      write.table(node, out_filename, sep = ',', row.names = FALSE, quote = FALSE)
      
      
      nodenum <- nrow(node)
      deprition <-cbind(deprition,nodenum)
      
      
      num.edges = length(E(g)) 
      num.edges
      num.vertices = length(V(g))
      num.vertices
      connectance = edge_density(g,loops=FALSE)
      connectance
      average.degree = mean(igraph::degree(g))
      average.degree
      average.path.length = average.path.length(g) 
      average.path.length
      diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
      diameter
      edge.connectivity = edge_connectivity(g)
      edge.connectivity
      clustering.coefficient = transitivity(g) 
      clustering.coefficient
      no.clusters = no.clusters(g)
      no.clusters
      centralization.betweenness = centralization.betweenness(g)$centralization 
      centralization.betweenness
      centralization.degree = centralization.degree(g)$centralization
      centralization.degree
      
      network_property <-c(num.edges=num.edges,num.vertices=num.vertices,connectance=connectance,
                           average.degree=average.degree,average.path.length=average.path.length,
                           diameter=diameter,edge.connectivity=edge.connectivity,clustering.coefficient=clustering.coefficient,
                           no.clusters=no.clusters,centralization.betweenness=centralization.betweenness,
                           centralization.degree=centralization.degree)
      network_property <- data.frame(network_property)
      out_filename<-paste0(ecosys,satruation, "_network_property.csv") 
      write.csv(network_property, out_filename, row.names = T)
      resultdeprition <- rbind(resultdeprition,deprition)
    }
    
  }
  
  
  resultdeprition
  write.csv(resultdeprition,file = "边点数量.csv")
  
  library(igraph)
  library(ggplot2)
  library(rsq)
  library(ggpubr)
  library(ggpmisc)
  setwd("path/to/your/file")
  
  otulist = list.files(pattern="*.csv")
  otulist<-as.data.frame(otulist)
  otulist$eco <- substr(otulist$otulist,1,1)
  otulist$satru <- substr(otulist$otulist,2,nchar(otulist$otulist)-34)
  for (ecosys in unique(otulist$eco)){
    otu_eco <- subset(otulist,eco==ecosys)
    record <- c()
    record <- cbind(record,substr(otu_eco[1,1],1,9))
    print(otu_eco$satru)
    savename <- substr(otu_eco[1,1],1,nchar(otu_eco[1,1])-34)
    result <- c()
    for (satrution in unique(otu_eco$satru)) {
      otu_satru <- subset(otu_eco,satru==satrution)
      a1=read.csv( otu_satru[,1],header = TRUE) 
      a1$source=as.factor(a1$source)
      a1$target=as.factor(a1$target)
      b1=a1[,c(1:2)]
      b1=t(b1)
      b1=t(b1)
      ig=graph_from_edgelist(b1,directed=FALSE)
      natcon <- function(ig) {
        N   <- vcount(ig)
        adj <- get.adjacency(ig)
        evals <- eigen(adj)$value
        nc  <- log(mean(exp(evals)))
        nc / (N - log(N))
      }
      nc.attack <- function(ig) {
        hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=TRUE)
        sapply(1:round(vcount(ig)*.8), function(i) {
          ind <- hubord[1:i]
          tmp <- delete_vertices(ig, V(ig)$name[ind])
          natcon(tmp)
        }) }
      nc<- nc.attack(ig)
      nc
      ncoutname <- paste0("TotalNC_",ecosys,"_",satrution,".csv")
      write.csv(nc,ncoutname)
      nc<- read.csv(ncoutname)
      colnames(nc) <-c("X","NC")
      nc$RM <- nc$X/vcount(ig)
      nc$ecosystem <- ecosys
      nc$satrution <- satrution
      result <- rbind(result,nc)
    }
    formula <-  NC ~ RM 
    ncoutname <- paste0("Total-nc",ecosys,".csv")
    write.csv(result,ncoutname)
    result$satrution <- factor(result$satrution,levels=c("10%","40%","70%","100%"))
    p <- ggplot(data=result,mapping=aes( x=RM, y=NC,color = satrution))+
      geom_point(aes(color =satrution), size = 3, alpha =2) + 
      geom_smooth(aes(color =satrution),method = 'lm', formula =  y ~  x , se = T,linewidth=1.1,alpha=0.5)+
      scale_color_manual(values = c("
      stat_poly_eq(
        aes(label =  paste(..eq.label.., ..adj.rr.label.., ..p.value.label..,sep = "~`,`~")),p.digits = 2,rr.digits=3,
        formula =  y ~  x, parse = TRUE
      )+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.x = element_text(size = 20,face = "plain"))+
      theme(axis.text.y = element_text(size = 20,face = "plain"))+
      theme(axis.title.x = element_text(size = 20, face = "plain"))+
      theme(axis.title.y = element_text(size = 20, face = "plain"))+
      theme(legend.position = "top")
    out_filename <- paste0("Roboutness",ecosys,".pdf")
    ggsave(file=out_filename,p,  width = 9, height = 9)
  }
}

{
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_16S_all.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  setwd("path/to/your/file")
  otu<-read.csv("OTU_rare_ITS.csv", header=TRUE, row.names=1)
  group <- read.csv('分组.csv', header=TRUE)
  
  otu <- otu/colSums(otu)[1]
  otu <- as.data.frame(t(otu))
  colcounter <- ncol(otu)+1
  otu$names <-rownames(otu)
  input = merge(otu,group,by.x = "names",by.y ="sample" )
  rownames(input) <- input$names
  
  multi_ras<- function(colcount) {
    mean_10 <- mean(temp_10[,colcount])
    mean_40 <- mean(temp_40[,colcount])
    mean_70 <- mean(temp_70[,colcount])
    mean_100 <- mean(temp_100[,colcount])
    
    differ10_40 <- mean_10-mean_40
    differ10_70 <- mean_10-mean_70
    differ10_100 <- mean_10-mean_100
    differ40_70 <- mean_40-mean_70
    differ40_100 <- mean_40-mean_100
    differ70_100 <- mean_70-mean_100
    result_mean <- cbind(differ10_40,differ10_70,differ10_100,differ40_70,differ40_100,differ70_100)
    result_mean <- as.data.frame(t(as.data.frame(result_mean)))
    colnames(result_mean)[1] <- "differ"
    result_mean$otuname <- colnames(tempdata)[colcount]
    result_mean$ecosystem <- ecosys
    result_mean$type <- c("10_40","10_70","10_100","40_70","40_100","70_100")
    return(result_mean)
  }
  
  library(snowfall) 
  library(parallel)
  sfInit(parallel = TRUE, 12)
  finresult <- c()
  for (ecosys in unique(input$ecosystem)) {
    tempdata <- subset(input,ecosystem==ecosys)
    temp_10 <- subset(tempdata,group=="10%re")
    temp_40 <- subset(tempdata,group=="40%re")
    temp_70 <- subset(tempdata,group=="70%re") 
    temp_100 <- subset(tempdata,group=="100%re") 
    sfExport("temp_10","temp_40","temp_70","temp_100","tempdata","ecosys")
    result <- do.call(rbind, sfLapply(2:colcounter, multi_ras))
    finresult <- rbind(finresult,result)
  }
  
  write.csv(finresult, file = "16S_OTU相对多度变化.csv", row.names = F)
  write.csv(finresult, file = "itS_OTU相对多度变化.csv", row.names = F)
  
  
  library(snowfall) 
  library(parallel)
  sfInit(parallel = TRUE, 12)
  finresult <- c()
  for (ecosys in unique(input$ecosystem)) {
    tempdata <- subset(input,ecosystem==ecosys)
    tempotu <- tempdata
    tempotu[tempotu>0] <- 1
    tempotu <- as.data.frame(apply(tempotu, 2, as.numeric))
    halfrownum <- nrow(tempotu)/2
    tempdata <- tempdata[,which(colSums(tempotu) >= halfrownum)]
    temp_10 <- subset(tempdata,group=="10%re")
    temp_40 <- subset(tempdata,group=="40%re")
    temp_70 <- subset(tempdata,group=="70%re") 
    temp_100 <- subset(tempdata,group=="100%re") 
    sfExport("temp_10","temp_40","temp_70","temp_100","tempdata","ecosys")
    result <- do.call(rbind, sfLapply(2:(ncol(tempdata)-4), multi_ras))
    finresult <- rbind(finresult,result)
  }
  
  write.csv(finresult, file = "16S_一半以上_OTU相对多度变化.csv", row.names = F)
  write.csv(finresult, file = "itS_一半以上_OTU相对多度变化.csv", row.names = F)

  library(ggplot2)
  setwd("path/to/your/file")
  differ<-read.csv("16S_OTU相对多度变化.csv", header=TRUE, )
  half_differ<-read.csv("16S_一半以上_OTU相对多度变化.csv", header=TRUE,)
  
  setwd("path/to/your/file")
  differ<-read.csv("itS_OTU相对多度变化.csv", header=TRUE, )
  half_differ<-read.csv("itS_一半以上_OTU相对多度变化.csv", header=TRUE,)
  
  differ[differ$differ==0,1] <- NA
  differ<- na.omit(differ)
  half_differ[half_differ$differ==0,1] <- NA
  half_differ<- na.omit(half_differ)
  
  for (ecosys in unique(differ$ecosystem)) {
    tempdata <- subset(differ,ecosystem==ecosys)
    tempdata$type = factor(tempdata$type, levels = c("10_100","40_100","70_100","10_40","10_70","40_70"))
    p1<- ggplot(data = tempdata) + 
      geom_point(mapping = aes(x = type, y = differ,  size=5,color=differ>0,alpha= 0.001),shape=16)+
      scale_color_manual(values=c("
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 8,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text=element_text(size=8, color="black"))+
      theme(axis.text.x = element_text(size = 8,color="black",angle = 90)) 
    outputname <- paste0("相对多度变化_",ecosys,"_ITS.pdf")
    ggsave(outputname,p1, width = 9, height =6, useDingbats = FALSE,units = "cm")
  }
  
  for (ecosys in unique(half_differ$ecosystem)) {
    tempdata <- subset(half_differ,ecosystem==ecosys)
    tempdata$type = factor(tempdata$type, levels = c("10_100","40_100","70_100","10_40","10_70","40_70"))
    p1<- ggplot(data = tempdata) + 
      geom_point(mapping = aes(x = type, y = differ,  size=5,color=differ>0,alpha= 0.001),shape=16)+
      scale_color_manual(values=c("
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 8,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text=element_text(size=8, color="black"))+
      theme(axis.text.x = element_text(size = 8,color="black",angle = 90)) 
    outputname <- paste0("相对多度变化_一半以上_",ecosys,"_ITS.pdf")
    ggsave(outputname,p1, width = 9, height = 6, useDingbats = FALSE,units = "cm")
  }
}  
