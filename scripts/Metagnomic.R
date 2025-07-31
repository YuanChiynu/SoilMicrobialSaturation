{
  setwd("path/to/your/file")
  quant_data <- read.delim("KEGG_count_add_addNif.txt",header = F,check.names = F,sep = " ")
  colnames(quant_data)[2] = "KO"
  length_data <- read.delim("ko_length_addNif.txt",sep = " ",header = F)
  colnames(length_data) <- c("ID","gene ID","KO","Length")
  
  library(tidyr)
  c = colnames(quant_data)[3:ncol(quant_data)]
  data <- gather(quant_data,sample,value,c)
  
  length_data1 = subset(length_data,select = c(ID,Length))
  data1 = merge(data,length_data1,by.x = "
  data1$v1 = data1$Length*3
  
  data1$Lseq = 150
  N16Sdata <- read.delim("metadata.txt",header = TRUE,check.names = F)
  colnames(N16Sdata)[1] = "Name"
  N16Sdata = subset(N16Sdata,select = c(Name,n16S))
  data1 = merge(data1,N16Sdata,by.x = "sample",by.y = "Name")
  data1$L16s = 1432
  attach(data1)
  data1$value = as.numeric(data1$value)
  write.csv(data1,"temp.csv")
  data1 <- read.csv("temp.csv",header = T,row.names = 1)
  data1$RA <-  ( data1$value* data1$Lseq/ data1$v1)/( data1$n16S* data1$Lseq/ data1$L16s)
  detach(data1)
  
  qpcr = read.csv("copynum.csv",header = T)
  data1 = merge(data1,qpcr,by.x = "sample",by.y = "sample")
  
  
  data1$gene_copies = data1$RA*data1$copies
  data2 = subset(data1,select = c(sample,KO,gene_copies))
  
  write.csv(data1,"gene_copies.csv")
  write.csv(data2,"gene_copies2.csv")
  
  
  data3 = subset(data1,select = c(sample,KO,value,RA))
  write.csv(data3,"reads_ra.csv")
  
}

{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  
  carbonKO = read.csv("新KO.csv")
  carbonKO = carbonKO[carbonKO$cycle == "碳循环",]
  vcarbonKO = carbonKO$KEGG_KO
  
  data = genecopies_sum[genecopies_sum$ko %in% vcarbonKO,]
  
  datako = unique(data$ko)
  datagene = carbonKO[carbonKO$KEGG_KO %in% datako,]
  
  sum = aggregate(data$gene_copies,by = list(data$sample),FUN = "sum")
  colnames(sum) = c("sample","gene_copies")
  input = merge(sum,group,by = "sample")
  
  input$value <- log10(input$gene_copies)
  
  input$satru <- factor(input$satru,levels=c("10%","100%"))
  input$treat = factor(input$treat, levels = c("CK","C","N","pH","T","W"))
  
  library(ggplot2)
  library(ggpubr)
  p_sub <- ggbarplot(input, x="satru", y="value", add = "mean_se",color = "black",fill = "treat",
                     palette =c("
    ylab("碳循环基因拷贝数")+
    xlab("")+
    ylim(0,(max(input$value)*1.2))+
    stat_compare_means(aes(group=treat), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
  p_sub
  ggsave("总碳循环基因拷贝数.pdf",p_sub, width = 6, height = 6, useDingbats = FALSE)
}

{
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  
  nitroKO = read.csv("新KO.csv")
  nitroKO = nitroKO[nitroKO$cycle == "氮循环",]
  vnitroKO = nitroKO$KEGG_KO
  
  data = genecopies_sum[genecopies_sum$ko %in% vnitroKO,]
  
  datako = unique(data$ko)
  datagene = nitroKO[nitroKO$KEGG_KO %in% datako,]
  
  sum = aggregate(data$gene_copies,by = list(data$sample),FUN = "sum")
  colnames(sum) = c("sample","gene_copies")
  input = merge(sum,group,by = "sample")
  
  input$value <- log10(input$gene_copies)
  
  input$satru <- factor(input$satru,levels=c("10%","100%"))
  input$treat = factor(input$treat, levels = c("CK","C","N","pH","T","W"))
  
  library(ggplot2)
  library(ggpubr)
  p_sub <- ggbarplot(input, x="satru", y="value", add = "mean_se",color = "black",fill = "treat",
                     palette =c("
    ylab("氮循环基因拷贝数")+
    xlab("")+
    ylim(0,(max(input$value)*1.2))+
    stat_compare_means(aes(group=treat), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
  p_sub
  ggsave("总氮循环基因拷贝数.pdf",p_sub, width = 6, height = 6, useDingbats = FALSE)
}

{
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  
  phospKO = read.csv("新KO.csv")
  phospKO = phospKO[phospKO$cycle == "磷循环",]
  vphospKO = phospKO$KEGG_KO
  
  data = genecopies_sum[genecopies_sum$ko %in% vphospKO,]
  
  datako = unique(data$ko)
  datagene = phospKO[phospKO$KEGG_KO %in% datako,]
  
  sum = aggregate(data$gene_copies,by = list(data$sample),FUN = "sum")
  colnames(sum) = c("sample","gene_copies")
  input = merge(sum,group,by = "sample")
  
  input$value <- log10(input$gene_copies)
  
  input$satru <- factor(input$satru,levels=c("10%","100%"))
  input$treat = factor(input$treat, levels = c("CK","C","N","pH","T","W"))
  
  library(ggplot2)
  library(ggpubr)
  p_sub <- ggbarplot(input, x="satru", y="value", add = "mean_se",color = "black",fill = "treat",
                     palette =c("
    ylab("磷循环基因拷贝数")+
    xlab("")+
    ylim(0,(max(input$value)*1.2))+
    stat_compare_means(aes(group=treat), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
  p_sub
  ggsave("总磷循环基因拷贝数.pdf",p_sub, width = 6, height = 6, useDingbats = FALSE)
}

{
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  
  sulfiKO = read.csv("新KO.csv")
  sulfiKO = sulfiKO[sulfiKO$cycle == "硫循环",]
  vsulfiKO = sulfiKO$KEGG_KO
  
  data = genecopies_sum[genecopies_sum$ko %in% vsulfiKO,]
  
  datako = unique(data$ko)
  datagene = sulfiKO[sulfiKO$KEGG_KO %in% datako,]
  
  sum = aggregate(data$gene_copies,by = list(data$sample),FUN = "sum")
  colnames(sum) = c("sample","gene_copies")
  input = merge(sum,group,by = "sample")
  
  input$value <- log10(input$gene_copies)
  
  input$satru <- factor(input$satru,levels=c("10%","100%"))
  input$treat = factor(input$treat, levels = c("CK","C","N","pH","T","W"))
  
  library(ggplot2)
  library(ggpubr)
  p_sub <- ggbarplot(input, x="satru", y="value", add = "mean_se",color = "black",fill = "treat",
                     palette =c("
    ylab("硫循环基因拷贝数")+
    xlab("")+
    ylim(0,(max(input$value)*1.2))+
    stat_compare_means(aes(group=treat), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
    theme(axis.text.y  = element_text(size = 20,color="black"))+
    theme(axis.title.y = element_text(size = 20, color="black"))+
    theme(axis.title.x = element_text(size = 20, color="black"))+
    theme(legend.text=element_text(size=20, color="black"))+
    theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
  p_sub
  ggsave("总硫循环基因拷贝数.pdf",p_sub, width = 6, height = 6, useDingbats = FALSE)
}

{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  
  carbonKO = read.csv("新KO.csv")
  carbonKO = carbonKO[carbonKO$cycle == "碳循环",]
  vcarbonKO = carbonKO$KEGG_KO
  data = genecopies_sum[genecopies_sum$ko %in% vcarbonKO,]
  data = merge(data,group,by = "sample")
  data <- subset(data,gene_copies>0)
  input = merge(data,carbonKO,by.x = "ko",by.y ="KEGG_KO" )
  
  input$value <- log10(input$gene_copies)
  for (cfunc in unique(input$class)) {
    par(mfrow=c(1,3)) 
    tempdata <- subset(input,class==cfunc)
    library(ggplot2)
    library(ggpubr)
    sort_index<- sort(unique(tempdata$Gene_name))
    tempdata$Gene_name = factor(tempdata$Gene_name, levels = sort_index)

    p_sub <- ggbarplot(tempdata, x="Gene_name", y="value", add = "mean_se",color = "black",fill = "satru",
                       palette =c("
      ylab(cfunc)+
      xlab("")+
      ylim(0,(max(input$value)*1.2))+
      stat_compare_means(aes(group=satru), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(legend.text=element_text(size=20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
    p_sub 
    savename <- paste0("碳循环",cfunc,".pdf")
    
    ggsave(savename,p_sub, width = 18, height = 6, useDingbats = FALSE)
    }
  }
  
{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  nitroKO = read.csv("新KO.csv")
  nitroKO = nitroKO[nitroKO$cycle == "氮循环",]
  vnitroKO = nitroKO$KEGG_KO
  data = genecopies_sum[genecopies_sum$ko %in% vnitroKO,]
  data = merge(data,group,by = "sample")
  data <- subset(data,gene_copies>0)
  input = merge(data,nitroKO,by.x = "ko",by.y ="KEGG_KO" )
  
  input$value <- log10(input$gene_copies)
  for (cfunc in unique(input$class)) {
    par(mfrow=c(1,3)) 
    tempdata <- subset(input,class==cfunc)
    library(ggplot2)
    library(ggpubr)
    sort_index<- sort(unique(tempdata$Gene_name))
    tempdata$Gene_name = factor(tempdata$Gene_name, levels = sort_index)
    
    p_sub <- ggbarplot(tempdata, x="Gene_name", y="value", add = "mean_se",color = "black",fill = "satru",
                       palette =c("
      ylab(cfunc)+
      xlab("")+
      ylim(0,(max(input$value)*1.2))+
      stat_compare_means(aes(group=satru), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(legend.text=element_text(size=20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
    p_sub 
    savename <- paste0("氮循环",cfunc,".pdf")
    if (cfunc=="反硝化") {
      w <- 12
    }else{
      w <- 6
    }
    ggsave(savename,p_sub, width = w, height = 6, useDingbats = FALSE)
  }
}

{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  phospKO = read.csv("新KO.csv")
  phospKO = phospKO[phospKO$cycle == "磷循环",]
  vphospKO = phospKO$KEGG_KO
  data = genecopies_sum[genecopies_sum$ko %in% vphospKO,]
  data = merge(data,group,by = "sample")
  data <- subset(data,gene_copies>0)
  input = merge(data,phospKO,by.x = "ko",by.y ="KEGG_KO" )
  
  input$value <- log10(input$gene_copies)
  for (cfunc in unique(input$class)) {
    par(mfrow=c(1,3)) 
    tempdata <- subset(input,class==cfunc)
    library(ggplot2)
    library(ggpubr)
    sort_index<- sort(unique(tempdata$Gene_name))
    tempdata$Gene_name = factor(tempdata$Gene_name, levels = sort_index)
    
    p_sub <- ggbarplot(tempdata, x="Gene_name", y="value", add = "mean_se",color = "black",fill = "satru",
                       palette =c("
      ylab(cfunc)+
      xlab("")+
      ylim(0,(max(input$value)*1.2))+
      stat_compare_means(aes(group=satru), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(legend.text=element_text(size=20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
    p_sub 
    savename <- paste0("磷循环",cfunc,".pdf")
    
    ggsave(savename,p_sub, width = 9, height = 6, useDingbats = FALSE)
  }
}

{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  sulfiKO = read.csv("新KO.csv")
  sulfiKO = sulfiKO[sulfiKO$cycle == "硫循环",]
  vsulfiKO = sulfiKO$KEGG_KO
  data = genecopies_sum[genecopies_sum$ko %in% vsulfiKO,]
  data = merge(data,group,by = "sample")
  data <- subset(data,gene_copies>0)
  input = merge(data,sulfiKO,by.x = "ko",by.y ="KEGG_KO" )
  
  input$value <- log10(input$gene_copies)
  for (cfunc in unique(input$class)) {
    par(mfrow=c(1,3)) 
    tempdata <- subset(input,class==cfunc)
    library(ggplot2)
    library(ggpubr)
    sort_index<- sort(unique(tempdata$Gene_name))
    tempdata$Gene_name = factor(tempdata$Gene_name, levels = sort_index)
    
    p_sub <- ggbarplot(tempdata, x="Gene_name", y="value", add = "mean_se",color = "black",fill = "satru",
                       palette =c("
      ylab(cfunc)+
      xlab("")+
      ylim(0,(max(input$value)*1.2))+
      stat_compare_means(aes(group=satru), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(legend.text=element_text(size=20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
    p_sub 
    savename <- paste0("硫循环",cfunc,".pdf")
    
    ggsave(savename,p_sub, width = 9, height = 6, useDingbats = FALSE)
  }
}

{
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  KO = read.csv("新KO.csv")
  
  
  reresult <- c()
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    
    carbonKO = KO[KO$cycle == cyc,]
    vcarbonKO = carbonKO$KEGG_KO
    
    data = genecopies_sum[genecopies_sum$ko %in% vcarbonKO,]
    
    datako = unique(data$ko)
    datagene = carbonKO[carbonKO$KEGG_KO %in% datako,]
    
    sum = aggregate(data$gene_copies,by = list(data$sample),FUN = "sum")
    colnames(sum) = c("sample","gene_copies")
    input = merge(sum,group,by = "sample")
    
    input$value <- log10(input$gene_copies)
    for (sat in unique(input$satru)){
      sub_phi <- subset(input,satru==sat)
      mod1 = aov( value~treat, data= sub_phi)
      print(paste0(cyc,sat))
      print(summary(mod1))
      re = LSD.test(mod1,"treat",alpha = 0.05)
      re1 = re$groups
      re1$cycle <- cyc
      re1$satruation <- sat
      re1$treat <-rownames(re1) 
      reresult <- rbind(reresult,re1)
    }
  }
  write.csv(reresult,file = "总宏基因组多重比较.csv")
}

{
  library(ggplot2)
  library(ggpubr)
  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  KO = read.csv("新KO.csv")


  for (cycle in c("碳循环","氮循环","磷循环","硫循环")) {
    subKO = KO[KO$cycle == cycle,]
    vsubKO = subKO$KEGG_KO
    
    data = genecopies_sum[genecopies_sum$ko %in% vsubKO,]
    data = merge(data,group,by = "sample")
    data <- subset(data,gene_copies>0)
    input = merge(data,subKO,by.x = "ko",by.y ="KEGG_KO" )
    input$value <- log10(input$gene_copies)
    for(cfunc in unique(input$class)){
      pretemp <- subset(input,class==cfunc)
      for (gene in  unique(input$Gene_name)) {
        tempdata <- subset(pretemp,Gene_name==gene)
        tempdata$satru <- factor(tempdata$satru,levels=c("10%","100%"))
        tempdata$treat = factor(tempdata$treat, levels = c("CK","C","N","pH","T","W"))
        
        p_sub <- ggbarplot(tempdata, x="satru", y="value", add = "mean_se",color = "black",fill = "treat",
                           palette =c("
          ylab(gene)+
          xlab("")+
          ylim(0,(max(input$value)*1.2))+
          stat_compare_means(aes(group=treat), label = "p.signif", label.y = (max(input$value)), size = 10,method ="anova",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
          theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
          theme(axis.text.y  = element_text(size = 20,color="black"))+
          theme(axis.title.y = element_text(size = 20, color="black"))+
          theme(axis.title.x = element_text(size = 20, color="black"))+
          theme(legend.text=element_text(size=20, color="black"))+
          theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 
        p_sub 
        savename <- paste0(cycle,"_",cfunc,"_",gene, ".pdf")
        
        ggsave(savename,p_sub, width = 6, height = 6, useDingbats = FALSE)
      }
    }
  }
  
}

{
  library(pheatmap)
  library(tidyr)
  
  setwd("path/to/your/file")
  genecopies = read.csv("reads_ra.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("新KO.csv")
  
  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  rownames(heat_matrix) <- heat_matrix[,1]
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name
  
  cyc <- c("碳循环","氮循环","磷循环","硫循环")
  tempheat <- subset(input,input$cycle=="碳循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
  heat_gap_c <- tempheat[,c(49:53)]
  heat_data <- tempheat[,c(2:49)]
  sorted_columns <- group$sample
  heat_data <- heat_data[, sorted_columns]

  pheatmap(
    heat_data,
    scale = "row",
    cellwidth = 10, cellheight = 10,
    gaps_row = c(3, 23),
    gaps_col = 24,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "black",
    filename = "碳循环.pdf"
  ) 
  tempheat <- subset(input,input$cycle=="氮循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
  heat_gap_n <- tempheat[,c(49:53)]
  heat_data <- tempheat[,c(2:49)]
  sorted_columns <- group$sample
  heat_data <- heat_data[, sorted_columns]
  heat_data <- heat_data[-20,]

  pheatmap(
    heat_data,
    scale = "row",
    cellwidth = 10, cellheight = 10,
    gaps_row = c(1,5,12,13,14,18,21),
    gaps_col = 24,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "black",
    filename = "氮循环.pdf"
  ) 
  tempheat <- subset(input,input$cycle=="磷循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
  heat_gap_p <- tempheat[,c(49:53)]
  heat_data <- tempheat[,c(2:49)]
  sorted_columns <- group$sample
  heat_data <- heat_data[, sorted_columns]

  
  pheatmap(
    heat_data,
    scale = "row",
    cellwidth = 10, cellheight = 10,
    gaps_row = 4,
    gaps_col = 24,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "black",
    filename = "磷循环.pdf"
  ) 
  tempheat <- subset(input,input$cycle=="硫循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
  heat_gap_s <- tempheat[,c(49:53)]
  heat_data <- tempheat[,c(2:49)]
  sorted_columns <-group$sample
  heat_data <- heat_data[, sorted_columns]

  
  pheatmap(
    heat_data,
    scale = "row",
    cellwidth = 10, cellheight = 10,
    gaps_row = c(5,6),
    gaps_col = 24,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "black",
    filename = "硫循环.pdf"
  ) 
}

{
  
  setwd("path/to/your/file")
  genecopies = read.csv("reads_ra.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("新KO.csv")
  
  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "mean")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  rownames(heat_matrix) <- heat_matrix[,1]
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name
  resultslis <- c()
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
    heat_gap_s <- tempheat[,c(49:53)]
    heat_data <- tempheat[,c(2:49)]
    sorted_columns <-group$sample
    heat_data <- heat_data[, sorted_columns]
    for (row in 1:nrow(heat_data)) {
      genecv_10 <-heat_data[row,c(1:24)]
      genecv_100 <-heat_data[row,c(25:48)]
      sd_10 <- sd(genecv_10)
      mean_10 <- rowMeans(genecv_10)
      cv_10 <- (sd_10 / mean_10) * 100
      resultcatcher <- c(rownames(heat_data)[row],"10%",cv_10,heat_gap_s$cycle[row],heat_gap_s$class[row])
      resultslis <- rbind(resultslis,resultcatcher)
      sd_100 <- sd(genecv_100)
      mean_100 <- rowMeans(genecv_100)
      cv_100 <- (sd_100 / mean_100) * 100
      resultcatcher <- c(rownames(heat_data)[row],"100%",cv_100,heat_gap_s$cycle[row],heat_gap_s$class[row])
      resultslis <- rbind(resultslis,resultcatcher)
    }
  }
  write.csv(resultslis,file = "不同基因分饱和度cv.csv")
  resultslis <- as.data.frame(resultslis)
  resultslis[,3] <- as.numeric(resultslis[,3] )
  
  
  
  library(dplyr)
  setwd("path/to/your/file")
  resultslis = read.csv("不同基因分饱和度cv.csv")
  resultslis <- resultslis[,-1]
  colnames(resultslis) <- c("gene_name","satru","genecv","cycle","class")
  library(ggplot2)
  library(ggpubr)

  resultslis <- na.omit(resultslis)
  my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    input <- subset(resultslis,cycle==cyc)
    input$satru <- factor(input$satru,levels = c("100%","10%"))
    p_sub <- ggbarplot(input, x="gene_name", y="genecv",fill = "satru",
                       palette =c("
      xlab("")+
      ylab("")+
      coord_flip()+
      ylim(0,(max(input$genecv)*1.2))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 8,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text=element_text(size=8, color="black"))+
      theme(legend.position ="none")+
      theme(axis.text.x = element_text(size = 8,color="black",vjust = 1.0, hjust = 1.0))
      my_theme
      
    p_sub
    ggsave(paste0(cyc,"基因cv.pdf"),p_sub,width = 6,height = length(unique(input$gene_name))/2,units = "cm" )
  }
}

{
  library(dplyr)
  library(tidyr)
  setwd("path/to/your/file")
  genecopies = read.csv("reads_ra.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("新KO.csv")
  
  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "mean")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  rownames(heat_matrix) <- heat_matrix[,1]
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name
  resultslis <- c()
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
    heat_gap_s <- tempheat[,c(49:53)]
    heat_data <- tempheat[,c(2:49)]
    
    for (row in 1:nrow(heat_data)) {
      gene <- heat_data[row,]
      gene <- as.data.frame(t(gene))
      gene$name <- rownames(gene)
      genecv <-  merge(gene,group,by.x = "name",by.y ="sample" )
      for (treat in c("C","N","pH","T","W")) {
        select_treat <- c("CK",treat)
        genecv_fin <- genecv[genecv_fin$treat %in% c("CK",treat),]
        genecv_fin_10 <- genecv_fin[genecv_fin$satru %in% c("10%"),]
        genecv_fin_100 <-  genecv_fin[genecv_fin$satru %in% c("100%"),]
        sd_10 <- sd(genecv_fin_10[,2])
        mean_10 <- mean(genecv_fin_10[,2])
        cv_10 <- (sd_10 / mean_10) * 100
        resultcatcher <- c(colnames(gene)[1],"10%",cv_10,heat_gap_s$cycle[row],heat_gap_s$class[row],treat)
        resultslis <- rbind(resultslis,resultcatcher)
        
        sd_100 <- sd(genecv_fin_100[,2])
        mean_100 <- mean(genecv_fin_10[,2])
        cv_100 <- (sd_100 / mean_100) * 100
        resultcatcher <- c(colnames(gene)[1],"100%",cv_100,heat_gap_s$cycle[row],heat_gap_s$class[row],treat)
        resultslis <- rbind(resultslis,resultcatcher)
        
      }
    }
  }
  resultslis <- as.data.frame(resultslis)
  resultslis <- na.omit(resultslis)
  resultslis[,3] <- as.numeric(resultslis[,3] )
  
  write.csv(resultslis,file = "不同基因ck与处理cv.csv")
  resultslis <- read.csv(file ="不同基因ck与处理cv.csv",header = T,row.names = 1 )
 
  library(ggplot2)
  library(ggpubr)
  colnames(resultslis) <- c("gene_name","satru","genecv","cycle","class","treat")

  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempdata_cyc <- subset(resultslis,cycle==cyc)
    for (type in unique(tempdata_cyc$class)) {
      tempdata_class <- subset(tempdata_cyc,class==type)
      for (gene in unique(tempdata_class$gene_name)) {
        data <-  subset(tempdata_class,gene_name==gene)
        data$satru=factor(data$satru,levels=c("10%","100%"))
        p_sub <- ggbarplot(data, x="satru", y="genecv", add = "mean_se",color = "black",fill = "satru",
                           palette =c("
          ylab("CV")+
          xlab(gene)+
          ylim(0,(max(data$genecv)*1.2))+
          stat_compare_means(aes(group=satru),label = "p.signif", label.y = (max(data$genecv)),label.x.npc = "middle", size = 7,method ="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, Inf), symbols = c("***", "**", "*", "ns")))+
          theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
          theme(axis.text.y  = element_text(size = 10,color="black"))+
          theme(axis.title.y = element_text(size = 10, color="black"))+
          theme(axis.title.x = element_text(size = 10, color="black"))+
          theme(legend.text=element_text(size=10, color="black"))+
          theme(axis.text.x = element_text(size = 10,color="black",vjust = 1.0, hjust = 1.0))
        
        p_sub
        ggsave(paste0(cyc,"_",type,"_",gene,"cv.pdf"),p_sub, width = 2, height = 3, useDingbats = FALSE)
      }
    }
  }
}

{
  library(esc)
  library(tidyr)
  

  setwd("path/to/your/file")
  genecopies = read.csv("gene_copies2.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("KO.csv")
  
  genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  rownames(heat_matrix) <- heat_matrix[,1]
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name
  resultslis <- c()
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
    heat_gap_s <- tempheat[,c(49:53)]
    heat_data <- tempheat[,c(2:49)]
    
    for (row in 1:nrow(heat_data)) {
      gene <- heat_data[row,]
      gene <- as.data.frame(t(gene))
      gene$name <- rownames(gene)
      genecv <-  merge(gene,group,by.x = "name",by.y ="sample" )
      for (treat in c("C","N","pH","T","W")) {
        select_treat <- c("CK",treat)
        genecv_fin <- genecv[genecv$treat %in% select_treat,]
        genecv_fin_10 <- genecv_fin[genecv_fin$satru %in% c("10%"),]
        genecv_fin_100 <-  genecv_fin[genecv_fin$satru %in% c("100%"),]
        
        sd_10_ck <- sd(genecv_fin_10[genecv_fin_10$treat %in% "CK",2])
        mean_10_ck <-  mean(genecv_fin_10[genecv_fin_10$treat %in% "CK",2])
        
        sd_10_treat <- sd(genecv_fin_10[genecv_fin_10$treat %in% treat,2])
        mean_10_treat <-  mean(genecv_fin_10[genecv_fin_10$treat %in% treat,2])

        hedges_g_10<- esc_mean_sd(grp1m = sd_10_ck, grp1sd =mean_10_ck, grp1n = 4,grp2m = sd_10_treat, grp2sd =mean_10_treat, grp2n = 4, es.type = "g")
        resultcatcher <- c(colnames(gene)[1],"10%",hedges_g_10$es,heat_gap_s$cycle[row],heat_gap_s$class[row],treat)
        resultslis <- rbind(resultslis,resultcatcher)
        
        
        sd_100_ck <- sd(genecv_fin_100[genecv_fin_100$treat %in% "CK",2])
        mean_100_ck <-  mean(genecv_fin_100[genecv_fin_100$treat %in% "CK",2])
        
        sd_100_treat <- sd(genecv_fin_100[genecv_fin_100$treat %in% treat,2])
        mean_100_treat <-  mean(genecv_fin_100[genecv_fin_100$treat %in% treat,2])
        
        hedges_g_100<- esc_mean_sd(grp1m = sd_100_ck, grp1sd =mean_100_ck, grp1n = 4,grp2m = sd_100_treat, grp2sd =mean_100_treat, grp2n = 4, es.type = "g")
        resultcatcher <- c(colnames(gene)[1],"100%",hedges_g_100$es,heat_gap_s$cycle[row],heat_gap_s$class[row],treat)
        resultslis <- rbind(resultslis,resultcatcher)

      }
    }
  }

  colnames(resultslis) <- c("Gene_ID","satru","es","cyc","func","treat")
  resultslis <- as.data.frame(resultslis)
  resultslis$type<- paste0(resultslis$treat,resultslis$satru)
  write.csv(resultslis,file = "不同基因分饱和度、处理效应值hedges_g.csv")
  
  
  library(tidyr)
  setwd("path/to/your/file")
  data <-read.csv(file = "不同基因分饱和度、处理效应值hedges_g.csv",header = T,row.names = 1)
  data <- data[,c("Gene_ID","es","type")]
  gmatrix<- spread (data,key="type",value ="es",fill = 0)
  gmatrix <- na.omit(gmatrix)
  library(pheatmap)
  library(tidyr)
  KO = read.csv("KO.csv")
  input = merge(gmatrix,KO,by.x = "Gene_ID",by.y ="Gene_name" )
  rownames(input) <- input$Gene_ID
  
  tempheat <- subset(input,input$cycle=="碳循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_ID),]
  heat_gap_c <- tempheat[,c(1,12:15)]
  heat_data <- tempheat[,c(2:11)]

  
  pheatmap(
    heat_data,
    scale = "none",
    cellwidth = 20, cellheight = 20,
    color = colorRampPalette(c("
    gaps_row = c(3,23),
    gaps_col =c(2,4,6,8) ,
    cluster_rows = F,
    cluster_cols = F,
    display_numbers = TRUE,  
    number_format = "%.2f",  
    number_color = "black",  
    fontsize_number = 7,    
    border_color = "black",
    filename = "碳基因效应值.pdf"
  ) 
 
  tempheat <- subset(input,input$cycle=="氮循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_ID),]
  heat_gap_c <- tempheat[,c(1,12:15)]
  heat_data <- tempheat[,c(2:11)]


  
  pheatmap(
    heat_data,
    scale = "none",
    cellwidth = 20, cellheight = 20,
    color = colorRampPalette(c("
    gaps_row = c(1,5,12,13,14,18,21),
    gaps_col = c(2,4,6,8),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,  
    number_format = "%.2f",  
    number_color = "black",  
    fontsize_number = 7,    
    border_color = "black",
    filename = "氮基因效应值.pdf"
  ) 
  tempheat <- subset(input,input$cycle=="磷循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_ID),]
  heat_gap_c <- tempheat[,c(1,12:15)]
  heat_data <- tempheat[,c(2:11)]

  
  pheatmap(
    heat_data,
    scale = "none",
    cellwidth = 20, cellheight = 20,
    color = colorRampPalette(c("
    gaps_row = 4,
    gaps_col = c(2,4,6,8),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,  
    number_format = "%.2f",  
    number_color = "black",  
    fontsize_number = 7,    
    border_color = "black",
    filename = "磷基因效应值.pdf"
  ) 
  tempheat <- subset(input,input$cycle=="硫循环")
  tempheat <- tempheat[order(tempheat$class,tempheat$Gene_ID),]
  heat_gap_c <- tempheat[,c(1,12:15)]
  heat_data <- tempheat[,c(2:11)]

  
  
  pheatmap(
    heat_data,
    scale = "none",
    cellwidth = 20, cellheight = 20,
    color = colorRampPalette(c("
    gaps_row = c(5,6),
    gaps_col =  c(2,4,6,8),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,  
    number_format = "%.2f",  
    number_color = "black",  
    fontsize_number = 7,    
    border_color = "black",
    filename = "硫基因效应值.pdf"
  ) 

}

{
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(agricolae)
  library(tidyr)
  setwd("path/to/your/file")
  genecopies = read.csv("reads_ra.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("KO.csv")
  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)

  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name


  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    heat_data <- tempheat[,c(2:49)]
    heat_data <- as.data.frame(t(heat_data))
    nmds1 <- metaMDS(heat_data, autotransform =FALSE)
    nmds1.stress <- nmds1$stress
    nmds1.point <- data.frame(nmds1$point)
    nmds1.species <- data.frame(nmds1$species)
    sample_site <- nmds1.point[1:2]
    sample_site$names <- rownames(sample_site)
    names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
    sample_site <- merge(sample_site, group, by.x = 'names',  by.y = 'sample',all.x = TRUE)
    sample_site$satru <- factor(sample_site$satru,levels=c("10%","100%"))
    p <- ggplot(sample_site, aes(x=NMDS1, y=NMDS2, group = satru)) +
      geom_point(aes(color = satru, shape = treat), size = 4, alpha =2) + 
      scale_shape_manual(values = c(15,3,16,17,18,4)) + 
      scale_color_manual(values = c("
      stat_ellipse(data=sample_site,
                   geom = "polygon",level=0.95,
                   linetype = 1,size=0.8,
                   aes(fill=satru),
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
      theme(legend.position = "none")
    p 

    ggsave(p, file=paste0("新_",cyc,"_相对多度_NMDS.pdf"), width = 6, height = 6)
  }
  group1 <- group
  colnames(group1) <-c("sample","satru","treat" ) 
  group2 <- group
  colnames(group2) <-c("sample2","satru2","treat2" ) 
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    heat_data <- tempheat[,c(2:49)]
    heat_data <- as.data.frame(t(heat_data))
    prok.dist <- vegdist(heat_data, method="bray")
    prok.dist <- as.matrix(prok.dist)
    prok.dist <- as.data.frame(prok.dist)
    prok.dist$X1 <- rownames(prok.dist)
    prok.matr <- gather(prok.dist, key = "X2", value = "dissimi",-X1)
  
    colnames(prok.matr) <-c("sample1","sample2","similarity" ) 
    prok.matr$similarity <- 1-prok.matr$similarity
    simi_dis <- merge(prok.matr, group1, by.x = 'sample1',by.y = 'sample', all.x = TRUE)
    simi_dis <- merge(simi_dis, group2, by.x = 'sample2',by.y = 'sample2', all.x = TRUE)
    fin_dis_zhu<-subset(simi_dis,satru==satru2)
    fin_dis_zhu<-subset(fin_dis_zhu,treat=="CK")
    fin_dis_zhu<-subset(fin_dis_zhu,treat2!="CK")
    
    results_summary <- fin_dis_zhu %>%
      group_by(satru) %>%  
      summarise(
        mean_similarity = mean(similarity, na.rm = TRUE),  
        sd_similarity   = sd(similarity, na.rm = TRUE)      
      ) %>%
      ungroup()  
    
    browser()
    print(cyc)
    print(results_summary)}
    
    p_16s <- ggbarplot(fin_dis_zhu, x="satru", y="similarity", add = "mean_se",color = "black",fill = "satru",
                       palette = c("
      ylab("Similarity")+
      xlab("Satruration")+
      ylim(0,1)+
      stat_compare_means(aes(group=satru), label = "p.signif", label.y = 0.95,label.x.npc = "middle" , size = 10,method = "t.test")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+ 
      theme(axis.text.y  = element_text(size = 20,color="black"))+
      theme(axis.title.y = element_text(size = 20, color="black"))+
      theme(axis.title.x = element_text(size = 20, color="black"))+
      theme(legend.text=element_text(size=20, color="black"))+
      theme(axis.text.x = element_text(size = 20,color="black",vjust = 1.0, hjust = 1.0)) 

    p_16s
    ggsave(filename = paste0(cyc,'similarity_16s.pdf'), p_16s, width = 3, height = 6)
  }

}

{
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(agricolae)
  library(tidyr)
  setwd("path/to/your/file")
  genecopies = read.csv("reads_ra.csv",row.names = 1)
  group = read.csv("group.csv")
  KO = read.csv("KO.csv")
  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")
  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )
  rownames(input) <- input$Gene_name
  
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    heat_data <- tempheat[,c(2:49)]
    heat_data <- as.data.frame(t(heat_data))
    group_permanova <- group
    rownames(group_permanova) <- group_permanova[,1]
    group_permanova <- group_permanova[,-1]
    results1<- adonis2(formula = heat_data ~ satru+treat,data = group, permutations = 999, method = "bray")
    print(cyc)
    print(results1)
  }
  group1 <- group
  colnames(group1) <-c("sample","satru","treat" ) 
  group2 <- group
  colnames(group2) <-c("sample2","satru2","treat2" ) 
  for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
    tempheat <- subset(input,input$cycle==cyc)
    heat_data <- tempheat[,c(2:49)]
    heat_data <- as.data.frame(t(heat_data))
    prok.dist <- vegdist(heat_data)
    prok.dist <- as.matrix(prok.dist)
    prok.dist <- as.data.frame(prok.dist)
    prok.dist$X1 <- rownames(prok.dist)
    prok.matr <- gather(prok.dist, key = "X2", value = "dissimi",-X1)
    
    colnames(prok.matr) <-c("sample1","sample2","similarity" ) 
    prok.matr$similarity <- 1-prok.matr$similarity
    simi_dis <- merge(prok.matr, group1, by.x = 'sample1',by.y = 'sample', all.x = TRUE)
    simi_dis <- merge(simi_dis, group2, by.x = 'sample2',by.y = 'sample2', all.x = TRUE)
    fin_dis_zhu<-subset(simi_dis,satru==satru2)
    fin_dis_zhu<-subset(fin_dis_zhu,treat=="CK")
    fin_dis_zhu<-subset(fin_dis_zhu,treat2!="CK")
    fin_dis_1<-subset(fin_dis_zhu,satru=="10%")
    fin_dis_1 <- fin_dis_1$similarity
    fin_dis_2<-subset(fin_dis_zhu,satru=="100%")
    fin_dis_2 <- fin_dis_2$similarity
    mod1 <- t.test( x=fin_dis_2 ,y=fin_dis_1 )
    print(paste0(cyc,"——","-------"))
    print(mod1) 
    
    }
  }
  
  
setwd("path/to/your/file")
genecopies = read.csv("gene_copies2.csv",row.names = 1)
group = read.csv("group.csv")
genecopies_sum = aggregate(genecopies$gene_copies,by = list(genecopies$sample,genecopies$KO),FUN = "sum")
colnames(genecopies_sum) = c("sample","ko","gene_copies")

KO = read.csv("新KO.csv")
vKO = KO$KEGG_KO
data = genecopies_sum[genecopies_sum$ko %in% vKO,]
data = merge(data,group,by = "sample")
data <- subset(data,gene_copies>0)
input = merge(data,carbonKO,by.x = "ko",by.y ="KEGG_KO" )
sum = aggregate(input$gene_copies,by = list(input$sample,input$cycle),FUN = "sum")
colnames(sum) = c("sample","cycle","gene_copies")

for (cyc in c("碳循环","氮循环","磷循环","硫循环")) {
  tempheat <- subset(sum,sum$cycle==cyc)
  tempheat = merge(tempheat,group,by = "sample")
  mod1 <- aov( gene_copies~satru+treat, data= tempheat)
  print(paste0(cyc,"——","-------"))
  print(summary(mod1)) 
  }


