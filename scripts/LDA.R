library(vegan)
library(reshape2)
library(dplyr)
library(MASS)

{

  wd = c(
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file"
    )
  for (taxa in c("16s","its")){
    print(paste0("现在处理",taxa))
    if (taxa=="16s"){
      setwd("path/to/your/file")
      tax_table <- read.delim('metadata.txt', check.names = FALSE,sep ="path/to/your/file")
      otu_wd = wd[1:2]
    }
    else{
      setwd("path/to/your/file")
      tax_table <- read.delim('metadata.txt', check.names = FALSE,sep ="path/to/your/file")
      otu_wd = wd[3:4]
    }
    tax_table[] <- lapply(tax_table, function(x) {
      x_char <- as.character(x)
      x_char[is.na(x_char) | x_char == ""] <- "Unclassified"
      return(x_char)
    })
    colnames(tax_table) <- c("OTU_ID", "kingdom", "phylum",  "class",   "order",   "family",  "genus",   "species")
    tax_table$kingdom[tax_table$kingdom == "Unclassified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "Unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "Unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "Unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "Unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "Unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "Unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unclassified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unidentified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unidentified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "d__unidentified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "p__unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "c__unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "o__unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "f__unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "g__unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "s__unidentified"] <- "s__Unclassified"
    tax_table$kingdom <- sub("^.", "k", tax_table$kingdom)
    
    
    for (target_wd in otu_wd ){
      setwd(target_wd)
      otulist = list.files(pattern="*.csv")
      print(paste0("当前文件夹",target_wd,"path/to/your/file"))
      print(otulist)
      for (i in 1:length(otulist)) {
        otu_table <- read.csv(otulist[i], header=TRUE, check.names=F)
        colnames(otu_table)[1] <- "OTU_ID"
        
        
        merged_data <- merge(otu_table, tax_table, by="OTU_ID")  
        merged_data$phylum =paste0(merged_data$kingdom,'.',merged_data$phylum)
        merged_data$class=paste0(merged_data$phylum,'.',merged_data$class)
        merged_data$order=paste0(merged_data$class,'.',merged_data$order)
        merged_data$family=paste0(merged_data$order,'.',merged_data$family)
        merged_data$genus=paste0(merged_data$family,'.',merged_data$genus)
        merged_data$species=paste0(merged_data$genus,'.',merged_data$species)
        filename <- paste0("qq",otulist[i])
        path <- file.path("path/to/your/file",filename)
        write.csv(merged_data,file=path,row.names=F)
        print(paste0(otulist[i],"处理完毕"))
      }
    }
  }
}


{
  library(permute)
  library(vegan)
  wd = c(
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file"
  )
  print("开始筛选相对多度在万分之一以上的otu")
  for (target_wd in wd){
    print("当前路径:")
    print(target_wd)
    setwd(target_wd)
    otu_list=otulist = list.files(pattern="*.csv")
    for(otu_file in otu_list){
      OTU_Tax16sCP <-read.csv(otu_file, header=TRUE,check.names = F)
      
      datacol = ncol(OTU_Tax16sCP)-7
      colSums(OTU_Tax16sCP[2:datacol])
      rare = sum(OTU_Tax16sCP[,2])
      
      Kingdom=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$kingdom),FUN='sum')
      
      
      Phylum=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$phylum),FUN='sum')
      
      colSums(Phylum[,2:datacol])
      
      all <-  Kingdom
      
      precount <- nrow(all)
      for (i in 1:nrow(Phylum)) {
        all[precount+i,] <- Phylum[i,]
      }
      
      
      
      Class=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$class),FUN='sum')
      
      precount <- nrow(all)
      for (i in 1:nrow(Class)) {
        all[precount+i,] <- Class[i,]
      }
      
      
      Order=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$order),FUN='sum')
      
      precount <- nrow(all)
      for (i in 1:nrow(Order)) {
        all[precount+i,] <- Order[i,]
      }
      
      
      Family=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$family),FUN='sum')
      
      precount <- nrow(all)
      for (i in 1:nrow(Family)) {
        all[precount+i,] <- Family[i,]
      }
      
      
      Genus=aggregate(x=OTU_Tax16sCP[,2:datacol]/rare,by=list(OTU_Tax16sCP$genus),FUN='sum')
      
      precount <- nrow(all)
      for (i in 1:nrow(Genus)) {
        all[precount+i,] <- Genus[i,]
      }
      sum_taxa_name <- paste0(substr(otu_file,1,nchar(otu_file)-4),"分类求和.csv")
      sum_taxa_path <- file.path(target_wd,"sum_taxa")
      dir.create(sum_taxa_path,showWarnings = F)
      write.csv(all, file =file.path(sum_taxa_path,sum_taxa_name),row.names = FALSE)
      
      colSums(Phylum[,2:datacol])
      colSums(Class[,2:datacol])
      colSums(Order[,2:datacol])
      colSums(Family[,2:datacol])
      colSums(Genus[,2:datacol])
      
      all$average = rowMeans(all[,2:datacol])
      screen0.0001 <-all
      result <- subset(screen0.0001[which(screen0.0001$average>0.0001),])
      filter_sum_name <- paste0(substr(otu_file,1,nchar(otu_file)-4),"筛选万分之一.csv")
      filter_sum_path <- file.path(target_wd,"filter_sum_taxa")
      dir.create(filter_sum_path,showWarnings = F)
      write.csv(result[,1:ncol(result)-1], file = file.path(filter_sum_path,filter_sum_name),row.names = FALSE)
      print(paste0("完成",otu_file,"的筛选"))
    }
  }
}


{
  library(reshape2)
  library(dplyr)
  library(multcomp)  
  library(car)       
  wd = c(
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file"
  )
  eco_wd <- wd[c(1,2)]
  treat_wd <- wd[c(3,4)]
  for (otu_type in c("eco","treat")){
    if  (otu_type=="eco"){
      temp_wd <- eco_wd
    }else{
      temp_wd <- treat_wd
    }
    print(paste0("处理",otu_type))
    for (target_wd in temp_wd){
      setwd(target_wd)
      print(paste0("文件路径",target_wd))
      otu_list <- list.files(pattern="*.csv")
      group <- read.csv(file = file.path(dirname(target_wd),'group','LEfSE分组.csv'),header = T)
      for (otu_file in otu_list) {
        data <-read.csv(otu_file, header=TRUE,row.names = 1)
        data_t <- as.data.frame(t(data))
        cols_count <- ncol(data_t)
        data_t$SampleID <- rownames(data_t)
        data_t_merged <- merge(data_t, group, by="SampleID")
        rownames(data_t_merged) <- data_t_merged[,1]
        data_t_merged <- data_t_merged[,-1]
        result <-  data.frame(
          taxa = character(0),
          clade_value = numeric(0),
          stringsAsFactors = FALSE
        )
        for (col_index in 1:cols_count) {
          temp_data <- data_t_merged[,c(col_index,cols_count+2)]
          result_taxa <- colnames(temp_data)[1]
          clade_value <- 6*log(sum(temp_data[,1])*10000)
          result <- rbind(result, data.frame(taxa = result_taxa, clade_value = clade_value))
        }
        clade_file_name <- paste0(substr(otu_file,1,nchar(otu_file)-4),"_clade_value.csv")
        clade_file_path <- file.path(dirname(target_wd),"clade_value")
        dir.create(clade_file_path,showWarnings = F)
        write.csv(result,file = file.path(clade_file_path,clade_file_name),row.names = F)
        
        print(paste0(otu_file,"处理完成"))
        
      }
      
    }
  }
  
}


{
  library(tidyverse)
  library(microeco)
  library(magrittr)
  
  library(snowfall) 
  library(parallel)
  
  
  wd <- c(
    "path/to/your/file",
    "path/to/your/file"
  )
  for (target_wd in wd){
    setwd(dirname(target_wd))
    tax_table <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
    tax_table[] <- lapply(tax_table, function(x) {
      x_char <- as.character(x)
      x_char[is.na(x_char) | x_char == ""] <- "Unclassified"
      return(x_char)
    })
    colnames(tax_table) <- c("kingdom", "phylum",  "class",   "order",   "family",  "genus",   "species")
    tax_table$kingdom[tax_table$kingdom == "Unclassified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "Unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "Unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "Unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "Unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "Unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "Unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unclassified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unidentified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unidentified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "d__unidentified"] <- "k__Unclassified"
    tax_table$phylum[tax_table$phylum == "p__unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "c__unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "o__unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "f__unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "g__unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "s__unidentified"] <- "s__Unclassified"
    tax_table$kingdom <- sub("^.", "k", tax_table$kingdom)
    
    setwd(target_wd)
    otulist = list.files(pattern="*.csv")
    for(otu_file in otulist){
      setwd(dirname(target_wd))
      sample_table <- read.csv('LEfSE分组.csv', row.names = 1, check.names = FALSE)
      setwd(target_wd)
      feature_table<- read.csv(otu_file, row.names = 1)
      name <- substr(otu_file,1,nchar(otu_file)-4)
      otu_colname <- colnames(feature_table)
      sample_table$rowselect <- rownames(sample_table)
      sample_table <- sample_table[sample_table$rowselect %in% otu_colname, ]
      
      head(feature_table)[1:6,1:6]; head(sample_table)[1:6, ]; head(tax_table)[,1:6]
      
      dataset <- microtable$new(sample_table = sample_table,
                                otu_table = feature_table, 
                                tax_table = tax_table)
      lefse <- trans_diff$new(dataset = dataset, 
                              method = "lefse", 
                              group = "group1", 
                              alpha =  0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")
      file_name <- paste0(name,'_lefes.csv')
      write.csv(as.data.frame(lefse$res_diff), file=file.path("path/to/your/file",file_name), row.names = F)
    }
  }
  
  wd <- c(
    "path/to/your/file",
    "path/to/your/file"
  )
  for (target_wd in wd){
    setwd(dirname(target_wd))
    tax_table <- read.delim('metadata.txt', row.name = 1, check.names = FALSE)
    tax_table[] <- lapply(tax_table, function(x) {
      x_char <- as.character(x)
      x_char[is.na(x_char) | x_char == ""] <- "Unclassified"
      return(x_char)
    })
    colnames(tax_table) <- c("kingdom", "phylum",  "class",   "order",   "family",  "genus",   "species")
    tax_table$kingdom[tax_table$kingdom == "Unclassified"] <- "d__Unclassified"
    tax_table$phylum[tax_table$phylum == "Unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "Unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "Unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "Unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "Unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "Unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unclassified"] <- "d__Unclassified"
    tax_table$phylum[tax_table$phylum == "unclassified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unclassified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unclassified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unclassified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unclassified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unclassified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "unidentified"] <- "d__Unclassified"
    tax_table$phylum[tax_table$phylum == "unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "unidentified"] <- "s__Unclassified"
    
    tax_table$kingdom[tax_table$kingdom == "d__unidentified"] <- "d__Unclassified"
    tax_table$phylum[tax_table$phylum == "p__unidentified"] <- "p__Unclassified"
    tax_table$class[tax_table$class == "c__unidentified"] <- "c__Unclassified"
    tax_table$order[tax_table$order == "o__unidentified"] <- "o__Unclassified"
    tax_table$family[tax_table$family == "f__unidentified"] <- "f__Unclassified"
    tax_table$genus[tax_table$genus == "g__unidentified"] <- "g__Unclassified"
    tax_table$species[tax_table$species == "s__unidentified"] <- "s__Unclassified"
    tax_table$kingdom <- sub("^.", "k", tax_table$kingdom)
    
    setwd(target_wd)
    otulist = list.files(pattern="*.csv")
    for(otu_file in otulist){
      setwd(dirname(target_wd))
      sample_table <- read.csv('LEfSE分组.csv', row.names = 1, check.names = FALSE)
      setwd(target_wd)
      feature_table<- read.csv(otu_file, row.names = 1)
      name <- substr(otu_file,1,nchar(otu_file)-4)
      otu_colname <- colnames(feature_table)
      sample_table$rowselect <- rownames(sample_table)
      sample_table <- sample_table[sample_table$rowselect %in% otu_colname, ]
      
      head(feature_table)[1:6,1:6]; head(sample_table)[1:6, ]; head(tax_table)[,1:6]
      
      dataset <- microtable$new(sample_table = sample_table,
                                otu_table = feature_table, 
                                tax_table = tax_table)
      lefse <- trans_diff$new(dataset = dataset, 
                              method = "lefse", 
                              group = "Group", 
                              alpha =  0.05, 
                              lefse_subgroup = NULL,
                              p_adjust_method = "none")
      file_name <- paste0(name,'_lefes.csv')
      write.csv(as.data.frame(lefse$res_diff), file=file.path("path/to/your/file",file_name), row.names = F)
    }
  }
}


{
  library(dplyr)
  library(stringr)
  wd = c(
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file"
  )
  for (target_wd in wd){
    setwd(target_wd)
    otu_list <- list.files(pattern="*.csv")
    for (otu_file in otu_list){
      tempdata <- read.csv(file=otu_file,header = T)
      tempdata$Taxa <- str_replace_all(tempdata$Taxa, fixed("|"), ".")
      tempdata <- tempdata %>%
        mutate(color_value = case_when(
          Group == "CK" ~ '
          Group == "C" ~ '
          Group == "N" ~ '
          Group == "pH" ~ '
          Group == "T" ~ '
          Group == "W" ~ '
          Group == "10%re" ~ '
          Group == "40%re" ~ '
          Group == "70%re" ~ '
          Group == "100%re" ~ '
          TRUE ~ '
        ))
      write.csv(tempdata,file =otu_file,row.names = F)
    }
  }
}


{
  library(stringr)
  wd = c(
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file",
    "path/to/your/file"
  )
  
  for (target_wd in wd){
    setwd(file.path(target_wd,'lda'))
    lefes_list <- list.files(pattern="*.csv")
    setwd(file.path(target_wd,'clade_value'))
    clade_list <- list.files(pattern="*.csv")
    for (lefes_file in lefes_list){
      lefes_data <- read.csv(file=file.path(target_wd,'lda',lefes_file),header = T)
      otu_name <- substr(lefes_file,1,nchar(lefes_file)-10)
      match_clade <- grep(otu_name, clade_list,value = TRUE)
      clade_data <- read.csv(file=file.path(target_wd,'clade_value',match_clade),header = T)
      colnames(clade_data) <- c('Taxa','clade_value')
      merge_data <- merge(lefes_data, clade_data, by="Taxa",all=T)
      file_name <- paste0(otu_name,"_lda_clade_merged.csv")
      write.csv(merge_data,file=file.path("path/to/your/file",file_name),row.names = F)
    }
  }
}


{
  wd <- "path/to/your/file"
  setwd(wd)
  otu_list <- list.files(pattern="*.csv")
  for (otu_file in otu_list){
    data <- read.csv(otu_file,header = T)
    data$Group[data$Group == "W"] <- "D"
    data$Group[data$Group == "T"] <- "W"
    write.csv(data,file =otu_file,row.names = F)
  }
  
  wd <- "path/to/your/file"
  setwd(wd)
  otu_list <- list.files(pattern="*.csv")
  for (otu_file in otu_list){
    data <- read.csv(otu_file,header = T)
    data <- data[!(is.na(data$clade_value) | data$clade_value == ""), ]
    rows_to_fix <- is.na(data$LDA) | data$LDA == ""
    data$LDA[rows_to_fix] <- 4
    data$color_value[rows_to_fix] <- "
    print(paste0(otu_file,"done"))
    write.csv(data,file =otu_file,row.names = F)
  }
  
}


{
  wd <- "path/to/your/file"
  setwd(wd)
  otu_list <- list.files(pattern="*.csv")
  process_taxa_string <- function(s) {
    if (is.na(s) || s == "") return(s)
    
    s_no_paren <- gsub("path/to/your/file", "", s)
    if (s!=s_no_paren){
      print(paste0("去除",s,"中的括号，去除后",s_no_paren))
    }
    
    parts <- strsplit(s_no_paren, "path/to/your/file")[[1]]
    if (length(parts) == 1) return(parts)  
    
    valid_parts <- character(0)
    for (i in seq_along(parts)) {
      if (i == 1) {
        valid_parts <- parts[i]
      } else {
        prev_part <- parts[i-1]
        curr_part <- parts[i]
        
        if (grepl("^[a-zA-Z]__[a-zA-Z].*", curr_part)) {
          valid_parts <- c(valid_parts, curr_part)
        } else {
          print(paste0(s,"检测到非法小数点"))
          valid_parts[length(valid_parts)] <- paste0(valid_parts[length(valid_parts)], curr_part)
        }
      }
    }
    
    paste(valid_parts, collapse = ".")
  }
  
  for (otu_file in otu_list){
    data <- read.csv(otu_file,header = T)
    data$Taxa <- sapply(data$Taxa, process_taxa_string)
    
    write.csv(data,file =otu_file,row.names = F)
  }
}

{
  wd <- "path/to/your/file"
  setwd(wd)
  otu_list <- list.files(pattern = "*.csv")
  
  if (!require("dplyr")) install.packages("dplyr")
  library(dplyr)
  
  for (otu_file in otu_list) {
    data <- read.csv(otu_file, header = TRUE, stringsAsFactors = FALSE)
    
    data$depth <- sapply(data$Taxa, function(x) lengths(regmatches(x, gregexpr("path/to/your/file", x))))
    
    data$is_placeholder <- ifelse(data$color_value == "
    
    father_taxa <- function(x) {
      parts <- unlist(strsplit(x, "path/to/your/file"))
      if (length(parts) > 1) {
        return(paste(parts[1:(length(parts)-1)], collapse = "."))
      } else {
        return(NA)
      }
    }
    data$parent_taxa <- sapply(data$Taxa,father_taxa )
        
    valid_nodes <- data$Taxa[!data$is_placeholder]
    
    find_ancestors <- function(node) {
      ancestors <- c()
      current <- node
      while (current %in% data$Taxa) {
        parent <- data$parent_taxa[data$Taxa == current][1]
        if (!is.na(parent)) {  
          ancestors <- c(ancestors, parent)
          current <- parent
        } else {
          break
        }
      }
      return(ancestors)
    }
    
    all_ancestors <- unique(unlist(sapply(valid_nodes, find_ancestors)))
    
    nodes_to_keep <- unique(c(valid_nodes, all_ancestors))
    
    filtered_data <- data[data$Taxa %in% nodes_to_keep, ]
    
    cols_to_remove <- c("depth", "is_placeholder", "parent_taxa")
    
    for (col in cols_to_remove) {
      if (col %in% names(filtered_data)) {
        filtered_data[[col]] <- NULL
      }
    }
    
    print(paste0(otu_file, " processed. Original rows: ", nrow(data), 
                 ", Filtered rows: ", nrow(filtered_data)))
    
    write.csv(filtered_data, file = otu_file, row.names = FALSE)
  }
}
