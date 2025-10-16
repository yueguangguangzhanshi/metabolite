#!/usr/bin/env Rscript

Sys.setenv(RSTUDIO_PANDOC="/usr/bin/pandoc")
Sys.setenv(OPENSSL_CONF="/dev/null")

args <- commandArgs(T)

project_id<-args[1]

metr_pkgs<-c('visNetwork', 'MicrobiomeProfiler', 'pathview', 'MetaboAnalystR', 'ggraph', 'igraph', 'grid', 'UpSetR', 'VennDiagram', 'corrplot', 'glue' ,'matrixStats','missForest','VIM','ggord','survminer','survival','highcharter','ropls','dplyr','reshape2','statnet','circlize','tibble','corrplot','psych','ComplexHeatmap','org.Hs.eg.db','limma','WebGestaltR','ggraph','igraph','plyr',"qvalue","digest", "BiocManager", "Cairo","ggrepel","rgl","knitr","mixOmics","pheatmap","RColorBrewer","gplots","ggplot2","openxlsx","stringr","scales","data.table","factoextra","FactoMineR","XML","downloader","ggsci","rJava","rChoiceDialogs","foreach","grid","seqinr","reshape","ontologyIndex","qvalue","VennDiagram","RSQLite","readr","tidyverse","readxl","gdata","svDialogs","ggseqlogo","rmotifx","rvest","rlist", "knitr", "rmarkdown", "bsselectR",'echarts4r', 'shiny', 'DT', 'plotly', 'ggbiplot','purrr')

for(i in 1:length(metr_pkgs)){
  
  library(metr_pkgs[i], character.only = TRUE)
  
}

project_db_path <- "/home/biognosis/shiny-server/BioAnalysis/db/projects_list.db"

# path_prefix <- '/home/biognosis/mnt/os2/BioAnalysis/'
path_prefix <- '/home/biognosis/shiny-server/BioAnalysis/'

source(paste0(path_prefix, '/etc/dep/bg_lib.r'), local = TRUE, echo = FALSE)

UpdateStatus <- function(status,id) {
  
  system("cp -rf /home/biognosis/mnt/os2/BioAnalysis/db/projects_list.db /home/biognosis/shiny-server/BioAnalysis/db/projects_list.db")
  
  projects_list <- dbConnect(RSQLite::SQLite(), project_db_path)
  
  dbSendStatement(projects_list, paste0("UPDATE Project SET status = '",status,"' WHERE id='",id,"'"))
  
  dbDisconnect(projects_list)
  
  system("cp -rf /home/biognosis/shiny-server/BioAnalysis/db/projects_list.db /home/biognosis/mnt/os2/BioAnalysis/db/projects_list.db")
  
}

#0:排队 1:完成 2:等待 3:错误
# project_id <- "174435965127oGpytUzjM7"

data_path <- paste0(path_prefix, '/data/usrdata/', project_id,'/data')

report_path <- paste0(path_prefix, '/data/usrdata/', project_id,'/report')

para_path <- paste0(data_path, '/para.Rdata')

if(file.exists(para_path)){
  
  load(para_path)
  
  input <- readRDS(paste0(data_path,"/input.rds"))
  
  rv <- readRDS(paste0(data_path,"/rv.rds"))
  
  tryCatch({
    
    if(para[["info"]][["type"]] == 'pre_process'){
      
      dir.create(report_path)
      
      setwd(report_path)
      
      protein.file <- list.files(path = data_path, pattern = ".*tsv$",recursive = T,full.names = T)
      
      protein.data <- read_tsv(protein.file)
      
      group.file <- list.files(path = data_path, pattern = "IdentificationsOverview.*xls$",recursive = T,full.names = T)
      
      sample.data <- read_xls(group.file)
      
      quant.data.raw <- protein.data[,c("PG.ProteinAccessions",grep("PG.Quantity$",colnames(protein.data),value = T))]%>%column_to_rownames("PG.ProteinAccessions")
      
      match.index <- match(str_replace_all(colnames(quant.data.raw),"\\[.*\\] |(.PG.Quantity)",""),sample.data$FileName)
      
      analysis.name <- sample.data$`Sample name(分析用的名称）`[match.index]
      
      colnames(quant.data.raw) <- analysis.name
      
      group <- sample.data$`Condition（分组）`[match.index]
      
      ## 1. 原始数据 pca ----
      dir.create("1.原始数据PCA")
      
      wb <- createWorkbook()
      addWorksheet(wb, "原始数据")
      writeData(wb, "原始数据", quant.data.raw, rowNames = TRUE)
      saveWorkbook(wb, "1.原始数据PCA/原始数据.xlsx", overwrite = TRUE)
      
      if (para[['condition']][['miss.value.handle.method']]=="none") {
        
        ### 不含缺失值
        quant.data <- na.omit(quant.data.raw)%>%as.data.frame()
        
        tryCatch({
          
          quant.data[,1:ncol(quant.data)]<-lapply(quant.data[,1:ncol(quant.data)],as.numeric)
          
          if (para[['condition']][['miss.value.handle.show']]=="normal.qc") {
            
            get.pca.plot(quant.data,group,name="1.原始数据PCA/不含缺失值的原始数据PCA.pdf")
            
          } else if (para[['condition']][['miss.value.handle.show']]=="normal") {
            
            get.pca.plot(quant.data[,which(group!="QC")],group[which(group!="QC")],name="1.原始数据PCA/不含缺失值的原始数据PCA.pdf")
            
          }
          
        }, error=function(e){
          
          ""
          
        })
        
      }else if (para[['condition']][['miss.value.handle.method']]=="zero") {
        
        ### 缺失值全部设为 0
        quant.data <- quant.data.raw
        
        quant.data[is.na(quant.data)]<-0
        quant.data[quant.data=="NaN"]<-0
        
        quant.data[,1:ncol(quant.data)]<-lapply(quant.data[,1:ncol(quant.data)],as.numeric)
        
        if (para[['condition']][['miss.value.handle.show']]=="normal.qc") {
          
          get.pca.plot(quant.data,group,name="1.原始数据PCA/缺失值全部设为0的原始数据PCA.pdf")
          
        } else if (para[['condition']][['miss.value.handle.show']]=="normal") {
          
          get.pca.plot(quant.data[,which(group!="QC")],group[which(group!="QC")],name="1.原始数据PCA/缺失值全部设为0的原始数据PCA.pdf")
          
        }
        
      }
      
      ## 2. 缺失值过滤 ----
      dir.create("2.缺失值过滤")
      
      if (para[['condition']][['miss.value.handle.batch']]=="none.batch") {
        
        if (para[['condition']][['miss.value.handle.group']]=="inter.group"){
          
          ### 组内缺失值过滤
          missing.index <- list()
          
          for (i in unique(group)) {
            
            missing.index[[i]] <- which(apply(quant.data.raw[,which(group==i)], 1, function(x){sum(is.na(x))>(length(x)*para[['condition']][['miss.value.handle.cutoff']]/100)}))
            
          }
          
          if (para[['condition']][["miss.value.handle.condition"]]=="none.inter.condition.all") {
            
            missing.index <- Reduce(intersect,missing.index)
            
          }else if(para[['condition']][["miss.value.handle.condition"]]=="none.inter.condition.any"){
            
            missing.index <- Reduce(union,missing.index)
            
          }
          
          quant.data2 <- quant.data.raw[-missing.index,]
          
          wb <- createWorkbook()
          addWorksheet(wb, "缺失值过滤数据")
          writeData(wb, "缺失值过滤数据", quant.data2, rowNames = TRUE)
          saveWorkbook(wb, "2.缺失值过滤/组内缺失值过滤数据.xlsx", overwrite = TRUE)
          
        }else if(para[['condition']][['miss.value.handle.group']]=="golbal.group"){
          
          ### 全局缺失值过滤
          quant.data2 <- quant.data.raw[-which(apply(quant.data.raw, 1, function(x){sum(is.na(x))>(length(x)*para[['condition']][['miss.value.handle.cutoff']]/100)})),]
          
          wb <- createWorkbook()
          addWorksheet(wb, "缺失值过滤数据")
          writeData(wb, "缺失值过滤数据", quant.data2, rowNames = TRUE)
          saveWorkbook(wb, "2.缺失值过滤/全局缺失值过滤数据.xlsx", overwrite = TRUE)
          
        }
        
      }else if(para[['condition']][['miss.value.handle.batch']]=="across.batch"){
        
        batch <- sample.data$Batch[match.index]
        
        missing.index <- list()
        
        for (i in unique(batch)) {
          
          missing.index[[i]] <- which(apply(quant.data.raw[,which(batch==i)], 1, function(x){sum(is.na(x))>(length(x)*para[['condition']][['miss.value.handle.cutoff']]/100)}))
          
        }
        
        if (para[['condition']][["miss.value.handle.condition"]]=="across.condition.all") {
          
          missing.index <- Reduce(intersect,missing.index)
          
        }else if(para[['condition']][["miss.value.handle.condition"]]=="across.condition.any"){
          
          missing.index <- Reduce(union,missing.index)
          
        }
        
        quant.data2 <- quant.data.raw[-missing.index,]
        
        wb <- createWorkbook()
        addWorksheet(wb, "缺失值过滤数据")
        writeData(wb, "缺失值过滤数据", quant.data2, rowNames = TRUE)
        saveWorkbook(wb, "2.缺失值过滤/缺失值过滤数据.xlsx", overwrite = TRUE)
        
      }
      
      ## 3. 剩余缺失值处理方法 ----
      ### 缺失值做填充 ----
      if (para[['condition']][['miss.value.handle.remain']]=="remain.fill") {
        
        #### 4. 填充方式 ----
        dir.create("4.缺失值填充")
        
        ##### 不填充 ----
        if (para[['condition']][['miss.value.fill']]=="none"){
          
          quant.data.fill <- quant.data2
          
        }
        
        ##### 根据组内有效数据占比 ----
        if (para[['condition']][['miss.value.fill']]=="inter.valid") {
          
          quant.data.fill <- list()
          
          for (i in unique(group)) {
            
            quant.data3 <- quant.data2[,which(group==i)]
            
            if (sum(!is.na(quant.data3),na.rm = T)<nrow(quant.data3)*ncol(quant.data3)*0.5) {
              
              quant.data.fill[[i]] <- quant.data3
              
              # 如果组内有效数据占比< 50%
              if (para[['condition']][['miss.value.valid.below50']]=="min.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- min(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.below50']]=="mean.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- mean(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.below50']]=="median.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- median(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.below50']]=="knn.inter") {
                
                knn.result <- kNN(quant.data.fill[[i]])
                rownames(knn.result) <- rownames(quant.data.fill[[i]])
                quant.data.fill[[i]] <- knn.result[,1:ncol(quant.data.fill[[i]])]
                
              }
              
              if (para[['condition']][['miss.value.valid.below50']]=="user.defined") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- para[['condition']][[para[['condition']][['miss.value.valid.below50']]]]
                
              }
              
            }else{
              
              quant.data.fill[[i]] <- quant.data3
              
              # 如果组内有效数据占比 ≥ 50%
              if (para[['condition']][['miss.value.valid.over50']]=="min.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- min(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.over50']]=="mean.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- mean(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.over50']]=="median.inter") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- median(as.matrix(quant.data.fill[[i]]),na.rm = T)
                
              }
              
              if (para[['condition']][['miss.value.valid.over50']]=="knn.inter") {
                
                knn.result <- kNN(quant.data.fill[[i]])
                rownames(knn.result) <- rownames(quant.data.fill[[i]])
                quant.data.fill[[i]] <- knn.result[,1:ncol(quant.data.fill[[i]])]
                
              }
              
              if (para[['condition']][['miss.value.valid.over50']]=="user.defined") {
                
                quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- para[['condition']][[para[['condition']][['miss.value.valid.over50']]]]
                
              }
              
            }
            
          }
          
          quant.data.fill <- do.call(bind_cols,quant.data.fill)
          
          quant.data.fill <- quant.data.fill[,match(colnames(quant.data2),colnames(quant.data.fill))]
          
        }
        
        ##### 根据全局有效数据占比 ----
        if (para[['condition']][['miss.value.fill']]=="global.valid") {
          
          if (sum(!is.na(quant.data2),na.rm = T)<nrow(quant.data2)*ncol(quant.data2)*0.5) {
            
            quant.data.fill <- quant.data2
            
            # 如果组内有效数据占比< 50%
            if (para[['condition']][['miss.value.valid.below50']]=="min.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- min(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.below50']]=="mean.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- mean(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.below50']]=="median.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- median(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.below50']]=="knn.inter") {
              
              knn.result <- kNN(quant.data.fill)
              rownames(knn.result) <- rownames(quant.data.fill)
              quant.data.fill <- knn.result[,1:ncol(quant.data.fill)]
              
            }
            
            if (para[['condition']][['miss.value.valid.below50']]=="user.defined") {
              
              quant.data.fill[is.na(quant.data.fill)] <- para[['condition']][[para[['condition']][['miss.value.valid.below50']]]]
              
            }
            
          }else{
            
            quant.data.fill <- quant.data2
            
            # 如果组内有效数据占比 ≥ 50%
            if (para[['condition']][['miss.value.valid.over50']]=="min.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- min(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.over50']]=="mean.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- mean(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.over50']]=="median.inter") {
              
              quant.data.fill[is.na(quant.data.fill)] <- median(as.matrix(quant.data.fill),na.rm = T)
              
            }
            
            if (para[['condition']][['miss.value.valid.over50']]=="knn.inter") {
              
              knn.result <- kNN(quant.data.fill)
              rownames(knn.result) <- rownames(quant.data.fill)
              quant.data.fill <- knn.result[,1:ncol(quant.data.fill)]
              
            }
            
            if (para[['condition']][['miss.value.valid.over50']]=="user.defined") {
              
              quant.data.fill[is.na(quant.data.fill)] <- para[['condition']][[para[['condition']][['miss.value.valid.over50']]]]
              
            }
            
          }
          
        }
        
        ##### 简洁填充 ----
        if (para[['condition']][['miss.value.fill']]=="min.inter") {
          
          quant.data.fill <- list()
          
          for (i in unique(group)) {
            
            quant.data.fill[[i]] <- quant.data2[,which(group==i)]
            
            quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- min(as.matrix(quant.data.fill[[i]]),na.rm = T)
            
          }
          
          quant.data.fill <- do.call(bind_cols,quant.data.fill)
          
          quant.data.fill <- quant.data.fill[,match(colnames(quant.data2),colnames(quant.data.fill))]
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="mean.inter") {
          
          quant.data.fill <- list()
          
          for (i in unique(group)) {
            
            quant.data.fill[[i]] <- quant.data2[,which(group==i)]
            
            quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- mean(as.matrix(quant.data.fill[[i]]),na.rm = T)
            
          }
          
          quant.data.fill <- do.call(bind_cols,quant.data.fill)
          
          quant.data.fill <- quant.data.fill[,match(colnames(quant.data2),colnames(quant.data.fill))]
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="median.inter") {
          
          quant.data.fill <- list()
          
          for (i in unique(group)) {
            
            quant.data.fill[[i]] <- quant.data2[,which(group==i)]
            
            quant.data.fill[[i]][is.na(quant.data.fill[[i]])] <- median(as.matrix(quant.data.fill[[i]]),na.rm = T)
            
          }
          
          quant.data.fill <- do.call(bind_cols,quant.data.fill)
          
          quant.data.fill <- quant.data.fill[,match(colnames(quant.data2),colnames(quant.data.fill))]
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="knn.inter") {
          
          quant.data.fill <- list()
          
          for (i in unique(group)) {
            
            quant.data.fill[[i]] <- quant.data2[,which(group==i)]
            
            knn.result <- kNN(quant.data.fill[[i]])
            rownames(knn.result) <- rownames(quant.data.fill[[i]])
            quant.data.fill[[i]] <- knn.result[,1:ncol(quant.data.fill[[i]])]
            
          }
          
          quant.data.fill <- do.call(bind_cols,quant.data.fill)
          
          quant.data.fill <- quant.data.fill[,match(colnames(quant.data2),colnames(quant.data.fill))]
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="min.global") {
          
          quant.data.fill <- quant.data2
          
          quant.data.fill[is.na(quant.data.fill)] <- min(as.matrix(quant.data.fill),na.rm = T)
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="mean.global") {
          
          quant.data.fill <- quant.data2
          
          quant.data.fill[is.na(quant.data.fill)] <- mean(as.matrix(quant.data.fill),na.rm = T)
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="median.global") {
          
          quant.data.fill <- quant.data2
          
          quant.data.fill[is.na(quant.data.fill)] <- median(as.matrix(quant.data.fill),na.rm = T)
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="knn.global") {
          
          knn.data.temp <- kNN(quant.data2)
          quant.data.fill <- knn.data.temp[,1:ncol(quant.data2)]
          rownames(quant.data.fill) <- rownames(quant.data2)
          
        }
        
        if (para[['condition']][['miss.value.fill']]=="user.defined") {
          
          quant.data.fill <- quant.data2
          
          quant.data.fill[is.na(quant.data.fill)] <- para[['condition']][["miss.value.fill.defined"]]
          
        }
        
        wb <- createWorkbook()
        addWorksheet(wb, "缺失值填充数据")
        writeData(wb, "缺失值填充数据", quant.data.fill, rowNames = TRUE)
        saveWorkbook(wb, "4.缺失值填充/缺失值填充数据.xlsx", overwrite = TRUE)
        
      }
      
      ### 保留包含缺失值的行，但在数据集中标记这些行 ----
      if(para[['condition']][['miss.value.handle.remain']]=="remain.reserve"){
        
        quant.data.fill <- quant.data2
        
      }
      
      ## 5. 中位数处理 ----
      if (para[['condition']][['normalized.handle.method']]=="none"){""}else{
        
        dir.create("5.中位数处理")
        
        ### 中位数平滑归一化
        if (para[['condition']][['normalized.handle.method']]=="median.polish") {
          
          polish.median<-apply(quant.data.fill,2,median,na.rm=T)
          
          for (i in 1:ncol(quant.data.fill)) {
            
            quant.data.fill[,i]<-(quant.data.fill[,i]-polish.median[i])/(median(abs(quant.data.fill[,i]-polish.median[i]),na.rm = T)+1)
            
          }
          
        }
        
        ### 中位数归一化
        if (para[['condition']][['normalized.handle.method']]=="median") {
          
          median.data<-apply(quant.data.fill,2,median,na.rm=T)
          
          for (i in 1:ncol(quant.data.fill)) {
            
            quant.data.fill[,i]<-quant.data.fill[,i]/median.data[i]
            
          }
          
        }
        
        wb <- createWorkbook()
        addWorksheet(wb, "中位数处理")
        writeData(wb, "中位数处理", quant.data.fill, rowNames = TRUE)
        saveWorkbook(wb, "5.中位数处理/中位数处理数据.xlsx", overwrite = TRUE)
        
      }
      
      ## 6. 批次矫正 ----
      if (para[['condition']][['miss.value.handle.batch']]=="none.batch"){""}else{
        
        if (para[['condition']][['batch.correct']]=="none"){""}else{
          
          dir.create("6.批次校正")
          
          ### combat校正
          if (para[['condition']][['batch.correct']]=="combat") {
            
            quant.data.fill<-sva::ComBat(quant.data.fill%>%as.matrix(),batch = batch)
            
          }
          
          ### 参考样本（指定一个已经上传的样本，做除法）
          if (para[['condition']][['batch.correct']]=="reference.sample") {
            
            reference.data <- quant.data.fill[,para[['condition']][['batch.correct.reference.sample']]]
            
            for (i in 1:ncol(quant.data.fill)) {
              
              quant.data.fill[,i]<-quant.data.fill[,i]/reference.data[i]
              
            }
            
          }
          
          wb <- createWorkbook()
          addWorksheet(wb, "批次矫正")
          writeData(wb, "批次矫正", quant.data.fill, rowNames = TRUE)
          saveWorkbook(wb, "6.批次矫正/批次矫正数据.xlsx", overwrite = TRUE)
          
        }
        
      }
      
      ## 7. 处理后的数据 PCA ----
      dir.create("7.处理后的数据PCA")
      
      if (para[['condition']][['pca.handle.method']]=="none") {
        
        ### 不含缺失值
        quant.data.fill <- na.omit(quant.data.fill)%>%as.data.frame()
        
        tryCatch({
          
          quant.data.fill[,1:ncol(quant.data.fill)]<-lapply(quant.data.fill[,1:ncol(quant.data.fill)],as.numeric)
          
          if (para[['condition']][['pca.handle.show']]=="normal.qc") {
            
            get.pca.plot(quant.data.fill,group,name="7.处理后的数据PCA/不含缺失值的处理后的数据PCA.pdf")
            
          } else if (para[['condition']][['pca.handle.show']]=="normal") {
            
            get.pca.plot(quant.data.fill[,which(group!="QC")],group[which(group!="QC")],name="7.处理后的数据PCA/不含缺失值的处理后的数据PCA.pdf")
            
          }
          
        }, error=function(e){
          
          ""
          
        })
        
      }else if (para[['condition']][['pca.handle.method']]=="zero") {
        
        ### 缺失值全部设为 0
        quant.data.fill[is.na(quant.data.fill)]<-0
        quant.data.fill[quant.data.fill=="NaN"]<-0
        
        quant.data.fill[,1:ncol(quant.data.fill)]<-lapply(quant.data.fill[,1:ncol(quant.data.fill)],as.numeric)
        
        if (para[['condition']][['pca.handle.show']]=="normal.qc") {
          
          get.pca.plot(quant.data.fill,group,name="7.处理后的数据PCA/缺失值全部设为0的处理后的数据PCA.pdf")
          
        } else if (para[['condition']][['pca.handle.show']]=="normal") {
          
          get.pca.plot(quant.data.fill[,which(group!="QC")],group[which(group!="QC")],name="7.处理后的数据PCA/缺失值全部设为0的处理后的数据PCA.pdf")
          
        }
        
      }
      
    }
    
    if(para[["info"]][["type"]] == 'metabonomics'){
      
      dir.create(report_path, recursive = T)
      
      # # hmdb 数据库
      # load(file=paste0(path_prefix, '/db/parse_hmdb.rda'))
      # hmdb_db2 <- tidyr::separate_rows(hmdb_db, synonyms, sep = "; ")
      # hmdb_db2 <- dplyr::distinct(hmdb_db2, synonyms, .keep_all = T)
      
      setwd(report_path)
      
      # if (grepl("zip", input[['metaboQuant']][['type']])) {
      #   
      #   system(paste0('unzip ', data_path, '/metabonomics.zip -d', data_path))
      #   
      # }else{
      #   
      #   system(paste0('7z -y x ', data_path, '/qc.7z -o', data_path))
      #   
      # }
      
      metabo.file <- c(paste0(data_path, "/metabonomics_NEG.csv"), paste0(data_path, "/metabonomics_POS.csv"))
      
      names(metabo.file) <- c("NEG","POS")
      
      sample.data <- list(POS = rv$condition.data.pos, NEG = rv$condition.data.neg)
      
      # group.data <- list(POS = input$metaboQuant_compare_select_pos, NEG = input$metaboQuant_compare_select_neg)
      group.data <- input$metaboQuant_compare_select
      
      # 存储当前比较中POS和NEG的结果
      summary_stats_list <- results_mode_diff_name <- results_mode_diff <- results_mode <- metabo.data.raw <- list()
      
      # 数据预处理与各组的差异筛选
      for (ion_type in c("NEG","POS")) {
        
        metabo.data.raw[[ion_type]] <- read_csv(metabo.file[[ion_type]])
        
        metabo.data <- metabo.data.raw[[ion_type]][,c("peak_name", sample.data[[ion_type]]$`File name`)] %>% remove_rownames() %>% column_to_rownames("peak_name")
        
        metabo.data[metabo.data==input$miss.value.handle.type] <- NA
        
        group <- sample.data[[ion_type]]$Condition
        
        # 数据预处理 ----
        dir.create("数据预处理")
        
        ## 缺失值过滤 ----
        dir.create("数据预处理/缺失值过滤")
        
        if (input[['miss.value.handle.group']]=="inter.group"){
          
          ### 组内缺失值过滤 ----
          missing.index <- list()
          
          for (i in unique(group)) {
            
            missing.index[[i]] <- which(apply(metabo.data[ ,which(group==i)], 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)}))
            
          }
          
          # if (para[['condition']][["miss.value.handle.condition"]]=="none.inter.condition.all") {
          #   
          #   missing.index <- Reduce(intersect,missing.index)
          #   
          # }else if(para[['condition']][["miss.value.handle.condition"]]=="none.inter.condition.any"){
          #   
          #   missing.index <- Reduce(union,missing.index)
          #   
          # }
          
          missing.index <- Reduce(union,missing.index)
          
          if (length(missing.index)>0) {
            
            metabo.data2 <- metabo.data[-missing.index,]
            
          }else{
            
            metabo.data2 <- metabo.data
            
          }
          
          wb <- createWorkbook()
          addWorksheet(wb, "缺失值过滤数据")
          writeData(wb, "缺失值过滤数据", metabo.data2, rowNames = TRUE)
          saveWorkbook(wb, paste0("数据预处理/缺失值过滤/组内缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
        }else if(input[['miss.value.handle.group']]=="golbal.group"){
          
          ### 全局缺失值过滤 ----
          metabo.data2 <- metabo.data[-which(apply(metabo.data, 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)})),]
          
          wb <- createWorkbook()
          addWorksheet(wb, "缺失值过滤数据")
          writeData(wb, "缺失值过滤数据", metabo.data2, rowNames = TRUE)
          saveWorkbook(wb, paste0("数据预处理/缺失值过滤/全局缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
        }
        
        ## 缺失值做填充 ----
        dir.create("数据预处理/缺失值填充")
        
        ### 不填充 ----
        if (input[['miss.value.fill']]=="none"){
          
          metabo.data.fill <- metabo.data2
          
        }

        ### 简洁填充 ----
        if (input[['miss.value.fill']]=="min.global") {
          
          metabo.data.fill <- metabo.data2
          
          metabo.data.fill[is.na(metabo.data.fill)] <- min(as.matrix(metabo.data.fill),na.rm = T)
          
        }
        
        if (input[['miss.value.fill']]=="mean.global") {
          
          metabo.data.fill <- metabo.data2
          
          metabo.data.fill[is.na(metabo.data.fill)] <- mean(as.matrix(metabo.data.fill),na.rm = T)
          
        }
        
        if (input[['miss.value.fill']]=="min2") {
          
          metabo.data.fill <- metabo.data2
          
          metabo.data.fill[is.na(metabo.data.fill)] <- mean(as.matrix(metabo.data.fill),na.rm = T)/2
          
        }
        
        if (input[['miss.value.fill']]=="median.global") {
          
          metabo.data.fill <- metabo.data2
          
          metabo.data.fill[is.na(metabo.data.fill)] <- median(as.matrix(metabo.data.fill),na.rm = T)
          
        }
        
        if (input[['miss.value.fill']]=="knn.global") {
          
          knn.data.temp <- kNN(metabo.data2)
          metabo.data.fill <- knn.data.temp[,1:ncol(metabo.data2)]
          rownames(metabo.data.fill) <- rownames(metabo.data2)
          
        }
        
        if (input[['miss.value.fill']]=="rf") {
          
          row_names <- rownames(metabo.data2)
          
          # 执行随机森林填充
          set.seed(123)
          rf_imp <- missForest(
            as.matrix(metabo.data2),
            maxiter = 5,
            ntree = 100,
            verbose = FALSE  # 关闭进度显示
          )
          
          # 重构数据框
          metabo.data.fill <- as.data.frame(rf_imp$ximp)
          rownames(metabo.data.fill) <- row_names
          
        }
        
        wb <- createWorkbook()
        addWorksheet(wb, "缺失值填充数据")
        writeData(wb, "缺失值填充数据", metabo.data.fill, rowNames = TRUE)
        saveWorkbook(wb, paste0("数据预处理/缺失值填充/缺失值填充数据_", ion_type, ".xlsx"), overwrite = TRUE)
        
        ## 归一化 ----
        if (input[['normalized.handle.method']]=="none"){""}else{
          
          dir.create("数据预处理/归一化")
          
          ### sum归一化
          if (input[['normalized.handle.method']]=="sum") {
            
            sum.median<-apply(metabo.data.fill,2,sum,na.rm=T)
            
            for (i in 1:ncol(metabo.data.fill)) {
              
              metabo.data.fill[,i]<-metabo.data.fill[,i]/sum.median[i]
              
            }
            
          }
          
          ### QC样本归一化
          if (input[['normalized.handle.method']]=="qc") {
            
            
            
          }
          
          ### 概率商归一化
          if (input[['normalized.handle.method']]=="prob_quot") {
            
            # 生成参考样本（中位数）
            ref_sample <- apply(metabo.data.fill, 2, median, na.rm = TRUE)
            
            # 计算商矩阵（每个特征值/参考样本对应特征值）
            quotient_matrix <- sweep(metabo.data.fill, 2, ref_sample, FUN = "/")
            
            # 提取归一化因子（每行商的中位数）
            norm_factors <- apply(quotient_matrix, 1, median, na.rm = TRUE)
            
            # 应用归一化因子校正数据
            metabo.data.fill <- sweep(metabo.data.fill, 1, norm_factors, FUN = "/")
            
          }
          
          wb <- createWorkbook()
          addWorksheet(wb, "归一化处理")
          writeData(wb, "归一化处理", metabo.data.fill, rowNames = TRUE)
          saveWorkbook(wb, paste0("数据预处理/归一化/归一化处理数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
        }
        
        ## 定量值取log ----
        if (input[["log.handle.method"]]=="log2") {
          
          metabo.data.fill <- log2(metabo.data.fill+1)
          
        }
       
        if (input[["log.handle.method"]]=="log10") {
          
          metabo.data.fill <- log10(metabo.data.fill+1)
          
        }
        
        ## RSD 过滤 ----
        dir.create("数据预处理/RSD过滤")
        
        if (input$log.rsd.method=="none") {
          
          metabo.data.fill <- metabo.data.fill %>% 
            dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)) %>% 
            dplyr::filter(qc_rsd <= input$rsd.cutoff)
          
        }
        
        if (input$log.rsd.method=="log2") {
          
          metabo.data.fill <- log2(metabo.data.fill+1)%>%
            dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)) %>% 
            dplyr::filter(qc_rsd <= input$rsd.cutoff)
          
        }
        
        if (input$log.rsd.method=="log10") {
          
          metabo.data.fill <- log10(metabo.data.fill+1) %>% 
            dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.fill[,grep("QC", group, ignore.case = T)]),na.rm = T)) %>% 
            dplyr::filter(qc_rsd <= input$rsd.cutoff)
          
        }
        
        wb <- createWorkbook()
        addWorksheet(wb, "RSD过滤")
        writeData(wb, "RSD过滤", metabo.data.fill, rowNames = TRUE)
        saveWorkbook(wb, paste0("数据预处理/RSD过滤/RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
        
        ## 批次矫正 ----
        if (input[['batch.correct']]=="none"){""}else{
          
          dir.create("数据预处理/批次校正")
          
          batch <- sample.data$batch
          
          ### combat校正
          if (input[['batch.correct']]=="combat") {
            
            metabo.data.fill<-sva::ComBat(metabo.data.fill%>%as.matrix(), batch = batch)
            
          }
          
          ### limma 校正
          if (input[['batch.correct']]=="limma") {
            
            # 获取实验分组变量（示例用"group"，需替换实际列名）
            group <- sample.data$group  
            
            # 构建设计矩阵
            design <- model.matrix(~group)
            
            # 数据log转换判断（示例条件，需自定义）
            if (input[['log.handle.method']]=="none") { 
              metabo.data.fill <- log2(metabo.data.fill + 1)
            }
            
            # 执行批次校正
            metabo.data.fill <- limma::removeBatchEffect(
              metabo.data.fill,
              batch = batch,
              design = design
            )
            
            # 数据还原
            if (input[['log.handle.method']]=="none") {
              metabo.data.fill <- 2^metabo.data.fill - 1
            }
            
          }
          
          wb <- createWorkbook()
          addWorksheet(wb, "批次矫正")
          writeData(wb, "批次矫正", metabo.data.fill, rowNames = TRUE)
          saveWorkbook(wb, paste0("数据预处理/批次矫正/批次矫正数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
        }
        
        assign(paste0("metabo.data.", ion_type), metabo.data.fill)
        
        # QC 质控 ----
        dir.create("QC")
        
        select_zoom<-"function (event) {
    var text,
    label;
    if (event.xAxis) {
      text = 'min: ' + Highcharts.numberFormat(event.xAxis[0].min, 2) + ', max: ' + Highcharts.numberFormat(event.xAxis[0].max, 2);
    } else {
      text = 'Selection reset';
    }
    label = this.renderer.label(text, 100, 120)
    .attr({
      fill: Highcharts.getOptions().colors[0],
      padding: 10,
      r: 5,
      zIndex: 8
    })
    .css({
      color: '#FFFFFF'
    })
    .add();
    setTimeout(function () {
      label.fadeOut();
    }, 1000);
  }"
        
        ## QC样本的TIC重叠图 ----
        dir.create("QC/TIC")
        
        QC.RT <- data.frame(EG.ApexRT=metabo.data.raw[[ion_type]]$rt[match(rownames(metabo.data.fill), metabo.data.raw[[ion_type]]$peak_name)], metabo.data.fill[,grep("QC", sample.data[[ion_type]]$Condition)]) %>% reshape2::melt(.,id=c("EG.ApexRT"),variable.name = "R.FileName",value.name = "FG.MS1Quantity")
        
        QC.RT2 <- dplyr::mutate(QC.RT,EG.ApexRT.bin=round(EG.ApexRT,digits = 0))
        
        QC.RT3 <- QC.RT2 %>% dplyr::group_by(R.FileName,EG.ApexRT.bin)%>%dplyr::summarize(Quantity.sum=sum(FG.MS1Quantity,na.rm=T))
        
        QC.RT3$type1<-"MS1"
        
        QC.TIC <- QC.RT3
        
        g_TIC_MS1 <- {
          
          highchart() %>% 
            hc_add_yAxis(lineWidth = 3,title = list(text = "Intensity"))%>%
            hc_xAxis(title = list(text = "Retention Time (min)")) %>%
            hc_add_series(data = QC.TIC[which(QC.TIC$type1=="MS1"),],type="spline",hcaes(x = EG.ApexRT.bin, y = Quantity.sum, group = R.FileName)) %>%
            hc_tooltip(split=T,valueDecimals=3)%>%
            hc_plotOptions(series = list(marker = list(symbol = "circle")))%>%
            hc_chart(zoomType="x",events=list(selection=select_zoom))%>%
            hc_exporting(enabled = TRUE,buttons = list( contextButton = list(menuItems = list('downloadPNG', 'downloadSVG',"downloadPDF","downloadJPEG","printChart","viewFullscreen"))))
          
        }
        
        htmlwidgets::saveWidget(widget = g_TIC_MS1, file = paste0("QC/TIC/TIC_", ion_type, ".html"))
        
        webshot::webshot(url = paste0("QC/TIC/TIC_", ion_type, ".html"), 
                         file = c(paste0("QC/TIC/TIC_", ion_type, ".png"),
                                  paste0("QC/TIC/TIC_", ion_type, ".pdf")),
                         delay = 3)
        
        write_csv(QC.TIC, paste0("QC/TIC/TIC_", ion_type, ".csv"))
        
        ## QC样本的相关性图 ----
        dir.create("QC/correlation")
        
        result<-corr.test(metabo.data.fill[,grep("QC", sample.data[[ion_type]]$Condition)], method = "pearson",adjust="none",alpha=.05)
        rmatrix<-result$r
        pmatrix<-result$p
        
        col <- colorRampPalette(c("darkblue", "white", "red"))(200)
        
        pdf(paste0("QC/correlation/correlation_", ion_type, ".pdf"))
        p1 <- corrplot::corrplot(rmatrix,method = "number",col = col,tl.col="black")
        dev.off()
        pdf(paste0("QC/correlation/correlation_", ion_type, ".png"))
        p1 <- corrplot::corrplot(rmatrix,method ="number",col = col,tl.col="black")
        dev.off()
        write.csv(rmatrix, paste0("QC/correlation/correlation_", ion_type, ".csv"))
        write.csv(pmatrix, paste0("QC/correlation/correlation_pvalue_", ion_type, ".csv"))
        
        ## QC样本的RSD分布图 ----
        dir.create("QC/RSD")
        
        dat <- data_to_boxplot(metabo.data.fill, qc_rsd, name="QC.RSD")
        
        g_rsd<-highchart() %>%
          hc_xAxis(type = "category") %>%
          hc_add_series_list(dat) %>%
          hc_chart(zoomType="x",events=list(selection=select_zoom))%>%
          hc_exporting(
            enabled = TRUE,
            buttons = list(
              contextButton = list(
                menuItems = list('downloadPNG', 'downloadSVG',"downloadPDF","downloadJPEG","printChart","viewFullscreen")
              ))
          )
        
        htmlwidgets::saveWidget(widget = g_rsd, file = paste0("QC/RSD/RSD_", ion_type, ".html"))
        
        webshot::webshot(url = paste0("QC/RSD/RSD_", ion_type, ".html"), 
                         file = c(paste0("QC/RSD/RSD_", ion_type, ".png"),
                                  paste0("QC/RSD/RSD_", ion_type, ".pdf")),
                         delay = 3)
        
        write.csv(metabo.data.fill$qc_rsd, paste0("QC/RSD/RSD_", ion_type, ".csv"))
        
        ## pca ----
        dir.create("QC/PCA")
        
        pca_qc <- metabo.data.fill[,-which(colnames(metabo.data.fill)=="qc_rsd")]
        
        pca_qc[is.na(pca_qc)]<-0
        
        pca_qc[pca_qc=="NaN"]<-0
        
        pca_qc[,1:ncol(pca_qc)]<-lapply(pca_qc[,1:ncol(pca_qc)],as.numeric)
        
        pc.cr <- prcomp(t(pca_qc))
        
        pca_group <- sample.data[[ion_type]]$Condition
        
        pca_group[-grep("qc", pca_group, ignore.case = T)] <- "Sample"
        
        g_pca <- {
          
          ggord::ggord(pc.cr, pca_group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, alpha=0.5, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
            geom_text(aes(label=lab),vjust=-0.5,hjust="inward",check_overlap = T,size=3,color="black",show.legend=T) +
            theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8)) +
            labs(title="PCA Plot",fill=NULL,colour=NULL,shape=NULL)
          
        } %>% ggplotly()
        
        g_pca$x$data[[3]]$showlegend<-F
        g_pca$x$data[[4]]$showlegend<-F
        
        htmlwidgets::saveWidget(widget = g_pca, file = paste0("QC/PCA/PCA_", ion_type, ".html"))
        
        webshot::webshot(url = paste0("QC/PCA/PCA_", ion_type, ".html"), 
                         file = c(paste0("QC/PCA/PCA_", ion_type, ".png"),
                                  paste0("QC/PCA/PCA_", ion_type, ".pdf")),
                         delay = 3)
        
        ## pls-da ----
        dir.create("QC/PLS-DA", recursive = TRUE, showWarnings = FALSE)
        
        # 移除qc_rsd列并处理缺失值
        plsda_qc <- metabo.data.fill[, -which(colnames(metabo.data.fill) == "qc_rsd")]
        plsda_qc[is.na(plsda_qc)] <- 0
        plsda_qc[plsda_qc == "NaN"] <- 0
        plsda_qc[, 1:ncol(plsda_qc)] <- lapply(plsda_qc[, 1:ncol(plsda_qc)], as.numeric)
        
        # 转置数据（样本在行，代谢物在列）
        X <- t(plsda_qc)
        
        # 获取分组信息（与PCA一致）
        plsda_group <- sample.data[[ion_type]]$Condition
        plsda_group[-grep("qc", plsda_group, ignore.case = T)] <- "Sample"
        
        # 运行PLS-DA模型
        plsda_result <- mixOmics::plsda(
          X = X, 
          Y = plsda_group, 
          ncomp = 2          # 选择主成分数
        )
        
        # 提取得分矩阵
        scores <- plsda_result$variates$X[, 1:2]
        
        # 将得分和分组合并为数据框
        plot_data <- data.frame(
          Comp1 = scores[, 1],
          Comp2 = scores[, 2],
          Group = plsda_group
        )%>%rownames_to_column("sample")
        
        # 绘制静态图
        g_plsda <- ggplot(plot_data, aes(x = Comp1, y = Comp2, color=Group)) +
          geom_point(size = 4) +
          stat_ellipse(level = 0.95, alpha = 0.2) +  # 添加95%置信椭圆
          geom_text_repel(aes(label=sample))+
          # scale_color_manual(values = group_colors) +
          labs(title = "PLS-DA Score Plot", 
               x = "Component 1", 
               y = "Component 2") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 12),
                text = element_text(size = 8))
        
        ggsave(paste0("QC/PLS-DA/PLS-DA_", ion_type, ".pdf"), g_plsda)
        ggsave(paste0("QC/PLS-DA/PLS-DA_", ion_type, ".png"), g_plsda)
        
        # 多元统计分析 ----
        dir.create("多元统计分析")
        
        ## pca ----
        dir.create("多元统计分析/PCA")
        
        NcrossvalI <- min(7,ncol(pca_qc))
        
        model_pca <- tryCatch(opls(t(pca_qc),scaleC = 'standard',
                                   fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI),
                              error = function(cnd) {opls(t(pca_qc),predI=3,scaleC = 'standard',
                                                          fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
                              })
        
        if(ncol(model_pca@scoreMN)==1){
          
          model_pca <- opls(t(pca_qc),predI=2,scaleC = 'standard',# center 就是UV。
                            fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
          
        }
        
        # 新版本ropls可能没有了warning和error提示。所以需要自行进行判断。
        if(ncol(model_pca@scoreMN)==0){
          model_pca <- opls(t(pca_qc),predI=3,scaleC = 'standard',# center 就是UV。
                            fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
        }
        
        df_pca <- model_pca@scoreMN %>% data.frame() %>% 
          rownames_to_column('sample') %>% 
          inner_join(sample.data[[ion_type]] %>% dplyr::select(`File name`, Condition), by=c("sample"="File name"))
        
        g <- ggpubr::ggscatter(df_pca,
                               x='p1',
                               y='p2',
                               fill = 'Condition',
                               color = 'Condition',
                               shape = 'Condition',
                               label = df_pca$sample,
                               show.legend.text = F,
                               ellipse = T)+
          # scale_color_manual(values = ngs_values[c(g1,g2) %>% sort]  %>% `+`(1))+
          # scale_fill_manual(values = ngs_values[c(g1,g2) %>% sort]  %>% `+`(1))+
          # scale_shape_manual(values = ngs_values[c(g1,g2) %>% sort] )+
          xlab(glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'))+
          ylab(glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)'))+
          ggtitle('PCA')+
          theme_bw()+
          # coord_fixed(ratio=1)+
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave(g,filename = paste0("多元统计分析/PCA/PCA_all_", ion_type, ".pdf"),width = 6.5,height = 6)
        ggsave(g,filename = paste0("多元统计分析/PCA/PCA_all_", ion_type, ".png"),width = 6.5,height = 6)
        
        # permutt_model <- matrix(ncol = 7) %>% data.frame()
        # colnames(permutt_model) <- c('title','type','A','N','R2X(cum)','R2Y(cum)','Q2(cum)')
        # permutt_model$N <- as.character(permutt_model$N)
        # 
        # permutt_model <- permutt_model %>% 
        #   add_row(title="PCA",type = 'PCA',A=length(samples),N=model_pca@summaryDF$pre %>% as.character(),
        #           `R2X(cum)`=model_pca@summaryDF$`R2X(cum)`)
        
        for (g in 1:length(group.data)) {
          
          experimental <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[1]
          control <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[2]
          
          pca_qc <- metabo.data.fill[, which(sample.data[[ion_type]]$Condition %in% c(experimental, control))]
          
          pca_qc[is.na(pca_qc)]<-0
          
          pca_qc[pca_qc=="NaN"]<-0
          
          pca_qc[,1:ncol(pca_qc)]<-lapply(pca_qc[,1:ncol(pca_qc)],as.numeric)
          
          pc.cr<-prcomp(t(pca_qc))
          
          # pca_group = sample.data$group[match(colnames(pca_qc), sample.data$sample)]
          
          NcrossvalI <- min(7,ncol(pca_qc))
          
          model_pca <- tryCatch(opls(t(pca_qc),scaleC = 'standard',
                                     fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI),
                                error = function(cnd) {opls(t(pca_qc),predI=3,scaleC = 'standard',
                                                            fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
                                })
          
          if(ncol(model_pca@scoreMN)==1){
            model_pca <- opls(t(pca_qc),predI=2,scaleC = 'standard',# center 就是UV。
                              fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
          }
          
          # 新版本ropls可能没有了warning和error提示。所以需要自行进行判断。
          if(ncol(model_pca@scoreMN)==0){
            model_pca <- opls(t(pca_qc),predI=3,scaleC = 'standard',# center 就是UV。
                              fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI)
          }
          
          df_pca <- model_pca@scoreMN %>% data.frame() %>% 
            rownames_to_column('sample') %>% 
            inner_join(sample.data[[ion_type]] %>% dplyr::select(`File name`, Condition), by=c("sample"="File name"))
          
          g1 <- ggpubr::ggscatter(df_pca,
                                 x='p1',
                                 y='p2',
                                 fill = 'Condition',
                                 color = 'Condition',
                                 shape = 'Condition',
                                 label = df_pca$sample,
                                 show.legend.text = F,
                                 ellipse = T)+
            # scale_color_manual(values = ngs_values[c(g1,g2) %>% sort]  %>% `+`(1))+
            # scale_fill_manual(values = ngs_values[c(g1,g2) %>% sort]  %>% `+`(1))+
            # scale_shape_manual(values = ngs_values[c(g1,g2) %>% sort] )+
            xlab(glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'))+
            ylab(glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)'))+
            ggtitle('PCA')+
            theme_bw()+
            # coord_fixed(ratio=1)+
            theme(plot.title = element_text(hjust = 0.5))
          
          ggsave(g1,filename = paste0('多元统计分析/PCA/PCA_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 6.5,height = 6)
          ggsave(g1,filename = paste0('多元统计分析/PCA/PCA_', experimental, "_vs_", control, "_", ion_type, ".png"))
          
        }
        
        ## pls-da/opls-da ----
        dir.create("多元统计分析/PLS-DA")
        dir.create("多元统计分析/OPLS-DA")
        
        ### 所有样本 ----
        # 移除qc_rsd列并处理缺失值
        plsda_qc <- metabo.data.fill[, -which(colnames(metabo.data.fill) == "qc_rsd")]
        plsda_qc[is.na(plsda_qc)] <- 0
        plsda_qc[plsda_qc == "NaN"] <- 0
        plsda_qc[, 1:ncol(plsda_qc)] <- lapply(plsda_qc[, 1:ncol(plsda_qc)], as.numeric)
        
        # 转置数据（样本在行，代谢物在列）
        X <- t(plsda_qc)
        
        # 获取分组信息（与PCA一致）
        plsda_group <- sample.data[[ion_type]]$Condition
        
        # 运行PLS-DA模型
        plsda_result <- mixOmics::plsda(
          X = X, 
          Y = plsda_group, 
          ncomp = 2          # 选择主成分数
        )
        
        # 提取得分矩阵
        scores <- plsda_result$variates$X[, 1:2]
        
        # 将得分和分组合并为数据框
        plot_data <- data.frame(
          Comp1 = scores[, 1],
          Comp2 = scores[, 2],
          Group = plsda_group
        )%>%rownames_to_column("sample")
        
        # 绘制静态图
        g_plsda <- ggplot(plot_data, aes(x = Comp1, y = Comp2, color = Group)) +
          geom_point(size = 4) +
          stat_ellipse(level = 0.95, alpha = 0.2) +  # 添加95%置信椭圆
          geom_text_repel(aes(label=sample))+
          # scale_color_manual(values = group_colors) +
          labs(title = "PLS-DA Score Plot", 
               x = "Component 1", 
               y = "Component 2") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5, size = 12),
                text = element_text(size = 8))
        
        ggsave(paste0("多元统计分析/PLS-DA/PLS-DA_all_", ion_type, ".pdf"), g_plsda)
        ggsave(paste0("多元统计分析/PLS-DA/PLS-DA_all_", ion_type, ".png"), g_plsda)
        
        # PLS-DA模型与置换检验
        tryCatch({
          
          plsda_model <- opls(t(plsda_qc), plsda_group, crossvalI=ifelse(ncol(plsda_qc)<7,ncol(plsda_qc),7), orthoI=0)
          
          pdf(paste0("多元统计分析/PLS-DA/PLS-DA_置换检验_", ion_type, ".pdf"), width = 5, height = 5)
          plot(plsda_model, typeVc = "permutation")  # 置换检验图
          dev.off()
          
          png(paste0("多元统计分析/PLS-DA/PLS-DA_置换检验_", ion_type, ".png"), width = 800, height = 600)
          plot(plsda_model, typeVc = "permutation")  # 置换检验图
          dev.off()
          
          pdf(paste0("多元统计分析/PLS-DA/PLS-DA_模型诊断_", ion_type, ".pdf"), width = 5, height = 5)
          plot(plsda_model, typeVc = "summary")  # 模型诊断图
          dev.off()
          
          png(paste0("多元统计分析/PLS-DA/PLS-DA_模型诊断_", ion_type, ".png"), width = 5, height = 600)
          plot(plsda_model, typeVc = "summary")  # 模型诊断图
          dev.off()
          
          opls_result_data<-plsda_model@modelDF
          colnames(opls_result_data)[which(colnames(opls_result_data)=="Signif.")]<-"Significance"
          write.csv(opls_result_data, paste0("多元统计分析/PLS-DA/PLS-DA_模型诊断_", ion_type, ".csv"))
          
        },error=function(e){
          
          ""
          
        })
        
        ### 分组样本（PLSDA和OPLS-DA） ----
        for (g in 1:length(group.data)) {
          
          experimental <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[1]
          control <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[2]
          
          plsda_qc <- metabo.data.fill[, which(sample.data[[ion_type]]$Condition %in% c(experimental, control))]
          
          plsda_qc[is.na(plsda_qc)] <- 0
          plsda_qc[plsda_qc == "NaN"] <- 0
          plsda_qc[, 1:ncol(plsda_qc)] <- lapply(plsda_qc[, 1:ncol(plsda_qc)], as.numeric)
          
          # 转置数据（样本在行，代谢物在列）
          X <- t(plsda_qc)
          
          # 获取分组信息
          plsda_group <- sample.data[[ion_type]]$Condition[match(colnames(plsda_qc), sample.data[[ion_type]]$`File name`)]
          
          # 运行PLS-DA模型
          plsda_result <- mixOmics::plsda(
            X = X, 
            Y = plsda_group, 
            ncomp = 2          # 选择主成分数
          )
          
          # 提取得分矩阵
          scores <- plsda_result$variates$X[, 1:2]
          
          # 将得分和分组合并为数据框
          plot_data <- data.frame(
            Comp1 = scores[, 1],
            Comp2 = scores[, 2],
            Group = plsda_group
          )%>%rownames_to_column("sample")
          
          # 绘制静态图
          g_plsda <- ggplot(plot_data, aes(x = Comp1, y = Comp2, color = Group)) +
            geom_point(size = 4) +
            stat_ellipse(level = 0.95, alpha = 0.2) +  # 添加95%置信椭圆
            geom_text_repel(aes(label=sample))+
            # scale_color_manual(values = group_colors) +
            labs(title = "PLS-DA Score Plot", 
                 x = "Component 1", 
                 y = "Component 2") +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 12),
                  text = element_text(size = 8))
          
          # ggsave("多元统计分析/PLS-DA/PLS-DA_all_", ion_type, ".pdf", g_plsda)
          # ggsave("多元统计分析/PLS-DA/PLS-DA_all_", ion_type, ".png", g_plsda)
          
          ggsave(g1,filename = paste0('多元统计分析/PLS-DA/PLS-DA_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 6.5,height = 6)
          ggsave(g1,filename = paste0('多元统计分析/PLS-DA/PLS-DA_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 6.5,height = 6)
          
          # PLS-DA模型与置换检验
          tryCatch({
            
            tryCatch({
              
              plsda_model <- opls(t(plsda_qc), plsda_group, crossvalI=ifelse(ncol(plsda_qc)<7,ncol(plsda_qc),7), permI = 20)
              
            },error=function(e){
              
              plsda_model <<- opls(t(plsda_qc), plsda_group, predI = 1, crossvalI=ifelse(ncol(plsda_qc)<7,ncol(plsda_qc),7), permI = 20)
              
            })
            
            pdf(paste0('多元统计分析/PLS-DA/PLS-DA_置换检验_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 5,height = 5)
            plot(plsda_model, typeVc = "permutation")  # 置换检验图
            dev.off()
            
            png(paste0('多元统计分析/PLS-DA/PLS-DA_置换检验_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 800,height = 600)
            plot(plsda_model, typeVc = "permutation")  # 置换检验图
            dev.off()
            
            pdf(paste0('多元统计分析/PLS-DA/PLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 5,height = 5)
            plot(plsda_model, typeVc = "summary")  # 模型诊断图
            dev.off()
            
            # png(paste0('多元统计分析/PLS-DA/PLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 800,height = 600)
            # plot(plsda_model, typeVc = "summary")  # 模型诊断图
            # dev.off()
            
            opls_result_data<-plsda_model@modelDF
            colnames(opls_result_data)[which(colnames(opls_result_data)=="Signif.")]<-"Significance"
            write.csv(opls_result_data,paste0('多元统计分析/PLS-DA/PLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".csv"))
            
            # OPLS-DA模型与置换检验
            plsda_model <- opls(t(plsda_qc), plsda_group, orthoI=NA, crossvalI=ifelse(ncol(plsda_qc)<7,ncol(plsda_qc),7), permI = 20)
            
            pdf(paste0('多元统计分析/OPLS-DA/OPLS-DA_置换检验_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 5,height = 5)
            plot(plsda_model, typeVc = "permutation")  # 置换检验图
            dev.off()
            
            png(paste0('多元统计分析/OPLS-DA/OPLS-DA_置换检验_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 800,height = 600)
            plot(plsda_model, typeVc = "permutation")  # 置换检验图
            dev.off()
            
            pdf(paste0('多元统计分析/OPLS-DA/OPLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 5,height = 5)
            plot(plsda_model, typeVc = "summary")  # 模型诊断图
            dev.off()
            
            # png(paste0('多元统计分析/OPLS-DA/OPLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 800,height = 600)
            # plot(plsda_model, typeVc = "summary")  # 模型诊断图
            # dev.off()
            
            opls_result_data<-plsda_model@modelDF
            colnames(opls_result_data)[which(colnames(opls_result_data)=="Signif.")]<-"Significance"
            write.csv(opls_result_data,paste0('多元统计分析/OPLS-DA/OPLS-DA_模型诊断_', experimental, "_vs_", control, "_", ion_type, ".csv"))
            
          },error=function(e){
            
            ""
            
          })
          
        }
        
        # 差异筛选 ----
        dir.create("差异筛选")
        
        dataMat <- metabo.data.fill
        
        for (g in 1:length(group.data)){
          
          cmp_name <- group.data[g]
          expGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
          ctrlGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
          
          # 创建文件夹（如果不存在）用于保存当前比较的结果
          dir.create(paste0("差异筛选/",cmp_name), recursive = TRUE, showWarnings = F)
          
          subData <- dataMat[, which(sample.data[[ion_type]]$Condition %in% c(expGroup, ctrlGroup))]
          
          # 分别提取实验组和对照组的样本名称（在subData中）
          expSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == expGroup]
          ctrlSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == ctrlGroup]
          
          # 计算每个代谢物的t检验p值（在当前比较样本上）
          pvals <- apply(subData, 1, function(x) {
            
            if (input$pvalue_type=="ttest") {
              
              t_res <- t.test(x[expSamples], x[ctrlSamples])
              
            }
            
            if (input$pvalue_type=="wilcox_test") {
              
              t_res <- wilcox.test(x[expSamples], x[ctrlSamples])
              
            }
            
            return(t_res$p.value)
            
          })
          
          # 计算Fold Change：实验组均值 / 对照组均值
          fc <- apply(subData, 1, function(x) {
            
            mean_exp <- mean(x[expSamples], na.rm = TRUE)
            mean_ctrl <- mean(x[ctrlSamples], na.rm = TRUE)
            
            if (input$log.handle.method=="none") {
              
              fc <- mean_exp / mean_ctrl
              
            }else{
              
              fc <- abs(mean_exp - mean_ctrl)
              
            }
            
            return(fc)
            
          })
          
          # 利用ropls包进行PLS-DA，计算VIP值（只取预测成分predI=1）
          selectedSamples <- c(expSamples, ctrlSamples)
          X <- t(subData[, selectedSamples])
          y <- factor(c(rep("exp", length(expSamples)), rep("ctrl", length(ctrlSamples))))
          model <- opls(X, y, predI = 1, orthoI = NA, crossvalI=ifelse(length(selectedSamples)<7,length(selectedSamples),7))
          png(paste0("差异筛选/", cmp_name, "/", cmp_name, "_plsda_", ion_type, ".png"))
          plot(model, typeVc = "x-score")  # 生成预测分类图
          dev.off()
          pdf(paste0("差异筛选/", cmp_name, "/", cmp_name, "_plsda_", ion_type, ".pdf"))
          plot(model, typeVc = "x-score")  # 生成预测分类图
          dev.off()
          vip <- getVipVn(model)
          vip <- vip[rownames(subData)]
          
          # 整合统计结果和当前比较的定量数据
          res_df <- data.frame(peak_name = rownames(subData),
                               name = metabo.data.raw[[ion_type]]$name[match(rownames(subData), metabo.data.raw[[ion_type]]$peak_name)],
                               pvalue = pvals,
                               FC = fc,
                               VIP = vip,
                               stringsAsFactors = FALSE)
          
          # 增加上下调标志
          res_df[["difference"]] <- "nodiff"
          
          if (input$fc_cutoff_judge & !input$vip_cutoff_judge) {
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC > input$fc_cutoff)] <- "up"
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC < (1/input$fc_cutoff))] <- "down"
          }
          
          if (!input$fc_cutoff_judge & input$vip_cutoff_judge) {
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC > 1 & res_df$VIP>input$vip_cutoff)] <- "up"
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC < 1 & res_df$VIP>input$vip_cutoff)] <- "down"
          }
          
          if (input$fc_cutoff_judge & input$vip_cutoff_judge) {
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC > input$fc_cutoff & res_df$VIP>input$vip_cutoff)] <- "up"
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC < (1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
          }
          
          if (!input$fc_cutoff_judge & !input$vip_cutoff_judge) {
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC > 1)] <- "up"
            res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$FC < 1)] <- "down"
          }
          
          # 匹配注释
          anno_metabo_data <- metabo.data.raw[[ion_type]][,c("peak_name", "id_kegg", "id_hmdb")]
          anno_metabo_data$peak_name <- as.character(anno_metabo_data$peak_name)
          res_df <- cbind(left_join(res_df, anno_metabo_data), subData)
          
          # 添加Mode标记
          res_df$peak_name <- paste0(res_df$peak_name, "_", ion_type)
          
          colnames(res_df) <- str_replace_all(colnames(res_df), "(POS)|(pos)|(NEG)|(neg)", "ion")
          
          results_mode[[cmp_name]][[ion_type]] <- res_df
          
          results_mode_diff[[cmp_name]][[ion_type]] <- res_df_diff <- res_df[which(res_df$difference != "nodiff"),]
          
          results_mode_diff_name[[cmp_name]][[ion_type]] <- res_df_diff_name  <- res_df[which(res_df$difference != "nodiff" & !is.na(res_df$name)),]
          
          summary_stats_list[[cmp_name]][[ion_type]] <- data.frame(compare=cmp_name,
                                                                   ion_type=ion_type,
                                                                   raw_feature_num = nrow(metabo.data),
                                                                   remove_missing_num = nrow(metabo.data2),
                                                                   preprocess_feature_num=nrow(res_df),
                                                                   diff_feature=nrow(res_df_diff),
                                                                   diff_feature_up=sum(res_df_diff$difference=="up"),
                                                                   diff_feature_down=sum(res_df_diff$difference=="down"),
                                                                   diff_feature_name=nrow(res_df_diff_name),
                                                                   diff_feature_name_up=sum(res_df_diff_name$difference=="up"),
                                                                   diff_feature_name_down=sum(res_df_diff_name$difference=="down"))
          
          # 保存结果到当前比较的文件夹
          wb <- createWorkbook()
          addWorksheet(wb, "组间对比")
          addWorksheet(wb, "组间对比差异")
          addWorksheet(wb, "组间对比差异name")
          writeData(wb, "组间对比", res_df)
          writeData(wb, "组间对比差异", res_df_diff)
          writeData(wb, "组间对比差异name", res_df_diff_name)
          saveWorkbook(wb, paste0("差异筛选/", cmp_name, "/", cmp_name, "_组间对比_", ion_type, ".xlsx"), overwrite = TRUE)
          
        }
        
      }
      
      # 各组NEG+POS合并后的生信
      # metabo.data <- rbind(metabo.data.NEG, metabo.data.POS)
      
      for (cmp_name in group.data) {
        
        experimental <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
        control <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
        
        # 合并当前比较中POS和NEG的结果（行合并，每行为一个代谢物的统计结果）
        results_mode_merged <- do.call(rbind, results_mode[[cmp_name]])
        results_mode_diff_merged <- do.call(rbind, results_mode_diff[[cmp_name]])
        results_mode_diff_name_merged <- do.call(rbind, results_mode_diff_name[[cmp_name]])
        
        expSamples <- sample.data[["POS"]]$`File name`[sample.data[["POS"]]$Condition == experimental]%>%str_replace_all("(POS)|(pos)","ion")
        ctrlSamples <- sample.data[["POS"]]$`File name`[sample.data[["POS"]]$Condition == control]%>%str_replace_all("(POS)|(pos)","ion")
        sample_cols <- c(expSamples, ctrlSamples)
        
        expr_data <- results_mode_diff_merged[, sample_cols]
        rownames(expr_data) <- results_mode_diff_merged$peak_name
        
        ### 1. HCA（层级聚类热图和表格） ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/HCA")
        dir.create(output_folder, recursive = T)
        
        # 对行（feature）进行标准化（Z-score），便于聚类比较
        expr_data_scaled <- t(scale(t(expr_data)))
        
        # 绘制热图
        pheatmap(expr_data_scaled, 
                 clustering_distance_rows = "euclidean", 
                 clustering_method = "complete", 
                 main = "Hierarchical Clustering Heatmap of Differential Features",
                 filename = file.path(output_folder, "HCA_heatmap.pdf"),
                 height = max(0.1*nrow(expr_data),30))
        
        # 计算层级聚类（行聚类）
        hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
        # 保存聚类顺序表格
        cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
        write.xlsx(cluster_order, file = file.path(output_folder, "HCA.xlsx"), rowNames = T)
        
        ### 2. Correlation Analysis ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/correlation")
        dir.create(output_folder, recursive = T)
        
        # 这里计算样本之间的相关性
        cor_mat <- cor(expr_data, use = "pairwise.complete.obs")
        
        # 绘制相关性热图
        pdf(file.path(output_folder, "correlation_sample.pdf"))
        corrplot(cor_mat, method = "color", tl.col = "black", 
                 # addCoef.col = "black",
                 col = colorRampPalette(c("blue", "white", "red"))(200),
                 title = "Correlation Among Samples", mar = c(0,0,1,0))
        dev.off()
        
        # 保存相关性矩阵
        write.xlsx(cor_mat, file = file.path(output_folder, "correlation_sample.xlsx"), rowNames = TRUE)
        
        # 这里计算feature之间的相关性
        cor_mat <- cor(t(expr_data), use = "pairwise.complete.obs")
        
        # 绘制相关性热图
        pdf(file.path(output_folder, "correlation_feature.pdf"))
        corrplot(cor_mat, method = "color", tl.col = "black", tl.cex = 0.2,
                 # addCoef.col = "black",
                 col = colorRampPalette(c("blue", "white", "red"))(200),
                 title = "Correlation Among Features", mar = c(0,0,1,0))
        dev.off()
        
        # 保存相关性矩阵
        write.xlsx(cor_mat, file = file.path(output_folder, "correlation_feature.xlsx"), rowNames = TRUE)
        
        ### 3. Volcano Plot ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/volcano")
        dir.create(output_folder, recursive = T)
        
        diff_features_all <- results_mode_merged %>%
          mutate(
            log2FC = log2(FC),
            negLog10P = -log10(pvalue),
            regulation = case_when(
              pvalue < input$p_cutoff & FC > input$fc_cutoff & VIP > input$vip_cutoff ~ "Up",
              pvalue < input$p_cutoff & FC < (1/input$fc_cutoff) & VIP > input$vip_cutoff ~ "Down",
              TRUE ~ "NS"  # 其他情况均为非显著
            )
          )
        
        # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
        max_val <- ceiling(max(abs(diff_features_all$log2FC[is.finite(diff_features_all$log2FC)]), na.rm = TRUE))
        
        volcano_plot <- ggplot(diff_features_all, aes(x = log2FC, y = negLog10P, color = regulation, size=VIP)) +
          geom_point(alpha = 0.6) +  # 统一点的大小和透明度
          scale_color_manual(
            values = c("Up" = "red", "Down" = "blue", "NS" = "grey70"),  # NS用灰色
            labels = c("Up", "Down", "Not Significant")
          ) +
          geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)), 
                     linetype = "dashed", color = "black") +
          geom_hline(yintercept = -log10(input$p_cutoff), linetype = "dashed", color = "black") +
          theme_minimal() +
          labs(x = "log2(Fold Change)", y = "-log10(p-value)", color = "Regulation") +
          scale_x_continuous(limits = c(-max_val, max_val))
        
        ggsave(filename = file.path(output_folder, "volcano.png"), volcano_plot, width = 7, height = 6)
        ggsave(filename = file.path(output_folder, "volcano.pdf"), volcano_plot, width = 7, height = 6)
        write.xlsx(diff_features_all, file.path(output_folder, "volcano.xlsx"))
        
        ## 差异的代谢物（有name的）---------------------------
        expr_data <- results_mode_diff_name_merged[,sample_cols]
        rownames(expr_data) <- paste0(results_mode_diff_name_merged$peak_name, "_", results_mode_diff_name_merged$name)
        
        ### 1. HCA（层级聚类热图和表格） ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/HCA")
        # 对行（feature）进行标准化（Z-score），便于聚类比较
        expr_data_scaled <- t(scale(t(expr_data)))
        
        # 绘制热图
        truncated_labels <- substr(
          rownames(expr_data_scaled), 
          1, 
          pmin(nchar(rownames(expr_data_scaled)), 20)  # 最多保留20个字符
        )
        # 添加省略号表示截断
        truncated_labels <- ifelse(
          nchar(rownames(expr_data_scaled)) > 20,
          paste0(truncated_labels, "..."),
          truncated_labels
        )
        
        pheatmap(
          expr_data_scaled,
          labels_row = truncated_labels,  # 使用截断后的标签
          clustering_distance_rows = "euclidean",
          clustering_method = "complete",
          main = "Hierarchical Clustering Heatmap of Differential Features",
          filename = file.path(output_folder, "HCA_heatmap_metabolite.pdf")
        )
        
        # 计算层级聚类（行聚类）
        hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
        # 保存聚类顺序表格
        cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
        write.xlsx(cluster_order, file = file.path(output_folder, "HCA_metabolite.xlsx"), rowNames = T)
        
        ### 2. Correlation Analysis ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/correlation")
        
        # 这里计算样本之间的相关性
        cor_mat <- cor(expr_data, use = "pairwise.complete.obs")
        
        # 绘制相关性热图
        pdf(file.path(output_folder, "correlation_metabolite_sample.pdf"))
        corrplot(cor_mat, method = "color", tl.col = "black", 
                 col = colorRampPalette(c("blue", "white", "red"))(200),
                 title = "Correlation Among Samples", mar = c(0,0,1,0))
        dev.off()
        
        # 保存相关性矩阵
        write.xlsx(cor_mat, file = file.path(output_folder, "correlation_metabolite_sample.xlsx"), rowNames = TRUE)
        
        # 这里计算feature之间的相关性
        cor_mat <- cor(t(expr_data), use = "pairwise.complete.obs")
        
        # 获取原始名称（假设行名是特征名）
        feature_names <- rownames(cor_mat)
        
        # 截断名称（超过20字符则保留前17字符 + ...）
        truncated_names <- ifelse(
          nchar(feature_names) > 20,
          paste0(substr(feature_names, 1, 17), "..."),
          feature_names
        )
        
        # 更新矩阵的行列名
        rownames(cor_mat) <- truncated_names
        colnames(cor_mat) <- truncated_names
        
        # 绘制相关性热图
        pdf(file.path(output_folder, "correlation_metabolite_feature.pdf"))
        corrplot(
          cor_mat,
          method = "color",
          tl.col = "black",
          col = colorRampPalette(c("blue", "white", "red"))(200),
          title = "Correlation Among Features",
          mar = c(0, 0, 1, 0)
        )
        dev.off()
        
        # 保存相关性矩阵
        write.xlsx(t(expr_data), file = file.path(output_folder, "correlation_metabolite_feature.xlsx"), rowNames = TRUE)
        
        ### 3. Volcano Plot ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/volcano")
        
        # 计算 log2(FC)
        diff_features_all <- results_mode_merged[which(!is.na(results_mode_merged$name)),] %>%
          mutate(
            log2FC = log2(FC),
            negLog10P = -log10(pvalue),
            regulation = case_when(
              pvalue < input$p_cutoff & FC > input$fc_cutoff & VIP > input$vip_cutoff ~ "Up",
              pvalue < input$p_cutoff & FC < (1/input$fc_cutoff) & VIP > input$vip_cutoff ~ "Down",
              TRUE ~ "NS"  # 其他情况均为非显著
            )
          )
        
        # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
        max_val <- ceiling(max(abs(diff_features_all$log2FC[is.finite(diff_features_all$log2FC)]), na.rm = TRUE))
        
        volcano_plot <- ggplot(diff_features_all, aes(x = log2FC, y = negLog10P, color = regulation, size=VIP)) +
          geom_point(alpha = 0.6) +  # 统一点的大小和透明度
          scale_color_manual(
            values = c("Up" = "red", "Down" = "blue", "NS" = "grey70"),  # NS用灰色
            labels = c("Up", "Down", "Not Significant")
          ) +
          geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)), 
                     linetype = "dashed", color = "black") +
          geom_hline(yintercept = -log10(input$p_cutoff), linetype = "dashed", color = "black") +
          theme_minimal() +
          labs(x = "log2(Fold Change)", y = "-log10(p-value)", color = "Regulation") +
          scale_x_continuous(limits = c(-max_val, max_val))
        
        ggsave(filename = file.path(output_folder, "volcano_metabolite.png"), volcano_plot, width = 7, height = 6)
        ggsave(filename = file.path(output_folder, "volcano_metabolite.pdf"), volcano_plot, width = 7, height = 6)
        write.xlsx(diff_features_all, file.path(output_folder, "volcano_metabolite.xlsx"))
        
        ### 4. VIP分析 ---------------------
        output_folder <- paste0("差异筛选/", cmp_name, "/VIP_analysis")
        dir.create(output_folder, recursive = T)
        
        # 绘制VIP分布直方图
        vip_hist <- ggplot(results_mode_diff_name_merged, aes(x = VIP)) + 
          geom_histogram(binwidth = 0.2, fill = "skyblue", color = "black", alpha = 0.8) +
          geom_vline(xintercept = input$vip_cutoff, linetype = "dashed", color = "red") +
          theme_minimal() +
          labs(title = "VIP Value Distribution", x = "VIP", y = "Count")
        
        ggsave(file.path(output_folder, "VIP_distribution.png"), vip_hist, width = 7, height = 6)
        ggsave(file.path(output_folder, "VIP_distribution.pdf"), vip_hist, width = 7, height = 6)
        
        # 筛选高VIP代谢物 (VIP ≥ 1.0)
        high_vip_metabolites <- results_mode_diff_name_merged %>% 
          filter(VIP >= input$vip_cutoff) %>%
          arrange(desc(VIP))
        
        write.xlsx(high_vip_metabolites, 
                   file.path(output_folder, "high_VIP_metabolite.xlsx"),
                   rowNames = FALSE)
        
        # ### 5. 差异代谢物分类 ------------------
        # output_folder <- paste0("差异筛选/", cmp_name, "/差异代谢物分类")
        # dir.create(output_folder, recursive = T)
        # 
        # # 定义分类层级列表
        # hmdb_levels <- c("kingdom", "super_class", "class", "sub_class")
        # 
        # # 遍历每个分类层级
        # for (level in hmdb_levels) {
        #   # 统计各分类的代谢物数量
        #   class_counts <- diff_features %>%
        #     dplyr::group_by(!!sym(level)) %>%
        #     dplyr::summarise(Count = n()) %>%
        #     dplyr::arrange(desc(Count)) %>%
        #     dplyr::mutate(Percentage = round(Count / sum(Count) * 100, 1))%>%.[which(!is.na(.[[1]])),]
        #   
        #   # 处理长分类名称（换行）
        #   class_counts[[level]] <- gsub(" ", "\n", class_counts[[level]])
        #   
        #   # 绘制柱状图
        #   p <- ggplot(class_counts, aes(x = reorder(!!sym(level), -Count), y = Count)) +
        #     geom_bar(stat = "identity", fill = "#4DBBD5B2", width = 0.8) +
        #     geom_text(aes(label = Count), vjust = -0.5, size = 4, color = "black") +
        #     labs(
        #       title = paste("HMDB", str_to_title(gsub("_", " ", level))),
        #       x = NULL, y = "Number of Metabolites"
        #     ) +
        #     theme_classic() +
        #     theme(
        #       axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        #       plot.title = element_text(hjust = 0.5, face = "bold"),
        #       panel.grid.major.y = element_line(color = "grey90")
        #     )
        #   
        #   # 保存结果
        #   ggsave(
        #     filename = file.path(output_folder, paste0("HMDB_", level, "_barplot.pdf")),
        #     plot = p, width = 8, height = 6
        #   )
        #   write.xlsx(
        #     class_counts,
        #     file = file.path(output_folder, paste0("HMDB_", level, "_counts.xlsx"))
        #   )
        # }
        
        ### 6.差异代谢物表达量的箱线图 ---------------
        output_folder <- paste0("差异筛选/", cmp_name, "/差异代谢物表达量")
        dir.create(output_folder, recursive = T, showWarnings = F)
        
        # 转换为长数据格式（代谢物ID + 样本 + 定量值 + 分组）
        expr_long <- results_mode_diff_name_merged[,c("name", sample_cols)] %>%
          tidyr::pivot_longer(cols = -name, 
                              names_to = "sample", 
                              values_to = "intensity") %>%
          dplyr::left_join(sample.data$POS,
                           by = c("sample" = "File name")) %>%
          dplyr::mutate(group = factor(Condition, levels = unique(sample.data$POS$Condition)))
        
        # (1) 分面箱线图（按代谢物分组）
        p_group <- ggplot(expr_long, aes(x = group, y = intensity, fill = group)) +
          geom_boxplot(alpha = 0.8, outlier.size = 1) +
          scale_fill_brewer(palette = "Set2") +
          labs(title = "代谢物定量值分布（按实验组）", x = NULL, y = "Normalized Intensity") +
          theme_classic() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 12)) 
        
        # 保存PDF（适合代谢物数量较少时）
        ggsave(file.path(output_folder, "Metabolite_Boxplot_by_Group.pdf"), p_group, width = 12, height = 8)
        ggsave(file.path(output_folder, "Metabolite_Boxplot_by_Group.png"), p_group, width = 12, height = 8)
        
        # (2) 分面箱线图（按样本）
        p_sample <- ggplot(expr_long, aes(x = sample, y = intensity, fill = group)) +
          geom_boxplot(width = 0.6) +
          scale_fill_manual(values = c("#1B9E77", "#D95F02")) + # 自定义分组颜色
          labs(title = "代谢物定量值分布（按样本）", x = "Sample", y = "Normalized Intensity") +
          theme_bw() +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12))
        
        # 保存高分辨率图片（适合样本数量多时）
        ggsave(file.path(output_folder, "Metabolite_Boxplot_by_Sample.pdf"), p_sample, width = 16, height = 10, dpi = 300)
        ggsave(file.path(output_folder, "Metabolite_Boxplot_by_Sample.png"), p_sample, width = 16, height = 10, dpi = 300)
        
        # (3) 所有代谢物全局分布
        if (length(unique(expr_long$name)) < 50) {
          
          # # 定义需要比较的组别（假设是两组比较）
          # my_comparisons <- list(unique(expr_long$group))  # 例如 list(c("Control", "Treatment"))
          
          if (input$pvalue_type=="ttest") {
            stat_compare_means_method <- "t.test"
          }else if  (input$pvalue_type=="wilcox_test") {
            stat_compare_means_method <- "wilcox.test"
          }
          
          p_global <- ggboxplot(expr_long, x = "group", y = "intensity",color = "group", palette = "npg")+
            stat_compare_means(method = stat_compare_means_method,label = "p.format")+
            labs(title = "Metabolites Comparison", x = NULL, y = "Normalized Intensity") +
            # theme_classic() +
            # 按代谢物分面
            facet_wrap(~name, scales = "free", ncol = 4)+
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 10, face = "bold"),  # 分面标签样式
              strip.background = element_blank()     # 分面标签背景
            )
          
          ggsave(
            file.path(output_folder, "Global_Metabolite_Comparison.pdf"), 
            p_global,
            width = 16, 
            height = ceiling(length(unique(expr_long$name))/4) * 3,  # 动态调整高度
            dpi = 300
          )
          
          ggsave(
            file.path(output_folder, "Global_Metabolite_Comparison.png"), 
            p_global,
            width = 16, 
            height = ceiling(length(unique(expr_long$name))/4) * 3,  # 动态调整高度
            dpi = 300
          )
          
        }
        
        ### 7.差异代谢物的网络图 ---------------
        output_folder <- paste0("差异筛选/", cmp_name, "/差异代谢物的网络图")
        dir.create(output_folder, recursive = T, showWarnings = F)
        
        # 计算代谢物间的相关系数矩阵
        cor_mat <- cor(t(expr_data), method = "spearman")  # 使用Spearman相关系数
        
        # 获取原始名称（假设行名是特征名）
        feature_names <- rownames(cor_mat)
        
        # 原始截断代码
        truncated_names <- ifelse(
          nchar(feature_names) > 20,
          paste0(substr(feature_names, 1, 17), "..."),
          feature_names
        )
        
        # 确保唯一性：添加序号
        truncated_names <- make.unique(truncated_names, sep = "-")
        rownames(cor_mat) <- truncated_names
        colnames(cor_mat) <- truncated_names
        
        cor_df <- as.data.frame(as.table(cor_mat)) %>% 
          dplyr::rename(Metabolite1 = Var1, Metabolite2 = Var2, Correlation = Freq) %>%
          dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff)  # 去除自相关
        
        diff_features <- results_mode_diff_name_merged %>%
          dplyr::mutate(
            log2FC = log2(FC),  # 计算log2转换的差异倍数
            regulation = case_when(
              FC > input$fc_cutoff & pvalue < input$p_cutoff ~ "Up",
              FC < (1/input$fc_cutoff) & pvalue < input$p_cutoff ~ "Down",
              TRUE ~ "NS"
            )
          )
        
        diff_features$regulation <- factor(diff_features$regulation, levels = c("Up", "Down"))
        
        diff_features$truncated_name <- truncated_names
        
        # 确保代谢物名称唯一
        diff_features_unique <- diff_features %>%
          dplyr::distinct(peak_name, .keep_all = TRUE)  # 去重
        
        # 创建igraph图对象
        g <- graph_from_data_frame(
          d = cor_df[, 1:3],  # 边数据（Metabolite1, Metabolite2, Correlation）
          directed = T,
          vertices = diff_features_unique %>% dplyr::select(truncated_name, log2FC, regulation)  # 节点属性
        )
        
        # 节点属性设置
        V(g)$size <- abs(V(g)$log2FC) * 3  # 大小反映log2FC绝对值
        V(g)$color <- ifelse(V(g)$regulation == "Up", "#E64B35FF", "#3C5488FF")  # 上调红，下调蓝
        
        # 边属性设置
        E(g)$width <- abs(E(g)$Correlation) * 0.2  # 边宽与相关性强度正相关
        E(g)$color <- ifelse(E(g)$Correlation > 0, "#D73027", "#4575B4")  # 正相关红，负相关蓝
        
        # 生成带图例的网络图
        set.seed(123)
        network_plot <- ggraph(g, layout = "fr") + 
          # 边图层（颜色和宽度）
          geom_edge_link(
            aes(edge_width = abs(Correlation),  # 边宽度映射到相关性绝对值
                edge_color = ifelse(Correlation > 0, "Positive", "Negative")),  # 边颜色映射到正负
            alpha = 0.6
          ) +
          # 节点图层（颜色和大小）
          geom_node_point(
            aes(color = regulation,    # 颜色映射到上下调
                size = abs(log2FC)),   # 大小映射到log2FC绝对值
            alpha = 0.8
          ) +
          # 文本标签（防重叠）
          geom_node_text(
            aes(label = name), 
            size = 3, 
            repel = TRUE,
            color = "black"
          ) +
          # 边颜色映射（正红负蓝）
          scale_edge_color_manual(
            name = "Correlation Type",
            values = c("Positive" = "#D73027", "Negative" = "#4575B4"),
            labels = c("Positive", "Negative")
          ) +
          # 边宽度映射（连续型）
          scale_edge_width_continuous(
            name = "|Correlation|", 
            range = c(0.5, 4),  # 边宽度范围
            breaks = c(0.2, 0.5, 0.8)  # 图例刻度
          ) +
          # 节点颜色映射（离散型）
          scale_color_manual(
            name = "Regulation",
            values = c("Up" = "#E64B35FF", "Down" = "#3C5488FF"),
            labels = c("Up", "Down")
          ) +
          # 节点大小映射（连续型）
          scale_size_continuous(
            name = "|log2FC|", 
            range = c(3, 10),  # 节点大小范围
            breaks = c(1, 2, 3)  # 图例刻度
          ) +
          theme_void() +
          labs(title = "差异代谢物相关性网络") +
          # 图例位置与样式调整
          theme(
            legend.position = "right",  # 图例在右侧
            legend.box = "vertical",    # 图例垂直排列
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        
        ggsave(file.path(output_folder, "Network_CytoscapeStyle.pdf"), network_plot, width = 12, height = 10)
        ggsave(file.path(output_folder, "Network_CytoscapeStyle.png"), network_plot, width = 12, height = 10, dpi = 300)
        
        # 导出网络数据供Cytoscape进一步编辑
        write.graph(g, file.path(output_folder, "metabolite_network.graphml"), format = "graphml")
        
        ### KEGG 富集 ---------------
        output_folder <- paste0("差异筛选/", cmp_name, "/KEGG富集")
        dir.create(output_folder, recursive = T, showWarnings = F)
        
        # 假设diff_features包含kegg_id和log2FC列
        # 提取差异代谢物KEGG信息
        kegg_data <- diff_features %>% separate_rows(id_kegg)%>%
          dplyr::filter(!is.na(name)) %>%
          dplyr::distinct(id_kegg, .keep_all = TRUE) %>%
          dplyr::filter(!is.na(id_kegg))%>%.[which(.[["id_kegg"]]!="NA"),]
        
        mb3 <- enrichMBKEGG(kegg_data$id_kegg)
        
        mb3@result$Description <- factor(mb3@result$Description,levels = mb3@result$Description)
        
        p_dot <- ggplot(mb3@result, 
                        aes(x = -log10(p.adjust),
                            y = reorder(Description, -log10(p.adjust)))) +  # 排序
          geom_point(aes(size = Count, color = FoldEnrichment), alpha = 0.8) +
          scale_color_gradient2(low = "blue", mid = "grey", high = "red", 
                                midpoint = median(mb3@result$FoldEnrichment),
                                name = "Fold Enrichment") +
          scale_size(range = c(3, 10), name = "Metabolite Count") +
          labs(x = "-log10(Adj. p-value)", y = "",
               title = "KEGG Pathway Enrichment (Sorted by Significance)") +  # 标题注明排序依据
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 12, hjust = 0, color = "black"),
            plot.margin = margin(l = 2.5, unit = "cm")
          )
        
        ggsave(file.path(output_folder, "KEGG_dotplot.pdf"), p_dot, width = 10, height = 8)
        ggsave(file.path(output_folder, "KEGG_dotplot.png"), p_dot, width = 10, height = 8)
        
        write.xlsx(mb3@result, file.path(output_folder, "KEGG_enrichment.xlsx"))
        
        cat("已保存", cmp_name, "的结果至文件。\n")
        
        # 统计各阶段的feature数量及差异筛选的上下调数量
        diff_count <- nrow(results_mode_diff_merged)
        diff_up <- sum(results_mode_diff_merged$difference == "up", na.rm = TRUE)
        diff_down <- sum(results_mode_diff_merged$difference == "down", na.rm = TRUE)
        diff_name_count <- nrow(results_mode_diff_name_merged)
        diff_name_up <- sum(results_mode_diff_name_merged$difference == "up", na.rm = TRUE)
        diff_name_down <- sum(results_mode_diff_name_merged$difference == "down", na.rm = TRUE)
        
        # 将当前比较的统计信息保存到列表中
        summary_stats_list[[cmp_name]][["Merge"]] <- data.frame(
          compare = cmp_name,
          ion_type = "POS&NEG",
          raw_feature_num = summary_stats_list[[cmp_name]][["POS"]]$raw_feature_num+summary_stats_list[[cmp_name]][["NEG"]]$raw_feature_num,
          remove_missing_num = summary_stats_list[[cmp_name]][["POS"]]$remove_missing_num+summary_stats_list[[cmp_name]][["NEG"]]$remove_missing_num,
          preprocess_feature_num = nrow(results_mode_merged),
          diff_feature = diff_count,
          diff_feature_up = diff_up,
          diff_feature_down = diff_down,
          diff_feature_name = diff_name_count,
          diff_feature_name_up = diff_name_up,
          diff_feature_name_down = diff_name_down,
          stringsAsFactors = FALSE
        )
        
      }
      
      # 合并所有分组比较的统计信息
      final_summary <- do.call(rbind, do.call(rbind, summary_stats_list))
      # 保存最终统计结果到 "差异筛选/差异数量统计.xlsx"
      write.xlsx(final_summary, file = "差异筛选/差异数量统计.xlsx", rowNames = FALSE)
      cat("已保存差异数量统计结果至文件：差异筛选/差异数量统计.xlsx\n")
      
      ## feature --------------------
      ### Venn ------------------
      output_folder <- paste0("差异筛选/venn")
      dir.create(output_folder, recursive = T)
      
      diff_features_list <- list()
      
      # 遍历所有比较组
      for(cmp_name in group.data) {
        
        diff_features_list[[cmp_name]] <- do.call(rbind, results_mode_diff[[cmp_name]])%>%.[,c("peak_name")]
        
      }
      
      # 绘制Venn图（当比较组数量2-5时适用）
      if(length(diff_features_list) >= 2 && length(diff_features_list) <= 5) {
        
        # 动态选择调色板（支持2-5组）
        if(length(diff_features_list) == 2) {
          fill_colors <- brewer.pal(3, "Paired")[1:2]  # 使用Paired调色板前两色
        } else {
          fill_colors <- brewer.pal(length(diff_features_list), "Set2")  # Set2支持3-8色
        }
        
        venn.plot <- venn.diagram(
          x = diff_features_list,
          filename = file.path(output_folder, "Venn_diagram.png"),
          imagetype = "png",
          fill = fill_colors,  # 修正后的颜色参数
          alpha = 0.6,
          cat.cex = 1.2,
          cat.fontface = "bold",
          margin = 0.1,
          cex = 1.5,
          fontfamily = "sans",
          cat.fontfamily = "sans"
        )
        
        ### Upset ------------------------
        output_folder <- paste0("差异筛选/upset")
        dir.create(output_folder, recursive = T)
        
        # 创建二进制矩阵
        all_features <- unique(unlist(diff_features_list))
        
        binary_matrix <- data.frame(
          feature = all_features,
          stringsAsFactors = FALSE
        )
        
        # 为每个比较组生成存在性标记
        for(cmp in names(diff_features_list)) {
          binary_matrix[[cmp]] <- as.integer(binary_matrix$feature %in% diff_features_list[[cmp]])
        }
        
        # 生成Upset图（带交互式HTML输出）
        upset_plot <- upset(
          binary_matrix,
          nsets = length(diff_features_list),
          nintersects = 20,
          sets = names(diff_features_list),
          mainbar.y.label = "Feature Intersections",
          sets.x.label = "Differential Features per Comparison",
          text.scale = 2
        )
        
        # 保存图像
        pdf(file.path(output_folder, "UpSet_plot.pdf"), width = 20, height = 12)
        print(upset_plot)
        dev.off()
        
        png(file.path(output_folder, "UpSet_plot.png"), width = 800, height = 600)
        print(upset_plot)
        dev.off()
        
        # 保存二进制矩阵
        write.xlsx(binary_matrix, file.path(output_folder, "UpSet_data.xlsx"))
        
      }
      
      ## metabolites ----------------
      ### Venn for names ---------------------
      output_folder_name <- paste0("差异筛选/venn")
      
      # 收集所有比较组的差异代谢物名称（基于name）
      diff_names_list <- list()
      
      # 遍历所有比较组
      for(cmp_name in group.data) {
        
        diff_names_list[[cmp_name]] <- do.call(rbind, results_mode_diff_name[[cmp_name]])%>%.[,c("peak_name")]
        
      }
      
      # 绘制Venn图（当比较组数量2-5时适用）
      if(length(diff_names_list) >= 2 && length(diff_names_list) <= 5) {
        
        # 动态选择调色板
        if(length(diff_names_list) == 2) {
          
          fill_colors <- brewer.pal(3, "Paired")[1:2]
          
        } else {
          
          fill_colors <- brewer.pal(length(diff_names_list), "Set2")
          
        }
        
        venn.plot <- venn.diagram(
          x = diff_names_list,
          filename = file.path(output_folder_name, "Venn_diagram_metabolites.png"),
          imagetype = "png",
          fill = fill_colors,
          alpha = 0.6,
          cat.cex = 1.2,
          cat.fontface = "bold",
          margin = 0.1,
          cex = 1.5,
          fontfamily = "sans",
          cat.fontfamily = "sans"
        )
        
        ### Upset for names ------------------------
        output_folder_name <- paste0("差异筛选/upset")
        
        # 创建二进制矩阵
        all_names <- unique(unlist(diff_names_list))
        
        binary_matrix_name <- data.frame(
          name = all_names,
          stringsAsFactors = FALSE
        )
        
        for(cmp in names(diff_names_list)) {
          
          binary_matrix_name[[cmp]] <- as.integer(binary_matrix_name$name %in% diff_names_list[[cmp]])
          
        }
        
        # 生成Upset图
        upset_plot_name <- upset(
          binary_matrix_name,
          nsets = length(diff_names_list),
          nintersects = 20,
          sets = names(diff_names_list),
          mainbar.y.label = "Metabolite Name Intersections",
          sets.x.label = "Differential Names per Comparison",
          text.scale = 2
        )
        
        pdf(file.path(output_folder_name, "UpSet_plot_metabolites.pdf"), width = 20, height = 12)
        print(upset_plot_name)
        dev.off()
        
        png(file.path(output_folder_name, "UpSet_plot_metabolites.png"), width = 800, height = 600)
        print(upset_plot_name)
        dev.off()
        
        write.xlsx(binary_matrix_name, file.path(output_folder_name, "UpSet_data_metabolites.xlsx"))
        
      }
      
      # 出报告 -----------------
      render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report.Rmd"), output_file = paste0(path_prefix, 'data/usrdata/', project_id, '/report/report.html'))
      
    }
    
    # 打包和通知 ----
    files_to_extract_path <- list.dirs(report_path, full.names = TRUE,recursive = F)
    
    zip::zip(zipfile = paste0(report_path, '/report.zip'), files = files_to_extract_path, include_directories = TRUE, mode = "cherry-pick")
    
    UpdateStatus(2,project_id)
    
    #发送邮件
    mail <- input$email
    
    message <- paste0('From: "Omicsolution" <marketing@omicsolution.com>
To: "',stringr::str_replace(mail,"\\@.*",""),'" <',mail,'>
Subject: Task finished inform!

Dear user,

Your project (', input$project_name, ') have finished.Please go to http://192.168.100.198:3838/BioAnalysis/?demo&id2=',project_id,' to check the result.')
    
    curl::send_mail(mail_from = "marketing@omicsolution.com",
                    mail_rcpt = mail, 
                    message, 
                    smtp_server = 'smtps://smtp.exmail.qq.com',
                    username = 'marketing@omicsolution.com', 
                    password  = 'ERPUuVncXNJjc8gT')
    
  },error=function(e){
    
    #发送邮件
    mail <- input$email
    
    message <- paste0('From: "Omicsolution" <marketing@omicsolution.com>
To: "',stringr::str_replace(mail,"\\@.*",""),'" <',mail,'>
Subject: Task finished inform!

Dear user,

Your project (', input$project_name, ') have failure.',e,' Please check the raw data and group file.')
    
    curl::send_mail(mail_from = "marketing@omicsolution.com",
                    mail_rcpt = mail, 
                    message, 
                    smtp_server = 'smtps://smtp.exmail.qq.com',
                    username = 'marketing@omicsolution.com', 
                    password  = 'ERPUuVncXNJjc8gT')
    
  })
  
}else{
  
  UpdateStatus(3,project_id)
  
}