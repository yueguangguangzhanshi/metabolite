#!/usr/bin/env Rscript
library(DBI)

Sys.setenv(RSTUDIO_PANDOC="/usr/bin/pandoc")
Sys.setenv(OPENSSL_CONF="/dev/null")

path_prefix <- '/home/saitoasuka/shiny-server/BioAnalysis'
path_prefix2<-"/home/saitoasuka/shiny-server/BioAnalysis"

source(paste0(path_prefix, '/etc/dep/bg_lib.r'), local = TRUE, echo = FALSE)

UpdateStatus <- function(status,id) {
  
  # system(paste0("cp -rf ", path_prefix, "/db/projects_list.db ", path_prefix2, "/db/projects_list.db"))
  
  projects_list <- dbConnect(RSQLite::SQLite(), paste0(path_prefix2, "/db/projects_list.db"))
  
  dbSendStatement(projects_list, paste0("UPDATE Project SET status = '",status,"' WHERE id='",id,"'"))
  
  dbDisconnect(projects_list)
  
  # system(paste0("cp -rf ", path_prefix2, "/db/projects_list.db ", path_prefix, "/db/projects_list.db"))
  
}

projects_list <- dbConnect(RSQLite::SQLite(), paste0(path_prefix, "/db/projects_list.db"))

project_status_table <- as.data.frame(dbGetQuery(projects_list, paste0("SELECT * FROM Project WHERE status=0")))

dbDisconnect(projects_list)

if (nrow(project_status_table)>0) {
  
  project_id <- project_status_table[which(project_status_table$time==min(project_status_table$time)),"id"]
  
  UpdateStatus(4,project_id)
  
  metr_pkgs<-c('statTarget', 'fs', 'MetNormalizer', 'visNetwork', 'MicrobiomeProfiler', 'pathview', 'MetaboAnalystR', 'ggraph', 'igraph', 'grid', 'UpSetR', 'VennDiagram', 'corrplot', 'glue' ,'matrixStats','missForest','VIM','ggord','survminer','survival','highcharter','ropls','dplyr','reshape2','statnet','circlize','tibble','corrplot','psych','ComplexHeatmap','org.Hs.eg.db','limma','WebGestaltR','ggraph','igraph','plyr',"qvalue","digest", "BiocManager", "Cairo","ggrepel","rgl","knitr","mixOmics","pheatmap","RColorBrewer","gplots","ggplot2","openxlsx","stringr","scales","data.table","factoextra","FactoMineR","XML","downloader","ggsci","rJava","rChoiceDialogs","foreach","grid","seqinr","reshape","ontologyIndex","qvalue","VennDiagram","RSQLite","readr","tidyverse","readxl","gdata","svDialogs","ggseqlogo","rmotifx","rvest","rlist", "knitr", "rmarkdown", "bsselectR",'echarts4r', 'shiny', 'DT', 'plotly', 'ggbiplot','purrr')
  
  for(i in 1:length(metr_pkgs)){
    
    library(metr_pkgs[i], character.only = TRUE)
    
  }
  
  data_path <- paste0(path_prefix, '/data/usrdata/', project_id, '/data')
  
  report_path <- paste0(path_prefix, '/data/usrdata/', project_id, '/report')
  
  url_path <- paste0("https://project.omicsolution.com/reports/", project_id, '/report')
  
  para_path <- paste0(data_path, '/para.Rdata')
  
  if(file.exists(para_path)){
    
    load(para_path)
    
    input <- readRDS(paste0(data_path,"/input.rds"))
    
    rv <- readRDS(paste0(data_path,"/rv.rds"))
    
    tryCatch({
      
      if(para[["info"]][["type"]] == 'sjtu'){
        
        dir.create(report_path, recursive = T)
        
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
        
        if (input$single_ion_judge) {
          
          metabo.file <- paste0(data_path, "/metabonomics.csv")
          
          ion_type_all <- "single"
          
          names(metabo.file) <- ion_type_all
          
          sample.data <- list(single = rv$condition.data.single)
          
        } else {
          
          metabo.file <- c(paste0(data_path, "/metabonomics_NEG.csv"), paste0(data_path, "/metabonomics_POS.csv"))
          
          ion_type_all <- c("NEG","POS")
          
          names(metabo.file) <- ion_type_all
          
          sample.data <- list(POS = rv$condition.data.pos, NEG = rv$condition.data.neg)
          
        }
        
        if (input[['batch.correct']]!="none") {
          
          sample.data$batch <- rv$condition.data.pos[[grep("batch",colnames(rv$condition.data.pos),value = T,ignore.case = T)]]
          
        }
        
        group.data <- input$metaboQuant_compare_select
        
        if (!is.na(input$key_column_kegg))  if (input$key_column_kegg  == 0) input$key_column_kegg  <- NA
        if (!is.na(input$key_column_hmdb))  if (input$key_column_hmdb  == 0) input$key_column_hmdb  <- NA
        if (!is.na(input$key_column_lmsd))  if (input$key_column_lmsd  == 0) input$key_column_lmsd  <- NA
        if (!is.na(input$key_column_name))  if (input$key_column_name  == 0) input$key_column_name  <- NA
        
        if (!is.na(input$key_column_name)) {
          
          if (is.na(input$key_column_kegg)) {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_name,
              input$key_column_sample_start:input$key_column_sample_end
            )
          } else {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_name,
              input$key_column_kegg,
              input$key_column_sample_start:input$key_column_sample_end
            )
          }
          
        }else{
          
          if (is.na(input$key_column_kegg)) {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_sample_start:input$key_column_sample_end
            )
          } else {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_kegg,
              input$key_column_sample_start:input$key_column_sample_end
            )
          }
          
        }
        
        # 存储当前比较中POS和NEG的结果
        metabo.data.raw.all.col <- summary_stats_list <- results_mode_diff_name <- results_mode_diff <- results_mode <- metabo.data.raw <- list()
        
        sum_pca <- list(); sum_plsda <- list(); sum_oplsda <- list()
        
        # kegg 数据库
        load(file=paste0(path_prefix, '/db/kegg_id_db.rda'))
        
        if (is.na(input$key_column_kegg)) {
          
          ## 构建 name → KEGG ID 索引
          hmdb2 <- tidyr::separate_rows(hmdb, synonyms2, sep = "; ")
          hmdb2[["synonyms_tolower"]] <- tolower(hmdb2$synonyms2)
          hmdb2 <- dplyr::distinct(hmdb2, synonyms_tolower, kegg_id, .keep_all = T)
          kegg2 <- tidyr::separate_rows(kegg, name, sep = "; ")
          kegg2[["name_tolower"]] <- tolower(kegg2$name)
          kegg2 <- dplyr::distinct(kegg2, name_tolower, kegg_id, .keep_all = T)
          # 预处理：过滤掉kegg_id为NA的行
          hmdb2 <- hmdb2[!is.na(hmdb2$kegg_id), ]
          kegg2 <- kegg2[!is.na(kegg2$kegg_id), ]
          # 创建快速查询索引
          hmdb_index <- split(hmdb2$kegg_id, hmdb2$synonyms_tolower)
          kegg_index <- split(kegg2$kegg_id, kegg2$name_tolower)
          
          lmsd_name_tbl <- lmsd %>%
            dplyr::select(LM_ID, KEGG_ID, NAME, SYSTEMATIC_NAME) %>%
            tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
            dplyr::filter(!is.na(nm), nm != "") %>%
            dplyr::mutate(nm_tolower = tolower(nm)) %>%
            dplyr::filter(!is.na(KEGG_ID), KEGG_ID != "") %>%
            dplyr::distinct(nm_tolower, KEGG_ID)
          lmsd_index <- split(lmsd_name_tbl$KEGG_ID, lmsd_name_tbl$nm_tolower)
          
          ## 构建 HMDB accession → KEGG ID 索引
          hmdb3 <- dplyr::distinct(hmdb2, accession, kegg_id, .keep_all = T)
          accession_index <- split(hmdb3$kegg_id, hmdb3$accession)
          
          ## 构建 LM_ID → KEGG ID 索引
          lmsd_lm_tbl <- lmsd %>%
            dplyr::filter(!is.na(LM_ID), LM_ID != "", !is.na(KEGG_ID), KEGG_ID != "") %>%
            dplyr::distinct(LM_ID, KEGG_ID)
          lmsd_lm_index <- split(lmsd_lm_tbl$KEGG_ID, lmsd_lm_tbl$LM_ID)
          
        }
        
        # 数据预处理与各组的差异筛选 ----
        for (ion_type in ion_type_all) {
          
          dir.create("1.搜库原始结果")
          
          # 读取数据
          raw_data <- metabo.data.raw.all.col[[ion_type]] <- readr::read_csv(metabo.file[[ion_type]])  # 或 read.csv()
          
          # 构造想保留的列的索引（按顺序）
          selected_columns_index <- c(
            input$key_column_peakname,
            input$key_column_mz,
            input$key_column_rt,
            input$key_column_name,
            which(colnames(raw_data) == "formula"),
            which(colnames(raw_data) == "smiles"),
            which(colnames(raw_data) == "inchikey"),
            which(colnames(raw_data) == "adduct"),
            input$key_column_kegg,
            input$key_column_hmdb,
            input$key_column_lmsd,
            input$key_column_sample_start:input$key_column_sample_end
          )
          
          # 去除不存在的列索引（防止有些列缺失）
          selected_columns_index <- selected_columns_index[selected_columns_index %in% seq_along(colnames(raw_data))]
          
          # 重新排序并提取列
          raw_data <- raw_data[, selected_columns_index]
          
          # 写入文件（可按需用 write_tsv 或 write.csv）
          readr::write_csv(raw_data, file = file.path("1.搜库原始结果", basename(metabo.file[[ion_type]])))
          
          metabo.data.raw[[ion_type]] <- metabo.data.raw.all.col[[ion_type]]%>%.[, selected_column]
          
          if (!is.na(input$key_column_name)) {
            
            if (is.na(input$key_column_kegg)) {
              
              colnames(metabo.data.raw[[ion_type]])[1:4] <- c("peak_name", "mz", "rt", "name")
              
              # 获取待匹配的名称向量
              query_names <- tolower(metabo.data.raw[[ion_type]][["name"]])
              
              matched_kegg <- lapply(query_names, function(name_str) {
                sub_names <- unlist(strsplit(name_str, ";"))
                sub_ids <- sapply(sub_names, function(name) {
                  name <- trimws(name)
                  
                  hmdb_ids <- if (!is.null(hmdb_index[[name]])) unname(hmdb_index[[name]]) else character(0)
                  kegg_ids <- if (!is.null(kegg_index[[name]])) unname(kegg_index[[name]]) else character(0)
                  lmsd_ids <- if (!is.null(lmsd_index[[name]])) unname(lmsd_index[[name]]) else character(0)
                  all_ids  <- unique(c(hmdb_ids, kegg_ids, lmsd_ids))
                  
                  # 去括号再试一次
                  if (length(all_ids) == 0 && grepl("\\([^()]*\\)$", name)) {
                    name2 <- sub("\\s*\\([^()]*\\)$", "", name)
                    hmdb2 <- if (!is.null(hmdb_index[[name2]])) unname(hmdb_index[[name2]]) else character(0)
                    kegg2 <- if (!is.null(kegg_index[[name2]])) unname(kegg_index[[name2]]) else character(0)
                    lmsd2 <- if (!is.null(lmsd_index[[name2]])) unname(lmsd_index[[name2]]) else character(0)
                    all_ids <- unique(c(hmdb2, kegg2, lmsd2))
                  }
                  
                  if (length(all_ids) > 0) paste(all_ids, collapse = "/") else "NA"
                }, USE.NAMES = FALSE)
                
                paste(sub_ids, collapse = ";")
              })
              
              # 4.2 HMDB accession → KEGG，生成 hmdb_kegg（若提供 HMDB 列）
              if (!is.na(input$key_column_hmdb)) {
                hmdb_kegg <- vapply(seq_len(nrow(metabo.data.raw.all.col[[ion_type]])), function(i) {
                  acc <- metabo.data.raw.all.col[[ion_type]][[input$key_column_hmdb]][i]
                  if (is.null(acc) || is.na(acc) || acc == "") return(NA_character_)
                  acc_list <- strsplit(toupper(as.character(acc)), ";")[[1]]
                  
                  ids <- sapply(acc_list, function(a) {
                    a <- trimws(a)
                    res <- accession_index[[a]]
                    if (!is.null(res) && length(res) > 0) paste(unique(res), collapse = "/") else NA_character_
                  }, USE.NAMES = FALSE)
                  
                  paste(ids, collapse = ";")
                }, character(1))
              }
              
              # 4.3 LM_ID → KEGG，生成 lmsd_kegg（若提供 LMSD 列）
              if (!is.na(input$key_column_lmsd)) {
                lm_col <- metabo.data.raw.all.col[[ion_type]][[ input$key_column_lmsd ]]
                lmsd_kegg <- vapply(seq_along(lm_col), function(i) {
                  val <- lm_col[i]
                  if (is.null(val) || is.na(val) || val == "") return(NA_character_)
                  lm_list <- strsplit(as.character(val), ";")[[1]]
                  ids <- sapply(lm_list, function(lm) {
                    key <- trimws(lm)
                    res <- lmsd_lm_index[[key]]
                    if (!is.null(res) && length(res) > 0) paste(unique(res), collapse = "/") else NA_character_
                  }, USE.NAMES = FALSE)
                  paste(ids, collapse = ";")
                }, character(1))
              }
              
              # 4.4 组装输出（按需附加）
              tmp <- metabo.data.raw[[ion_type]][, 1:4]
              tmp$id_kegg <- unlist(matched_kegg)
              
              if (!is.na(input$key_column_hmdb)) tmp$hmdb_kegg <- hmdb_kegg
              if (!is.na(input$key_column_lmsd)) tmp$lmsd_kegg <- lmsd_kegg
              
              metabo.data.raw[[ion_type]] <- data.frame(
                tmp,
                metabo.data.raw[[ion_type]][, -c(1:4)],
                check.names = FALSE
              )
              
            }else{
              
              colnames(metabo.data.raw[[ion_type]])[1:5] <- c("peak_name", "mz", "rt", "name", "id_kegg")
              
            }
            
          }else{
            
            # ========= 当没有 name 列时：先用 KEGG/HMDB/LMSD 的 ID 反查 name，再用 name 反查/补全 id_kegg =========
            message("未提供 name 列，基于 KEGG/HMDB/LMSD ID 自动回填 name，并补全 id_kegg...")
            
            # —— 1) 构建 ID -> name 的索引（不 reload，直接用上面已存在的 kegg/hmdb/lmsd 对象）——
            kegg_name_index <- kegg %>%
              dplyr::select(kegg_id, name) %>%
              dplyr::filter(!is.na(kegg_id), kegg_id != "", !is.na(name), name != "") %>%
              dplyr::distinct(kegg_id, .keep_all = TRUE) %>%
              tibble::deframe()
            
            hmdb_name_index <- hmdb %>%
              dplyr::select(accession, name) %>%
              dplyr::filter(!is.na(accession), accession != "", !is.na(name), name != "") %>%
              dplyr::distinct(accession, .keep_all = TRUE) %>%
              tibble::deframe()
            
            lmsd_name_index <- lmsd %>%
              dplyr::select(LM_ID, NAME) %>%
              dplyr::filter(!is.na(LM_ID), LM_ID != "", !is.na(NAME), NAME != "") %>%
              dplyr::distinct(LM_ID, .keep_all = TRUE) %>%
              tibble::deframe()
            
            # —— 2) 构建 name -> KEGG_ID 的索引（用于后续用 name 反推 id_kegg）——
            # 2.1 KEGG: name 可能是以分号分隔的多别名
            kegg2 <- tidyr::separate_rows(kegg, name, sep = "; ") %>%
              dplyr::mutate(name_tolower = tolower(name)) %>%
              dplyr::filter(!is.na(name_tolower), name_tolower != "", !is.na(kegg_id), kegg_id != "") %>%
              dplyr::distinct(name_tolower, kegg_id)
            
            kegg_index_by_name <- split(kegg2$kegg_id, kegg2$name_tolower)
            
            # 2.2 HMDB: 用 synonyms2 里的别名对 KEGG_ID 做映射（通过 hmdb$kegg_id）
            hmdb2 <- tidyr::separate_rows(hmdb, synonyms2, sep = "; ") %>%
              dplyr::mutate(synonyms_tolower = tolower(synonyms2)) %>%
              dplyr::filter(!is.na(synonyms_tolower), synonyms_tolower != "", !is.na(kegg_id), kegg_id != "") %>%
              dplyr::distinct(synonyms_tolower, kegg_id)
            
            hmdb_index_by_name <- split(hmdb2$kegg_id, hmdb2$synonyms_tolower)
            
            # 2.3 LMSD: NAME/Systematic name -> KEGG_ID
            lmsd_name_tbl <- lmsd %>%
              dplyr::select(LM_ID, KEGG_ID, NAME, SYSTEMATIC_NAME) %>%
              tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
              dplyr::filter(!is.na(nm), nm != "", !is.na(KEGG_ID), KEGG_ID != "") %>%
              dplyr::mutate(nm_tolower = tolower(nm)) %>%
              dplyr::distinct(nm_tolower, KEGG_ID)
            
            lmsd_index_by_name <- split(lmsd_name_tbl$KEGG_ID, lmsd_name_tbl$nm_tolower)
            
            # —— 3) 用 ID 反查得到 name —— 
            n_row <- nrow(metabo.data.raw.all.col[[ion_type]])
            
            name_from_kegg <- name_from_hmdb <- name_from_lmsd <- rep(NA_character_, n_row)
            
            # KEGG id -> name
            if (!is.na(input$key_column_kegg)) {
              name_from_kegg <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_kegg]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(kegg_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # HMDB accession -> name
            if (!is.na(input$key_column_hmdb)) {
              name_from_hmdb <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_hmdb]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(hmdb_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # LMSD LM_ID -> name
            if (!is.na(input$key_column_lmsd)) {
              name_from_lmsd <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_lmsd]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(lmsd_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # 优先级：KEGG > HMDB > LMSD
            filled_name <- dplyr::coalesce(name_from_kegg, name_from_hmdb, name_from_lmsd)
            
            # —— 4) 在得到 filled_name 的基础上，反查/补全 KEGG_ID —— 
            derive_kegg_from_name <- function(nm_str) {
              if (is.na(nm_str) || nm_str == "") return(NA_character_)
              toks <- strsplit(tolower(nm_str), "[;|,]")[[1]] %>% trimws()
              ids <- unique(c(
                unlist(kegg_index_by_name[toks]),
                unlist(hmdb_index_by_name[toks]),
                unlist(lmsd_index_by_name[toks])
              ))
              ids <- ids[!is.na(ids) & ids != ""]
              if (length(ids) == 0) return(NA_character_)
              paste(unique(ids), collapse = "/")     # 多命中用 / 连接
            }
            
            derived_kegg_by_name <- vapply(filled_name, derive_kegg_from_name, FUN.VALUE = character(1))
            
            # —— 5) 组装输出：插入 name，且创建/补全 id_kegg —— 
            #   目标列顺序：peak_name, mz, rt, name, id_kegg, <sample cols...>
            base_cols <- metabo.data.raw.all.col[[ion_type]][, c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt
            )]
            
            colnames(base_cols) <- c("peak_name","mz","rt")
            
            # 如果原始里有 id_kegg 列，就取出（以便只对 NA/空值做补全）；否则新建
            has_kegg_col <- !is.na(input$key_column_kegg) &&
              input$key_column_kegg %in% seq_along(colnames(metabo.data.raw.all.col[[ion_type]]))
            
            if (has_kegg_col) {
              old_kegg <- as.character(metabo.data.raw.all.col[[ion_type]][[input$key_column_kegg]])
              old_kegg[is.na(old_kegg) | trimws(old_kegg) == ""] <- NA_character_
              # 优先已存在的 KEGG，其次用 name 推断
              id_kegg_final <- ifelse(!is.na(old_kegg), old_kegg, derived_kegg_by_name)
            } else {
              id_kegg_final <- derived_kegg_by_name
            }
            
            tmp <- cbind(
              base_cols,
              name    = filled_name,
              id_kegg = id_kegg_final
            )
            
            metabo.data.raw[[ion_type]] <- data.frame(
              tmp,
              metabo.data.raw.all.col[[ion_type]][, input$key_column_sample_start:input$key_column_sample_end, drop = FALSE],
              check.names = FALSE
            )
            
            message("已自动生成 name，并据此补全 id_kegg。")
            
          }
          
          if (sum(is.na(metabo.data.raw[[ion_type]]$peak_name))>0) {
            
            metabo.data.raw[[ion_type]]$peak_name[which(is.na(metabo.data.raw[[ion_type]]$peak_name))] <- "NA"
            
          }
          
          metabo.data <- metabo.data.raw[[ion_type]][,c("peak_name", sample.data[[ion_type]]$`File name`)] %>% remove_rownames() %>% column_to_rownames("peak_name")
          
          metabo.data[metabo.data==input$miss.value.handle.type] <- NA
          
          group <- sample.data[[ion_type]]$Condition
          
          # 数据预处理 ----
          dir.create("2.数据预处理")
          
          ## 缺失值过滤 ----
          dir.create("2.数据预处理/1.缺失值过滤")
          
          if (input[['miss.value.handle.group']]=="inter.group"){
            
            ### 组内缺失值过滤 ----
            missing.index <- list()
            
            for (i in unique(group)) {
              
              missing.index[[i]] <- which(apply(metabo.data[ ,which(group==i)], 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)}))
              
            }
            
            missing.index <- Reduce(union,missing.index)
            
            if (length(missing.index)>0) {
              
              metabo.data2 <- metabo.data[-missing.index,]
              
            }else{
              
              metabo.data2 <- metabo.data
              
            }
            
            wb <- createWorkbook()
            addWorksheet(wb, "缺失值过滤数据")
            writeData(wb, "缺失值过滤数据", metabo.data2 %>% rownames_to_column(var = "peak_name"), rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/1.缺失值过滤/组内缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          if (input[['miss.value.handle.group']]=="global.group"){
            
            global.judge <- apply(metabo.data, 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)})
            
            if (sum(global.judge>0)) {
              
              metabo.data2 <- metabo.data[-which(global.judge),]
              
            }else{
              
              metabo.data2 <- metabo.data
              
            }
            
            wb <- createWorkbook()
            addWorksheet(wb, "缺失值过滤数据")
            writeData(wb, "缺失值过滤数据", metabo.data2 %>% rownames_to_column(var = "peak_name"), rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/1.缺失值过滤/全局缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          ## 缺失值做填充 ----
          dir.create("2.数据预处理/2.缺失值填充")
          
          ### 不填充 ----
          if (input[['miss.value.fill']]=="none"){
            
            metabo.data.fill <- metabo.data2
            
          }
          
          ### 填充 ----
          if (input[['miss.value.fill']]=="mean.group") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::mutate(group = group[match(sample, colnames(metabo.data2))]) %>%  # 添加分组信息
              dplyr::group_by(metabolite, group) %>%  # 按代谢物和组分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  # 关键修改：对全NA的分组返回0
                  coalesce(mean(value, na.rm = TRUE), 0), 
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="mean.global") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::group_by(metabolite) %>%  # 按代谢物分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(mean(value, na.rm = TRUE), 0),  # 如果全为 NA，填充为 0
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="median.group") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::mutate(group = group[match(sample, colnames(metabo.data2))]) %>%  # 添加分组信息
              dplyr::group_by(metabolite, group) %>%  # 按代谢物和组分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(median(value, na.rm = TRUE), 0), # 组内中位数填充
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="median.global") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::group_by(metabolite) %>%  # 按代谢物分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(median(value, na.rm = TRUE), 0),  # 如果全为 NA，填充为 0
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="min.global") {
            
            metabo.data.fill <- metabo.data2
            
            metabo.data.fill[is.na(metabo.data.fill)] <- min(as.matrix(metabo.data.fill),na.rm = T)
            
          }
          
          if (input[['miss.value.fill']]=="min2") {
            
            metabo.data.fill <- metabo.data2
            
            metabo.data.fill[is.na(metabo.data.fill)] <- mean(as.matrix(metabo.data.fill),na.rm = T)/2
            
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
          writeData(wb, "缺失值填充数据", metabo.data.fill %>% rownames_to_column(var = "peak_name"), rowNames = F)
          saveWorkbook(wb, paste0("2.数据预处理/2.缺失值填充/缺失值填充数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
          ## RSD 过滤 ----
          dir.create("2.数据预处理/3.RSD过滤")
          
          if (input$log.rsd.method=="none") {
            
            metabo.data.rsd <- metabo.data.fill %>% 
              dplyr::mutate(qc_rsd = matrixStats::rowSds(as.matrix(metabo.data.fill[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/matrixStats::rowMeans2(as.matrix(metabo.data.fill[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)) %>% 
              dplyr::filter(qc_rsd <= input$rsd.cutoff)
            
          }
          
          if (input$log.rsd.method=="log2") {
            
            metabo.data.rsd <- metabo.data.fill
            
            metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)] <- log2(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]+1)
            
            metabo.data.rsd <- metabo.data.rsd %>%
              dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)) %>% 
              dplyr::filter(qc_rsd <= input$rsd.cutoff)
            
          }
          
          if (input$log.rsd.method=="log10") {
            
            metabo.data.rsd <- metabo.data.fill
            
            metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)] <- log10(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]+1)
            
            metabo.data.rsd <- metabo.data.rsd %>% 
              dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)) %>% 
              dplyr::filter(qc_rsd <= input$rsd.cutoff)
            
          }
          
          wb <- createWorkbook()
          addWorksheet(wb, "RSD过滤")
          writeData(wb, "RSD过滤", metabo.data.rsd %>% rownames_to_column(var = "peak_name"), rowNames = F)
          saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
          metabo.data.rsd2 <- metabo.data.rsd%>%rownames_to_column("peak_name")%>%.[,c("peak_name","qc_rsd")]
          
          tryCatch({
            
            metabo.data.rsd2$peak_name <- as.character(metabo.data.rsd2$peak_name)
            metabo.data.rsd.whole <- left_join(metabo.data.rsd2, metabo.data.raw[[ion_type]],by="peak_name")
            
          },error=function(e){
            
            metabo.data.rsd2$peak_name <- as.double(metabo.data.rsd2$peak_name)
            metabo.data.rsd.whole <<- left_join(metabo.data.rsd2, metabo.data.raw[[ion_type]],by="peak_name")
            
          })
          
          write.xlsx(metabo.data.rsd.whole, paste0("2.数据预处理/过滤后的鉴定结果表_", ion_type, ".xlsx"))
          
          ## 归一化 ----
          metabo.data.fill <- metabo.data.fill.pre <- metabo.data.rsd[,-which(colnames(metabo.data.rsd)=="qc_rsd")]
          
          if (input[['normalized.handle.method']]=="none"){""}else{
            
            dir.create("2.数据预处理/4.归一化")
            
            ### sum归一化
            if (input[['normalized.handle.method']]=="sum") {
              
              sum.median <- apply(metabo.data.fill,2,sum,na.rm=T)
              
              for (i in 1:ncol(metabo.data.fill)) {
                
                metabo.data.fill[,i] <- input$sum_coef*metabo.data.fill[,i]/sum.median[i]
                
              }
              
            }
            
            ### QC样本归一化
            if (input[['normalized.handle.method']]=="qc") {
              
              data_csv <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)], metabo.data.fill)
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              
              sample_info <- sample.data[[ion_type]][,c("File name", "injection.order", "Condition")]%>%dplyr::rename(sample.name=`File name`, class=Condition)
              sample_info[["class"]][grep("^qc$", trimws(sample_info[["class"]]), ignore.case = T)] <- "QC"
              sample_info[["class"]][-grep("^qc$", trimws(sample_info[["class"]]), ignore.case = T)] <- "Subject"
              
              write_csv(data_csv, paste0("data.csv"))
              write_csv(sample_info, paste0("sample.info.csv"))
              
              metNor(
                ms1.data.name = "data.csv",
                sample.info.name = "sample.info.csv",
                minfrac.qc = 0,
                minfrac.sample = 0,
                optimization = TRUE,
                multiple = 5,
                threads = 3
              )
              
              metabo.data.rsd <- metabo.data.fill <- read_csv("svr_normalization_result/data_svr_normalization.csv")
              
              metabo.data.rsd[["QC.nor.rsd"]] <- metabo.data.rsd[["QC.nor.rsd"]]/100
              metabo.data.rsd[["sample.nor.rsd"]] <- metabo.data.rsd[["sample.nor.rsd"]]/100
              
              # 正确生成累积分布数据的方法
              cdf_data <- metabo.data.rsd %>% 
                # 确保使用正确的列名（假设实际列名是"qc_rsd"）
                dplyr::arrange(QC.nor.rsd) %>% 
                # 添加频次计数列（如果原数据没有）
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 5, by = 0.01)), # 创建1%间隔的区间
                  count = 1  # 每行代表一个峰
                ) %>% 
                dplyr::group_by(bin) %>% 
                dplyr::summarise(
                  qc_rsd = mean(QC.nor.rsd, na.rm = TRUE), # 取区间中值
                  count = sum(count) # 计算每个区间的峰数量
                ) %>% 
                dplyr::mutate(
                  cumulative = cumsum(count)/sum(count)*100 # 计算累积百分比
                ) %>% 
                dplyr::filter(!is.na(bin)) # 移除空区间
              
              metabo.data.rsd <- metabo.data.rsd %>% dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              wb <- createWorkbook()
              addWorksheet(wb, "RSD过滤")
              writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = F)
              saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              metabo.data.fill <- metabo.data.rsd
              colnames(metabo.data.fill) <- str_replace(colnames(metabo.data.fill), "^Sample","")
              metabo.data.fill <- metabo.data.fill[,c("name", sample.data[[ion_type]]$`File name`)]%>%column_to_rownames("name")
              
              system(paste0("mv svr_normalization_result 2.数据预处理/4.归一化/svr_normalization_result", ion_type))
              
            }
            
            ### QC-RLSC样本归一化
            if (input[['normalized.handle.method']]=="qc-rlsc") {
              
              # metabo.data.fill <- metabo.data.fill.pre
              ## 1) 仍然沿用你已写好的 data.csv / sample.info.csv 生成 ----------
              data_csv <- cbind(
                metabo.data.raw.all.col[[ion_type]][
                  match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]),
                  c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)
                ],
                metabo.data.fill
              )
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              data_csv$name <- as.character(data_csv$name)
              
              sample_info <- sample.data[[ion_type]][,c("File name", "injection.order", "Condition")] %>%
                dplyr::rename(sample.name=`File name`, class=Condition)
              sample_info[["class"]][grep("^qc$", trimws(sample_info[["class"]]), ignore.case = TRUE)] <- "QC"
              sample_info[["class"]][-grep("^qc$", trimws(sample_info[["class"]]), ignore.case = TRUE)] <- "Subject"
              
              readr::write_csv(data_csv,      "data.csv")
              readr::write_csv(sample_info,   "sample.info.csv")
              
              ## 2) 跑 QC-RLSC（qcrlscR::qc.rlsc.wrap），输出到 qcrlsc_normalization_result ----------
              if (!requireNamespace("qcrlscR", quietly = TRUE)) install.packages("qcrlscR")
              
              dir.create("qcrlsc_normalization_result", showWarnings = FALSE)
              
              # a) 读入
              dat_raw  <- readr::read_csv("data.csv",        show_col_types = FALSE)
              meta_raw <- readr::read_csv("sample.info.csv", show_col_types = FALSE)
              
              # b) 按注入顺序排序 & 构造 QC/批次标签
              meta <- meta_raw %>%
                dplyr::transmute(
                  sample = .data[["sample.name"]],
                  class  = ifelse(grepl("(?i)^qc$|\\bqc\\b|[._-]qc\\d*$", .data[["class"]]), "QC", "Sample"),
                  order  = as.numeric(.data[["injection.order"]]),
                  batch  = "1"
                ) %>%
                dplyr::arrange(order)
              
              stopifnot(any(meta$class=="QC"))
              
              # c) 把矩阵整理成 行=样本 × 列=特征
              stopifnot(all(meta$sample %in% colnames(dat_raw)))
              feat_id <- "name"  # 非样本列里，把 "name" 当特征ID
              mat <- as.matrix(dat_raw[, meta$sample, drop = FALSE])   # 列=样本
              rownames(mat) <- dat_raw[[feat_id]]
              X <- t(mat)  # 行=样本，列=特征（qcrlscR 需要）
              
              # d) 轻度缺失过滤（仅用于拟合，输出时会拼回）
              miss_rate <- colMeans(is.na(X))
              keep <- miss_rate < 0.5
              X_fit <- X[, keep, drop = FALSE]
              
              # e) 单批就关掉批内/批间步骤
              cls_qc <- factor(ifelse(meta$class=="QC","qc","sample"), levels=c("qc","sample"))
              cls_bl <- factor(meta$batch)
              intra_use <- FALSE
              shift_use <- FALSE
              
              # f) 跑 QC-RLSC（log10=TRUE + GCV 优化 span + 多项式 degree=2）
              Xc_fit <- qcrlscR::qc.rlsc.wrap(
                dat    = as.data.frame(X_fit),
                cls.qc = cls_qc,
                cls.bl = cls_bl,
                method = "divide",     # Dunn 2011 更推荐等比校正
                intra  = intra_use,    # 单批 FALSE
                opti   = TRUE,         # GCV 自动寻优 span
                log10  = TRUE,
                outl   = TRUE,
                shift  = shift_use,    # 单批 FALSE
                degree = 2
              )
              
              # g) 拼回所有特征
              Xc <- matrix(NA_real_, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))
              Xc[, keep] <- as.matrix(Xc_fit)
              
              # h) 写出 “样本×特征” 和 “特征×样本”
              out_s_by_f <- Xc %>% as.data.frame() %>% tibble::rownames_to_column("name")
              readr::write_csv(out_s_by_f, file.path("qcrlsc_normalization_result","qcrlsc_sample_by_feature.csv"))
              
              out_f_by_s <- t(Xc) %>% as.data.frame() %>% tibble::rownames_to_column("name")
              
              # 统一 key 类型：都用字符
              dat_raw  <- readr::read_csv("data.csv", show_col_types = FALSE) %>%
                dplyr::mutate(name = as.character(name))
              
              out_f_by_s <- out_f_by_s %>%
                dplyr::mutate(name = as.character(name))
              
              data_qcrlsc_norm <- dplyr::left_join(
                dat_raw[, c("name","mz","rt")],
                out_f_by_s, by = "name"
              )
              
              # 为了和你 SVR 的命名风格一致，做一个总表：特征基本信息 + 归一化后的样本强度
              data_qcrlsc_norm <- dplyr::left_join(
                dat_raw[, c("name","mz","rt")],   # 特征信息
                out_f_by_s, by="name"
              )
              readr::write_csv(data_qcrlsc_norm, file.path("qcrlsc_normalization_result","data_qcrlsc_normalization.csv"))
              
              ## 3) 计算 RSD（QC 与 Subject），并生成和你 SVR 分支一致的后续产物 ----------
              # 样本分组
              qc_samples     <- meta$sample[meta$class=="QC"]
              subject_samples<- meta$sample[meta$class!="QC"]
              
              # 取 “特征×样本” 矩阵（不含 name/mz/rt）
              M <- as.matrix(data_qcrlsc_norm[, subject_samples, drop=FALSE])
              M_qc <- as.matrix(data_qcrlsc_norm[, qc_samples, drop=FALSE])
              
              # RSD 函数（百分数）
              rsd_perc <- function(v) {
                v <- as.numeric(v)
                mu <- mean(v, na.rm=TRUE)
                if (!is.finite(mu) || mu<=0) return(NA_real_)
                100 * stats::sd(v, na.rm=TRUE) / mu
              }
              
              QC_rsd_perc     <- apply(M_qc,     1, rsd_perc)
              Sample_rsd_perc <- apply(M,        1, rsd_perc)
              
              metabo.data.rsd <- data_qcrlsc_norm %>%
                dplyr::mutate(
                  `QC.nor.rsd.perc`     = QC_rsd_perc,
                  `sample.nor.rsd.perc` = Sample_rsd_perc,
                  `QC.nor.rsd`          = QC_rsd_perc/100,       # ← 分数形式（0~1）
                  `sample.nor.rsd`      = Sample_rsd_perc/100
                )
              
              # CDF 数据（注意：这里用“分数制”的 QC.nor.rsd，1% = 0.01）
              cdf_data <- metabo.data.rsd %>%
                dplyr::arrange(QC.nor.rsd) %>%
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 0.50, by = 0.01), include.lowest = TRUE)
                ) %>%
                dplyr::filter(!is.na(bin)) %>%
                dplyr::group_by(bin) %>%
                dplyr::summarise(
                  qc_rsd_mid = mean(QC.nor.rsd, na.rm = TRUE),
                  count      = dplyr::n(),
                  .groups    = "drop"
                ) %>%
                dplyr::mutate(cumulative = cumsum(count)/sum(count)*100)
              
              # RSD 过滤（沿用你原来的阈值语义：input$rsd.cutoff 是“分数制”，比如 0.3 代表 30%）
              metabo.data.rsd <- metabo.data.rsd %>%
                dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              # 写 Excel
              wb <- openxlsx::createWorkbook()
              openxlsx::addWorksheet(wb, "RSD过滤")
              openxlsx::writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = FALSE)
              dir.create("2.数据预处理/3.RSD过滤", recursive = TRUE, showWarnings = FALSE)
              openxlsx::saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              # 回填 metabo.data.fill，列只保留样本列
              metabo.data.fill <- metabo.data.rsd
              colnames(metabo.data.fill) <- stringr::str_replace(colnames(metabo.data.fill), "^Sample","")
              metabo.data.fill <- metabo.data.fill[, c("name", sample.data[[ion_type]]$`File name`)] %>% tibble::column_to_rownames("name")
              
              # 挪动输出目录
              system(paste0("mv qcrlsc_normalization_result 2.数据预处理/4.归一化/", "qcrlsc_normalization_result_", ion_type))
              
            }
            
            ### NormAE样本归一化
            if (input[['normalized.handle.method']]=="normae") {
              
              # metabo.data.fill <- metabo.data.fill.pre
              
              data_csv <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)], metabo.data.fill)
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              
              sample_info <- sample.data[[ion_type]][,c("File name", "Condition", "injection.order", "batch")]%>%dplyr::rename(sample.name=`File name`, class=Condition)
              ## 假设你已经生成了 data_csv 和 sample_info（与之前相同）
              ## 1) 基础匹配检查：样本名是否一致
              sample_cols <- colnames(data_csv)[-(1:3)]
              stopifnot(!anyDuplicated(sample_cols))
              stopifnot(!anyDuplicated(sample_info$sample.name))
              
              # 映射并提醒不匹配
              not_in_data   <- setdiff(sample_info$sample.name, sample_cols)
              not_in_sample <- setdiff(sample_cols, sample_info$sample.name)
              if (length(not_in_data) || length(not_in_sample)) {
                cat("【警告】样本名不匹配：\n- 在 sample_info 不在 data 的：", paste(head(not_in_data,20), collapse=", "), "\n",
                    "- 在 data 不在 sample_info 的：", paste(head(not_in_sample,20), collapse=", "), "\n", sep = "")
                stop("请先修正样本名匹配再运行 NormAE。")
              }
              
              ## 2) 类型与取值检查：class/batch/injection.order
              sample_info$class <- ifelse(grepl("^qc$", trimws(sample_info$class), ignore.case = TRUE), "QC", "Subject")
              sample_info$batch <- suppressWarnings(as.integer(sample_info$batch))
              sample_info$batch[is.na(sample_info$batch)] <- 1L
              sample_info$injection.order <- suppressWarnings(as.integer(sample_info$injection.order))
              if (any(is.na(sample_info$injection.order))) {
                stop("injection.order 中存在非数字条目，请先清洗为整数。")
              }
              
              # cat("Class 计数：\n"); print(table(sample_info$class))
              # cat("各 batch×class 计数：\n"); print(table(sample_info$batch, sample_info$class))
              
              ## 3) 数值矩阵清洗：去掉非有限值、零方差行
              X <- as.matrix(data_csv[, -(1:3)])
              mode(X) <- "numeric"
              
              # 标记非有限值
              non_finite_rate <- mean(!is.finite(X))
              cat(sprintf("非有限值比例：%.4f\n", non_finite_rate))
              
              # 将非有限值先设为 NA，计算零方差行
              X[!is.finite(X)] <- NA_real_
              
              # 至少保留 >=3 个有效观测，且范围>0（避免零方差）
              keep <- rowSums(is.finite(X)) >= 3 & apply(X, 1, function(v) {
                rng <- range(v, na.rm = TRUE)
                is.finite(rng[1]) && is.finite(rng[2]) && (diff(rng) > 0)
              })
              
              cat(sprintf("原始特征数：%d，过滤后保留：%d（去掉全NA/零方差）\n", nrow(X), sum(keep)))
              
              # 可选：对剩余少量 NA 做填补（NormAE/Sklearn 一般不接受 NA）
              X <- X[keep, , drop = FALSE]
              X[is.na(X)] <- 0
              
              # 写回 data.csv（前三列为 name/mz/rt）
              data_csv_clean <- cbind(data_csv[keep, 1:3, drop = FALSE], as.data.frame(X, check.names = FALSE))
              readr::write_csv(data_csv_clean, "data.csv")  # 覆盖为清洗后的
              readr::write_csv(sample_info, "sample.info.csv")
              
              
              dir.create("normae_normalization_result", showWarnings = FALSE)
              
              # —— 用 conda 的 normae 环境运行 —— #
              run_normae <- function(meta = "data.csv", sample = "sample.info.csv", out = "./normae_normalization_result",
                                     env = "normae", log_file = "normae_run.log") {
                meta <- normalizePath(meta); sample <- normalizePath(sample); out <- normalizePath(out)
                log_path <- file.path(out, log_file)
                
                # 找 conda
                conda <- Sys.which("conda")
                if (conda == "") {
                  cand <- path.expand(c("~/anaconda3/bin/conda", "~/miniconda3/bin/conda", "/opt/anaconda3/bin/conda"))
                  conda <- cand[file.exists(cand)][1]
                  if (is.na(conda)) stop("找不到 conda，可临时把 anaconda3/miniconda3/bin 加入 PATH。")
                }
                
                run_cmd <- function(args, tag) {
                  out_lines <- system2(conda, args, stdout = TRUE, stderr = TRUE)
                  cat(out_lines, sep = "\n", file = log_path, append = TRUE)
                  status <- attr(out_lines, "status")
                  list(status = if (is.null(status)) 0L else status, out = out_lines, tag = tag)
                }
                
                base <- c("run", "--no-capture-output", "-n", env)
                
                # 方案1：python -m normae（最稳）
                r1 <- run_cmd(c(base, "python", "-m", "normae",
                                "--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                              "python -m normae")
                if (r1$status == 0) { message("NormAE 完成（python -m normae）。日志：", log_path); return(invisible(r1$out)) }
                
                # 方案2：CLI 脚本 normae
                r2 <- run_cmd(c(base, "normae",
                                "--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                              "normae CLI")
                if (r2$status == 0) { message("NormAE 完成（normae CLI）。日志：", log_path); return(invisible(r2$out)) }
                
                # 方案3：绝对路径兜底（按你的安装路径修改）
                normae_bin <- path.expand(file.path("~", "anaconda3", "envs", env, "bin", "normae"))
                if (file.exists(normae_bin)) {
                  r3 <- system2(normae_bin,
                                c("--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                                stdout = TRUE, stderr = TRUE)
                  cat(r3, sep = "\n", file = log_path, append = TRUE)
                  st3 <- attr(r3, "status"); if (is.null(st3)) st3 <- 0L
                  if (st3 == 0) { message("NormAE 完成（绝对路径）。日志：", log_path); return(invisible(r3)) }
                }
                
                # 汇总报错片段
                tail_msg <- paste(tail(c(r1$out, r2$out), 40), collapse = "\n")
                stop("NormAE 运行失败。请查看日志：", log_path, "\n关键信息：\n", tail_msg)
              }
              
              run_normae()
              
              metabo.data.rsd <- read_csv("normae_normalization_result/X_clean.csv")
              
              qc_cols <- sample.data[[ion_type]]$`File name`[grep("QC",sample.data[[ion_type]]$Condition,ignore.case = T)]
              rsd_perc <- function(v){ mu <- mean(v, na.rm=TRUE); if(!is.finite(mu)||mu<=0) return(NA_real_); 100*stats::sd(v,na.rm=TRUE)/mu }
              metabo.data.rsd$QC.nor.rsd <- if (length(qc_cols) >= 2) apply(metabo.data.rsd[, qc_cols, drop=FALSE], 1, rsd_perc) else rep(NA_real_, nrow(metabo.data.rsd))
              
              # CDF 数据（注意：这里用“分数制”的 QC.nor.rsd，1% = 0.01）
              cdf_data <- metabo.data.rsd %>%
                dplyr::arrange(QC.nor.rsd) %>%
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 0.50, by = 0.01), include.lowest = TRUE)
                ) %>%
                dplyr::filter(!is.na(bin)) %>%
                dplyr::group_by(bin) %>%
                dplyr::summarise(
                  qc_rsd_mid = mean(QC.nor.rsd, na.rm = TRUE),
                  count      = dplyr::n(),
                  .groups    = "drop"
                ) %>%
                dplyr::mutate(cumulative = cumsum(count)/sum(count)*100)
              
              # RSD 过滤（沿用你原来的阈值语义：input$rsd.cutoff 是“分数制”，比如 0.3 代表 30%）
              metabo.data.rsd <- metabo.data.rsd %>%
                dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              # 写 Excel
              wb <- openxlsx::createWorkbook()
              openxlsx::addWorksheet(wb, "RSD过滤")
              openxlsx::writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = FALSE)
              dir.create("2.数据预处理/3.RSD过滤", recursive = TRUE, showWarnings = FALSE)
              openxlsx::saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              # 回填 metabo.data.fill，列只保留样本列
              metabo.data.fill <- metabo.data.rsd
              metabo.data.fill <- metabo.data.fill[, c("name", sample.data[[ion_type]]$`File name`)] %>% tibble::column_to_rownames("name")
              
              system(paste0("mv normae_normalization_result 2.数据预处理/4.归一化/normae_normalization_result_", ion_type))
              
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
            
            ### 75分位数归一化
            if (input[['normalized.handle.method']]=="percent_0.75"){
              
              # 定义处理函数
              normalize_by_quantile <- function(column) {
                # 计算75%分位数
                q75 <- quantile(column, 0.75, na.rm = TRUE)
                # 处理特殊情况：如果分位数为0，则直接返回原始值
                if(q75 == 0) {
                  return(column * 1000)
                } else {
                  # 归一化处理：每个值除以75%分位数，再乘以1000
                  return(column / q75 * 1000)
                }
              }
              
              # 对每一列应用处理函数
              # 使用sapply会返回矩阵，这里保持原始数据框结构
              metabo.data.fill <- as.data.frame(lapply(metabo.data.fill, normalize_by_quantile), row.names = rownames(metabo.data.fill))
              
            }
            
            wb <- createWorkbook()
            addWorksheet(wb, "归一化处理")
            writeData(wb, "归一化处理", metabo.data.fill %>% rownames_to_column(var = "peak_name"), rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/4.归一化/归一化处理数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          ## 定量值取log ----
          if (input[["log.handle.method"]]=="log2") {
            
            metabo.data.fill <- log2(metabo.data.fill+1)
            
          }
          
          if (input[["log.handle.method"]]=="log10") {
            
            metabo.data.fill <- log10(metabo.data.fill+1)
            
          }
          
          ## 批次矫正 ----
          if (input[['batch.correct']]=="none"){""}else{
            
            dir.create("2.数据预处理/5.批次校正")
            
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
            writeData(wb, "批次矫正", metabo.data.fill %>% rownames_to_column(var = "peak_name"), rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/5.批次矫正/批次矫正数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          assign(paste0("metabo.data.", ion_type), metabo.data.fill)
          
          # QC 质控 ----
          dir.create("3.QC")
          
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
          dir.create("3.QC/1.TIC")
          
          QC.RT <- data.frame(EG.ApexRT=metabo.data.raw[[ion_type]]$rt[match(rownames(metabo.data.fill), metabo.data.raw[[ion_type]]$peak_name)]%>%as.numeric(), metabo.data.fill[,grepl("^QC$", trimws(sample.data[[ion_type]]$Condition), ignore.case = T)]) %>% reshape2::melt(.,id=c("EG.ApexRT"),variable.name = "R.FileName",value.name = "FG.MS1Quantity")
          
          QC.RT2 <- dplyr::mutate(QC.RT,EG.ApexRT.bin=round(EG.ApexRT,digits = 0))
          
          QC.RT3 <- QC.RT2 %>% dplyr::group_by(R.FileName,EG.ApexRT.bin)%>%dplyr::summarize(Quantity.sum=sum(FG.MS1Quantity,na.rm=T))
          
          QC.RT3$type1<-"MS1"
          
          QC.TIC <- QC.RT3
          
          g_TIC_MS1 <- {
            
            highchart() %>% 
              hc_add_yAxis(lineWidth = 3,title = list(text = "Intensity"))%>%
              hc_xAxis(title = list(text = "Retention Time")) %>%
              hc_add_series(data = QC.TIC[which(QC.TIC$type1=="MS1"),],type="spline",hcaes(x = EG.ApexRT.bin, y = Quantity.sum, group = R.FileName)) %>%
              hc_tooltip(split=T,valueDecimals=3)%>%
              hc_plotOptions(series = list(marker = list(symbol = "circle")))%>%
              hc_chart(zoomType="x",events=list(selection=select_zoom))%>%
              hc_exporting(enabled = TRUE,buttons = list( contextButton = list(menuItems = list('downloadPNG', 'downloadSVG',"downloadPDF","downloadJPEG","printChart","viewFullscreen"))))
            
          }
          
          htmlwidgets::saveWidget(widget = g_TIC_MS1, file = paste0("3.QC/1.TIC/TIC_", ion_type, ".html"))
          
          webshot::webshot(url = paste0("3.QC/1.TIC/TIC_", ion_type, ".html"), 
                           file = c(paste0("3.QC/1.TIC/TIC_", ion_type, ".png"),
                                    paste0("3.QC/1.TIC/TIC_", ion_type, ".pdf")),
                           delay = 3)
          
          system(paste0("rm -rf '3.QC/1.TIC/TIC_", ion_type, "_files'"))
          
          write_csv(QC.TIC, paste0("3.QC/1.TIC/TIC_", ion_type, ".csv"))
          
          ## QC样本的相关性图 ----
          dir.create("3.QC/2.correlation")
          
          result<-corr.test(metabo.data.fill[,grepl("^QC$", trimws(sample.data[[ion_type]]$Condition), ignore.case = T)], method = "pearson",adjust="none",alpha=.05)
          rmatrix<-result$r
          pmatrix<-result$p
          
          col <- colorRampPalette(c("darkblue", "white", "red"))(200)
          
          pdf(paste0("3.QC/2.correlation/correlation_", ion_type, ".pdf"))
          p1 <- corrplot::corrplot(rmatrix,method = "circle",col = col,tl.col="black",type="upper")
          dev.off()
          png(paste0("3.QC/2.correlation/correlation_", ion_type, ".png"))
          p1 <- corrplot::corrplot(rmatrix,method ="circle",col = col,tl.col="black",type="upper")
          dev.off()
          write.csv(rmatrix, paste0("3.QC/2.correlation/correlation_", ion_type, ".csv"))
          write.csv(pmatrix, paste0("3.QC/2.correlation/correlation_pvalue_", ion_type, ".csv"))
          
          ## QC样本的RSD分布图 ----
          dir.create("3.QC/3.RSD")
          
          dat <- data_to_boxplot(metabo.data.rsd, qc_rsd, name="QC.RSD")
          
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
          
          htmlwidgets::saveWidget(widget = g_rsd, file = paste0("3.QC/3.RSD/RSD_", ion_type, ".html"))
          
          webshot::webshot(url = paste0("3.QC/3.RSD/RSD_", ion_type, ".html"), 
                           file = c(paste0("3.QC/3.RSD/RSD_", ion_type, ".png"),
                                    paste0("3.QC/3.RSD/RSD_", ion_type, ".pdf")),
                           delay = 3)
          
          system(paste0("rm -rf '3.QC/3.RSD/RSD_", ion_type, "_files'"))
          
          write.csv(metabo.data.rsd$qc_rsd, paste0("3.QC/3.RSD/RSD_", ion_type, ".csv"))
          
          ## pca ----
          dir.create("3.QC/4.PCA")
          
          pca_qc <- metabo.data.fill
          
          pca_qc[is.na(pca_qc)]<-0
          
          pca_qc[pca_qc=="NaN"]<-0
          
          pca_qc[,1:ncol(pca_qc)]<-lapply(pca_qc[,1:ncol(pca_qc)],as.numeric)
          
          pc.cr <- prcomp(t(pca_qc), scale. = T)
          
          pca_group <- sample.data[[ion_type]]$Condition
          
          # 自动生成颜色映射规则（无需预知分组名）
          generate_color_mapping <- function(groups) {
            # 按字母顺序排序保证可重复性
            sorted_groups <- sort(unique(as.character(groups)))
            
            # 核心配色（参考图片中的黄金色系）
            base_colors <- c("#f8766d", "#1F77B4")  # 黄金 + 经典蓝
            
            # 动态扩展配色方案
            if(length(sorted_groups) > 2) {
              extended_colors <- viridis::viridis(
                n = length(sorted_groups) - 2,
                begin = 0.2,  # 避免过亮黄色
                end = 0.8,    # 避免过深紫色
                option = "A"  # viridis方案
              )
              color_palette <- c(base_colors, extended_colors)
            } else {
              color_palette <- base_colors
            }
            
            # 创建命名向量
            setNames(color_palette, sorted_groups)
          }
          
          dynamic_colors <- generate_color_mapping(pca_group)
          
          g <- ggord::ggord(pc.cr, pca_group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, alpha=0.5, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
            geom_text(aes(label=lab),vjust=-0.5,hjust="inward",check_overlap = T,size=3,color="black",show.legend=T) +
            theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8), panel.grid = element_blank()) +
            labs(title="PCA Plot",fill=NULL,colour=NULL,shape=NULL)+
            scale_color_brewer(palette="Set2") +   # 修改点和文本颜色
            scale_fill_brewer(palette="Set2")      # 修改椭圆填充色
          
          g_pca <- g %>% ggplotly()
          
          g_pca$x$data[[3]]$showlegend<-F
          g_pca$x$data[[4]]$showlegend<-F
          g_pca$x$data[[5]]$showlegend<-T
          g_pca$x$data[[5]]$name<-"label"
          
          htmlwidgets::saveWidget(widget = g_pca, file = paste0("3.QC/4.PCA/PCA_", ion_type, ".html"))
          
          system(paste0("rm -rf '3.QC/4.PCA/PCA_", ion_type, "_files'"))
          
          g <- ggord::ggord(pc.cr, pca_group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, alpha=0.5, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
            geom_text_repel(aes(label=lab))+
            theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8), panel.grid = element_blank()) +
            labs(title="PCA Plot", fill=NULL, colour=NULL, shape=NULL)+
            scale_color_brewer(palette="Set2") +   # 修改点和文本颜色
            scale_fill_brewer(palette="Set2")
          
          ggsave(paste0("3.QC/4.PCA/PCA_", ion_type, ".png"), g)
          ggsave(paste0("3.QC/4.PCA/PCA_", ion_type, ".pdf"), g)
          
          ## CV分布图 ----
          dir.create("3.QC/5.CV", recursive = TRUE, showWarnings = FALSE)
          
          metabo.data.cv <- metabo.data.fill
          cv_group <- data.frame(Sample=sample.data[[ion_type]]$`File name`, Group=sample.data[[ion_type]]$Condition)
          
          # 假设你的数据框名称为df，列为样本名，行为代谢物
          # 1. 数据预处理（长格式转换和分组）
          df_long <- metabo.data.cv %>%
            rownames_to_column("Metabolite") %>%
            pivot_longer(-Metabolite, names_to = "Sample", values_to = "Intensity") %>%
            left_join(.,cv_group)
          
          # 2. 计算分组CV（过滤零均值情况）
          cv_data <- df_long %>%
            dplyr::group_by(Metabolite, Group) %>%
            dplyr::summarise(
              CV = sd(Intensity) / mean(Intensity),
              .groups = "drop"
            ) %>%
            dplyr::filter(is.finite(CV))  # 移除无效计算结果
          
          # 3. 计算累积比例（阶梯型分布）
          cumulative_plot <- cv_data %>%
            dplyr::group_by(Group) %>%
            dplyr::arrange(CV) %>% 
            dplyr::mutate(Cumulative = seq_along(CV)/n()) %>%
            dplyr::ungroup()
          
          # 针对QC组添加延伸至CV=1.0的数据点
          if ("QC" %in% cumulative_plot$Group) {
            # 获取QC组的最大CV值
            qc_max <- cumulative_plot %>% 
              dplyr::filter(Group == "QC") %>% 
              dplyr::summarise(max_cv = max(CV)) %>% 
              dplyr::pull(max_cv)
            
            # 如果最大CV小于1.0则添加延伸点
            if (qc_max < 1) {
              cumulative_plot <- cumulative_plot %>%
                bind_rows(
                  tibble(
                    Metabolite = NA_character_,
                    Group = "QC",
                    CV = 1.0,
                    Cumulative = 1.0
                  )
                ) %>%
                dplyr::arrange(Group, CV)
            }
          }
          
          ref_lines <- data.frame(
            type = c("vline", "vline", "hline", "hline"),
            value = c(0.3, 0.5, 0.75, 0.85),
            label = c("CV=0.3", "CV=0.5", "75%", "85%"),
            x = c(0.3, 0.5, 0.95, 0.95), # 文本x坐标
            y = c(0.95, 0.95, 0.75, 0.85) # 文本y坐标
          )
          
          # 4. 绘制图形
          p1 <- ggplot(cumulative_plot, aes(x = CV, y = Cumulative, color = Group)) +
            geom_step(direction = "hv", linewidth = 0.8) +  # 阶梯状曲线
            geom_vline(xintercept = c(0.3, 0.5), linetype = "dashed", color = "gray40") +
            geom_hline(yintercept = c(0.75, 0.85), linetype = "dashed", color = "gray40") +
            scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0, 0.3, 0.5, 1)) +  # 固定x轴范围
            scale_y_continuous(labels = scales::percent_format(), breaks = c(0, 0.75, 0.85, 1)) +    # y轴百分比
            labs(x = "Coefficient of Variation (CV)", 
                 y = "Cumulative Proportion",
                 color = "Sample Group") +
            theme_minimal() +
            theme(legend.position = "top",
                  panel.grid.minor = element_blank())
          
          ggsave(paste0("3.QC/5.CV/cv_", ion_type, ".pdf"), p1)
          ggsave(paste0("3.QC/5.CV/cv_", ion_type, ".png"), p1)
          
          # 判断是否除QC外只有一组数据 ----
          if (length(unique(sample.data[[ion_type]][["Condition"]][!grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T)]))>1) {
            
            # 多元统计分析 ----
            dir.create("4.多元统计分析")
            
            ## pca ----
            dir.create(paste0("4.多元统计分析/1.PCA/all/", ion_type),recursive = T)
            
            NcrossvalI <- min(7,ncol(pca_qc))
            
            model_pca <- tryCatch(opls(t(pca_qc), scaleC = 'standard', fig.pdfC = 'none', info.txtC = 'none', crossvalI = NcrossvalI),
                                  error = function(cnd) {opls(t(pca_qc), predI=3, scaleC = 'standard', fig.pdfC = 'none', info.txtC = 'none', crossvalI = NcrossvalI)
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
              inner_join(sample.data[[ion_type]] %>% dplyr::select(`File name`, Condition), by=c("sample"="File name"))%>%
              .[!grepl("^qc$", trimws(.[["Condition"]]),ignore.case = T),]
            
            if (length(unique(df_pca$Condition))>1) {
              
              g <- ggpubr::ggscatter(df_pca,
                                     x='p1',
                                     y='p2',
                                     fill = 'Condition',
                                     color = 'Condition',
                                     shape = 'Condition',
                                     # label = df_pca$sample,
                                     show.legend.text = F,
                                     ellipse = T)+
                scale_color_brewer(palette="Set2") +   # 修改点和文本颜色
                scale_fill_brewer(palette="Set2")+
                geom_text_repel(aes(label=sample))+
                xlab(glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'))+
                ylab(glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)'))+
                ggtitle('PCA')+
                theme_bw()+
                # coord_fixed(ratio=1)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
              
              ggsave(g,filename = paste0("4.多元统计分析/1.PCA/all/", ion_type, "/PCA_all_", ion_type, ".pdf"),width = 10,height = 6)
              ggsave(g,filename = paste0("4.多元统计分析/1.PCA/all/", ion_type, "/PCA_all_", ion_type, ".png"),width = 10,height = 6)
              
            }
            
            sum_pca[[paste0("all sample ", ion_type)]] <- data.frame(title=paste0("all sample ", ion_type), type="PCA", A=nrow(df_pca), N=nrow(model_pca@modelDF), `R2X(cum)`=model_pca@summaryDF$`R2X(cum)`)
            
            for (g in 1:length(group.data)) {
              
              dir.create(paste0("4.多元统计分析/1.PCA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              
              experimental <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[1]
              control <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[2]
              
              pca_qc <- metabo.data.fill[, which(sample.data[[ion_type]]$Condition %in% c(experimental, control))]
              
              pca_qc[is.na(pca_qc)]<-0
              
              pca_qc[pca_qc=="NaN"]<-0
              
              pca_qc[,1:ncol(pca_qc)]<-lapply(pca_qc[,1:ncol(pca_qc)],as.numeric)
              
              pc.cr<-prcomp(t(pca_qc), scale. = T)
              
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
              
              sum_pca[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                            type="PCA", 
                                                                                            A=nrow(df_pca), 
                                                                                            N=nrow(model_pca@modelDF), 
                                                                                            `R2X(cum)`=model_pca@summaryDF$`R2X(cum)`)
              
              g1 <- ggpubr::ggscatter(df_pca,
                                      x='p1',
                                      y='p2',
                                      fill = 'Condition',
                                      color = 'Condition',
                                      shape = 'Condition',
                                      # label = df_pca$sample,
                                      show.legend.text = F,
                                      ellipse = T)+
                scale_color_brewer(palette="Set2") +   # 修改点和文本颜色
                scale_fill_brewer(palette="Set2")+
                geom_text_repel(aes(label=sample))+
                xlab(glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'))+
                ylab(glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)'))+
                ggtitle('PCA')+
                theme_bw()+
                # coord_fixed(ratio=1)+
                theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())
              
              ggsave(g1,filename = paste0('4.多元统计分析/1.PCA/',group.data[g],'/',ion_type,'/PCA_', experimental, "_vs_", control, "_", ion_type, ".pdf"),width = 10,height = 6)
              ggsave(g1,filename = paste0('4.多元统计分析/1.PCA/',group.data[g],'/',ion_type,'/PCA_', experimental, "_vs_", control, "_", ion_type, ".png"),width = 10,height = 6)
              
            }
            
            ## pls-da/opls-da ----
            dir.create(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type), recursive = T, showWarnings = F)
            dir.create("4.多元统计分析/3.OPLS-DA")
            
            ### 所有样本 ----
            # 获取plsda数据
            plsda_data <- metabo.data.fill[, !grepl("^qc$", trimws(group), ignore.case = T)]
            
            # 获取分组信息
            plsda_group <- group[!grepl("^qc$", trimws(group), ignore.case = T)]
            
            NcrossvalI <- min(7,ncol(plsda_data))
            
            # 运行PLS-DA模型
            plsda_model <- opls(t(plsda_data), plsda_group, crossvalI=NcrossvalI, permI = 200)
            
            if (nrow(plsda_model@modelDF)==0 | nrow(plsda_model@modelDF)==1) {
              
              plsda_model <- opls(t(plsda_data), plsda_group, predI = 2, crossvalI=NcrossvalI, permI = 200)
              
            }
            
            if (nrow(plsda_model@modelDF)>=2){
              
              # 从模型对象中提取得分矩阵
              pls_scores <- data.frame(
                p1 = plsda_model@scoreMN[, 1],
                p2 = plsda_model@scoreMN[, 2],
                Groups = plsda_group,
                sample = rownames(t(plsda_data))
              )
              
              x_lab <- paste0("P1(", round(plsda_model@modelDF$R2X[1]*100, 1), "*'%')")
              y_lab <- paste0("P2(", round(plsda_model@modelDF$R2X[2]*100, 1), "*'%')")
              
              g_plsda <- ggplot(pls_scores, aes(x = p1, y = p2, color = Groups)) +
                geom_point(size = 4) +
                stat_ellipse(level = 0.95, alpha = 0.8) +
                geom_text_repel(aes(label = sample)) +
                scale_color_manual(values = dynamic_colors) + 
                labs(
                  x = parse(text = x_lab),  # 解析数学表达式
                  y = parse(text = y_lab),
                  title = "PLS-DA Score Plot"
                ) +
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))
              
              ggsave(paste0("4.多元统计分析/2.PLS-DA/all/",ion_type,"/PLS-DA_all_", ion_type, ".pdf"), g_plsda, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/2.PLS-DA/all/",ion_type,"/PLS-DA_all_", ion_type, ".png"), g_plsda, width = 10, height = 6)
              
            }
            
            sum_plsda[[paste0("all sample ", ion_type)]] <- data.frame(title=paste0("all sample ", ion_type), 
                                                                       type="PLS-DA", 
                                                                       A=nrow(pls_scores), 
                                                                       N=nrow(plsda_model@modelDF), 
                                                                       `R2X(cum)`=plsda_model@summaryDF$`R2X(cum)`, 
                                                                       `R2Y(cum)`=plsda_model@summaryDF$`R2Y(cum)`, 
                                                                       `Q2(cum)`=plsda_model@summaryDF$`Q2(cum)`)
            
            # 从模型对象提取置换检验结果
            perm_data <- plsda_model@suppLs$permMN %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column("perm_id") %>% 
              tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                  names_to = "metric", 
                                  values_to = "value")
            
            # 计算原始模型指标
            orig_r2y <- plsda_model@summaryDF$`R2Y(cum)`
            orig_q2 <- plsda_model@summaryDF$`Q2(cum)`
            
            # 生成置换检验图
            perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
              # 绘制置换散点
              geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
              scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
              # 添加回归线
              geom_smooth(
                aes(group = metric, color = metric),
                method = "lm", 
                formula = y ~ x,
                se = FALSE, 
                fullrange = TRUE,  # 启用全范围延伸
                # linewidth = 1.2,   # 加粗线宽匹配图片
                linetype = "dashed"
              ) +
              # 标注原始模型值（精确坐标定位）
              annotate("point", x = orig_r2y, y = orig_r2y, 
                       color = "#00BFC4", size = 3, shape = 18) +  # R2Y菱形标记
              annotate("point", x = orig_q2, y = 0, 
                       color = "#7CAE00", size = 3, shape = 18) +  # Q2菱形标记
              # 坐标轴设置（匹配图片范围）
              scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
              scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
              # 标签系统（完全复现图片文字）
              labs(title = "permutation_test",
                   x = "200 permutations 1 components",
                   y = "value",
                   color = "Metric") +
              # 主题优化（精确字体和网格匹配）
              theme_bw() +
              theme(panel.grid.minor = element_blank(),
                    axis.text = element_text(size = 10, color = "black"),
                    legend.position = "right",plot.title = element_text(hjust = 0.5))
            
            # 输出高清图
            ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
            ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
            
            ### 分组样本（PLSDA和OPLS-DA） ----
            for (g in 1:length(group.data)) {
              
              dir.create(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              dir.create(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              
              experimental <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[1]
              control <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[2]
              
              plsda_data <- metabo.data.fill[, which(group %in% c(experimental, control))]
              
              # 获取分组信息
              plsda_group <- group[which(group %in% c(experimental, control))]
              
              NcrossvalI <- min(7,ncol(plsda_data))
              
              #### 运行PLS-DA模型 ----
              plsda_model <- opls(t(plsda_data), plsda_group, crossvalI=NcrossvalI, permI = 200)
              
              if (nrow(plsda_model@modelDF)==0 | nrow(plsda_model@modelDF)==1) {
                
                plsda_model <- opls(t(plsda_data), plsda_group, predI = 2, crossvalI=NcrossvalI, permI = 200)
                
              }
              
              # if (nrow(plsda_model@modelDF)==0) {
              #   
              #   plsda_model <- opls(t(plsda_data), plsda_group, predI = 1, crossvalI=NcrossvalI, permI = 200)
              #   
              # }
              
              if (nrow(plsda_model@modelDF)>=2) {
                
                # 从模型对象中提取得分矩阵
                pls_scores <- data.frame(
                  p1 = plsda_model@scoreMN[, 1],
                  p2 = plsda_model@scoreMN[, 2],
                  Groups = plsda_group,
                  sample = rownames(t(plsda_data))
                )
                
                x_lab <- paste0("P1(", round(plsda_model@modelDF$R2X[1]*100, 1), "*'%')")
                y_lab <- paste0("P2(", round(plsda_model@modelDF$R2X[2]*100, 1), "*'%')")
                
                g_plsda <- ggplot(pls_scores, aes(x = p1, y = p2, color = Groups)) +
                  geom_point(size = 4) +
                  stat_ellipse(level = 0.95, alpha = 0.8) +
                  scale_color_manual(values = dynamic_colors) + 
                  geom_text_repel(aes(label = sample)) +
                  labs(
                    x = parse(text = x_lab),  # 解析数学表达式
                    y = parse(text = y_lab),
                    title = "PLS-DA Score Plot"
                  ) +
                  theme_bw()+theme(plot.title = element_text(hjust = 0.5))
                
                ggsave(g_plsda,filename = paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, ".pdf"), width = 10, height = 6)
                ggsave(g_plsda,filename = paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, ".png"), width = 10, height = 6)
                
              }
              
              sum_plsda[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                              type="PLS-DA", 
                                                                                              A=nrow(pls_scores), 
                                                                                              N=nrow(plsda_model@modelDF), 
                                                                                              `R2X(cum)`=plsda_model@summaryDF$`R2X(cum)`, 
                                                                                              `R2Y(cum)`=plsda_model@summaryDF$`R2Y(cum)`, 
                                                                                              `Q2(cum)`=plsda_model@summaryDF$`Q2(cum)`)
              
              # 从模型对象提取置换检验结果
              perm_data <- plsda_model@suppLs$permMN %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("perm_id") %>% 
                tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                    names_to = "metric", 
                                    values_to = "value")
              
              # 计算原始模型指标
              orig_r2y <- plsda_model@summaryDF$`R2Y(cum)`
              orig_q2 <- plsda_model@summaryDF$`Q2(cum)`
              
              # 生成置换检验图
              perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
                # 绘制置换散点
                geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
                scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
                # 添加回归线
                geom_smooth(
                  aes(group = metric, color = metric),
                  method = "lm", 
                  formula = y ~ x,
                  se = FALSE, 
                  fullrange = TRUE,  # 启用全范围延伸
                  # linewidth = 1.2,   # 加粗线宽匹配图片
                  linetype = "dashed"
                ) +
                # 坐标轴设置（匹配图片范围）
                scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
                scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
                # 标签系统（完全复现图片文字）
                labs(title = "permutation_test",
                     subtitle = paste0("R2Y(cum)=", orig_r2y, "; Q2(cum)=", orig_q2),
                     x = "200 permutations 1 components",
                     y = "value",
                     color = "Metric") +
                # 主题优化（精确字体和网格匹配）
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 10, color = "black"),
                      legend.position = "right",
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
              
              # 输出高清图
              ggsave(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
              
              #### 运行OPLS-DA模型 ----
              oplsda_data <- metabo.data.fill[, which(group %in% c(experimental, control))]
              
              # 获取分组信息
              oplsda_group <- group[which(group %in% c(experimental, control))]
              
              NcrossvalI <- min(7,ncol(oplsda_data))
              
              oplsda_model <- opls(t(oplsda_data), oplsda_group, orthoI=NA, crossvalI=NcrossvalI, permI = 200)
              
              if (nrow(oplsda_model@modelDF)==0) {
                
                oplsda_model <- opls(t(oplsda_data), oplsda_group, predI = 1, orthoI=1, crossvalI=NcrossvalI, permI = 200)
                
              }
              
              # 从模型对象中提取得分矩阵
              opls_scores <- data.frame(
                p1 = oplsda_model@scoreMN[, 1],       # 预测主成分得分
                orth1 = oplsda_model@orthoScoreMN[,1], # 正交主成分得分
                Groups = oplsda_group,
                sample = rownames(t(oplsda_data))
              )
              
              x_lab <- paste0("t[1]~(", round(oplsda_model@modelDF$R2X[1]*100, 1), "*'%')")
              y_lab <- paste0("to[1]~(", round(oplsda_model@modelDF$R2X[2]*100, 1), "*'%')")
              
              g_oplsda <- ggplot(opls_scores, aes(x = p1, y = orth1, color = Groups)) +
                geom_point(size = 4) + 
                scale_color_manual(values = dynamic_colors) + 
                stat_ellipse(level = 0.95, alpha = 0.8) +
                geom_text_repel(aes(label = sample)) +
                labs(
                  x = parse(text = x_lab),  # 解析数学表达式
                  y = parse(text = y_lab),
                  title = "OPLS-DA Score Plot"
                ) +
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))
              
              ggsave(g_oplsda, filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, ".pdf"), width = 10, height = 6)
              ggsave(g_oplsda, filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, ".png"), width = 10, height = 6)
              
              sum_oplsda[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                               type="OPLS-DA", 
                                                                                               A=nrow(opls_scores), 
                                                                                               N=paste0(sum(grepl("^p",rownames(oplsda_model@modelDF))), "+", sum(grepl("^o",rownames(oplsda_model@modelDF)))), 
                                                                                               `R2X(cum)`=oplsda_model@summaryDF$`R2X(cum)`, 
                                                                                               `R2Y(cum)`=oplsda_model@summaryDF$`R2Y(cum)`, `Q2(cum)`=oplsda_model@summaryDF$`Q2(cum)`)
              
              # 从模型对象提取置换检验结果
              perm_data <- oplsda_model@suppLs$permMN %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("perm_id") %>% 
                tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                    names_to = "metric", 
                                    values_to = "value")
              
              # 计算原始模型指标
              orig_r2y <- oplsda_model@summaryDF$`R2Y(cum)`
              orig_q2 <- oplsda_model@summaryDF$`Q2(cum)`
              
              # 生成置换检验图
              perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
                # 绘制置换散点
                geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
                scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
                # 添加回归线
                geom_smooth(
                  aes(group = metric, color = metric),
                  method = "lm", 
                  formula = y ~ x,
                  se = FALSE, 
                  fullrange = TRUE,  # 启用全范围延伸
                  # linewidth = 1.2,   # 加粗线宽匹配图片
                  linetype = "dashed"
                ) +
                # 坐标轴设置（匹配图片范围）
                scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
                scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
                # 标签系统（完全复现图片文字）
                labs(title = "permutation_test",
                     subtitle = paste0("R2Y(cum)=", orig_r2y, "; Q2(cum)=", orig_q2),
                     x = "200 permutations 1 components",
                     y = "value",
                     color = "Metric") +
                # 主题优化（精确字体和网格匹配）
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 10, color = "black"),
                      legend.position = "right",
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
              
              # 输出高清图
              ggsave(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
              
              # S-Plot
              loadings <- oplsda_model@loadingMN[, 1]  # 第一预测成分的载荷（协方差）
              score_p1 <- oplsda_model@scoreMN[, 1]      # 第一预测成分得分
              vip_values <- oplsda_model@vipVn           # VIP值
              
              # 提取变量名称（确保与载荷、标准差、VIP的名称一致）
              variable_names <- names(loadings)
              
              p_cov <- cov(score_p1, t(oplsda_data[names(loadings),]))%>%as.numeric()
              p_corr <- cor(score_p1, t(oplsda_data[names(loadings),]))%>%as.numeric()
              
              # 构建S-Plot数据框
              splot_data <- data.frame(
                Variable = variable_names,
                Loading = p_cov,
                Correlation = p_corr,
                VIP = vip_values
              )
              
              # 标记VIP > 1.0的变量
              splot_data$Significant <- ifelse(splot_data$VIP > 1.0, "Yes", "No")
              
              # 绘制S-Plot
              splot <- ggplot(splot_data, aes(x = Loading, y = Correlation)) +
                # 绘制点，按VIP显著性着色
                geom_point(aes(color = Significant), size = 2.5, alpha = 0.7) +
                # 设置颜色方案（红色：显著；灰色：不显著）
                scale_color_manual(values = c("Yes" = "red", "No" = "gray50")) +
                # 添加参考线（x=0, y=0）
                geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
                # 设置主题和标签
                theme_minimal() +
                labs(
                  title = "OPLS-DA S-Plot",
                  x = "Loading [p1]",
                  y = "Correlation [p1]",
                  color = "VIP > 1.0"
                ) +
                theme(
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.title = element_text(size = 10),
                  legend.position = "bottom"
                )
              
              # 标记VIP值最高的前5个变量
              top_n <- 5  # 可调整为你需要的数量
              top_vars <- splot_data %>% 
                arrange(desc(VIP)) %>% 
                head(top_n)
              
              if (nrow(top_vars) > 0) {
                splot <- splot +
                  geom_text(data = top_vars, 
                            aes(label = Variable), 
                            hjust = 0.5, vjust = -0.7,  # 标签位置在点下方
                            size = 3, color = "black",
                            check_overlap = TRUE  # 避免标签重叠
                  )
              }
              
              # 保存图形
              ggsave(
                plot = splot,
                filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_", experimental, "_vs_", control, "_", ion_type, ".pdf"),
                width = 8, height = 6
              )
              ggsave(
                plot = splot,
                filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_", experimental, "_vs_", control, "_", ion_type, ".png"),
                width = 8, height = 6, dpi = 300
              )
              
              # 保存S-Plot数据
              write.csv(
                splot_data,
                paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_Data_", experimental, "_vs_", control, "_", ion_type, ".csv"),
                row.names = FALSE
              )
              
            }
            
            # 差异筛选 ----
            dir.create("5.差异代谢物分析")
            
            dataMat <- metabo.data.fill
            
            for (g in 1:length(group.data)){
              
              cmp_name <- group.data[g]
              expGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
              ctrlGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
              
              # 创建文件夹（如果不存在）用于保存当前比较的结果
              output_dir <- paste0("5.差异代谢物分析/1.分组对比/",cmp_name,"/0.差异筛选表格")
              dir.create(output_dir, recursive = TRUE, showWarnings = F)
              
              subData <- dataMat[, which(sample.data[[ion_type]]$Condition %in% c(expGroup, ctrlGroup))]
              
              # 分别提取实验组和对照组的样本名称（在subData中）
              expSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == expGroup]
              ctrlSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == ctrlGroup]
              
              # 计算每个代谢物的t检验p值（在当前比较样本上）
              pvals <- apply(subData, 1, function(x) {
                
                tryCatch({
                  
                  if (input$pvalue_type=="ttest") {
                    
                    t_res <- t.test(x[expSamples], x[ctrlSamples])
                    
                  }
                  
                  if (input$pvalue_type=="wilcox_test") {
                    
                    t_res <- wilcox.test(x[expSamples], x[ctrlSamples])
                    
                  }
                  
                  return(t_res$p.value)
                  
                },error=function(e){
                  
                  return(1)
                  
                })
                
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
              
              tryCatch({
                
                # 利用ropls包进行PLS-DA，计算VIP值（只取预测成分predI=1）
                selectedSamples <- c(expSamples, ctrlSamples)
                X <- t(subData[, selectedSamples])
                y <- factor(c(rep("exp", length(expSamples)), rep("ctrl", length(ctrlSamples))))
                model <- opls(X, y, predI = 1, orthoI = NA, crossvalI=ifelse(length(selectedSamples)<7,length(selectedSamples),7))
                
                if (nrow(model@modelDF)==0) {
                  
                  model <- opls(X, y, predI = 1, orthoI = 1, crossvalI=ifelse(length(selectedSamples)<7,length(selectedSamples),7))
                  
                }
                
                vip <- getVipVn(model)
                vip <- vip[rownames(subData)]
                
              },error=function(e){
                
                vip <<- 0
                
              })
              
              # 整合统计结果和当前比较的定量数据
              res_df <- data.frame(peak_name = rownames(subData),
                                   name = metabo.data.raw[[ion_type]]$name[match(rownames(subData), metabo.data.raw[[ion_type]]$peak_name)],
                                   pvalue = pvals,
                                   adj_p_value = p.adjust(pvals, method = input$padjust_method),
                                   FC = fc,
                                   Log2FC = log2(fc),
                                   VIP = vip,
                                   stringsAsFactors = FALSE)
              
              # 增加上下调标志
              res_df[["difference"]] <- "nodiff"
              
              # fc筛选
              if (input$fc_cutoff_judge & !input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
              }
              
              # vip筛选
              if (!input$fc_cutoff_judge & input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
              }
              
              # p筛选
              if (!input$fc_cutoff_judge & !input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge) {
                  res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$adj_p_value < input$p_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$adj_p_value < input$p_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$pvalue < input$p_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$pvalue < input$p_cutoff)] <- "down"
                }
              }
              
              # fc和p值筛选
              if (input$fc_cutoff_judge & !input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
                }
              }
              
              # vip和p值筛选
              if (!input$fc_cutoff_judge & input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
                }
              }
              
              # fc和vip筛选
              if (input$fc_cutoff_judge & input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
              }
              
              # fc，vip，p值筛选
              if (input$fc_cutoff_judge & input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
                }
              }
              
              # 匹配注释
              selected_columns_index <- c("peak_name", "mz", "rt", "id_kegg", "hmdb_kegg", "lmsd_kegg")
              # 去除不存在的列索引（防止有些列缺失）
              selected_columns_index <- selected_columns_index[selected_columns_index %in% colnames(metabo.data.raw[[ion_type]])]
              anno_metabo_data <- metabo.data.raw[[ion_type]][,selected_columns_index]
              anno_metabo_data$peak_name <- as.character(anno_metabo_data$peak_name)
              res_df <- cbind(left_join(res_df, anno_metabo_data), subData)
              
              # 添加Mode标记
              res_df$peak_name <- paste0(res_df$peak_name, "_", ion_type)
              
              if (!input$single_ion_judge) {
                
                colnames(res_df) <- str_replace_all(colnames(res_df), "(\\.POS)|(\\.pos)|(\\.NEG)|(\\.neg)", "")
                
              }
              
              results_mode[[cmp_name]][[ion_type]] <- res_df
              
              results_mode_diff[[cmp_name]][[ion_type]] <- res_df_diff <- res_df[which(res_df$difference != "nodiff"),]
              
              results_mode_diff_name[[cmp_name]][[ion_type]] <- res_df_diff_name  <- res_df[which(res_df$difference != "nodiff" & !is.na(res_df$name)),]
              
              summary_stats_list[[cmp_name]][[ion_type]] <- data.frame(compare=cmp_name,
                                                                       ion_type=ion_type,
                                                                       raw_feature_num = nrow(metabo.data),
                                                                       remove_missing_num = nrow(metabo.data2),
                                                                       preprocess_feature_num=nrow(res_df),
                                                                       preprocess_feature_num_name=sum(!is.na(res_df$name)),
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
              saveWorkbook(wb, paste0(output_dir, "/", cmp_name, "_组间对比_", ion_type, ".xlsx"), overwrite = TRUE)
              
            }
            
          } else {
            
            summary_stats_list[["no comapre"]][[ion_type]] <- data.frame(compare="no compare",
                                                                         ion_type=ion_type,
                                                                         raw_feature_num = nrow(metabo.data),
                                                                         remove_missing_num = nrow(metabo.data2),
                                                                         preprocess_feature_num=nrow(metabo.data.fill),
                                                                         preprocess_feature_num_name=sum(!is.na(metabo.data.raw[[ion_type]]$name[match(rownames(metabo.data.fill), metabo.data.raw[[ion_type]]$peak_name)])))
            
          }
          
        }
        
        # 保存结果rda数据
        save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw"), file = "result.rda")
        # load(file = "result.rda")
        # ion_type = "NEG"
        # ion_type = "single"
        # sum_pca <- read.xlsx(paste0("4.多元统计分析/1.PCA/PCA_summary.xlsx"))
        # sum_plsda <- read.xlsx(paste0("4.多元统计分析/2.PLS-DA/PLS-DA_summary.xlsx"))
        # sum_oplsda <- read.xlsx(paste0("4.多元统计分析/3.OPLS-DA/OPLS-DA_summary.xlsx"))
        
        # 判断是否除QC外有多组数据
        if (length(unique(sample.data[[ion_type]][["Condition"]][!grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T)]))>1){
          
          sum_pca <- do.call(rbind, sum_pca)
          write.xlsx(sum_pca, paste0("4.多元统计分析/1.PCA/PCA_summary.xlsx"))
          
          sum_plsda <- do.call(rbind, sum_plsda)
          write.xlsx(sum_plsda, paste0("4.多元统计分析/2.PLS-DA/PLS-DA_summary.xlsx"))
          
          sum_oplsda <- do.call(rbind, sum_oplsda)
          write.xlsx(sum_oplsda, paste0("4.多元统计分析/3.OPLS-DA/OPLS-DA_summary.xlsx"))
          
          # 各组NEG+POS合并后的生信 ----
          for (cmp_name in group.data) {
            
            experimental <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
            control <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
            
            # 合并当前比较中POS和NEG的结果（行合并，每行为一个代谢物的统计结果）
            if (!input$single_ion_judge) {
              
              for (i in 1:length(results_mode[[cmp_name]])) {
                
                colnames(results_mode[[cmp_name]][[i]]) <- colnames(results_mode[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                colnames(results_mode_diff[[cmp_name]][[i]]) <- colnames(results_mode_diff[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                colnames(results_mode_diff_name[[cmp_name]][[i]]) <- colnames(results_mode_diff_name[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                
              }
              
            }
            
            results_mode_merged <- do.call(rbind, results_mode[[cmp_name]])
            results_mode_diff_merged <- do.call(rbind, results_mode_diff[[cmp_name]])
            results_mode_diff_name_merged <- do.call(rbind, results_mode_diff_name[[cmp_name]])
            
            # 保存结果到当前比较的文件夹
            wb <- createWorkbook()
            addWorksheet(wb, "组间对比")
            addWorksheet(wb, "组间对比差异")
            addWorksheet(wb, "组间对比差异name")
            writeData(wb, "组间对比", results_mode_merged)
            writeData(wb, "组间对比差异", results_mode_diff_merged)
            writeData(wb, "组间对比差异name", results_mode_diff_name_merged)
            saveWorkbook(wb, paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/0.差异筛选表格/", cmp_name, "_组间对比_merged.xlsx"), overwrite = TRUE)
            
            if (input$single_ion_judge) {
              
              sample.data[["POS"]] <- sample.data[[ion_type]]
              
              sample.data[["POS"]]$`File name2` <- sample.data[["POS"]]$`File name`
              
            }else{
              
              sample.data[["POS"]]$`File name2` <- sample.data[["POS"]]$`File name`%>%str_replace_all("(POS)|(pos)","")%>%str_replace_all("\\.\\.","\\.")
              
            }
            
            expSamples <- sample.data[["POS"]]$`File name2`[sample.data[["POS"]]$Condition == experimental]
            ctrlSamples <- sample.data[["POS"]]$`File name2`[sample.data[["POS"]]$Condition == control]
            sample_cols <- c(expSamples, ctrlSamples)
            
            ## 差异的代谢物（feature）---------------------------
            if (nrow(results_mode_diff_merged)>0) {
              
              expr_data <- results_mode_diff_merged[, sample_cols]
              rownames(expr_data) <- results_mode_diff_merged$peak_name
              
              ### 1. HCA（层级聚类热图和表格） ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/1.HCA")
              dir.create(output_folder, recursive = T)
              
              # 对行（feature）进行标准化（Z-score），便于聚类比较
              # expr_data_scaled <- expr_data
              expr_data_scaled <- t(scale(t(expr_data)))
              
              hca_group <- c(sample.data[[ion_type]]$Condition[which(sample.data[[ion_type]]$Condition == experimental)],
                             sample.data[[ion_type]]$Condition[which(sample.data[[ion_type]]$Condition == control)])
              annotation_col <- data.frame(Group = hca_group)
              rownames(annotation_col) <- colnames(expr_data_scaled) # 确保样本名对齐
              
              # 获取数据范围（假设已标准化，范围对称）
              data_range <- range(expr_data_scaled)
              max_abs <- min(data_range[2],3,na.rm = T)
              min_abs <- max(data_range[1],-3,na.rm = T)
              
              # 定义对称的数值断点（覆盖数据范围，0居中）
              bk <- seq(min_abs, max_abs, length.out = 100)
              
              # 定义颜色渐变（100级）
              heatmap_colors <- colorRampPalette(c("#0072B2", "white", "#F08080"))(100)
              
              # 绘制热图
              pheatmap::pheatmap(expr_data_scaled, 
                                 # scale = "row",
                                 cluster_cols = F,
                                 color = heatmap_colors,
                                 border_color = NA,
                                 breaks = bk,
                                 clustering_distance_rows = "euclidean", 
                                 clustering_method = "complete", 
                                 main = "Heatmap of Differential Features",
                                 filename = file.path(output_folder, "heatmap.pdf"),
                                 height = max(0.1*nrow(expr_data),15),
                                 show_rownames=F, show_colnames=F,
                                 annotation_col = annotation_col, # 添加列分组注释
                                 annotation_names_col = TRUE)
              
              # 计算层级聚类（行聚类）
              hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
              # 保存聚类顺序表格
              cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
              write.xlsx(cluster_order, file = file.path(output_folder, "HCA.xlsx"), rowNames = T)
              
              ### 2. Correlation Analysis ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/2.correlation")
              
              dir.create(output_folder, recursive = T)
              
              # 这里计算feature之间的相关性
              cor_mat <- cor(t(expr_data), use = "pairwise.complete.obs")
              
              if (nrow(cor_mat)>20) {
                
                # 步骤1：计算特征的综合相关性评分（取绝对值均值）
                feature_scores <- rowMeans(abs(cor_mat), na.rm = TRUE)
                
                # 步骤2：按评分降序排列，提取前20个特征名
                top_features <- names(sort(feature_scores, decreasing = TRUE)[1:20])
                
                # 步骤3：从原始矩阵中筛选出这20个特征的相关性子矩阵
                cor_mat_top <- cor_mat[top_features, top_features]
                
                # 步骤4：绘制筛选后的热图（调整标签可读性）
                pdf(file.path(output_folder, "correlation_top20_features.pdf"))
                corrplot(cor_mat_top,
                         method = "circle",
                         type = "upper",
                         tl.col = "black",
                         tl.cex = 0.6,      # 增大标签字号方便查看
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Top 20 Features by Mean Absolute Correlation",
                         mar = c(0,0,1,0),
                         order = "hclust",  # 按层次聚类排序[2,3](@ref)
                         addrect = 2)        # 添加聚类分界线[2,3](@ref)
                dev.off()
                
              }else{
                
                # 绘制相关性热图
                pdf(file.path(output_folder, "correlation_features.pdf"))
                corrplot(cor_mat, method = "circle", tl.col = "black", tl.cex = 0.6,
                         type="upper",
                         # addCoef.col = "black",
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Correlation Among Features", mar = c(0,0,1,0))
                dev.off()
                
                cor_mat_top <- cor_mat
                
              }
              
              # 保存相关性矩阵
              write.xlsx(cor_mat_top, file = file.path(output_folder, "correlation_features.xlsx"), rowNames = TRUE)
              
              ### 3. Volcano Plot ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/3.volcano")
              dir.create(output_folder, recursive = T)
              
              # 计算 log2(FC)
              diff_features_all <- results_mode_merged %>% dplyr::mutate(negLog10P = -log10(pvalue))
              
              # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
              max_val <- ceiling(max(abs(diff_features_all$Log2FC[is.finite(diff_features_all$Log2FC)]), na.rm = TRUE))
              
              # volcano_plot <- ggplot(diff_features_all, 
              #                        aes(x = Log2FC, 
              #                            y = VIP, 
              #                            color = negLog10P,  # 将color映射到VIP值
              #                            shape = difference)) +  # shape映射调控方向
              #   geom_point(alpha = 0.6, size = 3) +  # 固定点大小，通过shape区分
              #   # 颜色梯度设置
              #   scale_color_gradientn(
              #     colors = c("blue", "pink", "red"),  # 蓝-白-红渐变[6,14](@ref)
              #     # values = scales::rescale(c(0, -log10(0.05), -log10(0.01))),  # 标准化梯度位置
              #     name = "-log10P",  # 颜色图例标题
              #     # limits = c(0, max(diff_features_all$negLog10P)),  # 设置颜色范围
              #     na.value = "grey50"  # 缺失值颜色
              #   ) +
              #   # 形状映射设置
              #   scale_shape_manual(
              #     values = c("up" = 24,   # 上三角
              #                "down" = 25, # 下三角
              #                "nodiff" = 21),  # 圆形[9,14](@ref)
              #     labels = c("Down-regulated", 
              #                "Not significant", 
              #                "Up-regulated"),
              #     name = "Regulation"  # 形状图例标题
              #   ) +
              #   # 参考线保持原样
              #   geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)), linetype = "dashed", color = ifelse(input$fc_cutoff_judge, "black", NA)) +
              #   geom_hline(yintercept = input$vip_cutoff, linetype = "dashed", color = ifelse(input$vip_cutoff_judge, "black", NA)) +
              #   theme_minimal() +
              #   labs(x = "log2(Fold Change)", y = "VIP") +
              #   scale_x_continuous(limits = c(-max_val, max_val)) +
              #   # 增强可视化效果
              #   guides(
              #     color = guide_colorbar(barwidth = 1, barheight = 10),  # 颜色图例垂直排列
              #     shape = guide_legend(override.aes = list(color = "black"))  # 形状图例统一颜色
              #   ) +
              #   theme(
              #     legend.position = "right",
              #     legend.box = "vertical"  # 图例垂直排列
              #   )
              
              if (input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference, size=VIP)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (input$p_cutoff_judge & !input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (!input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = VIP, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = input$vip_cutoff, 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "VIP") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              ggsave(filename = file.path(output_folder, "volcano.png"), volcano_plot, width = 7, height = 6)
              ggsave(filename = file.path(output_folder, "volcano.pdf"), volcano_plot, width = 7, height = 6)
              write.xlsx(diff_features_all, file.path(output_folder, "volcano.xlsx"))
              
              ### 4. VIP分析 ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/4.VIP_analysis")
              dir.create(output_folder, recursive = T)
              
              # 绘制VIP分布直方图
              vip_hist <- ggplot(results_mode_merged, aes(x = VIP)) + 
                geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black", alpha = 0.8) +
                geom_vline(xintercept = input$vip_cutoff, linetype = "dashed", color = "red") +
                theme_minimal() +
                labs(title = "VIP Value Distribution", x = "VIP", y = "Count")+
                theme(plot.title = element_text(hjust = 0.5))
              
              ggsave(file.path(output_folder, "VIP_distribution.png"), vip_hist, width = 7, height = 6)
              ggsave(file.path(output_folder, "VIP_distribution.pdf"), vip_hist, width = 7, height = 6)
              
              # 1. 设置 bin 宽度（和图一致）
              binwidth <- 0.2
              
              # 2. 创建 bin 列（分组）
              vip_binned <- results_mode_merged %>%
                dplyr::mutate(VIP_bin = cut(VIP,
                                            breaks = seq(floor(min(VIP, na.rm = TRUE)),
                                                         ceiling(max(VIP, na.rm = TRUE)),
                                                         by = binwidth),
                                            include.lowest = TRUE, right = FALSE)) %>%
                dplyr::count(VIP_bin, name = "Count")
              
              # 3. 保存为 Excel
              write.xlsx(vip_binned, file = file.path(output_folder, "VIP_distribution_table.xlsx"))
              
              # 绘制气泡图
              if (nrow(results_mode_diff_name_merged)>20) {
                vip_data <- dplyr::arrange(results_mode_diff_merged, desc(VIP))%>%.[1:20,]
              }else{
                vip_data <- dplyr::arrange(results_mode_diff_merged, desc(VIP))
              }
              
              custom_colors <- c("up" = "#F08080", "down" = "#0072B2") 
              vip_data$peak_name <- factor(vip_data$peak_name, levels = unique(vip_data$peak_name[order(vip_data$VIP)]))
              vip_bubble <- ggplot(vip_data, aes(x = VIP, y = peak_name, size = abs(Log2FC), color = difference)) +
                geom_point(alpha = 0.7) +
                scale_size_continuous(range = c(1, 10)) +
                scale_color_manual(values = custom_colors) +
                labs(x = 'VIP', y = 'peak_name', title = 'VIP_distribution', size = "Log2FC") +
                theme_minimal() +
                theme(
                  panel.grid = element_blank(),  # 去除网格线
                  panel.border = element_rect(    # 添加坐标轴四边的框线
                    color = "black", 
                    fill = NA, 
                    linewidth = 0.5
                  ),
                  plot.title = element_text(hjust = 0.5)  # 将标题设置为居中
                )
              
              ggsave(file.path(output_folder, "VIP_bubble.png"), vip_bubble, width = 7, height = 6)
              ggsave(file.path(output_folder, "VIP_bubble.pdf"), vip_bubble, width = 7, height = 6)
              write.xlsx(vip_data, file.path(output_folder, "VIP_bubble.xlsx"), rowNames = FALSE)
              
              ### 6.差异网络分析 ---------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/6.差异网络分析")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              if (nrow(expr_data)<5000) {
                  
                # 计算代谢物间的相关系数矩阵
                cor_mat2 <- cor_mat <- cor(t(expr_data), method = "spearman")  # 使用Spearman相关系数
                
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
                  dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
                
                cor_df2 <- as.data.frame(as.table(cor_mat2)) %>% 
                  dplyr::rename(Metabolite1 = Var1, Metabolite2 = Var2, Correlation = Freq) %>%
                  dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
                
                diff_features <- results_mode_diff_merged %>% dplyr::mutate(regulation = difference)
                
                diff_features$regulation <- factor(diff_features$regulation, levels = c("up", "down"))
                
                diff_features$truncated_name <- truncated_names
                
                if (nrow(diff_features)>100) {
                  
                  # 按VIP排序并取前100
                  top_features <- diff_features %>%
                    dplyr::arrange(desc(VIP)) %>%
                    head(100)
                  
                  # 筛选cor_df只包含这100个代谢物的边
                  cor_df <- cor_df %>%
                    dplyr::filter(Metabolite1 %in% top_features$truncated_name & 
                                    Metabolite2 %in% top_features$truncated_name)
                  
                  cor_df2 <- cor_df2 %>%
                    dplyr::filter(Metabolite1 %in% top_features$name2 & 
                                    Metabolite2 %in% top_features$name2)
                  
                  # 更新diff_features为筛选后的100个
                  diff_features <- top_features
                  
                }
                
                # 确保代谢物名称唯一
                diff_features_unique <- diff_features %>% dplyr::distinct(peak_name, .keep_all = TRUE)  # 去重
                
                # 创建igraph图对象
                g <- graph_from_data_frame(
                  d = cor_df[, 1:3],  # 边数据（Metabolite1, Metabolite2, Correlation）
                  directed = T,
                  vertices = diff_features_unique %>% dplyr::select(truncated_name, Log2FC, regulation)  # 节点属性
                )
                
                # 节点属性设置
                V(g)$size <- abs(V(g)$Log2FC) * 3  # 大小反映Log2FC绝对值
                V(g)$color <- ifelse(V(g)$regulation == "up", "#E64B35FF", "#3C5488FF")  # 上调红，下调蓝
                
                # 边属性设置
                E(g)$width <- abs(E(g)$Correlation) * 0.2  # 边宽与相关性强度正相关
                E(g)$color <- ifelse(E(g)$Correlation > 0, "#E64B35FF", "#3C5488FF")  # 正相关红，负相关蓝
                
                # 生成带图例的网络图
                set.seed(123)
                network_plot <- ggraph(g, layout = "fr") + 
                  # 边图层（颜色和宽度）
                  geom_edge_link(
                    aes(edge_width = abs(Correlation),  # 边宽度映射到相关性绝对值
                        edge_color = ifelse(Correlation > 0, "Positive", "Negative")),  # 边颜色映射到正负
                    alpha = 0.6) +
                  # 节点图层（颜色和大小）
                  geom_node_point(
                    aes(color = regulation,    # 颜色映射到上下调
                        size = abs(Log2FC)),   # 大小映射到Log2FC绝对值
                    alpha = 0.8) +
                  # 文本标签（防重叠）
                  geom_node_text(
                    aes(label = name), 
                    size = 3, 
                    repel = TRUE,
                    color = "black") +
                  # 边颜色映射（正红负蓝）
                  scale_edge_color_manual(
                    name = "Correlation Type",
                    values = c("Positive" = "#E64B35FF", "Negative" = "#3C5488FF")
                  ) +
                  # 边宽度映射（连续型）
                  scale_edge_width_continuous(
                    name = "|Correlation|", 
                    range = c(0.5, 4),  # 边宽度范围
                    breaks = c(0.1, 0.4, 0.7)  # 图例刻度
                  ) +
                  # 节点颜色映射（离散型）
                  scale_color_manual(
                    name = "Regulation",
                    values = c("up" = "#E64B35FF", "down" = "#3C5488FF")
                  ) +
                  # 节点大小映射（连续型）
                  scale_size_continuous(
                    name = "|Log2FC|", 
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
                write.graph(g, file.path(output_folder, "feature_network.graphml"), format = "graphml")
                
                # 正确生成节点属性表的方法
                node_table <- diff_features %>%
                  # 首先确保只选择实际存在的列
                  dplyr::select(
                    truncated_name,  # 截断后的代谢物名称（必须存在于diff_features中）
                    Log2FC,         # 差异倍数
                    VIP,            # VIP值
                    regulation      # 上下调标记
                  ) %>%
                  # 添加原始名称（如果feature_names是独立向量）
                  dplyr::mutate(
                    Original_Name = rownames(diff_features),  # 假设原始名称是行名
                    Direction = ifelse(regulation == "up", "Upregulated", "Downregulated"),
                    Metabolite = truncated_name  # 重命名列
                  ) %>%
                  # 重新排列列顺序
                  dplyr::select(
                    Metabolite,
                    Original_Name,
                    Log2FC,
                    VIP,
                    Regulation = regulation,
                    Direction
                  )
                
                # 生成边属性表（已筛选相关性强的边）
                edge_table <- cor_df %>%
                  dplyr::mutate(
                    Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative"),
                    Abs_Correlation = abs(Correlation)  # 相关性绝对值
                  ) %>%
                  dplyr::select(
                    Source = Metabolite1,        # 源节点
                    Target = Metabolite2,        # 目标节点
                    Correlation,                 # 相关系数
                    Correlation_Type,            # 正/负相关
                    Abs_Correlation              # 相关性强度
                  )
                
                # 计算网络拓扑指标
                node_topology <- data.frame(
                  Metabolite = V(g)$name,
                  Degree = igraph::degree(g),            # 连接度
                  Betweenness = igraph::betweenness(g),  # 中介中心性
                  Closeness = igraph::closeness(g)       # 接近中心性
                )
                
                # 合并节点属性与拓扑指标
                network_summary <- node_table %>%
                  dplyr::left_join(node_topology, by = c("Metabolite" = "Metabolite"))
                
                # 创建新的Excel工作簿
                wb <- createWorkbook()
                
                # 1. 添加节点属性表（Sheet1）
                addWorksheet(wb, sheetName = "节点属性")
                writeData(wb, sheet = "节点属性", x = node_table, startCol = 1, startRow = 1)
                
                # 2. 添加边属性表（Sheet2）
                addWorksheet(wb, sheetName = "边属性")
                writeData(wb, sheet = "边属性", x = edge_table, startCol = 1, startRow = 1)
                
                # 3. 添加网络拓扑表（Sheet3，可选）
                addWorksheet(wb, sheetName = "网络拓扑")
                writeData(wb, sheet = "网络拓扑", x = network_summary, startCol = 1, startRow = 1)
                
                # 保存Excel文件
                saveWorkbook(wb, file = file.path(output_folder, "feature_network.xlsx"), overwrite = TRUE)
                
              }
              
              ### 7.差异代谢物的FC Rank图 ----
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/7.差异代谢物的Rank图")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 计算Rank，即FC从小到大的排序值
              rank_data <- results_mode_merged
              rank_data <- rank_data %>%
                dplyr::arrange(FC) %>%
                dplyr::mutate(Rank = row_number())
              
              # 找出Rank排名最前面和最末尾的10个点
              top_bottom_10 <- rank_data %>%
                dplyr::arrange(Rank) %>%
                dplyr::filter(Rank %in% c(1:10, (nrow(rank_data) - 9):nrow(rank_data)))
              
              # 自定义颜色向量，这里假设difference列有不同的类别，你可以根据实际情况修改颜色
              custom_colors <- c("up" = "#F08080", "down" = "#0072B2", "nodiff" = "gray")
              
              # 绘制散点图
              p_rank <- ggplot(rank_data, aes(x = Rank, y = Log2FC, color = difference)) +
                geom_point() +
                geom_text_repel(data = top_bottom_10, aes(label = peak_name), 
                                segment.color = "black", segment.size = 0.5, size = 3,
                                max.overlaps = Inf,        # 不限制标注数量
                                force = 5,                 # 增大排斥力
                                box.padding = 0.5,         # 文字与文字之间的间距
                                point.padding = 0.3        # 文字与点的间距
                ) +
                scale_color_manual(values = custom_colors) +
                theme_minimal() +
                # 去掉背景网格线
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank()) +
                # 添加纵坐标为1和 -1 的参考线
                geom_hline(yintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)), linetype = "dashed", color = "gray") +
                # 为横纵坐标加上边框
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
              
              ggsave(file.path(output_folder, "feature_with_rank.png"), p_rank, width = 8, height = 6)
              ggsave(file.path(output_folder, "feature_with_rank.pdf"), p_rank, width = 8, height = 6)
              
              # 将结果保存为新的 CSV 文件
              write.csv(rank_data, file.path(output_folder, "feature_with_rank.csv"), row.names = FALSE)
              
            }
            
            ## 差异的代谢物（有name的）---------------------------
            if (nrow(results_mode_diff_name_merged)>0) {
              
              expr_data <- results_mode_diff_name_merged[,sample_cols]%>%as.data.frame()
              rownames(expr_data) <- paste0(results_mode_diff_name_merged$peak_name, "_", results_mode_diff_name_merged$name)
              
              ### 1. HCA（层级聚类热图和表格） ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/1.HCA")
              
              # 对行（feature）进行标准化（Z-score），便于聚类比较
              expr_data_scaled <- t(scale(t(expr_data)))
              rownames(expr_data_scaled) <- results_mode_diff_name_merged$name
              
              # 绘制热图
              truncated_labels <- substr(rownames(expr_data_scaled), 1, pmin(nchar(rownames(expr_data_scaled)), 40)) # 最多保留40个字符
              # 添加省略号表示截断
              truncated_labels <- ifelse(
                nchar(rownames(expr_data_scaled)) > 20,
                paste0(truncated_labels, "..."),
                truncated_labels
              )
              
              # 获取数据范围（假设已标准化，范围对称）
              data_range <- range(expr_data_scaled)
              max_abs <- min(data_range[2],3,na.rm = T)
              min_abs <- max(data_range[1],-3,na.rm = T)
              
              # 定义对称的数值断点（覆盖数据范围，0居中）
              bk <- seq(min_abs, max_abs, length.out = 100)
              
              # 定义颜色渐变（100级）
              heatmap_colors <- colorRampPalette(c("#0072B2", "white", "#F08080"))(100)
              
              pheatmap::pheatmap(
                expr_data_scaled,
                # scale = "row",
                cluster_cols = F,
                color = heatmap_colors,
                border_color = NA,
                breaks = bk,
                labels_row = truncated_labels,  # 使用截断后的标签
                show_colnames = F,
                clustering_distance_rows = "euclidean",
                clustering_method = "complete",
                main = "Heatmap of Differential Features",
                height = max(6, nrow(expr_data_scaled)*0.14),
                width = max(8, ncol(expr_data_scaled)*0.25),
                filename = file.path(output_folder, "heatmap_name.pdf"),
                annotation_col = annotation_col, # 添加列分组注释
                annotation_names_col = TRUE)
              
              # 计算层级聚类（行聚类）
              hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
              # 保存聚类顺序表格
              cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
              write.xlsx(cluster_order, file = file.path(output_folder, "heatmap_name.xlsx"), rowNames = T)
              
              if (nrow(expr_data) > 100) {
                
                # 上调：FC > 1，按 FC 降序
                up_features <- results_mode_diff_name_merged %>%
                  dplyr::filter(FC > 1) %>%
                  dplyr::arrange(desc(FC)) %>%
                  dplyr::slice_head(n = 25)
                
                # 下调：FC < 1，按 FC 升序
                down_features <- results_mode_diff_name_merged %>%
                  dplyr::filter(FC < 1) %>%
                  dplyr::arrange(FC) %>%
                  dplyr::slice_head(n = 25)
                
                # 合并上下调
                selected_diff <- bind_rows(up_features, down_features)
                
                # 提取子集表达矩阵
                expr_data_sub <- expr_data[rownames(expr_data) %in% paste0(selected_diff$peak_name, "_", selected_diff$name), ]
                
                # 标准化
                expr_data_scaled <- t(scale(t(expr_data_sub)))
                rownames(expr_data_scaled) <- selected_diff$name
                
                # 截断标签
                truncated_labels <- substr(rownames(expr_data_scaled), 1, pmin(nchar(rownames(expr_data_scaled)), 40))
                truncated_labels <- ifelse(
                  nchar(rownames(expr_data_scaled)) > 20,
                  paste0(truncated_labels, "..."),
                  truncated_labels
                )
                
                # 颜色断点和配色
                data_range <- range(expr_data_scaled)
                max_abs <- min(data_range[2], 3, na.rm = T)
                min_abs <- max(data_range[1], -3, na.rm = T)
                bk <- seq(min_abs, max_abs, length.out = 100)
                heatmap_colors <- colorRampPalette(c("#0072B2", "white", "#F08080"))(100)
                
                # 保存热图 PDF
                pheatmap::pheatmap(
                  expr_data_scaled,
                  cluster_cols = FALSE,
                  color = heatmap_colors,
                  border_color = NA,
                  breaks = bk,
                  labels_row = truncated_labels,
                  show_colnames = FALSE,
                  clustering_distance_rows = "euclidean",
                  clustering_method = "complete",
                  main = "Heatmap of Differential Features",
                  height = max(6, nrow(expr_data_scaled) * 0.14),
                  width = max(8, ncol(expr_data_scaled) * 0.25),
                  filename = file.path(output_folder, "heatmap_name_down&up_top25.pdf"),
                  annotation_col = annotation_col,
                  annotation_names_col = TRUE
                )
                
                # 保存排序数据表格
                hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
                cluster_order <- data.frame(expr_data_scaled[hc_rows$order, ])
                write.xlsx(cluster_order, file = file.path(output_folder, "heatmap_name_down&up_top25.xlsx"), rowNames = TRUE)
                
              }
              
              ### 2. Correlation Analysis ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/2.correlation")
              
              # 这里计算feature之间的相关性
              cor_mat <- cor(t(expr_data), use = "pairwise.complete.obs")
              
              # 获取原始名称（假设行名是特征名）
              feature_names <- results_mode_diff_name_merged$name
              
              # 截断名称（超过20字符则保留前17字符 + ...）
              truncated_names <- ifelse(
                nchar(feature_names) > 20,
                paste0(substr(feature_names, 1, 17), "..."),
                feature_names
              )
              
              # 更新矩阵的行列名
              rownames(cor_mat) <- truncated_names
              colnames(cor_mat) <- truncated_names
              
              if (ncol(cor_mat)>20) {
                
                # 步骤1：计算特征的综合相关性评分（取绝对值均值）
                feature_scores <- rowMeans(abs(cor_mat), na.rm = TRUE)
                
                # 步骤2：按评分降序排列，提取前20个特征名
                top_features <- names(sort(feature_scores, decreasing = TRUE)[1:20])
                
                # 步骤3：从原始矩阵中筛选出这20个特征的相关性子矩阵
                cor_mat_top <- cor_mat[top_features, top_features]
                
                # 步骤4：绘制筛选后的热图（调整标签可读性）
                pdf(file.path(output_folder, "correlation_top20_features_name.pdf"))
                corrplot(cor_mat_top,
                         method = "circle",
                         type = "upper",
                         tl.col = "black",
                         tl.cex = 0.6,      # 增大标签字号方便查看
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Top 20 Features by Mean Absolute Correlation",
                         mar = c(0,0,1,0),
                         order = "hclust",  # 按层次聚类排序[2,3](@ref)
                         addrect = 2)        # 添加聚类分界线[2,3](@ref)
                dev.off()
                
              }else{
                
                # 绘制相关性热图
                pdf(file.path(output_folder, "correlation_features_name.pdf"))
                corrplot(cor_mat, method = "circle", tl.col = "black", type = "upper",
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Correlation Among Samples", mar = c(0,0,1,0))
                dev.off()
                
              }
              
              # 保存相关性矩阵
              write.xlsx(t(expr_data), file = file.path(output_folder, "correlation_features_name.xlsx"), rowNames = TRUE)
              
              ### 3. Volcano Plot ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/3.volcano")
              
              # 计算 log2(FC)
              diff_features_all <- results_mode_merged[which(!is.na(results_mode_merged$name)),] %>% dplyr::mutate(negLog10P = -log10(pvalue))
              
              # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
              max_val <- ceiling(max(abs(diff_features_all$Log2FC[is.finite(diff_features_all$Log2FC)]), na.rm = TRUE))
              
              # volcano_plot <- ggplot(diff_features_all, 
              #                        aes(x = Log2FC, 
              #                            y = VIP, 
              #                            color = negLog10P,  # 将color映射到VIP值
              #                            shape = difference)) +  # shape映射调控方向
              #   geom_point(alpha = 0.6, size = 3) +  # 固定点大小，通过shape区分
              #   # 颜色梯度设置
              #   scale_color_gradientn(
              #     colors = c("blue", "pink", "red"),  # 蓝-白-红渐变[6,14](@ref)
              #     # values = scales::rescale(c(0, -log10(0.05), -log10(0.01))),  # 标准化梯度位置
              #     name = "-log10P",  # 颜色图例标题
              #     # limits = c(0, max(diff_features_all$negLog10P)),  # 设置颜色范围
              #     na.value = "grey50"  # 缺失值颜色
              #   ) +
              #   # 形状映射设置
              #   scale_shape_manual(
              #     values = c("up" = 24,   # 上三角
              #                "down" = 25, # 下三角
              #                "nodiff" = 21),  # 圆形[9,14](@ref)
              #     labels = c("Down-regulated", 
              #                "Not significant", 
              #                "Up-regulated"),
              #     name = "Regulation"  # 形状图例标题
              #   ) +
              #   # 参考线保持原样
              #   geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)), linetype = "dashed", color = ifelse(input$fc_cutoff_judge, "black", NA)) +
              #   geom_hline(yintercept = input$vip_cutoff, linetype = "dashed", color = ifelse(input$vip_cutoff_judge, "black", NA)) +
              #   theme_minimal() +
              #   labs(x = "log2(Fold Change)", y = "VIP") +
              #   scale_x_continuous(limits = c(-max_val, max_val)) +
              #   # 增强可视化效果
              #   guides(
              #     color = guide_colorbar(barwidth = 1, barheight = 10),  # 颜色图例垂直排列
              #     shape = guide_legend(override.aes = list(color = "black"))  # 形状图例统一颜色
              #   ) +
              #   theme(
              #     legend.position = "right",
              #     legend.box = "vertical"  # 图例垂直排列
              #   )
              
              if (input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference, size=VIP)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (input$p_cutoff_judge & !input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (!input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = VIP, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "difference"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = input$vip_cutoff, 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "VIP") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              ggsave(filename = file.path(output_folder, "volcano_name.png"), volcano_plot, width = 7, height = 6)
              ggsave(filename = file.path(output_folder, "volcano_name.pdf"), volcano_plot, width = 7, height = 6)
              write.xlsx(diff_features_all, file.path(output_folder, "volcano_name.xlsx"))
              
              ### 4. VIP分析 ---------------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/4.VIP_analysis")
              dir.create(output_folder, recursive = T)
              
              if (nrow(results_mode_diff_name_merged)>20) {
                vip_data <- dplyr::arrange(results_mode_diff_name_merged, desc(VIP))%>%.[1:20,]
              }else{
                vip_data <- dplyr::arrange(results_mode_diff_name_merged, desc(VIP))
              }
              
              # 截断名称（超过20字符则保留前17字符 + ...）
              truncated_names <- ifelse(
                nchar(vip_data$name) > 20,
                paste0(substr(vip_data$name, 1, 17), "..."),
                vip_data$name
              )
              
              # 更新矩阵的行列名
              vip_data$truncated_names <- truncated_names
              
              # 绘制气泡图
              custom_colors <- c("up" = "#F08080", "down" = "#0072B2") 
              vip_data$truncated_names <- factor(vip_data$truncated_names, levels = unique(vip_data$truncated_names[order(vip_data$VIP)]))
              vip_bubble <- ggplot(vip_data, aes(x = VIP, y = truncated_names, size = abs(Log2FC), color = difference)) +
                geom_point(alpha = 0.7) +
                scale_size_continuous(range = c(1, 10)) +
                scale_color_manual(values = custom_colors) +
                labs(x = 'VIP', y = 'name', title = 'VIP_distribution', size = "Log2FC") +
                theme_minimal() +
                theme(
                  panel.grid = element_blank(),  # 去除网格线
                  panel.border = element_rect(    # 添加坐标轴四边的框线
                    color = "black", 
                    fill = NA, 
                    linewidth = 0.5
                  ),
                  plot.title = element_text(hjust = 0.5)  # 将标题设置为居中
                )
              
              ggsave(file.path(output_folder, "VIP_bubble_name.png"), vip_bubble, width = 7, height = 6)
              ggsave(file.path(output_folder, "VIP_bubble_name.pdf"), vip_bubble, width = 7, height = 6)
              write.xlsx(vip_data, file.path(output_folder, "VIP_bubble_name.xlsx"), rowNames = FALSE)
              
              ### 5.差异代谢物表达量的箱线图 ---------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/5.差异代谢物表达量")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 转换为长数据格式（代谢物ID + 样本 + 定量值 + 分组）
              if (length(unique(results_mode_diff_name_merged$name)) < 24){
                
                expr_long <- results_mode_diff_name_merged[,c("name", sample_cols)] %>%
                  tidyr::pivot_longer(cols = -name, names_to = "sample", values_to = "intensity") %>%
                  dplyr::left_join(sample.data$POS, by = c("sample" = "File name2")) %>%
                  dplyr::mutate(group = factor(Condition, levels = unique(sample.data$POS$Condition)))
                
              }else{
                
                expr_long <- results_mode_diff_name_merged%>%dplyr::arrange(desc(VIP))%>%.[1:24, c("name", sample_cols)] %>%
                  tidyr::pivot_longer(cols = -name, names_to = "sample", values_to = "intensity") %>%
                  dplyr::left_join(sample.data$POS, by = c("sample" = "File name2")) %>%
                  dplyr::mutate(group = factor(Condition, levels = unique(sample.data$POS$Condition)))
                
              }
              
              # (3) 所有代谢物全局分布
              if (input$pvalue_type=="ttest") {
                stat_compare_means_method <- "t.test"
              }else if(input$pvalue_type=="wilcox_test") {
                stat_compare_means_method <- "wilcox.test"
              }
              
              truncate_and_force_wrap <- function(x, max_length = 30, wrap_width = 15) {
                long_idx <- str_length(x) > max_length
                
                # 超长名字先截断加省略号
                x[long_idx] <- paste0(str_sub(x[long_idx], 1, max_length), "...")
                
                # 超长名字才强制换行
                x[long_idx] <- sapply(x[long_idx], function(s) {
                  chars <- strsplit(s, "")[[1]]
                  paste(sapply(seq(1, length(chars), by = wrap_width),
                               function(i) paste0(chars[i:min(i+wrap_width-1, length(chars))], collapse = "")),
                        collapse = "\n")
                })
                
                return(x)
              }
              
              p_global <- ggboxplot(expr_long, x = "group", y = "intensity",color = "group", palette = "npg")+
                stat_compare_means(method = stat_compare_means_method, label = "p.format", vjust = 1, hjust=-0.5)+
                labs(title = "Metabolites Comparison", x = NULL, y = "Normalized Intensity") +
                # theme_classic() +
                # 按代谢物分面
                facet_wrap(~name, scales = "free", ncol = 6, nrow = 4, labeller = labeller(name = function(x) truncate_and_force_wrap(x, max_length = 30, wrap_width = 15)))+
                theme(
                  axis.text.x = element_text(angle = 45, hjust = 1),
                  strip.text = element_text(size = 10, face = "bold"),  # 分面标签样式
                  strip.background = element_blank()     # 分面标签背景
                )
              
              ggsave(
                file.path(output_folder, "Global_Metabolite_Comparison.pdf"), 
                p_global,
                width = 18, 
                height = ceiling(length(unique(expr_long$name))/4) * 2,  # 动态调整高度
                dpi = 300,
                limitsize = FALSE
              )
              
              ggsave(
                file.path(output_folder, "Global_Metabolite_Comparison.png"), 
                p_global,
                width = 18, 
                height = ceiling(length(unique(expr_long$name))/4) * 2,  # 动态调整高度
                dpi = 300,
                limitsize = FALSE
              )
              
              ### 6.差异网络分析 ---------------
              output_folder <- paste0("5.差异代谢物分析/1.分组对比/", cmp_name, "/6.差异网络分析")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 计算代谢物间的相关系数矩阵
              cor_mat2 <- cor_mat <- cor(t(expr_data), method = "spearman")  # 使用Spearman相关系数
              
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
              
              cor_df2 <- as.data.frame(as.table(cor_mat2)) %>% 
                dplyr::rename(Metabolite1 = Var1, Metabolite2 = Var2, Correlation = Freq) %>%
                dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
              
              diff_features <- results_mode_diff_name_merged %>% dplyr::mutate(regulation = difference)
              
              diff_features$regulation <- factor(diff_features$regulation, levels = c("up", "down"))
              
              diff_features$truncated_name <- truncated_names
              
              if (nrow(diff_features) > 100) {
                
                # 按VIP排序并取前100
                top_features <- diff_features %>%
                  dplyr::arrange(desc(VIP)) %>%
                  head(100)
                
                # 筛选cor_df只包含这100个代谢物的边
                cor_df <- cor_df %>%
                  dplyr::filter(Metabolite1 %in% top_features$truncated_name & 
                                  Metabolite2 %in% top_features$truncated_name)
                
                cor_df2 <- cor_df2 %>%
                  dplyr::filter(Metabolite1 %in% top_features$name2 & 
                                  Metabolite2 %in% top_features$name2)
                
                # 更新diff_features为筛选后的100个
                diff_features <- top_features
                
              }
              
              # 确保代谢物名称唯一
              diff_features_unique <- diff_features %>% dplyr::distinct(peak_name, .keep_all = TRUE)  # 去重
              
              # 创建igraph图对象
              g <- graph_from_data_frame(
                d = cor_df[, 1:3],  # 边数据（Metabolite1, Metabolite2, Correlation）
                directed = T,
                vertices = diff_features_unique %>% dplyr::select(truncated_name, Log2FC, regulation)  # 节点属性
              )
              
              # 节点属性设置
              V(g)$size <- abs(V(g)$Log2FC) * 3  # 大小反映Log2FC绝对值
              V(g)$color <- ifelse(V(g)$regulation == "up", "#E64B35FF", "#3C5488FF")  # 上调红，下调蓝
              
              # 边属性设置
              E(g)$width <- abs(E(g)$Correlation) * 0.2  # 边宽与相关性强度正相关
              E(g)$color <- ifelse(E(g)$Correlation > 0, "#E64B35FF", "#3C5488FF")  # 正相关红，负相关蓝
              
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
                      size = abs(Log2FC)),   # 大小映射到Log2FC绝对值
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
                  values = c("Positive" = "#E64B35FF", "Negative" = "#3C5488FF")
                ) +
                # 边宽度映射（连续型）
                scale_edge_width_continuous(
                  name = "|Correlation|", 
                  range = c(0.5, 4),  # 边宽度范围
                  breaks = c(0.1, 0.4, 0.7)  # 图例刻度
                ) +
                # 节点颜色映射（离散型）
                scale_color_manual(
                  name = "Regulation",
                  values = c("up" = "#E64B35FF", "down" = "#3C5488FF")
                ) +
                # 节点大小映射（连续型）
                scale_size_continuous(
                  name = "|Log2FC|", 
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
              
              ggsave(file.path(output_folder, "Network_CytoscapeStyle_name.pdf"), network_plot, width = 12, height = 10)
              
              ggsave(file.path(output_folder, "Network_CytoscapeStyle_name.png"), network_plot, width = 12, height = 10, dpi = 300)
              
              # 导出网络数据供Cytoscape进一步编辑
              write.graph(g, file.path(output_folder, "feature_network_name.graphml"), format = "graphml")
              
              # 正确生成节点属性表的方法
              node_table <- diff_features %>%
                # 首先确保只选择实际存在的列
                dplyr::select(
                  truncated_name,  # 截断后的代谢物名称（必须存在于diff_features中）
                  Log2FC,         # 差异倍数
                  VIP,            # VIP值
                  regulation      # 上下调标记
                ) %>%
                # 添加原始名称（如果feature_names是独立向量）
                dplyr::mutate(
                  Original_Name = rownames(diff_features),  # 假设原始名称是行名
                  Direction = ifelse(regulation == "up", "Upregulated", "Downregulated"),
                  Metabolite = truncated_name  # 重命名列
                ) %>%
                # 重新排列列顺序
                dplyr::select(
                  Metabolite,
                  Original_Name,
                  Log2FC,
                  VIP,
                  Regulation = regulation,
                  Direction
                )
              
              # 生成边属性表（已筛选相关性强的边）
              edge_table <- cor_df %>%
                dplyr::mutate(
                  Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative"),
                  Abs_Correlation = abs(Correlation)  # 相关性绝对值
                ) %>%
                dplyr::select(
                  Source = Metabolite1,        # 源节点
                  Target = Metabolite2,        # 目标节点
                  Correlation,                 # 相关系数
                  Correlation_Type,            # 正/负相关
                  Abs_Correlation              # 相关性强度
                )
              
              # 计算网络拓扑指标
              node_topology <- data.frame(
                Metabolite = V(g)$name,
                Degree = igraph::degree(g),            # 连接度
                Betweenness = igraph::betweenness(g),  # 中介中心性
                Closeness = igraph::closeness(g)       # 接近中心性
              )
              
              # 合并节点属性与拓扑指标
              network_summary <- node_table %>%
                dplyr::left_join(node_topology, by = c("Metabolite" = "Metabolite"))
              
              # 创建新的Excel工作簿
              wb <- createWorkbook()
              
              # 1. 添加节点属性表（Sheet1）
              addWorksheet(wb, sheetName = "节点属性")
              writeData(wb, sheet = "节点属性", x = node_table, startCol = 1, startRow = 1)
              
              # 2. 添加边属性表（Sheet2）
              addWorksheet(wb, sheetName = "边属性")
              writeData(wb, sheet = "边属性", x = edge_table, startCol = 1, startRow = 1)
              
              # 3. 添加网络拓扑表（Sheet3，可选）
              addWorksheet(wb, sheetName = "网络拓扑")
              writeData(wb, sheet = "网络拓扑", x = network_summary, startCol = 1, startRow = 1)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder, "feature_network_name.xlsx"), overwrite = TRUE)
              
              ### 7.KEGG 富集 ---------------
              output_folder <- paste0("6.KEGG富集/", cmp_name)
              dir.create(paste0(output_folder,"/KEGG_pathway"), recursive = T, showWarnings = F)
              
              # 一个通用的小工具：拆分 + 去噪
              .split_ids <- function(x) {
                if (is.null(x) || all(is.na(x))) return(character(0))
                x <- ifelse(is.na(x), "", x)
                # 先按 ; 再按 /
                out <- strsplit(x, ";") |> unlist(use.names = FALSE)
                out <- strsplit(out, "/") |> unlist(use.names = FALSE)
                out <- trimws(out)
                out[!(is.na(out) | out == "" | out == "NA")]
              }
              
              # 把若干列（可能缺失）合并为去重 id_kegg 字符串
              .combine_kegg_cols <- function(df, cols) {
                if (length(cols) == 0) {
                  # 没有任何可用来源，直接返回原 df 的 id_kegg 或 NA
                  return(df$id_kegg)
                }
                purrr::pmap_chr(df[, cols, drop = FALSE], function(...) {
                  vals <- list(...)
                  ids <- unlist(lapply(vals, .split_ids), use.names = FALSE)
                  ids <- unique(ids)
                  if (length(ids) == 0) NA_character_ else paste(ids, collapse = ";")
                })
              }
              
              if (is.na(input$key_column_kegg)) {
                # 没有原始 KEGG 列：尝试把（名称匹配得到的）id_kegg + hmdb_kegg + lmsd_kegg 合并去重
                # 针对差异特征
                cols_df  <- intersect(c("id_kegg", "hmdb_kegg", "lmsd_kegg"),
                                      names(results_mode_diff_name_merged))
                diff_features <- results_mode_diff_name_merged %>%
                  dplyr::mutate(id_kegg = .combine_kegg_cols(cur_data_all(), cols_df))
                
                # 针对全部特征
                cols_all <- intersect(c("id_kegg", "hmdb_kegg", "lmsd_kegg"),
                                      names(results_mode_merged))
                all_features <- results_mode_merged %>%
                  dplyr::mutate(id_kegg = .combine_kegg_cols(cur_data_all(), cols_all))
                # 展开到逐 ID 形式并去 NA
                all_features <- all_features %>%
                  tidyr::separate_rows(id_kegg, sep = ";") %>%
                  dplyr::filter(!is.na(id_kegg), id_kegg != "")
                
              }else{
                
                diff_features <- results_mode_diff_name_merged
                all_features <- results_mode_merged%>%separate_rows(id_kegg)
                
              }
              
              tryCatch({
                
                # 提取差异代谢物KEGG信息
                kegg_data <- diff_features %>% separate_rows(id_kegg)%>%
                  dplyr::filter(!is.na(name)) %>%
                  dplyr::distinct(id_kegg, .keep_all = TRUE) %>%
                  dplyr::filter(!is.na(id_kegg))%>%.[which(.[["id_kegg"]]!="NA"),]
                
                mb3 <- enrichMBKEGG(kegg_data$id_kegg)
                enrich_result <- mb3@result%>%dplyr::arrange(desc(FoldEnrichment))
                enrich_result$Description <- factor(enrich_result$Description,levels = enrich_result$Description)
                enrich_result <- dplyr::rename(enrich_result, MetaboRatio=GeneRatio, metaboID=geneID)
                
                ## kegg通路制作
                # 从差异分析结果提取代谢物ID和Log2FC
                kegg_fc_data <- diff_features %>% 
                  dplyr::select(id_kegg, Log2FC, difference) %>%
                  dplyr::filter(!is.na(id_kegg)) %>%
                  dplyr::distinct(id_kegg, .keep_all = TRUE)%>%
                  tidyr::separate_rows(id_kegg)%>%.[which(.[["id_kegg"]]!="NA"),]
                
                setwd(paste0(output_folder, "/KEGG_pathway"))
                
                # 更新颜色
                update_node_colors <- function(cpd_data, metabo_data) {
                  
                  cpd_data %>%
                    mutate(
                      # 1. 分割多值ID
                      mapped_ids = ifelse(all.mapped == "", 
                                          list(character(0)), 
                                          str_split(all.mapped, ",")),
                      
                      # 2. 获取每个节点的状态向量
                      status_vec = map(mapped_ids, ~ {
                        if (length(.x) == 0) return(character(0))
                        # 关键修正：匹配difference而非Log2FC
                        matched_status <- metabo_data$difference[metabo_data$id_kegg %in% .x]
                        if (length(matched_status) > 0) matched_status else character(0)
                      }),
                      
                      # 3. 应用着色规则
                      mol.col = map_chr(status_vec, function(states) {
                        if (length(states) == 0) return("#FFFFFF")  # 无匹配保持白色
                        
                        # 核心着色逻辑（基于difference状态）
                        if (any(states == "up") & !any(states == "down")) {
                          "#FF0000"  # 红：仅含up
                        } else if (any(states == "down") & !any(states == "up")) {
                          "#0000FF"  # 蓝：仅含down
                        } else if (any(states == "up") & any(states == "down")) {
                          "#FFA500"  # 橙：混合状态
                        }
                      })
                    ) %>%
                    select(-mapped_ids, -status_vec)
                  
                }
                
                enrich_result$match <- "FALSE"
                
                # 目标匹配数量
                target_n <- nrow(enrich_result)
                matched_n <- 0
                
                # 遍历所有富集到的通路ID
                for (pathway_id in enrich_result$ID) {
                  
                  if (matched_n >= target_n) break
                  
                  # 提取当前通路的代谢物ID和Log2FC
                  pathway_metabo <- kegg_fc_data %>% dplyr::filter(id_kegg %in% (str_split(enrich_result$metaboID[which(enrich_result$ID==pathway_id)], "/", simplify = T)%>%as.character()))
                  
                  tryCatch({
                    
                    # 生成通路图（自动下载通路PNG和XML文件）
                    pathview.data <- pathview(
                      cpd.data = setNames(pathway_metabo$Log2FC, pathway_metabo$id_kegg),
                      pathway.id = gsub("map", "", pathway_id), # 去除前缀
                      species = input$species, # 物种缩写（如人类hsa）
                      kegg.native = T, # 生成KEGG原生风格图
                      same.layer = T, 
                      new.signature = F,
                      discrete = list(gene=FALSE, cpd=TRUE),
                      limit = list(gene=1, cpd=1),
                      bins = list(gene=2, cpd=2),   # 颜色分阶数
                      low = "blue", mid = "orange", high = "red", # 颜色梯度
                      out.suffix = "metabo_plot" # 输出文件后缀
                    )
                    
                    pathview.data$plot.data.cpd$mol.col[match(all_features$id_kegg, pathview.data$plot.data.cpd$kegg.names)%>%na.omit%>%as.numeric()] <- "#D3D3D3"
                      pathview.data$plot.data.cpd <- update_node_colors(pathview.data$plot.data.cpd, pathway_metabo)
                      
                      keggview.native(plot.data.cpd = pathview.data$plot.data.cpd[,-10], 
                                      pathway.name = gsub("map", input$species, pathway_id), 
                                      cols.ts.cpd = pathview.data$plot.data.cpd[,10], 
                                      limit = list(gene=1, cpd=1), 
                                      bins = list(gene=2, cpd=2), 
                                      low = list(gene = "green", cpd = "blue"), 
                                      mid = list(gene = "gray", cpd = "orange"), 
                                      high = list(gene = "red", cpd = "red"),
                                      discrete = list(gene=T, cpd=T),
                                      out.suffix = "metabo_plot")
                      
                  },error=function(e){})
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".metabo_plot.png"))) {
                    
                    if (file.info(paste0(gsub("map", input$species, pathway_id), ".metabo_plot.png"))%>%.$size>0) {
                      
                      enrich_result$match[which(enrich_result$ID == pathway_id)] <- "TRUE"
                      
                      matched_n <- matched_n + 1
                      
                    }
                    
                  }
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".png"))) {
                    
                    file.remove(paste0(gsub("map", input$species, pathway_id), ".png"))
                    
                  }
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".xml"))) {
                    
                    file.remove(paste0(gsub("map", input$species, pathway_id), ".xml"))
                    
                  }
                  
                }
                
                setwd(report_path)
                
                enrich_result <- enrich_result[which(enrich_result$match=="TRUE"), ]
                
                ## 气泡图
                # 创建标签截断函数
                truncate_labels <- function(x) {
                  ifelse(nchar(x) > 30, paste0(substr(x, 1, 30), "..."), x)
                }
                
                if (nrow(enrich_result)>20) {
                  enrich_result2 <- enrich_result[1:20,]
                }else{
                  enrich_result2 <- enrich_result
                }
                
                p_dot <- ggplot(enrich_result2, aes(x = FoldEnrichment,y = forcats::fct_rev(Description))) +  # 排序
                  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
                  scale_color_gradient(
                    low = "#56B1F7",  # 深蓝
                    high = "#132B43",  # 亮蓝
                    space = "Lab"  # 使用Lab色彩空间增强渐变平滑度[6](@ref)
                  )+
                  scale_size(name = "Counts", limits = c(0, min(max(enrich_result2$Count),20))) +
                  scale_y_discrete(
                    labels = truncate_labels,  # 应用截断函数[7,8](@ref)
                    expand = expansion(add = 0.2)  # 增加标签间距避免截断重叠
                  )+
                  labs(x = "Fold Enrichment", y = "", title = "KEGG Pathway Enrichment") +  # 标题注明排序依据
                  theme_bw() +
                  theme(
                    axis.text.y = element_text(
                      size = 12, 
                      color = "black",
                      hjust = 1,        # 右对齐核心参数[7,8](@ref)
                      margin = margin(r = 10)  # 右侧留出10单位边距[1](@ref)
                    ), plot.margin = margin(l = 3.5, r = 1, unit = "cm")  # 增大左侧绘图边距
                  )
                
                ggsave(file.path(output_folder, "KEGG_dotplot.pdf"), p_dot, width = 10, height = 8)
                ggsave(file.path(output_folder, "KEGG_dotplot.png"), p_dot, width = 10, height = 8)
                
                # 表格导出
                enrich_long <- enrich_result %>%
                  dplyr::mutate(ID_split = str_split(metaboID, "/")) %>%  # 拆分斜杠分隔的ID
                  unnest(ID_split)
                
                enrich_with_names <- enrich_long %>%
                  dplyr::left_join(
                    kegg_data %>% dplyr::distinct(id_kegg, name),  # 去重以避免重复合并
                    by = c("ID_split" = "id_kegg")          # 将enrich的ID_split与kegg的id_kegg匹配
                  )
                
                # 聚合名称到原始行格式
                enrich_final <- enrich_with_names %>%
                  dplyr::group_by(ID) %>%                            # 按富集条目分组
                  dplyr::summarise(
                    name = paste(unique(na.omit(name)), collapse = "/")  # 去重后合并名称
                  ) %>%
                  dplyr::right_join(enrich_result, by = "ID") %>%    # 合并回原始数据框
                  dplyr::relocate(name, .after = "metaboID") 
                
                write.xlsx(enrich_final, file.path(output_folder, "KEGG_enrichment.xlsx"))
                
              },error=function(e){
                
                setwd(report_path)
                
              })
              
              setwd(report_path)
              
              ### 8.MSEA -----------------
              output_folder <- paste0("7.MSEA/", cmp_name)
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 过滤无效ID（如空值或NA）
              msea_input <- separate_rows(results_mode_merged, "id_kegg")%>%as.data.frame()%>%.[!is.na(.$id_kegg) & .$id_kegg != "" & .$id_kegg != "NA", ]%>%dplyr::distinct(id_kegg,.keep_all = T)
              
              tryCatch({
                
                if (nrow(msea_input)>0) {
                  
                  if(exists("mSet")) rm(mSet)
                  
                  mSet <- InitDataObjects("conc", "msetora", FALSE, default.dpi = 72)
                  
                  system(paste0("cp -rf ", path_prefix, "/db/*.qs ", "./"))
                  
                  # Create vector consisting of compounds 
                  tmp.vec <- msea_input$id_kegg
                  
                  mSet<-Setup.MapData(mSet, tmp.vec)
                  
                  # Cross referencing
                  mSet<-CrossReferencing(mSet, "kegg")
                  
                  # Create the mapping results table
                  mSet<-CreateMappingResultTable(mSet)
                  
                  # Create list of candidates to replace the compound
                  mSet <- GetCandidateList(mSet)
                  
                  # Set the metabolite filter
                  mSet<-SetMetabolomeFilter(mSet, F);
                  
                  # Select metabolite set library, refer to 
                  mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0);
                  
                  # Calculate hypergeometric score, results table generated in your working directory
                  mSet<-CalculateHyperScore(mSet)
                  
                  # Plot the ORA, bar-graph
                  mSet<-PlotORA(mSet, paste0(output_folder, "/MSEA_"), "bar", "png", 300, width=NA)
                  
                  # export report
                  msea_table <- mSet$analSet$ora.mat
                  
                  write.xlsx(msea_table%>%as.data.frame%>%rownames_to_column("Metabolite pathway"), paste0(output_folder, "/MSEA_result.xlsx"))
                  
                  # 绘制enrichment score
                  # 自定义函数生成ES曲线
                  plotMetabEnrichment <- function(mSet, pathway, output_folder){
                    # 获取排序后的代谢物列表（基于logFC）
                    ranked_metabs <- msea_input[order(-msea_input$FC), ]
                    
                    # 获取当前通路代谢物集合
                    pathway_metabs <- mSet[["dataSet"]][["map.table"]][,"Query"][match(mSet$analSet$ora.hits[[pathway]], mSet[["dataSet"]][["map.table"]][,"Match"])]
                    
                    # 计算累计ES值
                    es_profile <- nrow(ranked_metabs)
                    hit_indices <- which((ranked_metabs$id_kegg) %in% pathway_metabs)
                    nh <- length(hit_indices)
                    nm <- nrow(ranked_metabs)
                    
                    if(nh > 0){
                      
                      phit <- rep(0, nm)
                      phit[hit_indices] <- abs(ranked_metabs$FC[hit_indices])^1 # 权重指数可调
                      phit <- cumsum(phit/sum(phit))
                      
                      pmiss <- cumsum(!(1:nm %in% hit_indices))
                      pmiss <- pmiss/(nm - nh)
                      
                      es_profile <- phit - pmiss
                      
                    }
                    
                    # 绘制ES曲线
                    df <- data.frame(
                      Position = 1:nm,
                      ES = es_profile,
                      Type = ifelse(1:nm %in% hit_indices, "Hit", "Background")
                    )
                    
                    p <- ggplot(df, aes(x=Position)) +
                      geom_line(aes(y=ES), color="#2c7bb6", size=1) +
                      geom_rug(aes(color=Type), sides="b") +
                      geom_hline(yintercept=0, linetype="dashed") +
                      labs(title=paste("Enrichment Profile:", pathway),
                           y="Enrichment Score") +
                      theme_minimal()
                    
                    ggsave(paste0(output_folder, "/", pathway, "_ES_plot.png"), p, width=8, height=6)
                    
                  }
                  
                  # 生成所有显著通路的ES图
                  sig_pathways <- rownames(mSet$analSet$ora.mat)[mSet$analSet$ora.mat%>%as.data.frame()%>%.$`Raw p` < 0.05]
                  
                  lapply(sig_pathways, function(pw){
                    
                    plotMetabEnrichment(mSet, pw, output_folder)
                    
                  })
                  
                  # 生成综合GSEA风格热图（参考clusterProfiler的gseaplot2）
                  metab_ranks <- msea_input$FC
                  names(metab_ranks) <- msea_input$id_kegg
                  metab_ranks <- sort(metab_ranks, decreasing=TRUE)
                  
                  top_pathways <- head(sig_pathways, 10)
                  plotlist <- lapply(top_pathways, function(pw){
                    fgsea::plotEnrichment(
                      pathway = mSet[["dataSet"]][["map.table"]][,"Query"][match(mSet$analSet$ora.hits[[pw]], mSet[["dataSet"]][["map.table"]][,"Match"])],
                      stats = metab_ranks
                    ) + labs(title=pw)
                  })
                  
                  pdf(paste0(output_folder, "/Top10Pathways_ES.pdf"), width=11, height=8)
                  gridExtra::grid.arrange(grobs = plotlist[1:min(6, length(plotlist))], nrow=3, ncol=2)
                  dev.off()
                  
                }
                
              },error=function(e){
                
                
                
              })
              
              ### 9.数据库分类（差异代谢物-有 name） ----------------
              output_folder <- paste0("8.数据库分类/", cmp_name)
              dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
              
              ## 小工具：清洗名称（去两端空格、转小写、去掉末尾括号内容）
              .clean_nm <- function(x) {
                if (is.null(x)) return(NA_character_)
                x <- tolower(trimws(x))
                x <- sub("\\s*\\([^()]*\\)$", "", x)  # 去掉末尾括号
                trimws(x)
              }
              
              ## 1) 构建 HMDB / LMSD 的“名称 → 分类”映射表
              # HMDB: name + synonyms2 → class / sub_class
              hmdb_name_map <- hmdb %>%
                dplyr::select(name, synonyms2, class, sub_class) %>%
                tidyr::pivot_longer(cols = c(name, synonyms2), names_to = "src", values_to = "nm") %>%
                tidyr::separate_rows(nm, sep = ";\\s*") %>%
                dplyr::mutate(nm = .clean_nm(nm)) %>%
                dplyr::filter(!is.na(nm), nm != "") %>%
                dplyr::distinct(nm, class, sub_class)
              
              # LMSD: NAME + SYSTEMATIC_NAME → CATEGORY / MAIN_CLASS / SUB_CLASS
              lmsd_name_map <- lmsd %>%
                dplyr::select(NAME, SYSTEMATIC_NAME, CATEGORY, MAIN_CLASS, SUB_CLASS) %>%
                tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
                dplyr::mutate(nm = .clean_nm(nm)) %>%
                dplyr::filter(!is.na(nm), nm != "") %>%
                dplyr::distinct(nm, CATEGORY, MAIN_CLASS, SUB_CLASS)
              
              ## 2) 取差异代谢物（必须有 name），拆分多名称并清洗
              diff_with_name <- results_mode_diff_name_merged %>%
                dplyr::filter(!is.na(name), trimws(name) != "") %>%
                tidyr::separate_rows(name, sep = ";") %>%
                dplyr::mutate(nm = .clean_nm(name)) %>%
                dplyr::filter(!is.na(nm), nm != "")
              
              ## 3) 关联到 HMDB / LMSD 分类
              diff_anno <- diff_with_name %>%
                dplyr::left_join(hmdb_name_map, by = "nm") %>%
                dplyr::left_join(lmsd_name_map, by = "nm") %>%
                dplyr::relocate(nm, .before = 1)
              
              # 导出明细
              readr::write_csv(
                diff_anno,
                file.path(output_folder, sprintf("差异代谢物_数据库分类明细_%s.csv", cmp_name))
              )
              
              ## 4) 通用画饼工具
              plot_cat_and_save <- function(df, col, title, out_prefix) {
                if (!col %in% names(df)) return(invisible(NULL))
                stat <- df %>%
                  dplyr::filter(!is.na(.data[[col]]), trimws(.data[[col]]) != "") %>%
                  dplyr::count(.data[[col]], name = "n") %>%
                  dplyr::arrange(dplyr::desc(n)) %>%
                  dplyr::mutate(pct = n / sum(n))
                
                if (nrow(stat) == 0) return(invisible(NULL))
                
                out_dir <- output_folder
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                
                if (nrow(stat) <= 10) {
                  ## 先把非有限/非正的 pct 去掉，避免 geom_col() 报错
                  stat <- stat %>%
                    dplyr::mutate(pct = as.numeric(pct)) %>%
                    dplyr::filter(is.finite(pct), pct > 0)
                  
                  # 计算角度、侧别、坐标（在你的基础上微调为更安全的半径/范围）
                  stat <- stat %>%
                    dplyr::mutate(
                      ymax  = cumsum(pct),
                      ymin  = dplyr::lag(ymax, default = 0),
                      ymid  = (ymax + ymin) / 2,                 # 扇区中点(0-1)
                      angle = ymid * 360 - 90,
                      side  = dplyr::if_else(angle > 90 & angle < 270, "left", "right"),
                      label = paste0(.data[[col]], " (", scales::percent(pct), ")"),
                      # 径向段终点（略在饼图外）
                      x_rad = dplyr::if_else(side == "left", 0.93, 1.07),
                      # 肘形水平段终点（更外侧，便于排布标签）
                      x_elb = dplyr::if_else(side == "left", 0.75, 1.25),
                      # 标签放在肘形终点稍外一点
                      x_lab = dplyr::if_else(side == "left", 0.72, 1.28),
                      hjust = dplyr::if_else(side == "left", 1, 0)
                    )
                  
                  # 扇区（关键：limits 覆盖 [0.5,1.5] 并再留余量；width < 1 防边界擦边）
                  p <- ggplot(stat, aes(x = 1, y = pct, fill = .data[[col]])) +
                    geom_col(width = 0.92, color = "white", na.rm = TRUE) +  # na.rm 进一步静默无效行
                    coord_polar(theta = "y", clip = "off") +
                    theme_void() +
                    labs(title = title, fill = NULL) +
                    # 覆盖 [0.45,1.55]，确保整条环（中心1、厚度~0.92）在范围内，并给肘线/标签留空
                    scale_x_continuous(limits = c(0.45, 1.55), expand = expansion(mult = c(0, 0))) +
                    theme(
                      legend.position = "none",
                      plot.margin = margin(10, 110, 10, 90)   # 两侧更大的外边距，防截断
                    )
                  
                  # 径向引线（从 1.00 到 x_rad）
                  p <- p +
                    geom_segment(
                      data = dplyr::filter(stat, side == "right"),
                      aes(x = 1.00, xend = x_rad, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    ) +
                    geom_segment(
                      data = dplyr::filter(stat, side == "left"),
                      aes(x = 1.00, xend = x_rad, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    )
                  
                  # 肘形水平引线（x_rad -> x_elb）
                  p <- p +
                    geom_segment(
                      data = dplyr::filter(stat, side == "right"),
                      aes(x = x_rad, xend = x_elb, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    ) +
                    geom_segment(
                      data = dplyr::filter(stat, side == "left"),
                      aes(x = x_rad, xend = x_elb, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    )
                  
                  # 标签（优先 ggrepel，沿 y 方向自动避让；放在外圈 x_lab）
                  if (requireNamespace("ggrepel", quietly = TRUE)) {
                    p <- p +
                      ggrepel::geom_text_repel(
                        data = dplyr::filter(stat, side == "right"),
                        aes(x = x_lab, y = ymid, label = label),
                        inherit.aes = FALSE,
                        size = 3, direction = "y",
                        nudge_x = 0, box.padding = 0.3, point.padding = 0.25,
                        min.segment.length = 0, segment.size = 0.3, max.overlaps = Inf
                      ) +
                      ggrepel::geom_text_repel(
                        data = dplyr::filter(stat, side == "left"),
                        aes(x = x_lab, y = ymid, label = label),
                        inherit.aes = FALSE,
                        size = 3, direction = "y",
                        nudge_x = 0, box.padding = 0.3, point.padding = 0.25,
                        min.segment.length = 0, segment.size = 0.3, max.overlaps = Inf
                      )
                  } else {
                    p <- p +
                      geom_text(
                        data = dplyr::filter(stat, side == "right"),
                        aes(x = x_lab, y = ymid, label = label, hjust = 0),
                        inherit.aes = FALSE, size = 3
                      ) +
                      geom_text(
                        data = dplyr::filter(stat, side == "left"),
                        aes(x = x_lab, y = ymid, label = label, hjust = 1),
                        inherit.aes = FALSE, size = 3
                      )
                  }
                  
                } else {
                  stat <- stat %>% dplyr::mutate(label = scales::percent(pct))
                  max_n <- max(stat$n, na.rm = TRUE)
                  pad   <- max(1, max_n * 0.35)
                  
                  p <- ggplot(stat, aes(x = reorder(.data[[col]], n), y = n, fill = .data[[col]])) +
                    geom_col() +
                    coord_flip(clip = "off") +
                    scale_y_continuous(limits = c(0, max_n + pad),
                                       expand = expansion(mult = c(0, 0.05))) +
                    labs(title = paste0(title), x = NULL, y = "Count") +
                    theme_minimal() +
                    theme(
                      legend.position = "none",
                      plot.margin = margin(10, 90, 10, 10)
                    ) +
                    geom_text(aes(label = label), hjust = -0.05, size = 3)
                }
                
                ggsave(file.path(out_dir, paste0(out_prefix, ".png")), p, width = 8, height = 6, dpi = 300)
                ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), p, width = 8, height = 6)
                
                readr::write_csv(stat, file.path(out_dir, paste0(out_prefix, "_统计表.csv")))
              }
              
              ## 5) 分别绘制 5 类饼图
              # HMDB: class
              plot_cat_and_save(
                diff_anno, "class",
                title = sprintf("HMDB class - %s", cmp_name),
                out_prefix = sprintf("HMDB_class_%s", cmp_name)
              )
              
              # HMDB: sub_class
              plot_cat_and_save(
                diff_anno, "sub_class",
                title = sprintf("HMDB sub_class - %s", cmp_name),
                out_prefix = sprintf("HMDB_sub_class_%s", cmp_name)
              )
              
              # LMSD: CATEGORY
              plot_cat_and_save(
                diff_anno, "CATEGORY",
                title = sprintf("LMSD CATEGORY - %s", cmp_name),
                out_prefix = sprintf("LMSD_CATEGORY_%s", cmp_name)
              )
              
              # LMSD: MAIN_CLASS
              plot_cat_and_save(
                diff_anno, "MAIN_CLASS",
                title = sprintf("LMSD MAIN_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_MAIN_CLASS_%s", cmp_name)
              )
              
              # LMSD: SUB_CLASS
              plot_cat_and_save(
                diff_anno, "SUB_CLASS",
                title = sprintf("LMSD SUB_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_SUB_CLASS_%s", cmp_name)
              )
              
              ### 10.数据库分类（总代谢物-有 name） ----------------
              output_folder <- paste0("8.数据库分类/", cmp_name)
              dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
              
              ## 2) 取差异代谢物（必须有 name），拆分多名称并清洗
              total_with_name <- results_mode_merged %>%
                dplyr::filter(!is.na(name), trimws(name) != "") %>%
                tidyr::separate_rows(name, sep = ";") %>%
                dplyr::mutate(nm = .clean_nm(name)) %>%
                dplyr::filter(!is.na(nm), nm != "")
              
              ## 3) 关联到 HMDB / LMSD 分类
              total_anno <- total_with_name %>%
                dplyr::left_join(hmdb_name_map, by = "nm") %>%
                dplyr::left_join(lmsd_name_map, by = "nm") %>%
                dplyr::relocate(nm, .before = 1)
              
              # 导出明细
              readr::write_csv(
                total_anno,
                file.path(output_folder, sprintf("总代谢物_数据库分类明细_%s.csv", cmp_name))
              )
              
              ## 4) 分别绘制 5 类饼图
              # HMDB: class
              plot_cat_and_save(
                total_anno, "class",
                title = sprintf("HMDB class - %s", cmp_name),
                out_prefix = sprintf("HMDB_class_all_%s", cmp_name)
              )
              
              # HMDB: sub_class
              plot_cat_and_save(
                total_anno, "sub_class",
                title = sprintf("HMDB sub_class - %s", cmp_name),
                out_prefix = sprintf("HMDB_sub_class_all_%s", cmp_name)
              )
              
              # LMSD: CATEGORY
              plot_cat_and_save(
                total_anno, "CATEGORY",
                title = sprintf("LMSD CATEGORY - %s", cmp_name),
                out_prefix = sprintf("LMSD_CATEGORY_all_%s", cmp_name)
              )
              
              # LMSD: MAIN_CLASS
              plot_cat_and_save(
                total_anno, "MAIN_CLASS",
                title = sprintf("LMSD MAIN_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_MAIN_CLASS_all_%s", cmp_name)
              )
              
              # LMSD: SUB_CLASS
              plot_cat_and_save(
                total_anno, "SUB_CLASS",
                title = sprintf("LMSD SUB_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_SUB_CLASS_all_%s", cmp_name)
              )
              
            }
            
            cat("已保存", cmp_name, "的结果至文件。\n")
            
            # 统计各阶段的feature数量及差异筛选的上下调数量
            diff_count <- nrow(results_mode_diff_merged)
            diff_up <- sum(results_mode_diff_merged$difference == "up", na.rm = TRUE)
            diff_down <- sum(results_mode_diff_merged$difference == "down", na.rm = TRUE)
            diff_name_count <- nrow(results_mode_diff_name_merged)
            diff_name_up <- sum(results_mode_diff_name_merged$difference == "up", na.rm = TRUE)
            diff_name_down <- sum(results_mode_diff_name_merged$difference == "down", na.rm = TRUE)
            
            # 将当前比较的统计信息保存到列表中
            if (!input$single_ion_judge) {
              
              summary_stats_list[[cmp_name]][["Merge"]] <- data.frame(
                compare = cmp_name,
                ion_type = "POS&NEG",
                raw_feature_num = summary_stats_list[[cmp_name]][["POS"]]$raw_feature_num+summary_stats_list[[cmp_name]][["NEG"]]$raw_feature_num,
                remove_missing_num = summary_stats_list[[cmp_name]][["POS"]]$remove_missing_num+summary_stats_list[[cmp_name]][["NEG"]]$remove_missing_num,
                preprocess_feature_num = nrow(results_mode_merged),
                preprocess_feature_num_name = sum(!is.na(results_mode_merged$name)),
                diff_feature = diff_count,
                diff_feature_up = diff_up,
                diff_feature_down = diff_down,
                diff_feature_name = diff_name_count,
                diff_feature_name_up = diff_name_up,
                diff_feature_name_down = diff_name_down,
                stringsAsFactors = FALSE
              )
              
            }
            
          }
          
          # 合并所有分组比较的统计信息
          final_summary <- do.call(rbind, do.call(rbind, summary_stats_list))
          
          # 保存结果rda数据
          save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw", "final_summary"), file = "result.rda")
          
          # 保存最终统计结果到 "差异筛选/差异数量统计.xlsx"
          write.xlsx(final_summary, file = "5.差异代谢物分析/差异数量统计.xlsx", rowNames = FALSE)
          cat("已保存差异数量统计结果至文件：5.差异代谢物分析/差异数量统计.xlsx\n")
          
          # 组间差异结果 ---------------
          tryCatch({
            
            ## feature --------------------
            diff_features_list <- list()
            
            # 遍历所有比较组
            for(cmp_name in group.data) {
              
              diff_features_list[[cmp_name]] <- do.call(rbind, results_mode_diff[[cmp_name]])%>%.[,c("peak_name")]
              
            }
            
            ### Venn for names ---------------------
            if(length(group.data) >= 2 & length(group.data) <=5) {
              
              output_folder <- paste0("5.差异代谢物分析/2.组间对比/1.venn")
              dir.create(output_folder, recursive = T)
              
              # 动态选择调色板（支持2-5组）
              if(length(diff_features_list) == 2) {
                
                fill_colors <- brewer.pal(3, "Paired")[1:2]  # 使用Paired调色板前两色
                
              } else {
                
                fill_colors <- brewer.pal(length(diff_features_list), "Set2")  # Set2支持3-8色
                
              }
              
              venn.plot <- venn.diagram(
                x = diff_features_list,
                filename = file.path(output_folder, "Venn_diagram.png"),
                disable.logging = T,
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
              
              # 生成Venn交集信息（修正列名问题）
              venn_data <- get.venn.partitions(diff_features_list) %>%
                rename(
                  Set_Combination = ..set..,  # 根据网页1的列名修正
                  Peak_Names = ..values..,
                  Count = ..count..
                )
              
              # 合并所有差异代谢物完整信息
              full_data <- do.call(rbind, lapply(names(results_mode_diff), function(cmp_name){
                cbind(Comparison = cmp_name, 
                      do.call(rbind, results_mode_diff[[cmp_name]]) %>% 
                        # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb))
                        dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg))
              })) %>% dplyr::distinct(peak_name, .keep_all = TRUE)
              
              # 初始化Excel工作簿
              wb <- createWorkbook()
              
              # 初始化一个空数据框存储所有数据
              combined_data <- data.frame()
              
              # 遍历每个Venn交集区域
              for(i in 1:nrow(venn_data)) {
                # 获取当前交集区域的peak_name列表
                current_peaks <- unlist(venn_data$Peak_Names[i])
                
                # 提取代谢物信息并添加标识列
                sheet_data <- full_data %>%
                  dplyr::filter(peak_name %in% current_peaks) %>%
                  # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb) %>%
                  dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg) %>%
                  dplyr::mutate(Set_Combination = venn_data$Set_Combination[i])  # 关键：添加标识列[2,6](@ref)
                
                # 合并到总数据框
                combined_data <- dplyr::bind_rows(combined_data, sheet_data)  # 使用dplyr的合并方法[6](@ref)
              }
              
              # 创建单个工作表并写入数据
              addWorksheet(wb, sheetName = "Combined_Venn_Data")
              writeData(wb, sheet = "Combined_Venn_Data", x = combined_data)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder, "Venn_diagram.xlsx"), overwrite = TRUE)
              
            }
            
            ### Upset ------------------------
            if(length(group.data) >= 2) {
              
              output_folder <- paste0("5.差异代谢物分析/2.组间对比/2.upset")
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
            # 收集所有比较组的差异代谢物名称（基于name）
            diff_names_list <- list()
            
            # 遍历所有比较组
            for(cmp_name in group.data) {
              
              diff_names_list[[cmp_name]] <- do.call(rbind, results_mode_diff_name[[cmp_name]])%>%.[,c("peak_name")]
              
            }
            
            ### Venn for names ---------------------
            if(length(diff_names_list) >= 2 & length(diff_names_list) <= 5) {
              
              output_folder_name <- paste0("5.差异代谢物分析/2.组间对比/1.venn")
              
              # 动态选择调色板
              if(length(diff_names_list) == 2) {
                
                fill_colors <- brewer.pal(3, "Paired")[1:2]
                
              } else {
                
                fill_colors <- brewer.pal(length(diff_names_list), "Set2")
                
              }
              
              venn.plot <- venn.diagram(
                x = diff_names_list,
                filename = file.path(output_folder_name, "Venn_diagram_name.png"),
                disable.logging = T,
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
              
              # 生成Venn交集信息（修正列名问题）
              venn_data <- get.venn.partitions(diff_names_list) %>%
                rename(
                  Set_Combination = ..set..,  # 根据网页1的列名修正
                  Peak_Names = ..values..,
                  Count = ..count..
                )
              
              # 合并所有差异代谢物完整信息
              full_data <- do.call(rbind, lapply(names(results_mode_diff), function(cmp_name){
                cbind(Comparison = cmp_name, 
                      do.call(rbind, results_mode_diff[[cmp_name]]) %>% 
                        # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb))
                        dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg))
              })) %>% dplyr::distinct(peak_name, .keep_all = TRUE)
              
              # 初始化Excel工作簿
              wb <- createWorkbook()
              
              # 初始化一个空数据框存储所有数据
              combined_data <- data.frame()
              
              # 遍历每个Venn交集区域
              for(i in 1:nrow(venn_data)) {
                # 获取当前交集区域的peak_name列表
                current_peaks <- unlist(venn_data$Peak_Names[i])
                
                # 提取代谢物信息并添加标识列
                sheet_data <- full_data %>%
                  dplyr::filter(peak_name %in% current_peaks) %>%
                  # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb) %>%
                  dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg) %>%
                  dplyr::mutate(Set_Combination = venn_data$Set_Combination[i])  # 关键：添加标识列[2,6](@ref)
                
                # 合并到总数据框
                combined_data <- dplyr::bind_rows(combined_data, sheet_data)  # 使用dplyr的合并方法[6](@ref)
              }
              
              # 创建单个工作表并写入数据
              addWorksheet(wb, sheetName = "Combined_Venn_Data")
              writeData(wb, sheet = "Combined_Venn_Data", x = combined_data)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder_name, "Venn_diagram_name.xlsx"), overwrite = TRUE)
              
            }
            
            ### Upset for names ------------------------
            if(length(diff_names_list) >= 2){
              
              output_folder_name <- paste0("5.差异代谢物分析/2.组间对比/2.upset")
              
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
              
              pdf(file.path(output_folder_name, "UpSet_plot_name.pdf"), width = 20, height = 12)
              print(upset_plot_name)
              dev.off()
              
              png(file.path(output_folder_name, "UpSet_plot_name.png"), width = 800, height = 600)
              print(upset_plot_name)
              dev.off()
              
              write.xlsx(binary_matrix_name, file.path(output_folder_name, "UpSet_data_name.xlsx"))
              
            }
            
          },error=function(e){
            
            
            
          })
          
        }else{
          
          if (!input$single_ion_judge){
            
            final_summary <- data.frame(
              compare = "no compare",
              ion_type = "POS&NEG",
              raw_feature_num = summary_stats_list[["POS"]]$raw_feature_num+summary_stats_list[["NEG"]]$raw_feature_num,
              remove_missing_num = summary_stats_list[["POS"]]$remove_missing_num+summary_stats_list[["NEG"]]$remove_missing_num,
              preprocess_feature_num = summary_stats_list[["POS"]]$preprocess_feature_num+summary_stats_list[["NEG"]]$preprocess_feature_num,
              preprocess_feature_num_name = summary_stats_list[["POS"]]$preprocess_feature_num_name+summary_stats_list[["NEG"]]$preprocess_feature_num_name,
              stringsAsFactors = FALSE
            )
            
          }
          
          # 合并所有分组比较的统计信息
          final_summary <- do.call(rbind, do.call(rbind, summary_stats_list))
          
          # 保存结果rda数据
          save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw", "final_summary"), file = "result.rda")
          
        }
        
        # 出报告 -----------------
        para_sum <- list(miss=c("Within the group" = "inter.group",
                                "Global filter" = "global.group"),
                         fill=c("Intra-group mean" = "mean.group", 
                                "Global mean" = "mean.global", 
                                "Median within the group" = "median.group",
                                "Global median" = "median.global", 
                                "None" = "none",
                                "Min" = "min.global",
                                "Min/2" = "min2",
                                "KNN" = "knn.global",
                                "Randomforest" = "rf"),
                         normalize=c("Sum" = "sum",
                                     "None" = "none",
                                     "QC svr"= "qc",
                                     "QC RLSC"= "qc-rlsc",
                                     "QC NormAE"= "normae",
                                     "Probability quotient" = "prob_quot",
                                     "0.75 quantile" = "percent_0.75")
        )
        
        para_data <- data.frame("缺失值类型" = input$miss.value.handle.type, 
                                "缺失值过滤方法" = names(para_sum$miss)[which(para_sum$miss==input$miss.value.handle.group)], 
                                "缺失值过滤比率" = input$miss.value.handle.cutoff,
                                "填充方式" = names(para_sum$fill)[which(para_sum$fill==input$miss.value.fill)],
                                "归一化方法" = names(para_sum$normalize)[which(para_sum$normalize==input$normalized.handle.method)],
                                "SUM的系数" = input$sum_coef,
                                "定量值取log" = input$log.handle.method,
                                "RSD过滤阈值" = input$rsd.cutoff,
                                "RSD取log" = input$log.rsd.method,
                                "批次矫正" = input$batch.correct,
                                "相关性阈值" = input$correlation.cutoff,
                                "是否使用FC筛选" = input$fc_cutoff_judge,
                                "Fold Change Cutoff" = input$fc_cutoff,
                                "是否使用pvalue筛选" = input$p_cutoff_judge,
                                "检验类型" = input$pvalue_type,
                                "P-value Cutoff" = input$p_cutoff,
                                "是否使用p-adjust筛选" = input$padjust_judge,
                                "padjust方法" = input$padjust_method,
                                "是否使用VIP筛选" = input$vip_cutoff_judge,
                                "VIP Cutoff" = input$vip_cutoff,
                                "KEGG物种" = input$species,
                                "选择的列数" = paste0(selected_column, collapse = ";"))%>%t%>%data.frame()
        
        if (input$single_ion_judge) {
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_single.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_server_single.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report_server.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_single_word.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.docx'))
        }else{
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_server.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report_server.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_word.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.docx'))
        }
        
      }
      
      if(para[["info"]][["type"]] == 'omicsolution'){
        
        dir.create(report_path, recursive = T)
        
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
        
        if (input$single_ion_judge) {
          
          metabo.file <- paste0(data_path, "/metabonomics.csv")
          
          ion_type_all <- "single"
          
          names(metabo.file) <- ion_type_all
          
          sample.data <- list(single = rv$condition.data.single)
          
        } else {
          
          metabo.file <- c(paste0(data_path, "/metabonomics_NEG.csv"), paste0(data_path, "/metabonomics_POS.csv"))
          
          ion_type_all <- c("NEG","POS")
          
          names(metabo.file) <- ion_type_all
          
          sample.data <- list(POS = rv$condition.data.pos, NEG = rv$condition.data.neg)
          
        }
        
        if (input[['batch.correct']]!="none") {
          
          sample.data$batch <- rv$condition.data.pos[[grep("batch",colnames(rv$condition.data.pos),value = T,ignore.case = T)]]
          
        }
        
        group.data <- input$metaboQuant_compare_select
        
        if (!is.na(input$key_column_kegg))  if (input$key_column_kegg  == 0) input$key_column_kegg  <- NA
        if (!is.na(input$key_column_hmdb))  if (input$key_column_hmdb  == 0) input$key_column_hmdb  <- NA
        if (!is.na(input$key_column_lmsd))  if (input$key_column_lmsd  == 0) input$key_column_lmsd  <- NA
        if (!is.na(input$key_column_name))  if (input$key_column_name  == 0) input$key_column_name  <- NA
        
        if (!is.na(input$key_column_name)) {
          
          if (is.na(input$key_column_kegg)) {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_name,
              input$key_column_sample_start:input$key_column_sample_end
            )
          } else {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_name,
              input$key_column_kegg,
              input$key_column_sample_start:input$key_column_sample_end
            )
          }
          
        }else{
          
          if (is.na(input$key_column_kegg)) {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_sample_start:input$key_column_sample_end
            )
          } else {
            selected_column <- c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt,
              input$key_column_kegg,
              input$key_column_sample_start:input$key_column_sample_end
            )
          }
          
        }
        
        # 存储当前比较中POS和NEG的结果
        metabo.data.raw.all.col <- summary_stats_list <- results_mode_diff_name <- results_mode_diff <- results_mode <- metabo.data.raw <- list()
        
        sum_pca <- list(); sum_plsda <- list(); sum_oplsda <- list()
        
        # kegg 数据库
        load(file=paste0(path_prefix, '/db/kegg_id_db.rda'))
        
        if (is.na(input$key_column_kegg)) {
          
          ## 构建 name → KEGG ID 索引
          hmdb2 <- tidyr::separate_rows(hmdb, synonyms2, sep = "; ")
          hmdb2[["synonyms_tolower"]] <- tolower(hmdb2$synonyms2)
          hmdb2 <- dplyr::distinct(hmdb2, synonyms_tolower, kegg_id, .keep_all = T)
          kegg2 <- tidyr::separate_rows(kegg, name, sep = "; ")
          kegg2[["name_tolower"]] <- tolower(kegg2$name)
          kegg2 <- dplyr::distinct(kegg2, name_tolower, kegg_id, .keep_all = T)
          # 预处理：过滤掉kegg_id为NA的行
          hmdb2 <- hmdb2[!is.na(hmdb2$kegg_id), ]
          kegg2 <- kegg2[!is.na(kegg2$kegg_id), ]
          # 创建快速查询索引
          hmdb_index <- split(hmdb2$kegg_id, hmdb2$synonyms_tolower)
          kegg_index <- split(kegg2$kegg_id, kegg2$name_tolower)
          
          lmsd_name_tbl <- lmsd %>%
            dplyr::select(LM_ID, KEGG_ID, NAME, SYSTEMATIC_NAME) %>%
            tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
            dplyr::filter(!is.na(nm), nm != "") %>%
            dplyr::mutate(nm_tolower = tolower(nm)) %>%
            dplyr::filter(!is.na(KEGG_ID), KEGG_ID != "") %>%
            dplyr::distinct(nm_tolower, KEGG_ID)
          lmsd_index <- split(lmsd_name_tbl$KEGG_ID, lmsd_name_tbl$nm_tolower)
          
          ## 构建 HMDB accession → KEGG ID 索引
          hmdb3 <- dplyr::distinct(hmdb2, accession, kegg_id, .keep_all = T)
          accession_index <- split(hmdb3$kegg_id, hmdb3$accession)
          
          ## 构建 LM_ID → KEGG ID 索引
          lmsd_lm_tbl <- lmsd %>%
            dplyr::filter(!is.na(LM_ID), LM_ID != "", !is.na(KEGG_ID), KEGG_ID != "") %>%
            dplyr::distinct(LM_ID, KEGG_ID)
          lmsd_lm_index <- split(lmsd_lm_tbl$KEGG_ID, lmsd_lm_tbl$LM_ID)
          
        }
        
        # 数据预处理与各组的差异筛选 ----
        for (ion_type in ion_type_all) {
          
          # dir.create("1.Raw_data")
          dir.create("1.搜库原始结果")
          
          # system(paste0("cp -rf ", metabo.file[[ion_type]], " '1.Raw_data'"))
          system(paste0("cp -rf ", metabo.file[[ion_type]], " '1.搜库原始结果'"))
          
          metabo.data.raw.all.col[[ion_type]] <- read_csv(metabo.file[[ion_type]], locale = readr::locale(encoding = "UTF-8"))
          
          metabo.data.raw[[ion_type]] <- metabo.data.raw.all.col[[ion_type]][,selected_column]
          
          if (!is.na(input$key_column_name)) {
            
            if (is.na(input$key_column_kegg)) {
              
              colnames(metabo.data.raw[[ion_type]])[1:4] <- c("peak_name", "mz", "rt", "name")
              
              # 获取待匹配的名称向量
              query_names <- tolower(metabo.data.raw[[ion_type]][["name"]])
              
              matched_kegg <- lapply(query_names, function(name_str) {
                sub_names <- unlist(strsplit(name_str, ";"))
                sub_ids <- sapply(sub_names, function(name) {
                  name <- trimws(name)
                  
                  hmdb_ids <- if (!is.null(hmdb_index[[name]])) unname(hmdb_index[[name]]) else character(0)
                  kegg_ids <- if (!is.null(kegg_index[[name]])) unname(kegg_index[[name]]) else character(0)
                  lmsd_ids <- if (!is.null(lmsd_index[[name]])) unname(lmsd_index[[name]]) else character(0)
                  all_ids  <- unique(c(hmdb_ids, kegg_ids, lmsd_ids))
                  
                  # 去括号再试一次
                  if (length(all_ids) == 0 && grepl("\\([^()]*\\)$", name)) {
                    name2 <- sub("\\s*\\([^()]*\\)$", "", name)
                    hmdb2 <- if (!is.null(hmdb_index[[name2]])) unname(hmdb_index[[name2]]) else character(0)
                    kegg2 <- if (!is.null(kegg_index[[name2]])) unname(kegg_index[[name2]]) else character(0)
                    lmsd2 <- if (!is.null(lmsd_index[[name2]])) unname(lmsd_index[[name2]]) else character(0)
                    all_ids <- unique(c(hmdb2, kegg2, lmsd2))
                  }
                  
                  if (length(all_ids) > 0) paste(all_ids, collapse = "/") else "NA"
                }, USE.NAMES = FALSE)
                
                paste(sub_ids, collapse = ";")
              })
              
              # 4.2 HMDB accession → KEGG，生成 hmdb_kegg（若提供 HMDB 列）
              if (!is.na(input$key_column_hmdb)) {
                hmdb_kegg <- vapply(seq_len(nrow(metabo.data.raw.all.col[[ion_type]])), function(i) {
                  acc <- metabo.data.raw.all.col[[ion_type]][[input$key_column_hmdb]][i]
                  if (is.null(acc) || is.na(acc) || acc == "") return(NA_character_)
                  acc_list <- strsplit(toupper(as.character(acc)), ";")[[1]]
                  
                  ids <- sapply(acc_list, function(a) {
                    a <- trimws(a)
                    res <- accession_index[[a]]
                    if (!is.null(res) && length(res) > 0) paste(unique(res), collapse = "/") else NA_character_
                  }, USE.NAMES = FALSE)
                  
                  paste(ids, collapse = ";")
                }, character(1))
              }
              
              # 4.3 LM_ID → KEGG，生成 lmsd_kegg（若提供 LMSD 列）
              if (!is.na(input$key_column_lmsd)) {
                lm_col <- metabo.data.raw.all.col[[ion_type]][[ input$key_column_lmsd ]]
                lmsd_kegg <- vapply(seq_along(lm_col), function(i) {
                  val <- lm_col[i]
                  if (is.null(val) || is.na(val) || val == "") return(NA_character_)
                  lm_list <- strsplit(as.character(val), ";")[[1]]
                  ids <- sapply(lm_list, function(lm) {
                    key <- trimws(lm)
                    res <- lmsd_lm_index[[key]]
                    if (!is.null(res) && length(res) > 0) paste(unique(res), collapse = "/") else NA_character_
                  }, USE.NAMES = FALSE)
                  paste(ids, collapse = ";")
                }, character(1))
              }
              
              # 4.4 组装输出（按需附加）
              tmp <- metabo.data.raw[[ion_type]][, 1:4]
              tmp$id_kegg <- unlist(matched_kegg)
              
              if (!is.na(input$key_column_hmdb)) tmp$hmdb_kegg <- hmdb_kegg
              if (!is.na(input$key_column_lmsd)) tmp$lmsd_kegg <- lmsd_kegg
              
              metabo.data.raw[[ion_type]] <- data.frame(
                tmp,
                metabo.data.raw[[ion_type]][, -c(1:4)],
                check.names = FALSE
              )
              
            }else{
              
              colnames(metabo.data.raw[[ion_type]])[1:5] <- c("peak_name", "mz", "rt", "name", "id_kegg")
              
            }
            
          } else {
            
            # ========= 当没有 name 列时：先用 KEGG/HMDB/LMSD 的 ID 反查 name，再用 name 反查/补全 id_kegg =========
            message("未提供 name 列，基于 KEGG/HMDB/LMSD ID 自动回填 name，并补全 id_kegg...")
            
            # —— 1) 构建 ID -> name 的索引（不 reload，直接用上面已存在的 kegg/hmdb/lmsd 对象）——
            kegg_name_index <- kegg %>%
              dplyr::select(kegg_id, name) %>%
              dplyr::filter(!is.na(kegg_id), kegg_id != "", !is.na(name), name != "") %>%
              dplyr::distinct(kegg_id, .keep_all = TRUE) %>%
              tibble::deframe()
            
            hmdb_name_index <- hmdb %>%
              dplyr::select(accession, name) %>%
              dplyr::filter(!is.na(accession), accession != "", !is.na(name), name != "") %>%
              dplyr::distinct(accession, .keep_all = TRUE) %>%
              tibble::deframe()
            
            lmsd_name_index <- lmsd %>%
              dplyr::select(LM_ID, NAME) %>%
              dplyr::filter(!is.na(LM_ID), LM_ID != "", !is.na(NAME), NAME != "") %>%
              dplyr::distinct(LM_ID, .keep_all = TRUE) %>%
              tibble::deframe()
            
            # —— 2) 构建 name -> KEGG_ID 的索引（用于后续用 name 反推 id_kegg）——
            # 2.1 KEGG: name 可能是以分号分隔的多别名
            kegg2 <- tidyr::separate_rows(kegg, name, sep = "; ") %>%
              dplyr::mutate(name_tolower = tolower(name)) %>%
              dplyr::filter(!is.na(name_tolower), name_tolower != "", !is.na(kegg_id), kegg_id != "") %>%
              dplyr::distinct(name_tolower, kegg_id)
            
            kegg_index_by_name <- split(kegg2$kegg_id, kegg2$name_tolower)
            
            # 2.2 HMDB: 用 synonyms2 里的别名对 KEGG_ID 做映射（通过 hmdb$kegg_id）
            hmdb2 <- tidyr::separate_rows(hmdb, synonyms2, sep = "; ") %>%
              dplyr::mutate(synonyms_tolower = tolower(synonyms2)) %>%
              dplyr::filter(!is.na(synonyms_tolower), synonyms_tolower != "", !is.na(kegg_id), kegg_id != "") %>%
              dplyr::distinct(synonyms_tolower, kegg_id)
            
            hmdb_index_by_name <- split(hmdb2$kegg_id, hmdb2$synonyms_tolower)
            
            # 2.3 LMSD: NAME/Systematic name -> KEGG_ID
            lmsd_name_tbl <- lmsd %>%
              dplyr::select(LM_ID, KEGG_ID, NAME, SYSTEMATIC_NAME) %>%
              tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
              dplyr::filter(!is.na(nm), nm != "", !is.na(KEGG_ID), KEGG_ID != "") %>%
              dplyr::mutate(nm_tolower = tolower(nm)) %>%
              dplyr::distinct(nm_tolower, KEGG_ID)
            
            lmsd_index_by_name <- split(lmsd_name_tbl$KEGG_ID, lmsd_name_tbl$nm_tolower)
            
            # —— 3) 用 ID 反查得到 name —— 
            n_row <- nrow(metabo.data.raw.all.col[[ion_type]])
            
            name_from_kegg <- name_from_hmdb <- name_from_lmsd <- rep(NA_character_, n_row)
            
            # KEGG id -> name
            if (!is.na(input$key_column_kegg)) {
              name_from_kegg <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_kegg]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(kegg_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # HMDB accession -> name
            if (!is.na(input$key_column_hmdb)) {
              name_from_hmdb <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_hmdb]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(hmdb_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # LMSD LM_ID -> name
            if (!is.na(input$key_column_lmsd)) {
              name_from_lmsd <- vapply(
                seq_len(n_row),
                function(i) {
                  val <- metabo.data.raw.all.col[[ion_type]][[input$key_column_lmsd]][i]
                  if (is.na(val) || val == "") return(NA_character_)
                  ids <- strsplit(as.character(val), "[;|,]")[[1]] %>% trimws()
                  hits <- unique(na.omit(lmsd_name_index[ids]))
                  if (length(hits) > 0) paste(hits, collapse = ";") else NA_character_
                },
                character(1)
              )
            }
            
            # 优先级：KEGG > HMDB > LMSD
            filled_name <- dplyr::coalesce(name_from_kegg, name_from_hmdb, name_from_lmsd)
            
            # —— 4) 在得到 filled_name 的基础上，反查/补全 KEGG_ID —— 
            derive_kegg_from_name <- function(nm_str) {
              if (is.na(nm_str) || nm_str == "") return(NA_character_)
              toks <- strsplit(tolower(nm_str), "[;|,]")[[1]] %>% trimws()
              ids <- unique(c(
                unlist(kegg_index_by_name[toks]),
                unlist(hmdb_index_by_name[toks]),
                unlist(lmsd_index_by_name[toks])
              ))
              ids <- ids[!is.na(ids) & ids != ""]
              if (length(ids) == 0) return(NA_character_)
              paste(unique(ids), collapse = "/")     # 多命中用 / 连接
            }
            
            derived_kegg_by_name <- vapply(filled_name, derive_kegg_from_name, FUN.VALUE = character(1))
            
            # —— 5) 组装输出：插入 name，且创建/补全 id_kegg —— 
            #   目标列顺序：peak_name, mz, rt, name, id_kegg, <sample cols...>
            base_cols <- metabo.data.raw.all.col[[ion_type]][, c(
              input$key_column_peakname,
              input$key_column_mz,
              input$key_column_rt
            )]
            
            colnames(base_cols) <- c("peak_name","mz","rt")
            
            # 如果原始里有 id_kegg 列，就取出（以便只对 NA/空值做补全）；否则新建
            has_kegg_col <- !is.na(input$key_column_kegg) &&
              input$key_column_kegg %in% seq_along(colnames(metabo.data.raw.all.col[[ion_type]]))
            
            if (has_kegg_col) {
              old_kegg <- as.character(metabo.data.raw.all.col[[ion_type]][[input$key_column_kegg]])
              old_kegg[is.na(old_kegg) | trimws(old_kegg) == ""] <- NA_character_
              # 优先已存在的 KEGG，其次用 name 推断
              id_kegg_final <- ifelse(!is.na(old_kegg), old_kegg, derived_kegg_by_name)
            } else {
              id_kegg_final <- derived_kegg_by_name
            }
            
            tmp <- cbind(
              base_cols,
              name    = filled_name,
              id_kegg = id_kegg_final
            )
            
            metabo.data.raw[[ion_type]] <- data.frame(
              tmp,
              metabo.data.raw.all.col[[ion_type]][, input$key_column_sample_start:input$key_column_sample_end, drop = FALSE],
              check.names = FALSE
            )
            
            message("已自动生成 name，并据此补全 id_kegg。")
            
          }
          
          if (sum(is.na(metabo.data.raw[[ion_type]]$peak_name))>0) {
            
            metabo.data.raw[[ion_type]]$peak_name[which(is.na(metabo.data.raw[[ion_type]]$peak_name))] <- "NA"
            
          }
          
          metabo.data <- metabo.data.raw[[ion_type]][, c("peak_name", sample.data[[ion_type]]$`File name`)] %>% remove_rownames() %>% column_to_rownames("peak_name")
          
          metabo.data[metabo.data==input$miss.value.handle.type] <- NA
          
          group <- sample.data[[ion_type]]$Condition
          
          # 数据预处理 ----
          # dir.create("2.Data_preprocessing")
          dir.create("2.数据预处理")
          
          ## 缺失值过滤 ----
          # dir.create("2.Data_preprocessing/1.Missing_value_filter")
          dir.create("2.数据预处理/1.缺失值过滤")
          
          if (input[['miss.value.handle.group']]=="inter.group"){
            
            ### 组内缺失值过滤 ----
            missing.index <- list()
            
            for (i in unique(group)) {
              
              missing.index[[i]] <- which(apply(metabo.data[ ,which(group==i)], 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)}))
              
            }
            
            missing.index <- Reduce(union,missing.index)
            
            if (length(missing.index)>0) {
              
              metabo.data2 <- metabo.data[-missing.index,]
              
            }else{
              
              metabo.data2 <- metabo.data
              
            }
            
            metabo.data.filter <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data2), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data2)
            
            char_cols <- sapply(metabo.data.filter, is.character)
            metabo.data.filter[,char_cols] <- lapply(metabo.data.filter[,char_cols], function(col) {
              stringi::stri_enc_toutf8(col, validate = TRUE)
            })
            
            wb <- createWorkbook()
            addWorksheet(wb, "缺失值过滤数据")
            writeData(wb, "缺失值过滤数据", metabo.data.filter, rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/1.缺失值过滤/组内缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          if(input[['miss.value.handle.group']]=="global.group"){
            
            global.judge <- apply(metabo.data, 1, function(x){sum(is.na(x))>(length(x)*input[['miss.value.handle.cutoff']]/100)})
            
            if (sum(global.judge>0)) {
              
              metabo.data2 <- metabo.data[-which(global.judge),]
              
            }else{
              
              metabo.data2 <- metabo.data
              
            }
            
            metabo.data.filter <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data2), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data2)
            
            char_cols <- sapply(metabo.data.filter, is.character)
            metabo.data.filter[,char_cols] <- lapply(metabo.data.filter[,char_cols], function(col) {
              stringi::stri_enc_toutf8(col, validate = TRUE)
            })
            
            wb <- createWorkbook()
            addWorksheet(wb, "缺失值过滤数据")
            writeData(wb, "缺失值过滤数据", metabo.data.filter, rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/1.缺失值过滤/全局缺失值过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          ## 缺失值做填充 ----
          dir.create("2.数据预处理/2.缺失值填充")
          
          ### 不填充 ----
          if (input[['miss.value.fill']]=="none"){
            
            metabo.data.fill <- metabo.data2
            
          }
          
          ### 填充 ----
          if (input[['miss.value.fill']]=="mean.group") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::mutate(group = group[match(sample, colnames(metabo.data2))]) %>%  # 添加分组信息
              dplyr::group_by(metabolite, group) %>%  # 按代谢物和组分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  # 关键修改：对全NA的分组返回0
                  coalesce(mean(value, na.rm = TRUE), 0), 
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="mean.global") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::group_by(metabolite) %>%  # 按代谢物分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(mean(value, na.rm = TRUE), 0),  # 如果全为 NA，填充为 0
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="median.group") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::mutate(group = group[match(sample, colnames(metabo.data2))]) %>%  # 添加分组信息
              dplyr::group_by(metabolite, group) %>%  # 按代谢物和组分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(median(value, na.rm = TRUE), 0), # 组内中位数填充
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="median.global") {
            
            metabo.data.fill <- metabo.data2 %>%
              rownames_to_column(var = "metabolite") %>%  # 保留代谢物名称
              pivot_longer(cols = -metabolite, names_to = "sample", values_to = "value") %>%  # 转换为长格式
              dplyr::group_by(metabolite) %>%  # 按代谢物分组
              dplyr::mutate(
                value = ifelse(
                  is.na(value), 
                  coalesce(median(value, na.rm = TRUE), 0),  # 如果全为 NA，填充为 0
                  value
                )
              ) %>% 
              .[,c("metabolite", "sample", "value")] %>% 
              pivot_wider(names_from = sample, values_from = value) %>%  # 转回宽格式
              dplyr::ungroup() %>%
              column_to_rownames(var = "metabolite")  # 恢复代谢物为行名
            
          }
          
          if (input[['miss.value.fill']]=="min.global") {
            
            metabo.data.fill <- metabo.data2
            
            metabo.data.fill[is.na(metabo.data.fill)] <- min(as.matrix(metabo.data.fill),na.rm = T)
            
          }
          
          if (input[['miss.value.fill']]=="min2") {
            
            metabo.data.fill <- metabo.data2
            
            metabo.data.fill[is.na(metabo.data.fill)] <- mean(as.matrix(metabo.data.fill),na.rm = T)/2
            
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
          
          metabo.data.fill2 <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data.fill)
          
          char_cols <- sapply(metabo.data.fill2, is.character)
          metabo.data.fill2[,char_cols] <- lapply(metabo.data.fill2[,char_cols], function(col) {
            stringi::stri_enc_toutf8(col, validate = TRUE)
          })
          
          wb <- createWorkbook()
          addWorksheet(wb, "缺失值填充数据")
          writeData(wb, "缺失值填充数据", metabo.data.fill2, rowNames = F)
          saveWorkbook(wb, paste0("2.数据预处理/2.缺失值填充/缺失值填充数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
          ## RSD 过滤 ----
          dir.create("2.数据预处理/3.RSD过滤")
          
          if (input$log.rsd.method=="none") {
            
            metabo.data.rsd <- metabo.data.fill %>% dplyr::mutate(qc_rsd = matrixStats::rowSds(as.matrix(metabo.data.fill[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/matrixStats::rowMeans2(as.matrix(metabo.data.fill[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T))
            
          }
          
          if (input$log.rsd.method=="log2") {
            
            metabo.data.rsd <- metabo.data.fill
            
            metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)] <- log2(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]+1)
            
            metabo.data.rsd <- metabo.data.rsd %>%
              dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T))
            
          }
          
          if (input$log.rsd.method=="log10") {
            
            metabo.data.rsd <- metabo.data.fill
            
            metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)] <- log10(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]+1)
            
            metabo.data.rsd <- metabo.data.rsd %>%
              dplyr::mutate(qc_rsd = rowSds(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T)/rowMeans2(as.matrix(metabo.data.rsd[,grepl("^QC$", trimws(group), ignore.case = T)]),na.rm = T))
            
          }
          
          # 正确生成累积分布数据的方法
          cdf_data <- metabo.data.rsd %>% 
            # 确保使用正确的列名（假设实际列名是"qc_rsd"）
            dplyr::arrange(qc_rsd) %>% 
            # 添加频次计数列（如果原数据没有）
            dplyr::mutate(
              bin = cut(qc_rsd, breaks = seq(0, 5, by = 0.01)), # 创建1%间隔的区间
              count = 1  # 每行代表一个峰
            ) %>% 
            dplyr::group_by(bin) %>% 
            dplyr::summarise(
              qc_rsd = mean(qc_rsd, na.rm = TRUE), # 取区间中值
              count = sum(count) # 计算每个区间的峰数量
            ) %>% 
            dplyr::mutate(
              cumulative = cumsum(count)/sum(count)*100 # 计算累积百分比
            ) %>% 
            dplyr::filter(!is.na(bin)) # 移除空区间
          
          metabo.data.rsd <- metabo.data.rsd %>% dplyr::filter(qc_rsd <= input$rsd.cutoff)
          
          metabo.data.rsd2 <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.rsd), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data.rsd)
          
          char_cols <- sapply(metabo.data.rsd2, is.character)
          metabo.data.rsd2[,char_cols] <- lapply(metabo.data.rsd2[,char_cols], function(col) {
            stringi::stri_enc_toutf8(col, validate = TRUE)
          })
          
          wb <- createWorkbook()
          addWorksheet(wb, "RSD过滤")
          writeData(wb, "RSD过滤", metabo.data.rsd2, rowNames = F)
          saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
          
          ## 归一化 ----
          metabo.data.fill <- metabo.data.rsd[,-which(colnames(metabo.data.rsd)=="qc_rsd")]
          
          if (input[['normalized.handle.method']]=="none"){""}else{
            
            dir.create("2.数据预处理/4.归一化")
            
            ### sum归一化
            if (input[['normalized.handle.method']]=="sum") {
              
              sum.median <- apply(metabo.data.fill,2,sum,na.rm=T)
              
              for (i in 1:ncol(metabo.data.fill)) {
                
                metabo.data.fill[,i] <- input$sum_coef*metabo.data.fill[,i]/sum.median[i]
                
              }
              
            }
            
            ### QC样本归一化
            if (input[['normalized.handle.method']]=="qc") {
              
              data_csv <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)], metabo.data.fill)
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              
              sample_info <- sample.data[[ion_type]][,c("File name", "injection.order", "Condition")]%>%dplyr::rename(sample.name=`File name`, class=Condition)
              sample_info[["class"]][grep("^qc$", trimws(sample_info[["class"]]), ignore.case = T)] <- "QC"
              sample_info[["class"]][-grep("^qc$", trimws(sample_info[["class"]]), ignore.case = T)] <- "Subject"
              
              write_csv(data_csv, paste0("data.csv"))
              write_csv(sample_info, paste0("sample.info.csv"))
              
              metNor(
                ms1.data.name = "data.csv",
                sample.info.name = "sample.info.csv",
                minfrac.qc = 0,
                minfrac.sample = 0,
                optimization = TRUE,
                multiple = 5,
                threads = 3
              )
              
              metabo.data.rsd <- metabo.data.fill <- read_csv("svr_normalization_result/data_svr_normalization.csv")
              
              metabo.data.rsd[["QC.nor.rsd"]] <- metabo.data.rsd[["QC.nor.rsd"]]/100
              metabo.data.rsd[["sample.nor.rsd"]] <- metabo.data.rsd[["sample.nor.rsd"]]/100
              
              # 正确生成累积分布数据的方法
              cdf_data <- metabo.data.rsd %>% 
                # 确保使用正确的列名（假设实际列名是"qc_rsd"）
                dplyr::arrange(QC.nor.rsd) %>% 
                # 添加频次计数列（如果原数据没有）
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 5, by = 0.01)), # 创建1%间隔的区间
                  count = 1  # 每行代表一个峰
                ) %>% 
                dplyr::group_by(bin) %>% 
                dplyr::summarise(
                  qc_rsd = mean(QC.nor.rsd, na.rm = TRUE), # 取区间中值
                  count = sum(count) # 计算每个区间的峰数量
                ) %>% 
                dplyr::mutate(
                  cumulative = cumsum(count)/sum(count)*100 # 计算累积百分比
                ) %>% 
                dplyr::filter(!is.na(bin)) # 移除空区间
              
              metabo.data.rsd <- metabo.data.rsd %>% dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              wb <- createWorkbook()
              addWorksheet(wb, "RSD过滤")
              writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = F)
              saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              metabo.data.fill <- metabo.data.rsd
              colnames(metabo.data.fill) <- str_replace(colnames(metabo.data.fill), "^Sample","")
              metabo.data.fill <- metabo.data.fill[,c("name", sample.data[[ion_type]]$`File name`)]%>%column_to_rownames("name")

              system(paste0("mv svr_normalization_result 2.数据预处理/4.归一化/svr_normalization_result_", ion_type))
              
            }
            
            ### QC-RLSC样本归一化
            if (input[['normalized.handle.method']]=="qc-rlsc") {
              
              # metabo.data.fill <- metabo.data.fill.pre
              ## 1) 仍然沿用你已写好的 data.csv / sample.info.csv 生成 ----------
              data_csv <- cbind(
                metabo.data.raw.all.col[[ion_type]][
                  match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]),
                  c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)
                ],
                metabo.data.fill
              )
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              data_csv$name <- as.character(data_csv$name)
              
              sample_info <- sample.data[[ion_type]][,c("File name", "injection.order", "Condition")] %>%
                dplyr::rename(sample.name=`File name`, class=Condition)
              sample_info[["class"]][grep("^qc$", trimws(sample_info[["class"]]), ignore.case = TRUE)] <- "QC"
              sample_info[["class"]][-grep("^qc$", trimws(sample_info[["class"]]), ignore.case = TRUE)] <- "Subject"
              
              readr::write_csv(data_csv,      "data.csv")
              readr::write_csv(sample_info,   "sample.info.csv")
              
              ## 2) 跑 QC-RLSC（qcrlscR::qc.rlsc.wrap），输出到 qcrlsc_normalization_result ----------
              if (!requireNamespace("qcrlscR", quietly = TRUE)) install.packages("qcrlscR")
              
              dir.create("qcrlsc_normalization_result", showWarnings = FALSE)
              
              # a) 读入
              dat_raw  <- readr::read_csv("data.csv",        show_col_types = FALSE)
              meta_raw <- readr::read_csv("sample.info.csv", show_col_types = FALSE)
              
              # b) 按注入顺序排序 & 构造 QC/批次标签
              meta <- meta_raw %>%
                dplyr::transmute(
                  sample = .data[["sample.name"]],
                  class  = ifelse(grepl("(?i)^qc$|\\bqc\\b|[._-]qc\\d*$", .data[["class"]]), "QC", "Sample"),
                  order  = as.numeric(.data[["injection.order"]]),
                  batch  = "1"
                ) %>%
                dplyr::arrange(order)
              
              stopifnot(any(meta$class=="QC"))
              
              # c) 把矩阵整理成 行=样本 × 列=特征
              stopifnot(all(meta$sample %in% colnames(dat_raw)))
              feat_id <- "name"  # 非样本列里，把 "name" 当特征ID
              mat <- as.matrix(dat_raw[, meta$sample, drop = FALSE])   # 列=样本
              rownames(mat) <- dat_raw[[feat_id]]
              X <- t(mat)  # 行=样本，列=特征（qcrlscR 需要）
              
              # d) 轻度缺失过滤（仅用于拟合，输出时会拼回）
              miss_rate <- colMeans(is.na(X))
              keep <- miss_rate < 0.5
              X_fit <- X[, keep, drop = FALSE]
              
              # e) 单批就关掉批内/批间步骤
              cls_qc <- factor(ifelse(meta$class=="QC","qc","sample"), levels=c("qc","sample"))
              cls_bl <- factor(meta$batch)
              intra_use <- FALSE
              shift_use <- FALSE
              
              # f) 跑 QC-RLSC（log10=TRUE + GCV 优化 span + 多项式 degree=2）
              Xc_fit <- qcrlscR::qc.rlsc.wrap(
                dat    = as.data.frame(X_fit),
                cls.qc = cls_qc,
                cls.bl = cls_bl,
                method = "divide",     # Dunn 2011 更推荐等比校正
                intra  = intra_use,    # 单批 FALSE
                opti   = TRUE,         # GCV 自动寻优 span
                log10  = TRUE,
                outl   = TRUE,
                shift  = shift_use,    # 单批 FALSE
                degree = 2
              )
              
              # g) 拼回所有特征
              Xc <- matrix(NA_real_, nrow=nrow(X), ncol=ncol(X), dimnames=dimnames(X))
              Xc[, keep] <- as.matrix(Xc_fit)
              
              # h) 写出 “样本×特征” 和 “特征×样本”
              out_s_by_f <- Xc %>% as.data.frame() %>% tibble::rownames_to_column("name")
              readr::write_csv(out_s_by_f, file.path("qcrlsc_normalization_result","qcrlsc_sample_by_feature.csv"))
              
              out_f_by_s <- t(Xc) %>% as.data.frame() %>% tibble::rownames_to_column("name")
              
              # 统一 key 类型：都用字符
              dat_raw  <- readr::read_csv("data.csv", show_col_types = FALSE) %>%
                dplyr::mutate(name = as.character(name))
              
              out_f_by_s <- out_f_by_s %>%
                dplyr::mutate(name = as.character(name))
              
              data_qcrlsc_norm <- dplyr::left_join(
                dat_raw[, c("name","mz","rt")],
                out_f_by_s, by = "name"
              )
              
              # 为了和你 SVR 的命名风格一致，做一个总表：特征基本信息 + 归一化后的样本强度
              data_qcrlsc_norm <- dplyr::left_join(
                dat_raw[, c("name","mz","rt")],   # 特征信息
                out_f_by_s, by="name"
              )
              readr::write_csv(data_qcrlsc_norm, file.path("qcrlsc_normalization_result","data_qcrlsc_normalization.csv"))
              
              ## 3) 计算 RSD（QC 与 Subject），并生成和你 SVR 分支一致的后续产物 ----------
              # 样本分组
              qc_samples     <- meta$sample[meta$class=="QC"]
              subject_samples<- meta$sample[meta$class!="QC"]
              
              # 取 “特征×样本” 矩阵（不含 name/mz/rt）
              M <- as.matrix(data_qcrlsc_norm[, subject_samples, drop=FALSE])
              M_qc <- as.matrix(data_qcrlsc_norm[, qc_samples, drop=FALSE])
              
              # RSD 函数（百分数）
              rsd_perc <- function(v) {
                v <- as.numeric(v)
                mu <- mean(v, na.rm=TRUE)
                if (!is.finite(mu) || mu<=0) return(NA_real_)
                100 * stats::sd(v, na.rm=TRUE) / mu
              }
              
              QC_rsd_perc     <- apply(M_qc,     1, rsd_perc)
              Sample_rsd_perc <- apply(M,        1, rsd_perc)
              
              metabo.data.rsd <- data_qcrlsc_norm %>%
                dplyr::mutate(
                  `QC.nor.rsd.perc`     = QC_rsd_perc,
                  `sample.nor.rsd.perc` = Sample_rsd_perc,
                  `QC.nor.rsd`          = QC_rsd_perc/100,       # ← 分数形式（0~1）
                  `sample.nor.rsd`      = Sample_rsd_perc/100
                )
              
              # CDF 数据（注意：这里用“分数制”的 QC.nor.rsd，1% = 0.01）
              cdf_data <- metabo.data.rsd %>%
                dplyr::arrange(QC.nor.rsd) %>%
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 0.50, by = 0.01), include.lowest = TRUE)
                ) %>%
                dplyr::filter(!is.na(bin)) %>%
                dplyr::group_by(bin) %>%
                dplyr::summarise(
                  qc_rsd_mid = mean(QC.nor.rsd, na.rm = TRUE),
                  count      = dplyr::n(),
                  .groups    = "drop"
                ) %>%
                dplyr::mutate(cumulative = cumsum(count)/sum(count)*100)
              
              # RSD 过滤（沿用你原来的阈值语义：input$rsd.cutoff 是“分数制”，比如 0.3 代表 30%）
              metabo.data.rsd <- metabo.data.rsd %>%
                dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              # 写 Excel
              wb <- openxlsx::createWorkbook()
              openxlsx::addWorksheet(wb, "RSD过滤")
              openxlsx::writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = FALSE)
              dir.create("2.数据预处理/3.RSD过滤", recursive = TRUE, showWarnings = FALSE)
              openxlsx::saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              # 回填 metabo.data.fill，列只保留样本列
              metabo.data.fill <- metabo.data.rsd
              colnames(metabo.data.fill) <- stringr::str_replace(colnames(metabo.data.fill), "^Sample","")
              metabo.data.fill <- metabo.data.fill[, c("name", sample.data[[ion_type]]$`File name`)] %>% tibble::column_to_rownames("name")
              
              # 挪动输出目录
              system(paste0("mv qcrlsc_normalization_result 2.数据预处理/4.归一化/", "qcrlsc_normalization_result_", ion_type))
              
            }
            
            ### NormAE样本归一化
            if (input[['normalized.handle.method']]=="normae") {
              
              # metabo.data.fill <- metabo.data.fill.pre
              
              data_csv <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), c(input$key_column_peakname, input$key_column_mz, input$key_column_rt)], metabo.data.fill)
              colnames(data_csv)[1:3] <- c("name","mz","rt")
              
              sample_info <- sample.data[[ion_type]][,c("File name", "Condition", "injection.order", "batch")]%>%dplyr::rename(sample.name=`File name`, class=Condition)
              ## 假设你已经生成了 data_csv 和 sample_info（与之前相同）
              ## 1) 基础匹配检查：样本名是否一致
              sample_cols <- colnames(data_csv)[-(1:3)]
              stopifnot(!anyDuplicated(sample_cols))
              stopifnot(!anyDuplicated(sample_info$sample.name))
              
              # 映射并提醒不匹配
              not_in_data   <- setdiff(sample_info$sample.name, sample_cols)
              not_in_sample <- setdiff(sample_cols, sample_info$sample.name)
              if (length(not_in_data) || length(not_in_sample)) {
                cat("【警告】样本名不匹配：\n- 在 sample_info 不在 data 的：", paste(head(not_in_data,20), collapse=", "), "\n",
                    "- 在 data 不在 sample_info 的：", paste(head(not_in_sample,20), collapse=", "), "\n", sep = "")
                stop("请先修正样本名匹配再运行 NormAE。")
              }
              
              ## 2) 类型与取值检查：class/batch/injection.order
              sample_info$class <- ifelse(grepl("^qc$", trimws(sample_info$class), ignore.case = TRUE), "QC", "Subject")
              sample_info$batch <- suppressWarnings(as.integer(sample_info$batch))
              sample_info$batch[is.na(sample_info$batch)] <- 1L
              sample_info$injection.order <- suppressWarnings(as.integer(sample_info$injection.order))
              if (any(is.na(sample_info$injection.order))) {
                stop("injection.order 中存在非数字条目，请先清洗为整数。")
              }
              
              # cat("Class 计数：\n"); print(table(sample_info$class))
              # cat("各 batch×class 计数：\n"); print(table(sample_info$batch, sample_info$class))
              
              ## 3) 数值矩阵清洗：去掉非有限值、零方差行
              X <- as.matrix(data_csv[, -(1:3)])
              mode(X) <- "numeric"
              
              # 标记非有限值
              non_finite_rate <- mean(!is.finite(X))
              cat(sprintf("非有限值比例：%.4f\n", non_finite_rate))
              
              # 将非有限值先设为 NA，计算零方差行
              X[!is.finite(X)] <- NA_real_
              
              # 至少保留 >=3 个有效观测，且范围>0（避免零方差）
              keep <- rowSums(is.finite(X)) >= 3 & apply(X, 1, function(v) {
                rng <- range(v, na.rm = TRUE)
                is.finite(rng[1]) && is.finite(rng[2]) && (diff(rng) > 0)
              })
              
              cat(sprintf("原始特征数：%d，过滤后保留：%d（去掉全NA/零方差）\n", nrow(X), sum(keep)))
              
              # 可选：对剩余少量 NA 做填补（NormAE/Sklearn 一般不接受 NA）
              X <- X[keep, , drop = FALSE]
              X[is.na(X)] <- 0
              
              # 写回 data.csv（前三列为 name/mz/rt）
              data_csv_clean <- cbind(data_csv[keep, 1:3, drop = FALSE], as.data.frame(X, check.names = FALSE))
              readr::write_csv(data_csv_clean, "data.csv")  # 覆盖为清洗后的
              readr::write_csv(sample_info, "sample.info.csv")
              
              dir.create("normae_normalization_result", showWarnings = FALSE)
              
              # —— 用 conda 的 normae 环境运行 —— #
              run_normae <- function(meta = "data.csv", sample = "sample.info.csv", out = "./normae_normalization_result",
                                     env = "normae", log_file = "normae_run.log") {
                meta <- normalizePath(meta); sample <- normalizePath(sample); out <- normalizePath(out)
                log_path <- file.path(out, log_file)
                
                # 找 conda
                conda <- Sys.which("conda")
                if (conda == "") {
                  cand <- path.expand(c("~/anaconda3/bin/conda", "~/miniconda3/bin/conda", "/opt/anaconda3/bin/conda"))
                  conda <- cand[file.exists(cand)][1]
                  if (is.na(conda)) stop("找不到 conda，可临时把 anaconda3/miniconda3/bin 加入 PATH。")
                }
                
                run_cmd <- function(args, tag) {
                  out_lines <- system2(conda, args, stdout = TRUE, stderr = TRUE)
                  cat(out_lines, sep = "\n", file = log_path, append = TRUE)
                  status <- attr(out_lines, "status")
                  list(status = if (is.null(status)) 0L else status, out = out_lines, tag = tag)
                }
                
                base <- c("run", "--no-capture-output", "-n", env)
                
                # 方案1：python -m normae（最稳）
                r1 <- run_cmd(c(base, "python", "-m", "normae",
                                "--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                              "python -m normae")
                if (r1$status == 0) { message("NormAE 完成（python -m normae）。日志：", log_path); return(invisible(r1$out)) }
                
                # 方案2：CLI 脚本 normae
                r2 <- run_cmd(c(base, "normae",
                                "--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                              "normae CLI")
                if (r2$status == 0) { message("NormAE 完成（normae CLI）。日志：", log_path); return(invisible(r2$out)) }
                
                # 方案3：绝对路径兜底（按你的安装路径修改）
                normae_bin <- path.expand(file.path("~", "anaconda3", "envs", env, "bin", "normae"))
                if (file.exists(normae_bin)) {
                  r3 <- system2(normae_bin,
                                c("--meta_csv", meta, "--sample_csv", sample, "--output_dir", out),
                                stdout = TRUE, stderr = TRUE)
                  cat(r3, sep = "\n", file = log_path, append = TRUE)
                  st3 <- attr(r3, "status"); if (is.null(st3)) st3 <- 0L
                  if (st3 == 0) { message("NormAE 完成（绝对路径）。日志：", log_path); return(invisible(r3)) }
                }
                
                # 汇总报错片段
                tail_msg <- paste(tail(c(r1$out, r2$out), 40), collapse = "\n")
                stop("NormAE 运行失败。请查看日志：", log_path, "\n关键信息：\n", tail_msg)
              }
              
              run_normae()
              
              metabo.data.rsd <- read_csv("normae_normalization_result/X_clean.csv")
              
              qc_cols <- sample.data[[ion_type]]$`File name`[grep("QC",sample.data[[ion_type]]$Condition,ignore.case = T)]
              rsd_perc <- function(v){ mu <- mean(v, na.rm=TRUE); if(!is.finite(mu)||mu<=0) return(NA_real_); 100*stats::sd(v,na.rm=TRUE)/mu }
              metabo.data.rsd$QC.nor.rsd <- if (length(qc_cols) >= 2) apply(metabo.data.rsd[, qc_cols, drop=FALSE], 1, rsd_perc) else rep(NA_real_, nrow(metabo.data.rsd))
              
              # CDF 数据（注意：这里用“分数制”的 QC.nor.rsd，1% = 0.01）
              cdf_data <- metabo.data.rsd %>%
                dplyr::arrange(QC.nor.rsd) %>%
                dplyr::mutate(
                  bin = cut(QC.nor.rsd, breaks = seq(0, 0.50, by = 0.01), include.lowest = TRUE)
                ) %>%
                dplyr::filter(!is.na(bin)) %>%
                dplyr::group_by(bin) %>%
                dplyr::summarise(
                  qc_rsd_mid = mean(QC.nor.rsd, na.rm = TRUE),
                  count      = dplyr::n(),
                  .groups    = "drop"
                ) %>%
                dplyr::mutate(cumulative = cumsum(count)/sum(count)*100)
              
              # RSD 过滤（沿用你原来的阈值语义：input$rsd.cutoff 是“分数制”，比如 0.3 代表 30%）
              metabo.data.rsd <- metabo.data.rsd %>%
                dplyr::filter(QC.nor.rsd <= input$rsd.cutoff)
              
              # 写 Excel
              wb <- openxlsx::createWorkbook()
              openxlsx::addWorksheet(wb, "RSD过滤")
              openxlsx::writeData(wb, "RSD过滤", metabo.data.rsd, rowNames = FALSE)
              dir.create("2.数据预处理/3.RSD过滤", recursive = TRUE, showWarnings = FALSE)
              openxlsx::saveWorkbook(wb, paste0("2.数据预处理/3.RSD过滤/qc样本归一化后的RSD过滤数据_", ion_type, ".xlsx"), overwrite = TRUE)
              
              metabo.data.rsd$qc_rsd <- metabo.data.rsd$QC.nor.rsd
              
              # 回填 metabo.data.fill，列只保留样本列
              metabo.data.fill <- metabo.data.rsd
              metabo.data.fill <- metabo.data.fill[, c("name", sample.data[[ion_type]]$`File name`)] %>% tibble::column_to_rownames("name")
              
              system(paste0("mv normae_normalization_result 2.数据预处理/4.归一化/normae_normalization_result_", ion_type))
              
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
            
            ### 75分位数归一化
            if (input[['normalized.handle.method']]=="percent_0.75"){
              
              # 定义处理函数
              normalize_by_quantile <- function(column) {
                # 计算75%分位数
                q75 <- quantile(column, 0.75, na.rm = TRUE)
                # 处理特殊情况：如果分位数为0，则直接返回原始值
                if(q75 == 0) {
                  return(column * 1000)
                } else {
                  # 归一化处理：每个值除以75%分位数，再乘以1000
                  return(column / q75 * 1000)
                }
              }
              
              # 对每一列应用处理函数
              # 使用sapply会返回矩阵，这里保持原始数据框结构
              metabo.data.fill <- as.data.frame(lapply(metabo.data.fill, normalize_by_quantile), row.names = rownames(metabo.data.fill))
              
            }
            
            metabo.data.fill2 <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data.fill)
            
            char_cols <- sapply(metabo.data.fill2, is.character)
            metabo.data.fill2[,char_cols] <- lapply(metabo.data.fill2[,char_cols], function(col) {
              stringi::stri_enc_toutf8(col, validate = TRUE)
            })
            
            wb <- createWorkbook()
            addWorksheet(wb, "归一化处理")
            writeData(wb, "归一化处理", metabo.data.fill2, rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/4.归一化/归一化处理数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          ## 定量值取log ----
          if (input[["log.handle.method"]]=="log2") {
            
            metabo.data.fill <- log2(metabo.data.fill+1)
            
          }
          
          if (input[["log.handle.method"]]=="log10") {
            
            metabo.data.fill <- log10(metabo.data.fill+1)
            
          }
          
          ## 批次矫正 ----
          if (input[['batch.correct']]=="none"){""}else{
            
            dir.create("2.数据预处理/5.批次校正")
            
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
            
            metabo.data.fill2 <- cbind(metabo.data.raw.all.col[[ion_type]][match(rownames(metabo.data.fill), metabo.data.raw.all.col[[ion_type]][[input$key_column_peakname]]), 1:(input$key_column_sample_start-1)], metabo.data.fill)
            
            char_cols <- sapply(metabo.data.fill2, is.character)
            metabo.data.fill2[,char_cols] <- lapply(metabo.data.fill2[,char_cols], function(col) {
              stringi::stri_enc_toutf8(col, validate = TRUE)
            })
            
            wb <- createWorkbook()
            addWorksheet(wb, "批次矫正")
            writeData(wb, "批次矫正", metabo.data.fill2, rowNames = F)
            saveWorkbook(wb, paste0("2.数据预处理/5.批次矫正/批次矫正数据_", ion_type, ".xlsx"), overwrite = TRUE)
            
          }
          
          ## 保存变量 ----
          assign(paste0("metabo.data.", ion_type), metabo.data.fill)
          
          # QC 质控 ----
          dir.create("3.QC")
          
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
          dir.create("3.QC/1.TIC")
          
          QC.RT <- data.frame(EG.ApexRT=metabo.data.raw[[ion_type]]$rt[match(rownames(metabo.data.fill), metabo.data.raw[[ion_type]]$peak_name)]%>%as.numeric(), metabo.data.fill[, grepl("^qc$", trimws(group), ignore.case = T)]) %>% reshape2::melt(.,id=c("EG.ApexRT"),variable.name = "R.FileName",value.name = "FG.MS1Quantity")
          
          QC.RT2 <- dplyr::mutate(QC.RT,EG.ApexRT.bin=round(EG.ApexRT,digits = 0))
          
          QC.RT3 <- QC.RT2 %>% dplyr::group_by(R.FileName,EG.ApexRT.bin)%>%dplyr::summarize(Quantity.sum=sum(FG.MS1Quantity,na.rm=T))
          
          QC.RT3$type1<-"MS1"
          
          QC.TIC <- QC.RT3
          
          g_TIC_MS1 <- {
            
            highchart() %>% 
              hc_add_yAxis(lineWidth = 3,title = list(text = "Intensity"))%>%
              hc_xAxis(title = list(text = "Retention Time")) %>%
              hc_add_series(data = QC.TIC[which(QC.TIC$type1=="MS1"),],type="spline",hcaes(x = EG.ApexRT.bin, y = Quantity.sum, group = R.FileName)) %>%
              hc_tooltip(split=T,valueDecimals=3)%>%
              hc_plotOptions(series = list(marker = list(symbol = "circle")))%>%
              hc_chart(zoomType="x",events=list(selection=select_zoom))%>%
              hc_exporting(enabled = TRUE,buttons = list( contextButton = list(menuItems = list('downloadPNG', 'downloadSVG',"downloadPDF","downloadJPEG","printChart","viewFullscreen"))))
            
          }
          
          htmlwidgets::saveWidget(widget = g_TIC_MS1, file = paste0("3.QC/1.TIC/TIC_", ion_type, ".html"))
          
          webshot::webshot(url = paste0("3.QC/1.TIC/TIC_", ion_type, ".html"), 
                           file = c(paste0("3.QC/1.TIC/TIC_", ion_type, ".png"),
                                    paste0("3.QC/1.TIC/TIC_", ion_type, ".pdf")),
                           delay = 3)
          
          system(paste0("rm -rf '3.QC/1.TIC/TIC_", ion_type, "_files'"))
          
          write_csv(QC.TIC, paste0("3.QC/1.TIC/TIC_", ion_type, ".csv"))
          
          ## QC样本的相关性图 ----
          dir.create(paste0("3.QC/2.correlation/", ion_type),recursive = T,showWarnings = F)
          
          result<-corr.test(metabo.data.fill[, grepl("^qc$", trimws(group), ignore.case = T)], method = "pearson",adjust="none",alpha=.05)
          rmatrix<-result$r
          pmatrix<-result$p
          
          col <- colorRampPalette(c("darkblue", "white", "red"))(200)
          
          pdf(paste0("3.QC/2.correlation/", ion_type, "/correlation_", ion_type, ".pdf"))
          p1 <- corrplot::corrplot(rmatrix,method = "circle",col = col,tl.col="black",type="upper")
          dev.off()
          png(paste0("3.QC/2.correlation/", ion_type, "/correlation_", ion_type, ".png"))
          p1 <- corrplot::corrplot(rmatrix,method ="circle",col = col,tl.col="black",type="upper")
          dev.off()
          write.csv(rmatrix, paste0("3.QC/2.correlation/", ion_type, "/correlation_", ion_type, ".csv"))
          write.csv(pmatrix, paste0("3.QC/2.correlation/", ion_type, "/correlation_pvalue_", ion_type, ".csv"))
          
          ## QC样本的RSD分布图 ----
          dir.create(paste0("3.QC/3.RSD/", ion_type),recursive = T,showWarnings = F)
          
          # 添加30%参考线配置
          vline_30 <- list(
            label = list(text = "Acceptance criteria (30%)"),
            color = "#666666",
            dashStyle = "Dash",
            width = 2,
            value = 30
          )
          
          # 生成最终图表
          g_rsd <- highchart() %>%
            hc_add_series(cdf_data, hcaes(x = qc_rsd*100, y = cumulative), 
                          type = "line", name = "Cumulative %", 
                          lineWidth = 2,
                          color = "#FFD700") %>%
            hc_xAxis(title = list(text = "RSD (%)"),
                     plotLines = list(vline_30),
                     min = 0, max = 100,
                     tickInterval = 10) %>%  # 添加刻度间隔
            hc_yAxis(title = list(text = "% of peaks"),
                     min = 0, max = 100,
                     labels = list(format = "{value}%"),
                     tickInterval = 25) %>%  # 匹配图片刻度
            hc_tooltip(
              headerFormat = "RSD: {point.x:.1f}%<br/>",
              pointFormat = "累积比例: {point.y:.1f}%",
              valueDecimals = 1
            ) %>% hc_plotOptions(line = list(marker = list(enabled = FALSE))) %>%
            hc_chart(
              style = list(
                # fontFamily = "Arial",
                fontSize = "18px"  # 全局基准字号[4,6](@ref)
              )
            )
          
          htmlwidgets::saveWidget(widget = g_rsd, file = paste0("3.QC/3.RSD/", ion_type, "/RSD_", ion_type, ".html"))
          
          webshot::webshot(url = paste0("3.QC/3.RSD/", ion_type, "/RSD_", ion_type, ".html"), 
                           file = c(paste0("3.QC/3.RSD/", ion_type, "/RSD_", ion_type, ".png"),
                                    paste0("3.QC/3.RSD/", ion_type, "/RSD_", ion_type, ".pdf")),
                           delay = 3)
          
          system(paste0("rm -rf '3.QC/3.RSD/", ion_type, "/RSD_", ion_type, "_files'"))
          
          write.csv(cdf_data, paste0("3.QC/3.RSD/", ion_type, "/RSD_", ion_type, ".csv"))
          
          ## PCA ----
          dir.create(paste0("3.QC/4.PCA/", ion_type),recursive = T,showWarnings = F)
          
          metabo.data.fill2 <- metabo.data.fill
          
          metabo.data.fill2[is.na(metabo.data.fill2)]<-0
          
          metabo.data.fill2[metabo.data.fill2=="NaN"]<-0
          
          metabo.data.fill2[,1:ncol(metabo.data.fill2)]<-lapply(metabo.data.fill2[,1:ncol(metabo.data.fill2)],as.numeric)
          
          pc.cr <- prcomp(t(metabo.data.fill2), scale. = T)
          
          pca_group <- group
          
          # 自动生成颜色映射规则（无需预知分组名）
          generate_color_mapping <- function(groups) {
            # 按字母顺序排序保证可重复性
            sorted_groups <- sort(unique(as.character(groups)))
            
            # 核心配色（参考图片中的黄金色系）
            base_colors <- c("#f8766d", "#1F77B4")  # 黄金 + 经典蓝
            
            # 动态扩展配色方案
            if(length(sorted_groups) > 2) {
              extended_colors <- viridis::viridis(
                n = length(sorted_groups) - 2,
                begin = 0.2,  # 避免过亮黄色
                end = 0.8,    # 避免过深紫色
                option = "A"  # viridis方案
              )
              color_palette <- c(base_colors, extended_colors)
            } else {
              color_palette <- base_colors
            }
            
            # 创建命名向量
            setNames(color_palette, sorted_groups)
          }
          
          dynamic_colors <- generate_color_mapping(pca_group)
          
          g <- ggord::ggord(pc.cr, pca_group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
            scale_color_manual(values = dynamic_colors) + 
            geom_text(aes(label=lab),vjust=-0.5,hjust="inward",check_overlap = T,size=3,color="black",show.legend=T) +
            theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8)) +
            labs(title="PCA Plot",fill=NULL,colour=NULL,shape=NULL)
          
          g_pca <- g %>% ggplotly()
          
          g_pca$x$data[[3]]$showlegend<-F
          g_pca$x$data[[4]]$showlegend<-F
          g_pca$x$data[[5]]$showlegend<-T
          g_pca$x$data[[5]]$name<-"label"
          
          htmlwidgets::saveWidget(widget = g_pca, file = paste0("3.QC/4.PCA/", ion_type, "/PCA_", ion_type, ".html"))
          
          system(paste0("rm -rf '3.QC/4.PCA/", ion_type, "/PCA_", ion_type, "_files'"))
          
          g_pca <- ggord::ggord(pc.cr, pca_group, obslab=F, arrow=NULL, txt=NULL, ellipse =T, poly=F, coord_fix=F, facet = F, size=4, repel=T) +
            # geom_text(aes(label=lab),vjust=-0.5,hjust="inward",check_overlap = T,size=3,color="black",show.legend=T) +
            geom_text_repel(aes(label=lab))+
            scale_color_manual(values = dynamic_colors) + 
            theme(plot.title = element_text(hjust = 0.5, size=12), text = element_text(size = 8)) +
            labs(title="PCA Plot",fill=NULL,colour=NULL,shape=NULL)
          
          ggsave(paste0("3.QC/4.PCA/", ion_type, "/PCA_", ion_type, ".png"), g_pca, width = 10, height = 6)
          ggsave(paste0("3.QC/4.PCA/", ion_type, "/PCA_", ion_type, ".pdf"), g_pca, width = 10, height = 6)
          
          names(pca_group) <- colnames(metabo.data.fill2)
          
          # 取PC1/PC2坐标 + 分组 + 样本名
          scores <- as.data.frame(pc.cr$x[, 1:2])
          colnames(scores) <- c("PC1", "PC2")
          scores$sample <- rownames(scores)
          scores$group  <- as.character(pca_group[rownames(scores)])
          
          # 按组标记是否在level置信椭圆之外（默认95%）
          flag_outside_by_group <- function(df, level = 0.95) {
            # 当前组的中心和协方差
            ctr <- colMeans(df[, c("PC1", "PC2")], na.rm = TRUE)
            cv  <- stats::cov(df[, c("PC1", "PC2")], use = "complete.obs")
            # 马氏距离平方
            md2 <- stats::mahalanobis(df[, c("PC1","PC2")], center = ctr, cov = cv)
            # 椭圆阈值（df=2）
            thr <- stats::qchisq(level, df = 2)
            df$md2 <- md2
            df$outside <- md2 > thr
            df
          }
          
          # 分组计算；需要至少3个点才能求协方差，少于3个的组全部标记为NA
          scores_flagged <- scores %>%
            dplyr::group_by(group) %>%
            dplyr::group_modify(~ if (nrow(.x) >= 3) {flag_outside_by_group(.x, level = 0.95)} else dplyr::mutate(.x, md2 = NA_real_, outside = NA)) %>%
            ungroup()
          
          # 取置信椭圆外的样本
          outside_samples <- scores_flagged %>% filter(outside %in% TRUE)
          
          # 查看或保存
          # print(outside_samples[, c("sample", "group", "PC1", "PC2", "md2")])
          write.csv(outside_samples, file = paste0("3.QC/4.PCA/", ion_type, "/PCA_outside_95.csv"), row.names = FALSE)
          
          # 判断是否除QC外有多组数据 ----
          if (length(unique(sample.data[[ion_type]][["Condition"]][!grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T)]))>1) {
            
            # 多元统计分析 ----
            dir.create("4.多元统计分析")
            
            ## pca ----
            dir.create(paste0("4.多元统计分析/1.PCA/all/", ion_type), recursive = T, showWarnings = F)
            
            ### 整体样本 ----
            pca_data <- metabo.data.fill2[, !grepl("^qc$", trimws(group), ignore.case = T)]
            
            NcrossvalI <- min(7,ncol(pca_data))
            
            model_pca <- opls(t(pca_data), scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
            
            if(ncol(model_pca@scoreMN)==1){
              
              model_pca <- opls(t(pca_data), predI=2, scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
              
            }
            
            if(ncol(model_pca@scoreMN)==0){
              
              model_pca <- opls(t(pca_data), predI=3, scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
              
            }
            
            df_pca <- model_pca@scoreMN %>% data.frame() %>% 
              rownames_to_column('sample') %>% 
              inner_join(sample.data[[ion_type]] %>% dplyr::select(`File name`, Condition), by=c("sample"="File name")) %>% dplyr::rename("Groups"="Condition")
            
            g_pca <- ggplot(df_pca, aes(x = p1, y = p2, color = Groups)) +
              geom_point(size = 4) +
              stat_ellipse(level = 0.95, alpha = 0.8) +  # 添加95%置信椭圆
              geom_text_repel(aes(label=sample),show.legend = F)+
              scale_color_manual(values = dynamic_colors) + 
              labs(title = "PCA Score Plot", 
                   x = glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'),
                   y = glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)')) +
              theme_bw() +
              theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 12))
            
            ggsave(paste0("4.多元统计分析/1.PCA/all/", ion_type, "/PCA_all_", ion_type, ".png"), g_pca, width = 10, height = 6)
            ggsave(paste0("4.多元统计分析/1.PCA/all/", ion_type, "/PCA_all_", ion_type, ".pdf"), g_pca, width = 10, height = 6)
            
            sum_pca[[paste0("all sample ", ion_type)]] <- data.frame(title=paste0("all sample ", ion_type), type="PCA", A=nrow(df_pca), N=nrow(model_pca@modelDF), `R2X(cum)`=model_pca@summaryDF$`R2X(cum)`)
            
            ### 分组样本 ----
            for (g in 1:length(group.data)) {
              
              dir.create(paste0("4.多元统计分析/1.PCA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              
              experimental <- str_split(group.data[g], pattern = "_vs_") %>% unlist %>% .[1]
              control <- str_split(group.data[g], pattern = "_vs_") %>% unlist %>% .[2]
              
              pca_data <- metabo.data.fill2[, which(group %in% c(experimental, control))]
              
              NcrossvalI <- min(7,ncol(pca_data))
              
              model_pca <- opls(t(pca_data), scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
              
              if(ncol(model_pca@scoreMN)==1){
                
                model_pca <- opls(t(pca_data), predI=2, scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
                
              }
              
              if(ncol(model_pca@scoreMN)==0){
                
                model_pca <- opls(t(pca_data), predI=3, scaleC = 'standard', fig.pdfC = 'none',info.txtC = 'none',crossvalI = NcrossvalI, permI = 200)
                
              }
              
              sum_pca[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                            type="PCA", 
                                                                                            A=nrow(df_pca), 
                                                                                            N=nrow(model_pca@modelDF), 
                                                                                            `R2X(cum)`=model_pca@summaryDF$`R2X(cum)`)
              
              df_pca <- model_pca@scoreMN %>% data.frame() %>% 
                rownames_to_column('sample') %>% 
                inner_join(sample.data[[ion_type]] %>% dplyr::select(`File name`, Condition), by=c("sample"="File name"))%>%dplyr::rename("Groups"="Condition")
              
              g_pca <- ggplot(df_pca, aes(x = p1, y = p2, color = Groups)) +
                geom_point(size = 4) +
                stat_ellipse(level = 0.95, alpha = 0.8) +  # 添加95%置信椭圆
                geom_text_repel(aes(label=sample),show.legend = F)+
                scale_color_manual(values = dynamic_colors) + 
                labs(title = "PCA Score Plot", 
                     x = glue('PCA 1 ({round(model_pca@modelDF$R2X[1]*100)}%)'),
                     y = glue('PCA 2 ({round(model_pca@modelDF$R2X[2]*100)}%)')) +
                theme_bw() +
                theme(plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 12))
              
              ggsave(g_pca,filename = paste0("4.多元统计分析/1.PCA/", group.data[g], "/", ion_type, "/PCA_", experimental, "_vs_", control, "_", ion_type, ".pdf"), width = 10, height = 6)
              ggsave(g_pca,filename = paste0("4.多元统计分析/1.PCA/", group.data[g], "/", ion_type, "/PCA_", experimental, "_vs_", control, "_", ion_type, ".png"), width = 10, height = 6)
              
            }
            
            ## pls-da/opls-da ----
            dir.create(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type), recursive = T, showWarnings = F)
            dir.create("4.多元统计分析/3.OPLS-DA")
            
            ### 所有样本 ----
            # 获取plsda数据
            plsda_data <- metabo.data.fill2[, !grepl("^qc$", trimws(group), ignore.case = T)]
            
            # 获取分组信息
            plsda_group <- group[!grepl("^qc$", trimws(group), ignore.case = T)]
            
            NcrossvalI <- min(7,ncol(plsda_data))
            
            # 运行PLS-DA模型
            plsda_model <- opls(t(plsda_data), plsda_group, crossvalI=NcrossvalI, permI = 200)
            
            if (nrow(plsda_model@modelDF)==0 | nrow(plsda_model@modelDF)==1) {
              
              plsda_model <- opls(t(plsda_data), plsda_group, predI = 2, crossvalI=NcrossvalI, permI = 200)
              
            }
            
            # if (nrow(plsda_model@modelDF)==0) {
            #   
            #   plsda_model <- opls(t(plsda_data), plsda_group, predI = 1, crossvalI=NcrossvalI, permI = 200)
            #   
            # }
            
            if (nrow(plsda_model@modelDF)>=2){
              
              # 从模型对象中提取得分矩阵
              pls_scores <- data.frame(
                p1 = plsda_model@scoreMN[, 1],
                p2 = plsda_model@scoreMN[, 2],
                Groups = plsda_group,
                sample = rownames(t(plsda_data))
              )
              
              x_lab <- paste0("P1(", round(plsda_model@modelDF$R2X[1]*100, 1), "*'%')")
              y_lab <- paste0("P2(", round(plsda_model@modelDF$R2X[2]*100, 1), "*'%')")
              
              g_plsda <- ggplot(pls_scores, aes(x = p1, y = p2, color = Groups)) +
                geom_point(size = 4) +
                stat_ellipse(level = 0.95, alpha = 0.8) +
                geom_text_repel(aes(label = sample)) +
                scale_color_manual(values = dynamic_colors) + 
                labs(
                  x = parse(text = x_lab),  # 解析数学表达式
                  y = parse(text = y_lab),
                  title = "PLS-DA Score Plot"
                ) +
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))
              
              ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, ".png"), g_plsda, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, ".pdf"), g_plsda, width = 10, height = 6)
              
            }
            
            sum_plsda[[paste0("all sample ", ion_type)]] <- data.frame(title=paste0("all sample ", ion_type), 
                                                                       type="PLS-DA", 
                                                                       A=nrow(pls_scores), 
                                                                       N=nrow(plsda_model@modelDF), 
                                                                       `R2X(cum)`=plsda_model@summaryDF$`R2X(cum)`, 
                                                                       `R2Y(cum)`=plsda_model@summaryDF$`R2Y(cum)`, 
                                                                       `Q2(cum)`=plsda_model@summaryDF$`Q2(cum)`)
            
            # 从模型对象提取置换检验结果
            perm_data <- plsda_model@suppLs$permMN %>% 
              as.data.frame() %>% 
              tibble::rownames_to_column("perm_id") %>% 
              tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                  names_to = "metric", 
                                  values_to = "value")
            
            # 计算原始模型指标
            orig_r2y <- plsda_model@summaryDF$`R2Y(cum)`
            orig_q2 <- plsda_model@summaryDF$`Q2(cum)`
            
            # 生成置换检验图
            perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
              # 绘制置换散点
              geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
              scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
              # 添加回归线
              geom_smooth(
                aes(group = metric, color = metric),
                method = "lm", 
                formula = y ~ x,
                se = FALSE, 
                fullrange = TRUE,  # 启用全范围延伸
                # linewidth = 1.2,   # 加粗线宽匹配图片
                linetype = "dashed"
              ) +
              # 标注原始模型值（精确坐标定位）
              annotate("point", x = orig_r2y, y = orig_r2y, 
                       color = "#00BFC4", size = 3, shape = 18) +  # R2Y菱形标记
              annotate("point", x = orig_q2, y = 0, 
                       color = "#7CAE00", size = 3, shape = 18) +  # Q2菱形标记
              # 坐标轴设置（匹配图片范围）
              scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
              scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
              # 标签系统（完全复现图片文字）
              labs(title = "permutation_test",
                   x = "200 permutations 1 components",
                   y = "value",
                   color = "Metric") +
              # 主题优化（精确字体和网格匹配）
              theme_bw() +
              theme(panel.grid.minor = element_blank(),
                    axis.text = element_text(size = 10, color = "black"),
                    legend.position = "right",plot.title = element_text(hjust = 0.5))
            
            # 输出高清图
            ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
            ggsave(paste0("4.多元统计分析/2.PLS-DA/all/", ion_type, "/PLS-DA_all_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
            
            ### 分组样本（PLSDA和OPLS-DA） ----
            for (g in 1:length(group.data)) {
              
              dir.create(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              dir.create(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type), recursive = T, showWarnings = F)
              
              experimental <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[1]
              control <- str_split(group.data[g], pattern = "_vs_")%>%unlist%>%.[2]
              
              plsda_data <- metabo.data.fill2[, which(group %in% c(experimental, control))]
              
              # 获取分组信息
              plsda_group <- group[which(group %in% c(experimental, control))]
              
              NcrossvalI <- min(7,ncol(plsda_data))
              
              #### 运行PLS-DA模型 ----
              plsda_model <- opls(t(plsda_data), plsda_group, crossvalI=NcrossvalI, permI = 200)
              
              if (nrow(plsda_model@modelDF)==0 | nrow(plsda_model@modelDF)==1) {
                
                plsda_model <- opls(t(plsda_data), plsda_group, predI = 2, crossvalI=NcrossvalI, permI = 200)
                
              }
              
              # if (nrow(plsda_model@modelDF)==0) {
              #   
              #   plsda_model <- opls(t(plsda_data), plsda_group, predI = 1, crossvalI=NcrossvalI, permI = 200)
              #   
              # }
              
              if (nrow(plsda_model@modelDF)>=2) {
                
                # 从模型对象中提取得分矩阵
                pls_scores <- data.frame(
                  p1 = plsda_model@scoreMN[, 1],
                  p2 = plsda_model@scoreMN[, 2],
                  Groups = plsda_group,
                  sample = rownames(t(plsda_data))
                )
                
                x_lab <- paste0("P1(", round(plsda_model@modelDF$R2X[1]*100, 1), "*'%')")
                y_lab <- paste0("P2(", round(plsda_model@modelDF$R2X[2]*100, 1), "*'%')")
                
                g_plsda <- ggplot(pls_scores, aes(x = p1, y = p2, color = Groups)) +
                  geom_point(size = 4) +
                  stat_ellipse(level = 0.95, alpha = 0.8) +
                  scale_color_manual(values = dynamic_colors) + 
                  geom_text_repel(aes(label = sample)) +
                  labs(
                    x = parse(text = x_lab),  # 解析数学表达式
                    y = parse(text = y_lab),
                    title = "PLS-DA Score Plot"
                  ) +
                  theme_bw()+theme(plot.title = element_text(hjust = 0.5))
                
                ggsave(g_plsda,filename = paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, ".pdf"), width = 10, height = 6)
                ggsave(g_plsda,filename = paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, ".png"), width = 10, height = 6)
                
              }
              
              sum_plsda[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                              type="PLS-DA", 
                                                                                              A=nrow(pls_scores), 
                                                                                              N=nrow(plsda_model@modelDF), 
                                                                                              `R2X(cum)`=plsda_model@summaryDF$`R2X(cum)`, 
                                                                                              `R2Y(cum)`=plsda_model@summaryDF$`R2Y(cum)`, 
                                                                                              `Q2(cum)`=plsda_model@summaryDF$`Q2(cum)`)
              
              # 从模型对象提取置换检验结果
              perm_data <- plsda_model@suppLs$permMN %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("perm_id") %>% 
                tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                    names_to = "metric", 
                                    values_to = "value")
              
              # 计算原始模型指标
              orig_r2y <- plsda_model@summaryDF$`R2Y(cum)`
              orig_q2 <- plsda_model@summaryDF$`Q2(cum)`
              
              # 生成置换检验图
              perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
                # 绘制置换散点
                geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
                scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
                # 添加回归线
                geom_smooth(
                  aes(group = metric, color = metric),
                  method = "lm", 
                  formula = y ~ x,
                  se = FALSE, 
                  fullrange = TRUE,  # 启用全范围延伸
                  # linewidth = 1.2,   # 加粗线宽匹配图片
                  linetype = "dashed"
                ) +
                # 坐标轴设置（匹配图片范围）
                scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
                scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
                # 标签系统（完全复现图片文字）
                labs(title = "permutation_test",
                     subtitle = paste0("R2Y(cum)=", orig_r2y, "; Q2(cum)=", orig_q2),
                     x = "200 permutations 1 components",
                     y = "value",
                     color = "Metric") +
                # 主题优化（精确字体和网格匹配）
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 10, color = "black"),
                      legend.position = "right",
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
              
              # 输出高清图
              ggsave(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/2.PLS-DA/", group.data[g], "/", ion_type, "/PLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
              
              #### 运行OPLS-DA模型 ----
              oplsda_data <- metabo.data.fill2[, which(group %in% c(experimental, control))]
              
              # 获取分组信息
              oplsda_group <- group[which(group %in% c(experimental, control))]
              
              NcrossvalI <- min(7,ncol(oplsda_data))
              
              oplsda_model <- opls(t(oplsda_data), oplsda_group, orthoI=NA, crossvalI=NcrossvalI, permI = 200)
              
              if (nrow(oplsda_model@modelDF)==0) {
                
                oplsda_model <- opls(t(oplsda_data), oplsda_group, predI = 1, orthoI=1, crossvalI=NcrossvalI, permI = 200)
                
              }
              
              # 从模型对象中提取得分矩阵
              opls_scores <- data.frame(
                p1 = oplsda_model@scoreMN[, 1],       # 预测主成分得分
                orth1 = oplsda_model@orthoScoreMN[,1], # 正交主成分得分
                Groups = oplsda_group,
                sample = rownames(t(oplsda_data))
              )
              
              x_lab <- paste0("t[1]~(", round(oplsda_model@modelDF$R2X[1]*100, 1), "*'%')")
              y_lab <- paste0("to[1]~(", round(oplsda_model@modelDF$R2X[2]*100, 1), "*'%')")
              
              g_oplsda <- ggplot(opls_scores, aes(x = p1, y = orth1, color = Groups)) +
                geom_point(size = 4) + 
                scale_color_manual(values = dynamic_colors) + 
                stat_ellipse(level = 0.95, alpha = 0.8) +
                geom_text_repel(aes(label = sample)) +
                labs(
                  x = parse(text = x_lab),  # 解析数学表达式
                  y = parse(text = y_lab),
                  title = "OPLS-DA Score Plot"
                ) +
                theme_bw()+theme(plot.title = element_text(hjust = 0.5))
              
              ggsave(g_oplsda, filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, ".pdf"), width = 10, height = 6)
              ggsave(g_oplsda, filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, ".png"), width = 10, height = 6)
              
              sum_oplsda[[paste0(experimental, "_vs_", control, "_", ion_type)]] <- data.frame(title=paste0(experimental, "_vs_", control, "_", ion_type), 
                                                                                               type="OPLS-DA", 
                                                                                               A=nrow(opls_scores), 
                                                                                               N=paste0(sum(grepl("^p",rownames(oplsda_model@modelDF))), "+", sum(grepl("^o",rownames(oplsda_model@modelDF)))), 
                                                                                               `R2X(cum)`=oplsda_model@summaryDF$`R2X(cum)`, 
                                                                                               `R2Y(cum)`=oplsda_model@summaryDF$`R2Y(cum)`, `Q2(cum)`=oplsda_model@summaryDF$`Q2(cum)`)
              
              # 从模型对象提取置换检验结果
              perm_data <- oplsda_model@suppLs$permMN %>% 
                as.data.frame() %>% 
                tibble::rownames_to_column("perm_id") %>% 
                tidyr::pivot_longer(cols = c(`R2Y(cum)`, `Q2(cum)`),
                                    names_to = "metric", 
                                    values_to = "value")
              
              # 计算原始模型指标
              orig_r2y <- oplsda_model@summaryDF$`R2Y(cum)`
              orig_q2 <- oplsda_model@summaryDF$`Q2(cum)`
              
              # 生成置换检验图
              perm_plot <- ggplot(perm_data, aes(x = sim, y = value)) +
                # 绘制置换散点
                geom_point(aes(color = metric), alpha = 0.6, size = 2.5) +
                scale_color_manual(values = c("R2Y(cum)" = "#00BFC4", "Q2(cum)" = "#7CAE00")) +
                # 添加回归线
                geom_smooth(
                  aes(group = metric, color = metric),
                  method = "lm", 
                  formula = y ~ x,
                  se = FALSE, 
                  fullrange = TRUE,  # 启用全范围延伸
                  # linewidth = 1.2,   # 加粗线宽匹配图片
                  linetype = "dashed"
                ) +
                # 坐标轴设置（匹配图片范围）
                scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
                scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
                # 标签系统（完全复现图片文字）
                labs(title = "permutation_test",
                     subtitle = paste0("R2Y(cum)=", orig_r2y, "; Q2(cum)=", orig_q2),
                     x = "200 permutations 1 components",
                     y = "value",
                     color = "Metric") +
                # 主题优化（精确字体和网格匹配）
                theme_bw() +
                theme(panel.grid.minor = element_blank(),
                      axis.text = element_text(size = 10, color = "black"),
                      legend.position = "right",
                      plot.title = element_text(hjust = 0.5),
                      plot.subtitle = element_text(hjust = 0.5))
              
              # 输出高清图
              ggsave(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.png"), perm_plot, dpi = 300, width = 10, height = 6)
              ggsave(paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/OPLS-DA_", experimental, "_vs_", control, "_", ion_type, "_permutation_plot.pdf"), perm_plot, dpi = 300, width = 10, height = 6)
              
              # S-Plot
              loadings <- oplsda_model@loadingMN[, 1]  # 第一预测成分的载荷（协方差）
              score_p1 <- oplsda_model@scoreMN[, 1]      # 第一预测成分得分
              vip_values <- oplsda_model@vipVn           # VIP值
              
              # 提取变量名称（确保与载荷、标准差、VIP的名称一致）
              variable_names <- names(loadings)
              
              p_cov <- cov(score_p1, t(oplsda_data[names(loadings),]))%>%as.numeric()
              p_corr <- cor(score_p1, t(oplsda_data[names(loadings),]))%>%as.numeric()
              
              # 构建S-Plot数据框
              splot_data <- data.frame(
                Variable = variable_names,
                Loading = p_cov,
                Correlation = p_corr,
                VIP = vip_values
              )
              
              # 标记VIP > 1.0的变量
              splot_data$Significant <- ifelse(splot_data$VIP > 1.0, "Yes", "No")
              
              # 绘制S-Plot
              splot <- ggplot(splot_data, aes(x = Loading, y = Correlation)) +
                # 绘制点，按VIP显著性着色
                geom_point(aes(color = Significant), size = 2.5, alpha = 0.7) +
                # 设置颜色方案（红色：显著；灰色：不显著）
                scale_color_manual(values = c("Yes" = "red", "No" = "gray50")) +
                # 添加参考线（x=0, y=0）
                geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
                geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
                # 设置主题和标签
                theme_minimal() +
                labs(
                  title = "OPLS-DA S-Plot",
                  x = "Loading [p1]",
                  y = "Correlation [p1]",
                  color = "VIP > 1.0"
                ) +
                theme(
                  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.title = element_text(size = 10),
                  legend.position = "bottom"
                )
              
              # 标记VIP值最高的前5个变量
              top_n <- 5  # 可调整为你需要的数量
              top_vars <- splot_data %>% 
                arrange(desc(VIP)) %>% 
                head(top_n)
              
              if (nrow(top_vars) > 0) {
                splot <- splot +
                  geom_text(data = top_vars, 
                            aes(label = Variable), 
                            hjust = 0.5, vjust = -0.7,  # 标签位置在点下方
                            size = 3, color = "black",
                            check_overlap = TRUE  # 避免标签重叠
                  )
              }
              
              # 保存图形
              ggsave(
                plot = splot,
                filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_", experimental, "_vs_", control, "_", ion_type, ".pdf"),
                width = 8, height = 6
              )
              ggsave(
                plot = splot,
                filename = paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_", experimental, "_vs_", control, "_", ion_type, ".png"),
                width = 8, height = 6, dpi = 300
              )
              
              # 保存S-Plot数据
              write.csv(
                splot_data,
                paste0("4.多元统计分析/3.OPLS-DA/", group.data[g], "/", ion_type, "/S-Plot_Data_", experimental, "_vs_", control, "_", ion_type, ".csv"),
                row.names = FALSE
              )
              
            }
            
            # 差异筛选 ----
            dir.create("5.差异筛选")
            
            dataMat <- metabo.data.fill
            
            for (g in 1:length(group.data)){
              
              cmp_name <- group.data[g]
              expGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
              ctrlGroup <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
              
              # 创建文件夹（如果不存在）用于保存当前比较的结果
              dir.create(paste0("5.差异筛选/",cmp_name), recursive = TRUE, showWarnings = F)
              
              subData <- dataMat[, which(sample.data[[ion_type]]$Condition %in% c(expGroup, ctrlGroup))]
              
              # 分别提取实验组和对照组的样本名称（在subData中）
              expSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == expGroup]
              ctrlSamples <- sample.data[[ion_type]]$`File name`[sample.data[[ion_type]]$Condition == ctrlGroup]
              
              # 计算每个代谢物的t检验p值（在当前比较样本上）
              pvals <- apply(subData, 1, function(x) {
                
                tryCatch({
                  
                  if (input$pvalue_type=="ttest") {
                    
                    t_res <- t.test(x[expSamples], x[ctrlSamples])
                    
                  }
                  
                  if (input$pvalue_type=="wilcox_test") {
                    
                    t_res <- wilcox.test(x[expSamples], x[ctrlSamples])
                    
                  }
                  
                  return(t_res$p.value)
                  
                },error=function(e){
                  
                  return(1)
                  
                })
                
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
              
              tryCatch({
                
                # 利用ropls包进行PLS-DA，计算VIP值（只取预测成分predI=1）
                selectedSamples <- c(expSamples, ctrlSamples)
                X <- t(subData[, selectedSamples])
                y <- factor(c(rep("exp", length(expSamples)), rep("ctrl", length(ctrlSamples))))
                model <- opls(X, y, predI = 1, orthoI = NA, crossvalI=ifelse(length(selectedSamples)<7,length(selectedSamples),7))
                if (nrow(model@modelDF)==0) {
                  model <- opls(X, y, predI = 1, orthoI = 1, crossvalI=ifelse(length(selectedSamples)<7,length(selectedSamples),7))
                }
                # png(paste0("5.差异筛选/", cmp_name, "/", cmp_name, "_plsda_", ion_type, ".png"))
                # plot(model, typeVc = "x-score")  # 生成预测分类图
                # dev.off()
                # pdf(paste0("5.差异筛选/", cmp_name, "/", cmp_name, "_plsda_", ion_type, ".pdf"))
                # plot(model, typeVc = "x-score")  # 生成预测分类图
                # dev.off()
                vip <- getVipVn(model)
                vip <- vip[rownames(subData)]
                
              },error=function(e){
                
                vip <<- 0
                
              })
              
              # 整合统计结果和当前比较的定量数据
              res_df <- data.frame(peak_name = rownames(subData),
                                   name = metabo.data.raw[[ion_type]]$name[match(rownames(subData), metabo.data.raw[[ion_type]]$peak_name)],
                                   pvalue = pvals,
                                   adj_p_value = p.adjust(pvals, method = ifelse(input$padjust_judge, input$padjust_method, "none")),
                                   FC = fc,
                                   Log2FC = log2(fc),
                                   VIP = vip,
                                   rsd = metabo.data.rsd$qc_rsd,
                                   stringsAsFactors = FALSE)
              
              # 1. 定义需要添加的目标列名列表
              target_cols <- c("formula", "confidence_level", "adduct", "total_score", "mz_error", "rt_error_abs", "rt_error_rela", "ms2_score", "iden_score")
              
              # 2. 获取当前ion_type下的列名
              current_cols <- colnames(metabo.data.raw.all.col[[ion_type]])
              
              # 3. 筛选实际存在的目标列
              cols_to_add <- target_cols[target_cols %in% current_cols]
              
              # 4. 添加匹配的列到res_df
              if (length(cols_to_add) > 0) {
                # 创建包含目标列的临时数据框
                temp_df <- metabo.data.raw.all.col[[ion_type]][
                  match(res_df$peak_name, metabo.data.raw.all.col[[ion_type]]$peak_name),
                  cols_to_add,
                  drop = FALSE  # 确保单列时仍保持数据框结构
                ]
                
                # 将目标列添加到结果数据框
                res_df <- cbind(res_df, temp_df)
                
                # 可选：打印添加的列信息
                message("已成功添加以下列到结果数据框：", paste(cols_to_add, collapse = ", "))
              } else {
                message("未找到目标列：", paste(target_cols, collapse = ", "))
              }
              
              # 增加上下调标志
              res_df[["difference"]] <- "nodiff"
              
              # fc筛选
              if (input$fc_cutoff_judge & !input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
              }
              
              # vip筛选
              if (!input$fc_cutoff_judge & input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
              }
              
              # p筛选
              if (!input$fc_cutoff_judge & !input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge) {
                  res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$adj_p_value < input$p_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$adj_p_value < input$p_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$Log2FC > 0 & res_df$pvalue < input$p_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$Log2FC < 0 & res_df$pvalue < input$p_cutoff)] <- "down"
                }
              }
              
              # fc和p值筛选
              if (input$fc_cutoff_judge & !input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff))] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff))] <- "down"
                }
              }
              
              # vip和p值筛选
              if (!input$fc_cutoff_judge & input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > 0 & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < 0 & res_df$VIP>input$vip_cutoff)] <- "down"
                }
              }
              
              # fc和vip筛选
              if (input$fc_cutoff_judge & input$vip_cutoff_judge & !input$p_cutoff_judge) {
                res_df[["difference"]][which(res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                res_df[["difference"]][which(res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
              }
              
              # fc，vip，p值筛选
              if (input$fc_cutoff_judge & input$vip_cutoff_judge & input$p_cutoff_judge) {
                if (input$padjust_judge){
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$adj_p_value < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
                }else{
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC > log2(input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "up"
                  res_df[["difference"]][which(res_df$pvalue < input$p_cutoff & res_df$Log2FC < log2(1/input$fc_cutoff) & res_df$VIP>input$vip_cutoff)] <- "down"
                }
              }
              
              # 匹配注释
              selected_columns_index <- c("peak_name", "mz", "rt", "id_kegg", "hmdb_kegg", "lmsd_kegg")
              # 去除不存在的列索引（防止有些列缺失）
              selected_columns_index <- selected_columns_index[selected_columns_index %in% colnames(metabo.data.raw[[ion_type]])]
              anno_metabo_data <- metabo.data.raw[[ion_type]][,selected_columns_index]
              anno_metabo_data$peak_name <- as.character(anno_metabo_data$peak_name)
              res_df <- cbind(left_join(res_df, anno_metabo_data), subData)
              
              # 添加Mode标记
              res_df$peak_name <- paste0(res_df$peak_name, "_", ion_type)
              
              if (!input$single_ion_judge) {
                
                colnames(res_df) <- str_replace_all(colnames(res_df), "(\\.POS)|(\\.pos)|(\\.NEG)|(\\.neg)", "")
                
              }
              
              results_mode[[cmp_name]][[ion_type]] <- res_df
              
              results_mode_diff[[cmp_name]][[ion_type]] <- res_df_diff <- res_df[which(res_df$difference != "nodiff"),]
              
              results_mode_diff_name[[cmp_name]][[ion_type]] <- res_df_diff_name  <- res_df[which(res_df$difference != "nodiff" & !is.na(res_df$name)),]
              
              summary_stats_list[[cmp_name]][[ion_type]] <- data.frame(compare=cmp_name,
                                                                       ion_type=ion_type,
                                                                       raw_feature_num = nrow(metabo.data),
                                                                       remove_missing_num = nrow(metabo.data2),
                                                                       preprocess_feature_num=nrow(res_df),
                                                                       preprocess_feature_num_name=sum(!is.na(res_df$name)),
                                                                       diff_feature=0,
                                                                       diff_feature_up=0,
                                                                       diff_feature_down=0,
                                                                       diff_feature_name=0,
                                                                       diff_feature_name_up=0,
                                                                       diff_feature_name_down=0)
              
            }
            
          } else {
            
            summary_stats_list[["no comapre"]][[ion_type]] <- data.frame(compare="no compare",
                                                                         ion_type=ion_type,
                                                                         raw_feature_num = nrow(metabo.data),
                                                                         remove_missing_num = nrow(metabo.data2),
                                                                         preprocess_feature_num=nrow(metabo.data.fill),
                                                                         preprocess_feature_num_name=sum(!is.na(metabo.data.raw[[ion_type]]$name[match(rownames(metabo.data.fill), metabo.data.raw[[ion_type]]$peak_name)])))
            
          }
          
        }
        
        # 保存结果rda数据
        save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw"), file = "result.rda")
        # load(file = "result.rda")
        # ion_type = "NEG"
        # ion_type = "single"
        # sum_pca <- read.xlsx(paste0("4.多元统计分析/1.PCA/PCA_summary.xlsx"))
        # sum_plsda <- read.xlsx(paste0("4.多元统计分析/2.PLS-DA/PLS-DA_summary.xlsx"))
        # sum_oplsda <- read.xlsx(paste0("4.多元统计分析/3.OPLS-DA/OPLS-DA_summary.xlsx"))
        
        if (length(unique(sample.data[[ion_type]][["Condition"]][!grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T)]))>1) {
          
          sum_pca <- do.call(rbind, sum_pca)
          write.xlsx(sum_pca, paste0("4.多元统计分析/1.PCA/PCA_summary.xlsx"))
          
          sum_plsda <- do.call(rbind, sum_plsda)
          write.xlsx(sum_plsda, paste0("4.多元统计分析/2.PLS-DA/PLS-DA_summary.xlsx"))
          
          sum_oplsda <- do.call(rbind, sum_oplsda)
          write.xlsx(sum_oplsda, paste0("4.多元统计分析/3.OPLS-DA/OPLS-DA_summary.xlsx"))
          
          # 各组NEG+POS合并后的生信 ----
          for (cmp_name in group.data) {
            
            experimental <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[1]
            control <- str_split(cmp_name, pattern = "_vs_")%>%unlist%>%.[2]
            
            # 合并当前比较中POS和NEG的结果（行合并，每行为一个代谢物的统计结果）
            if (!input$single_ion_judge) {
              
              for (i in 1:length(results_mode[[cmp_name]])) {
                
                colnames(results_mode[[cmp_name]][[i]]) <- colnames(results_mode[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                colnames(results_mode_diff[[cmp_name]][[i]]) <- colnames(results_mode_diff[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                colnames(results_mode_diff_name[[cmp_name]][[i]]) <- colnames(results_mode_diff_name[[cmp_name]][[2]])%>%str_replace_all("(POS)|(pos)","")
                
              }
              
            }
            
            if (input$single_ion_judge) {
              
              sample.data[["POS"]] <- sample.data[[ion_type]]
              
              sample.data[["POS"]]$`File name2` <- sample.data[["POS"]]$`File name`
              
            }else{
              
              sample.data[["POS"]]$`File name2` <- sample.data[["POS"]]$`File name`%>%str_replace_all("(\\.POS)|(\\.pos)","")%>%str_replace_all("\\.\\.","\\.")
              
            }
            
            expSamples <- sample.data[["POS"]]$`File name2`[sample.data[["POS"]]$Condition == experimental]
            ctrlSamples <- sample.data[["POS"]]$`File name2`[sample.data[["POS"]]$Condition == control]
            sample_cols <- c(expSamples, ctrlSamples)
            
            group <- sample.data[["POS"]]$Condition
            
            # 合并当前比较中POS和NEG的结果（行合并，每行为一个代谢物的统计结果）
            results_mode_merged <- do.call(rbind, results_mode[[cmp_name]])
            
            # 数据的拆分和去冗余
            wb <- createWorkbook()
            addWorksheet(wb, "合并")
            addWorksheet(wb, "拆分")
            addWorksheet(wb, "去冗余")
            addWorksheet(wb, "差异筛选")
            addWorksheet(wb, "差异筛选_name")
            writeData(wb, "合并", results_mode_merged)
            
            # 筛选出数据框中实际存在的列
            if (is.na(input$key_column_kegg)){
              
              existing_cols <- c("name", "formula", "confidence_level", "adduct", "total_score", "mz_error", "rt_error_abs", "rt_error_rela", "ms2_score", "iden_score", "id_kegg") %>% intersect(colnames(results_mode_merged))
              
            }else{
              
              existing_cols <- c("name", "formula", "confidence_level", "adduct", "total_score", "mz_error", "rt_error_abs", "rt_error_rela", "ms2_score", "iden_score") %>% intersect(colnames(results_mode_merged))
              
            }
            
            # 执行拆分操作
            results_mode_merged <- results_mode_merged %>% dplyr::mutate(across(all_of(existing_cols), as.character))
            results_mode_merged <- results_mode_merged %>% separate_rows(all_of(existing_cols), sep = ";")
            writeData(wb, "拆分", results_mode_merged)
            
            results_mode_merged$name2 <- paste0(results_mode_merged$peak_name, "_", results_mode_merged$name)
            results_mode_merged$sum_intensity <- apply(results_mode_merged[, sample_cols], 1, sum)
            results_mode_merged$normalized_name <- tolower(results_mode_merged$name)
            
            # 执行去冗余操作
            if ("total_score" %in% colnames(results_mode_merged)) {
              
              if (sum(grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T))>0) {
                
                results_mode_merged <- rbind(dplyr::arrange(results_mode_merged[which(!is.na(results_mode_merged$name)), ], rsd, desc(total_score), desc(sum_intensity)) %>% dplyr::distinct(normalized_name, .keep_all = T), results_mode_merged[which(is.na(results_mode_merged$name)), ])
                
              }else{
                
                results_mode_merged <- rbind(dplyr::arrange(results_mode_merged[which(!is.na(results_mode_merged$name)), ], desc(total_score), desc(sum_intensity)) %>% dplyr::distinct(normalized_name, .keep_all = T), results_mode_merged[which(is.na(results_mode_merged$name)), ])
                
              }
              
            }else{
              
              if (sum(grepl("^qc$", trimws(sample.data[[ion_type]][["Condition"]]), ignore.case = T))>0) {
                
                results_mode_merged <- rbind(dplyr::arrange(results_mode_merged[which(!is.na(results_mode_merged$name)), ], rsd, desc(sum_intensity)) %>% dplyr::distinct(normalized_name, .keep_all = T), results_mode_merged[which(is.na(results_mode_merged$name)), ])
                
              }else{
                
                results_mode_merged <- rbind(dplyr::arrange(results_mode_merged[which(!is.na(results_mode_merged$name)), ], desc(sum_intensity)) %>% dplyr::distinct(normalized_name, .keep_all = T), results_mode_merged[which(is.na(results_mode_merged$name)), ])
                
              }
              
            }
            
            results_mode_merged <- results_mode_merged[, c(colnames(results_mode_merged)[which(!colnames(results_mode_merged) %in% sample_cols)], sample_cols)]
            writeData(wb, "去冗余", results_mode_merged)
            
            results_mode_diff_merged <- results_mode_merged[which(results_mode_merged$difference != "nodiff"),]
            writeData(wb, "差异筛选", results_mode_diff_merged)
            
            results_mode_diff_name_merged <- results_mode_merged[which(results_mode_merged$difference != "nodiff" & !is.na(results_mode_merged$name)),]
            writeData(wb, "差异筛选_name", results_mode_diff_name_merged)
            
            saveWorkbook(wb, paste0("5.差异筛选/", cmp_name, "/", cmp_name, "_表格拆分与去重.xlsx"), overwrite = TRUE)
            
            ## 差异的代谢物（feature）---------------------------
            if (nrow(results_mode_diff_merged)>0) {
              
              expr_data <- results_mode_diff_merged[, sample_cols]
              rownames(expr_data) <- results_mode_diff_merged$name2
              
              ### 1. HCA（层级聚类热图和表格） ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/1.HCA")
              dir.create(output_folder, recursive = T)
              
              # 对行（feature）进行标准化（Z-score），便于聚类比较
              # expr_data_scaled <- expr_data
              expr_data_scaled <- t(scale(t(expr_data)))
              
              hca_group <- c(group[which(group == experimental)], group[which(group == control)])
              annotation_col <- data.frame(Group = hca_group)
              rownames(annotation_col) <- colnames(expr_data_scaled) # 确保样本名对齐
              
              # 构建注释颜色列表
              annotation_colors <- list(
                Group = dynamic_colors[which(names(dynamic_colors) %in% c(experimental, control))]  # 对应annotation_col的分组列名
              )
              
              # 获取数据范围（假设已标准化，范围对称）
              data_range <- range(expr_data_scaled)
              max_abs <- min(data_range[2],3,na.rm = T)
              min_abs <- max(data_range[1],-3,na.rm = T)
              
              # 定义对称的数值断点（覆盖数据范围，0居中）
              bk <- seq(min_abs, max_abs, length.out = 100)
              
              # 定义颜色渐变（100级）
              heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
              
              # 绘制热图
              pheatmap::pheatmap(expr_data_scaled, 
                                 # scale = "row",
                                 cluster_cols = input$hca_col_cluster,
                                 color = heatmap_colors,
                                 breaks = bk,
                                 clustering_distance_rows = "euclidean", 
                                 clustering_method = "complete", 
                                 main = "Heatmap of Differential Features",
                                 filename = file.path(output_folder, "heatmap.pdf"),
                                 height = max(0.1*nrow(expr_data),15),
                                 show_rownames=F, show_colnames=F,
                                 annotation_col = annotation_col, # 添加列分组注释
                                 annotation_names_col = TRUE,
                                 annotation_colors = annotation_colors)
              
              # 计算层级聚类（行聚类）
              hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
              # 保存聚类顺序表格
              cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
              write.xlsx(cluster_order, file = file.path(output_folder, "HCA.xlsx"), rowNames = T)
              
              ### 2. Correlation Analysis ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/2.correlation")
              dir.create(output_folder, recursive = T)
              
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
              
              if (nrow(cor_mat)>20) {
                
                # 步骤1：计算特征的综合相关性评分（取绝对值均值）
                feature_scores <- rowMeans(abs(cor_mat), na.rm = TRUE)
                
                # 步骤2：按评分降序排列，提取前20个特征名
                top_features <- names(sort(feature_scores, decreasing = TRUE)[1:20])
                
                # 步骤3：从原始矩阵中筛选出这20个特征的相关性子矩阵
                cor_mat_top <- cor_mat[top_features, top_features]
                
                # 步骤4：绘制筛选后的热图（调整标签可读性）
                pdf(file.path(output_folder, "correlation_top20_features.pdf"))
                corrplot(cor_mat_top,
                         method = "circle",
                         type = "upper",
                         tl.col = "black",
                         tl.cex = 0.6,      # 增大标签字号方便查看
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Top 20 Features by Mean Absolute Correlation",
                         mar = c(0,0,1,0),
                         order = "hclust",  # 按层次聚类排序[2,3](@ref)
                         addrect = 2)        # 添加聚类分界线[2,3](@ref)
                dev.off()
                
              }else{
                
                # 绘制相关性热图
                pdf(file.path(output_folder, "correlation_features.pdf"))
                corrplot(cor_mat, method = "circle", tl.col = "black", tl.cex = 0.6,
                         type="upper",
                         # addCoef.col = "black",
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Correlation Among Features", mar = c(0,0,1,0))
                dev.off()
                
                cor_mat_top <- cor_mat
                
              }
              
              # 保存相关性矩阵
              write.xlsx(cor_mat_top, file = file.path(output_folder, "correlation_features.xlsx"), rowNames = TRUE)
              
              ### 3. Volcano Plot ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/3.volcano")
              dir.create(output_folder, recursive = T)
              
              # 计算 log2(FC)
              diff_features_all <- results_mode_merged %>% dplyr::mutate(negLog10P = -log10(pvalue))
              
              # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
              max_val <- ceiling(max(abs(diff_features_all$Log2FC[is.finite(diff_features_all$Log2FC)]), na.rm = TRUE))
              
              if (input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference, size=VIP)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (input$p_cutoff_judge & !input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (!input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = VIP, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = input$vip_cutoff, 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "VIP") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              ggsave(filename = file.path(output_folder, "volcano.png"), volcano_plot, width = 7, height = 6)
              ggsave(filename = file.path(output_folder, "volcano.pdf"), volcano_plot, width = 7, height = 6)
              write.xlsx(diff_features_all, file.path(output_folder, "volcano.xlsx"))
              
              ### 4. VIP分析 ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/4.VIP_analysis")
              dir.create(output_folder, recursive = T)
              
              # 数据预处理（确保代谢物名称按VIP降序排列）
              vip_bubble_data <- results_mode_diff_merged %>% dplyr::arrange(desc(VIP))
              
              if (nrow(vip_bubble_data)>20) {
                
                vip_bubble_data2 <- vip_bubble_data[1:20,]
                
              }else{
                
                vip_bubble_data2 <- vip_bubble_data
                
              }
              
              vip_bubble_data2$peak_name <- factor(vip_bubble_data2$peak_name, levels = rev(unique(vip_bubble_data2$peak_name)))
              
              # 气泡图绘制
              vip_bubble <- ggplot(vip_bubble_data2, aes(x = VIP, y = peak_name, color = VIP, size = VIP)) + 
                scale_color_gradient(
                  low = "#56B1F7",  # 深蓝
                  high = "#132B43",  # 亮蓝
                  name = "VIP_Score",
                  space = "Lab"  # 使用Lab色彩空间增强渐变平滑度[6](@ref)
                )+
                geom_point() +
                # 坐标轴与主题
                labs(x = "VIP Value", y = "Metabolites") +
                theme_bw() +
                theme(
                  axis.text.y = element_text(size = 8),  # 代谢物名称字号
                  legend.position = "right"
                )
              
              ggsave(file.path(output_folder, "VIP_distribution.png"), vip_bubble, width = 8, height = 6)
              ggsave(file.path(output_folder, "VIP_distribution.pdf"), vip_bubble, width = 8, height = 6)
              
              write.xlsx(vip_bubble_data, file.path(output_folder, "VIP_feature.xlsx"), rowNames = FALSE)
              
              ### 6.差异网络分析 ---------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/6.差异网络分析")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              if(nrow(expr_data)<5000){
                
                # 计算代谢物间的相关系数矩阵
                cor_mat2 <- cor_mat <- cor(t(expr_data), method = "spearman")  # 使用Spearman相关系数
                
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
                  dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
                
                cor_df2 <- as.data.frame(as.table(cor_mat2)) %>% 
                  dplyr::rename(Metabolite1 = Var1, Metabolite2 = Var2, Correlation = Freq) %>%
                  dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
                
                diff_features <- results_mode_diff_merged %>% dplyr::mutate(regulation = difference)
                
                diff_features$regulation <- factor(diff_features$regulation, levels = c("up", "down"))
                
                diff_features$truncated_name <- truncated_names
                
                if (nrow(diff_features)>50) {
                  
                  # 按VIP排序并取前50
                  top_features <- diff_features %>%
                    dplyr::arrange(desc(VIP)) %>%
                    head(50)
                  
                  # 筛选cor_df只包含这50个代谢物的边
                  cor_df <- cor_df %>% dplyr::filter(Metabolite1 %in% top_features$truncated_name & Metabolite2 %in% top_features$truncated_name)
                  
                  # 更新diff_features为筛选后的50个
                  diff_features <- top_features
                  
                }
                
                # 创建igraph图对象
                g <- graph_from_data_frame(
                  d = cor_df[, 1:3],  # 边数据（Metabolite1, Metabolite2, Correlation）
                  directed = T,
                  vertices = diff_features %>% dplyr::select(truncated_name, Log2FC, regulation)  # 节点属性
                )
                
                # 删除孤立节点
                deg <- igraph::degree(g, mode = "all")
                g <- delete_vertices(g, which(deg == 0))
                
                # 节点属性设置
                V(g)$size <- abs(V(g)$Log2FC) * 3  # 大小反映Log2FC绝对值
                V(g)$color <- ifelse(V(g)$regulation == "up", "#E64B35FF", "#3C5488FF")  # 上调红，下调蓝
                
                # 边属性设置
                E(g)$width <- abs(E(g)$Correlation) * 0.2  # 边宽与相关性强度正相关
                E(g)$color <- ifelse(E(g)$Correlation > 0, "#E64B35FF", "#3C5488FF")  # 正相关红，负相关蓝
                
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
                        size = abs(Log2FC)),   # 大小映射到Log2FC绝对值
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
                    values = c("Positive" = "#E64B35FF", "Negative" = "#3C5488FF")
                  ) +
                  # 边宽度映射（连续型）
                  scale_edge_width_continuous(
                    name = "|Correlation|", 
                    range = c(0.5, 4),  # 边宽度范围
                    breaks = c(0.1, 0.4, 0.7)  # 图例刻度
                  ) +
                  # 节点颜色映射（离散型）
                  scale_color_manual(
                    name = "Regulation",
                    values = c("up" = "#E64B35FF", "down" = "#3C5488FF")
                  ) +
                  # 节点大小映射（连续型）
                  scale_size_continuous(
                    name = "|Log2FC|", 
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
                write.graph(g, file.path(output_folder, "feature_network.graphml"), format = "graphml")
                
                # 正确生成节点属性表的方法
                node_table <- diff_features %>%
                  # 首先确保只选择实际存在的列
                  dplyr::select(
                    name2,  # 截断后的代谢物名称（必须存在于diff_features中）
                    Log2FC,         # 差异倍数
                    VIP,            # VIP值
                    regulation      # 上下调标记
                  ) %>%
                  # 添加原始名称（如果feature_names是独立向量）
                  dplyr::mutate(
                    Original_Name = rownames(diff_features),  # 假设原始名称是行名
                    Direction = ifelse(regulation == "up", "Upregulated", "Downregulated"),
                    Metabolite = name2  # 重命名列
                  ) %>%
                  # 重新排列列顺序
                  dplyr::select(
                    Metabolite,
                    Original_Name,
                    Log2FC,
                    VIP,
                    Regulation = regulation,
                    Direction
                  )
                
                # 生成边属性表（已筛选相关性强的边）
                edge_table <- cor_df %>%
                  dplyr::mutate(
                    Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative"),
                    Abs_Correlation = abs(Correlation)  # 相关性绝对值
                  ) %>%
                  dplyr::select(
                    Source = Metabolite1,        # 源节点
                    Target = Metabolite2,        # 目标节点
                    Correlation,                 # 相关系数
                    Correlation_Type,            # 正/负相关
                    Abs_Correlation              # 相关性强度
                  )
                
                # 计算网络拓扑指标
                node_topology <- data.frame(
                  Metabolite = V(g)$name,
                  Degree = igraph::degree(g),            # 连接度
                  Betweenness = igraph::betweenness(g),  # 中介中心性
                  Closeness = igraph::closeness(g)       # 接近中心性
                )
                
                # 合并节点属性与拓扑指标
                network_summary <- node_table %>%
                  dplyr::left_join(node_topology, by = c("Metabolite" = "Metabolite"))
                
                # 创建新的Excel工作簿
                wb <- createWorkbook()
                
                # 1. 添加节点属性表（Sheet1）
                addWorksheet(wb, sheetName = "节点属性")
                writeData(wb, sheet = "节点属性", x = node_table, startCol = 1, startRow = 1)
                
                # 2. 添加边属性表（Sheet2）
                addWorksheet(wb, sheetName = "边属性")
                writeData(wb, sheet = "边属性", x = edge_table, startCol = 1, startRow = 1)
                
                # 3. 添加网络拓扑表（Sheet3，可选）
                addWorksheet(wb, sheetName = "网络拓扑")
                writeData(wb, sheet = "网络拓扑", x = network_summary, startCol = 1, startRow = 1)
                
                # 保存Excel文件
                saveWorkbook(wb, file = file.path(output_folder, "feature_network.xlsx"), overwrite = TRUE)
                
              }
              
            }
            
            ## 差异的代谢物（有name的）---------------------------
            if (nrow(results_mode_diff_name_merged)>0) {
              
              expr_data <- results_mode_diff_name_merged[,sample_cols]
              rownames(expr_data) <- paste0(results_mode_diff_name_merged$name, "_", results_mode_diff_name_merged$peak_name)
              
              ### 1. HCA（层级聚类热图和表格） ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/1.HCA")
              
              # 对行（feature）进行标准化（Z-score），便于聚类比较
              # expr_data_scaled <- expr_data
              expr_data_scaled <- t(scale(t(expr_data)))
              
              # 绘制热图
              truncated_labels <- substr(rownames(expr_data_scaled), 1, pmin(nchar(rownames(expr_data_scaled)), 40)) # 最多保留40个字符
              # 添加省略号表示截断
              truncated_labels <- ifelse(
                nchar(rownames(expr_data_scaled)) > 20,
                paste0(truncated_labels, "..."),
                truncated_labels
              )
              
              # 获取数据范围（假设已标准化，范围对称）
              data_range <- range(expr_data_scaled)
              max_abs <- min(data_range[2],3,na.rm = T)
              min_abs <- max(data_range[1],-3,na.rm = T)
              
              # 定义对称的数值断点（覆盖数据范围，0居中）
              bk <- seq(min_abs, max_abs, length.out = 100)
              
              # 定义颜色渐变（100级）
              heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
              
              pheatmap::pheatmap(
                expr_data_scaled,
                # scale = "row",
                cluster_cols = input$hca_col_cluster,
                color = heatmap_colors,
                breaks = bk,
                labels_row = truncated_labels,  # 使用截断后的标签
                show_colnames = F,
                clustering_distance_rows = "euclidean",
                clustering_method = "complete",
                main = "Heatmap of Differential Features",
                height = max(6, nrow(expr_data_scaled)*0.14),
                width = max(8, ncol(expr_data_scaled)*0.25),
                filename = file.path(output_folder, "heatmap_name.pdf"),
                annotation_col = annotation_col, # 添加列分组注释
                annotation_names_col = TRUE,
                annotation_colors = annotation_colors)
              
              # 计算层级聚类（行聚类）
              hc_rows <- hclust(dist(expr_data_scaled), method = "complete")
              # 保存聚类顺序表格
              cluster_order <- data.frame(expr_data_scaled[hc_rows$order,])
              write.xlsx(cluster_order, file = file.path(output_folder, "heatmap_name.xlsx"), rowNames = T)
              
              ### 2. Correlation Analysis ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/2.correlation")
              
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
              
              if (ncol(cor_mat)>20) {
                
                # 步骤1：计算特征的综合相关性评分（取绝对值均值）
                feature_scores <- rowMeans(abs(cor_mat), na.rm = TRUE)
                
                # 步骤2：按评分降序排列，提取前20个特征名
                top_features <- names(sort(feature_scores, decreasing = TRUE)[1:20])
                
                # 步骤3：从原始矩阵中筛选出这20个特征的相关性子矩阵
                cor_mat_top <- cor_mat[top_features, top_features]
                
                # 步骤4：绘制筛选后的热图（调整标签可读性）
                pdf(file.path(output_folder, "correlation_top20_features_name.pdf"))
                corrplot(cor_mat_top,
                         method = "circle",
                         type = "upper",
                         tl.col = "black",
                         tl.cex = 0.6,      # 增大标签字号方便查看
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Top 20 Features by Mean Absolute Correlation",
                         mar = c(0,0,1,0),
                         order = "hclust",  # 按层次聚类排序[2,3](@ref)
                         addrect = 2)        # 添加聚类分界线[2,3](@ref)
                dev.off()
                
              }else{
                
                # 绘制相关性热图
                pdf(file.path(output_folder, "correlation_features_name.pdf"))
                corrplot(cor_mat, method = "circle", tl.col = "black", type = "upper",
                         col = colorRampPalette(c("blue", "white", "red"))(200),
                         title = "Correlation Among Samples", mar = c(0,0,1,0))
                dev.off()
                
              }
              
              # 保存相关性矩阵
              write.xlsx(t(expr_data), file = file.path(output_folder, "correlation_features_name.xlsx"), rowNames = TRUE)
              
              ### 3. Volcano Plot ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/3.volcano")
              
              # 计算 log2(FC)
              diff_features_all <- results_mode_merged[which(!is.na(results_mode_merged$name)),] %>% dplyr::mutate(negLog10P = -log10(pvalue))
              
              # 计算 log2(FC) 的最大绝对值（用于设置对称的横坐标范围）
              max_val <- ceiling(max(abs(diff_features_all$Log2FC[is.finite(diff_features_all$Log2FC)]), na.rm = TRUE))
              
              if (input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference, size=VIP)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (input$p_cutoff_judge & !input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = negLog10P, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = -log10(input$p_cutoff), 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "-log10(p-value)") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              if (!input$p_cutoff_judge & input$vip_cutoff_judge) {
                
                volcano_plot <- ggplot(diff_features_all, aes(x = Log2FC, y = VIP, color = difference)) + 
                  geom_point(alpha = 0.6) +
                  # 颜色映射设置（红蓝灰）
                  scale_color_manual(
                    values = c(up = "red", down = "blue", nodiff = "grey50"),
                    labels = c(up = "Up", down = "Down", nodiff = "nodiff"),
                    name = "Status"
                  ) +
                  # 参考线保持原样
                  geom_vline(xintercept = c(log2(input$fc_cutoff), -log2(input$fc_cutoff)),
                             linetype = "dashed", color = "black") +
                  geom_hline(yintercept = input$vip_cutoff, 
                             linetype = "dashed", color = "black") +
                  theme_minimal() +
                  labs(x = "log2(Fold Change)", y = "VIP") +
                  scale_x_continuous(limits = c(-max_val, max_val)) +
                  # 图例排版优化
                  theme(
                    legend.position = "right",
                    legend.box = "vertical",
                    legend.spacing.y = unit(0.2, "cm")
                  )
                
              }
              
              ggsave(filename = file.path(output_folder, "volcano_name.png"), volcano_plot, width = 7, height = 6)
              ggsave(filename = file.path(output_folder, "volcano_name.pdf"), volcano_plot, width = 7, height = 6)
              write.xlsx(diff_features_all, file.path(output_folder, "volcano_name.xlsx"))
              
              ### 4. VIP分析 ---------------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/4.VIP_analysis")
              dir.create(output_folder, recursive = T)
              
              # 数据预处理（确保代谢物名称按VIP降序排列）
              vip_bubble_data <- results_mode_diff_merged %>% dplyr::arrange(desc(VIP)) %>% dplyr::filter(!is.na(name))
              
              if (nrow(vip_bubble_data)>20) {
                
                vip_bubble_data2 <- vip_bubble_data[1:20,]
                
              }else{
                
                vip_bubble_data2 <- vip_bubble_data
                
              }
              
              vip_bubble_data2$name <- factor(vip_bubble_data2$name, levels = rev(unique(vip_bubble_data2$name)))
              
              # 创建标签截断函数
              truncate_labels <- function(x) {
                
                ifelse(nchar(x) > 30, paste0(substr(x, 1, 30), "..."), x)
                
              }
              
              # 气泡图绘制
              vip_bubble <- ggplot(vip_bubble_data2, aes(x = VIP, y = name, color = VIP, size = VIP)) + 
                scale_color_gradient(
                  low = "#56B1F7",  # 深蓝
                  high = "#132B43",  # 亮蓝
                  name = "VIP_Score",
                  space = "Lab"  # 使用Lab色彩空间增强渐变平滑度[6](@ref)
                )+
                geom_point() +
                scale_y_discrete(
                  labels = truncate_labels,  # 应用截断函数[7,8](@ref)
                  expand = expansion(add = 0.2)  # 增加标签间距避免截断重叠
                )+
                # 坐标轴与主题
                labs(x = "VIP Value", y = "Metabolites") +
                theme_bw() +
                theme(
                  axis.text.y = element_text(size = 8),  # 代谢物名称字号
                  legend.position = "right"
                )
              
              ggsave(file.path(output_folder, "VIP_distribution_metabolite.png"), vip_bubble, width = 8, height = 6)
              ggsave(file.path(output_folder, "VIP_distribution_metabolite.pdf"), vip_bubble, width = 8, height = 6)
              
              write.xlsx(vip_bubble_data, file.path(output_folder, "VIP_metabolite.xlsx"), rowNames = FALSE)
              
              ### 5.差异代谢物表达量的箱线图 ---------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/5.差异代谢物表达量")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 转换为长数据格式（代谢物ID + 样本 + 定量值 + 分组）
              if (length(unique(results_mode_diff_name_merged$name)) < 50){
                
                expr_long <- results_mode_diff_name_merged[,c("name", sample_cols)] %>%
                  tidyr::pivot_longer(cols = -name, names_to = "sample", values_to = "intensity") %>%
                  dplyr::left_join(sample.data$POS, by = c("sample" = "File name2")) %>%
                  dplyr::mutate(group = factor(Condition, levels = unique(sample.data$POS$Condition)))
                
              }else{
                
                expr_long <- results_mode_diff_name_merged%>%dplyr::arrange(desc(VIP))%>%.[1:50, c("name", sample_cols)] %>%
                  tidyr::pivot_longer(cols = -name, names_to = "sample", values_to = "intensity") %>%
                  dplyr::left_join(sample.data$POS, by = c("sample" = "File name2")) %>%
                  dplyr::mutate(group = factor(Condition, levels = unique(sample.data$POS$Condition)))
                
              }
              
              # 所有代谢物全局分布
              if (input$pvalue_type=="ttest") {
                
                stat_compare_means_method <- "t.test"
                
              }else if(input$pvalue_type=="wilcox_test") {
                
                stat_compare_means_method <- "wilcox.test"
                
              }
              
              truncate_and_force_wrap <- function(x, max_length = 30, wrap_width = 15) {
                long_idx <- str_length(x) > max_length
                
                # 超长名字先截断加省略号
                x[long_idx] <- paste0(str_sub(x[long_idx], 1, max_length), "...")
                
                # 超长名字才强制换行
                x[long_idx] <- sapply(x[long_idx], function(s) {
                  chars <- strsplit(s, "")[[1]]
                  paste(sapply(seq(1, length(chars), by = wrap_width),
                               function(i) paste0(chars[i:min(i+wrap_width-1, length(chars))], collapse = "")),
                        collapse = "\n")
                })
                
                return(x)
              }
              
              p_global <- ggboxplot(expr_long, x = "group", y = "intensity",color = "group", palette = "npg")+
                stat_compare_means(method = stat_compare_means_method,label = "p.format", vjust = 1, hjust=-0.5)+
                labs(title = "Metabolites Comparison", x = NULL, y = "Normalized Intensity") +
                scale_color_manual(values = dynamic_colors) + 
                # theme_classic() +
                # 按代谢物分面
                facet_wrap(~name, scales = "free", ncol = 6, labeller = labeller(name = function(x) truncate_and_force_wrap(x, max_length = 30, wrap_width = 15)))+
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
                dpi = 300,
                limitsize = FALSE
              )
              
              ggsave(
                file.path(output_folder, "Global_Metabolite_Comparison.png"), 
                p_global,
                width = 16, 
                height = ceiling(length(unique(expr_long$name))/4) * 3,  # 动态调整高度
                dpi = 300,
                limitsize = FALSE
              )
              
              ### 6.差异网络分析 ---------------
              output_folder <- paste0("6.差异代谢物可视化/", cmp_name, "/6.差异网络分析")
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 计算代谢物间的相关系数矩阵
              cor_mat2 <- cor_mat <- cor(t(expr_data), method = "spearman")  # 使用Spearman相关系数
              
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
              
              cor_df2 <- as.data.frame(as.table(cor_mat2)) %>% 
                dplyr::rename(Metabolite1 = Var1, Metabolite2 = Var2, Correlation = Freq) %>%
                dplyr::filter(Metabolite1 != Metabolite2, abs(Correlation) > input$correlation.cutoff) # 去除自相关
              
              diff_features <- results_mode_diff_name_merged %>% dplyr::mutate(regulation = difference)
              
              diff_features$regulation <- factor(diff_features$regulation, levels = c("up", "down"))
              
              diff_features$truncated_name <- truncated_names
              
              if (nrow(diff_features)>50) {
                
                # 按VIP排序并取前50
                top_features <- diff_features %>%
                  dplyr::arrange(desc(VIP)) %>%
                  head(50)
                
                # 筛选cor_df只包含这50个代谢物的边
                cor_df <- cor_df %>%
                  dplyr::filter(Metabolite1 %in% top_features$truncated_name & 
                                  Metabolite2 %in% top_features$truncated_name)
                
                # 更新diff_features为筛选后的50个
                diff_features <- top_features
                
              }
              
              # 创建igraph图对象
              g <- graph_from_data_frame(
                d = cor_df[, 1:3],  # 边数据（Metabolite1, Metabolite2, Correlation）
                directed = T,
                vertices = diff_features %>% dplyr::select(truncated_name, Log2FC, regulation)  # 节点属性
              )
              
              # 删除孤立节点
              deg <- igraph::degree(g, mode = "all")
              g <- delete_vertices(g, which(deg == 0))
              
              # 节点属性设置
              V(g)$size <- abs(V(g)$Log2FC) * 3  # 大小反映Log2FC绝对值
              V(g)$color <- ifelse(V(g)$regulation == "up", "#E64B35FF", "#3C5488FF")  # 上调红，下调蓝
              
              # 边属性设置
              E(g)$width <- abs(E(g)$Correlation) * 0.2  # 边宽与相关性强度正相关
              E(g)$color <- ifelse(E(g)$Correlation > 0, "#E64B35FF", "#3C5488FF")  # 正相关红，负相关蓝
              
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
                      size = abs(Log2FC)),   # 大小映射到Log2FC绝对值
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
                  values = c("Positive" = "#E64B35FF", "Negative" = "#3C5488FF")
                ) +
                # 边宽度映射（连续型）
                scale_edge_width_continuous(
                  name = "|Correlation|", 
                  range = c(0.5, 4),  # 边宽度范围
                  breaks = c(0.1, 0.4, 0.7)  # 图例刻度
                ) +
                # 节点颜色映射（离散型）
                scale_color_manual(
                  name = "Regulation",
                  values = c("up" = "#E64B35FF", "down" = "#3C5488FF")
                ) +
                # 节点大小映射（连续型）
                scale_size_continuous(
                  name = "|Log2FC|", 
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
              
              ggsave(file.path(output_folder, "Network_CytoscapeStyle_name.pdf"), network_plot, width = 12, height = 10)
              ggsave(file.path(output_folder, "Network_CytoscapeStyle_name.png"), network_plot, width = 12, height = 10, dpi = 300)
              
              # 导出网络数据供Cytoscape进一步编辑
              write.graph(g, file.path(output_folder, "feature_network_name.graphml"), format = "graphml")
              
              # 正确生成节点属性表的方法
              node_table <- diff_features %>%
                # 首先确保只选择实际存在的列
                dplyr::select(
                  name2,  # 截断后的代谢物名称（必须存在于diff_features中）
                  Log2FC,         # 差异倍数
                  VIP,            # VIP值
                  regulation      # 上下调标记
                ) %>%
                # 添加原始名称（如果feature_names是独立向量）
                dplyr::mutate(
                  Original_Name = rownames(diff_features),  # 假设原始名称是行名
                  Direction = ifelse(regulation == "up", "Upregulated", "Downregulated"),
                  Metabolite = name2  # 重命名列
                ) %>%
                # 重新排列列顺序
                dplyr::select(
                  Metabolite,
                  Original_Name,
                  Log2FC,
                  VIP,
                  Regulation = regulation,
                  Direction
                )
              
              # 生成边属性表（已筛选相关性强的边）
              edge_table <- cor_df %>%
                dplyr::mutate(
                  Correlation_Type = ifelse(Correlation > 0, "Positive", "Negative"),
                  Abs_Correlation = abs(Correlation)  # 相关性绝对值
                ) %>%
                dplyr::select(
                  Source = Metabolite1,        # 源节点
                  Target = Metabolite2,        # 目标节点
                  Correlation,                 # 相关系数
                  Correlation_Type,            # 正/负相关
                  Abs_Correlation              # 相关性强度
                )
              
              # 计算网络拓扑指标
              node_topology <- data.frame(
                Metabolite = V(g)$name,
                Degree = igraph::degree(g),            # 连接度
                Betweenness = igraph::betweenness(g),  # 中介中心性
                Closeness = igraph::closeness(g)       # 接近中心性
              )
              
              # 合并节点属性与拓扑指标
              network_summary <- node_table %>%
                dplyr::left_join(node_topology, by = c("Metabolite" = "Metabolite"))
              
              # 创建新的Excel工作簿
              wb <- createWorkbook()
              
              # 1. 添加节点属性表（Sheet1）
              addWorksheet(wb, sheetName = "节点属性")
              writeData(wb, sheet = "节点属性", x = node_table, startCol = 1, startRow = 1)
              
              # 2. 添加边属性表（Sheet2）
              addWorksheet(wb, sheetName = "边属性")
              writeData(wb, sheet = "边属性", x = edge_table, startCol = 1, startRow = 1)
              
              # 3. 添加网络拓扑表（Sheet3，可选）
              addWorksheet(wb, sheetName = "网络拓扑")
              writeData(wb, sheet = "网络拓扑", x = network_summary, startCol = 1, startRow = 1)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder, "feature_network_name.xlsx"), overwrite = TRUE)
              
              ### 7.KEGG 富集 ---------------
              output_folder <- paste0("7.KEGG富集/", cmp_name)
              dir.create(paste0(output_folder, "/KEGG_pathway"), recursive = T, showWarnings = F)
              
              # 一个通用的小工具：拆分 + 去噪
              .split_ids <- function(x) {
                if (is.null(x) || all(is.na(x))) return(character(0))
                x <- ifelse(is.na(x), "", x)
                # 先按 ; 再按 /
                out <- strsplit(x, ";") |> unlist(use.names = FALSE)
                out <- strsplit(out, "/") |> unlist(use.names = FALSE)
                out <- trimws(out)
                out[!(is.na(out) | out == "" | out == "NA")]
              }
              
              # 把若干列（可能缺失）合并为去重 id_kegg 字符串
              .combine_kegg_cols <- function(df, cols) {
                if (length(cols) == 0) {
                  # 没有任何可用来源，直接返回原 df 的 id_kegg 或 NA
                  return(df$id_kegg)
                }
                purrr::pmap_chr(df[, cols, drop = FALSE], function(...) {
                  vals <- list(...)
                  ids <- unlist(lapply(vals, .split_ids), use.names = FALSE)
                  ids <- unique(ids)
                  if (length(ids) == 0) NA_character_ else paste(ids, collapse = ";")
                })
              }
              
              if (is.na(input$key_column_kegg)) {
                # 没有原始 KEGG 列：尝试把（名称匹配得到的）id_kegg + hmdb_kegg + lmsd_kegg 合并去重
                # 针对差异特征
                cols_df  <- intersect(c("id_kegg", "hmdb_kegg", "lmsd_kegg"),
                                      names(results_mode_diff_name_merged))
                diff_features <- results_mode_diff_name_merged %>%
                  dplyr::mutate(id_kegg = .combine_kegg_cols(cur_data_all(), cols_df))
                
                # 针对全部特征
                cols_all <- intersect(c("id_kegg", "hmdb_kegg", "lmsd_kegg"),
                                      names(results_mode_merged))
                all_features <- results_mode_merged %>%
                  dplyr::mutate(id_kegg = .combine_kegg_cols(cur_data_all(), cols_all))
                # 展开到逐 ID 形式并去 NA
                all_features <- all_features %>%
                  tidyr::separate_rows(id_kegg, sep = ";") %>%
                  dplyr::filter(!is.na(id_kegg), id_kegg != "")
                
              }else{
                
                diff_features <- results_mode_diff_name_merged
                all_features <- results_mode_merged%>%separate_rows(id_kegg)
                
              }
              
              tryCatch({
                
                # 提取差异代谢物KEGG信息
                kegg_data <- diff_features %>% separate_rows(id_kegg)%>%
                  dplyr::filter(!is.na(name)) %>%
                  dplyr::distinct(id_kegg, .keep_all = TRUE) %>%
                  dplyr::filter(!is.na(id_kegg))%>%.[which(.[["id_kegg"]]!="NA"),]
                
                mb3 <- enrichMBKEGG(kegg_data$id_kegg)
                enrich_result <- mb3@result%>%dplyr::arrange(desc(FoldEnrichment))
                enrich_result$Description <- factor(enrich_result$Description, levels = enrich_result$Description)
                enrich_result <- dplyr::rename(enrich_result, MetaboRatio=GeneRatio, metaboID=geneID)
                
                ## kegg通路制作
                # 从差异分析结果提取代谢物ID和Log2FC
                kegg_fc_data <- diff_features %>% 
                  dplyr::select(id_kegg, Log2FC, difference) %>%
                  dplyr::filter(!is.na(id_kegg)) %>%
                  dplyr::distinct(id_kegg, .keep_all = TRUE)%>%
                  tidyr::separate_rows(id_kegg)%>%.[which(.[["id_kegg"]]!="NA"),]
                
                setwd(paste0(output_folder, "/KEGG_pathway"))
                
                # 更新颜色
                update_node_colors <- function(cpd_data, metabo_data) {
                  
                  cpd_data %>%
                    mutate(
                      # 1. 分割多值ID
                      mapped_ids = ifelse(all.mapped == "", 
                                          list(character(0)), 
                                          str_split(all.mapped, ",")),
                      
                      # 2. 获取每个节点的状态向量
                      status_vec = map(mapped_ids, ~ {
                        if (length(.x) == 0) return(character(0))
                        # 关键修正：匹配difference而非Log2FC
                        matched_status <- metabo_data$difference[metabo_data$id_kegg %in% .x]
                        if (length(matched_status) > 0) matched_status else character(0)
                      }),
                      
                      # 3. 应用着色规则
                      mol.col = map_chr(status_vec, function(states) {
                        if (length(states) == 0) return("#FFFFFF")  # 无匹配保持白色
                        
                        # 核心着色逻辑（基于difference状态）
                        if (any(states == "up") & !any(states == "down")) {
                          "#FF0000"  # 红：仅含up
                        } else if (any(states == "down") & !any(states == "up")) {
                          "#0000FF"  # 蓝：仅含down
                        } else if (any(states == "up") & any(states == "down")) {
                          "#FFA500"  # 橙：混合状态
                        }
                      })
                    ) %>%
                    select(-mapped_ids, -status_vec)
                  
                }
                
                enrich_result$match <- "FALSE"
                
                # 目标匹配数量
                target_n <- nrow(enrich_result)
                matched_n <- 0
                
                # 遍历所有富集到的通路ID
                for (pathway_id in enrich_result$ID) {
                  
                  if (matched_n >= target_n) break
                  
                  # 提取当前通路的代谢物ID和Log2FC
                  pathway_metabo <- kegg_fc_data %>% dplyr::filter(id_kegg %in% (str_split(enrich_result$metaboID[which(enrich_result$ID==pathway_id)], "/", simplify = T)%>%as.character()))
                  
                  tryCatch({
                    
                    # 生成通路图（自动下载通路PNG和XML文件）
                    pathview.data <- pathview(
                      cpd.data = setNames(pathway_metabo$Log2FC, pathway_metabo$id_kegg),
                      pathway.id = gsub("map", "", pathway_id), # 去除前缀
                      species = input$species, # 物种缩写（如人类hsa）
                      kegg.native = T, # 生成KEGG原生风格图
                      same.layer = T, 
                      new.signature = F,
                      discrete = list(gene=FALSE, cpd=TRUE),
                      limit = list(gene=1, cpd=1),
                      bins = list(gene=2, cpd=2),   # 颜色分阶数
                      low = "blue", mid = "orange", high = "red", # 颜色梯度
                      out.suffix = "metabo_plot" # 输出文件后缀
                    )
                    
                    pathview.data$plot.data.cpd$mol.col[match(all_features$id_kegg, pathview.data$plot.data.cpd$kegg.names)%>%na.omit%>%as.numeric()] <- "#D3D3D3"
                      pathview.data$plot.data.cpd <- update_node_colors(pathview.data$plot.data.cpd, pathway_metabo)
                      
                      keggview.native(plot.data.cpd = pathview.data$plot.data.cpd[,-10], 
                                      pathway.name = gsub("map", input$species, pathway_id), 
                                      cols.ts.cpd = pathview.data$plot.data.cpd[,10], 
                                      limit = list(gene=1, cpd=1), 
                                      bins = list(gene=2, cpd=2), 
                                      low = list(gene = "green", cpd = "blue"), 
                                      mid = list(gene = "gray", cpd = "orange"), 
                                      high = list(gene = "red", cpd = "red"),
                                      discrete = list(gene=T, cpd=T),
                                      plot.col.key = FALSE,
                                      out.suffix = "metabo_plot")
                      
                  },error=function(e){})
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".metabo_plot.png"))) {
                    
                    if (file.info(paste0(gsub("map", input$species, pathway_id), ".metabo_plot.png"))%>%.$size>0) {
                      
                      enrich_result$match[which(enrich_result$ID == pathway_id)] <- "TRUE"
                      
                      matched_n <- matched_n + 1
                      
                    }else{
                      
                      file.remove(paste0(gsub("map", input$species, pathway_id), ".metabo_plot.png"))
                      
                    }
                    
                  }
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".png"))) {
                    
                    file.remove(paste0(gsub("map", input$species, pathway_id), ".png"))
                    
                  }
                  
                  if (file.exists(paste0(gsub("map", input$species, pathway_id), ".xml"))) {
                    
                    file.remove(paste0(gsub("map", input$species, pathway_id), ".xml"))
                    
                  }
                  
                }
                
                setwd(report_path)
                
                enrich_result <- enrich_result[which(enrich_result$match=="TRUE"), ]
                
                ## 气泡图
                # 创建标签截断函数
                truncate_labels <- function(x) {
                  ifelse(nchar(x) > 30, paste0(substr(x, 1, 30), "..."), x)
                }
                
                if (nrow(enrich_result)>20) {
                  enrich_result2 <- enrich_result[1:20,]
                }else{
                  enrich_result2 <- enrich_result
                }
                
                p_dot <- ggplot(enrich_result2, aes(x = FoldEnrichment,y = forcats::fct_rev(Description))) +  # 排序
                  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.8) +
                  scale_color_gradient(
                    low = "#56B1F7",  # 深蓝
                    high = "#132B43",  # 亮蓝
                    space = "Lab"  # 使用Lab色彩空间增强渐变平滑度[6](@ref)
                  )+
                  scale_size(name = "Metabolite Count", limits = c(0, min(max(enrich_result2$Count),20))) +
                  scale_y_discrete(
                    labels = truncate_labels,  # 应用截断函数[7,8](@ref)
                    expand = expansion(add = 0.2)  # 增加标签间距避免截断重叠
                  )+
                  labs(x = "Fold Enrichment", y = "", title = "KEGG Pathway Enrichment (Sorted by Significance)") +  # 标题注明排序依据
                  theme_bw() +
                  theme(
                    axis.text.y = element_text(
                      size = 12, 
                      color = "black",
                      hjust = 1,        # 右对齐核心参数[7,8](@ref)
                      margin = margin(r = 10)  # 右侧留出10单位边距[1](@ref)
                    ), plot.margin = margin(l = 3.5, r = 1, unit = "cm")  # 增大左侧绘图边距
                  )
                
                ggsave(file.path(output_folder, "KEGG_dotplot.pdf"), p_dot, width = 10, height = 8)
                ggsave(file.path(output_folder, "KEGG_dotplot.png"), p_dot, width = 10, height = 8)
                
                # 表格导出
                enrich_long <- enrich_result %>%
                  dplyr::mutate(ID_split = str_split(metaboID, "/")) %>%  # 拆分斜杠分隔的ID
                  unnest(ID_split)
                
                enrich_with_names <- enrich_long %>%
                  dplyr::left_join(
                    kegg_data %>% dplyr::distinct(id_kegg, name),  # 去重以避免重复合并
                    by = c("ID_split" = "id_kegg")          # 将enrich的ID_split与kegg的id_kegg匹配
                  )
                
                # 聚合名称到原始行格式
                enrich_final <- enrich_with_names %>%
                  dplyr::group_by(ID) %>%                            # 按富集条目分组
                  dplyr::summarise(
                    name = paste(unique(na.omit(name)), collapse = "/")  # 去重后合并名称
                  ) %>%
                  dplyr::right_join(enrich_result, by = "ID") %>%    # 合并回原始数据框
                  dplyr::relocate(name, .after = "metaboID") 
                
                enrich_final[["ID"]] <- str_replace(enrich_final[["ID"]], "map", input$species)
                
                write.xlsx(enrich_final, file.path(output_folder, "KEGG_enrichment.xlsx"))
                
              },error=function(e){
                
                setwd(report_path)
                
              })
              
              setwd(report_path)
              
              ### 8.MSEA -----------------
              output_folder <- paste0("8.MSEA/", cmp_name)
              dir.create(output_folder, recursive = T, showWarnings = F)
              
              # 过滤无效ID（如空值或NA）
              msea_input <- separate_rows(results_mode_merged, "id_kegg")%>%as.data.frame()%>%.[!is.na(.$id_kegg) & .$id_kegg != "" & .$id_kegg != "NA", ]%>%dplyr::distinct(id_kegg,.keep_all = T)
              
              tryCatch({
                
                if (nrow(msea_input)>0) {
                  
                  if(exists("mSet")) rm(mSet)
                  
                  mSet <- InitDataObjects("conc", "msetora", FALSE, default.dpi = 72)
                  
                  system(paste0("cp -rf ", path_prefix, "/db/*.qs ", "./"))
                  
                  # Create vector consisting of compounds 
                  tmp.vec <- msea_input$id_kegg
                  
                  mSet<-Setup.MapData(mSet, tmp.vec)
                  
                  # Cross referencing
                  mSet<-CrossReferencing(mSet, "kegg")
                  
                  # Create the mapping results table
                  mSet<-CreateMappingResultTable(mSet)
                  
                  # Create list of candidates to replace the compound
                  mSet <- GetCandidateList(mSet)
                  
                  # Set the metabolite filter
                  mSet<-SetMetabolomeFilter(mSet, F);
                  
                  # Select metabolite set library, refer to 
                  mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 0);
                  
                  # Calculate hypergeometric score, results table generated in your working directory
                  mSet<-CalculateHyperScore(mSet)
                  
                  # Plot the ORA, bar-graph
                  mSet<-PlotORA(mSet, paste0(output_folder, "/MSEA_"), "bar", "png", 300, width=NA)
                  
                  # export report
                  msea_table <- mSet$analSet$ora.mat
                  
                  write.xlsx(msea_table%>%as.data.frame%>%rownames_to_column("Metabolite pathway"), paste0(output_folder, "/MSEA_result.xlsx"))
                  
                  # 绘制enrichment score
                  # 自定义函数生成ES曲线
                  plotMetabEnrichment <- function(mSet, pathway, output_folder){
                    # 获取排序后的代谢物列表（基于logFC）
                    ranked_metabs <- msea_input[order(-msea_input$FC), ]
                    
                    # 获取当前通路代谢物集合
                    pathway_metabs <- mSet[["dataSet"]][["map.table"]][,"Query"][match(mSet$analSet$ora.hits[[pathway]], mSet[["dataSet"]][["map.table"]][,"Match"])]
                    
                    # 计算累计ES值
                    es_profile <- nrow(ranked_metabs)
                    hit_indices <- which((ranked_metabs$id_kegg) %in% pathway_metabs)
                    nh <- length(hit_indices)
                    nm <- nrow(ranked_metabs)
                    
                    if(nh > 0){
                      
                      phit <- rep(0, nm)
                      phit[hit_indices] <- abs(ranked_metabs$FC[hit_indices])^1 # 权重指数可调
                      phit <- cumsum(phit/sum(phit))
                      
                      pmiss <- cumsum(!(1:nm %in% hit_indices))
                      pmiss <- pmiss/(nm - nh)
                      
                      es_profile <- phit - pmiss
                      
                    }
                    
                    # 绘制ES曲线
                    df <- data.frame(
                      Position = 1:nm,
                      ES = es_profile,
                      Type = ifelse(1:nm %in% hit_indices, "Hit", "Background")
                    )
                    
                    p <- ggplot(df, aes(x=Position)) +
                      geom_line(aes(y=ES), color="#2c7bb6", size=1) +
                      geom_rug(aes(color=Type), sides="b") +
                      geom_hline(yintercept=0, linetype="dashed") +
                      labs(title=paste("Enrichment Profile:", pathway),
                           y="Enrichment Score") +
                      theme_minimal()
                    
                    ggsave(paste0(output_folder, "/", pathway, "_ES_plot.png"), p, width=8, height=6)
                    
                  }
                  
                  # 生成所有显著通路的ES图
                  sig_pathways <- rownames(mSet$analSet$ora.mat)[mSet$analSet$ora.mat%>%as.data.frame()%>%.$`Raw p` < 0.05]
                  
                  lapply(sig_pathways, function(pw){
                    
                    plotMetabEnrichment(mSet, pw, output_folder)
                    
                  })
                  
                  # 生成综合GSEA风格热图（参考clusterProfiler的gseaplot2）
                  metab_ranks <- msea_input$FC
                  names(metab_ranks) <- msea_input$id_kegg
                  metab_ranks <- sort(metab_ranks, decreasing=TRUE)
                  
                  top_pathways <- head(sig_pathways, 10)
                  plotlist <- lapply(top_pathways, function(pw){
                    fgsea::plotEnrichment(
                      pathway = mSet[["dataSet"]][["map.table"]][,"Query"][match(mSet$analSet$ora.hits[[pw]], mSet[["dataSet"]][["map.table"]][,"Match"])],
                      stats = metab_ranks
                    ) + labs(title=pw)
                  })
                  
                  pdf(paste0(output_folder, "/Top10Pathways_ES.pdf"), width=11, height=8)
                  gridExtra::grid.arrange(grobs = plotlist[1:min(6, length(plotlist))], nrow=3, ncol=2)
                  dev.off()
                  
                }
                
              },error=function(e){
                
                
                
              })
              
              ### 9.数据库分类（差异代谢物-有 name） ----------------
              output_folder <- paste0("9.数据库分类/", cmp_name)
              dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
              
              ## 小工具：清洗名称（去两端空格、转小写、去掉末尾括号内容）
              .clean_nm <- function(x) {
                if (is.null(x)) return(NA_character_)
                x <- tolower(trimws(x))
                x <- sub("\\s*\\([^()]*\\)$", "", x)  # 去掉末尾括号
                trimws(x)
              }
              
              ## 1) 构建 HMDB / LMSD 的“名称 → 分类”映射表
              # HMDB: name + synonyms2 → class / sub_class
              hmdb_name_map <- hmdb %>%
                dplyr::select(name, synonyms2, class, sub_class) %>%
                tidyr::pivot_longer(cols = c(name, synonyms2), names_to = "src", values_to = "nm") %>%
                tidyr::separate_rows(nm, sep = ";\\s*") %>%
                dplyr::mutate(nm = .clean_nm(nm)) %>%
                dplyr::filter(!is.na(nm), nm != "") %>%
                dplyr::distinct(nm, class, sub_class)
              
              # LMSD: NAME + SYSTEMATIC_NAME → CATEGORY / MAIN_CLASS / SUB_CLASS
              lmsd_name_map <- lmsd %>%
                dplyr::select(NAME, SYSTEMATIC_NAME, CATEGORY, MAIN_CLASS, SUB_CLASS) %>%
                tidyr::pivot_longer(cols = c(NAME, SYSTEMATIC_NAME), names_to = "src", values_to = "nm") %>%
                dplyr::mutate(nm = .clean_nm(nm)) %>%
                dplyr::filter(!is.na(nm), nm != "") %>%
                dplyr::distinct(nm, CATEGORY, MAIN_CLASS, SUB_CLASS)
              
              ## 2) 取差异代谢物（必须有 name），拆分多名称并清洗
              diff_with_name <- results_mode_diff_name_merged %>%
                dplyr::filter(!is.na(name), trimws(name) != "") %>%
                tidyr::separate_rows(name, sep = ";") %>%
                dplyr::mutate(nm = .clean_nm(name)) %>%
                dplyr::filter(!is.na(nm), nm != "")
              
              ## 3) 关联到 HMDB / LMSD 分类
              diff_anno <- diff_with_name %>%
                dplyr::left_join(hmdb_name_map, by = "nm") %>%
                dplyr::left_join(lmsd_name_map, by = "nm") %>%
                dplyr::relocate(nm, .before = 1)
              
              # 导出明细
              readr::write_csv(
                diff_anno,
                file.path(output_folder, sprintf("差异代谢物_数据库分类明细_%s.csv", cmp_name))
              )
              
              ## 4) 通用画饼工具
              plot_cat_and_save <- function(df, col, title, out_prefix) {
                if (!col %in% names(df)) return(invisible(NULL))
                stat <- df %>%
                  dplyr::filter(!is.na(.data[[col]]), trimws(.data[[col]]) != "") %>%
                  dplyr::count(.data[[col]], name = "n") %>%
                  dplyr::arrange(dplyr::desc(n)) %>%
                  dplyr::mutate(pct = n / sum(n))
                
                if (nrow(stat) == 0) return(invisible(NULL))
                
                out_dir <- output_folder
                dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
                
                if (nrow(stat) <= 10) {
                  ## 先把非有限/非正的 pct 去掉，避免 geom_col() 报错
                  stat <- stat %>%
                    dplyr::mutate(pct = as.numeric(pct)) %>%
                    dplyr::filter(is.finite(pct), pct > 0)
                  
                  # 计算角度、侧别、坐标（在你的基础上微调为更安全的半径/范围）
                  stat <- stat %>%
                    dplyr::mutate(
                      ymax  = cumsum(pct),
                      ymin  = dplyr::lag(ymax, default = 0),
                      ymid  = (ymax + ymin) / 2,                 # 扇区中点(0-1)
                      angle = ymid * 360 - 90,
                      side  = dplyr::if_else(angle > 90 & angle < 270, "left", "right"),
                      label = paste0(.data[[col]], " (", scales::percent(pct), ")"),
                      # 径向段终点（略在饼图外）
                      x_rad = dplyr::if_else(side == "left", 0.93, 1.07),
                      # 肘形水平段终点（更外侧，便于排布标签）
                      x_elb = dplyr::if_else(side == "left", 0.75, 1.25),
                      # 标签放在肘形终点稍外一点
                      x_lab = dplyr::if_else(side == "left", 0.72, 1.28),
                      hjust = dplyr::if_else(side == "left", 1, 0)
                    )
                  
                  # 扇区（关键：limits 覆盖 [0.5,1.5] 并再留余量；width < 1 防边界擦边）
                  p <- ggplot(stat, aes(x = 1, y = pct, fill = .data[[col]])) +
                    geom_col(width = 0.92, color = "white", na.rm = TRUE) +  # na.rm 进一步静默无效行
                    coord_polar(theta = "y", clip = "off") +
                    theme_void() +
                    labs(title = title, fill = NULL) +
                    # 覆盖 [0.45,1.55]，确保整条环（中心1、厚度~0.92）在范围内，并给肘线/标签留空
                    scale_x_continuous(limits = c(0.45, 1.55), expand = expansion(mult = c(0, 0))) +
                    theme(
                      legend.position = "none",
                      plot.margin = margin(10, 110, 10, 90)   # 两侧更大的外边距，防截断
                    )
                  
                  # 径向引线（从 1.00 到 x_rad）
                  p <- p +
                    geom_segment(
                      data = dplyr::filter(stat, side == "right"),
                      aes(x = 1.00, xend = x_rad, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    ) +
                    geom_segment(
                      data = dplyr::filter(stat, side == "left"),
                      aes(x = 1.00, xend = x_rad, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    )
                  
                  # 肘形水平引线（x_rad -> x_elb）
                  p <- p +
                    geom_segment(
                      data = dplyr::filter(stat, side == "right"),
                      aes(x = x_rad, xend = x_elb, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    ) +
                    geom_segment(
                      data = dplyr::filter(stat, side == "left"),
                      aes(x = x_rad, xend = x_elb, y = ymid, yend = ymid),
                      inherit.aes = FALSE, linewidth = 0.4
                    )
                  
                  # 标签（优先 ggrepel，沿 y 方向自动避让；放在外圈 x_lab）
                  if (requireNamespace("ggrepel", quietly = TRUE)) {
                    p <- p +
                      ggrepel::geom_text_repel(
                        data = dplyr::filter(stat, side == "right"),
                        aes(x = x_lab, y = ymid, label = label),
                        inherit.aes = FALSE,
                        size = 3, direction = "y",
                        nudge_x = 0, box.padding = 0.3, point.padding = 0.25,
                        min.segment.length = 0, segment.size = 0.3, max.overlaps = Inf
                      ) +
                      ggrepel::geom_text_repel(
                        data = dplyr::filter(stat, side == "left"),
                        aes(x = x_lab, y = ymid, label = label),
                        inherit.aes = FALSE,
                        size = 3, direction = "y",
                        nudge_x = 0, box.padding = 0.3, point.padding = 0.25,
                        min.segment.length = 0, segment.size = 0.3, max.overlaps = Inf
                      )
                  } else {
                    p <- p +
                      geom_text(
                        data = dplyr::filter(stat, side == "right"),
                        aes(x = x_lab, y = ymid, label = label, hjust = 0),
                        inherit.aes = FALSE, size = 3
                      ) +
                      geom_text(
                        data = dplyr::filter(stat, side == "left"),
                        aes(x = x_lab, y = ymid, label = label, hjust = 1),
                        inherit.aes = FALSE, size = 3
                      )
                  }
                  
                } else {
                  stat <- stat %>% dplyr::mutate(label = scales::percent(pct))
                  max_n <- max(stat$n, na.rm = TRUE)
                  pad   <- max(1, max_n * 0.35)
                  
                  p <- ggplot(stat, aes(x = reorder(.data[[col]], n), y = n, fill = .data[[col]])) +
                    geom_col() +
                    coord_flip(clip = "off") +
                    scale_y_continuous(limits = c(0, max_n + pad),
                                       expand = expansion(mult = c(0, 0.05))) +
                    labs(title = paste0(title), x = NULL, y = "Count") +
                    theme_minimal() +
                    theme(
                      legend.position = "none",
                      plot.margin = margin(10, 90, 10, 10)
                    ) +
                    geom_text(aes(label = label), hjust = -0.05, size = 3)
                }
                
                ggsave(file.path(out_dir, paste0(out_prefix, ".png")), p, width = 8, height = 6, dpi = 300)
                ggsave(file.path(out_dir, paste0(out_prefix, ".pdf")), p, width = 8, height = 6)
                
                readr::write_csv(stat, file.path(out_dir, paste0(out_prefix, "_统计表.csv")))
              }
              
              ## 5) 分别绘制 5 类饼图
              # HMDB: class
              plot_cat_and_save(
                diff_anno, "class",
                title = sprintf("HMDB class - %s", cmp_name),
                out_prefix = sprintf("HMDB_class_%s", cmp_name)
              )
              
              # HMDB: sub_class
              plot_cat_and_save(
                diff_anno, "sub_class",
                title = sprintf("HMDB sub_class - %s", cmp_name),
                out_prefix = sprintf("HMDB_sub_class_%s", cmp_name)
              )
              
              # LMSD: CATEGORY
              plot_cat_and_save(
                diff_anno, "CATEGORY",
                title = sprintf("LMSD CATEGORY - %s", cmp_name),
                out_prefix = sprintf("LMSD_CATEGORY_%s", cmp_name)
              )
              
              # LMSD: MAIN_CLASS
              plot_cat_and_save(
                diff_anno, "MAIN_CLASS",
                title = sprintf("LMSD MAIN_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_MAIN_CLASS_%s", cmp_name)
              )
              
              # LMSD: SUB_CLASS
              plot_cat_and_save(
                diff_anno, "SUB_CLASS",
                title = sprintf("LMSD SUB_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_SUB_CLASS_%s", cmp_name)
              )
              
              ### 10.数据库分类（总代谢物-有 name） ----------------
              output_folder <- paste0("9.数据库分类/", cmp_name)
              dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
              
              ## 2) 取差异代谢物（必须有 name），拆分多名称并清洗
              total_with_name <- results_mode_merged %>%
                dplyr::filter(!is.na(name), trimws(name) != "") %>%
                tidyr::separate_rows(name, sep = ";") %>%
                dplyr::mutate(nm = .clean_nm(name)) %>%
                dplyr::filter(!is.na(nm), nm != "")
              
              ## 3) 关联到 HMDB / LMSD 分类
              total_anno <- total_with_name %>%
                dplyr::left_join(hmdb_name_map, by = "nm") %>%
                dplyr::left_join(lmsd_name_map, by = "nm") %>%
                dplyr::relocate(nm, .before = 1)
              
              # 导出明细
              readr::write_csv(
                total_anno,
                file.path(output_folder, sprintf("总代谢物_数据库分类明细_%s.csv", cmp_name))
              )
              
              ## 4) 分别绘制 5 类饼图
              # HMDB: class
              plot_cat_and_save(
                total_anno, "class",
                title = sprintf("HMDB class - %s", cmp_name),
                out_prefix = sprintf("HMDB_class_all_%s", cmp_name)
              )
              
              # HMDB: sub_class
              plot_cat_and_save(
                total_anno, "sub_class",
                title = sprintf("HMDB sub_class - %s", cmp_name),
                out_prefix = sprintf("HMDB_sub_class_all_%s", cmp_name)
              )
              
              # LMSD: CATEGORY
              plot_cat_and_save(
                total_anno, "CATEGORY",
                title = sprintf("LMSD CATEGORY - %s", cmp_name),
                out_prefix = sprintf("LMSD_CATEGORY_all_%s", cmp_name)
              )
              
              # LMSD: MAIN_CLASS
              plot_cat_and_save(
                total_anno, "MAIN_CLASS",
                title = sprintf("LMSD MAIN_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_MAIN_CLASS_all_%s", cmp_name)
              )
              
              # LMSD: SUB_CLASS
              plot_cat_and_save(
                total_anno, "SUB_CLASS",
                title = sprintf("LMSD SUB_CLASS - %s", cmp_name),
                out_prefix = sprintf("LMSD_SUB_CLASS_all_%s", cmp_name)
              )
              
            }
            
            cat("已保存", cmp_name, "的结果至文件。\n")
            
            # 统计各阶段的feature数量及差异筛选的上下调数量 ----------
            diff_count <- nrow(results_mode_diff_merged)
            diff_up <- sum(results_mode_diff_merged$difference == "up", na.rm = TRUE)
            diff_down <- sum(results_mode_diff_merged$difference == "down", na.rm = TRUE)
            diff_name_count <- nrow(results_mode_diff_name_merged)
            diff_name_up <- sum(results_mode_diff_name_merged$difference == "up", na.rm = TRUE)
            diff_name_down <- sum(results_mode_diff_name_merged$difference == "down", na.rm = TRUE)
            
            if (!input$single_ion_judge){
              
              results_mode_diff_merged_POS <- results_mode_diff_merged[grepl("POS", results_mode_diff_merged$peak_name), ]
              results_mode_diff_name_merged_POS <- results_mode_diff_name_merged[grepl("POS", results_mode_diff_name_merged$peak_name), ]
              results_mode_diff_merged_NEG <- results_mode_diff_merged[grepl("NEG", results_mode_diff_merged$peak_name), ]
              results_mode_diff_name_merged_NEG <- results_mode_diff_name_merged[grepl("NEG", results_mode_diff_name_merged$peak_name), ]
              
              # 将当前比较的统计信息保存到列表中
              summary_stats_list[[cmp_name]][["POS"]] <- data.frame(
                compare = cmp_name,
                ion_type = "POS",
                raw_feature_num = summary_stats_list[[cmp_name]][["POS"]]$raw_feature_num,
                remove_missing_num = summary_stats_list[[cmp_name]][["POS"]]$remove_missing_num,
                preprocess_feature_num = summary_stats_list[[cmp_name]][["POS"]]$preprocess_feature_num,
                preprocess_feature_num_name = summary_stats_list[[cmp_name]][["POS"]]$preprocess_feature_num_name,
                diff_feature = nrow(results_mode_diff_merged_POS),
                diff_feature_up = sum(results_mode_diff_merged_POS$difference == "up", na.rm = TRUE),
                diff_feature_down = sum(results_mode_diff_merged_POS$difference == "down", na.rm = TRUE),
                diff_feature_name = nrow(results_mode_diff_name_merged_POS),
                diff_feature_name_up = sum(results_mode_diff_name_merged_POS$difference == "up", na.rm = TRUE),
                diff_feature_name_down = sum(results_mode_diff_name_merged_POS$difference == "down", na.rm = TRUE),
                stringsAsFactors = FALSE
              )
              
              summary_stats_list[[cmp_name]][["NEG"]] <- data.frame(
                compare = cmp_name,
                ion_type = "NEG",
                raw_feature_num = summary_stats_list[[cmp_name]][["NEG"]]$raw_feature_num,
                remove_missing_num = summary_stats_list[[cmp_name]][["NEG"]]$remove_missing_num,
                preprocess_feature_num = summary_stats_list[[cmp_name]][["NEG"]]$preprocess_feature_num,
                preprocess_feature_num_name = summary_stats_list[[cmp_name]][["NEG"]]$preprocess_feature_num_name,
                diff_feature = nrow(results_mode_diff_merged_NEG),
                diff_feature_up = sum(results_mode_diff_merged_NEG$difference == "up", na.rm = TRUE),
                diff_feature_down = sum(results_mode_diff_merged_NEG$difference == "down", na.rm = TRUE),
                diff_feature_name = nrow(results_mode_diff_name_merged_NEG),
                diff_feature_name_up = sum(results_mode_diff_name_merged_NEG$difference == "up", na.rm = TRUE),
                diff_feature_name_down = sum(results_mode_diff_name_merged_NEG$difference == "down", na.rm = TRUE),
                stringsAsFactors = FALSE
              )
              
              summary_stats_list[[cmp_name]][["Merge"]] <- data.frame(
                compare = cmp_name,
                ion_type = "POS&NEG",
                raw_feature_num = summary_stats_list[[cmp_name]][["POS"]]$raw_feature_num+summary_stats_list[[cmp_name]][["NEG"]]$raw_feature_num,
                remove_missing_num = summary_stats_list[[cmp_name]][["POS"]]$remove_missing_num+summary_stats_list[[cmp_name]][["NEG"]]$remove_missing_num,
                preprocess_feature_num = nrow(results_mode_merged),
                preprocess_feature_num_name = sum(!is.na(results_mode_merged$name)),
                diff_feature = nrow(results_mode_diff_merged),
                diff_feature_up = sum(results_mode_diff_merged$difference == "up", na.rm = TRUE),
                diff_feature_down = sum(results_mode_diff_merged$difference == "down", na.rm = TRUE),
                diff_feature_name = nrow(results_mode_diff_name_merged),
                diff_feature_name_up = sum(results_mode_diff_name_merged$difference == "up", na.rm = TRUE),
                diff_feature_name_down = sum(results_mode_diff_name_merged$difference == "down", na.rm = TRUE),
                stringsAsFactors = FALSE
              )
              
            }else{
              
              summary_stats_list[[cmp_name]][["single"]] <- data.frame(
                compare = cmp_name,
                ion_type = "single",
                raw_feature_num = summary_stats_list[[cmp_name]][["single"]]$raw_feature_num,
                remove_missing_num = summary_stats_list[[cmp_name]][["single"]]$remove_missing_num,
                preprocess_feature_num = nrow(results_mode_merged),
                preprocess_feature_num_name = sum(!is.na(results_mode_merged$name)),
                diff_feature = nrow(results_mode_diff_merged),
                diff_feature_up = sum(results_mode_diff_merged$difference == "up", na.rm = TRUE),
                diff_feature_down = sum(results_mode_diff_merged$difference == "down", na.rm = TRUE),
                diff_feature_name = nrow(results_mode_diff_name_merged),
                diff_feature_name_up = sum(results_mode_diff_name_merged$difference == "up", na.rm = TRUE),
                diff_feature_name_down = sum(results_mode_diff_name_merged$difference == "down", na.rm = TRUE),
                stringsAsFactors = FALSE
              )
              
            }
            
          }
          
          # 合并所有分组比较的统计信息
          final_summary <- do.call(rbind, do.call(rbind, summary_stats_list))
          
          # 保存结果rda数据
          save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw", "final_summary"), file = "result.rda")
          
          # 保存最终统计结果到 "差异筛选/差异数量统计.xlsx"
          write.xlsx(final_summary, file = "5.差异筛选/差异数量统计.xlsx", rowNames = FALSE)
          cat("已保存差异数量统计结果至文件：5.差异筛选/差异数量统计.xlsx\n")
          
          # 组间差异结果 ---------------
          tryCatch({
            
            ## feature --------------------
            diff_features_list <- list()
            
            # 遍历所有比较组
            for(cmp_name in group.data) {
              
              diff_features_list[[cmp_name]] <- do.call(rbind, results_mode_diff[[cmp_name]])%>%.[,c("peak_name")]
              
            }
            
            ### Venn for names ---------------------
            if(length(group.data) >= 2 & length(group.data) <=5) {
              
              output_folder <- paste0("6.差异代谢物可视化/1.venn")
              dir.create(output_folder, recursive = T)
              
              # 动态选择调色板（支持2-5组）
              if(length(diff_features_list) == 2) {
                
                fill_colors <- brewer.pal(3, "Paired")[1:2]  # 使用Paired调色板前两色
                
              } else {
                
                fill_colors <- brewer.pal(length(diff_features_list), "Set2")  # Set2支持3-8色
                
              }
              
              venn.plot <- venn.diagram(
                x = diff_features_list,
                filename = file.path(output_folder, "Venn_diagram.png"),
                disable.logging = T,
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
              
              # 生成Venn交集信息（修正列名问题）
              venn_data <- get.venn.partitions(diff_features_list) %>%
                rename(
                  Set_Combination = ..set..,  # 根据网页1的列名修正
                  Peak_Names = ..values..,
                  Count = ..count..
                )
              
              # 合并所有差异代谢物完整信息
              full_data <- do.call(rbind, lapply(names(results_mode_diff), function(cmp_name){
                cbind(Comparison = cmp_name, 
                      do.call(rbind, results_mode_diff[[cmp_name]]) %>% 
                        # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb))
                        dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg))
              })) %>% dplyr::distinct(peak_name, .keep_all = TRUE)
              
              # 初始化Excel工作簿
              wb <- createWorkbook()
              
              # 初始化一个空数据框存储所有数据
              combined_data <- data.frame()
              
              # 遍历每个Venn交集区域
              for(i in 1:nrow(venn_data)) {
                # 获取当前交集区域的peak_name列表
                current_peaks <- unlist(venn_data$Peak_Names[i])
                
                # 提取代谢物信息并添加标识列
                sheet_data <- full_data %>%
                  dplyr::filter(peak_name %in% current_peaks) %>%
                  # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb) %>%
                  dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg) %>%
                  dplyr::mutate(Set_Combination = venn_data$Set_Combination[i])  # 关键：添加标识列[2,6](@ref)
                
                # 合并到总数据框
                combined_data <- dplyr::bind_rows(combined_data, sheet_data)  # 使用dplyr的合并方法[6](@ref)
              }
              
              # 创建单个工作表并写入数据
              addWorksheet(wb, sheetName = "Combined_Venn_Data")
              writeData(wb, sheet = "Combined_Venn_Data", x = combined_data)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder, "Venn_diagram.xlsx"), overwrite = TRUE)
              
            }
            
            ### Upset ------------------------
            if(length(group.data) >= 2) {
              
              output_folder <- paste0("6.差异代谢物可视化/2.upset")
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
            # 收集所有比较组的差异代谢物名称（基于name）
            diff_names_list <- list()
            
            # 遍历所有比较组
            for(cmp_name in group.data) {
              
              diff_names_list[[cmp_name]] <- do.call(rbind, results_mode_diff_name[[cmp_name]])%>%.[,c("peak_name")]
              
            }
            
            ### Venn for names ---------------------
            if(length(diff_names_list) >= 2 & length(diff_names_list) <= 5) {
              
              output_folder_name <- paste0("6.差异代谢物可视化/1.venn")
              
              # 动态选择调色板
              if(length(diff_names_list) == 2) {
                
                fill_colors <- brewer.pal(3, "Paired")[1:2]
                
              } else {
                
                fill_colors <- brewer.pal(length(diff_names_list), "Set2")
                
              }
              
              venn.plot <- venn.diagram(
                x = diff_names_list,
                filename = file.path(output_folder_name, "Venn_diagram_name.png"),
                disable.logging = T,
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
              
              # 生成Venn交集信息（修正列名问题）
              venn_data <- get.venn.partitions(diff_names_list) %>%
                rename(
                  Set_Combination = ..set..,  # 根据网页1的列名修正
                  Peak_Names = ..values..,
                  Count = ..count..
                )
              
              # 合并所有差异代谢物完整信息
              full_data <- do.call(rbind, lapply(names(results_mode_diff), function(cmp_name){
                cbind(Comparison = cmp_name, 
                      do.call(rbind, results_mode_diff[[cmp_name]]) %>% 
                        # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb))
                        dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg))
              })) %>% dplyr::distinct(peak_name, .keep_all = TRUE)
              
              # 初始化Excel工作簿
              wb <- createWorkbook()
              
              # 初始化一个空数据框存储所有数据
              combined_data <- data.frame()
              
              # 遍历每个Venn交集区域
              for(i in 1:nrow(venn_data)) {
                # 获取当前交集区域的peak_name列表
                current_peaks <- unlist(venn_data$Peak_Names[i])
                
                # 提取代谢物信息并添加标识列
                sheet_data <- full_data %>%
                  dplyr::filter(peak_name %in% current_peaks) %>%
                  # dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg, id_hmdb) %>%
                  dplyr::select(peak_name, name, pvalue, adj_p_value, Log2FC, VIP, difference, id_kegg) %>%
                  dplyr::mutate(Set_Combination = venn_data$Set_Combination[i])  # 关键：添加标识列[2,6](@ref)
                
                # 合并到总数据框
                combined_data <- dplyr::bind_rows(combined_data, sheet_data)  # 使用dplyr的合并方法[6](@ref)
              }
              
              # 创建单个工作表并写入数据
              addWorksheet(wb, sheetName = "Combined_Venn_Data")
              writeData(wb, sheet = "Combined_Venn_Data", x = combined_data)
              
              # 保存Excel文件
              saveWorkbook(wb, file = file.path(output_folder_name, "Venn_diagram_name.xlsx"), overwrite = TRUE)
              
            }
            
            ### Upset for names ------------------------
            if(length(diff_names_list) >= 2) {
              
              output_folder_name <- paste0("6.差异代谢物可视化/2.upset")
              
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
              
              pdf(file.path(output_folder_name, "UpSet_plot_name.pdf"), width = 20, height = 12)
              print(upset_plot_name)
              dev.off()
              
              png(file.path(output_folder_name, "UpSet_plot_name.png"), width = 800, height = 600)
              print(upset_plot_name)
              dev.off()
              
              write.xlsx(binary_matrix_name, file.path(output_folder_name, "UpSet_data_name.xlsx"))
              
            }
            
          },error=function(e){
            
            
            
          })
          
        }else{
          
          if (!input$single_ion_judge){
            
            final_summary <- data.frame(
              compare = "no compare",
              ion_type = "POS&NEG",
              raw_feature_num = summary_stats_list[["POS"]]$raw_feature_num+summary_stats_list[["NEG"]]$raw_feature_num,
              remove_missing_num = summary_stats_list[["POS"]]$remove_missing_num+summary_stats_list[["NEG"]]$remove_missing_num,
              preprocess_feature_num = summary_stats_list[["POS"]]$preprocess_feature_num+summary_stats_list[["NEG"]]$preprocess_feature_num,
              preprocess_feature_num_name = summary_stats_list[["POS"]]$preprocess_feature_num_name+summary_stats_list[["NEG"]]$preprocess_feature_num_name,
              stringsAsFactors = FALSE
            )
            
          }
          
          # 合并所有分组比较的统计信息
          final_summary <- do.call(rbind, do.call(rbind, summary_stats_list))
          
          # 保存结果rda数据
          save(list = c("summary_stats_list", "results_mode_diff_name", "results_mode_diff", "results_mode", "metabo.data.raw", "final_summary"), file = "result.rda")
          
        }
        
        # 出报告 -----------------
        para_sum <- list(miss=c("Within the group" = "inter.group",
                                "Global filter" = "global.group"),
                         fill=c("Intra-group mean" = "mean.group", 
                                "Global mean" = "mean.global", 
                                "Median within the group" = "median.group",
                                "Global median" = "median.global", 
                                "None" = "none",
                                "Min" = "min.global",
                                "Min/2" = "min2",
                                "KNN" = "knn.global",
                                "Randomforest" = "rf"),
                         normalize=c("Sum" = "sum",
                                     "None" = "none",
                                     "QC svr"= "qc",
                                     "QC RLSC"= "qc-rlsc",
                                     "QC NormAE"= "normae",
                                     "Probability quotient" = "prob_quot",
                                     "0.75 quantile" = "percent_0.75")
        )
        
        para_data <- data.frame("缺失值类型" = input$miss.value.handle.type, 
                                "缺失值过滤方法" = names(para_sum$miss)[which(para_sum$miss==input$miss.value.handle.group)], 
                                "缺失值过滤比率" = input$miss.value.handle.cutoff,
                                "填充方式" = names(para_sum$fill)[which(para_sum$fill==input$miss.value.fill)],
                                "归一化方法" = names(para_sum$normalize)[which(para_sum$normalize==input$normalized.handle.method)],
                                "SUM的系数" = input$sum_coef,
                                "定量值取log" = input$log.handle.method,
                                "RSD过滤阈值" = input$rsd.cutoff,
                                "RSD取log" = input$log.rsd.method,
                                "批次矫正" = input$batch.correct,
                                "相关性阈值" = input$correlation.cutoff,
                                "是否使用FC筛选" = input$fc_cutoff_judge,
                                "Fold Change Cutoff" = input$fc_cutoff,
                                "是否使用pvalue筛选" = input$p_cutoff_judge,
                                "检验类型" = input$pvalue_type,
                                "P-value Cutoff" = input$p_cutoff,
                                "是否使用p-adjust筛选" = input$padjust_judge,
                                "padjust方法" = input$padjust_method,
                                "是否使用VIP筛选" = input$vip_cutoff_judge,
                                "VIP Cutoff" = input$vip_cutoff,
                                "聚类热图是否使用样本聚类" = input$hca_col_cluster,
                                "KEGG物种" = input$species,
                                "选择的列数" = paste0(selected_column, collapse = ";"))%>%t%>%data.frame()
        
        if (input$single_ion_judge) {
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution_single.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution_server_single.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report_server.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution_single_word.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.docx'))
        }else{
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution_server.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report_server.html'))
          rmarkdown::render(paste0(path_prefix, "/etc/dep/report_templete/quant/metabolite_report_omicsolution_word.Rmd"), output_file = paste0(path_prefix, '/data/usrdata/', project_id, '/report/report.docx'))
        }
        
      }
      
      # 打包和通知 ----
      files_to_extract_path <- c(
        list.dirs(report_path, full.names = TRUE, recursive = FALSE),
        file.path(report_path, "report.html")
      )
      
      if (sum(grepl("svr_normalization_result", files_to_extract_path))>0) {
        files_to_extract_path <- files_to_extract_path[-grep("svr_normalization_result", files_to_extract_path)]
      }
      
      # # 创建临时目录
      # temp_dir <- tempfile("report_temp_")
      # dir_create(temp_dir)
      # 
      # # 把每个待压缩的项目复制到临时目录下
      # for (f in files_to_extract_path) {
      #   system(paste("cp -rf", f, temp_dir))
      # }
      # 
      # files_to_extract_path2 <- c(
      #   list.dirs(temp_dir, full.names = TRUE, recursive = FALSE),
      #   file.path(temp_dir, "report.html")
      # )
      
      # # 压缩临时目录下的所有内容
      # zipfile_path <- file.path(temp_dir, "report.zip")
      zipfile_path <- file.path(report_path, "report.zip")
      
      zip::zip(
        zipfile = zipfile_path,
        files = files_to_extract_path,
        include_directories = TRUE,
        mode = "cherry-pick"
      )
      
      # # 复制压缩包回原目录
      # system(paste("cp -f", zipfile_path, file.path(report_path, "report.zip")))
      
      UpdateStatus(2,project_id)
      
      #发送邮件
      mail <- input$email
      
      message <- paste0('From: "Omicsolution" <marketing@omicsolution.com>
To: "',stringr::str_replace(mail,"\\@.*",""),'" <',mail,'>
Subject: Task finished inform!

Dear user,

Your project (', input$project_name, ') have finished.Please go to https://project.omicsolution.com/BioAnalysis/?demo&id2=',project_id,' to check the result.')
      
      curl::send_mail(mail_from = "marketing@omicsolution.com",
                      mail_rcpt = mail, 
                      message, 
                      smtp_server = 'smtps://smtp.exmail.qq.com',
                      username = 'marketing@omicsolution.com', 
                      password  = 'ERPUuVncXNJjc8gT')
      
    },error=function(e){
      
      UpdateStatus(3,project_id)
      
      #发送邮件
      mail <- input$email
      
      message <- paste0('From: "Omicsolution" <marketing@omicsolution.com>
To: "',stringr::str_replace(mail,"\\@.*",""),'" <',mail,'>
Subject: Task finished inform!

Dear user,

Your project ', input$project_name, '(', project_id ,') have failure.',e,' Please check the raw data and group file.')
      
      curl::send_mail(mail_from = "marketing@omicsolution.com",
                      mail_rcpt = mail, 
                      message, 
                      smtp_server = 'smtps://smtp.exmail.qq.com',
                      username = 'marketing@omicsolution.com', 
                      password  = 'ERPUuVncXNJjc8gT')
      
    })
    
  }
  
}