# 主界面 ----
output$ident_panel <- renderUI({
  
  box(solidHeader = TRUE, status = "primary", width = 12, style = 'margin: 0 auto; min-height: 500px;',
      
      uiOutput("ident_species_selector"),
      
      # uiOutput("ident_project_confirm"),
      
      shinyjs::useShinyjs(),
      
      useShinyalert()
      
  )
  
})

observeEvent(input$ident_selected_species, {
  
  species_sci_name <<- input$ident_selected_species %>% str_replace_all('( \\().*', '')
  
  species_index <<- which(species_list$sci_names == species_sci_name)
  
  species_db_name <<- species_list$name[species_index]
  
  species_ch_names <<- species_list$ch_names[species_index]
  
  species_kegg <<- species_list$kegg[species_index]
  
  species_taxid <<- species_list$taxid[species_index]
  
  species_updatetime <<- species_list$updatetime[species_index]
  
})

output$ident_species_selector <- renderUI({
  
  speciesdb <- dbConnect(RSQLite::SQLite(), "/home/biognosis/mnt/os2/BioAnalysis/db/speciesdb.db")	
  
  species_list <<- dbGetQuery(speciesdb, paste0("SELECT * FROM Species")) %>% as.data.frame
  
  dbDisconnect(speciesdb)
  
  choices <<- list()
  
  species_L1 <- species_list$species_L1 %>% unique
  
  for(i in 1:length(species_L1)){
    
    species_list_sub <- species_list[which(species_list$species_L1 == species_L1[i]),]
    
    species_names <<- species_list_sub$sci_names
    
    species_ch_names <<- species_list_sub$ch_names
    
    choices_string <- paste0(species_names, ' (', species_ch_names, ')')
    
    choices[[species_L1[i]]] <<- choices_string %>% as.list
    
  }
  
  fluidRow(class = 'para_panel', style = 'margin-right: auto; margin-left: auto;',
           
           column(2, '' ),
           
           column(8,
                  
                  box(title = "Upload File", solidHeader = TRUE, status = "primary", width = 12, style = '',
                      
                      div(style = ' max-width: 78%;  margin-top: 1em;  margin-right: auto; margin-left: auto; ',
                          
                          fileInput("glycoRaw", paste0("Upload data file:"),
                                    multiple = FALSE,
                                    accept = c(".zip"))),
                      
                      div(style = ' max-width: 78%;  margin-top: 1em;  margin-right: auto; margin-left: auto; ',
                          
                          selectizeInput(inputId = 'ident_selected_species', label = NULL, choices = choices, selected = NULL, multiple = FALSE, options = list(openOnFocus = FALSE))),
                      
                      div(style = 'max-width: 70%; margin-bottom: 3em; margin-right: auto; margin-left: auto;',
                          
                          shiny::downloadLink(outputId = 'ident_downloadFasta', label = shinysky::actionButton(inputId = 'dl_btn_fasta', label = 'FASTA Download', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 1em; margin-top: 2em;') %>% paste0 %>% HTML())
                          
                      )),
                  
                  fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
                           
                           column(offset = 0,width = 12,
                                  
                                  div(style = 'max-width: 600px; margin-right: auto; margin-left: auto;',
                                      
                                      box(title = "Parameter Info", solidHeader = TRUE, status = "primary", width = 12, style = '',
                                          
                                          div(style = 'max-width: 85%; margin: 3em auto;',
                                              
                                              selectInput(inputId = "miss.value.handle.method",
                                                          "PCA缺失值处理方式:",
                                                          c("不含缺失值" = "none",
                                                            "缺失值全部设为0" = "zero")
                                                          ),
                                              
                                              selectInput(inputId = "miss.value.handle.show",
                                                          "PCA展示选项：",
                                                          c("展示组别 + QC + 批次" = "normal.qc",
                                                            "展示组别 + 批次" = "normal")
                                              ),
                                              
                                              selectInput(inputId = "miss.value.handle.batch",
                                                          "缺失值过滤范围:",
                                                          c("批次内" = "none.batch",
                                                            "跨批次" = "across.batch")
                                              ),
                                              
                                              selectInput(inputId = "miss.value.handle.group",
                                                          "过滤方式:",
                                                          c("组内缺失值占比过滤" = "inter.group",
                                                            "全局缺失值占比过滤" = "global.group")
                                              ),
                                              
                                              numericInput(inputId = "miss.value.handle.cutoff",
                                                          "缺失值阈值过滤:",
                                                          value=50,min=0,max=90,step = 10
                                              ),
                                              
                                              uiOutput("miss.value.handle.o"),
                                              
                                              selectInput(inputId = "miss.value.handle.remain",
                                                          "剩余缺失值处理方法:",
                                                          c("缺失值做填充" = "remain.fill",
                                                            "保留包含缺失值的行，但在数据集中标记这些行" = "remain.reserve")
                                              ),
                                              
                                              uiOutput("miss.value.fill.o"),
                                              
                                              uiOutput("miss.value.valid.o"),
                                              
                                              selectInput(inputId = "normalized.handle.method",
                                                          "中位数处理:",
                                                          c("不处理" = "none",
                                                            "中位数平滑" = "median.polish",
                                                            "中位数归一化"= "median")
                                              ),
                                              
                                              uiOutput("batch.correct.o"),
                                              
                                              selectInput(inputId = "pca.handle.method",
                                                          "处理后的数据PCA:",
                                                          c("不含缺失值" = "none",
                                                            "缺失值全部设为0" = "zero")
                                              ),
                                              
                                              selectInput(inputId = "pca.handle.show",
                                                          "处理后的PCA展示选项：",
                                                          c("展示组别 + QC + 批次" = "normal.qc",
                                                            "展示组别 + 批次" = "normal")
                                              ),
                                              
                                              p('Project Name:', style = 'font-size: 14px; font-weight: bold;'),
                                              
                                              div(style = '',textInput('project_name', label = NULL, value = paste0('Project_',format(Sys.time(),"%Y%m%d%H%M")))),
                                              
                                              p('Organization Information (And fill in contact information):', style = 'font-size: 14px; font-weight: bold;'),
                                              
                                              div(style = '',textInput('organ_info', label = NULL)),
                                              
                                              p('E-mail Address:', style = 'font-size: 14px; font-weight: bold;'),
                                              
                                              div(style = '',textInput('email', label = NULL, value = ''))),
                                          
                                          div(style = 'width: 100%; text-align:center;',shinysky::actionButton(inputId = 'confirm_ident1', label = 'Confirm', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-bottom: 0; margin-top: 3em;'))
                                          
                                      )
                                      
                                  )
                                  
                           ))
                  
           ),
           
           column(2, '' )
           
  )
  
})

observe({
  
  if (!is.null(input$miss.value.handle.batch) & !is.null(input$miss.value.handle.group)) {
    
    if (input$miss.value.handle.batch=="none.batch") {
      
      if (input$miss.value.handle.group=="inter.group") {
        
        output$miss.value.handle.o <- renderUI({
          
          selectInput(inputId = "miss.value.handle.condition",
                      "过滤条件：",
                      c("所有组别都满足过滤条件即删除该行" = "none.inter.condition.all",
                        "任何一个组别满足过滤条件即删除该行" = "none.inter.condition.any"))
          
        })
        
      }else if (input$miss.value.handle.group=="global.group") {
        
        output$miss.value.handle.o <- renderUI({
          
          selectInput(inputId = "miss.value.handle.condition",
                      "过滤条件：",
                      c("满足过滤条件即删除该行" = "none.global.condition"))
          
        })
        
      }
      
    }else if (input$miss.value.handle.batch=="across.batch") {
        
        output$miss.value.handle.o <- renderUI({
          
          selectInput(inputId = "miss.value.handle.condition",
                      "过滤条件：",
                      c("所有批次都满足过滤条件即删除该行" = "across.condition.all",
                        "任何一个批次满足过滤条件即删除该行" = "across.condition.any"))
          
        })
        
      }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.handle.remain)) {
    
    if (input$miss.value.handle.remain=="remain.fill") {
      
      output$miss.value.fill.o <- renderUI({
        
        fluidPage(
          
          selectInput(inputId = "miss.value.fill",
                      "填充方式：",
                      c("不填充" = "none",
                        "根据组内有效数据占比"="inter.valid",
                        "根据全局有效数据占比"="global.valid",
                        "简洁填充:组内KNN" = "knn.inter",
                        "简洁填充:组内最小值" = "min.inter",
                        "简洁填充:组内均值" = "mean.inter",
                        "简洁填充:组内中位数" = "median.inter",
                        "简洁填充:全局最小值" = "min.global",
                        "简洁填充:全局均值" = "mean.global",
                        "简洁填充:全局中位数" = "median.global",
                        "简洁填充:全局KNN" = "knn.global",
                        "简洁填充:用户自定义的值" = "user.defined")
          ),
          
          uiOutput("miss.value.fill.defined.o")
          
        )
        
      })
      
    }else{
      
      output$miss.value.fill.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.fill)) {
    
    if (input$miss.value.fill=="user.defined") {
      
      output$miss.value.fill.defined.o <- renderUI({
        
        numericInput(inputId = "miss.value.fill.defined","用户自定义的数值：",value = 0)
        
      })
      
    }else{
      
      output$miss.value.fill.defined.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.fill)) {
    
    if (input$miss.value.fill == "inter.valid") {
      
      output$miss.value.valid.o <- renderUI({
        
        fluidPage(
          
          selectInput(inputId = "miss.value.valid.below50",
                      "如果组内有效数据占比< 50%：",
                      c("组内最小值" = "min.inter",
                        "组内均值" = "mean.inter",
                        "组内中位数" = "median.inter",
                        "组内KNN" = "knn.inter",
                        "用户自定义的数值" = "user.defined")),
          
          uiOutput("miss.value.valid.below50.defined.o"),
          
          selectInput(inputId = "miss.value.valid.over50",
                      "如果组内有效数据占比 ≥ 50%：",
                      c("组内最小值" = "min.inter",
                        "组内均值" = "mean.inter",
                        "组内中位数" = "median.inter",
                        "组内KNN" = "knn.inter",
                        "用户自定义的数值" = "user.defined")),
          
          uiOutput("miss.value.valid.over50.defined.o")
          
        )
        
      })
      
    }else if (input$miss.value.fill == "global.valid") {
      
      output$miss.value.valid.o <- renderUI({
        
        fluidPage(
          
          selectInput(inputId = "miss.value.valid.below50",
                      "如果全局有效数据占比< 50%：",
                      c("全局最小值" = "min.inter",
                        "全局均值" = "mean.inter",
                        "全局中位数" = "median.inter",
                        "全局KNN" = "knn.inter",
                        "用户自定义的数值" = "user.defined")),
          
          uiOutput("miss.value.valid.below50.defined.o"),
          
          selectInput(inputId = "miss.value.valid.over50",
                      "如果全局有效数据占比 ≥ 50%：",
                      c("全局最小值" = "min.inter",
                        "全局均值" = "mean.inter",
                        "全局中位数" = "median.inter",
                        "全局KNN" = "knn.inter",
                        "用户自定义的数值" = "user.defined")),
          
          uiOutput("miss.value.valid.over50.defined.o")
          
        )
        
      })
      
    }else{
      
      output$miss.value.valid.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.valid.below50)) {
    
    if (input$miss.value.valid.below50=="user.defined") {
      
      output$miss.value.valid.below50.defined.o <- renderUI({
        
        numericInput(inputId = "miss.value.valid.below50.defined","用户自定义的数值：",value = 0)
        
      })
      
    }else{
      
      output$miss.value.valid.below50.defined.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.valid.over50)) {
    
    if (input$miss.value.valid.over50=="user.defined") {
      
      output$miss.value.valid.over50.defined.o <- renderUI({
        
        numericInput(inputId = "miss.value.valid.over50.defined","用户自定义的数值：",value = 0)
        
      })
      
    }else{
      
      output$miss.value.valid.over50.defined.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$miss.value.handle.batch)) {
    
    if (input$miss.value.handle.batch=="across.batch") {
      
      output$batch.correct.o <- renderUI({
        
        fluidPage(
          
          selectInput(inputId = "batch.correct",
                      "批次矫正：",
                      c("不矫正" = "none",
                        "COMBAT" = "combat",
                        "参考样本（指定一个已经上传的样本，做除法）" = "reference.sample")
          ),
          
          uiOutput("batch.correct.reference.sample.o")
          
        )
        
      })
      
    }else{
      
      output$batch.correct.o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$batch.correct)) {
    
    if (input$batch.correct=="reference.sample") {
      
      output$batch.correct.reference.sample.o <- renderUI({
        
        textInput(inputId = "batch.correct.reference.sample",
                  label = "指定一个已经上传的样本名称（请使用IdentificationsOverview.xls表格中的Sample name(分析用的名称）列里的名称）：")
        
      })
      
    }else{
      
      output$batch.correct.reference.sample.o <- NULL
      
    }
    
  }
  
})

# download ----
output$ident_downloadFasta <- downloadHandler(
  
  filename = function() {
    
    paste0(species_db_name, '.fasta', sep='')
    
  },
  
  content = function(file) {
    
    file.copy(paste0('/home/biognosis/shiny-server/GAP/db/fasta/', species_db_name, '.fasta'), file)
    
  }
  
)

# confirm ----
observeEvent(input$confirm_ident1, {
  
  if (trimws(input$organ_info)==""|trimws(input$email)=="") {
    
    showModal(modalDialog(
      title = "Important message",
      paste0("Please fill in blank!"),
      easyClose = T
    ))
    
  }else{
    
    process_item<-c()
    
    for (i in c("ident_stat_check")) {
      
      if (input[[i]]) {
        
        process_item<-c(process_item,i)
        
      }
      
    }
    
    project_id <- generate_jobid()
    
    path_prefix <<- paste0('/home/biognosis/mnt/os2/BioAnalysis/')
    
    data_path <<- paste0(path_prefix,'data/usrdata/', project_id, '/data')
    
    para <<- list();para[['info']]<<-list(para[['info']]);para[['condition']]<<-list(para[['condition']])
    
    para[['info']][['type']] <<- 'pre_process';
    
    para[['info']][['data_path']] <<- data_path
    
    dir.create(path = para[['info']][['data_path']], recursive = TRUE)
    
    para[['info']][['project_id']] <<- project_id
    
    para[['info']][['db_sci_name']] <<- species_names
    
    para[['info']][['glycoType']] <<- input$glycoType
    
    para[['info']][['project_name']] <<- input$project_name
    
    para[['info']][['user']] <<- input$organ_info
    
    para[["info"]][["email"]] <<- input$email
    
    para[["info"]][["upload_fasta"]] <<- input$upload_fasta
    
    para[['info']][['db_name']] <<- species_db_name
    
    para[['info']][['ch_names']] <<- species_ch_names
    
    para[['info']][['process_item']] <<- process_item
    
    ## 1.  原始数据PCA
    para[['condition']][['miss.value.handle.method']] <<- input$miss.value.handle.method
    
    para[['condition']][['miss.value.handle.show']] <<- input$miss.value.handle.show
    
    ## 2.  缺失值过滤
    para[['condition']][['miss.value.handle.batch']] <<- input$miss.value.handle.batch
    
    para[['condition']][['miss.value.handle.group']] <<- input$miss.value.handle.group
    
    para[['condition']][['miss.value.handle.cutoff']] <<- input$miss.value.handle.cutoff
    
    if (!is.null(input$miss.value.handle.condition)) {
      
      para[['condition']][['miss.value.handle.condition']] <<- input$miss.value.handle.condition
      
    }
    
    ## 3. 剩余缺失值处理方法
    para[['condition']][['miss.value.handle.remain']] <<- input$miss.value.handle.remain
    
    ## 4. 填充方式
    para[['condition']][['miss.value.fill']] <<- input$miss.value.fill
    
    if (!is.null(input$miss.value.fill.defined)) {
      
      para[['condition']][['miss.value.fill.defined']] <<- input$miss.value.fill.defined
      
    }
    
    ### 如果有效数据占比< 50%
    if (!is.null(input$miss.value.valid.below50)) {
      
      para[['condition']][['miss.value.valid.below50']] <<- input$miss.value.valid.below50
      
    }
    
    if (!is.null(input$miss.value.valid.below50.defined)) {
      
      para[['condition']][['miss.value.valid.below50.defined']] <<- input$miss.value.valid.below50.defined
      
    }
    
    ### 如果有效数据占比 ≥ 50%（自动排除上一选择）
    if (!is.null(input$miss.value.valid.over50)) {
      
      para[['condition']][['miss.value.valid.over50']] <<- input$miss.value.valid.over50
      
    }
    
    if (!is.null(input$miss.value.valid.over50.defined)) {
      
      para[['condition']][['miss.value.valid.over50.defined']] <<- input$miss.value.valid.over50.defined
      
    }
    
    ## 5. 中位数处理
    para[['condition']][['normalized.handle.method']] <<- input$normalized.handle.method
    
    ## 6. 批次矫正
    if (!is.null(input$batch.correct)) {
      
      para[['condition']][['batch.correct']] <<- input$batch.correct
      
    }
    
    if (!is.null(input$batch.correct.reference.sample)) {
      
      para[['condition']][['batch.correct.reference.sample']] <<- input$batch.correct.reference.sample
      
    }
    
    ## 7. 处理后的数据PCA
    para[['condition']][['pca.handle.method']] <<- input$pca.handle.method
    
    para[['condition']][['pca.handle.show']] <<- input$pca.handle.show
    
    # ----
    file.copy(from = input[["glycoRaw"]][['datapath']], to = paste0(para[['info']][['data_path']], '/data.zip'))
    
    para[['info']][['data_type']]<<-"zip"
    
    system(paste0('unzip ', data_path,'/data.zip -d ', data_path))
    
    save(list = c('para'), file = paste0(para[['info']][['data_path']], '/para.Rdata'), envir = .GlobalEnv)
    
    system(paste0("cp -rf /home/biognosis/mnt/os2/BioAnalysis/db/projects_list.db /home/biognosis/shiny-server/BioAnalysis/db/projects_list.db"))
    
    projects_list <- dbConnect(RSQLite::SQLite(), "/home/biognosis/shiny-server/BioAnalysis/db/projects_list.db", flags = SQLITE_RWC)
    
    dbSendQuery(projects_list, paste0("INSERT INTO ",dbListTables(projects_list)," VALUES ('", para[['info']]['project_id'], "','", para[['info']]['project_name'], "', '", input$organ_info, "', '1', '", as.character(Sys.time()),"', 'pre_process', '",para[["info"]]["email"],"')"))
    
    dbDisconnect(projects_list)
    
    system(paste0("cp -rf /home/biognosis/shiny-server/BioAnalysis/db/projects_list.db /home/biognosis/mnt/os2/BioAnalysis/db/projects_list.db"))
    
    shinyalert::shinyalert(
      title = "Important message",
      text=paste0("Your task is created. We will informed your email when it finished!"),
      type="success",
      callbackR = function(x) { if(x != FALSE)
        runjs(paste0(
          "(function() {
		location.href='';
		}());"
        ))
      }
    )
    
    paste0("'/usr/local/lib/R/bin/Rscript' '", path_prefix, "/etc/dep/bg_script.r' ",para[['info']][['project_id']]) %>% system(wait = F)
    
  }
  
})