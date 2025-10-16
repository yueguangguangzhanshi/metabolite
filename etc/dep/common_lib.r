#common

#点击左侧边展开
observeEvent(input$lsb_status, {
  
  js$getsidebarclass()
  
  if(input$tabs %in% c("overview")){
    
    if(grep(input$sidebarclass, pattern = 'control-sidebar-open') %>% length > 0){
      
      # runjs(
      # "$('#my_maintab').animate({width: '70%'}, 300);
      # setTimeout(function(){
      # $('#my_widget').slideDown('fast').css('display','inline-block');
      # },350)"
      # )
      
    } else {
      
      runjs(
        "
				setTimeout(function(){
					$('#my_maintab').animate({width: '70%'}, 200);;
				},180)"
      )
      
    }
    
  }
  
})

observeEvent(input$selected_species, {
  
  species_sci_name <<- input$selected_species %>% str_replace_all('( \\().*', '')
  
  species_index <<- which(species_list$sci_names == species_sci_name)
  
  species_db_name <<- species_list$name[species_index]
  
  species_ch_names <<- species_list$ch_names[species_index]
  
  species_kegg <<- species_list$kegg[species_index]
  
  species_taxid <<- species_list$taxid[species_index]
  
  species_updatetime <<- species_list$updatetime[species_index]
  
})

output$downloadFasta <- downloadHandler(
  
  filename = function() {
    
    paste0(species_db_name, '.fasta', sep='')
    
  },
  
  content = function(file) {
    
    file.copy(paste0('../db/fasta/', species_db_name, '.fasta'), file)
    
  }
  
)

output$dl_report <- downloadHandler(
  
  filename = function() {
    
    'report.zip'
    
  },
  
  content = function(file) {
    
    file.copy(paste0('../data/usrdata/', DataRow_report(), '/', 'report.zip'), file)
    
  }
  
)

output$match_panel <- renderUI({
  
  box(title = 'Spectrum Match', solidHeader = TRUE, status = "primary", width = 12, style = 'margin: 0 auto; min-height: 500px;',
      
      h3(tags$a("pGlyco3",href="https://github.com/pFindStudio/pGlyco3/releases",target="_blank")),
      h4("pGlyco is a software tool designed for the analysis of intact glycopeptides by using mass spectrometry.\n\n"),
      h3(tags$a("pFind",href="http://pfind.ict.ac.cn/software/pFind/index.html",target="_blank")),
      h4("pFind® 3 is developed as an upgrade of pFind 2.8 .\n

Shotgun proteomics has grown rapidly in recent decades, but a large fraction of tandem mass spectrometry data in shotgun proteomics are not successfully identified. Thus we developed a novel database search algorithm, Open-pFind, to efficiently identify peptides even in an ultra-large search space which takes into account unexpected modifications, amino acid mutations, semi- or non-specific digestion and co-eluting peptides. Open-pFind has now been integrated into pFind 3 as the default workflow, and pFind 3 also supports the restricted search mode, which can be switched in the interface.\n\n")
      
  )
  
})

SelectedRow_summary <- eventReactive(input$summary_click,{
  
  index_temp <- as.numeric(strsplit(input$summary_click, "_")[[1]][2])
  
  return(index_temp)
  
})

SelectedRow_pep <- eventReactive(input$pep_click,{
  
  index_temp <- as.numeric(strsplit(input$pep_click, "_")[[1]][2])
  
  return(index_temp)
  
})

SelectedRow_report <- eventReactive(input$report_click,{
  
  index_temp <- as.numeric(strsplit(input$report_click, "_")[[1]][2])
  
  return(index_temp)
  
})

SelectedRow_download <- eventReactive(input$download_click,{
  
  index_temp <- as.numeric(strsplit(input$download_click, "_")[[1]][2])
  
  return(index_temp)
  
})

SelectedRow_manage <- eventReactive(input$manage_click,{
  
  index_temp <- as.numeric(strsplit(input$manage_click, "_")[[1]][2])
  
  return(index_temp)
  
})

SelectedRow_msg <- eventReactive(input$msg_click,{
  
  index_temp <- as.numeric(strsplit(input$msg_click, "_")[[1]][2])
  
  return(index_temp)
  
})

DataRow_summary <- eventReactive(input$summary_click,{
  
  return(generate_project_status_table_output(pri_level, username, status_box)[SelectedRow_summary(), 'id'])
  
})

DataRow_report <- eventReactive(input$report_click,{
  
  return(generate_project_status_table_output(pri_level, username, status_box)[SelectedRow_report(), 'id'])
  
})

DataRow_manage <- eventReactive(input$manage_click,{
  
  return(generate_project_status_table_output(pri_level, username, status_box)[SelectedRow_manage(), 'id'])
  
})

DataRow_pep <- eventReactive(input$pep_click,{
  
  return(generate_project_status_table_output(pri_level, username, status_box)[SelectedRow_pep(), 'id'])
  
})

DataRow_download <- eventReactive(input$download_click,{
  
  Row_id <- generate_project_status_table_output(pri_level, username, status_box)[SelectedRow_download(), 'id']
  
  return(Row_id)
  
})

#关闭modal后清除data并且重置input$select_button的值

observeEvent(input$last_modal_close, {
  
  js$resetClicksummary()
  
  js$resetClickmanage()
  
  js$resetClickmsg()
  
  js$resetClickreport()
  
  js$resetClickspdetail()
  
  runjs(
    "x = new Date().getTime().toLocaleString(); Shiny.onInputChange('refresh_notify', x);"
  )
  
})

observeEvent(input$spdetail_click, {
  
  if (input$spdetail_click != 'null'){
    
    toggleModal(session, "modal_spdetail", "open")
    
  }
  
})

observeEvent(input$summary_click, {
  
  if (input$summary_click != 'null'){
    
    toggleModal(session, "modal_summary", "open")
    
  }
  
})

observeEvent(input$report_click, {
  
  if (input$report_click != 'null'){
    
    toggleModal(session, "modal_report", "open")
    
    # output$tablerReport<-renderUI(tablerDash::tablerCard(
    #   # title = "REPORT",
    #   width = 8,
    #   overflow = T,
    #   shiny::downloadLink(outputId = 'dl_report', label = shinysky::actionButton(inputId = 'dl_btn_report', label = 'DOWNLOAD', styleclass = 'primary', style = '') %>% paste0 %>% HTML()),
    #   
    #   generate_report_content(DataRow_report())
    # ))
    
  }
  
})

observeEvent(input$manage_click, {
  
  if (input$manage_click != 'null'){
    
    toggleModal(session, "modal_manage", "open")
    
  }
  
})

observeEvent(input$msg_click, {
  
  if (input$msg_click != 'null'){
    
    toggleModal(session, "modal_msg", "open")
    
  }
  
})

observeEvent(input$create_click, {
  
  if (input$create_click != 'null'){
    
    toggleModal(session, "modal_create", "open")
    
  }
  
})

output$popup_spdetail <- renderUI({
  
  shiny::req(input$selected_species)
  
  shinyBS::bsModal("modal_spdetail", title = "DATABASE INFO", trigger = "spdetail_click", size = "small",
                   
                   generate_spdetail_content(input$selected_species, species_list)
                   
  )
  
})

output$popup_msg <- renderUI({
  
  shinyBS::bsModal("modal_msg", title = "MESSAGE", trigger = "msg_click", size = "middle",
                   
                   generate_msg_content(row_msg = SelectedRow_msg(), msg_am = msg_am)
                   
  )
  
})

output$popup_summary <- renderUI({
  
  shinyBS::bsModal("modal_summary", title = "SUMMARY", trigger = "summary_click", size = "middle",
                   
                   generate_summary_content(DataRow_summary(), pri_level, username)
                   
  )
  
})

output$popup_report <- renderUI({
  
  if(DataRow_report() %>% is.na == FALSE){
    
    shiny::addResourcePath(prefix = paste0(DataRow_report(), '_report_dir'), directoryPath = paste0(path_prefix, '/data/usrdata/', DataRow_report()))
    
    shinyBS::bsModal("modal_report", title = "REPORT", trigger = "report_click", size = "large",

                     shiny::downloadLink(outputId = 'dl_report', label = shinysky::actionButton(inputId = 'dl_btn_report', label = 'DOWNLOAD', styleclass = 'primary', style = '') %>% paste0 %>% HTML()),

                     generate_report_content(DataRow_report())

    )
    
    # uiOutput("tablerReport")

  }
  
})

output$popup_manage <- renderUI({
  
  shinyBS::bsModal("modal_manage", title = "MANAGEMENT", trigger = "manage_click", size = "large",
                   
                   generate_manage_content(DataRow_manage())
                   
  )
  
})

output$project_table_panel <- renderUI({
  
  check <- input$tabs
  
  shiny::div(style="font-size: 14px;", dataTableOutput(outputId="project_table"), uiOutput("hidden_downloads"))
  
})

output$project_table <- {DT::renderDataTable(
  
  generate_project_status_table_output(pri_level, username, status_box) %>% select(., -id) %>% arrange(., desc(DATE)),
  
  escape = FALSE,
  
  extensions = c('Buttons'),
  
  options = list(
    
    dom = 'Bftp',
    
    scroller = TRUE,
    
    #paging = FALSE,
    
    rowCallback = JS("function(r,d) {$(r).attr('height', '40px')}"),
    
    columnDefs = list(list(className = 'dt-center', targets = c(2,3)), list(className = 'dt-head-center', targets = 4)),
    
    buttons = list(
      
      list(extend = 'collection', text = 'Start New', className = 'myDTbutton',
           
           buttons = project_tb_btn
           
      )
      
    )
    
  ),
  
  selection = "none"
  
)}

output$usertext <- renderText({
  
  username
  
})

observeEvent(input$refresh_notify, {
  
  msg <<- get_msg(pri_level, username)
  
  msg_am <<- rbind(msg$anc, msg$msg)
  
  if(nrow(msg$notif) > 0){
    
    for(notif_index in 1:nrow(msg$notif)){
      
      shiny::showNotification(paste0(msg$notif[notif_index, 'content']), id = paste0('mynotify_', notif_index), type = "default", duration = NULL, closeButton = FALSE, action = shiny::a(style = 'text-decoration: none; font-weight: 400; color: #888888; margin-top: 1.5em;', href = "#", onmousedown = paste0("Shiny.onInputChange('close_notif', ", notif_index, ");"), "Close"))
      
    }
    
  }
  
})

observeEvent(input$close_notif, {
  
  removeNotification(paste0('mynotify_', input$close_notif))
  
  #标记已读
  
  message_db <- dbConnect(RSQLite::SQLite(), paste0(path_prefix, "/db/msg.db"))
  
  dbSendStatement(message_db, paste0("UPDATE Message SET rdtag = 1 WHERE msg_index = '", msg$notif[input$close_notif, 'msg_index'] %>% as.numeric, "'"))
  
  dbDisconnect(message_db)
  
})

output$server_announcement <- renderUI({
  
  refresh <- input$refresh_notify
  
  #1. 全局通知(anc) 2. 一次性提醒(notif) 3. 一般消息(msg)
  
  if(msg$anc %>% nrow >= 1){
    
    box(title = "Announcement",  solidHeader = TRUE, status = "primary", width = 12, style = 'font-weight: 500',
        
        lapply(1:nrow(msg$anc), FUN = generate_msg_button, msg = msg_am) %>% unlist %>% paste0("<div style = 'margin-top: 0.3em; max-height: 10em; '>", ., '</div>', collapse = '') %>% shiny::HTML() %>% div(style = 'max-height:100px; overflow-y: auto;')
        
    )
    
  } else {
    
    box(title = "Announcement",  solidHeader = TRUE, status = "primary", width = 12, style = 'font-weight: 500',
        
        'No announcement yet.'
        
    )
    
  }
  
})

output$server_message <- renderUI({
  
  refresh <- input$refresh_notify
  
  #1. 全局通知(anc) 2. 一次性提醒(notif) 3. 一般消息(msg)
  
  if(msg$msg %>% nrow >= 1){
    
    box(title = "Message",  solidHeader = TRUE, status = "primary", width = 12, style = 'font-weight: 500',
        
        lapply((nrow(msg$anc) + 1):nrow(msg_am), FUN = generate_msg_button, msg = msg_am) %>% unlist %>% paste0("<div style = 'margin-top: 0.3em; max-height: 10em; '>", ., '</div>', collapse = '') %>% shiny::HTML() %>% div(style = 'max-height:100px; overflow-y: auto;')
        
    )
    
  } else {
    
    box(title = "Message",  solidHeader = TRUE, status = "primary", width = 12, style = 'font-weight: 500',
        
        'No message yet.'
        
    )
    
  }
  
})

get_controlbar <- function(pri_level, condition, extra_vars = ''){
  
  if('Master' %in% pri_conf[[pri_level %>% as.character]] | 'Admin' %in% pri_conf[[pri_level %>% as.character]]){
    
    if(file.exists(paste0(path_prefix, '/.stopped'))){
      
      init_css <- ""
      
      stop_css <- "non-active"
      
      status_str <- 'Stopped'
      
      status_color <- '#FF9933'
      
      
    } else {
      
      init_css <- "non-active"
      
      stop_css <- ""
      
      status_str <- 'Running'
      
      status_color <- '#009933'
      
    }
    
    if(condition == 'overview') {
      
      output$my_controlbar <- renderUI({
        
        #master, admin
        
        controlbarMenu(
          
          id = "control_menu_id",
          
          controlbarItem(
            
            "Config",
            
            shiny::div(
              
              shiny::div(style = 'display: block; margin-top: 1em; margin-left: 1.5em',
                         
                         selectInput(inputId = 'config_theme', choices = c('default'), selected = 'default', label = 'Plot Theme:', width = '85%'),
                         
                         shinysky::actionButton(inputId = 'apply_config', label = 'Apply', styleclass = 'info', style = 'margin-right: 0.5em;', size = "mini")
                         
              )
              
            )
            
          ),
          
          controlbarItem(
            
            "Server",
            
            shiny::div(
              
              shiny::div(style = 'display: inline-block; margin-top: 1em;',
                         
                         shiny::a(style = 'margin-left: 1.3em; font-weight: bold; color: #FFFFFF;', "Background Status:"),
                         
                         shiny::a(style = paste0('text-align: left; font-weight: bold; color: ', status_color, ';'), status_str)
                         
              ),
              
              shiny::div(style = 'text-align: center; display: block; margin-top: 1em;',
                         
                         shinysky::actionButton(inputId = 'init_bgserve', label = 'Start', styleclass = 'info', style = 'margin-right: 0.5em;', class = init_css, size = "mini"),
                         
                         shinysky::actionButton(inputId = 'stop_bgserve', label = 'Stop', styleclass = 'info', style = '', class = stop_css, size = "mini")
                         
              )
              
            )
            
          )
          
        )
        
      })
      
    } else if(condition == 'quant') {
      
      rsb_quant_peptide_method_class <- 'myhide'
      
      rsb_quant_p_convert_class <- 'myhide'
      
      if(extra_vars$custom_pvalue == TRUE | is.null(extra_vars$custom_pvalue)){
        
        rsb_quant_p_convert_class <- 'myshow'
        
      }
      
      #肽段定量及仅肽段时更新test
      
      if(extra_vars$peptide_quant_check == 'peptide_extract'){
        
        rsb_quant_peptide_method_class <- 'myshow'
        
      }
      
      if(extra_vars$peptide_quant_check == 'peptide_quant'){
        
        rsb_quant_peptide_method_class <- 'myshow'
        
      }
      
      output$my_controlbar <- renderUI({
        
        #master, admin
        
        controlbarMenu(
          
          id = "control_menu_id",
          
          controlbarItem(
            
            "Method",
            
            shiny::div(
              
              shiny::div(style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em',
                         
                         selectInput(inputId = 'fill_method', label = "NA filling Method", choices = c('default', 'Fill 0', 'KNN'), selected = 'default')
                         
              ),
              
              shiny::div(id = 'rsb_quant_test_method', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em',
                         
                         selectInput(inputId = 'test_method', label = "Test Method:", choices = c('One-way ANOVA', 'T-Test'), selected = 'One-way ANOVA')
                         
              ),
              
              shiny::div(id = 'rsb_quant_p_convert', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em', class = rsb_quant_p_convert_class,
                         
                         selectInput(inputId = 'p_convert', label = "P-Value Type:", choices = c('Significance', 'P-Value'), selected = 'P-Value')
                         
              ),
              
              shiny::div(id = 'rsb_quant_peptide_method', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em', class = rsb_quant_peptide_method_class,
                         
                         selectInput(inputId = 'pep_method', label = "Peptide Quantative Method：", choices = c('top-3 peptide', 'Mean'), selected = 'top-3 peptide')
                         
              )
              
            )
            
          ),
          
          controlbarItem(
            
            "Config",
            
            shiny::div(
              
              shiny::div(style = 'display: block; margin-top: 1em; margin-left: 1.5em',
                         
                         selectInput(inputId = 'config_theme', choices = c('default'), selected = 'default', label = 'Plot Theme:', width = '85%'),
                         
                         shinysky::actionButton(inputId = 'apply_config', label = 'Apply', styleclass = 'info', style = 'margin-right: 0.5em;', size = "mini")
                         
              )
              
            )
            
          ),
          
          controlbarItem(
            
            "Server",
            
            shiny::div(
              
              shiny::div(style = 'display: inline-block; margin-top: 1em;',
                         
                         shiny::a(style = 'margin-left: 1.3em; font-weight: bold; color: #FFFFFF;', "Background Status:"),
                         
                         shiny::a(style = paste0('text-align: left; font-weight: bold; color: ', status_color, ';'), status_str)
                         
              ),
              
              shiny::div(style = 'text-align: center; display: block; margin-top: 1em;',
                         
                         shinysky::actionButton(inputId = 'init_bgserve', label = 'Start', styleclass = 'info', style = 'margin-right: 0.5em;', class = init_css, size = "mini"),
                         
                         shinysky::actionButton(inputId = 'stop_bgserve', label = 'Stop', styleclass = 'info', style = '', class = stop_css, size = "mini")
                         
              )
              
            )
            
          )
          
        )
        
      })
      
    }
    
    
  } else {
    
    if(condition == 'overview') {
      
      output$my_controlbar <- renderUI({
        
        #other users
        
        controlbarMenu(
          
          id = "control_menu_id",
          
          controlbarItem(
            
            "Config",
            
            shiny::div(
              
              selectInput(inputId = 'config_theme', choices = c('default'), selected = 'default', label = 'Plot Theme:'),
              
              shinysky::actionButton(inputId = 'apply_config', label = 'Apply', styleclass = 'info', style = 'margin-right: 0.5em;')
              
            )
            
          )
          
        )
        
      })
      
    } else if(condition == 'quant') {
      
      rsb_quant_peptide_method_class <- 'myhide'
      
      rsb_quant_p_convert_class <- 'myhide'
      
      if(extra_vars$custom_pvalue == TRUE | is.null(extra_vars$custom_pvalue)){
        
        rsb_quant_p_convert_class <- 'myshow'
        
      }
      
      #肽段定量及仅肽段时更新test
      
      if(extra_vars$peptide_quant_check == 'peptide_extract'){
        
        rsb_quant_peptide_method_class <- 'myshow'
        
      }
      
      if(extra_vars$peptide_quant_check == 'peptide_quant'){
        
        rsb_quant_peptide_method_class <- 'myshow'
        
      }
      
      output$my_controlbar <- renderUI({
        
        #other users
        
        controlbarMenu(
          
          id = "control_menu_id",
          
          controlbarItem(
            
            "Method",
            
            shiny::div(
              
              shiny::div(style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em',
                         
                         selectInput(inputId = 'fill_method', label = "NA filling Method", choices = c('default', 'Fill 0', 'KNN'), selected = 'default')
                         
              ),
              
              shiny::div(id = 'rsb_quant_test_method', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em',
                         
                         selectInput(inputId = 'test_method', label = "Test Method:", choices = c('One-way ANOVA', 'T-Test'), selected = 'One-way ANOVA')
                         
              ),
              
              shiny::div(id = 'rsb_quant_p_convert', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em', class = rsb_quant_p_convert_class,
                         
                         selectInput(inputId = 'p_convert', label = "P-Value Type:", choices = c('Significance', 'P-Value'), selected = 'P-Value')
                         
              ),
              
              shiny::div(id = 'rsb_quant_peptide_method', style = 'display: block; margin-top: 1em; margin-left: 1.5em; margin-right: 1.5em', class = rsb_quant_peptide_method_class,
                         
                         selectInput(inputId = 'pep_method', label = "Peptide Quantative Method：", choices = c('top-3 peptide', 'Mean'), selected = 'top-3 peptide')
                         
              )
              
            )
            
          ),
          
          controlbarItem(
            
            "Config",
            
            shiny::div(
              
              selectInput(inputId = 'config_theme', choices = c('default'), selected = 'default', label = 'Plot Theme:'),
              
              shinysky::actionButton(inputId = 'apply_config', label = 'Apply', styleclass = 'info', style = 'margin-right: 0.5em;')
              
            )
            
          )
          
        )
        
      })
      
    }
    
  }
  
}