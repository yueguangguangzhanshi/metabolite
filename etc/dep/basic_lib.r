shinyInput <- function(FUN, len, id, ...) {
  inputs <- character(len)
  for (i in seq_len(len)) {
    inputs[i] <- as.character(FUN(paste0(id, i), ...))}
  inputs
}

path_prefix <- '/home/biognosis/shiny-server/GAP'

username<-"GAP"

makeElement <- function(data, name)
{
  shiny::div(style = "margin: 0.5em; border-width:1px;border-style:solid;color:grey;border-radius:2px;min-height: 50px;",
             drag = name,
             shiny::div(class = "active", style = "margin: 0.5em auto 0 auto; font-weight:bold; font-size: 12px; color:black; text-align:center", name),
             shiny::div(class = "active", style = "margin: 0 auto 0.5em auto;font-size: 10px; color:black; text-align:center", sprintf("Class: %s", class(data[[name]])))
  )
}

generate_jobid <- function() {
  
  a <- paste0(collapse = '', sample(x = c(letters, LETTERS, 0:9), size = 12, replace = TRUE))
  
  paste0(round(as.numeric(Sys.time()),0), a)
  
}

grey_screen <- function (status){
  
  if(status == 'on'){
    
    shinyjs::removeCssClass(id = "loading-content", class = "hide")
    
    shinyjs::addCssClass(id = "loading-content", class = "show")
    
  } else {
    
    shinyjs::removeCssClass(id = "loading-content", class = "show")
    
    shinyjs::addCssClass(id = "loading-content", class = "hide")
    
  }
  
}

myhide <- function (id){
  
  for(i in id){
    
    shinyjs::removeCssClass(id = i, class = "myshow")
    
    shinyjs::addCssClass(id = i, class = "myhide")
    
  }
  
}

myshow <- function (id){
  
  for(i in id){
    
    shinyjs::removeCssClass(id = i, class = "myhide")
    
    shinyjs::addCssClass(id = i, class = "myshow")
    
  }
  
  
}

pri_function <- c("Spectrum identification", "Glycosylation analysis", "Quantitative analysis", "Functional analysis", "Clinical analysis")

jsCode <- '
  shinyjs.getcookie = function(params) {
    var cookie = Cookies.get("uid");
    if (typeof cookie !== "undefined") {
      Shiny.onInputChange("jscookie", cookie);
    } else {
      var cookie = "";
      Shiny.onInputChange("jscookie", cookie);
    }
  }
  shinyjs.setcookie = function(params) {
    Cookies.set("uid", escape(params), { expires: 7 });  
    Shiny.onInputChange("jscookie", params);
  }
  shinyjs.rmcookie = function(params) {
    Cookies.remove("uid");
    Shiny.onInputChange("jscookie", "");
  }
'

status_box <- list(
  
  "code" = 0:4,
  
  "status_color" = c("info", "info", "success", "warning", "danger"),
  
  "bkcolor" = c("#3366CC","#3366CC","#006633","#FF9933","#FF0000"),
  
  "status_msg" = c("Waiting","Processing","Finished","Evaluating","Error"),
  
  "status_icon" = c("clock","clock","check-circle","exclamation-circle","times")
  
)

accid_process <- function(id){
  
  return(id %>% unlist %>% unname %>% str_replace_all("(\t).*", "") %>% str_replace_all("[[:punct:]]", "_") %>% str_replace_all("^(tr\\||sp\\|)", "") %>% str_replace_all("\\|.*", ""))
  
}

generate_project_status_table_output <- function(pri_level, username, status_box){
  
  #获取Project Info:
  projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")	
  
  if(pri_level <= 3){
    
    project_status <- dbGetQuery(projects_list, paste0("SELECT * FROM Project")) %>% as.data.frame
    
  } else {
    
    project_status <- dbGetQuery(projects_list, paste0("SELECT * FROM Project WHERE user='",username,"'")) %>% as.data.frame
    
  }
  
  dbDisconnect(projects_list)
  
  if(project_status %>% nrow == 0){
    
    project_status <- data.frame('id' = '', 'PROJECT' = '', 'STATUS' = '', 'DATE' = '', 'OPTION' = '', stringsAsFactors = F) %>% .[-1, ]
    
  } else {
    
    for(i in 1:length(status_box[[1]])){
      
      #project_status[which(project_status$status == status_box$code[i]), '项目名称'] <- lapply(project_status$name[which(project_status$status == status_box$code[i])], shiny::span, style = paste0('font-weight: bold;')) %>% lapply(., paste0, '') %>% paste0(shiny::span(icon(status_box$status_icon[i]), style = paste0('margin-right: 0.5em; font-weight: bold; color:', status_box$bkcolor[i])) %>% paste0, .)
      if(which(project_status$status == status_box$code[i]) %>% length > 0){
        
        project_status[which(project_status$status == status_box$code[i]), 'PROJECT'] <- lapply(project_status$name[which(project_status$status == status_box$code[i])], shiny::span, style = paste0('font-weight: bold;')) %>% lapply(., paste0, '') %>% paste0(shiny::span(style = paste0('margin-right: 0.5em; font-weight: bold; color:', status_box$bkcolor[i])) %>% paste0, .)
        
        project_status[which(project_status$status == status_box$code[i]), 'STATUS'] <- status_box$status_msg[i] %>% lapply(FUN = shiny::a, style = paste0('color: #FFFFFF; font-weight:bold; font-size: 10px; text-decoration: none;')) %>% lapply(FUN = shiny::div, style = paste0('background-color:', status_box$bkcolor[i], '; width: 5.5em; height: 1.5em; text-align: center; margin: auto auto;')) %>% lapply(., paste0, '')
        
        project_status[which(project_status$status == status_box$code[i]), 'OPTION'] <- lapply(which(project_status$status == status_box$code[i]), FUN = generate_option_button, project_id = project_status$id, color_i = i, status_box = status_box, pri_level = pri_level, current_window_size = input$current_window_size, project_type = project_status$type) %>% unlist
        
      }
      
    }
    
    project_status[, 'DATE'] <- project_status$time %>% as.Date %>% as.character
    
    if(pri_level <= 3){
      
      project_status[, 'ID'] <- project_status$id
      
      project_status[, 'USER'] <- project_status$user
      
      project_status <- project_status[, c('id', 'PROJECT', 'STATUS', 'USER', 'DATE', 'OPTION')]
      
    } else {
      
      project_status <- project_status[, c('id', 'PROJECT', 'STATUS', 'DATE', 'OPTION')]
      
    }
    
  }
  
  lapply(1:nrow(project_status), function(i) {
    
    output[[paste0('dl_button_', project_status[i, 'id'])]] <- downloadHandler(
      
      filename = function() {
        
        'report.zip'
        
      },
      
      content = function(file) {
        
        file.copy(paste0('../data/usrdata/', project_status[i, 'id'], '/report.zip'), file)
        
      }
      
    )
    
  })
  
  output$hidden_downloads <- renderUI(
    
    lapply(1:nrow(project_status), function(i) {
      
      downloadLink(paste0("dl_button_", project_status[i, 'id']), "download", class = "hiddenLink")
      
    })
    
  )
  
  return(project_status)
  
}

generate_msg_button <- function (row_index, msg){
  
  time_elipse <- (Sys.time() %>% as.numeric - msg[row_index, 'valid_from'] %>% as.numeric)/60
  
  if(time_elipse <= 2) {
    
    msg_time <- paste0(floor(time_elipse), ' minute ago')
    
  }
  
  if(time_elipse > 2 & time_elipse < 60) {
    
    msg_time <- paste0(floor(time_elipse), ' minutes ago')
    
  }
  
  if(time_elipse >= 60 & time_elipse < 60*2) {
    
    msg_time <- paste0(floor(time_elipse/60), ' hour ago')
    
  }
  
  if(time_elipse >= 60*2 & time_elipse < 60*24) {
    
    msg_time <- paste0(floor(time_elipse/60), ' hours ago')
    
  }
  
  if(time_elipse >= 60*24 & time_elipse < 60*24*2) {
    
    msg_time <- paste0(floor(time_elipse/(60*24)), ' day ago')
    
  }
  
  if(time_elipse >= 60*24*2 & time_elipse < 60*24*7) {
    
    msg_time <- paste0(floor(time_elipse/(60*24)), ' days ago')
    
  }
  
  if(time_elipse >= 60*24*7 & time_elipse < 60*24*7*2) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7)), ' week ago')
    
  }
  
  if(time_elipse >= 60*24*7*2 & time_elipse < 60*24*7*4) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7)), ' weeks ago')
    
  }
  
  if(time_elipse >= 60*24*7*4 & time_elipse < 60*24*7*4*2) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7*4)), ' month ago')
    
  }
  
  if(time_elipse >= 60*24*7*4*2 & time_elipse < 60*24*7*4*12) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7*4)), ' months ago')
    
  }
  
  if(time_elipse >= 60*24*7*4*12 & time_elipse < 60*24*7*4*12*2) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7*4*12)), ' year ago')
    
  }
  
  if(time_elipse >= 60*24*7*4*12*2) {
    
    msg_time <- paste0(floor(time_elipse/(60*24*7*4*12)), ' years ago')
    
  }
  
  if(msg[row_index, 'typ'] == 3){
    
    if(msg[row_index, 'subject'] %>% str_length > 5){
      
      msg_title <- msg[row_index, 'subject'] %>% substr(start = 1, stop = 5) %>% paste0(., '...')
      
    } else {
      
      msg_title <- msg[row_index, 'subject']
      
    }
    
  } else {
    
    if(msg[row_index, 'subject'] %>% str_length > 8){
      
      msg_title <- msg[row_index, 'subject'] %>% substr(start = 1, stop = 8) %>% paste0(., '...')
      
    } else {
      
      msg_title <- msg[row_index, 'subject']
      
    }
    
  }
  
  
  # msg-database subject, content, sender, typ, rcvr, rdtag, rcpt_group,  valid_from, valid_thru
  #type: 1. 公告 2. 提醒 3.消息
  
  if(msg[row_index, 'typ'] == 3 & msg[row_index, 'rdtag'] == 0){
    
    paste0(
      
      "<div style = 'display:inline-block; vertical-align: top; margin-bottom: 0.3em; font-style: italic; color:#999999;'>",
      
      msg_time,
      
      "</div><div  style = 'color:#999999; margin-right:0.5em; width: 6em; float:left'>",
      
      shiny::a(href = '#', id = paste0('button_', row_index), style = paste0('margin-left:0.5em; color:#555555; text-decoration: none;'), onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"msg_click\",  this.id);', msg_title) %>% paste0,
      
      "</div><span style = 'color: #777777; font-weight: 600; display:-moz-inline-box; display:inline-block; width: 40px; margin-left: 1em; font-style: italic;'>NEW</span></div><div style = 'clear:both;'><div style = 'clear:both;'>",
      
      collapse = ''
      
    ) %>% return
    
    
  } else if(msg[row_index, 'typ'] == 3 & msg[row_index, 'rdtag'] == 1){
    
    paste0(
      
      "<div style = 'display:inline-block; vertical-align: top; margin-bottom: 0.3em; font-style: italic; color:#999999;'>",
      
      msg_time,
      
      "</div><div  style = 'color:#999999; margin-right:1.5em; width: 6em; float:left'>",
      
      shiny::a(href = '#', id = paste0('button_', row_index), style = paste0('margin-left:0.5em; color:#555555; text-decoration: none;'), onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"msg_click\",  this.id);', msg_title) %>% paste0,
      
      "</div><div style = 'clear:both;'>",
      
      collapse = ''
      
    ) %>% return
    
  } else {
    
    paste0(
      
      "<div style = 'display:inline-block; vertical-align: top; margin-bottom: 0.3em; '>",
      
      shiny::a(href = '#', id = paste0('button_', row_index), style = paste0('margin-left:0.5em; color:#555555; text-decoration: none;'), onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"msg_click\",  this.id);', msg_title) %>% paste0,
      
      "</div><div  style = 'color:#999999; font-style: italic; margin-right:0.5em; width: 8em; float:right'>",
      
      msg_time,
      
      "</div><div style = 'clear:both;'>",
      
      collapse = ''
      
    ) %>% return
    
  }
  
  
}

#载入权限设置

if(file.exists(paste0(path_prefix, "/etc/config/pri_level.csv"))){
  
  pri_table <- read_csv(paste0(path_prefix, "/etc/config/pri_level.csv"))
  
  pri_conf <- as.list(pri_table$`function` %>% sapply(., strsplit, split = ','))
  
  names(pri_conf) <- pri_table$pri_level %>% unlist
  
}

get_project_tb_btn <- function(pri_level){
  
  for(i in 1:length(pri_conf)){
    
    if(pri_level == names(pri_conf)[i]){
      
      if('guest' %in% pri_conf[[i]]){
        
        project_tb_btn <- list(
          
          list(extend = 'collection', text = 'Please Login', className = 'myDTbutton',
               
               action = DT::JS("")
               
          )
          
        )
        
      } else {
        
        project_tb_btn <- list()
        
        j <- 1
        
        for(k in 1:length(pri_conf[[i]])){
          
          if(pri_conf[[i]][k] %in% pri_function){
            
            project_tb_btn[[j]] <- list(extend = 'collection', text = pri_conf[[i]][k], className = 'myDTbutton',
                                        
                                        action = DT::JS(paste0("function ( e, dt, node, config ) {x = new Date().toLocaleString(); Shiny.setInputValue('start_", pri_conf[[i]][k] %>% tolower(),"', x, {priority: 'event'});}"))
                                        
            )
            
            j = j + 1
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(project_tb_btn)
  
}

get_left_menu <- function(pri_level){
  
  if(pri_level == 98){
    
    left_menu <- sidebarMenu(
      
      id = "tabs",
      
      menuItem("Overview", tabName = "overview", icon = icon("dashboard"), selected=TRUE)
      
      # menuItem("Tools", tabName = "tool_tab", icon = icon("th"))
      
    )
    
  } else if(pri_level <= 3){
    
    left_menu <- sidebarMenu(
      
      id = "tabs",
      
      menuItem("Introduction", tabName = "overview", icon = icon("dashboard"), selected=TRUE),
      
      menuItem("Start Analysis", icon = icon("th"),
               
               menuSubItem("Spectrum identification", tabName = "match_tab", icon = icon("database")),
               
               menuSubItem("Glycosylation analysis", tabName = "ident_tab", icon = icon("database")),
               
               menuSubItem("Quantitative analysis", tabName = "quant_tab", icon = icon("database")),
               
               menuSubItem("Functional analysis", tabName = "function_tab", icon = icon("database")),
               
               menuSubItem("Clinical analysis", tabName = "clinic_tab", icon = icon("database")),
               
               menuSubItem("Databases", tabName = "database_tab", icon = icon("database")),
               
               menuSubItem("Relative Websites", tabName = "rw_tab", icon = icon("database"))
               
      )
      
      # menuItem("Tools", tabName = "tool_tab", icon = icon("th"),
      #          
      #          menuSubItem("Peptide Properties", tabName = "peptides_property_tab", icon = icon("database")),
      #          
      #          menuSubItem("Protein Properties", tabName = "protein_property_tab", icon = icon("database"))
      #          
      # )
      
    )
    
  }
  
  return(left_menu)
  
}

generate_msg_content <- function(row_msg = NULL, msg_am = NULL){
  
  if(row_msg %>% is.null == FALSE & row_msg %>% is.na == FALSE){
    
    if(msg_am[row_msg, 'typ'] == 3){
      
      message_db <- dbConnect(RSQLite::SQLite(), paste0(path_prefix, "/db/msg.db"))
      
      # notif标记已读
      
      dbSendStatement(message_db, paste0("UPDATE Message SET rdtag = 1 WHERE msg_index = '", msg_am[row_msg, 'msg_index'] %>% as.numeric, "'"))
      
      dbDisconnect(message_db)
      
    }
    
  }
  
  paste0(
    
    h4(msg_am[row_msg, 'subject']),
    
    h6(style = 'color: #888888', paste0('Date:', ceiling(msg_am[row_msg, 'valid_from']/(60*60*24)) %>% as.Date(origin = "1970-01-01") %>% as.character())),
    
    hr(style = 'height:1px; border:none; border-top:1px solid #555555; margin-top: 10px;'),
    
    h5(style = 'min-height: 10em', msg_am[row_msg, 'content']),
    
    hr(style = 'height:1px;border:none;border-top:1px solid #555555;')
    
  ) %>% shiny::HTML()
  
}

generate_spdetail_content <- function (species_db_name, species_list){
  
  species_sci_name <<- species_db_name %>% str_replace_all('( \\().*', '')
  
  species_index <- which(species_list$sci_names == species_sci_name)
  
  species_ch_names <- species_list$ch_names[species_index]
  
  species_kegg <- species_list$kegg[species_index]
  
  species_taxid <- species_list$taxid[species_index]
  
  species_updatetime <- ceiling(species_list$updatetime[species_index]/(60*60*24)) %>% as.Date(origin = "1970-01-01") %>% as.character()
  
  paste0(
    
    #shiny::span(icon("info"), href = '#', style = 'margin-left:0.5em; color:#999999; text-decoration: none; display: inline-block; font-size: 14px; width: 2em; margin:0 auto'),
    
    p(paste0('Species：', species_db_name)),
    
    p(paste0('NCBI TAX ID：', species_taxid)),
    
    p(paste0('KEGG：', species_kegg)),
    
    p(paste0('Update Date：', species_updatetime))
    
  ) %>% shiny::HTML() %>% return
  
}

generate_summary_content <- function (project_id = NA, pri_level = NA, username = NA){
  
  if(is.na(project_id) == FALSE){
    
    projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")	
    
    project_table <- dbGetQuery(projects_list, paste0("SELECT * FROM Project WHERE id='", project_id, "'")) %>% as.data.frame
    
    dbDisconnect(projects_list)
    
    load(paste0('../data/usrdata/', project_id, '/data/para.Rdata'))
    
    process_item_word <- list(
      
      'var' = c('go_enrich', 'kegg_enrich', 'string', 'qc'), 
      
      'name' = c("GO Enrichment", 'KEGG Enrichment', 'STRING PPI', 'QC Report')
      
    )
    
    process_item_msg <- process_item_word$name[which(process_item_word$var %in% para[["info"]][["process_item"]])] %>% paste0(collapse = ', ')
    
    project_status <- project_table$status[which(project_table$id == project_id)]
    
    project_status_msg <- status_box$status_msg[which(status_box$code == project_status)]
    
    project_status_color <- status_box$bkcolor[which(status_box$code == project_status)]
    
    if(para[["info"]][["type"]] == 'function'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Overview:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name:', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Species:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(paste0(para[["info"]][["db_sci_name"]], ' (', para[["info"]][["ch_names"]], ')'), style = 'display: inline-block;')
                               
                               #shiny::a(icon("info"), href = '#', style = 'margin-left:0.5em; color:#999999; text-decoration: none; display: inline-block; font-size: 14px; width: 0.5em;', onmousedown = 'event.preventDefault();event.stopPropagation();Shiny.onInputChange(\"spdetail_click\", new Date().toLocaleString());', title = '点击查看详细信息') %>% paste0 %>% HTML()
                               
                           ),
                           
                           p(paste0('(Totally applied:', length(para[["data"]][["protein_list"]]), ' Proteins。)'), style = 'font-size: 12px; color: #555555; font-style: italic;'),
                           
                           p('Analysis:', style = 'font-size: 14px; font-weight: bold; margin-top: 2em'),
                           
                           div(style = 'margin-top: 1em; padding-left: 0;',
                               
                               p(process_item_msg)
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'IPA'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Info", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name:', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'ident'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name: ', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Species: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(paste0(para[["info"]][["db_sci_name"]], ' (', para[["info"]][["ch_names"]], ')'), style = 'display: inline-block;')
                               
                           ),
                           
                           p(paste0('(Totally ', length(para[['data']]), ' group(s).)'), style = 'font-size: 12px; color: #555555; font-style: italic;'),
                           
                           p('Applied Proteins：', style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'margin-top: 0.5em; padding-left: 0;',
                               
                               lapply(seq_len(as.numeric(para[["data"]] %>% length)), function(i) {
                                 
                                 count <- para[["data"]][[i]] %>% unlist %>% length
                                 
                                 group_name <- para[["data"]][['protein_list']][i] %>% names
                                 
                                 if(i == as.numeric(para[["data"]][['protein_list']] %>% length)){
                                   
                                   p('Group：', group_name, count, '.') %>% paste0 %>% shiny::HTML()
                                   
                                 } else {
                                   
                                   p('Group：', group_name, count, ';') %>% paste0 %>% shiny::HTML()
                                   
                                 }
                                 
                               })
                               
                           ),
                           
                           p('Analysis: ', style = 'font-size: 14px; font-weight: bold; margin-top: 2em'),
                           
                           div(style = 'margin-top: 1em; padding-left: 0;',
                               
                               p(process_item_msg)
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'qc_only'){
      
      #save(list = c('project_status_msg', 'para', 'project_table', 'process_item_msg', 'project_id', 'project_status_color'), file = '/home/biognosis/shiny-server/omicscloud/temp.Rdata')
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Report Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Report Name: ', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('QC Method: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(para[["info"]][["qc"]], style = 'display: inline-block;')
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'test_only'){
      
      #save(list = c('project_status_msg', 'para', 'project_table', 'process_item_msg', 'project_id', 'project_status_color'), file = '/home/biognosis/shiny-server/omicscloud/temp.Rdata')
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Report Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p('Test Report', style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p('Finished', style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:;'))
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    }
    
  }
  
}

generate_report_content <- function (project_id){
  
  if(is.na(project_id) == FALSE){
    
    # renderUI({
    #   HTML(read_file(paste0('../data/usrdata/',project_id, '/report.html')))
    # })
    renderUI({

      tags$div(style = '',

               tags$iframe(
                 
                 id = 'reportframe',
                 
                 # onload = "this.height=this.contentWindow.document.body.scrollHeight;onMyFrameLoad(this)",
                 
                 seamless = 'seamless',
                 
                 width = '100%', 
                 
                 height= '1000px',
                 
                 scrolling = 'yes',
                 
                 frameborder = '0',
                 
                 framespacing = '0',
                 
                 src = paste0(project_id, '_report_dir/report.html')
                 
               )

      )

    })
    
  }
  
}

# 'id'
# 'name'
# 'user'
# 'status'
# 'time'
generate_manage_content <- function (project_id){
  
  if(is.na(project_id) == FALSE){
    
    projects_list <- dbConnect(RSQLite::SQLite(), "../db/projects_list.db")	
    
    project_table <- dbGetQuery(projects_list, paste0("SELECT * FROM Project WHERE id='", project_id, "'")) %>% as.data.frame
    
    dbDisconnect(projects_list)
    
    load(paste0('../data/usrdata/', project_id, '/data/para.Rdata'))
    
    process_item_word <- list(
      
      'var' = c('go_enrich', 'kegg_enrich', 'string', 'qc'), 
      
      'name' = c("GO Enrichment", 'KEGG Enrichment', 'STRING PPI', 'QC Report')
      
    )
    
    process_item_msg <- process_item_word$name[which(process_item_word$var %in% para[["info"]][["process_item"]])] %>% paste0(collapse = ', ')
    
    project_status <- project_table$status[which(project_table$id == project_id)]
    
    project_status_msg <- status_box$status_msg[which(status_box$code == project_status)]
    
    project_status_color <- status_box$bkcolor[which(status_box$code == project_status)]
    
    if(para[["info"]][["type"]] == 'function'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Management:", solidHeader = TRUE, status = "primary", width = 3, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name:', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em; margin-left: 1.5em;',
                               
                               shinysky::actionButton(inputId = 'edit_project_name', label = '确认提交', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-bottom: 0; margin-top: 3em;')
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Species:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(paste0(para[["info"]][["db_sci_name"]], ' (', para[["info"]][["ch_names"]], ')'), style = 'display: inline-block;')
                               
                               #shiny::a(icon("info"), href = '#', style = 'margin-left:0.5em; color:#999999; text-decoration: none; display: inline-block; font-size: 14px; width: 0.5em;', onmousedown = 'event.preventDefault();event.stopPropagation();Shiny.onInputChange(\"spdetail_click\", new Date().toLocaleString());', title = '点击查看详细信息') %>% paste0 %>% HTML()
                               
                           ),
                           
                           p(paste0('(Totally applied:', length(para[["data"]][["protein_list"]]), ' Proteins。)'), style = 'font-size: 12px; color: #555555; font-style: italic;'),
                           
                           p('Analysis:', style = 'font-size: 14px; font-weight: bold; margin-top: 2em'),
                           
                           div(style = 'margin-top: 1em; padding-left: 0;',
                               
                               p(process_item_msg)
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'IPA'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Info", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name:', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'ident'){
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Project Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Project Name: ', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Project Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Species: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(paste0(para[["info"]][["db_sci_name"]], ' (', para[["info"]][["ch_names"]], ')'), style = 'display: inline-block;')
                               
                           ),
                           
                           p(paste0('(Totally ', length(para[['data']]), ' group(s).)'), style = 'font-size: 12px; color: #555555; font-style: italic;'),
                           
                           p('Applied Proteins：', style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'margin-top: 0.5em; padding-left: 0;',
                               
                               lapply(seq_len(as.numeric(para[["data"]] %>% length)), function(i) {
                                 
                                 count <- para[["data"]][[i]] %>% unlist %>% length
                                 
                                 group_name <- para[["data"]][['protein_list']][i] %>% names
                                 
                                 if(i == as.numeric(para[["data"]][['protein_list']] %>% length)){
                                   
                                   p('Group：', group_name, count, '.') %>% paste0 %>% shiny::HTML()
                                   
                                 } else {
                                   
                                   p('Group：', group_name, count, ';') %>% paste0 %>% shiny::HTML()
                                   
                                 }
                                 
                               })
                               
                           ),
                           
                           p('Analysis: ', style = 'font-size: 14px; font-weight: bold; margin-top: 2em'),
                           
                           div(style = 'margin-top: 1em; padding-left: 0;',
                               
                               p(process_item_msg)
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'qc_only'){
      
      #save(list = c('project_status_msg', 'para', 'project_table', 'process_item_msg', 'project_id', 'project_status_color'), file = '/home/biognosis/shiny-server/omicscloud/temp.Rdata')
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Report Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p(paste0('Report Name: ', para[["info"]][["project_name"]]), style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(project_status_msg, style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:', project_status_color, ';'))
                               
                           ),
                           
                           br(),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('QC Method: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p(para[["info"]][["qc"]], style = 'display: inline-block;')
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    } else if(para[["info"]][["type"]] == 'test_only'){
      
      #save(list = c('project_status_msg', 'para', 'project_table', 'process_item_msg', 'project_id', 'project_status_color'), file = '/home/biognosis/shiny-server/omicscloud/temp.Rdata')
      
      fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
               
               div(style = 'margin-right: auto; margin-left: auto;',
                   
                   box(title = "Report Info:", solidHeader = TRUE, status = "primary", width = 12, style = '',
                       
                       div(style = 'margin: 3em auto;',
                           
                           p('Test Report', style = 'font-size: 14px; font-weight: bold;'),
                           
                           div(style = 'display: inline-block; margin-top: 1.5em;',
                               
                               p('Status: ', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
                               
                               p('Finished', style = paste0('display: inline-block;; font-size: 14px; font-weight: bold; color:;'))
                               
                           )
                           
                       )
                       
                   )
                   
               )
               
      )
      
    }
    
  }
  
}

get_users <- function (){
  
  usrdb <- dbConnect(RSQLite::SQLite(), "../db/usrdb.db")
  
  usrdb_table <- dbGetQuery(usrdb, "SELECT * FROM User")
  
  dbDisconnect(usrdb)
  
  users <- paste0('uid:', usrdb_table$uid, '; username:', usrdb_table$username, '; email:', usrdb_table$email)
  
  return(users)
  
}

generate_server_status <- function(){
  
  dr <- read_csv('../.runtime')
  
  if(file.exists(paste0(path_prefix, '/.running'))){
    
    server_status <- shiny::span('Running', style = 'font-weight: bold; color: #009933; margin: 0.3em 0 0 0; display:block;')
    
  } else {
    
    server_status <- shiny::span('Maintainance', style = 'font-weight: bold; color: #FF9933; margin: 0.3em 0 0 0; display:block;')
    
  }
  
  if(dr$queue < 5) {
    
    server_load <- shiny::span('Low', style = 'font-weight: bold; color: #009933; font-size: 14px;')
    
  }
  
  if(dr$queue >= 5 & dr$queue < 20) {
    
    server_load <- shiny::span('Mid', style = 'font-weight: bold; color: #3366CC; font-size: 14px;')
    
  }
  
  if(dr$queue >= 20) {
    
    server_load <- shiny::span('High', style = 'font-weight: bold; color: #CC0000; font-size: 14px;')
    
  }
  
  server_status_var <- list()
  
  server_status_var[['HTML']] <- paste0(server_status, shiny::span('Server Load: ', style = 'font-weight: bold; color: #605CA8; font-size: 14px;'), server_load)
  
  server_status_var[['data']] <- dr
  
  return(server_status_var)
  
}

read_numeric_input <- function(x){
  
  if(grep(x, pattern = '-', fixed = TRUE) %>% length > 0){
    
    nums <- x %>% strsplit('-') %>% unlist
    
    if(length(nums) >= 2){
      
      return((nums[1]:nums[2]))
      
    } else {
      
      return(NULL)
      
    }
    
    
  } else if(grep(x, pattern = ',|( )|;|(\\.)') %>% length > 0){
    
    nums <- x %>% strsplit('[[:punct:]]') %>% unlist %>% as.numeric
    
    return(nums)
    
  } else {
    
    if(is.null(x)){
      
      return(NULL)
      
    } else {
      
      return(x)
      
    }
    
  }
  
}

#"status_color" = c("info","success","warning","danger"),

#"bkcolor" = c("#3399CC","#006633","#FF9933","#FF0000"),

generate_option_button <- function (row_index, project_id, color_i, status_box, pri_level, current_window_size, project_type){
  
  pri_level_content <- pri_conf[[pri_level %>% as.character]]
  
  if('Master' %in% pri_level_content | 'Admin' %in% pri_level_content){
    
    if(current_window_size >= 1800){
      
      button_manage <- shiny::tags$button('Manage', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #393A35; font-size:10px; cursor:default; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"manage_click\",  this.id);') %>% paste0
      
    } else {
      
      button_manage <- shiny::tags$button('MNG', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #393A35; font-size:10px; cursor:default; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"manage_click\",  this.id);') %>% paste0
      
    }
    
  } else {
    
    button_manage <- ''
    
  }
  
  if(project_type[row_index] == 'quant_pep'){
    
    if(current_window_size >= 1800){
      
      button_pepview <- shiny::tags$button('Peptide', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #006633; font-size:10px; cursor:default; border-color: #009933; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"pep_click\",  this.id);') %>% paste0
      
    } else {
      
      button_pepview <- shiny::tags$button('PEP', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #006633; font-size:10px; cursor:default; border-color: #009933; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"pep_click\",  this.id);') %>% paste0
      
    }
    
  } else {
    
    button_pepview <- ''
    
  }
  
  if(status_box$status_color[color_i] == 'success'){
    
    if(current_window_size >= 1800){
      
      button_report <- shiny::tags$button('Report', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #006633; font-size:10px; cursor:default; border-color: #009933; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"report_click\",  this.id);') %>% paste0
      
      button_download <- shiny::a(href = '#', onClick = paste0('document.getElementById("dl_button_', project_id[row_index],'").click()'), shiny::tags$button('Download', id = paste0('button_', project_id[row_index]), style = 'color: #FFFFFF; background-color: #393A35; font-size:10px; cursor:default; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'Shiny.onInputChange(\"download_click\",  this.id);')) %>% paste0 %>% HTML()
      
    } else {
      
      button_report <- shiny::tags$button('RPT', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #006633; font-size:10px; cursor:default; border-color: #009933; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = paste0('event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"report_click\",  this.id);')) %>% paste0
      
      button_download <- shiny::a(href = '#', onClick = paste0('document.getElementById("dl_button_', project_id[row_index],'").click()'), shiny::tags$button('DL', id = paste0('button_', project_id[row_index]), style = 'color: #FFFFFF; background-color: #393A35; font-size:10px; cursor:default; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'Shiny.onInputChange(\"download_click\",  this.id);')) %>% paste0 %>% HTML()
      
    }
    
  } else {
    
    button_report <- ''
    
    button_download <- ''
    
  }
  
  if(current_window_size >= 1800){
    
    paste0(
      
      shiny::tags$button('Summary', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #605CA8; font-size:10px; cursor:default; border-color: #669999; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"summary_click\",  this.id);') %>% paste0, 
      
      button_manage,
      
      button_report,
      
      button_download,
      
      button_pepview
      
    ) %>% return
    
  } else {
    
    paste0(
      
      shiny::tags$button('SUM', id = paste0('button_', row_index), style = 'color: #FFFFFF; background-color: #605CA8; font-size:10px; cursor:default; border-color: #669999; border-radius: 0.2em; margin-right: 0.5em; margin-bottom: 0.3em;  margin-top: 0.3em;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"summary_click\",  this.id);') %>% paste0, 
      
      button_manage,
      
      button_report,
      
      button_download,
      
      button_pepview
      
    ) %>% return
    
  }
  
}

checkpswd <- function(username, password_hash) {
  usrdb <- dbConnect(RSQLite::SQLite(), "../db/usrdb.db")	
  db_password_hash <- as.character(dbGetQuery(usrdb, paste("SELECT password FROM User WHERE username='",username,"'",sep="")))
  check <- identical(password_hash, db_password_hash)	
  dbDisconnect(usrdb)
  return(check)	
}

updatesessionid <- function(username, sessionid) {
  usrdb <- dbConnect(RSQLite::SQLite(), "../db/usrdb.db", flags = SQLITE_RW)	
  dbSendStatement(usrdb, paste("UPDATE User SET sessionid = '",sessionid,"' WHERE username='",username,"'",sep=""))
  dbSendStatement(usrdb, paste("UPDATE User SET logintime = '",as.character(Sys.time()),"' WHERE username='",username,"'",sep=""))
  dbDisconnect(usrdb)
}

run_time <- function(){
  
  dr <- data.frame(
    
    'maintain' = 0,
    
    'loop' = 1,
    
    'queue' = 10,
    
    'finished_project' = 20,
    
    'project_count' = 30,
    
    'user_count' = 5,
    
    'species_count' = 20,
    
    'version' = 1.01
    
  )
  
  write.csv(dr, file = paste0(path_prefix, '/.runtime'), row.names = FALSE)
  
}


get_msg <- function(pri_level, username){
  
  # msg-database subject, content, sender, type, rcvr, rdtag, rcpt_group,  valid_from, valid_thru
  message_db <- dbConnect(RSQLite::SQLite(), paste0(path_prefix, "/db/msg.db"))
  
  msg_list <- list()
  
  #提取anc typ 1 & rcpt_group >= pri_level
  msg_list[['anc']] <- dbGetQuery(message_db, paste0("SELECT * FROM Message WHERE valid_thru >= '", Sys.time() %>% as.numeric,"' and typ = 1 and rcpt_group >= ", pri_level)) %>% as.data.frame %>% arrange(., desc(valid_from))
  
  #提取notif typ 2 & rcvr = username & rdtag = 0
  msg_list[['notif']] <- dbGetQuery(message_db, paste0("SELECT * FROM Message WHERE valid_thru >= '", Sys.time() %>% as.numeric,"' and typ = 2 and rcvr = '", username, "' and rdtag = 0")) %>% as.data.frame %>% arrange(., desc(valid_from))
  
  #notif标记已读
  
  #dbSendQuery(message_db, paste0("UPDATE Message SET rdtag = '1' WHERE valid_thru >= '", Sys.time() %>% as.numeric,"' and typ = 2 and rcvr = '", username, "' and rdtag = 0"))
  
  #提取msg typ 3 & rcvr = username
  msg_list[['msg']] <- dbGetQuery(message_db, paste0("SELECT * FROM Message WHERE valid_thru >= '", Sys.time() %>% as.numeric,"' and typ = 3 and rcvr = '", username, "'")) %>% as.data.frame %>% arrange(., desc(valid_from))
  
  dbDisconnect(message_db)
  
  return(msg_list)
  
}
