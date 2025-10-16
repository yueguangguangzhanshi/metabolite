
observeEvent(input$start_clinic, {
  
  myhide(c('my_widget'))
  
  runjs("$('#my_maintab').animate({width: '100%'}, 300)")
  
  shinydashboard::updateTabItems(session = session, inputId = 'tabs', selected = 'clinic_tab')
  
})

# 主界面 ----

output$clinic_panel <- renderUI({
  
  box(title = 'Clinical Analysis', solidHeader = TRUE, status = "primary", width = 12, style = 'margin: 0 auto; min-height: 500px;',
      
      uiOutput("clinic_nav"),
      
      uiOutput("clinic_project_selector"),
      
      uiOutput("clinic_project_confirm"),
      
      shinyjs::useShinyjs(),
      
      useShinyalert()
      
  )
  
})

sizelimit<-30

output$clinic_project_selector <- renderUI({
  
  fluidRow(class = 'para_panel', style = 'margin-right: auto; margin-left: auto;',
           
           column(3, '' ),
           
           column(6,
                  
                  div(style = 'width: 400px; margin-right: auto; margin-left: auto;',
                      
                      box(title = "Clinical Analysis", solidHeader = TRUE, status = "primary", width = 12, style = '',

                          div(id = 'inst_container', style = 'margin-top: 1.5em; text-align:center; max-height: 2em; margin-bottom: -1em;', 
                              
                              div(style = 'display: inline-block; vertical-align: top; padding-top: 0.5em;',
                                  
                                  shiny::span(icon("info"), href = '', style = 'margin-left:0.5em; color:#393A35; text-decoration: none;font-size: 14px; width: 0.5em; margin-right: 0.8em;', title = ''),
                                  
                                  paste0('Click to download test ', shiny::downloadLink(outputId = 'download_clinic_testdata', label = 'Demo Data Set') %>% paste0 %>% HTML(), '.') %>% shiny::HTML()
                                  
                              )
                              
                          ) %>% paste0 %>% shiny::HTML(),
                          
                          runjs("Shiny.onInputChange('inst_msg', new Date().getTime().toLocaleString())"),
                          
                          div(style = ' max-width: 78%;  margin-top: 3em;  margin-right: auto; margin-left: auto; ',
                              
                              fileInput("glycoClinic", paste0("Upload the list file of quantitative results of glycopeptide library search and control it within ",sizelimit," MB:"),
                                        multiple = FALSE,
                                        accept = c(".list")),
                              
                              fileInput("ClinicGroup", paste0("Upload glycopeptide group file:"),
                                        multiple = FALSE,
                                        accept = c(".txt")),
                              
                              uiOutput("ClinicSurvival_o")

                          ),
                          
                          div(style = 'max-width: 70%; margin-bottom: 3em; margin-right: auto; margin-left: auto;',
                              
                              shinysky::actionButton(inputId = 'clinic_apply_parameter', label = 'Next', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin: 2em 0;')
                              
                          )
                          
                      )
                      
                  )
                  
           ), 
           
           column(3, '' )
           
  )
  
})

observe({
  
  if (!is.null(input$clinic_survival_check)) {
    
    if (input$clinic_survival_check) {
      
      output$ClinicSurvival_o<-renderUI(
        
        fileInput("ClinicSurvival", paste0("Upload survival information file (see demo for details):"),
                  multiple = FALSE,
                  accept = c(".xlsx"))
        
      )
      
    }else{
      
      output$ClinicSurvival_o<-NULL
      
    }
    
  }
  
})

# download ----

output$download_clinic_testdata <- downloadHandler(
  
  filename = function() {
    
    'GlycoClinic.zip'
    
  },
  
  content = function(file) {
    
    file.copy(paste0('/home/biognosis/shiny-server/GAP/data/testdata/GlycoClinic.zip'), file)
    
  }
  
)

# 提交数据 ----

observeEvent(input$clinic_apply_parameter, {
  
  myhide(c('clinic_project_selector'))
  
  myshow(c('clinic_project_confirm'))
  
  if(input[["glycoClinic"]][['datapath']] %>% is.null){
    
    shinyalert::shinyalert("Missing files !", "Please confirm data has been uploaded successfully.", type = "error")
    
    myhide(c('clinic_project_confirm'))
    
    myshow(c('clinic_project_selector'))
    
  }
  
  shinyjs::removeCssClass(id = "nav-clinic-select-project", class = "active")
  
  shinyjs::addCssClass(id = "nav-clinic-apply", class = "active")
  
  project_id <- generate_jobid()
  
  path_prefix <<- paste0('/home/biognosis/shiny-server/GAP')
  
  data_path <<- paste0('/home/biognosis/shiny-server/GAP/data/usrdata/', project_id, '/data')
  
  para <<- list();para[['info']]<<-list(para[['info']])
  
  para[['info']]['type'] <<- 'clinic'
  
  para[['info']]['data_path'] <<- data_path
  
  para[['info']]['project_id'] <<- project_id
  
  output$clinic_project_confirm <- renderUI({
    
    fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
             
             column(3, '' ),
             
             column(6,
                    
                    div(style = 'max-width: 600px; margin-right: auto; margin-left: auto;',
                        
                        box(title = "Report Info", solidHeader = TRUE, status = "primary", width = 12, style = '',
                            
                            div(style = 'max-width: 85%; margin: 3em auto;',
                                
                                p('Project Name:', style = 'font-size: 14px; font-weight: bold;'),
                                
                                div(style = '',textInput('project_name', label = NULL, value = paste0('Project_',format(Sys.time(),"%Y%m%d%H%M")))),
                                
                                p('Organization Information (And fill in contact information):', style = 'font-size: 14px; font-weight: bold;'),
                                
                                div(style = '',textInput('organ_info', label = NULL)),
                                
                                p('E-mail Address:', style = 'font-size: 14px; font-weight: bold;'),
                                
                                div(style = '',textInput('email', label = NULL, value = ''))
                                
                            ),
                            
                            div(style = 'width: 100%; text-align:center;',
                                
                                shinysky::actionButton(inputId = 'confirm_clinic1', label = 'Confirm', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-bottom: 0; margin-top: 3em;')
                                
                            ),
                            
                            div(style = 'width: 100%; text-align:center;',
                                
                                shinysky::actionButton(inputId = 'return_clinic_data_upload', label = 'Return', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-top: 0.5em; margin-bottom: 30px;')
                                
                            )
                            
                        )
                        
                    )
                    
             ), 
             
             column(3, '' )
             
    )
    
  })
  
})

# confirm ----

observeEvent(input$confirm_clinic1, {
  
  if (trimws(input$organ_info)==""|trimws(input$email)=="") {
    
    showModal(modalDialog(
      title = "Important message",
      paste0("Please fill in blank!"),
      easyClose = T
    ))
    
  }else{
    
    para[['info']][['project_name']] <<- input$project_name
    
    para[['info']][['user']] <<- input$organ_info
    
    para[["info"]][["email"]] <<- input$email
    
    process_item<-c()
    
    for (i in c("clinic_roc_check","clinic_survival_check")) {
      
      if (input[[i]]) {
        
        process_item<-c(process_item,i)
        
      }
      
    }
    
    para[['info']][['process_item']] <<- process_item
    
    dir.create(path = para[['info']][['data_path']], recursive = TRUE)
    
    file.copy(from = input[["glycoClinic"]][['datapath']], to = paste0(para[['info']][['data_path']], '/clinic.txt'))
    
    file.copy(from = input[["ClinicGroup"]][['datapath']], to = paste0(para[['info']][['data_path']], '/clinicGroup.txt'))
    
    if (input$clinic_survival_check) {
      
      file.copy(from = input[["ClinicSurvival"]][['datapath']], to = paste0(para[['info']][['data_path']], '/survival.xlsx'))
      
    }
    
    save(list = c('para'), file = paste0(para[['info']][['data_path']], '/para.Rdata'), envir = .GlobalEnv)
    
    system(paste0("cp -rf /home/biognosis/shiny-server/GAP/db/projects_list.db /home/biognosis/shiny-server/GAP/db1/projects_list.db"))
    
    projects_list <- dbConnect(RSQLite::SQLite(), "/home/biognosis/shiny-server/GAP/db1/projects_list.db", flags = SQLITE_RWC)
    
    dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", para[['info']][['project_id']], "','", para[['info']][['project_name']], "', '", input$organ_info, "', '1', '", as.character(Sys.time()),"', 'clinic', '",para[["info"]]["email"],"')"))
    
    dbDisconnect(projects_list)
    
    system(paste0("cp -rf /home/biognosis/shiny-server/GAP/db1/projects_list.db /home/biognosis/shiny-server/GAP/db/projects_list.db"))
    
    # paste0("'/usr/lib/R/bin/Rscript' '", path_prefix, "/etc/dep/bg_script.r' ",para[['info']]['project_id']) %>% system(wait = F)
    
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
    
  }
  
})

# 从结果确认返回选择 ----

observeEvent(input$return_clinic_data_upload, {
  
  myshow(c('clinic_species_selector'))
  
  myhide(c('clinic_project_confirm'))
  
  shinyjs::removeCssClass(id = "nav-clinic-apply", class = "active")
  
  shinyjs::addCssClass(id = "nav-clinic-select-project", class = "active")
  
})
