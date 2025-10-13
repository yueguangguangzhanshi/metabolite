metr_pkgs<-c('lubridate','shinyauthr','rJava','PKI','digest','shinydashboard','shinycssloaders','RSQLite','jsonlite','DT','V8','shinyBS','shiny','shinyjs','shinyalert','bslib','shinysky','shinythemes','shinydashboardPlus','shinyWidgets','readxl','readr','tidyverse','mixOmics','gplots','pheatmap','plotly','dragulaR','factoextra','FactoMineR','ggsci','seqinr','ggrepel','ggplot2','openxlsx','tools','echarts4r')

for(i in 1:length(metr_pkgs)){

library(metr_pkgs[i], character.only = TRUE)

}

path_prefix1<-paste0(Sys.getenv("HOME"),"/shiny-server/BioAnalysis")

path_prefix2 <- paste0(Sys.getenv("HOME"),"/shiny-server/BioAnalysis")

# connect to, or setup and connect to local SQLite db
if (file.exists("db/sessionids_db")) {
  db <- dbConnect(SQLite(), "db/sessionids_db")
} else {
  db <- dbConnect(SQLite(), "db/sessionids_db")
  dbCreateTable(db, "sessionids", c(user = "TEXT", sessionid = "TEXT", login_time = "TEXT"))
}

# a user who has not visited the app for this many days
# will be asked to login with user name and password again
cookie_expiry <- 7 # Days until session expires

# This function must accept two parameters: user and sessionid. It will be called whenever the user
# successfully logs in with a password.  This function saves to your database.

add_sessionid_to_db <- function(user, sessionid, conn = db) {
  tibble(user = user, sessionid = sessionid, login_time = as.character(now())) %>%
    dbWriteTable(conn, "sessionids", ., append = TRUE)
}

# This function must return a data.frame with columns user and sessionid  Other columns are also okay
# and will be made available to the app after log in as columns in credentials()$user_auth

get_sessionids_from_db <- function(conn = db, expiry = cookie_expiry) {
  dbReadTable(conn, "sessionids") %>%
    dplyr::mutate(login_time = lubridate::ymd_hms(login_time)) %>%
    as_tibble() %>%
    dplyr::filter(login_time > now() - days(expiry))
}

user_base <- readRDS(paste0(path_prefix2,"/db/user_base.rds"))

ui<-navbarPage(title = "Metabonomics Analysis",id = "inTabset",
               
               header = list(tags$link(rel = "stylesheet", 
                                       type = "text/css", 
                                       href = paste0(Sys.getenv("HOME"),"/shiny-server/BioAnalysis/www/css/custom.css")),
                             useShinyjs(),
                             tags$script(HTML('
$(document).on("show.bs.modal", "#analysis_quant", function() {
  $(this).data("bs.modal").options.backdrop = "static";
  $(this).data("bs.modal").options.keyboard = false;
});
')),
                             tags$style(HTML("
    .dt-autofill-select, .dt-autofill-list, .DTED_Lightbox_Wrapper, .dteditor-modal {
      z-index: 999999 !important;
    }
  ")),
                             tags$script(HTML("
    $(document).on('DOMNodeInserted', '.dt-autofill-select', function(e){
        setTimeout(function(){
            var autofill = $('.dt-autofill-select');
            if(autofill.length>0){
                $('body').append(autofill);
            }
        }, 10);
    });
  ")),
                             tags$style(HTML("
      .modal-xl .modal-dialog {
        width: 90% !important;
        max-width: 1200px;
      }
      .modal-xl .modal-body {
        max-height: 80vh;
        overflow-y: auto;
      }
    ")),
                             useShinydashboard(),
                             useShinyalert(force = TRUE),
                             div(shinybusy::add_busy_spinner(
                               spin = "fulfilling-bouncing-circle",
                               color = "#112446",
                               timeout = 100,
                               position = c("bottom-right"),
                               onstart = TRUE,
                               margins = c(10, 10),
                               height = "50px",
                               width = "50px"
                             ),
                             div(class = "d-flex w-100 justify-content-between",
                                 div(class = "ml-auto",
                                     textOutput("username",inline = T), 
                                     shinyauthr::logoutUI(id = "logout",style = "color: green;")
                                     )
                                 )
                             )),
               
               selected = "Metabonomics analysis",
               
               theme = shinythemes::shinytheme("journal"),
               
               tabPanel(title = "Home",icon = icon(name="home"),
                        
                        div(h1(span("Omics",style="color: #2e2e66;margin-right:-10px;"),
                               span("Bioinformation",style="color: #0a0a66;margin-right:-10px;"),
                               span("Analysis",style="color: #2e2e66;margin-right:-10px;"),style="font-family: Georgia; font-size: 60px;"),
                            h4('Omics Bioinformation Analysis Platform' ,style="font-family: Georgia; color:#F5F5F5;"),
                            br(),
                            br(),
                            div(style = "margin: 0 47% 0 47%;",
                                pickerInput(
                                  inputId = "gap_start",
                                  label = NULL, 
                                  choices = c("Data Pre-processing"),
                                  options = list(
                                    style = "btn-primary",
                                    `selected-text-format`= "journal",
                                    title = "Get Start"),
                                  width = "fit"
                                )),
                            style = "text-align: center; background-image: url('texturebg.png');padding: 60px"),
                        
                        br(),
                        br(),
                        
                        div(style = "margin-left:15px;",
                            h3("Introduction",style = "font-weight:bold;"),
                            h4("Omics Bioinformation Analysis Platform."),
                            br(),
                            
                            h3("Function",style = "font-weight:bold;"),
                            br(),
                            br(),
                            fluidRow(
                              
                              column(width = 3,div(span(class="fas fa-file-medical-alt fa-8x",style="color: #4682B4"),
                                                   h3("Spectrum Identification"),
                                                   p("Omics Bioinformation Analysis Platform"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-chart-pie fa-8x",style="color: #4682B4"),
                                                   h3("Data Pre-processing"),
                                                   p("Omics Bioinformation Analysis Platform"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-chart-line fa-8x",style="color: #4682B4"),
                                                   h3("Biomarker"),
                                                   p("Omics Bioinformation Analysis Platform"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-dna fa-8x",style="color: #4682B4"),
                                                   h3("Functional Analysis"),
                                                   p("Proteomics analysis with GO and KEGG"),
                                                   style = "text-align: center; background-repeat:no-repeat;"))
                              
                            ),
                            
                            br(),
                            br(),
                            
                            fluidRow(
                              
                              column(width = 3,div(span(class="fas fa-microscope fa-8x",style="color: #4682B4"),
                                                   h3("Clinical Analysis"),
                                                   p("Disease and prognosis prediction"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-database fa-8x",style="color: #4682B4"),
                                                   h3("Databases"),
                                                   p("Omics Bioinformation Analysis Platform"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-school fa-8x",style="color: #4682B4"),
                                                   h3("Relative Websites"),
                                                   p("Omics Bioinformation Analysis Platform"),
                                                   style = "text-align: center; background-repeat:no-repeat;")),
                              
                              column(width = 3,div(span(class="fas fa-people-carry fa-8x",style="color: #4682B4"),
                                                   h3("Question & Answers"),
                                                   p("The Q&As about the Omics Bioinformation Analysis operation guide"),
                                                   style = "text-align: center; background-repeat:no-repeat;"))
                              
                            ),
                            br(),
                            
                            h3("Reference",style = "font-weight:bold;"),
                            h4("[1] Omics Bioinformation Analysis Platform"),
                            h4("[2] Omics Bioinformation Analysis Platform"),
                            br(),
                            br()
                        )
                        
               ),
               
               navbarMenu(title = "Analysis",icon = icon(name="chart-bar"),
                          
                          tabPanel(title = "Metabonomics analysis", 
                                   shinyauthr::loginUI(id = "login", cookie_expiry = cookie_expiry),
                                   uiOutput("metabo_tag"),
                                   
                          ),
                          
                          tabPanel(title = "Data Pre-processing",
                                   
                                   # ident panel head ----
                                   div(style = "text-align: center; background-image: url('texturebg.png');padding: 60px; s", display='block',
                                       h1("Data Pre-processing",style="font-family: Georgia; color: #FFFFFF"),
                                       h5("Data Pre-processing",style="font-family: Georgia; color:#F5F5F5;"),
                                       br(),
                                       br(),
                                       actionButton("confirm_ident",label="Get Start",styleclass = "primary",size = 'large')),
                                   br(),
                                   br(),
                                   
                                   # ident analysis div ----
                                   fluidRow(
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border: 5px solid steelblue; height: 300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_stat_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="50%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing"))),
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_motif_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="100%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing"))),
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_gta_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="50%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing"))),
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_tg_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="100%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing")))
                                     
                                   ),
                                   br(),
                                   br(),
                                   
                                   fluidRow(
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border: 5px solid steelblue; height: 300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_gpgs_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="50%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing"))),
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_gpgp_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="50%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing"))),
                                     
                                     column(width = 3,
                                            div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                                                div(style="text-align: left; margin: 5% 0 0 5%;",
                                                    prettyCheckbox(
                                                      inputId = "ident_uniprot_check",
                                                      label = "", 
                                                      value = TRUE,
                                                      icon = icon("check"),
                                                      status = "primary",
                                                      bigger = T,
                                                      animation = "jelly"
                                                    )),
                                                img(src="pre_pca_check.png",height="50%",width="35%"),
                                                h4("Data Pre-processing"),
                                                p("Data Pre-processing")))
                                     
                                   ),
                                   br(),
                                   br(),
                                   bsModal("analysis_ident", "Data Pre-processing", 
                                           "confirm_ident",
                                           size="large",
                                           uiOutput("ident_panel"))
                                   
                          )
                          
               ),
               
               tabPanel(title = "Demos",uiOutput("report_o"),value = "demo"),
               
               tabPanelBody(type = "hidden",value="Appendix",uiOutput("appendix_o"))
               
)

server <- function(input, output, session) {
  
  # call the logout module with reactive trigger to hide/show
  logout_init <- shinyauthr::logoutServer(
    id = "logout",
    active = reactive(credentials()$user_auth)
  )
  
  # call login module supplying data frame, user and password cols
  # and reactive trigger
  credentials <- shinyauthr::loginServer(
    id = "login",
    data = user_base,
    user_col = user,
    pwd_col = password,
    cookie_logins = TRUE,
    sessionid_col = sessionid,
    cookie_getter = get_sessionids_from_db,
    cookie_setter = add_sessionid_to_db,
    sodium_hashed = TRUE,
    reload_on_logout = TRUE,
    log_out = reactive(logout_init())
  )
  
  # pulls out the user information returned from login module
  user_info <- reactive({
    credentials()$info
  })
  
  output$username <- renderText({
    
    paste("Logged in as:", user_info()$name)
    
  })
  
  sizelimit<-5000
  
  options(shiny.maxRequestSize = sizelimit * 1024 ^ 2)
  
  useShinyalert(force = TRUE)
  
  output$metabo_tag <- renderUI({
    
    req(credentials()$user_auth)
    
    tagList(
      
      # quant panel head ----
      div(style = "text-align: center; background-image: url('texturebg.png');padding: 60px; s", display='block',
          h1("Metabonomics analysis",style="font-family: Georgia; color: #FFFFFF"),
          h5("Focus on site specific quantitative analysis",style="font-family: Georgia; color:#F5F5F5;"),
          br(),
          br(),
          actionButton("confirm_quant",label="Get Start",styleclass = "primary",size = 'large')),
      br(),
      br(),
      
      # quant analysis div ----
      fluidRow(
        
        column(width = 3,
               div(style = "margin: 0 0 0 3%; text-align: center; border: 5px solid steelblue; height: 300px;",
                   div(style="text-align: left; margin: 5% 0 0 5%;",
                       prettyCheckbox(
                         inputId = "quant_pca_check",
                         label = "", 
                         value = TRUE,
                         icon = icon("check"),
                         status = "primary",
                         bigger = T,
                         animation = "jelly"
                       )),
                   img(src="quant_pca_check.png",height="50%",width="50%"),
                   h4("PCA"),
                   p("Principal components analysis"))),
        
        column(width = 3,
               div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                   div(style="text-align: left; margin: 5% 0 0 5%;",
                       prettyCheckbox(
                         inputId = "quant_oplsda_check",
                         label = "", 
                         value = TRUE,
                         icon = icon("check"),
                         status = "primary",
                         bigger = T,
                         animation = "jelly"
                       )),
                   img(src="quant_oplsda_check.png",height="50%",width="50%"),
                   h4("OPLS-DA"),
                   p("Orthogonal partial least-squares discrimination analysis"))),
        
        column(width = 3,
               div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                   div(style="text-align: left; margin: 5% 0 0 5%;",
                       prettyCheckbox(
                         inputId = "quant_heatmap_check",
                         label = "", 
                         value = TRUE,
                         icon = icon("check"),
                         status = "primary",
                         bigger = T,
                         animation = "jelly"
                       )),
                   img(src="quant_heatmap_check.png",height="50%",width="50%"),
                   h4("Heatmap"),
                   p("Cluster analysis of different groups"))),
        
        column(width = 3,
               div(style = "margin: 0 0 0 3%; text-align: center; border:5px solid steelblue; height:300px;",
                   div(style="text-align: left; margin: 5% 0 0 5%;",
                       prettyCheckbox(
                         inputId = "quant_volcano_check",
                         label = "", 
                         value = TRUE,
                         icon = icon("check"),
                         status = "primary",
                         bigger = T,
                         animation = "jelly"
                       )),
                   img(src="quant_volcano_check.png",height="50%",width="50%"),
                   h4("Variance analysis"),
                   p("Differential glycopeptides between two groups of samples (Please make sure there are at least 2 replicates in each group.)")))
        
      ),
      br(),
      br(),
      
      column(h3("Release Version:"),
             verbatimTextOutput("ReleaseVersion"),
             offset=2,width = 8),
      
      br(),
      br()
      
    )
  })
  
  dep_content <- 'quant_lib.r'
  
  for(i in dep_content){
    
    if(file.exists(paste0(path_prefix2,'/etc/dep/', i))) {
      
      source(paste0(path_prefix2,'/etc/dep/', i), local = TRUE)
      
    }
    
  }
  
  observeEvent(input$confirm_quant, {
    showModal(modalDialog(
      title = "Metabonomics analysis",
      size = "l",  # 必须设这个才能激活 modal-lg 样式
      easyClose = TRUE,
      footer = modalButton("Close"),
      div(class = "modal-xl",  # 👈 用这个 class 包住
          uiOutput("quant_panel")
      )
    ))
  })
  
  observeEvent(input$gap_start,{
    
    updateTabsetPanel(
      session = getDefaultReactiveDomain(),
      inputId= "inTabset",
      selected = input$gap_start
    )
    
    updatePickerInput(
      session,
      "gap_start",
      choices = c("Data Pre-processing"),
      options = list(
        style = "btn-primary",
        `selected-text-format`= "static",
        title = "Get Start")
    )
    
  })
  
  observe({
    
    updateTabsetPanel(session, 'inTabset', names(parseQueryString(session$clientData$url_search))[1]%>%str_replace("\\?",""))
    
  })
  
  observe({
    
    if (!is.null(names(parseQueryString(session$clientData$url_search)))) {
      
      if (names(parseQueryString(session$clientData$url_search))[1]%>%str_replace("\\?","")=="demo") {
        
        report_id<<-parseQueryString(session$clientData$url_search)%>%.$id2
        
        if (file.exists(paste0(path_prefix2, '/data/usrdata/', report_id, '/report/report_server.html'))) {
          
          system(paste0('cp -rf ',path_prefix2,'/data/usrdata/', report_id, '/report/report_server.html ',path_prefix1,'/www/report_tmp.html'))
          
          output$report_o<-renderUI(
            
            fluidPage(
              
              fluidRow(downloadButton("c_plot_download", "Result Download")),
              
              fluidRow(
                
                tags$iframe(
                  
                  id = 'reportframe',
                  
                  seamless = 'seamless',
                  
                  width = '100%',
                  
                  height= '1000px',
                  
                  scrolling = 'yes',
                  
                  frameborder = '0',
                  
                  framespacing = '0',
                  
                  src = 'report_tmp.html'))
              
            )
            
          )
          
        } else {
          
          shinyalert::shinyalert("ID Wrong !", "Please confirm report url is true.", type = "error")
          
        }
        
      }
      
    }
    
  })
  
  output$c_plot_download<-downloadHandler(
    
    filename = function() {
      
      # input1 <- readRDS(paste0(path_prefix1,'/data/usrdata/',report_id,'/data/input.rds'))
      
      paste0("BioAnalysis_",format(Sys.time(),"%Y%m%d%H%M"),".zip")
      
    },
    
    content = function(file) {
      
      file.copy(paste0(path_prefix2,'/data/usrdata/',report_id,'/report/report.zip'), file)
      
    }
    
  )
  
  output$ReleaseVersion<-renderPrint(
    
    cat(
      "2025.10.10 version 1.0.5 Add the function of filling in the name based on the annotation ID.",
      "2025.09.22 version 1.0.4 Fix kegg enrichment doesn't produce all pathways bug.",
      "2025.09.15 version 1.0.3 Add NormAE QC normalization function.",
      "2025.09.08 version 1.0.2 Add RLSC QC normalization function.",
      "2025.09.01 version 1.0.1 Add LMSD and HMDB database class function.",
      "2025.08.15 version 1.0.0 Metabolite analysis website launch",
      sep="\n")
    
  )
  
}

shinyApp(ui = ui, server = server)