#function

observeEvent(input$start_function, {

	myhide(c('my_widget'))

	runjs("$('#my_maintab').animate({width: '100%'}, 300);")
	
	shinydashboard::updateTabItems(session = session, inputId = 'tabs', selected = 'function_tab')
	
})

output$fun_species_selector <- renderUI({
				
	speciesdb <- dbConnect(RSQLite::SQLite(), "/home/biognosis/shiny-server/GAP/db/speciesdb.db")
	
	species_list <<- dbGetQuery(speciesdb, paste0("SELECT * FROM Species")) %>% as.data.frame
	
	dbDisconnect(speciesdb)

	choices <- list()
	
	species_L1 <- species_list$species_L1 %>% unique

	for(i in 1:length(species_L1)){
	
		species_list_sub <- species_list[which(species_list$species_L1 == species_L1[i]),]
		
		species_names <- species_list_sub$sci_names
	
		species_ch_names <- species_list_sub$ch_names
	
		choices_string <- paste0(species_names, ' (', species_ch_names, ')')
	
		choices[[species_L1[i]]] <- choices_string %>% as.list

	
	}
	
	fluidRow(class = 'para_panel', style = 'margin-right: auto; margin-left: auto;',

		column(3, '' ),
		
		column(6,
		
			div(style = 'width: 400px; margin-right: auto; margin-left: auto;',
			
				box(title = "Select species", solidHeader = TRUE, status = "primary", width = 12, style = '',

					div(style = ' max-width: 78%;  margin-top: 3em;  margin-right: auto; margin-left: auto; ',
					
						div(style = 'display: inline-block; width: 94%; padding-left:1em;',
						
							selectizeInput(inputId = 'fun_selected_species', label = NULL, choices = choices, selected = NULL, multiple = FALSE, options = list(openOnFocus = FALSE))
							
						),
						
						div(style = 'display: inline-block; vertical-align: top; padding-top: 0.5em;',
						
							shiny::a(icon("info"), href = '#', style = 'margin-left:0.5em; color:#999999; text-decoration: none;font-size: 14px; width: 0.5em;', onmousedown = 'event.preventDefault();event.stopPropagation();Shiny.onInputChange(\"spdetail_click\", new Date().toLocaleString());', title = 'Click to view details')
							
						),
						
					),

					div(style = 'max-width: 70%; margin-bottom: 3em; margin-right: auto; margin-left: auto;',
					
						shiny::downloadLink(outputId = 'fun_downloadFasta', label = shinysky::actionButton(inputId = 'dl_btn_fasta', label = 'FASTA Download', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 1em; margin-top: 2em;') %>% paste0 %>% HTML()),

						shinysky::actionButton(inputId = 'fun_apply_species', label = 'Next', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 0;')
					
					)

				)
				
			)
		
		), 
		
		column(3, '' )
	
	)
  
})

observeEvent(input$fun_selected_species, {
  
  species_sci_name <<- input$fun_selected_species %>% str_replace_all('( \\().*', '')
  
  species_index <<- which(species_list$sci_names == species_sci_name)
  
  species_db_name <<- species_list$name[species_index]
  
  species_ch_names <<- species_list$ch_names[species_index]
  
  species_kegg <<- species_list$kegg[species_index]
  
  species_taxid <<- species_list$taxid[species_index]
  
  species_updatetime <<- species_list$updatetime[species_index]
  
})

#download
output$fun_downloadFasta <- downloadHandler(

	filename = function() {

		paste0(species_db_name, '.fasta', sep='')

	},

	content = function(file) {

		file.copy(paste0('/home/biognosis/shiny-server/GAP/db/fasta/', species_db_name, '.fasta'), file)

	}
  
)

#从匹配返回输入
observeEvent(input$return_fun_input, {
	
	myhide(c('fun_mapping_output'))
	
	myshow(c('fun_para_input'))
	
	shinyjs::removeCssClass(id = "nav-fun-match-data", class = "active")

	shinyjs::addCssClass(id = "nav-fun-upload-data", class = "active")
	
})

#从结果确认返回匹配
observeEvent(input$return_mapping_input, {
	
	myshow(c('fun_mapping_output'))
	
	myhide(c('fun_project_confirm'))
	
	shinyjs::removeCssClass(id = "nav-fun-apply", class = "active")

	shinyjs::addCssClass(id = "nav-fun-match-data", class = "active")
	
})

#从输入返回物种选择
observeEvent(input$return_fun_species_selector, {
	
	myhide(c('fun_para_input'))
	
	myshow(c('fun_species_selector'))
	
	shinyjs::removeCssClass(id = "nav-fun-upload-data", class = "active")

	shinyjs::addCssClass(id = "nav-fun-select-species", class = "active")
	
})

#主界面 ----
output$fun_panel <- renderUI({

	box(title = 'FUNCTION', solidHeader = TRUE, status = "primary", width = 12, style = 'margin: 0 auto; min-height: 500px;',

		uiOutput("fun_nav"),
		
		uiOutput("fun_species_selector"),
	
		uiOutput("fun_para_input"),
		
		uiOutput("fun_mapping_output"),
		
		uiOutput("fun_project_confirm"),
		
		shinyjs::useShinyjs(),
		
		useShinyalert()

	)
	
})

observeEvent(input$fun_apply_species, {

	myhide('fun_species_selector')
	
	myshow('fun_para_input')
	
	shinyjs::removeCssClass(id = "nav-fun-select-species", class = "active")

	shinyjs::addCssClass(id = "nav-fun-upload-data", class = "active")

	output$fun_para_input <- renderUI({

		fluidRow(class = 'para_panel', style='margin-right:auto; margin-left:auto;',
		
			column(3, '' ),
			
			column(6,
			
				div(style = 'max-width: 400px; margin-right: auto; margin-left: auto;',
				
					box(title = "Upload File", solidHeader = TRUE, status = "primary", width = 12, style = '',
					
						div(style = 'max-width: 70%; margin: 3em auto;',
						
							h5('Paste ID or upload file:'),
						
							textAreaInput(inputId = 'protein_list', label = '', height = 80),

							div(style = 'margin-top:1em;',
								shiny::a(href = '#', id = 'sample_1', style = 'color:#555555; text-decoration: none;', onmousedown = 'event.preventDefault(); event.stopPropagation(); Shiny.onInputChange(\"fun_sample\", new Date().toLocaleString());', '#SAMPLES SET')
							),

							fileInput(inputId = 'fun_file_upload',  label = '', accept = c('.csv', '.txt', '.tsv', '.xls')),

							shinysky::actionButton(inputId = 'fun_mapping_input', label = 'Match database', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 0;'),
							
							shinysky::actionButton(inputId = 'return_fun_species_selector', label = 'Back', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-top: 1em;')
						
						)

					)
				
				)
				
			), 
		
			column(3, '' )
		
		)

	})
	
})

observeEvent(input$fun_sample, {

	sample_protein_list <- generate_protein_sample(species_db_name)
	
	updateTextAreaInput(session, inputId = 'protein_list', value = sample_protein_list)

})

observeEvent(input$fun_mapping_input, {

	if(is.null(input$fun_file_upload) & (input$protein_list == '')){
	
		shinyalert::shinyalert('Please enter ID or upload file', NULL, type = "error")
	
	} else {

		#updateTextInput(session, inputId = 'box_title', value = '数据库匹配')
		
		myhide(c('fun_para_input', 'fun_species_selector'))
		
		#shinyjs::addCssClass(id = 'fun_mapping_output', class = "show")
		
		myshow(c('fun_mapping_output'))
		
		shinyjs::removeCssClass(id = "nav-fun-upload-data", class = "active")

		shinyjs::addCssClass(id = "nav-fun-match-data", class = "active")

		match_result <- generate_fun_mapping_table(species_db_name, input$fun_file_upload, input$protein_list)
		
		mapping_table <- match_result[['annot_table']] %>% as.data.frame 

		if(nrow(mapping_table) > 0){
		
			mapping_table_out <<- cbind(Pick = lapply(1:nrow(mapping_table), FUN = generate_checkbox_fun) %>% unlist, mapping_table)

			output$fun_mapping_output <- renderUI({

				fluidRow(class = '', style = 'width: 90%; margin: 1em auto; padding: 0px;',
					
					DT::renderDataTable(
					
						mapping_table_out,
						
						escape = FALSE,	
						
						options = list(
						
							dom = 'lfrtp',
							
							scroller = TRUE, scrollX = T,
							
							rowCallback = JS("function(r,d) {$(r).attr('height', '30px;')}"),
							
							drawCallback= JS(
								'function(settings) {Shiny.bindAll(this.api().table().node());}'
							)

						),
						
						selection = "none"
					
					),


					shinysky::actionButton(inputId = 'apply_fun_input', label = 'Submit', styleclass = 'primary', style = 'padding: 6px 8px; width: 8em; float: right; margin-right: 1em; margin-top: 2em;'),
					
					shinysky::actionButton(inputId = 'return_fun_input', label = 'Back', styleclass = 'primary', style = 'padding: 6px 8px; width: 8em; float: right; margin-right: 1em; margin-top: 2em;')
				
				)

			})
			
		} else {
		
			shinyalert::shinyalert("No records were matched in the database.", type = "error", confirmButtonText = "Back", 
			
				callbackR = function(x) { if(x != FALSE)
				
					runjs(paste0(
						"(function() {
						location.href='",isolate(session$clientData$url_protocol),"//",isolate(session$clientData$url_hostname),":",isolate(session$clientData$url_port),"/GlycAP/' ;
						}());"
					))
					
				}
			
			)
		
		}

	}
	
})

observeEvent(input$apply_fun_input, {

	shinyjs::removeCssClass(id = "nav-fun-match-data", class = "active")

	shinyjs::addCssClass(id = "nav-fun-apply", class = "active")

	check_list <- list()

	for(i in 1:nrow(mapping_table_out)){

		check_list[[i]] <- ifelse(input[[paste0('checkbox_', i)]] %>% is.null, '', input[[paste0('checkbox_', i)]])
	
	}
	
	check_list <- check_list %>% unlist

	protein_list <- mapping_table_out[which(check_list == TRUE | check_list ==''), 'Accession']
	
	project_id <- generate_jobid()
	
	data_path <<- paste0('/home/biognosis/shiny-server/GAP/data/usrdata/', project_id, '/data')

	para <<- list();para[['info']]<<-list(para[['info']])
	
	para[['info']]['type'] <<- 'function'
	
	para[['info']]['data_path'] <<- data_path
	
	para[['info']]['project_id'] <<- project_id
	
	para[['info']]['db_sci_name'] <<- species_sci_name
	
	para[['info']]['db_name'] <<- species_db_name
	
	para[['info']]['ch_names'] <<- species_ch_names
	
	para[['info']]['kegg'] <<- species_kegg
	
	para[['info']]['taxid'] <<- species_taxid
	
	para[['info']]['updatetime'] <<- species_updatetime
	
	para[['data']][['protein_list']][['FUNC']] <<- protein_list

	output$fun_project_confirm <- renderUI({

		fluidRow(class = 'para_panel', style = 'margin-right:auto; margin-left:auto;',
		
			column(3, '' ),
			
			column(6,
			
				div(style = 'max-width: 400px; margin-right: auto; margin-left: auto;',
				
					box(title = "Project information", solidHeader = TRUE, status = "primary", width = 12, style = '',
					
						div(style = 'max-width: 70%; margin: 3em auto;',
						
						    p('Project Name:', style = 'font-size: 14px; font-weight: bold;'),
						    
						    div(style = '',textInput('project_name', label = NULL, value = paste0('Project_',format(Sys.time(),"%Y%m%d%H%M")))),
						    
						    p('Organization Information (And fill in contact information):', style = 'font-size: 14px; font-weight: bold;'),
						    
						    div(style = '',textInput('organ_info', label = NULL)),
						    
						    p('E-mail Address:', style = 'font-size: 14px; font-weight: bold;'),
						    
						    div(style = '',textInput('email', label = NULL, value = '')),

							div(style = 'display: inline-block; margin-top: 1.5em;',
							
								p('Species:', style = 'display: inline-block; font-size: 14px; font-weight: bold;'),
								
								p(input$fun_selected_species, style = 'display: inline-block;'),
							
								shiny::a(icon("info"), href = '#', style = 'margin-left:0.5em; color:#999999; text-decoration: none; display: inline-block; font-size: 14px; width: 0.5em;', onmousedown = 'event.preventDefault();event.stopPropagation();Shiny.onInputChange(\"spdetail_click\", new Date().toLocaleString());', title = 'Click to view details') %>% paste0 %>% HTML()
								
							),
							
							p(paste0('(Total selected ', length(protein_list), ' proteins)'), style = 'font-size: 12px; color: #555555; font-style: italic;'),
							
							p('Select chart style:', style = 'font-size: 14px; font-weight: bold; margin-top: 2em'),
							
							div(style = 'margin-top: 1em; padding-left: 0;',
							
								radioButtons(inputId = 'plot_theme', label = NULL,
								
									choices = list('default'),
									
									selected = list('default')
									
								)
								
							),
							
							shinysky::actionButton(inputId = 'confirm_fun1', label = 'Confirm', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 0; margin-top: 3em;'),
							
							shinysky::actionButton(inputId = 'return_mapping_input', label = 'Back', styleclass = 'primary', style = 'padding: 6px 8px; width: 100%; margin-bottom: 0; margin-top: 1em;')
						
						)

					)
				
				)
				
			), 
		
			column(3, '' )
		
		)

	})

	myhide('fun_mapping_output')
	
	myshow('fun_project_confirm')

})

observeEvent(input$confirm_fun1, {
  
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
    
    for (i in c("fun_go_check","fun_kegg_check","fun_string_check")) {
      
      if (input[[i]]) {
        
        process_item<-c(process_item,i)
        
      }
      
    }
    
    para[['info']][['process_item']] <<- process_item
    
    para[['info']]['plot_theme'] <<- input$plot_theme
    
    dir.create(path = para[['info']][['data_path']], recursive = TRUE)
    
    save(list = c('para'), file = paste0(para[['info']][['data_path']], '/para.Rdata'))
    
    path_prefix <<- paste0('/home/biognosis/shiny-server/GAP')
    
    # paste0("'/usr/lib/R/bin/Rscript' '", path_prefix, "/etc/dep/bg_script.r' ",para[['info']][['project_id']]) %>% system(wait = F)
    
    system("cp -rf /home/biognosis/mnt/os2/GlycAP/db/projects_list.db /home/biognosis/shiny-server/GAP/db/projects_list.db")
    
    projects_list <- dbConnect(RSQLite::SQLite(), "/home/biognosis/shiny-server/GAP/db/projects_list.db", flags = SQLITE_RWC)
    
    dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", para[['info']][['project_id']], "','", para[['info']][['project_name']], "', '", input$organ_info, "', '1', '", as.character(Sys.time()),"', 'function', '",para[["info"]]["email"],"')"))
    
    dbDisconnect(projects_list)
    
    system("cp -rf /home/biognosis/shiny-server/GAP/db/projects_list.db /home/biognosis/mnt/os2/GlycAP/db/projects_list.db")
    
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

#生成内容

generate_fun_mapping_table <- function(species_db_name, fun_file_upload, protein_list){

	if(fun_file_upload %>% is.null){
	
		input_acc <- protein_list %>% strsplit(split = ";\\s*|\\s+|,\\s*|-\\s*") %>% unlist
	
	} else {
	
		file_path <- fun_file_upload$datapath
		
		ext <- tools::file_ext(file_path)
	
		if(ext %in% c('csv', 'txt', 'tsv')) {
		
			input_acc <- read.csv(fun_file_upload$datapath, header = FALSE, stringsAsFactors = F) %>% unlist %>% unname

		}
		
		if(ext == 'xls') {
		
			input_acc <- read_xls(fun_file_upload$datapath, col_names = FALSE) %>% unlist %>% unname
		
		}
	
	}

	input_acc <- input_acc %>% accid_process
	
	db_acc <- paste0('/home/biognosis/shiny-server/GAP/db/acc_list/', species_db_name, '.txt') %>% read_csv(col_names = FALSE) %>% accid_process
	
	match_acc <- db_acc[which(db_acc %in% input_acc)]
	
	if(match_acc %>% length > 0){
	
		load(paste0('/home/biognosis/shiny-server/GAP/db/id_mapping/', species_db_name, '.Rdata'))
		
		if("Accession" %in% names(idmapping) == FALSE){
		
			idmapping$Accession <- idmapping$UniProtEntry
			
		}else{
		  
		  idmapping$UniProtEntry <- idmapping$Accession
		  
		}
	  
	  if("Protein names" %in% names(idmapping)){
	    
	    idmapping$Protein.Name <- idmapping$`Protein names`
	    
	  }
	
		match_result <- idmapping[which(idmapping$UniProtEntry %in% match_acc), ] %>% as.data.frame
		
		match_result <- match_result[, c('Accession', 'UniProtEntry', 'Protein.Name')]
		
	} else {
	
		match_result <- data.frame('Pick' = '', 'Accession' = '', 'UniProtEntry' = '', 'Protein.Name' = '') %>% .[-1, ]
	
	}
	
	return_value <- list()
	
	return_value[['input_list']] <- input_acc
	
	return_value[['annot_table']] <- match_result
	
	return_value[['matched']] <- length(match_acc)
	
	return_value[['bkground']] <- length(db_acc)
	
	return_value[['input']] <- length(input_acc)

	return(return_value)

}

generate_protein_sample <- function(species_db_name){

	db_acc <- paste0('/home/biognosis/shiny-server/GAP/db/acc_list/', species_db_name, '.txt') %>% read_csv(col_names = FALSE) %>% unlist
	
	db_acc[sample.int(length(db_acc), 100, replace = TRUE)] %>% accid_process %>% paste0(collapse = '\n') %>% return

}

Annot <- function (dr, species){

	load(file=paste("/home/biognosis/database/processed/",species,".Rdata",sep=""))
	
	if('Accession' %in% names(idmapping) == FALSE){
	
		idmapping$Accession <- idmapping$UniProtEntry
	
	}
	
	math_index <- sapply(dr, match, idmapping$Accession)
	
	dbs <- names(idmapping)[which(names(idmapping) %in% c("Accession", "UniProtEntry.Newest", "Status", "Length", "Mass", "Accession") == FALSE)]
	
	dr <- idmapping[, c('Accession', dbs), with = F] %>% as.data.frame

	return(dr)

}

generate_checkbox_fun <- function (row_index){

	shiny::checkboxInput(inputId = paste0('checkbox_', row_index), label = NULL, value = TRUE) %>% paste0 %>% return

}