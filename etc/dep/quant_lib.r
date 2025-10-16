# 主界面 ----
output$quant_panel <- renderUI({
  
  box(title = 'Metabonomics analysis', solidHeader = TRUE, status = "primary", width = 12, style = 'margin: 0 auto; min-height: 500px;',
      
      box(title = "Data upload", solidHeader = TRUE, status = "primary", width = 12, style = '', collapsible = T, collapsed = F, id = "upload_box", 
          
          div(style = 'max-width: 85%; margin: 3em auto;',
              
              checkboxInput("single_ion_judge", "Whether to use single-ion mode",value = FALSE),
              
              uiOutput("upload_quant_o")
              
          )
          
      ),
      
      box(title = "Data pre-processing", solidHeader = TRUE, status = "primary", width = 12, style = '', collapsible = T, collapsed = T, id = "preprocess_box", 
          
          div(style = 'max-width: 85%; margin: 3em auto;',
              
              selectInput(inputId = "miss.value.handle.type", "Types of missing values:", c(0, "NaN", "filtered", "NA")),
              
              selectInput(inputId = "miss.value.handle.group",
                          "Filtering missing values method:",
                          c("Within the group" = "inter.group",
                            "Global filter" = "global.group")
              ),
              
              numericInput(inputId = "miss.value.handle.cutoff",
                           "Missing value filtering threshold:",
                           value=50,min=0,max=90,step = 10
              ),
              
              selectInput(inputId = "miss.value.fill",
                          "Filling method：",
                          c("Intra-group mean" = "mean.group", 
                            "Global mean" = "mean.global", 
                            "Median within the group" = "median.group",
                            "Global median" = "median.global",
                            "None" = "none",
                            "Min" = "min.global",
                            "Min/2" = "min2",
                            "KNN" = "knn.global",
                            "Randomforest" = "rf")
              ),
              
              selectInput(inputId = "normalized.handle.method",
                          "Normalization method:",
                          c("Sum" = "sum",
                            "None" = "none",
                            "QC svr"= "qc",
                            "QC RLSC"= "qc-rlsc",
                            "QC NormAE"= "normae",
                            "Probability quotient" = "prob_quot",
                            "0.75 quantile" = "percent_0.75")
              ),
              
              uiOutput("sum_coef_o"),
              
              selectInput(inputId = "log.handle.method",
                          "The quantitative value is taken as log:",
                          c("None" = "none",
                            "log2" = "log2",
                            "log10" = "log10")
              ),
              
              numericInput(inputId = "rsd.cutoff", "QC RSD filter cutoff:", value=0.3, min=0, max=1, step = 0.1),
              
              selectInput(inputId = "log.rsd.method",
                          "Take the log from RSD:",
                          c("None" = "none",
                            "log2" = "log2",
                            "log10" = "log10")
              ),
              
              selectInput(inputId = "batch.correct",
                          "Batch correction：",
                          c("None" = "none",
                            "COMBAT" = "combat",
                            "Limma" = "limma")
              ),
              
              numericInput(inputId = "correlation.cutoff", "Correlation cutoff:", value=0.7,min=0,max=1,step = 0.1)
              
          ),
          
          div(style = 'width: 100%; text-align:center;',
              
              shinysky::actionButton(inputId = 'quant_apply_preprocess', label = 'Next', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-bottom: 0; margin-top: 3em;')
              
          ),
          
          div(style = 'width: 100%; text-align:center;',
              
              shinysky::actionButton(inputId = 'return_quant_data_upload', label = 'Return', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-top: 0.5em; margin-bottom: 30px;')
              
          )
          
      ),
      
      box(title = "Difference Cutoff", solidHeader = TRUE, status = "primary", width = 12, style = '', collapsible = T, collapsed = T, id = "diff_box", 
          
          div(style = 'max-width: 85%; margin: 3em auto;',
              
              checkboxInput("fc_cutoff_judge", "Whether to use Fold Change filtering",value = TRUE),
              uiOutput("fc_cutoff_o"),
              checkboxInput("p_cutoff_judge", "Whether to use pvalue filtering",value = TRUE),
              uiOutput("p_cutoff_o"),
              checkboxInput("padjust_judge", "Whether to use padjust filtering",value = TRUE),
              uiOutput("padjust_method_o"),
              checkboxInput("vip_cutoff_judge", "Whether to use VIP filtering",value = TRUE),
              uiOutput("vip_cutoff_o"),
              checkboxInput("hca_col_cluster", "Whether the samples are clustered in the clustering heatmap",value = FALSE),
              selectInput(inputId = "species", 
                          "Select species:", 
                          c("Human" = "hsa",
                            "Mouse" = "mmu",
                            "Rat" = "rno")),
              p('Project Name:', style = 'font-size: 14px; font-weight: bold;'),
              
              div(style = '',textInput('project_name', label = NULL, value = paste0('Project_',format(Sys.time(),"%Y%m%d%H%M")))),
              
              p('Organization Information:', style = 'font-size: 14px; font-weight: bold;'),
              
              div(style = '',textInput('organ_info', label = NULL)),
              
              p('E-mail Address:', style = 'font-size: 14px; font-weight: bold;'),
              
              div(style = '',textInput('email', label = NULL, value = '')),
              
              selectInput(inputId = "analysis_type", "Select the analysis type:", c("sjtu", "omicsolution"), selected = "omicsolution")
              
          ),
          
          div(style = 'width: 100%; text-align:center;', shinysky::actionButton(inputId = 'confirm_quant1', label = 'Confirm', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin-bottom: 0; margin-top: 3em;'))
          
      ),
      
      shinyjs::useShinyjs(),
      
      tags$script(HTML(
        "Shiny.addCustomMessageHandler('triggerResize', function(_) {
     setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 0);
   });"
      ))
      
  )
  
})

observe({
  
  if (!is.null(input$single_ion_judge)) {
    
    if (input$single_ion_judge) {
      
      output$upload_quant_o <- renderUI(
        tagList(
          h5(style = "color: red;", "Put up the table of the ionic mode on the flyer:"),
          fileInput("metaboQuant_single", paste0("Upload metabolites results. Control them within ",sizelimit," MB:"), multiple = FALSE, accept = c(".csv")),
          br(),
          DTOutput("metaboQuant_single_column"),
          uiOutput("metaboQuant_single_column_recognition"),
          br(),
          DTOutput("metaboQuant_single_group"),
          br(),
          uiOutput("metaboQuant_single_download_o"),
          br(),
          uiOutput("metaboQuant_compare")
        )
      )
      
    }else{
      
      output$upload_quant_o <- renderUI(
        tagList(
          h5(style = "color: red;", "Upload the positive ion mode table (After uploading and setting the groups, the negative ion table upload button will appear):"),
          fileInput("metaboQuant_pos", paste0("Upload metabolites POS results. Control them within ",sizelimit," MB:"), multiple = FALSE, accept = c(".csv")),
          br(),
          shinyBS::bsCollapse(
            id = "pos_sections",
            multiple = TRUE,
            open = c("pos_head_panel", "pos_key_panel"),
            shinyBS::bsCollapsePanel(
              title  = tagList(shiny::icon("info-circle"), "Positive ion meter header information"),
              value  = "pos_head_panel",
              style  = "primary",
              DTOutput("metaboQuant_pos_column")
            ),
            shinyBS::bsCollapsePanel(
              title  = tagList(shiny::icon("sliders-h"), "Set POS key columns"),
              value  = "pos_key_panel",
              style  = "danger",                      # ← 面板头变红
              div(class = "pos-key",                  # ← 给内容加类，方便定制
                  uiOutput("metaboQuant_pos_column_recognition")
              )
            )
          ),
          br(),
          DTOutput("metaboQuant_pos_group"),
          br(),
          uiOutput("metaboQuant_pos_download_o"),
          br(),
          uiOutput("neg_upload_panel"),
          br(),
          uiOutput("metaboQuant_compare")
        )
      )
      
    }
    
  }
  
})

observeEvent(rv$pos_uploaded, {
  if (isTRUE(rv$pos_uploaded)) {
    shinyBS::updateCollapse(session, "pos_sections", open = c("pos_head_panel","pos_key_panel"))
  }
}, ignoreInit = TRUE)

# CSS样式
tags$style("
  #quant_species_selector_box, #quant_project_confirm_box {
    transition: all 0.3s ease;
    overflow: hidden;
  }
")

tags$style(HTML("
  /* POS关键列：整体红色点缀 + 紧凑 */
  .pos-key .control-label { color: #d9534f; font-weight: 600; font-size: 12px; }
  .pos-key .form-group { margin-bottom: 8px; }
  .pos-key .form-control { height: 30px; padding: 2px 8px; font-size: 12px; }
  .pos-key .btn { padding: 6px 10px; }
"))

rv <- reactiveValues(
  pos_data = NULL,
  pos_uploaded = FALSE,
  neg_data = NULL,
  neg_uploaded = FALSE,
  neg_valid = FALSE,
  neg_alert_active = FALSE,
  neg_last_alert = NULL
)

callback <- c(
  "var tbl = $(table.table().node());",
  "var id = tbl.closest('.datatables').attr('id');",
  "table.on('autoFill', function(e, datatable, cells){",
  "  var out = [];",
  "  for(var i=0; i<cells.length; ++i){",
  "    var cells_i = cells[i];",
  "    for(var j=0; j < cells_i.length; ++j){",
  "      var c = cells_i[j];",
  "      var value = c.set === null ? '' : c.set;", # null causes problem in R
  "      out.push({row: c.index.row+1, col: c.index.column, value: value});",
  # if you want to color the autofilled cells, uncomment the the two lines below  
  #  "      $(table.cell(c.index.row, c.index.column).node())",
  #  "        .css('background-color', 'yellow');",
  "    }",
  "  }",
  "  Shiny.setInputValue(id + '_cells_filled:DT.cellInfo', out);",
  "  table.rows().invalidate();", # this updates the column type
  "});"
)

# 标签后加红星：用于“必填列”
req_label <- function(txt) {
  htmltools::HTML(sprintf('%s <span style="color:#d9534f;">*</span>', txt))
}

# 带最小值 1 的必填 numericInput（避免 0）
req_num <- function(inputId, label, value) {
  numericInput(inputId, label = req_label(label), value = value, min = 1)
}

# === 工具：安全取名 ===
get_colname_safe <- function(cols, idx) {
  if (is.null(idx) || is.na(idx) || idx <= 0 || idx > length(cols)) return(NA_character_)
  cols[idx]
}

# === 放宽规则的校验：仅校验样本列数量 + 关键列列名一致（可选列：仅在 POS 指定时要求一致） ===
validate_pos_neg_and_prepare <- function() {
  # 防抖用到的状态位（若未初始化，这里兜底）
  if (is.null(rv$neg_alert_active)) rv$neg_alert_active <- FALSE
  if (is.null(rv$neg_last_alert))   rv$neg_last_alert   <- NULL
  
  # 前置：两端都已上传 + POS 已设置样本起止
  if (!isTRUE(rv$pos_uploaded) || !isTRUE(rv$neg_uploaded)) return(FALSE)
  s1 <- input$key_column_sample_start
  e1 <- input$key_column_sample_end
  if (is.null(s1) || is.null(e1) || is.na(s1) || is.na(e1)) return(FALSE)
  
  posN <- ncol(rv$pos_data); negN <- ncol(rv$neg_data)
  
  # 1) 样本区间合法（在 POS）且 NEG 列数足够覆盖该区间
  bad_range <- function(s, e, N) is.na(s) || is.na(e) || s < 1 || e > N || s > e
  if (bad_range(s1, e1, posN)) {
    alert_once("Invalid sample range (POS)",
               sprintf("POS sample range must satisfy 1 ≤ start ≤ end ≤ %d.", posN),
               key = sprintf("pos-range-%s-%s", s1, e1))
    rv$neg_valid <- FALSE; return(FALSE)
  }
  if (e1 > negN) {
    alert_once("NEG has fewer columns",
               sprintf("NEG has only %d columns but POS sample end is %d.", negN, e1),
               key = sprintf("neg-too-short-%d-%d", negN, e1))
    rv$neg_valid <- FALSE; return(FALSE)
  }
  
  # 2) 样本列数量一致（使用 POS 的 [s1:e1] 作为标准）
  n_sample <- e1 - s1 + 1
  # 由于我们强制 s2=e1、s1 同步，此处只需保证 NEG 足够长，上面已检查；这里给信息提示用
  if (n_sample <= 0) {
    alert_once("Empty sample range", "Sample range length is 0. Please check.", key = "empty-sample")
    rv$neg_valid <- FALSE; return(FALSE)
  }
  
  # 3) 关键列（必填）列名一致：peak_name / mz / rt
  pos_cols <- colnames(rv$pos_data)
  neg_cols <- colnames(rv$neg_data)
  
  req_ids <- c(
    peak = input$key_column_peakname,
    mz   = input$key_column_mz,
    rt   = input$key_column_rt
  )
  # 可选列：只在 POS 指定时要求一致
  opt_ids <- c(
    name = input$key_column_name,
    kegg = input$key_column_kegg,
    hmdb = input$key_column_hmdb,
    lmsd = input$key_column_lmsd
  )
  
  # 检查必填列名一致
  bad_required <- character(0)
  for (nm in names(req_ids)) {
    idx <- req_ids[[nm]]
    pos_nm <- get_colname_safe(pos_cols, idx)
    neg_nm <- get_colname_safe(neg_cols, idx)
    if (is.na(idx) || idx == 0 || is.na(pos_nm) || is.na(neg_nm) || !identical(pos_nm, neg_nm)) {
      bad_required <- c(bad_required, sprintf("%s (POS: %s | NEG: %s)",
                                              nm, ifelse(is.na(pos_nm),"NA",pos_nm), ifelse(is.na(neg_nm),"NA",neg_nm)))
    }
  }
  if (length(bad_required) > 0) {
    alert_once("Key-column header mismatch",
               paste(c("Required key columns must match by header name at the same positions:",
                       paste0("- ", bad_required)), collapse = "\n"),
               key = paste("req-mis", paste(bad_required, collapse=";")))
    rv$neg_valid <- FALSE; return(FALSE)
  }
  
  # 检查可选列：仅对 POS 已指定（非 NA/非 0）的项目要求一致
  bad_optional <- character(0)
  for (nm in names(opt_ids)) {
    idx <- opt_ids[[nm]]
    if (!is.null(idx) && !is.na(idx) && idx > 0) {
      pos_nm <- get_colname_safe(pos_cols, idx)
      neg_nm <- get_colname_safe(neg_cols, idx)
      if (is.na(neg_nm) || !identical(pos_nm, neg_nm)) {
        bad_optional <- c(bad_optional, sprintf("%s (POS: %s | NEG: %s)",
                                                nm, ifelse(is.na(pos_nm),"NA",pos_nm), ifelse(is.na(neg_nm),"NA",neg_nm)))
      }
    }
  }
  if (length(bad_optional) > 0) {
    alert_once("Optional key header mismatch",
               paste(c("Optional keys that are set in POS must match in NEG:",
                       paste0("- ", bad_optional)), collapse = "\n"),
               key = paste("opt-mis", paste(bad_optional, collapse=";")))
    rv$neg_valid <- FALSE; return(FALSE)
  }
  
  # 4) ✅ 通过：构建 NEG 分组（Condition 同步 POS；否则 NA）
  s2 <- s1; e2 <- e1
  n_samp <- e2 - s2 + 1
  copyable <- !is.null(rv$condition.data.pos) &&
    nrow(rv$condition.data.pos) == n_samp &&
    "Condition" %in% colnames(rv$condition.data.pos)
  cond_vec <- if (copyable) rv$condition.data.pos$Condition else rep(NA, n_samp)
  if (!copyable) {
    alert_once("Condition not copied",
               "POS Condition is not ready or length mismatched, NEG Condition left as NA.",
               key = paste("cond-na", n_samp))
  }
  
  rv$condition.data.neg <- data.frame(
    `File name`      = neg_cols[s2:e2],
    Condition        = cond_vec,
    batch            = 1L,
    injection.order  = seq_len(n_samp),
    check.names      = FALSE
  )
  
  rv$neg_valid      <- TRUE
  rv$neg_last_alert <- NULL
  TRUE
}

# 防重复弹窗（同错不重复；正在显示时不叠加）
alert_once <- function(title, text, type = "error", key = NULL) {
  if (isTRUE(rv$neg_alert_active)) return(invisible(FALSE))
  if (!is.null(key) && identical(rv$neg_last_alert, key)) return(invisible(FALSE))
  rv$neg_alert_active <- TRUE
  rv$neg_last_alert   <- key
  shinyalert::shinyalert(
    title = title, text = text, type = type, closeOnClickOutside = TRUE,
    callbackR = function(x) rv$neg_alert_active <- FALSE
  )
  invisible(TRUE)
}

observeEvent(input$pos_sections, {
  session$sendCustomMessage("triggerResize", list())
}, ignoreInit = TRUE)

# pos ----
observeEvent(input$metaboQuant_pos, {
  req(input$metaboQuant_pos$datapath)
  rv$pos_data     <- readr::read_csv(input$metaboQuant_pos$datapath, show_col_types = FALSE)
  rv$pos_uploaded <- TRUE
  
  output$metaboQuant_pos_column <- DT::renderDT({
    DT::datatable(
      data.frame(`column name` = colnames(rv$pos_data),
                 `column number` = seq_along(colnames(rv$pos_data))),
      rownames = FALSE,
      caption = "Positive ion meter header information",
      options = list(paging = FALSE, scrollY = "500px", scrollCollapse = TRUE, fixedHeader = TRUE)
    )
  })
  outputOptions(output, "metaboQuant_pos_column", suspendWhenHidden = FALSE)
  
  `%or%` <- function(x, d) if (length(x) && !is.na(x[1])) x[1] else d
  
  output$metaboQuant_pos_column_recognition <- renderUI({
    req(rv$pos_uploaded, rv$pos_data)
    cols <- colnames(rv$pos_data)
    num_idx <- which(vapply(rv$pos_data, is.numeric, logical(1)))
    end_default <- if (length(num_idx)) max(num_idx) else ncol(rv$pos_data)
    
    `%or%` <- function(x, d) if (length(x) && !is.na(x[1])) x[1] else d
    
    tagList(
      fluidRow(
        column(4, req_num("key_column_peakname", "peak_name (required)",
                          which(grepl("^peak_name", cols, TRUE)) %or% 1)),
        column(4, req_num("key_column_mz", "mz (required)",
                          which(grepl("^mz$", cols, TRUE)) %or% 1)),
        column(4, req_num("key_column_rt", "rt (required)",
                          which(grepl("^rt$", cols, TRUE)) %or% 1))
      ),
      fluidRow(
        column(4, numericInput("key_column_name",  "name (optional)",
                               which(grepl("^name$", cols, TRUE)) %or% NA_integer_)),
        column(4, numericInput("key_column_kegg",  "kegg (optional)",
                               which(grepl("id_kegg|kegg", cols, TRUE)) %or% NA_integer_)),
        column(4, numericInput("key_column_hmdb",  "hmdb (optional)",
                               which(grepl("id_hmdb|hmdb", cols, TRUE)) %or% NA_integer_))
      ),
      fluidRow(
        column(4, numericInput("key_column_lmsd",  "lmsd (optional)",
                               which(grepl("Accepted Compound ID", cols, fixed = TRUE)) %or% NA_integer_)),
        column(4, req_num("key_column_sample_start", "sample start (required)", 1)),
        column(4, req_num("key_column_sample_end",   "sample end (required)", end_default))
      ),
      div(style = "margin-top:10px;",
          actionButton("key_column_confirm", "Confirm the key column",
                       class = "btn btn-danger"))  # ← 红色按钮
    )
  })
  
  outputOptions(output, "metaboQuant_pos_column_recognition", suspendWhenHidden = FALSE)
  
  shinyBS::updateCollapse(session, "pos_sections", open = c("pos_head_panel", "pos_key_panel"))
  session$sendCustomMessage("triggerResize", list())
})

observeEvent(input$key_column_confirm, {
  
  # 必填列集合
  required_vals <- c(input$key_column_peakname,
                     input$key_column_mz,
                     input$key_column_rt,
                     input$key_column_sample_start,
                     input$key_column_sample_end)
  
  # 判空/为0/NA
  if (any(is.na(required_vals) | required_vals == 0)) {
    shinyalert::shinyalert(
      title = "Required fields missing",
      text  = "Please fill in the mandatory columns marked with red asterisks.",
      type  = "warning",
      closeOnClickOutside = F
    )
    return()
  }
  
  if ((is.na(input$key_column_name)|input$key_column_name==0)&(is.na(input$key_column_kegg)&is.na(input$key_column_hmdb)&is.na(input$key_column_lmsd))) {
    
    shinyalert::shinyalert(
      title = "Important message",
      text = "If name column is unassigned, one of hmdb, lmsd, and kegg column must be assigned!",
      type = "warning",
      closeOnClickOutside = F
    )
    
  }else{
    
    rv$condition.data.pos <- data.frame(`File name`=colnames(rv$pos_data)[input$key_column_sample_start:input$key_column_sample_end],
                                        Condition=NA,
                                        batch = 1,
                                        injection.order = 1:(input$key_column_sample_end-input$key_column_sample_start+1),
                                        check.names = F)
    
    output$metaboQuant_pos_group <- renderDT(datatable(rv$condition.data.pos, 
                                                       selection = 'none', 
                                                       editable = list(target = "cell",disable = list(columns = c(0))),
                                                       extensions = 'AutoFill', 
                                                       callback = DT::JS(callback), 
                                                       rownames = F,
                                                       caption = "Set the grouping information",
                                                       options = list(
                                                         columnDefs = list(list(className = 'dt-left', targets = '_all')), 
                                                         pageLength = 20,
                                                         # lengthMenu = c(30,50,100),
                                                         lengthChange = F,
                                                         autoFill = list(columns = c(1), focus = 'click')
                                                       )),server = F)
    
    event_trigger <- reactive({
      
      list(input[["metaboQuant_pos_group_cells_filled"]], input[["metaboQuant_pos_group_cell_edit"]])
      
    })
    
    observeEvent(event_trigger(), ignoreInit = TRUE,{
      
      info <- rbind(input[["metaboQuant_pos_group_cells_filled"]],
                    input[["metaboQuant_pos_group_cell_edit"]])
      
      if(!is.null(info)){
        
        info <- unique(info)
        info$value[info$value==""] <- NA
        
        # 确保为赋值操作使用唯一的行索引
        info <- info[!duplicated(info$row), ]
        
        rv$condition.data.pos <- editData(rv$condition.data.pos,info,rownames = F)
        
      }
      
    })
    
    output$metaboQuant_pos_download_o <- renderUI(
      tagList(
        h5("You can choose to download the grouped data, edit it and then upload it"),
        
        fluidRow(
          downloadButton("condition.table.download.pos", "Download POS Condition Data")
        ),
        
        br(),
        
        checkboxInput("show_upload_pos", "Upload the edited group data file", value = FALSE),
        
        conditionalPanel(
          condition = "input.show_upload_pos == true",
          fluidRow(
            fileInput(
              "upload.condition.table.pos",
              label = 'Upload POS Condition file (only receive .tsv/.csv file)',
              accept = c('.tsv','.csv')
            )
          )
        )
      )
    )
    
    # 点击完正离子的确认按钮后出现负离子文件上传界面
    output$neg_upload_panel <- renderUI({
      # 只有在 POS 已完成关键列确认后才会进入到这个 renderUI（你原逻辑如此）
      # 本面板里也要避免“未上传 NEG 时显示提示文案”
      tagList(
        h5(style = "color: red;", "Upload the table of the negative ion mode:"),
        fileInput(
          "metaboQuant_neg", 
          paste0("Upload metabolites NEG results. Control them within ", sizelimit, " MB:"),
          multiple = FALSE, accept = c(".csv")
        ),
        
        # 若未上传 NEG：不展示任何标题/分组 UI/表头折叠（避免“字样”先出现）
        if (is.null(input$metaboQuant_neg) || is.null(input$metaboQuant_neg$datapath)) {
          NULL
        } else {
          tagList(
            # === 可折叠的负离子表头信息（默认折叠） ===
            shinyBS::bsCollapse(
              id = "neg_sections",
              multiple = TRUE,
              open = NULL,   # ⬅️ 默认不展开
              shinyBS::bsCollapsePanel(
                title = tagList(shiny::icon("info-circle"), "Negative ion meter head information"),
                value = "neg_head_panel",
                style = "info",
                DTOutput("metaboQuant_neg_column")   # 负离子表头表格（上传后再渲染）
              )
            ),
            br(),
            uiOutput("neg_grouping_panel")
          )
        }
      )
    })
    
  }
  
})

output$condition.table.download.pos <- downloadHandler(
  
  filename = function() {
    
    paste0("Metabolite_Group_POS_", Sys.Date(), ".csv")
    
  },
  
  content = function(file) {
    
    write_csv(rv$condition.data.pos,file)
    
  }
  
)

observe({
  
  if (!is.null(input$upload.condition.table.pos$datapath)) {
    
    # rv$condition.data.pos <- read_tsv(input[["upload.condition.table.pos"]][['datapath']],quoted_na = FALSE)
    
    # 获取文件扩展名并转换为小写
    file_ext <- tolower(tools::file_ext(input$upload.condition.table.pos$name))
    
    # 根据扩展名选择读取函数
    condition.data.pos <- if (file_ext == "csv") {
      read_csv(
        input[["upload.condition.table.pos"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 关闭列类型提示 
      )
    } else if (file_ext == "tsv") {
      read_tsv(
        input[["upload.condition.table.pos"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 保持与 CSV 一致的参数 
      )
    } else {
      stop("Unsupported file format. Only CSV/TSV allowed")
    }
    
    colnames(condition.data.pos) <- c("File name", "Condition", "batch", "injection.order")
    
    rv$condition.data.pos <- condition.data.pos
    
    event_trigger <- reactive({
      
      list(input[["metaboQuant_pos_group_cells_filled"]], input[["metaboQuant_pos_group_cell_edit"]])
      
    })
    
    observeEvent(event_trigger(), ignoreInit = TRUE,{
      
      info <- rbind(input[["metaboQuant_pos_group_cells_filled"]],
                    input[["metaboQuant_pos_group_cell_edit"]])
      
      if(!is.null(info)){
        
        info <- unique(info)
        info$value[info$value==""] <- NA
        rv$condition.data.pos <- editData(rv$condition.data.pos,info,rownames = F)
        
        # Create a new Replicate column without modifying the original data
        new_replicate <- rv$condition.data.pos$Replicate
        
        if (is.null(new_replicate)) new_replicate <- integer(nrow(rv$condition.data.pos))
        
        # identify which rows had their condition changed
        changed_rows <- info$row[info$col == match("Condition", colnames(rv$condition.data.pos))]
        
        # for each of those rows, update the Replicate column
        if (length(changed_rows) > 0) {
          conds <- unique(rv$condition.data.pos[changed_rows, "Condition"])
          for(cond in conds) {
            indices <- which(rv$condition.data.pos$Condition == cond)
            new_replicate[indices] <- seq_along(indices)
          }
        }
        
        # update the Replicate column in the original data
        rv$condition.data.pos$Replicate <- new_replicate
      }
      
    })
    
  }
  
})

# neg ----
output$neg_grouping_panel         <- renderUI(NULL)  # ������ 整块（标题+内容）的插槽
output$metaboQuant_neg_group      <- DT::renderDT(NULL)
output$metaboQuant_neg_download_o <- renderUI(NULL)
outputOptions(output, "neg_grouping_panel",         suspendWhenHidden = FALSE)
outputOptions(output, "metaboQuant_neg_group", suspendWhenHidden = FALSE)
outputOptions(output, "metaboQuant_neg_download_o", suspendWhenHidden = FALSE)

output$neg_grouping_panel <- renderUI({
  if (!isTRUE(rv$neg_uploaded)) return(NULL)
  tagList(
    h5("Set the grouping information:"),
    if (isTRUE(rv$neg_valid)) DTOutput("metaboQuant_neg_group") else div(em("Waiting for POS/NEG headers to align...")),
    br(),
    if (isTRUE(rv$neg_valid)) uiOutput("metaboQuant_neg_download_o") else NULL
  )
})

# 读取 NEG（只触发一次，不要包其它 req 导致不触发）
observeEvent(input$metaboQuant_neg, {
  req(input$metaboQuant_neg$datapath)
  rv$neg_data     <- readr::read_csv(input$metaboQuant_neg$datapath, show_col_types = FALSE)
  rv$neg_uploaded <- TRUE
  rv$neg_valid    <- FALSE      # 重新上传后需要重新校验
  
  output$metaboQuant_neg_column <- DT::renderDT({
    req(rv$neg_uploaded, rv$neg_data)
    cols <- colnames(rv$neg_data)
    DT::datatable(
      data.frame(
        `column name`   = cols,               # 先列名
        `column number` = seq_along(cols),    # 再列号（第2列）
        check.names     = FALSE
      ),
      rownames = FALSE,
      caption = "Negative ion meter header information",
      options = list(
        paging = FALSE,
        scrollY = "500px",
        scrollCollapse = TRUE,
        fixedHeader = TRUE,
        searching = FALSE,
        lengthChange = FALSE
      )
    )
  })
  
  outputOptions(output, "metaboQuant_neg_column", suspendWhenHidden = FALSE)
  
  # 上传后：强制展开负离子表头折叠 + 触发一次 resize（防 0 宽）
  shinyBS::updateCollapse(session, "neg_sections", open = "neg_head_panel")
  session$sendCustomMessage("triggerResize", list())
  
  # 紧接着做一次与 POS 的一致性校验；通过则准备分组数据
  rv$neg_valid <- validate_pos_neg_and_prepare()   # 用你当前无 *2 的精简版
})

observeEvent(input$neg_sections, { session$sendCustomMessage("triggerResize", list()) }, ignoreInit = TRUE)

valid_trigger <- reactive({
  list(input$key_column_sample_start, input$key_column_sample_end,
       input$metaboQuant_pos, input$metaboQuant_neg)
})

valid_trigger_debounced <- shiny::debounce(valid_trigger, millis = 300)

observeEvent(valid_trigger_debounced(), {
  if (isTRUE(rv$pos_uploaded) && isTRUE(rv$neg_uploaded)) {
    rv$neg_valid <- validate_pos_neg_and_prepare()
  }
}, ignoreInit = TRUE)

output$condition.table.download.neg <- downloadHandler(
  
  filename = function() {
    
    paste0("Metabolite_Group_neg_", Sys.Date(), ".csv")
    
  },
  
  content = function(file) {
    
    write_csv(rv$condition.data.neg,file)
    
  }
  
)

observe({
  
  if (!is.null(input$upload.condition.table.neg$datapath)) {
    
    # 获取文件扩展名并转换为小写
    file_ext <- tolower(tools::file_ext(input$upload.condition.table.neg$name))
    
    # 根据扩展名选择读取函数
    condition.data.neg <- if (file_ext == "csv") {
      read_csv(
        input[["upload.condition.table.neg"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 关闭列类型提示 
      )
    } else if (file_ext == "tsv") {
      read_tsv(
        input[["upload.condition.table.neg"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 保持与 CSV 一致的参数 
      )
    } else {
      stop("Unsupported file format. Only CSV/TSV allowed")
    }
    
    colnames(condition.data.neg) <- c("File name", "Condition", "batch", "injection.order")
    
    rv$condition.data.neg <- condition.data.neg
    
    event_trigger <- reactive({
      
      list(input[["metaboQuant_neg_group_cells_filled"]], input[["metaboQuant_neg_group_cell_edit"]])
      
    })
    
    observeEvent(event_trigger(), ignoreInit = TRUE,{
      
      info <- rbind(input[["metaboQuant_neg_group_cells_filled"]],
                    input[["metaboQuant_neg_group_cell_edit"]])
      
      if(!is.null(info)){
        
        info <- unique(info)
        info$value[info$value==""] <- NA
        rv$condition.data.neg <- editData(rv$condition.data.neg,info,rownames = F)
        
        # Create a new Replicate column without modifying the original data
        new_replicate <- rv$condition.data.neg$Replicate
        
        # identify which rows had their condition changed
        changed_rows <- info$row[info$col == match("Condition", colnames(rv$condition.data.neg))]
        
        # for each of those rows, update the Replicate column
        if (length(changed_rows) > 0) {
          conds <- unique(rv$condition.data.neg[changed_rows, "Condition"])
          for(cond in conds) {
            indices <- which(rv$condition.data.neg$Condition == cond)
            new_replicate[indices] <- seq_along(indices)
          }
        }
        
        # update the Replicate column in the original data
        rv$condition.data.neg$Replicate <- new_replicate
      }
      
    })
    
  }
  
})

# 仅当校验通过时渲染“分组表 + 下载/上传面板”；失败/回退则清空内容
observeEvent(rv$neg_valid, {
  if (isTRUE(rv$neg_valid)) {
    req(rv$condition.data.neg, nrow(rv$condition.data.neg) > 0)
    
    output$metaboQuant_neg_group <- DT::renderDT(
      DT::datatable(
        rv$condition.data.neg,
        selection = 'none',
        editable  = list(target = "cell", disable = list(columns = c(0))),
        extensions = 'AutoFill',
        rownames  = FALSE,
        caption   = "Set the grouping information",
        options   = list(
          columnDefs   = list(list(className = 'dt-left', targets = '_all')),
          pageLength   = 20, lengthChange = FALSE,
          autoFill     = list(columns = c(1), focus = 'click')
        )
      ),
      server = FALSE
    )
    
    output$metaboQuant_neg_download_o <- renderUI({
      tagList(
        h5("You can download the grouped data, edit it and then upload it"),
        fluidRow(downloadButton("condition.table.download.neg", "Download NEG Condition Data")),
        br(),
        checkboxInput("show_upload_neg", "Upload the edited group data file", value = FALSE),
        conditionalPanel(
          condition = "input.show_upload_neg == true",
          fluidRow(fileInput("upload.condition.table.neg",
                             label = 'Upload NEG Condition file (only .tsv/.csv)',
                             accept = c('.tsv','.csv','.TSV','.CSV')))
        )
      )
    })
    
  } else {
    # 未通过/回退：清空内容（标题仍由 neg_grouping_panel 维持）
    output$metaboQuant_neg_group      <- DT::renderDT(NULL)
    output$metaboQuant_neg_download_o <- renderUI(NULL)
  }
}, ignoreInit = TRUE)

# single ----
observe({
  
  if (!is.null(input$metaboQuant_single)) {
    
    metaboQuant_single_data <<- read_csv(input[["metaboQuant_single"]][['datapath']])
    
    output$metaboQuant_single_column <- renderDT({
      
      datatable(
        data.frame(`column name`=colnames(metaboQuant_single_data), 
                   `column number`=1:ncol(metaboQuant_single_data)),
        rownames = F,
        caption = "Single-ion meter head information",
        # 添加垂直滚动条配置
        options = list(
          paging = FALSE,
          scrollY = "500px",       # 设置固定高度（300px）
          scrollCollapse = TRUE,    # 内容不足时自动隐藏滚动条
          fixedHeader = TRUE
        )
      )
      
    })
    
    output$metaboQuant_single_column_recognition <- renderUI({
      
      tagList(
        h5("选择关键列"),
        column(width = 4, req_num("key_column_peakname", "The identified peak_name column", grep("^peak_name", colnames(metaboQuant_single_data))[1])),
        column(width = 4, req_num("key_column_mz", "The identified mz column", grep("^mz", colnames(metaboQuant_single_data))[1])),
        column(width = 4, req_num("key_column_rt", "The identified rt column", grep("^rt", colnames(metaboQuant_single_data))[1])),
        column(width = 4, numericInput("key_column_name", "The identified name column (optional)", value = grep("^name", colnames(metaboQuant_single_data))[1])),
        column(width = 4, numericInput("key_column_kegg", "The identified id_kegg column（optional）", value = grep("^id_kegg", colnames(metaboQuant_single_data))[1])),
        column(width = 4, numericInput("key_column_hmdb", "The identified id_hmdb column (optional)", value = grep("^id_hmdb", colnames(metaboQuant_single_data))[1])),
        column(width = 4, numericInput("key_column_lmsd", "The identified lmsd column (optional)", value = grep("Accepted Compound ID", colnames(metaboQuant_single_data))[1])),
        column(width = 4, req_num("key_column_sample_start", "Please select the sample quantitative value to start the column", 1)),
        column(width = 4, req_num("key_column_sample_end", "Please select the end column of the sample quantitative value", which(sapply(metaboQuant_single_data, is.numeric)) %>% max)),
        column(width = 12, actionButton("key_column_confirm_single", "Confirm the key column"))
      )
      
    })
    
  }
  
})

observeEvent(input$key_column_confirm_single, {
  
  # 必填列集合
  required_vals <- c(input$key_column_peakname,
                     input$key_column_mz,
                     input$key_column_rt,
                     input$key_column_sample_start,
                     input$key_column_sample_end)
  
  # 判空/为0/NA
  if (any(is.na(required_vals) | required_vals == 0)) {
    shinyalert::shinyalert(
      title = "Required fields missing",
      text  = "Please fill in the mandatory columns marked with red asterisks.",
      type  = "warning",
      closeOnClickOutside = F
    )
    return()
  }
  
  if ((is.na(input$key_column_name)|input$key_column_name==0)&(is.na(input$key_column_kegg)&is.na(input$key_column_hmdb)&is.na(input$key_column_lmsd))) {
    
    shinyalert::shinyalert(
      title = "Important message",
      text = "If name column is unassigned, one of hmdb, lmsd, and kegg column must be assigned!",
      type = "warning",
      closeOnClickOutside = F
    )
    
  }else{
    
    rv$condition.data.single <- data.frame(`File name`=colnames(metaboQuant_single_data)[input$key_column_sample_start:input$key_column_sample_end],
                                           Condition=NA,
                                           batch = 1,
                                           injection.order = 1:(input$key_column_sample_end-input$key_column_sample_start+1),
                                           check.names = F)
    
    output$metaboQuant_single_group <- renderDT(datatable(rv$condition.data.single, 
                                                          selection = 'none', 
                                                          editable = list(target = "cell",disable = list(columns = c(0))),
                                                          extensions = 'AutoFill', 
                                                          callback = DT::JS(callback), 
                                                          rownames = F,
                                                          options = list(
                                                            columnDefs = list(list(className = 'dt-left', targets = '_all')), 
                                                            pageLength = 20,
                                                            # lengthMenu = c(30,50,100),
                                                            lengthChange = F,
                                                            autoFill = list(columns = c(1), focus = 'click')
                                                          )),server = F)
    
    event_trigger <- reactive({
      
      list(input[["metaboQuant_single_group_cells_filled"]], input[["metaboQuant_single_group_cell_edit"]])
      
    })
    
    observeEvent(event_trigger(), ignoreInit = TRUE,{
      
      info <- rbind(input[["metaboQuant_single_group_cells_filled"]],
                    input[["metaboQuant_single_group_cell_edit"]])
      
      if(!is.null(info)){
        
        info <- unique(info)
        info$value[info$value==""] <- NA
        
        # 确保为赋值操作使用唯一的行索引
        info <- info[!duplicated(info$row), ]
        
        rv$condition.data.single <- editData(rv$condition.data.single,info,rownames = F)
        
      }
      
    })
    
    output$metaboQuant_single_download_o <- renderUI(
      tagList(
        h5("You can choose to download the grouped data, edit it and then upload it"),
        fluidRow(downloadButton("condition.table.download.single", "Download Contition Data")),
        br(),
        fluidRow(fileInput("upload.condition.table.single", label = 'Upload Condition file (only recevice .tsv/.csv file)', accept = c('.tsv','.csv')))
      )
    )
    
  }
  
})

output$condition.table.download.single <- downloadHandler(
  
  filename = function() {
    
    paste0("Metabolite_Group_", Sys.Date(), ".csv")
    
  },
  
  content = function(file) {
    
    write_csv(rv$condition.data.single,file)
    
  }
  
)

observe({
  
  if (!is.null(input$upload.condition.table.single$datapath)) {
    
    # 获取文件扩展名并转换为小写
    file_ext <- tolower(tools::file_ext(input$upload.condition.table.single$name))
    
    # 根据扩展名选择读取函数
    condition.data.single <- if (file_ext == "csv") {
      read_csv(
        input[["upload.condition.table.single"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 关闭列类型提示 
      )
    } else if (file_ext == "tsv") {
      read_tsv(
        input[["upload.condition.table.single"]][['datapath']],
        quoted_na = FALSE,
        show_col_types = FALSE  # 保持与 CSV 一致的参数 
      )
    } else {
      stop("Unsupported file format. Only CSV/TSV allowed")
    }
    
    colnames(condition.data.single) <- c("File name", "Condition", "batch", "injection.order")
    
    rv$condition.data.single <- condition.data.single
    
    event_trigger <- reactive({
      
      list(input[["metaboQuant_single_group_cells_filled"]], input[["metaboQuant_single_group_cell_edit"]])
      
    })
    
    observeEvent(event_trigger(), ignoreInit = TRUE,{
      
      info <- rbind(input[["metaboQuant_single_group_cells_filled"]],
                    input[["metaboQuant_single_group_cell_edit"]])
      
      if(!is.null(info)){
        
        info <- unique(info)
        info$value[info$value==""] <- NA
        rv$condition.data.single <- editData(rv$condition.data.single,info,rownames = F)
        
        # Create a new Replicate column without modifying the original data
        new_replicate <- rv$condition.data.single$Replicate
        
        # identify which rows had their condition changed
        changed_rows <- info$row[info$col == match("Condition", colnames(rv$condition.data.single))]
        
        # for each of those rows, update the Replicate column
        if (length(changed_rows) > 0) {
          conds <- unique(rv$condition.data.single[changed_rows, "Condition"])
          for(cond in conds) {
            indices <- which(rv$condition.data.single$Condition == cond)
            new_replicate[indices] <- seq_along(indices)
          }
        }
        
        # update the Replicate column in the original data
        rv$condition.data.single$Replicate <- new_replicate
      }
      
    })
    
  }
  
})

# compare select ----
observe({
  
  if (!is.null(input$single_ion_judge)) {
    
    if (input$single_ion_judge){
      
      req(rv$condition.data.single)
      
      unique_conditions <- na.omit(unique(rv$condition.data.single$Condition))
      
      # Create a data frame with required columns
      rv$comparison_data <- data.frame(
        Condition_Numerator = rep(unique_conditions, length(unique_conditions)),
        Condition_Denominator = rep(unique_conditions, each = length(unique_conditions)),
        stringsAsFactors = FALSE
      )
      
      output$metaboQuant_compare <- renderUI(
        tagList(
          selectInput(inputId = "metaboQuant_compare_select",
                      label = "Select the group for comparison",
                      choices = apply(rv$comparison_data[which(rv$comparison_data$Condition_Numerator != rv$comparison_data$Condition_Denominator),], 1, function(x) paste0(x, collapse = "_vs_")) %>% as.character(),
                      multiple = T),
          div(style = 'width: 100%; text-align:center;',
              shinysky::actionButton(inputId = 'quant_apply_upload', label = 'Next', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin: 2em 0;')
          )
        ))
      
    }else{
      
      req(rv$condition.data.pos)
      req(isTRUE(rv$neg_valid))
      req(rv$condition.data.neg)
      
      unique_conditions <- na.omit(unique(rv$condition.data.pos$Condition))
      
      # Create a data frame with required columns
      rv$comparison_data <- data.frame(
        Condition_Numerator = rep(unique_conditions, length(unique_conditions)),
        Condition_Denominator = rep(unique_conditions, each = length(unique_conditions)),
        stringsAsFactors = FALSE
      )
      
      output$metaboQuant_compare <- renderUI(
        tagList(
          selectInput(inputId = "metaboQuant_compare_select",
                      label = "Select the group for comparison",
                      choices = apply(rv$comparison_data[which(rv$comparison_data$Condition_Numerator != rv$comparison_data$Condition_Denominator),], 1, function(x) paste0(x, collapse = "_vs_")) %>% as.character(),
                      multiple = T),
          div(style = 'width: 100%; text-align:center;',
              shinysky::actionButton(inputId = 'quant_apply_upload', label = 'Next', styleclass = 'primary', style = 'padding: 6px 8px; width: 250px; margin: 2em 0;')
          )
        ))
      
    }
    
  }
  
})

# 条件panel ----
observe({
  
  if (!is.null(input$normalized.handle.method)) {
    
    if (input$normalized.handle.method=='sum') {
      
      output$sum_coef_o <- renderUI(
        tagList(
          numericInput("sum_coef", "The coefficient of Sum", value = 1, step=1)
        )
      )
      
    }else{
      
      output$sum_coef_o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$fc_cutoff_judge)) {
    
    if (input$fc_cutoff_judge) {
      
      output$fc_cutoff_o <- renderUI(
        tagList(
          numericInput("fc_cutoff","Fold Change Cutoff",value = 2, step=0.1)
        )
      )
      
    }else{
      
      output$fc_cutoff_o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$p_cutoff_judge)) {
    
    if (input$p_cutoff_judge) {
      
      output$p_cutoff_o <- renderUI(
        tagList(
          selectInput(inputId = "pvalue_type", "Select the type of test:", c("ttest", "wilcox_test")),
          numericInput("p_cutoff","P-value Cutoff",value = 0.05, step=0.01)
        )
      )
      
    }else{
      
      output$p_cutoff_o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$padjust_judge)) {
    
    if (input$padjust_judge) {
      
      output$padjust_method_o <- renderUI(
        tagList(
          selectInput(inputId = "padjust_method", "choose padjust method:", c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"), selected = "BH")
        )
      )
      
    }else{
      
      output$padjust_method_o <- NULL
      
    }
    
  }
  
})

observe({
  
  if (!is.null(input$vip_cutoff_judge)) {
    
    if (input$vip_cutoff_judge) {
      
      output$vip_cutoff_o <- renderUI(
        tagList(
          numericInput("vip_cutoff","VIP Cutoff",value = 1, step=0.01)
        )
      )
      
    }else{
      
      output$vip_cutoff_o <- NULL
      
    }
    
  }
  
})

# 提交数据 ----
observeEvent(input$quant_apply_upload, {
  
  if (input$single_ion_judge) {
    
    if (is.null(input[["metaboQuant_single"]][['datapath']])){
      
      shinyalert::shinyalert("Missing files !", "Please confirm data has been uploaded successfully.", 
                             type = "error",
                             callbackR = function(x) {
                               if(x != FALSE)
                                 runjs(paste0(
                                   "(function() {location.href='';}());"
                                 ))
                             })
      
    }
    
    if (sum(is.na(rv$condition.data.single$Condition))>0) {
      
      shinyalert::shinyalert("Missing files !", "Please confirm condition has been filled.", 
                             type = "error",
                             callbackR = function(x) {
                               if(x != FALSE)
                                 runjs(paste0(
                                   "(function() {location.href='';}());"
                                 ))
                             })
      
    }
    
  }else{
    
    if (!isTRUE(rv$neg_valid)) {
      shinyalert::shinyalert("NEG headers not aligned",
                             "Please make sure POS/NEG columns and non-sample headers are consistent.",
                             type = "error"
      )
      return()
    }
    
    if (is.null(input[["metaboQuant_pos"]][['datapath']]) | is.null(input[["metaboQuant_neg"]][['datapath']])){
      
      shinyalert::shinyalert("Missing files !", "Please confirm data has been uploaded successfully.", 
                             type = "error",
                             callbackR = function(x) {
                               if(x != FALSE)
                                 runjs(paste0(
                                   "(function() {location.href='';}());"
                                 ))
                             })
      
    }
    
    if (sum(is.na(rv$condition.data.pos$Condition))>0 | sum(is.na(rv$condition.data.neg$Condition))>0) {
      
      shinyalert::shinyalert("Missing files !", "Please confirm condition has been filled.", 
                             type = "error",
                             callbackR = function(x) {
                               if(x != FALSE)
                                 runjs(paste0(
                                   "(function() {location.href='';}());"
                                 ))
                             })
      
    }
    
  }
  
  # 通过 JavaScript 控制折叠
  runjs('
    // 折叠上传面板
    $("#upload_box").closest(".box").find(\'[data-widget="collapse"]\').click();
    // 展开预处理面板 
    $("#preprocess_box").closest(".box").find(\'[data-widget="collapse"]\').click();
  ')
  
})

# 提交预处理 ----
observeEvent(input$quant_apply_preprocess, {
  
  # 通过 JavaScript 控制折叠
  runjs('
    // 折叠上传面板
    $("#preprocess_box").closest(".box").find(\'[data-widget="collapse"]\').click();
    // 展开预处理面板 
    $("#diff_box").closest(".box").find(\'[data-widget="collapse"]\').click();
  ')
  
})

# confirm ----
observeEvent(input$confirm_quant1, {
  
  if (trimws(input$organ_info)==""|trimws(input$email)=="") {
    
    showModal(modalDialog(
      title = "Important message",
      paste0("Please fill in blank!"),
      easyClose = T
    ))
    
  }else{
    
    generate_jobid <- function() {
      
      a <- paste0(collapse = '', sample(x = c(letters, LETTERS, 0:9), size = 12, replace = TRUE))
      
      paste0(round(as.numeric(Sys.time()),0), a)
      
    }
    
    project_id <- generate_jobid()
    
    data_path <<- paste0(path_prefix2, '/data/usrdata/', project_id, '/data')
    
    para <<- list(); para[['info']] <<- list(para[['info']])
    
    para[['info']]['type'] <<- input$analysis_type
    
    para[['info']]['data_path'] <<- data_path
    
    para[['info']]['project_id'] <<- project_id
    
    dir.create(path = data_path, recursive = TRUE)
    
    if (input$single_ion_judge){
      
      file.copy(from = input$metaboQuant_single$datapath, to = paste0(data_path, '/metabonomics.csv'))
      
      write_csv(rv$condition.data.single, paste0(data_path, '/condition.csv'))
      
    }else{
      
      file.copy(from = input$metaboQuant_pos$datapath, to = paste0(data_path, '/metabonomics_POS.csv'))
      file.copy(from = input$metaboQuant_neg$datapath, to = paste0(data_path, '/metabonomics_NEG.csv'))
      
      write_csv(rv$condition.data.pos, paste0(data_path, '/condition_POS.csv'))
      write_csv(rv$condition.data.neg, paste0(data_path, '/condition_NEG.csv'))
      
    }
    
    inputs_list <- reactiveValuesToList(input)
    inputs_path <- paste0(data_path, '/input.rds')
    saveRDS(inputs_list, inputs_path)
    
    rv_list <- reactiveValuesToList(rv)
    rv_path <- paste0(data_path, '/rv.rds')
    saveRDS(rv_list, rv_path)
    
    save(list = c('para'), file = paste0(data_path, '/para.Rdata'), envir = .GlobalEnv)
    
    # system(paste0("cp -rf ", path_prefix2, "/db/projects_list.db ", path_prefix1, "/db/projects_list.db"))
    
    projects_list <- dbConnect(RSQLite::SQLite(), paste0(path_prefix1, "/db/projects_list.db"), flags = SQLITE_RWC)
    
    dbSendQuery(projects_list, paste0("INSERT INTO Project VALUES ('", project_id, "','", input$project_name, "', '", paste0(input$organ_info,"_",user_info()$name), "', '0', '", as.character(Sys.time()),"', 'metabonomics', '",input$email,"')"))
    
    dbDisconnect(projects_list)
    
    # system(paste0("cp -rf ", path_prefix1, "/db/projects_list.db ", path_prefix2, "/db/projects_list.db"))
    
    # paste0("'/usr/lib/R/bin/Rscript' '", path_prefix1, "/etc/dep/bg_script.r' ", project_id) %>% system(wait = F)
    
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

# 返回按钮逻辑
observeEvent(input$return_quant_data_upload, {
  # 通过 JavaScript 控制折叠
  runjs('
    // 展开上传面板
    $("#upload_box").closest(".box").find(\'[data-widget="collapse"]\').click();
    // 折叠预处理面板
    $("#preprocess_box").closest(".box").find(\'[data-widget="collapse"]\').click();
  ')
})