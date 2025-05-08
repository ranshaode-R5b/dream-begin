library("shiny")
#install.packages("shiny")
#install.packages("bslib")
#install.packages("htmltools")
#install.packages("rlang")
#install.packages("fastmap")
library(bslib)
#update.packages("shiny")
#update.packages("promises")
install.packages("htmlwidgets")
a <- read.table("/home/renshida/eqtm/data/data_cpg/29914364-GTP.txt",header=T)
a <- a[1:100,]


addResourcePath("manpage", "/home/renshida/eqtm/code")
ui <- page_navbar(
  title = strong("eQTMbase",style = "font-size: 1.3em;"),
  navbar_options = navbar_options(bg = "#2D89C8"),
  header = tags$head(
    tags$style(HTML(paste0(
      "/* 强制表头不换行 */",
      ".dataTable thead th {",
      "white-space: nowrap !important;", 
      "min-width: 120px !important;",    # 最小宽度
      "text-align: center !important;",
      "}",
      "/* 单元格内容不换行 */",
      ".dataTable tbody td {",
      "white-space: nowrap !important;",
      "text-align: center !important;",
      "}",
      "/* 调整行高 */",
      ".dataTable tbody tr {",
      "height: 43px !important;",
      "}"
    )))
    ),
  nav_panel(
    tags$b("Home"),
    fluidPage(
      # 主容器设置
      div(class = "main-container mt-3",
          # 标题区
          fluidRow(
            column(10,offset=1,div(
              tags$h3("Welcome to eQTMbase!", 
                      style = "font-weight: 700; color: #2D89C8;")
            ))
          ),
          
          # 文字
          fluidRow(
            column(10, offset = 1,
                   div(h4("Introduction", class = "mt-3"),
                       p("eQTMbase is a comprehensive website that aggregates human eQTM information, comprising data from 24 datasets and covering 11 distinct tissues. The platform supports searches by both CpG sites and genes, enabling users to selectively download data from various datasets. Additionally, eQTMbase provides multiple visualization outputs for these datasets and annotates their biological functions using extensive genomic features and functional evidence.",
                         class = "lead text-justify")
                   )
            )
          ),
          
          # 主示意图
          fluidRow(
            column(10, offset = 1,
                   div(class = "mt-5 img-container",
                       img(src = "https://via.placeholder.com/1200x400",
                           class = "img-fluid rounded-lg shadow")
                   )
            )
          ),
          hr(class = "my-5", style = "border-top: 2px solid #2D89C8;width: 75%; margin: 0 auto;"),
          # 功能卡片区
          fluidRow(
            column(12,
                   fluidRow(
                     class = "justify-content-center",  # 水平居中
                     lapply(1:3, function(i) {
                       column(4, class = "mb-4",
                              div(class = "feature-card",
                                  img(src = paste0("https://via.placeholder.com/300x200?text=Feature",i),
                                      class = "img-fluid rounded-top"),
                                  div(class = "card-body",
                                      p(switch(i,
                                               "1" = "eQTMbase 相关说明",
                                               "2" = "功能说明",
                                               "3" = "交互分析"),
                                        class = "card-text text-center")
                                  )
                              )
                       )
                     })
                   )
            )
          ),
          
          # 参考文献卡片
          fluidRow(
            column(10, offset = 1,
                   div(class = "mt-5",
                       card(class = "citation-card",
                            card_header("Reference", class = "bg-primary text-white", style = "font-size: 1.3em;"),
                            card_body(
                              markdown("If you have used eQTM data from this website, please cite the following papers:"),
                              tags$blockquote(
                                class = "blockquote",
                                p("[1] Nature Genetics...",style = "font-size: 0.9em;")
                              )
                            )
                       )
                   )
            )
          ),
          hr(class = "my-5", style = "border-top: 2px solid #2D89C8;width: 75%; margin: 0 auto;"),
          #contact us
          fluidRow(
            column(10, offset = 1,
                   div(class = "mt-3",
                       card(
                         class = "citation-card",
                         card_header("Contact Us", class = "bg-primary text-white", style = "font-size: 1.3em;"),
                         card_body(
                           tags$blockquote(class = "blockquote",p("Our platform serves as a comprehensive platform for eQTM data. If your research has yielded new eQTM datasets and you're interested in sharing them with fellow researchers, we encourage you to reach out to us. Please feel free to contact us at xxx@.com to discuss how we can assist in disseminating your valuable research.",style = "font-size: 0.9em;")))
                       )
                   )
            )
          )
      )
    )
  ),
  
  
  # 页面2：Search
  nav_panel(
    tags$b("Search"),
    style = "margin-left: 40px;",
    fluidPage(
      titlePanel(h3("Search",style = "font-weight: 600;")),
      sidebarLayout(
        sidebarPanel(
          width=3,
          style = "background-color: #f8f9fa; border-right: 2px solid #dee2e6;",
          
          # 物种选择
          div(
            style = "margin-bottom: 25px;",
            selectInput(
              "species",
              label = tags$span(icon("person"), " Species"),  # 修正标签格式
              choices = c("Human"), 
              selected = "Human",
              width = "100%"
            )
          ),
          
          # 搜索模式
          div(
            style = "margin-bottom: 25px;",
            # 搜索模式
            selectInput("model", 
                        tags$span(icon("search"), " Search Mode"), 
                        choices = c("CpG", "Gene"), 
                        selected = "CpG",
                        width = "100%"
            )
          ),
          
          conditionalPanel(
            condition = "input.model == 'CpG'",
            div(
              style = "margin-bottom: 25px;",
              selectInput("feature", 
                          tags$span(icon("glasses"), " Methylation Type"), 
                          choices = c("CpG", "DMR"), 
                          selected = "CpG",
                          width = "100%")
            )
          ),
          
          uiOutput("search_input_ui"),
          
          div(
            style = "margin-top: 30px; text-align: center;",  # 按钮居中
            actionButton("search_btn", "Search", 
                         class = "btn-primary",
                         style = "width: 100%; 
                                background-color: #2D89C8;
                                border-color: #2478ab;
                                padding: 10px 20px;
                                font-weight: 600;"
            )
          )
        ),
        
        mainPanel(
          # 下载按钮
          div(
            style = "text-align: right;",
            downloadButton("download_btn", "Download CSV")
          ),
          
          # 动态显示结果表格
          conditionalPanel(
            condition = "input.model == 'CpG'",
            div(
              style = "box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
              DTOutput("search_result")
            )
          ),
          conditionalPanel(
            condition = "input.model == 'Gene'",
            tabsetPanel(
              tabPanel(
                title = tags$span(
                  style = "font-weight: 650;",
                  "CpG-related eQTMs"
                ),
                div(
                  style = "margin-top: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
                  DTOutput("search_result_cpg_gene")
                )
              ),
              tabPanel(
                title = tags$span(
                  style = "font-weight: 650;",
                  "DMR-related eQTMs"
                ),
                div(
                  style = "margin-top: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
                  DTOutput("search_result_dmr_gene")
                )
              )
            )
          )
        )
      )
    )
  ),
  
  nav_panel(
    tags$b("Manhattan"),
    style = "margin-left: 40px;",
    fluidPage(
      titlePanel(h3("Manhattan Plot")),
      plotOutput("manhattan_plot", height = "400px"),
      tags$iframe(
        src = "manpage/man.html",
        height = "600px",
        width = "100%",
        frameborder = "0"
      )
    )
  ),
  
  #Download page
  nav_panel(tags$b("Download"),
            style = "margin-left: 40px;",
            div(class = "main-container mt-3",
                fluidRow(
                  column(10,offset=1,h3("Download eQTMs"))
                ),
                fluidRow(
                  column(10, offset = 1,
                         div(h4("On this page you can download the eQTM data available in different datasets.",
                                class = "lead text-justify mt-3")
                         )
                  )
                ),
            )
  ),
  nav_spacer(),
  nav_menu(
    title = tags$b("Links"),
    align = "right",
    nav_item(tags$a("Help", href = "https://shiny.posit.co")),
    nav_item(tags$a("News", href = "https://shiny.posit.co")),
    nav_item(tags$a("Contact Us", href = "http://bioailab.com/people.html"))
  )
)





server <- function(input, output, session) {
  # 读取数据
  raw_data  <- data1[1:100, ]
  raw_data1 <- data2[1:100, ]
  raw_data2 <- data11[1:100, ]
  
  # 统一处理 "Tissue" 为因子列
  if ("Tissue" %in% colnames(raw_data)) {
    raw_data$Tissue <- as.factor(raw_data$Tissue)
  }
  if ("Tissue" %in% colnames(raw_data1)) {
    raw_data1$Tissue <- as.factor(raw_data1$Tissue)
  }
  if ("Tissue" %in% colnames(raw_data2)) {
    raw_data2$Tissue <- as.factor(raw_data2$Tissue)
  }
  
  # 合并 CpG 数据集（raw_data + raw_data1）
  combined_data1 <- reactiveVal(rbind(raw_data, raw_data1))
  
  # 动态搜索输入框
  output$search_input_ui <- renderUI({
    if (input$model == "CpG") {
      if (!is.null(input$feature) && input$feature == "CpG") {
        textInput("search_text", label = tags$span(icon("file"), " Enter CpG:"), placeholder = "e.g. cg00000029")
      } else {
        textInput("search_text", label = tags$span(icon("file"), " Enter DMR:"), placeholder = "e.g. chr1:12345-67890")
      }
    } else {
      textInput("search_text", label = tags$span(icon("file"), " Enter Gene:"), placeholder = "e.g. MCOLN2")
    }
  })
  
  # 搜索过滤逻辑
  filtered_data_cpg <- reactiveVal(NULL)
  filtered_data_dmr <- reactiveVal(NULL)
  
  observeEvent(input$search_btn, {
    if (input$model == "CpG") {
      if (input$feature == "CpG") {
        df <- combined_data1()
        if (!is.null(input$search_text) && nzchar(input$search_text)) {
          df <- df[grepl(input$search_text, df$`CpG site`, ignore.case = TRUE), ]
        }
        filtered_data_cpg(df)
        filtered_data_dmr(NULL)
      } else {
        df <- raw_data2
        if (!is.null(input$search_text) && nzchar(input$search_text)) {
          df <- df[grepl(input$search_text, df$DMR, ignore.case = TRUE), ]
        }
        filtered_data_dmr(df)
        filtered_data_cpg(NULL)
      }
    } else {
      df_cpg <- rbind(raw_data, raw_data1)
      df_dmr <- raw_data2
      if (!is.null(input$search_text) && nzchar(input$search_text)) {
        df_cpg <- df_cpg[grepl(input$search_text, df_cpg$`Gene Symbol`, ignore.case = TRUE), ]
        df_dmr <- df_dmr[grepl(input$search_text, df_dmr$`Gene Symbol`, ignore.case = TRUE), ]
      }
      filtered_data_cpg(df_cpg)
      filtered_data_dmr(df_dmr)
    }
  })
  
  # 添加超链接处理函数
  format_study_column <- function(df) {
    if ("Study" %in% colnames(df) && "Links" %in% colnames(df)) {
      df$Study <- paste0("<a href='", df$Links, "' target='_blank'>", df$Study, "</a>")
      df <- df[, !(names(df) %in% "Links")]
    }
    df
  }
  
  # 主表格输出（只在 CpG 模式下显示）
  output$search_result <- renderDT({
    req(input$model == "CpG")
    df <- if (input$feature == "CpG") filtered_data_cpg() else filtered_data_dmr()
    req(df)
    df <- format_study_column(df)
    DT::datatable(df,
                  filter = "top",
                  options = list(
                    pageLength = 10,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    scrollY = "444px",
                    paging = TRUE,
                    dom = 'ltp'
                  ),
                  escape = FALSE,
                  rownames = FALSE)
  })
  
  # Gene 模式下的两个表格输出
  output$search_result_cpg_gene <- renderDT({
    req(input$model == "Gene", filtered_data_cpg())
    df <- format_study_column(filtered_data_cpg())
    DT::datatable(df,
                  filter = "top",
                  options = list(
                    pageLength = 10,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    scrollY = "444px",
                    paging = TRUE,
                    dom = 'ltp'
                  ),
                  escape = FALSE,
                  rownames = FALSE)
  })
  
  output$search_result_dmr_gene <- renderDT({
    req(input$model == "Gene", filtered_data_dmr())
    df <- format_study_column(filtered_data_dmr())
    DT::datatable(df,
                  filter = "top",
                  options = list(
                    pageLength = 10,
                    autoWidth = TRUE,
                    scrollX = TRUE,
                    scrollY = "444px",
                    paging = TRUE,
                    dom = 'ltp'
                  ),
                  escape = FALSE,
                  rownames = FALSE)
  })
  
  # 下载功能（导出当前模式下合并结果）
  output$download_btn <- downloadHandler(
    filename = function() {
      paste0("eqtm_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- if (input$model == "CpG") {
        if (input$feature == "CpG") filtered_data_cpg() else filtered_data_dmr()
      } else {
        rbind(
          filtered_data_cpg() %||% data.frame(),
          filtered_data_dmr() %||% data.frame()
        )
      }
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)




