# Load tidyverse
library(tidyverse)

tab_two_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Result Summary",
    
    sidebarLayout(
      sidebarPanel(
        h3("Upload Your Data"),
        textOutput(ns("data_source")),
        fileInput(ns("csv_file"), "Choose CSV File", accept = ".csv"),
        checkboxInput(ns("header"), "Header", TRUE),
        radioButtons(ns("sep"), "Separator",
                     choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                     selected = ","),
        numericInput(ns("threshold"), "Read Depth Threshold", value = 10, min = 0),
        br(),
        sliderInput(ns("success"), "Genotyping success threshold:",
                    min = .60, max = 1.00, value = .80, step = .05),
        
      ),
      
      mainPanel(
        fluidRow(
          column(3,
                 tags$div(class = "card text-white bg-secondary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Samples"),
                                   textOutput(ns("sample_number"))
                          )
                 )
          ),
          column(3,
                 tags$div(class = "card text-white bg-primary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Loci"),
                                   textOutput(ns("locus_number"))
                          )
                 )
          ),
          column(3,
                 tags$div(class = "card text-white bg-secondary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Success Rate"),
                                   textOutput(ns("sample_success"))
                          )
                 )
          ),
          column(3,
                 tags$div(class = "card text-white bg-primary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Read Depth"),
                                   textOutput(ns("overall_mean_read_depth"))
                          )
                 )
          )
        ),
        
        helpText("Samples containing `negative` are automatically removed"),
        br(),
        
        
        h3("Preview"),
        tableOutput(ns("preview")),
        br(),
        
        h3("Genotyping Success"),
        fluidRow(
          column(6, plotOutput(ns("sample_success_plot"))),
          column(6, plotOutput(ns("locus_success_plot")))
        ),
        
        fluidRow(
          column(6,tags$div(class = "card text-white bg-primary mb-3",
                            tags$div(class = "card-body",
                                     textOutput(ns("sample_success_proportion"))))),
          column(6, tags$div(class = "card text-white bg-primary mb-3",
                             tags$div(class = "card-body",
                                      textOutput(ns("locus_success_proportion")))))
        ),
        br(),
        
        h3("Read Depth Distribution"),
        fluidRow(
          column(6, plotOutput(ns("sample_readdepth_plot"))),
          column(6, plotOutput(ns("locus_readdepth_plot")))
        ),
        
        fluidRow(
          column(6,tags$div(class = "card text-white bg-primary mb-3",
                            tags$div(class = "card-body",
                                     textOutput(ns("sample_read_depth_summary"))))),
          column(6, tags$div(class = "card text-white bg-primary mb-3",
                             tags$div(class = "card-body",
                                      textOutput(ns("locus_read_depth_summary")))))
        ),
        br(),
        
        
        h3("Highly expressed samples and loci"),
        fluidRow(
          column(12, tableOutput(ns("sample_stats_max"))),
          column(12, tableOutput(ns("locus_stats_max")))
        ), 
        
        br(),
        
        h3("Lowly expressed samples and loci"),
        fluidRow(
          column(12, tableOutput(ns("sample_stats_min"))),
          column(12, tableOutput(ns("locus_stats_min")))
        ),
        br()
        
      )
    )
  )
}





tab_two_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    output$data_source <- renderText({
      if (is.null(input$csv_file)) {
        "Using default dataset"
      } else {
        paste("Using uploaded file:", input$csv_file$name)
      }
    })
    
    # Reactive: read uploaded CSV
    data <- reactive({
      
      # If user uploaded a file → use it
      if (!is.null(input$csv_file)) {
        
        df <- read.csv(
          input$csv_file$datapath,
          header = input$header,
          sep = input$sep
        )
        
      } else {
        
        # Otherwise load default CSV
        df <- read.csv(
          "data/Library_Genotypes_Counts_GTseq_Run4.csv",
          header = TRUE,
          sep = ","
        )
        
      }
      
      # Standardise column names safely
      df <- df %>%
        rename(
          On_Target = any_of("X.On.Target"),
          On_Target_Reads = any_of("On.Target.Reads"),
          GT = any_of("X.GT"),
          Raw_Reads = any_of("Raw.Reads")
        ) %>%
        filter(!str_detect(Sample, regex("negative", ignore_case = TRUE)))
      
      df
    })
    
    
    # Reactive: identify locus columns
    locus_cols <- reactive({
      df <- data() %>% select(starts_with("X"))
      # Return locus columns (those starting with 'X')
      grep("^X", names(df), value = TRUE)
    })
    
    
    # Reactive: calculate sample genotyping success
    sample_success <- reactive({
      df <- data()
      loci <- locus_cols()
      threshold <- input$threshold
      df$SuccessRate <- rowSums(df[loci] >= threshold) / ncol(df[loci])
      df[, c("Sample", "SuccessRate")]
    })
    
    # Reactive: calculate locus genotyping success
    locus_success <- reactive({
      df <- data()
      loci <- locus_cols()
      threshold <- input$threshold
      success <- colSums(df[loci] >= threshold) / nrow(df[loci])
      data.frame(Locus = loci, SuccessRate = success)
    })
    
    # Reactive: calculate mean and variance per sample
    sample_stats <- reactive({
      df <- data()
      loci <- locus_cols()
      data.frame(
        Sample = df$Sample,
        Sum_Reads = rowSums(df[loci]),
        Mean_ReadDepth = rowMeans(df[loci]),
        Variance_ReadDepth = apply(df[loci], 1, var),
        SD_ReadDepth = apply(df[loci], 1, sd),
        SE_ReadDepth = apply(df[loci], 1, sd) / sqrt(length(loci))
      )
    })
    
    
    locus_stats <- reactive({
      df <- data()
      loci <- locus_cols()
      data.frame(
        Locus = colnames(df[loci]),
        Sum_Reads = apply(df[loci], 2, sum),
        Mean_ReadDepth = apply(df[loci], 2, mean),
        Variance_ReadDepth = apply(df[loci], 2, var),
        SD_ReadDepth = apply(df[loci], 2, sd),
        SE_ReadDepth = apply(df[loci], 2, sd) / sqrt(length(loci))
      )
    })
    
    
    #########
    # Outputs
    #########
    
    #Sample number
    output$sample_number <- renderText({
      nrow(data())
    })
    
    #Locus number
    output$locus_number <- renderText({
      length(locus_cols())
    })
    
    #Sample genotyping success (Same when calculated by locus)
    output$sample_success <- renderText({
      df <- sample_success()
      paste(round(mean(df$SuccessRate, na.rm = TRUE)*100,1), "%")
    })
    
    #Read depth
    output$overall_mean_read_depth <- renderText({
      df <- data()
      loci <- locus_cols()
      paste(round(mean(unlist(df[loci])), 0), "x")
    })
    
    #Preview
    output$preview <- renderTable({
      head(data(), 3)[,1:10]
    })
    
    #Plots
    output$sample_success_plot <- renderPlot({
      df <- sample_success()
      ggplot(df, aes(x = SuccessRate)) +
        geom_histogram() +
        geom_vline(xintercept = input$success, col = "#d9534f", linetype = "dashed") +
        theme_minimal(base_size = 14) +
        labs(y = "Number of Samples", 
             x = "Sample Genotyping Success") +
        theme(panel.grid = element_blank(), 
              legend.position = "none",
              plot.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.grid.major = element_line(color = "#ebebeb"),
              #panel.grid.minor = element_line(color = "#ebebeb"),
              axis.text = element_text(color = "#20374c"),
              axis.title = element_text(color = "#20374c")
        )
    })
    
    output$locus_success_plot <- renderPlot({
      df <- locus_success()
      ggplot(df, aes(x = SuccessRate)) +
        geom_histogram() +
        geom_vline(xintercept = input$success, col = "#d9534f", linetype = "dashed") +
        theme_minimal(base_size = 14) +
        labs(y = "Number of loci", 
             x = "Locus Genotyping Success") +
        theme(panel.grid = element_blank(), 
              legend.position = "none",
              plot.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.grid.major = element_line(color = "#ebebeb"),
              #panel.grid.minor = element_line(color = "#ebebeb"),
              axis.text = element_text(color = "#20374c"),
              axis.title = element_text(color = "#20374c")
        )
    })
    
    
    #Summary numbers
    output$sample_success_proportion <- renderText({
      df <- sample_success()
      paste(sum(df$SuccessRate >= input$success), "samples with >", input$success*100, "% genotype")
    })
    
    output$locus_success_proportion <- renderText({
      df <- locus_success()
      paste(sum(df$SuccessRate >= input$success), "loci with >", input$success*100, "% genotype")
    })
    
    
    # Distribution plots
    
    output$sample_readdepth_plot <- renderPlot({
      stats <- sample_stats()
      
      ggplot(stats, aes(x = reorder(Sample, Mean_ReadDepth), y = Mean_ReadDepth)) +
        geom_errorbar(aes(ymin=Mean_ReadDepth-SE_ReadDepth, 
                          ymax=Mean_ReadDepth+SE_ReadDepth), width=0) +
        geom_point(col = "#d9534f") +
        theme_minimal(base_size = 14) +
        labs(title = "Read Depth by Sample",
             x = "Sample", y = expression(Mean~Read~Depth %+-% SE)) +
        theme(axis.text.x = element_blank(),
              panel.grid = element_blank(), 
              legend.position = "none",
              plot.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.background = element_rect(fill = "#dcdfe2", color = NA),
              #panel.border = element_rect(fill = "transparent", colour = "#20374c"),
              panel.grid.major.y = element_line(color = "#ebebeb"),
              #panel.grid.minor = element_line(color = "#ebebeb"),
              axis.text = element_text(color = "#20374c"),
              axis.title = element_text(color = "#20374c")
        )
      
    })
    
    
    output$locus_readdepth_plot <- renderPlot({
      stats <- locus_stats()
      
      ggplot(stats, aes(x = reorder(Locus, Mean_ReadDepth), y = Mean_ReadDepth)) +
        geom_errorbar(aes(ymin=Mean_ReadDepth-SE_ReadDepth, 
                          ymax=Mean_ReadDepth+SE_ReadDepth), width=0) +
        geom_point(col = "#d9534f") +
        scale_y_log10() +
        theme_minimal(base_size = 14) +
        labs(title = "Read Depth by Locus",
             x = "Locus", y = expression(Mean~Read~Depth %+-% SE)) +
        theme(axis.text.x = element_blank(),
              panel.grid = element_blank(), 
              legend.position = "none",
              plot.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.background = element_rect(fill = "#dcdfe2", color = NA),
              #panel.border = element_rect(fill = "transparent", colour = "#20374c"),
              panel.grid.major.y = element_line(color = "#ebebeb"),
              #panel.grid.minor = element_line(color = "#ebebeb"),
              axis.text = element_text(color = "#20374c"),
              axis.title = element_text(color = "#20374c")
        )
      
    })
    
    #Summary numbers
    output$sample_read_depth_summary <- renderText({
      df <- sample_stats()
      paste(sum(df$Mean_ReadDepth >= 100), "samples > 100x")
    })
    
    output$locus_read_depth_summary <- renderText({
      df <- locus_stats()
      paste(sum(df$Mean_ReadDepth >= 100), "loci > 100x")
    })    
    
    #High and low loci
    output$locus_stats_max <- renderTable({
      locus_stats() %>% 
        select(Locus, Mean_ReadDepth, SE_ReadDepth, Sum_Reads) %>% 
        mutate(PercentOflReads = Sum_Reads / sum(Sum_Reads) * 100) %>%
        slice_max(PercentOflReads, n = 10)
    })
    output$locus_stats_min <- renderTable({
      locus_stats() %>% 
        select(Locus, Mean_ReadDepth, SE_ReadDepth, Sum_Reads) %>% 
        mutate(PercentOfReads = Sum_Reads / sum(Sum_Reads) * 100) %>%
        slice_min(PercentOfReads, n = 10)
    })
    
    #High and low samples
    output$sample_stats_max <- renderTable({
      sample_stats() %>% 
        select(Sample, Mean_ReadDepth, SE_ReadDepth, Sum_Reads) %>% 
        mutate(PercentOflReads = Sum_Reads / sum(Sum_Reads) * 100) %>%
        slice_max(PercentOflReads, n = 10)
    })
    output$sample_stats_min <- renderTable({
      sample_stats() %>% 
        select(Sample, Mean_ReadDepth, SE_ReadDepth, Sum_Reads) %>% 
        mutate(PercentOfReads = Sum_Reads / sum(Sum_Reads) * 100) %>%
        slice_min(PercentOfReads, n = 10)
    })
    
    
    
  })
}
