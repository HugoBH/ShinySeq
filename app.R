library(shiny)
library(ggplot2)
library(readxl)
library(bslib)

# Load run kit information from Excel
kits <- read_excel("run_kits.xlsx")
#Reads column must not have identical values.

# UI ----

ui <- fluidPage(
  
  theme = bs_theme(bootswatch = "superhero"),
  titlePanel("Welcome to ShinySeq"),
  #helpText("Calculate read depth and sample costs for Genotyping by Sequencing projects"),
  tags$p(style = "font-size: 16px; font-style: color: #ccc;",
         "Estimate sequencing depth and cost based on run and library parameters for Genotyping by Sequencing."),
  #br(),
  
  sidebarLayout(
    sidebarPanel(
      h3("Select Run Parameters"),
      # Choose run kit dynamically from Excel file
      selectInput("kit", "Sequencing kit:",
                  choices = setNames(kits$Reads, kits$KitName),
                  selected = kits$Reads[1]),
      
      # Toggle paired end
      input_switch("switch", "Paired-end"), 

      # PhiX spiking slider
      sliderInput("phix", "PhiX:",
                  min = 5, max = 40, value = 20, step = 5),
      
      # Minimum read depth
      sliderInput("min_depth", "Minimum Read Depth:",
                  min = 50, max = 1000, value = 300, step = 50),
      
      br(),
      br(),
      h3("Select Library Parameters"),
      
      # Number of species slider
      sliderInput("n_species", "Number of Species:",
                  min = 1, max = 5, value = 1, step = 1),
      
      # Dynamic UI for species
      uiOutput("species_ui"),
      
      br(),
      br(),
      
      # Toggle log scale on plot
      radioButtons(
        "y_scale", "Y-axis scale:",
        choices = c("Linear" = "linear", "Logarithmic" = "log"),
        selected = "log"
      )
    ),
    
    mainPanel(
      h3("Run Summary"),
      tableOutput("summary"),
      br(),
      
      textOutput("total_reads"),
      textOutput("available_reads"),
      textOutput("reads_required"),
      #textOutput("cost_info"),
      br(),
      
      
      fluidRow(
        column(4,
               tags$div(class = "card text-white bg-secondary mb-3",
                        tags$div(class = "card-body",
                                 tags$h5(class = "card-title", "Read depth"),
                                 textOutput("read_depth")
                        )
               )
        ),
        
        column(4,
               tags$div(class = "card text-white bg-primary mb-3",
                        tags$div(class = "card-body",
                                 tags$h5(class = "card-title", "Sample cost"),
                                 textOutput("cost_per_genotype")
                        )
               )
        ),
        
        column(4,
               tags$div(class = "card text-white bg-secondary mb-3",
                        tags$div(class = "card-body",
                                 tags$h5(class = "card-title", "Run cost"),
                                 textOutput("cost_info")
                        )
               )
        )
      ),
      br(),     
      
      
      h3("Cost Scaling"),
      plotOutput("cost_plot"),
      br(),
      
      h3("Cost Breakdown"),
      helpText("Costs are approximate"),
      tableOutput("cost_summary"),
      br(), 
      
      h3("More Information"),
      uiOutput("kit_link"),
      textOutput("email"),
      br()
      
      
    )
  )
)


# Server ----

server <- function(input, output, session) {
  
  # Dynamic sliders for species
  output$species_ui <- renderUI({
    lapply(1:input$n_species, function(i) {
      tagList(
        fluidRow(
          column(6,
                 sliderInput(paste0("samples_", i), paste("Samples (Sp", i, ")", sep = ""),
                             min = 0, max = 5000, value = 3000, step = 100)
          ),
          column(6,
                 sliderInput(paste0("loci_", i), paste("Loci (Sp", i, ")", sep = ""),
                             min = 0, max = 1000, value = 300, step = 50)
          )
        ),
        hr()
      )
    })
  })
  
  
  
  # Kit info reactive
  kit_info <- reactive({
    kits[kits$Reads == as.numeric(input$kit), ]
  })
  
  # Sequencing stats
  sequencing_data <- reactive({
    phix <- as.numeric(input$phix)
    total_reads <- as.numeric(input$kit)
    
    # Adjust for paired-end switch
    if (input$switch) {
      total_reads <- total_reads / 2
    }
    
    available_reads <- total_reads * (1 - phix/100)
    
    data.frame(
      PhiX = phix,
      TotalReads = round(total_reads / 1e6, 1),
      TotalAvailableReads = round(available_reads / 1e6, 1)
    )
  })
  
  # Species data
  species_data <- reactive({
    n <- input$n_species
    samples <- sapply(1:n, function(i) input[[paste0("samples_", i)]])
    loci <- sapply(1:n, function(i) input[[paste0("loci_", i)]])
    
    data.frame(
      Species = 1:n,
      Samples = samples,
      Loci = loci,
      ReadsPerSpecies = samples * loci,
      PropRun = (samples * loci) / sum(samples * loci)
    )
  })
  
  # Summary table
  output$summary <- renderTable({
    species_data()
  })
  
  # Stats outputs
  output$total_reads <- renderText({
    df <- sequencing_data()
    paste(df$TotalReads[1], "Million single reads in this kit")
  })
  
  output$available_reads <- renderText({
    df <- sequencing_data()
    paste(df$TotalAvailableReads[1], "Million reads available with", df$PhiX, "% spiking")
  })
  
  output$reads_required <- renderText({
    df <- species_data()
    paste(sum(df$ReadsPerSpecies), "unique reads")
  })
  
  output$read_depth <- renderText({
    df <- species_data()
    reads <- as.numeric(input$kit) * (1 - as.numeric(input$phix)/100)
    paste(round(reads / sum(df$ReadsPerSpecies), 0),
          "reads per sample per locus")
  })
  
  
  output$read_depth <- renderText({
    df <- species_data()
    
    # Total available reads
    reads <- as.numeric(input$kit) * (1 - as.numeric(input$phix) / 100)
    
    # Calculate read depth
    depth <- round(reads / sum(df$ReadsPerSpecies), 0)
    
    # Return warning if depth too low
    if (depth < input$min_depth) {
      paste0(depth, " ⚠️ (<", input$min_depth, ")")
    } else {
      paste(depth, "")
    }
  })
  
  # Show cost
  output$cost_info <- renderText({
    ki <- kit_info()
    paste0("£", ki$Cost)
  })
  
  # Add weblink
  output$kit_link <- renderUI({
    tags$a(href = "https://youtu.be/93lrosBEW-Q?t=27", 
           "Click here for inspiration", 
           target = "_blank")
  })
  
  # Add email :
  output$email <- renderText({
    paste("hugo.harrison@bristol.ac.uk")
  })
  
  
  # Cost per genotype
  output$cost_per_genotype <- renderText({
    df <- species_data()
    ki <- kit_info()
    
    total_genotypes <- sum(df$Samples)
    
    if (total_genotypes > 0) {
      cost_pg <- ki$Cost / total_genotypes
      paste0("£", round(cost_pg, 4))
    } else {
      "Invalid # samples"
    }
  })
  
  
  #Cost plot
  output$cost_plot <- renderPlot({
    
    df_species <- species_data()
    sample_range_max <- 10000
    sum_loci <- sum(df_species$Loci)
    n_species <- as.numeric(input$n_species)
    phix <- as.numeric(input$phix)
    
    line_list <- list()
    
    for (i in seq_len(nrow(kits))) {
      kit_name <- kits$KitName[i]
      kit_cost <- kits$Cost[i]
      kit_reads <- kits$Reads[i]
      
      # Adjust for paired-end switch
      if (input$switch) {
        kit_reads <- kit_reads / 2
      }
      
      available_reads <- kit_reads * (1 - phix / 100)
      
      # cutoff S_max for read depth ≥ min_depth
      s_max <- if (sum_loci > 0) {
        floor(available_reads * n_species / (input$min_depth * sum_loci))
      } else 0L
      s_max <- min(s_max, sample_range_max)
      
      # valid segment (depth ≥ min_depth)
      if (s_max >= 1) {
        samples_seq <- seq_len(s_max)
        line_list[[length(line_list) + 1]] <- data.frame(
          Samples = samples_seq,
          CostPerSample = kit_cost / samples_seq,
          Kit = kit_name,
          Segment = "valid"
        )
      }
      
      # invalid segment (depth < 100)
      if (s_max < sample_range_max) {
        invalid_samples <- seq(s_max + 1, sample_range_max)
        if (length(invalid_samples) > 0) {
          line_list[[length(line_list) + 1]] <- data.frame(
            Samples = invalid_samples,
            CostPerSample = kit_cost / invalid_samples,
            Kit = kit_name,
            Segment = "invalid"
          )
        }
      }
    }
    
    plot_df <- do.call(rbind, line_list)
    
    # Build plot
    p <- ggplot(plot_df, aes(x = Samples, y = CostPerSample,
                             group = interaction(Kit, Segment))) +
      geom_line(aes(colour = Kit), linetype = "solid",
                data = subset(plot_df, Segment == "valid"), linewidth = 1) +
      #scale_linetype_manual(values = c(rep("solid", 5), rep("dotted", 3))) +
      geom_line(aes(colour = Kit), linetype = "dotted",
                data = subset(plot_df, Segment == "invalid"), linewidth = .8) +
      scale_colour_manual(values = c("#a63603", "#e6550d", "#fd8d3c", "#fdae6b", 
                                     "#a6cee3", "#1b7837", "#5aae61", "#a6dba0"))
    
    # Chosen parameters (red marker)
    chosen_samples <- sum(df_species$Samples)
    ki <- kit_info()
    chosen_cost_ps <- if (chosen_samples > 0) ki$Cost / chosen_samples else NA_real_
    available_reads_chosen <- as.numeric(ki$Reads) * (1 - phix / 100)
    s_max_chosen <- if (sum_loci > 0) floor(available_reads_chosen * n_species / (input$min_depth * sum_loci)) else 0
    x_point <- pmin(chosen_samples, sample_range_max)
    
    if (is.finite(chosen_cost_ps) && chosen_samples >= 1) {
      if (chosen_samples <= s_max_chosen) {
        p <- p + geom_point(aes(x = x_point, y = chosen_cost_ps),
                            colour = "#d9534f", size = 3)
      } else {
        p <- p + geom_point(aes(x = x_point, y = chosen_cost_ps),
                            colour = "#d9534f", size = 3, shape = 4) +
          annotate("text", x = x_point, y = chosen_cost_ps,
                   label = "⚠️ Low Read Depth", colour = "#d9534f", vjust = -1)
      }
    }
    
    # Labels and theme
    p <- p +
      labs(x = "Number of Samples", y = "Cost per Sample") +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(), 
            legend.position = c(.7,.7),
            plot.background = element_rect(fill = "#dcdfe2", color = NA),
            panel.background = element_rect(fill = "#dcdfe2", color = NA),
            panel.grid.major = element_line(color = "#ebebeb"),
            #panel.grid.minor = element_line(color = "#ebebeb"),
            axis.text = element_text(color = "#20374c"),
            axis.title = element_text(color = "#20374c")
      )
    
    # Axis scale
    if (input$y_scale == "log") {
      p <- p + scale_y_log10()
    } else {
      p <- p + scale_y_continuous(limits = c(0,1))
    }
    
    print(p)
  })
  
  
  # Cost summary table
  output$cost_summary.temp <- renderTable({
    df <- species_data()
    total_genotypes <- sum(df$Samples)
    ki <- kit_info()
    
    #per sample costs
    DNA_extraction = 0.25 
    PCR1 = 0.25
    PCR2 = 0.34
    QC =  121.00 
    TRAC = 654
    
    data.frame(
      Item = c("Extraction", "PCR1", "PCR2", "Sequencing", "QC", "Total"),
      PerSampleCost = c(paste("£", DNA_extraction, sep = ""),
                   paste("£", PCR1, sep = ""),
                   paste("£", PCR2, sep = ""),
                   paste("£", round(ki$Cost / total_genotypes, 2), sep = ""),
                   paste("£", QC, sep = ""), 
                   paste("£", round(sum(DNA_extraction, PCR1, PCR2) + 
                                      ki$Cost / total_genotypes + QC / total_genotypes, 2))),
      TotalCost = c(paste("£", DNA_extraction * total_genotypes, sep = ""),
                  paste("£", PCR1 * total_genotypes, sep = ""),
                  paste("£", PCR2 * total_genotypes, sep = ""),
                  paste("£", ki$Cost, sep = ""),
                  paste("£", QC, sep = ""), 
                  paste("£", round(sum(DNA_extraction, PCR1, PCR2) * total_genotypes + ki$Cost + QC, 2)))
    )
    
  })
  
  
  
  output$cost_summary <- renderTable({
    df <- species_data()
    total_genotypes <- sum(df$Samples)
    ki <- kit_info()
    
    # Define per-sample costs
    day_rate <- 258.25
    consumables <- c(Extraction = 0.25, PCR1 = 0.25, PCR2 = 0.34)
    trac <- c(Extraction = day_rate / 400,   #how many extractions one can do in a day
              PCR1 = day_rate / 1200,        #how many PCRs (3 x 4 plates) one can do in a day
              PCR2 = day_rate / 1200)       #how many PCRs (3 x 4 plates) one can do in a day
    sequencing <- ki$Cost / total_genotypes
    seq_trac <- day_rate / 2
    qc <- 121.00
    
    # Total cost calculation
    total_per_sample <- sum(consumables) + sequencing + qc / total_genotypes
    total_trac <- sum(trac) * total_genotypes + seq_trac + qc
    total_run_cost <- sum(consumables, trac) * total_genotypes + ki$Cost + seq_trac + qc
    
    # Create data frame
    data.frame(
      Item = c(names(consumables), "Sequencing", "QC", "Total"),
      Consumables = paste0("£", round(c(consumables, sequencing, qc / total_genotypes, total_per_sample), 2)),
      TRAC = paste0("£", round(c(trac * total_genotypes, seq_trac, qc, total_trac))), 
      TotalCost = paste0("£", round(c((consumables + trac) * total_genotypes, ki$Cost + seq_trac, qc, total_run_cost), 2))
    )
  })
  
  
  
  
  
}

shinyApp(ui = ui, server = server)