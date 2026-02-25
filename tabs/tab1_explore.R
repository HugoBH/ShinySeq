

# Load run kit information from Excel
kits <- read_excel("run_kits.xlsx") %>% 
  mutate(
    KitName   = str_squish(str_replace_all(KitName, "\u00A0", " ")),
    Instrument = str_squish(str_replace_all(Instrument, "\u00A0", " "))
  )

#Reads column must not have identical values.


tab_one_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Plan Your Run",
    
    sidebarLayout(
      sidebarPanel(
        
        h3("Select Run Parameters"),
        
        # Checkbox group input
        checkboxGroupInput(
          inputId = ns("Instrument"),
          label = "Select Instrument",
          choices = unique(kits$Instrument),
          selected = intersect("NextSeq 2000", unique(kits$Instrument))
        ),
        
        
        selectInput(ns("kit"), "Sequencing kit:",
                    choices = NULL,
                    selected = NULL),
        
        #bslib::input_switch(ns("switch"), "Paired-end"), 
        
        sliderInput(ns("phix"), "PhiX:",
                    min = 5, max = 40, value = 20, step = 5),
        
        sliderInput(ns("min_depth"), "Minimum Read Depth:",
                    min = 50, max = 1000, value = 500, step = 50),
        
        br(), br(),
        h3("Select Library Parameters"),
        
        sliderInput(ns("n_species"), "Number of Species:",
                    min = 1, max = 5, value = 1, step = 1),
        
        uiOutput(ns("species_ui")),
        
        br(), br(),
        
        radioButtons(
          ns("y_scale"), "Y-axis scale:",
          choices = c("Linear" = "linear", "Logarithmic" = "log"),
          selected = "log"
        )
      ),
      
      mainPanel(
        h3("Run Summary"),
        tableOutput(ns("summary")),
        br(),
        
        textOutput(ns("total_reads")),
        textOutput(ns("available_reads")),
        textOutput(ns("reads_required")),
        helpText("Note: Expect 80-90% On target reads"),
        br(), br(),
        
        fluidRow(
          column(4,
                 tags$div(class = "card text-white bg-secondary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Read depth"),
                                   textOutput(ns("read_depth"))
                          )
                 )
          ),
          column(4,
                 tags$div(class = "card text-white bg-primary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Sample cost"),
                                   textOutput(ns("cost_per_genotype"))
                          )
                 )
          ),
          column(4,
                 tags$div(class = "card text-white bg-secondary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Run cost"),
                                   textOutput(ns("cost_info"))
                          )
                 )
          )
        ),
        
        br(),
        h3("Cost Scaling"),
        plotOutput(ns("cost_plot")),
        br(),
        
        h3("Cost Breakdown"),
        helpText("Approximate pricing are purely indicative and vary by project"),
        
        fluidRow(
          column(8, 
                 tableOutput(ns("cost_summary"))
          ),
          column(4,
                 tags$div(class = "card text-white bg-primary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Per sample"),
                                   textOutput(ns("per_sample_cost"))
                          )
                 ),
                 tags$div(class = "card text-white bg-secondary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "Consumables only"),
                                   textOutput(ns("total_consumable_cost"))
                          )
                 ),
                 tags$div(class = "card text-white bg-primary mb-3",
                          tags$div(class = "card-body",
                                   tags$h5(class = "card-title", "TRAC costed project"),
                                   textOutput(ns("total_project_cost"))
                          )
                 )
          )
        ),
        
        br(),
        h3("More Information"),
        uiOutput(ns("kit_link")),
        textOutput(ns("email")),
        br()
      )
    )
  )
}




#Server -----
tab_one_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    output$species_ui <- renderUI({
      lapply(1:input$n_species, function(i) {
        tagList(
          fluidRow(
            column(6,
                   sliderInput(session$ns(paste0("samples_", i)), paste("Samples (Sp", i, ")"),
                               min = 0, max = 5000, value = 2000, step = 100)
            ),
            column(6,
                   sliderInput(session$ns(paste0("loci_", i)), paste("Loci (Sp", i, ")"),
                               min = 0, max = 500, value = 250, step = 25)
            )
          ),
          hr()
        )
      })
    })
    
    
    
    
    # Reactive filtered kits
    filtered_kits <- reactive({
      req(input$Instrument)                    # now works because Instrument is namespaced
      subset(kits, Instrument %in% input$Instrument)
    })
    
    # Populate/update kit dropdown whenever instruments change (including at init)
    observeEvent(input$Instrument, ignoreInit = FALSE, {
      df <- filtered_kits()
      
      # Build named choices: labels = KitName, values = KitName (consistent!)
      choices <- setNames(df$KitName, df$KitName)
      
      # Your desired default: 3rd row of the full kits table
      default_kit <- if (nrow(kits) >= 3) kits$KitName[3] else NULL
      
      # Try to keep the current selection (now also a KitName); else use default if present in filtered df; else first
      current <- isolate(input$kit)
      selected <-
        if (!is.null(current) && current %in% df$KitName) {
          current
        } else if (!is.null(default_kit) && default_kit %in% df$KitName) {
          default_kit
        } else if (nrow(df)) {
          df$KitName[3]
        } else {
          NULL
        }
      
      updateSelectInput(session, "kit",
                        choices = choices,
                        selected = selected
      )
    })
    
    # Require a kit before using it
    kit_info <- reactive({
      req(input$kit)
      kits[kits$KitName == input$kit, , drop = FALSE]
    })
  
    
    kit_info <- reactive({
      kits[kits$KitName == input$kit, ]
    })
    
    sequencing_data <- reactive({
      phix <- as.numeric(input$phix)
      total_reads <- kit_info()$Reads
      #if (input$switch) total_reads <- total_reads * 2
      available_reads <- total_reads * (1 - phix / 100)
      
      data.frame(
        PhiX = phix,
        TotalReads = round(total_reads / 1e6, 1),
        TotalAvailableReads = round(available_reads / 1e6, 1)
      )
    })
    
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
    
    output$summary <- renderTable({
      species_data()
    })
    
    output$total_reads <- renderText({
      df <- sequencing_data()
      paste(df$TotalReads[1], "Million single-end reads in this kit")
    })
    
    output$available_reads <- renderText({
      df <- sequencing_data()
      paste(df$TotalAvailableReads[1], "Million single-end reads available with", df$PhiX, "% spiking")
    })
    
    output$reads_required <- renderText({
      df <- species_data()
      paste(sum(df$ReadsPerSpecies), "unique reads")
    })
    
    output$read_depth <- renderText({
      df <- species_data()
      seq <- sequencing_data()
      depth <- round((seq$TotalAvailableReads[1] * 1e6) / sum(df$ReadsPerSpecies), 0)
      if (depth < input$min_depth) {
        paste0(depth, " ⚠️ (<", input$min_depth, ")")
      } else {
        paste(depth)
      }
    })
    
    output$cost_info <- renderText({
      ki <- kit_info()
      paste0("£", ki$Cost)
    })
    
    output$kit_link <- renderUI({
      tags$a(href = "https://youtu.be/93lrosBEW-Q?t=27", 
             "Click here for inspiration", 
             target = "_blank")
    })
    
    output$email <- renderText({
      "hugo.harrison@bristol.ac.uk"
    })
    
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
    
    
    #Cost plot ----
    output$cost_plot <- renderPlot({
      
      # Use the filtered kits based on selected Instruments
      df_kits <- filtered_kits()
      
      # If no kits available for the selected instruments, draw an empty message plot
      if (nrow(df_kits) == 0) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No kits available for the selected instrument(s)",
                   size = 5, colour = "#20374c") +
          theme_void() +
          theme(plot.background = element_rect(fill = "#dcdfe2", color = NA)) |> print()
        return(invisible(NULL))
      }
      
      # Species and sequencing params (unchanged)
      df_species <- species_data()
      sample_range_max <- 10000
      sum_loci <- sum(df_species$Loci)
      n_species <- as.numeric(input$n_species)
      phix <- as.numeric(input$phix)
      
      line_list <- list()
      
      # Build lines ONLY from the filtered kits
      for (i in seq_len(nrow(df_kits))) {
        kit_name <- df_kits$KitName[i]
        kit_cost <- as.numeric(df_kits$Cost[i])
        kit_reads <- as.numeric(df_kits$Reads[i])
        
        # Adjust for paired-end switch if you ever re-enable it
        # if (input$switch) { kit_reads <- kit_reads / 2 }
        
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
        
        # invalid segment (depth < min_depth)
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
      
      # If nothing to plot (e.g., all lines empty), show message
      if (length(line_list) == 0) {
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = "No plot available for current parameters",
                   size = 5, colour = "#20374c") +
          theme_void() +
          theme(plot.background = element_rect(fill = "#dcdfe2", color = NA)) |> print()
        return(invisible(NULL))
      }
      
      plot_df <- do.call(rbind, line_list)
      
      # Build plot (unchanged aesthetics)
      p <- ggplot(plot_df, aes(x = Samples, y = CostPerSample,
                               group = interaction(Kit, Segment))) +
        geom_line(aes(colour = Kit), linetype = "solid",
                  data = subset(plot_df, Segment == "valid"), linewidth = 1) +
        geom_line(aes(colour = Kit), linetype = "dotted",
                  data = subset(plot_df, Segment == "invalid"), linewidth = .8) +
        scale_color_viridis_d()
      
      # Chosen parameters (red marker)
      chosen_samples <- sum(df_species$Samples)
      ki <- kit_info()  # <- uses selected KitName → row lookup
      chosen_cost_ps <- if (chosen_samples > 0) as.numeric(ki$Cost) / chosen_samples else NA_real_
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
      
      # Labels, theme, scales (same as before)
      p <- p +
        labs(x = "Number of Samples", y = "Cost per Sample") +
        theme_minimal(base_size = 14) +
        theme(panel.grid = element_blank(),
              legend.position = c(.7,.7),
              plot.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.background = element_rect(fill = "#dcdfe2", color = NA),
              panel.grid.major = element_line(color = "#ebebeb"),
              axis.text = element_text(color = "#20374c"),
              axis.title = element_text(color = "#20374c"))
      
      if (input$y_scale == "log") {
        p <- p + scale_y_log10()
      } else {
        p <- p + scale_y_continuous(limits = c(0,1))
      }
      
      print(p)
    })
    
    
    # Cost summary table
    
    day_rate <- 258.25 #wetlab
    bioinformatics <- 357.90 #day rate bioinformatic
    
    consumables <- c(Extraction = 0.25, PCR1 = 0.25, PCR2 = 0.34)
    trac <- c(Extraction = day_rate / 400,   #how many extractions one can do in a day
              PCR1 = day_rate / 1200,        #how many PCRs (3 x 4 plates) one can do in a day
              PCR2 = day_rate / 1200)       #how many PCRs (3 x 4 plates) one can do in a day
    seq_trac <- day_rate / 2
    qc <- 121.00
    
    
    output$cost_summary <- renderTable({
      df <- species_data()
      total_genotypes <- sum(df$Samples)
      ki <- kit_info()
      
      # Define per-sample sequencing costs
      sequencing <- ki$Cost / total_genotypes
      
      # Total cost calculation
      total_per_sample <- sum(consumables) + sequencing
      total_trac <- sum(trac) * total_genotypes + seq_trac + qc + bioinformatics / 2
      total_run_cost <- sum(consumables, trac) * total_genotypes + ki$Cost + seq_trac + qc + bioinformatics / 2
      
      # Create data frame
      data.frame(
        Item = c(names(consumables), "Sequencing", "QC", "Bioinformatics", "Total"),
        Consumables = paste0("£", round(c(consumables, sequencing, 0, 0, total_per_sample), 2)),
        TRAC = paste0("£", round(c(trac * total_genotypes, seq_trac, qc, bioinformatics / 2, total_trac))), 
        TotalCost = paste0("£", round(c((consumables + trac) * total_genotypes, ki$Cost + seq_trac, qc, bioinformatics / 2, total_run_cost)))
      )
    })
    
    
    # Per-sample costs
    output$per_sample_cost <- renderText({
      df <- species_data()
      total_genotypes <- sum(df$Samples)
      ki <- kit_info()
      
      consumable_cost <- sum(consumables) + ki$Cost / total_genotypes
      paste0("£", round(consumable_cost, 2))
    })  
    
    # Total consumables
    output$total_consumable_cost <- renderText({
      df <- species_data()
      total_genotypes <- sum(df$Samples)
      ki <- kit_info()
      
      consumable_cost <- sum(consumables) * total_genotypes + ki$Cost 
      paste0("£", round(consumable_cost, 0))
    })  
    
    output$total_project_cost <- renderText({
      df <- species_data()
      total_genotypes <- sum(df$Samples)
      ki <- kit_info()
      
      total_cost <- sum(consumables, trac) * total_genotypes + ki$Cost + seq_trac + qc + bioinformatics / 2
      paste0("£", round(total_cost, 0))
    })
    
    # Add other outputs as needed...
  })
}

