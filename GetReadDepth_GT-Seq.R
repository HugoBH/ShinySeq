ui <- fluidPage(
  titlePanel("Sequencing Read Depth Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      # Choose available reads
      selectInput("reads", "Select run kit:",
                  choices = c("NextSeq500 1x75 High Output" = 400e6, #25-30Gb
                              "NextSeq500 2x75 Mid Output" = 260e6,  #50-60Gb
                              "NextSeq500 2x75 High Output" = 800e6),#50-60Gb
                  selected = 400e6),
      
      # Spiking slider
      sliderInput("phix", "PhiX:",
                  min = 10, max = 60, value = 30, step = 10),
      
      # Number of species slider
      sliderInput("n_species", "Number of Species:",
                  min = 1, max = 10, value = 1, step = 1),
      
      # Dynamic UI for species names, samples, and loci
      uiOutput("species_ui")
    ),
    
    mainPanel(
      h3("Species Summary"),
      tableOutput("summary"),
      #h3("Summary Stats"),
      textOutput("available_reads"),
      textOutput("total_reads"),
      h3("Read Depth"),
      textOutput("read_depth")
    )
  )
)

server <- function(input, output, session) {
  
  # Dynamic UI for each species: name, samples, loci
  output$species_ui <- renderUI({
    lapply(1:input$n_species, function(i) {
      tagList(
        textInput(paste0("species_name_", i), paste("Name of Species", i),
                   value = paste("Species", i)),
        fluidRow(
          column(6,
                 numericInput(paste0("samples_", i), 
                              paste("Samples for Species", i), 
                              value = 1, min = 1)
          ),
          column(6,
                 numericInput(paste0("loci_", i), 
                              paste("Loci for Species", i), 
                              value = 1, min = 1)
          )
        ),
        hr() # horizontal line between species for clarity
      )
    })
  })
  
  # Generate sequencing data frame
  sequencing_data <- reactive({
    total_reads <- as.numeric(input$reads)
    available_reads <- as.numeric(input$reads) * (1 - as.numeric(input$phix)/100)
    
    data.frame(
      TotalReads = round(total_reads / 1e6, 1),
      TotalAvailableReads = round(available_reads  / 1e6, 1)
    )
  })
  
  # Generate species data frame
  species_data <- reactive({
    n <- input$n_species
    names <- sapply(1:n, function(i) input[[paste0("species_name_", i)]])
    samples <- sapply(1:n, function(i) input[[paste0("samples_", i)]])
    loci <- sapply(1:n, function(i) input[[paste0("loci_", i)]])

    data.frame(
      Species = names,
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
  
  # Overall available reads
  output$available_reads <- renderText({
    df <- sequencing_data()
    if (!is.na(df$TotalAvailableReads[1])) {
      paste(df$TotalAvailableReads[1], "Million available reads")
    } else {
      "Please enter valid numbers for PhiX"
    }
  })
  
  # Overall total reads 
  output$total_reads <- renderText({
    df <- species_data()
    if (!is.na(sum(df$ReadsPerSpecies))) {
      paste(sum(df$ReadsPerSpecies), "reads required")
    } else {
      "Please enter valid numbers for samples and loci."
    }
  })
  
  # Overall read depths
  output$read_depth <- renderText({
    df <- species_data()
    if (!is.na(sum(df$ReadsPerSpecies))) {
      paste(round( (as.numeric(input$reads) * (1 - as.numeric(input$phix)/100)) / sum(df$ReadsPerSpecies) , 0), 
            "reads per sample per locus")
    } else {
      "Please enter valid numbers for samples and loci."
    }
  })
}

shinyApp(ui = ui, server = server)
