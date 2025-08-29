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
                  min = 10, max = 40, value = 20, step = 5),
      
      # Number of species slider
      sliderInput("n_species", "Number of Species:",
                  min = 1, max = 10, value = 2, step = 1),
      
      # Dynamic UI for species names, samples, and loci
      uiOutput("species_ui")
    ),
    
    mainPanel(
      h3("Species Summary"),
      tableOutput("summary"),
      h3("Stats Summary"),
      textOutput("total_reads"),
      textOutput("available_reads"),
      textOutput("reads_required"),
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
        fluidRow(
          column(6,
                 sliderInput(paste0("samples_", i), paste("Samples (Sp", i, ")", sep = ""),
                             min = 0, max = 5000, value = 100, step = 100)
          ),
          column(6,
                 sliderInput(paste0("loci_", i), paste("Loci (Sp", i, ")", sep = ""),
                             min = 0, max = 1000, value = 50, step = 50)
          )
        ),
        hr() # horizontal line between species for clarity
      )
    })
  })
  
  # Generate sequencing data frame
  sequencing_data <- reactive({
    phix <- as.numeric(input$phix)
    total_reads <- as.numeric(input$reads)
    available_reads <- as.numeric(input$reads) * (1 - as.numeric(input$phix)/100)
    
    data.frame(
      PhiX = phix,
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
      Species = 1:input$n_species,
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
  output$total_reads <- renderText({
    df <- sequencing_data()
    if (!is.na(df$TotalReads[1])) {
      paste(df$TotalReads[1], "Million reads in this kit")
    } else {
      "Please enter valid numbers for PhiX"
    }
  })
    
  # Overall available reads
  output$available_reads <- renderText({
    df <- sequencing_data()
    if (!is.na(df$TotalAvailableReads[1])) {
      paste(df$TotalAvailableReads[1], "Million reads available with", df$PhiX, "% spiking")
    } else {
      "Please enter valid numbers for PhiX"
    }
  })
  
  # Overall total reads 
  output$reads_required <- renderText({
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
