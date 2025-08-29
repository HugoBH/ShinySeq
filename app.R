library(shiny)
library(readxl)

# Load run kit information from Excel
kits <- read_excel("run_kits.xlsx")
#Reads column must not have identical values.

ui <- fluidPage(
  titlePanel("Sequencing Read Depth Calculator"),
  
  sidebarLayout(
    sidebarPanel(
      # Choose run kit dynamically from Excel file
      selectInput("kit", "Select run kit:",
                  choices = setNames(kits$Reads, kits$KitName),
                  selected = kits$Reads[1]),
      
      # PhiX spiking slider
      sliderInput("phix", "PhiX:",
                  min = 10, max = 40, value = 20, step = 5),
      
      # Number of species slider
      sliderInput("n_species", "Number of Species:",
                  min = 1, max = 10, value = 2, step = 1),
      
      # Dynamic UI for species
      uiOutput("species_ui")
    ),
    
    mainPanel(
      h3("Species Summary"),
      tableOutput("summary"),
      
      h3("Stats Summary"),
      textOutput("total_reads"),
      textOutput("available_reads"),
      textOutput("reads_required"),
      textOutput("cost_info"),
      textOutput("cost_per_genotype"),  
      
      h3("Read Depth"),
      textOutput("read_depth"),
      
      h3("More Information"),
      uiOutput("kit_link")
    )
  )
)

server <- function(input, output, session) {
  
  # Dynamic sliders for species
  output$species_ui <- renderUI({
    lapply(1:input$n_species, function(i) {
      tagList(
        fluidRow(
          column(6,
                 sliderInput(paste0("samples_", i), paste("Samples (Sp", i, ")"),
                             min = 0, max = 5000, value = 100, step = 100)
          ),
          column(6,
                 sliderInput(paste0("loci_", i), paste("Loci (Sp", i, ")"),
                             min = 0, max = 1000, value = 50, step = 50)
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
    paste(df$TotalReads[1], "Million reads in this kit")
  })
  
  output$available_reads <- renderText({
    df <- sequencing_data()
    paste(df$TotalAvailableReads[1], "Million reads available with", df$PhiX, "% spiking")
  })
  
  output$reads_required <- renderText({
    df <- species_data()
    paste(sum(df$ReadsPerSpecies), "reads required")
  })
  
  output$read_depth <- renderText({
    df <- species_data()
    reads <- as.numeric(input$kit) * (1 - as.numeric(input$phix)/100)
    paste(round(reads / sum(df$ReadsPerSpecies), 0),
          "reads per sample per locus")
  })
  
  # Show cost
  output$cost_info <- renderText({
    ki <- kit_info()
    paste("Cost of this kit:", paste0("$", ki$Cost))
  })
  
  # Add weblink
  output$kit_link <- renderUI({
    #ki <- kit_info()
    tags$a(href = "https://youtu.be/93lrosBEW-Q?t=27", 
           "Click here for more information", 
           target = "_blank")
  })
  
  # Cost per genotype
  output$cost_per_genotype <- renderText({
    df <- species_data()
    ki <- kit_info()
    
    total_genotypes <- sum(df$Samples)
    
    if (total_genotypes > 0) {
      cost_pg <- ki$Cost / total_genotypes
      paste("Cost per genotype:", paste0("$", round(cost_pg, 4)))
    } else {
      "Please enter valid numbers for samples"
    }
  })

  
}

shinyApp(ui = ui, server = server)