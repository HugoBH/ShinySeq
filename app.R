library(shiny)
library(bslib)

# Source tab scripts
source("tabs/tab1_explore.R")
source("tabs/tab2_results.R")

ui <- page_fillable(
  theme = bs_theme(bootswatch = "superhero"),
  titlePanel("Welcome to ShinySeq"),
  tags$p(style = "font-size: 16px; color: #ccc;",
         "A simple tool to explore run parameters on sequencing depth and costs for Genotyping by Sequencing"),
  
  navset_card_tab( 
    tab_one_ui("tab1"),
    tab_two_ui("tab2"), 
    nav_menu( 
      "Other links", 
      nav_panel("C", "Panel C content"), 
      "----", 
      "Description:", 
      nav_item( 
        a("Shiny", href = "https://shiny.posit.co", target = "_blank") 
      ) 
    ) 
  ), 
  id = "tab" 
)

server <- function(input, output, session) {
  tab_one_server("tab1")
  tab_two_server("tab2")
}

shinyApp(ui = ui, server = server)
