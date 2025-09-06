
tab_one_ui <- function(id) {
  ns <- NS(id)
  nav_panel(
    title = "Tab One",
    h3("Welcome to Tab One"),
    textOutput(ns("text_one"))
  )
}

tab_one_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    output$text_one <- renderText({
      "Hello from Tab One!"
    })
  })
}
