library(shiny)
runGitHub( "ShinySeq", "HugoBH", ref = "main")


#connect to shiny.io
install.packages('rsconnect')

library(rsconnect)
rsconnect::deployApp("app.R")

runApp("app.R")
