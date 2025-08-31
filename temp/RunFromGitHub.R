library(shiny)
runGitHub( "ShinySeq", "HugoBH", ref = "main")


#connect to shiny.io
install.packages('rsconnect')

library(rsconnect)
rsconnect::deployApp()

runApp("")


#Additional ideas
#Tab for user specific values (entry boxes instead of sliders)