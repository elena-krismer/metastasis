library(shiny)
library(mixOmics)
# use subsystem library created in python
library(reticulate)
library(dbplyr)
library(dplyr)
library(tidyverse)
source_python("code_chunks.py")

# get sampling data from all three samples an combine it into dataframe
brain <- read.csv('sampling_brain_metastasis_100.csv') 
brain_l <- nrow(brain)
brain$group <- c(rep('brain', brain_l))


breast <- read.csv('sampling_breast_cancer_100.csv')  
breast_l <- nrow(breast)
breast$group <- c(rep('breast', breast_l))

lung <- read.csv('sampling_lung_metastasis_100.csv') 
lung_l <- nrow(lung)
lung$group <- c(rep('lung', lung_l))

common_react <- Reduce(intersect, list(colnames(lung),colnames(brain),colnames(breast)))
# combine all matrices
flux <- rbind(brain[common_react], breast[common_react], lung[common_react])

# get list of all documented subsystems in human1
list_subsystem = list_subsystems()

# subset flux dataframe into reactions from subsystem
plsda_for_subsystem <- function(subsystem, flux, ellipse_logical = TRUE){
  # get reactions from subsystem using python function
  reactions = subsystem_reactions(subsystem)
  subsystem_data <- flux %>% 
    as_data_frame() %>%
    select_if(names(.) %in% reactions)
  # create vector with groups
  group <- flux$group
  
  # normalize data
  subsystem_data <- scale(subsystem_data, center=FALSE, scale=colSums(subsystem_data))
  
  pls_result <- mixOmics::splsda(subsystem_data, group) 
  #plot <- plotIndiv(pls_result, legend = TRUE, ellipse = ellipse_logical, title = subsystem)
  return(pls_result)
}

groups <- c("Brain metastasis", "Lung metastasis", "Breast cancer", 
            "Brain tissue", "Lung tissue", "Breast tissue")


ui <- fluidPage(
  titlePanel("PLS_DA"),
  
  sidebarPanel(
    checkboxGroupInput("group", "Select groups", groups)
    selectInput("subsystem", "Select subsystem", list_subsystem),
  ),
  mainPanel(
    plotOutput("plsda_samples"),
    plotOutput("plsda_var")
  )
)

server <- function(input, output){
  output$plsda_samples <- renderPlot({
    x <- input$subsystem
    plsda_result <- plsda_for_subsystem(x, flux)
    plotIndiv(plsda_result, ellipse = TRUE, title = x, legend = TRUE)
  })
  output$plsda_var <- renderPlot({
    x <- input$subsystem
    plsda_result <- plsda_for_subsystem(x, flux)
    plotVar(plsda_result)
  })
}
shinyApp(ui, server)

