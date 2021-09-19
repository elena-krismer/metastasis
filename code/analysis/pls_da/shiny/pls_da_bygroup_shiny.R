library(shiny)
library(mixOmics)
# use subsystem library created in python
library(reticulate)
library(dbplyr)
library(dplyr)
library(tidyverse)

source_python("code_chunks.py")

###############################################################################
# Extract sampling data
# create column specifying group
# combine to df
###############################################################################

brain_met <- read.csv("sampling_brain_metastasis_100.csv")
# create column with group
brain_l <- nrow(brain_met)
brain_met$group <- c(rep("Brain metastasis", brain_l))

breast_cancer <- read.csv("sampling_breast_cancer_100.csv")
breast_l <- nrow(breast_cancer)
breast_cancer$group <- c(rep("Breast cancer", breast_l))

lung_met <- read.csv("sampling_lung_metastasis_100.csv")
lung_l <- nrow(lung_met)
lung_met$group <- c(rep("Lung metastasis", lung_l))

brain_tissue <- read.csv("sampling_brain_tissue_100.csv")
brain_l <- nrow(brain_tissue)
brain_tissue$group <- c(rep("Brain tissue", brain_l))

breast_tissue <- read.csv("sampling_breast_tissue_100.csv")
breast_l <- nrow(breast_tissue)
breast_tissue$group <- c(rep("Breast tissue", breast_l))

lung_tissue <- read.csv("sampling_lung_tissue_100.csv")
lung_l <- nrow(lung_tissue)
lung_tissue$group <- c(rep("Lung tissue", lung_l))

ependymoma <- read.csv("sampling_ependymoma_100.csv")
l <- nrow(ependymoma)
ependymoma$group <- c(rep("Ependymoma", l))

glioblastoma <- read.csv("sampling_glioblastoma_100.csv")
l <- nrow(glioblastoma)
glioblastoma$group <- c(rep("Glioblastoma", l))

gbm_sur <- read.csv("sampling_gbm_surrounding_tissue_100.csv")
l <- nrow(gbm_sur)
gbm_sur$group <- c(rep("GBM surrounding tissue", l))

medullablastoma <- read.csv("sampling_medullablastoma_100.csv")
l <- nrow(medullablastoma)
medullablastoma$group <- c(rep("Medullablastoma", l))

pilocyticastrocytoma <- read.csv("sampling_pilocyticastrocytoma_100.csv")
l <- nrow(pilocyticastrocytoma)
pilocyticastrocytoma$group <- c(rep("Pilocyticastrocytoma", l))


# summarize reactions from all models
all_reactions <- unique(c(
  colnames(brain_met), colnames(breast_cancer), colnames(lung_met),
  colnames(brain_tissue), colnames(breast_tissue), colnames(lung_tissue),
  colnames(ependymoma), colnames(gbm_sur), colnames(glioblastoma),
  colnames(medullablastoma), colnames(pilocyticastrocytoma)
))

# fill in missing reactions as NA
brain_met[, setdiff(all_reactions, colnames(brain_met))] <- NA
brain_tissue[, setdiff(all_reactions, colnames(brain_tissue))] <- NA
breast_tissue[, setdiff(all_reactions, colnames(breast_tissue))] <- NA
breast_cancer[, setdiff(all_reactions, colnames(breast_cancer))] <- NA
lung_met[, setdiff(all_reactions, colnames(lung_met))] <- NA
lung_tissue[, setdiff(all_reactions, colnames(lung_tissue))] <- NA
ependymoma[, setdiff(all_reactions, colnames(ependymoma))] <- NA
glioblastoma[, setdiff(all_reactions, colnames(glioblastoma))] <- NA
gbm_sur[, setdiff(all_reactions, colnames(gbm_sur))] <- NA
medullablastoma[, setdiff(all_reactions, colnames(medullablastoma))] <- NA
pilocyticastrocytoma[, setdiff(all_reactions, colnames(pilocyticastrocytoma))] <- NA

# combine dataframes
flux <- rbind(brain_met, breast_cancer, lung_met, brain_tissue, breast_tissue,
  lung_tissue, glioblastoma, gbm_sur, pilocyticastrocytoma, ependymoma, medullablastoma,
  fill = TRUE
)
# set NA to zero flux
flux[is.na(flux)] <- 0



###############################################################################
# PLS-DA plot
###############################################################################

# subset flux dataframe into reactions from subsystem
plsda_for_subsystem <- function(subsystem, flux, groups_list, ellipse_logical = TRUE) {
  # get reactions from subsystem using python function
  reactions <- subsystem_reactions(subsystem)
  # extract rows with specified groups
  flux <- subset(flux, group %in% groups_list)
  # create vector with groups
  group <- flux$group
  # remove group column
  flux$group <- NULL
  # subset into subsystem dataframe
  subsystem_data <- flux[, colnames(flux) %in% reactions]
  # normalize data
  subsystem_data <- scale(subsystem_data, center = FALSE, scale = colSums(subsystem_data))

  pls_result <- mixOmics::splsda(subsystem_data, group)
  return(pls_result)
}


###############################################################################
# Shiny
###############################################################################

# get list of all documented subsystems in human1 for dropdownmenu
list_subsystem <- list_subsystems()

# groups for dropdown menu
groups <- c(
  "Brain metastasis", "Lung metastasis", "Breast cancer",
  "Brain tissue", "Lung tissue", "Breast tissue",
  "Glioblastoma","GBM surrounding tissue", "Ependymoma", 
  "Medullablastoma", "Pilocyticastrocytoma"
)


ui <- fluidPage(
  titlePanel("PLS_DA"),
  sidebarLayout(
  sidebarPanel(
    tags$form(
    checkboxGroupInput("group", "Select groups", groups),
    selectInput("subsystem", "Select subsystem", list_subsystem),
    sliderInput("cutoff", "Cut off  Variables with correlations below this cutoff in absolute value are not plotted ", value = 0.9, min = 0, max = 1),
    #textInput("reaction", "Get reaction description"),
    ),
    tags$form(
      textInput("reaction", "Get reaction description", "MAR06493"),
      actionButton("button", strong("Submit"))
    ),
    tags$form(
      h5(textOutput("reac")),
      h5(textOutput("description"))
    )
  ),
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Plot", plotOutput("plsda_samples"),
                         plotOutput("plsda_var")),
                tabPanel("Metabolic Atlas", 
                         htmlOutput("link")))
    )
  )
)


 
server <- function(input, output) {
  output$plsda_samples <- renderPlot({
    x <- input$subsystem
    g_list <- input$group
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = g_list)
    plotIndiv(plsda_result, ellipse = TRUE, title = x, legend = TRUE)
  })
  output$plsda_var <- renderPlot({
    x <- input$subsystem
    g_list <- input$group
    cutoff_defined <- input$cutoff
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = g_list)
    plotVar(plsda_result, cutoff = cutoff_defined)
  })
  string <- eventReactive(input$button, {
    x <-  get_reaction_description(input$reaction)
    if(length(x) == 0){
      x <- glue("https://metabolicatlas.org/explore/Human-GEM/gem-browser/reaction/{input$reaction}")
    }
    return(x)
  })
  string_reac <- eventReactive(input$button, {
    return(input$reaction)
  })
  output$reac <- renderText({
    string_reac()
  })   
  output$description <- renderText({
    string()
  }) 
  getPage<-function() {
    return(tags$iframe(src = "https://www.genome.jp/pathway/map00190"))
  }
  output$link<-renderUI({
    getPage()
  })
}
shinyApp(ui, server)
