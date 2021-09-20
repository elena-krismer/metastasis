library(shiny)
library(mixOmics)
# use subsystem library created in python
library(reticulate)
library(dbplyr)
library(dplyr)
library(glue)
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
# ttest
###############################################################################

# return tt-test results in html format
ttest_results <- function(groups_list, reaction, flux) {
  flux <- flux[, colnames(flux) %in% c(reaction, "group")] %>% as.data.frame()
  html_output <- glue("<title>{reaction} T-test comparison </title><br>")
  for (group1 in groups_list) {
    for (group2 in groups_list) {
      comparisongroups <- glue("{group1} vs. {group2}<br>")
      sample1 <- flux[flux$group == group1, ]
      sample2 <- flux[flux$group == group2, ]
      # ttest remove group names
      results <- t.test(as.vector(sample1[1]), as.vector(sample2[1]))
      statistic_string <- toString(results$statistic)
      results_string <- glue("<p>P-value: {results$p.value} <span class='tab'></span> Statistic:  {statistic_string}<p> <br><br>")
      html_output <- paste(html_output, comparisongroups, results_string)
    }
  }
  return(HTML(html_output))
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
  "Glioblastoma", "GBM surrounding tissue", "Ependymoma",
  "Medullablastoma", "Pilocyticastrocytoma"
)



ui <- fluidPage(
  titlePanel("PLS_DA"),
  sidebarLayout(
    # ----------------------------------------------------------------------------
    # Side bar Panel
    # ----------------------------------------------------------------------------
    sidebarPanel(
      tags$form(
        checkboxGroupInput("group", "Select groups", groups),
        selectInput("subsystem", "Select subsystem", list_subsystem),
        sliderInput("cutoff", "Cut off  Variables with correlations below this cutoff in absolute value are not plotted ", value = 0.9, min = 0, max = 1)
      ),
      tags$form(
        textInput("reaction", "Get reaction description", "MAR06493"),
        actionButton("button", strong("Submit"))
      ),
      tags$form(
        h4(textOutput("reac")),
        h4(htmlOutput("description"))
      )
    ),
    # ----------------------------------------------------------------------------
    # Main Panel
    # ----------------------------------------------------------------------------
    mainPanel(
      # subsetting main panel into tabs
      tabsetPanel(
        type = "tabs",
        # plots
        tabPanel("Plot", plotOutput("plsda_samples"), plotOutput("plsda_var")),
        # reaction description
        tabPanel("Reaction Description", htmlOutput("table", height = 400, width = 800)),
        # link to kegg database
        tabPanel(
          "KEGG", h5(strong(textOutput("keggheader"))), textOutput("keggids"),
          textInput("reaction_search", "Reaction", "R02736"),
          htmlOutput("link")
        ),
        tabPanel(
          "Paired T-test", h5(strong(textOutput("ttestheader"))), textOutput("humanoneids"),
          textInput("reaction_ttest", "Reaction", "MAR00193"),
          htmlOutput("ttest_result", height = 800, width = 1600)
        )
      )
    )
  )
)



server <- function(input, output, session) {
  # ----------------------------------------------------------------------------
  # PLS-DA plot
  # ----------------------------------------------------------------------------
  output$plsda_samples <- renderPlot({
    # get subsystem inpit
    x <- input$subsystem
    # use selected groups for subsetting df
    g_list <- input$group
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = g_list)
    plotIndiv(plsda_result, ellipse = TRUE, title = x, legend = TRUE)
  })
  # PLS- DA PLots with variable labels
  output$plsda_var <- renderPlot({
    x <- input$subsystem
    g_list <- input$group
    cutoff_defined <- input$cutoff
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = g_list)
    plotVar(plsda_result, cutoff = cutoff_defined)
  })
  # ----------------------------------------------------------------------------
  # reaction description look up
  # ----------------------------------------------------------------------------

  string_reac <- eventReactive(input$button, {
    return(input$reaction)
  })
  output$reac <- renderText({
    string_reac()
  })
  string <- eventReactive(input$button, {
    description <- get_reaction_description(input$reaction)
    x <- glue('{description}<br><br><a href="https://metabolicatlas.org/explore/Human-GEM/gem-browser/reaction/{input$reaction}"> "View Reaction on Metabolic atlas"</a>')
    return(x)
  })
  output$description <- renderUI({
    HTML(string())
  })

  # ----------------------------------------------------------------------------
  # description table
  # ----------------------------------------------------------------------------
  html_table_output <- function() {
    x <- input$subsystem
    cutoff_defined <- input$cutoff
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = input$group)
    plotVar_data <- plotVar(plsda_result, cutoff = cutoff_defined)
    reaction_names <- rownames(plotVar_data)
    return(HTML(html_table(reaction_names)))
  }
  reaction_search_input <- output$table <- renderUI({
    html_table_output()
  })


  # ----------------------------------------------------------------------------
  # for KEGG look up
  # ----------------------------------------------------------------------------
  keggheader <- function() {
    cutoff <- input$cutoff
    header <- glue("KEGG IDs for Reaction with a correlation above: {cutoff}")
    return(header)
  }

  output$keggheader <- renderText({
    keggheader()
  })
  reaction_list <- function() {
    x <- input$subsystem
    cutoff_defined <- input$cutoff
    plsda_result <- plsda_for_subsystem(x, flux, groups_list = input$group)
    plotVar_data <- plotVar(plsda_result, cutoff = cutoff_defined)
    reaction_names <- rownames(plotVar_data)
    return(reaction_names)
  }
  output$keggids <- renderText(get_keggID(reaction_list()))
  getPage <- function() {
    kegg_link <- glue("https://www.genome.jp/dbget-bin/www_bget?rn:{input$reaction_search}")
    return(tags$iframe(
      src = kegg_link, style = "width:100%;", frameborder = "0",
      id = "iframe", height = "500px"
    ))
  }
  output$link <- renderUI({
    getPage()
  })

  # ----------------------------------------------------------------------------
  # T-test
  # ----------------------------------------------------------------------------
  ttestheader <- function() {
    cutoff <- input$cutoff
    header <- glue("HumanOne IDs for Reaction with a correlation above: {cutoff}")
    return(header)
  }
  output$ttestheader <- renderText({
    ttestheader()
  })
  output$humanoneids <- renderText(reaction_list())
  output$ttest_result <- renderUI({
    ttest_results(groups_list = input$group, reaction = input$reaction_ttest, flux = flux)
  })
}
shinyApp(ui, server)
