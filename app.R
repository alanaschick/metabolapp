## Alana Schick, March 2019
## Metabolomics data shiny app
## People: Ian Lewis, Ryan Groves

library(tidyverse)
library(shiny)
library(shinythemes)
library(pheatmap)
library(factoextra)
source("helper_functions.R")

## Read in data
## Remove this later so user can input their own data
dat <- read.csv("data/Maven_output_edit.csv", head = TRUE, stringsAsFactors = FALSE)

## User interface
ui <- fluidPage(theme = shinytheme("spacelab"),

	titlePanel("Exploring metabolomics data"),
	
	mainPanel(
	  
	  tabsetPanel(type = "tabs",
	              
	    ## Tab 1 - Input Data
	    tabPanel("Input Data",
	             br(),
	             fileInput(inputId = "file1", label = "Choose CSV File",
	                       multiple = FALSE,
	                       accept = c("text/csv","text/comma-separated-values","text/plain",".csv")),
	             checkboxInput("sandbox", "Use built-in example dataset instead", value = FALSE),
	             radioButtons("disp", "Display",
	                          choices = c(Head = "head", All = "all"),
	                          selected = "head"),
	             hr(),
	             tableOutput("contents"),
	             tableOutput("contents2")
	    ),
	    
	    ## Tab 2 - Plots
	    tabPanel("Plots",
	             br(),
	             
	             ## Choose scale
	             selectInput(inputId = "scale", label = "Choose scale:", choices = c("log", "durbin", "row", "column", "fold", "ratio", "zScore", "none"), selected = "log"),
	             
	             br(),
	             ## If fold, ratio, or zScore, scale by:
	             conditionalPanel(
	               condition = "input.scale == 'fold' | input.scale == 'ratio' | input.scale == 'zScore'",
	               "If scaling by 'fold', 'ratio', or 'zScore', choose sample to scale by.",
	               br(),
	               uiOutput("get_sample_names")
	             ),  
	             hr(),
	             
	             ## PCA plot
	             h2("PCA"),
	             fluidRow(
	               column(3,
	                 checkboxInput("individuals", "Plot individuals", value = TRUE)
	               ),
	               column(3,
	                 checkboxInput("variables", "Plot variables", value = FALSE)
	               ),
	               column(6,
	                 checkboxInput("colourbygroup", "Colour by sample type", value = FALSE)
	                 #checkboxInput("addellipses", "Add ellipses for sample type", value = FALSE)
	               )
	             ),  
	             plotOutput("pca"),
	             hr(),
	             
	             ## Heatmap
	             h2("Heatmap"),
	             br(),
	             plotOutput("heatmap", height = 500, width = 900),
	             hr(),
	             
	             ## Barplots
	             h2("Barplots"),
	             br(),
	             uiOutput("get_compound_names"),
	             uiOutput("get_group_names"),
	             plotOutput("barplot"),
	             hr(),
	             
	             ##Other plots
	             h2("Violin plots"),
	             br(), 
	             uiOutput("get_compounds_violin"), 
	             plotOutput("violinplots"),
	             br()
	    ),
	    
	    ## Tab 3 - Statistics?
	    tabPanel("Statistics",
	             br()
	             
	    )
	    
	    ## Tab 4 - Report? 
	    #tabPanel("Write report?")
	    
	  )
	)
)	
	


## Server logic
server <- function(input, output, session) {
  
  ## Determine which dataset to use
  df <- reactive({
    if (input$sandbox){dat} else {read.csv(input$file1$datapath, header = TRUE, stringsAsFactors = FALSE)}
  })  
  
  #df <- reactive({read.csv(input$file1$datapath, header = TRUE, stringsAsFactors = FALSE)})
  
  ## Get sample names from data
  output$get_sample_names <- renderUI({
    sample_names <- names(metaSep(df())[[2]])
    selectInput("scale_by","Samples:", sample_names, multiple = T)
  })
  
  ## Get compound names for barplots from data
  output$get_compound_names <- renderUI({
    compound_names <- metaSep(df())[[1]]$compound
    selectInput("compound_list", "Choose compounds to plot:", compound_names, multiple = T)
  })  
  
  ## Get sample group names from data
  output$get_group_names <- renderUI({
    group_names <- unique(sapply(strsplit(names(metaSep(df())[[2]]), split = "\\."), function (x) x[1]))
    selectInput("group_list", "Choose sample groups:", group_names, multiple = T)
  })  
  
  ## Get compound names for violin plots
  output$get_compounds_violin <- renderUI({
    compound_names <- metaSep(df())[[1]]$compound
    selectInput("compound_list_violin", "Choose compounds to plot:", compound_names, multiple = T)
  })
  
  ## Display data contents
  output$contents <- renderTable({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows, if selected, will be shown.
    
    req(input$file1)
    if (input$disp == "head"){
      return(head(df()))
    }
    else {
      return(df())
    }
  }) 
  
  output$contents2 <- renderTable({
    
    req(input$sandbox)
    if (input$disp == "head"){
      return(head(df()))
    }
    else {
      return(df())
    }
  })
	
  ## PCA plot
	output$pca <- renderPlot({
	  #req(input$file1)
	  mavPlot(df(), pca = T, scale = input$scale, scaleby = input$scale_by, individuals = input$individuals, variables = input$variables, pcacolour = input$colourbygroup)
	})
	
	## Heatmap
	output$heatmap <- renderPlot({
	  #req(input$file1)
	  #mavPlot(df(), heat = T, rCst = T, scale = input$scale, scaleby = input$scale_by)
	  mavPlot(df(), pheatmap = T, scale = input$scale, scaleby = input$scale_by)
	})
	
	## Barplots
	output$barplot <- renderPlot({
	  req(input$compound_list, input$group_list)
	  mavPlot(df(), barPlot = T, scale = input$scale, scaleby = input$scale_by, compound_list = input$compound_list, group_list = input$group_list)
	}) 
	
	## Violin plots
	output$violinplots <- renderPlot({
	  req(input$compound_list_violin)
	  mavPlot(df(), vioPlot = T, scale = input$scale, scaleby = input$scale_by, compound_list = input$compound_list_violin)
	})
	
}

shinyApp(ui = ui, server = server)
