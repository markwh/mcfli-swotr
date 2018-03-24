#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(leaflet)
library(dygraphs)
library(plotly)
library(shinydashboard)
library(DT)

load("cache/reachdata.RData")

header <- dashboardHeader(
  title = "Pepsi Manning Evaluation"
)

body <- dashboardBody(
  fluidRow(
    column(width = 2,
           box(width = NULL, 
               checkboxGroupInput("selcases", label = "Cases", 
                              choices = names(reachdata), 
                              selected = names(reachdata)))
    ),
    column(width = 5,
           box(title = "Manning Prediction Scatterplot",
               width = NULL, solidHeader = TRUE,
               plotOutput("allCasesScatter", height = 500))
    ),
    column(width = 5,
           box(title = "Manning Prediction Stats Boxplots",
               width = NULL, solidHeader = TRUE,
               plotOutput("allCasesBoxplot", height = 500))
    )
  ),
  fluidRow(
    column(width = 2,
           selectInput("selcase", label = "Case", choices = names(reachdata),
                       selected = names(reachdata)[1]),
           uiOutput("selReachesUI")
    ),
    column(width = 10,
      tabsetPanel(
        tabPanel("Hydrograph",
                 plotOutput("selcaseHydrograph")
        ),
        tabPanel("Stats Boxplots"),
        tabPanel("Partial Residuals",
                 column(width = 6,
                        plotOutput("selCasePartResid_W"),
                        plotOutput("selCasePartResid_A")
                 ),
                 column(width = 6,
                        plotOutput("selCasePartResid_S"),
                        plotOutput("selCasePartResid_SW")
                 )
        )
      )
    )
  )
)

dashboardPage(
  header,
  dashboardSidebar(disable = TRUE),
  body
)
