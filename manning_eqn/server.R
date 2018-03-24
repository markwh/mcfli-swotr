#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(plotly)
library(purrr)
library(dplyr)

load("cache/reachdata.RData")
load("cache/Qmats.RData")
source("lib/utils.R")

# Data prep

manNlist <- map(reachdata, manningN_list)
constNlist <- map(manNlist, median, na.rm = TRUE)

manQlist <- map2(reachdata, constNlist, function(x, y) manningQ_list(x, y))

trueQlist <- map(reachdata, `[[`, "Q")

manQdf <- manQlist %>% 
  map(~data.frame(manQ = as.vector(.), 
                  time = rep(1:ncol(.), each = nrow(.)),
                  loc = rep(1:nrow(.), ncol(.)))) %>% 
  bind_rows(.id = "case")

trueQdf <- trueQlist %>% 
  map(~data.frame(trueQ = as.vector(.), 
                  time = rep(1:ncol(.), each = nrow(.)),
                  loc = rep(1:nrow(.), ncol(.)))) %>% 
  bind_rows(.id = "case")

Qdf <- left_join(manQdf, trueQdf, by = c("time", "loc", "case"))

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$allCasesScatter <- renderPlot({
    
    scatter_gg <- Qdf %>% 
      filter(case %in% input$selcases) %>% 
      ggplot(aes(x = trueQ, y = manQ)) +
      geom_point(aes(color = case), alpha = 0.4) +
      scale_x_log10() + scale_y_log10() +
      theme_minimal() +
      geom_abline(slope = 1, intercept = 0)
    # ggplotly(scatter_gg)
    scatter_gg
    
  })
  
  output$allCasesBoxplot <- renderPlot({
    evalstatdf <- Qdf %>% 
      filter(case %in% input$selcases) %>% 
      group_by(case, loc) %>% 
      summarize(RRMSE = RRMSE(pred = manQ, meas = trueQ),
                NRMSE = NRMSE(pred = manQ, meas = trueQ),
                rBIAS = rBIAS(pred = manQ, meas = trueQ),
                NSE = NSE(pred = manQ, meas = trueQ)) %>% 
      ungroup()
    
    boxplots_gg <- evalstatdf %>%
      gather(key = stat, value = value, -case, -loc) %>% 
      ggplot(aes(x = stat, y = value)) +
      geom_boxplot() + 
      theme_bw()
    
    boxplots_gg
    # ggplotly(boxplots_gg)
  })
  
  selCaseData <- reactive({
    swot_tidy(reachdata[[input$selcase]]) %>% 
      mutate(n = W^(-2/3) * A^(5/3) * S^(1/2) / Q,
             n_const = median(n),
             Qhat = W^(-2/3) * A^(5/3) * S^(1/2) / n_const,
             Qpart_W = Q / (A^(5/3) * S^(1/2) / n_const),
             Qpart_S = Q / (A^(5/3) * W^(-2/3) / n_const),
             Qpart_A = Q / (W^(-2/3) * S^(1/2) / n_const),
             Qpart_SW = Q / (A^(5/3) / n_const),
             loc = as.character(loc))
  })
  
  selCaseData_filt <- reactive({
    if (is.null(input$selReaches)) {
      return(selCaseData())
    }
    selCaseData() %>% 
      filter(loc %in% input$selReaches)
  })
  
  output$selReachesUI <- renderUI({
    checkboxGroupInput(inputId = "selReaches", label = "Reaches",
                       choices = unique(selCaseData()$loc), 
                       selected = unique(selCaseData()$loc))
  })
  
  output$selcaseHydrograph <- renderPlot({
    selCaseData_filt() %>% 
      ggplot(aes(x = time, y = Q)) + 
      geom_line(aes(group = loc)) +
      geom_line(aes(group = loc, y = Qhat), linetype = 2)
  })
  
  output$selCasePartResid_W <- renderPlot({
    selCaseData_filt() %>% 
      ggplot(aes(x = W, y = Qpart_W)) +
      geom_point(aes(color = loc)) +
      stat_smooth(aes(group = loc), method = "lm",  
                  # formula = log10(y) ~ log10(x), 
                  se = FALSE) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw()
  })
  
  output$selCasePartResid_S <- renderPlot({
    selCaseData_filt() %>% 
      ggplot(aes(x = S, y = Qpart_S)) +
      geom_point(aes(color = loc)) +
      stat_smooth(aes(group = loc), method = "lm", 
                  # formula = log10(y) ~ log10(x), 
                  se = FALSE) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw()
  })
  
  output$selCasePartResid_A <- renderPlot({
    selCaseData_filt() %>% 
      ggplot(aes(x = A, y = Qpart_A)) +
      geom_point(aes(color = loc)) +
      stat_smooth(aes(group = loc, color = loc), method = "lm",  
                  # formula = log10(y) ~ log10(x), 
                  se = FALSE) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw()
  })
  
  output$selCasePartResid_SW <- renderPlot({
    selCaseData_filt() %>% 
      ggplot(aes(x = S^(1/2) * W^(-2/3), y = Qpart_SW)) +
      geom_point(aes(color = loc)) +
      scale_x_log10() +
      scale_y_log10() +
      geom_abline(aes(slope = 1, intercept = 0)) +
      theme_bw()
  })
  
})
