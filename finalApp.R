#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

## Author: Mano Ranaweera
## mranawee@bu.edu
## BU BF591
## Final Project

# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(ggplot2)
library(colourpicker) # you might need to install this.
library(dplyr)
library('tidyverse')
library(PCAtools)
library('RColorBrewer')
library('readr')
library(grid)
library(gridExtra)
library(WGCNA)
library(ggbeeswarm)


# Define UI for application that draws a histogram
ui <- fluidPage(
  tabsetPanel(
    tabPanel('Samples',
      fileInput("sample_data", "Please input a sample information file: ", accept = ".csv"),
      tabsetPanel(
        tabPanel('Summary', DT::dataTableOutput('sum_table')),
        tabPanel('Table', DT::dataTableOutput('std_table')),
        tabPanel('Plots') #radioButtons("sampleParams", "Choose which column to plot:"))
                 
    )
    ),
    tabPanel('Counts',
      fileInput("norm_count_data", "Please input a normalized count matrix:", accept = ".csv"),
      sliderInput("var_percent", "Percentage of variance to filter genes by:", 0, 100, 25),
      sliderInput("non_zero", "Threshold for number of genes to be non-zero:", 0, 5000, 50),
      colourInput('col1', 'Choose a color for points below your threshold', 'red'), #for plotting
      colourInput('col2', 'Choose color for points above your threshold', 'black'),
      submitButton('show plot', icon('sync'), '400px'),
      
      tabsetPanel(
        tabPanel('Filtering summary', DT::dataTableOutput('filt_summary')),
        tabPanel('Diagnostic Variance Scatter Plot', plotOutput('norm_count_var_scatter')),
        tabPanel('Diagnostic Zeros scatter Plot', plotOutput('norm_count_zero_scatter')),
        tabPanel('Heatmap', plotOutput('heat')),
        tabPanel('PCA_plot', 
                 radioButtons('xpca', 'Please select a PC for the x-axis',
                              c('PC1',
                                'PC2',
                                'PC3',
                                'PC4',
                                'PC5',
                                'PC6',
                                'PC7',
                                'PC8',
                                'PC9')),
                 radioButtons('ypca', 'Please select a PC for the y-axis',
                              c('PC1',
                                'PC2',
                                'PC3',
                                'PC4',
                                'PC5',
                                'PC6',
                                'PC7',
                                'PC8',
                                'PC9')),
                 plotOutput('PCA_plt'))
      )
    ),
    tabPanel('Differential Expression',
             fileInput("de_results", "Please input a differential expression matrix", accept = ".csv"),
             radioButtons('button1', 'Please select a variable to plot on your x-axis', 
                          c('logFC',
                            'logCPM',
                            'LR',
                            'PValue',
                            'FDR')),
             radioButtons('button2', 'Please select a variable to plot on your y-axis', 
                          c('logFC',
                            'logCPM',
                            'LR',
                            'PValue',
                            'FDR')),
             colourInput('pval_col1', 'Choose a color for points below your adjusted p-value threshold', 'red'), #for plotting
             colourInput('pval_col2', 'Choose color for points above your adjusted p-value threshold', 'blue'),#for what meets p-value threshold
             sliderInput('slider', 'Select the magnitude of adjusted p-value you would like to color on the graph', -300, 0, -10),
             submitButton('show plot', icon('sync'), '400px'),
      tabsetPanel(
        tabPanel('Differential Expression Data Table',
                 
                 DT::dataTableOutput('table')
        ),
        tabPanel('Differential Expression Data Plot',
                 plotOutput('volcano')
        )
      )
    ),
    tabPanel('Individual Gene Expression',
             fileInput("cts", "Please input a normalized count matrix", accept = '.csv'),
             fileInput("samp_info", "Please input a sample information file", accept = '.csv'),
             radioButtons('categories', "Please select a categorical field",
                          c('Stage',
                            'Samples')),
             textInput("search", "Please search a gene of your choice(e.g Glyma06g12670.1, Glyma12g02076.11)", value = "", placeholder = "gene selection"),
             verbatimTextOutput("value"),
             radioButtons("plot_choice", "Please choose a type of plot",
                          c('bar',
                            'box',
                            'violin',
                            'beeswarm')),
             submitButton('show plot', icon('sync'), '400px'),
             tabsetPanel(
                tabPanel('Resulting Plot',
                      plotOutput('result'))
                )

      )
  )
)








  



# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #' load_Data
  #'
  #' @details Okay this one is a little weird but bear with me here. This is 
  #' still a "function", but it will take no arguments. The `reactive({})` bit 
  #' says "if any of my inputs (as in, input$...) are changed, run me again". 
  #' This is useful when a user clicks a new button or loads a new file. In 
  #' our case, look for the uploaded file's datapath argument and load it with 
  #' read.csv. Return this data frame in the normal return() style.
  
  #Loading data for all tabs
  
  #for loading sample data
  load_sample_data <- reactive({
    req(input$sample_data)
    return(read.csv(input$sample_data$datapath))
  })
  
  #for loading norm counts data
  load_norm_counts <- reactive({
    req(input$norm_count_data)
    return(read.csv(input$norm_count_data$datapath))
  })
  
  #for loading DE data
  load_DE_data <- reactive({
    req(input$de_results)
    return(read.csv(input$de_results$datapath))
  })
  
  #for loading individual gene expression data
  load_ind_counts <- reactive({
    req(input$cts)
    return(read.csv(input$cts$datapath))
  })
  
  load_samp_info <- reactive({
    req(input$samp_info)
    return(read.csv(input$samp_info$datapath))
  })
  
  
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param slider A negative integer value representing the magnitude of
  #' p-adjusted values to color. Most of our data will be between -1 and -300.
  #' @param color1 One of the colors for the points.
  #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  #' @details I bet you're tired of these plots by now. Me too, don't worry.
  #' This is _just_ a normal function. No reactivity, no bells, no whistles. 
  #' Write a normal volcano plot using geom_point, and integrate all the above 
  #' values into it as shown in the example app. The testing script will treat 
  #' this as a normal function.
  #' 
  #' !!sym() may be required to access column names in ggplot aes().
  #'
  #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    dataf$cutoff[dataf$PValue < 10^slider] <- 'BELOW PADJ MAGNITUDE'
    dataf$cutoff[dataf$PValue > 10^slider] <- 'ABOVE PADJ MAGNITUDE'
    plt <- dataf %>%
      ggplot() +
      geom_point(aes(x=!!sym(x_name), y=-log10(!!sym(y_name)), color = cutoff)) + #volcano points
      xlab(x_name) +
      ylab(paste('-log10',y_name, sep = ' ')) +
      scale_color_manual(values=c(color1, color2)) +
      ggtitle('Volcano Plot of DeSeq Differential Expression Results') +
      theme_minimal()
    
    return(plt)
  }
  
  diag_var_plot <- function(dataf, var_slider, zero_slider, color1, color2) { #x_name input from radion button
    #new columns for statistics
    dataf$rowVar <- apply(dataf[,-1], 1, var) #variance column
    dataf$num_zeros <- rowSums(dataf[,-1]==0) #number of zeros column
    dataf$median_count <- apply(dataf[,-1],1, median) #median column
    
    init_thresh <- var_slider/ 100
    final_thresh <- quantile(dataf$rowVar, probs = init_thresh)
    
    dataf$threshold[dataf$rowVar <= final_thresh | rowSums(dataf[,-1]==0) <= zero_slider] <- 'Passed'
    dataf$threshold[dataf$rowVar > final_thresh | rowSums(dataf[,-1]==0) > zero_slider] <- 'Failed' 
    
    #plot
    plt_variance <- dataf %>%
      ggplot() +
      geom_point(aes(x=!!sym('median_count'), y=!!sym('rowVar'), color = threshold)) + 
      xlab('median count') +
      ylab('variance') +
      scale_color_manual(values=c(color1, color2)) +
      ggtitle('Diagnostic Scatter Plot: Variance per gene') +
      theme_minimal()
    
    
    return(plt_variance)
  }
  
  diag_zero_plot <- function(dataf, var_slider, zero_slider, color1, color2) {
    dataf$rowVar <- apply(dataf[,-1], 1, var) #variance column
    dataf$num_zeros <- rowSums(dataf[,-1]==0) #number of zeros column
    dataf$median_count <- apply(dataf[,-1],1, median) #median column
    
    init_thresh <- var_slider/ 100
    final_thresh <- quantile(dataf$rowVar, probs = init_thresh)
    
    dataf$threshold[dataf$rowVar <= final_thresh | rowSums(dataf[,-1]==0) <= zero_slider] <- 'Passed'
    dataf$threshold[dataf$rowVar > final_thresh | rowSums(dataf[,-1]==0) > zero_slider] <- 'Failed' 
    
    #plot
    plt_zeros <- dataf %>%
      ggplot() +
      geom_point(aes(x=!!sym('median_count'), y=!!sym('num_zeros'), color = threshold)) + 
      xlab('median count') +
      ylab('Zeros') +
      scale_color_manual(values=c(color1, color2)) +
      ggtitle('Diagnostic Scatter Plot: Number of zeros per gene') +
      theme_minimal()
    return(plt_zeros)
  }
  
  #Tab with summary of table
  sample_summary <- function(dataf){
    ColumnName <- colnames(dataf)
    Type <- sapply(dataf, class)
    
    tib <- data.frame(ColumnName, Type) %>%
      as_tibble() %>%
      return()
    
  }
  
  std_tab <- function(dataf){
    return(read.csv(dataf))
  }
  
  dim <- function(dataf) {
    return(dim(dataf))
  }
  
  #Count summary
  counts_summary <- function(dataf, var_thresh, zeros_thresh) {
    num_samples <- ncol(dataf[,-1]) #number of samples
    num_genes <- nrow(dataf) #number of genes(rows)
    
    
    #variance <- apply(dataf[,-1], 1, var)
    
    dataf$rowVar<- apply(dataf[,-1], 1, var) #variance of each row
    #dataf$rowVar <- quantile(variance, probs = var_thresh) #percent variance of each row
    
    init_thresh <- var_thresh/ 100
    final_thresh <- quantile(dataf$rowVar, probs = init_thresh)
    
    dataf$threshold[dataf$rowVar <= final_thresh | rowSums(dataf[,-1]==0) <= zeros_thresh] <- 'Passed'
    dataf$threshold[dataf$rowVar > final_thresh | rowSums(dataf[,-1]==0) > zeros_thresh] <- 'Failed' 
    
    num_passed <- sum(dataf$threshold=="Passed") #number of passed genes
    portion_passed <- num_passed / num_genes
    percent_passed <- portion_passed * 100 #percentage of passed genes
    num_failed <- num_genes - num_passed #number of failed genes
    percent_failed <- (1 - portion_passed) * 100 #percent of failed genes
    
    sum_tib <- data.frame(num_samples, num_genes, num_passed, percent_passed, num_failed, percent_failed) %>%
      as_tibble() %>%
      return()

    
    #quantile(counts_matrix, probs = c(0.0, input$percentile))

  }
  
  heat_map <- function(dataf){
    #dataf.long <- pivot_longer(data = dataf,
    #                           cols = -c(1),
    #                           names_to = "Sample",
    #                           values_to = "Counts")
    #map <- ggplot(data = dataf.long, mapping = aes(x=Sample,
    #                                               y=X,
    #                                               fill = Counts))+
    #  geom_tile()+
    #  xlab(label = "Sample") +
    #  theme(axis.text.y = element_blank(), axis.ticks = element_blank()) %>%
    #  return()
    col.pal <- brewer.pal(8, 'RdBu')
    map <- heatmap(as.matrix(dataf[,-1]), col = col.pal)
    legend(x='left', legend = c('high','medium high','medium','medium low','low'),
           cex = 0.8, fill = col.pal)
    return(map)
  }
  
  plot_pca <- function(dataf, xpca, ypca) {
    pca <- pca(dataf[,-1])
    biplot(pca, x = xpca, y=ypca) %>%
      return()
  }
  

  
  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will 
  #' evaluate it normally. Not only does this function filter the data frame to 
  #' rows that are above the slider magnitude, it should also change the format 
  #' of the p-value columns to display more digits. This is so that it looks 
  #' better when displayed on the web page. I would suggest the function 
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_table <- function(dataf, slider) {
    filt_dataf <- filter(dataf, PValue< 10^slider) %>%  #dataf[!(dataf$padj > 10^slider),]
      mutate(PValue = formatC(.$PValue, digits = 4, format = 'e'))
             #padj = formatC(.$padj, digits = 4, format = 'e'))
    return(filt_dataf)
  }
  
  #merged data for individual gene expression
  ind_expression <- function(dataf, sampleInfo) {
    norm_data <- column_to_rownames(dataf, var = "X")
    t_norm_data <- t(norm_data)
    sample_info <- as.matrix(sampleInfo)
    
    merge <- cbind(sample_info, t_norm_data) %>%
      as_tibble() %>% mutate_at(-c(1,2), as.numeric) %>%
      return()
  }
  
  #Violin Plot
  violin_plot <- function(merged_data, group, gene) {
    plt <- merged_data %>% 
      ggplot(mapping = aes(x = !!sym(group), y=!!sym(gene))) +
      geom_violin() %>%
      return()
  }
  
  #Boxplot
  box_plot <- function(merged_data, group, gene) {
    plt <- merged_data %>%
      ggplot(mapping = aes(x=!!sym(group), y=!!sym(gene))) +
      geom_boxplot() %>%
      return()
  }
  
  #beeswarm plot
  bee_plot <- function(merged_data, group, gene) {
    plt <- merged_data %>%
      ggplot(mapping = aes(x=!!sym(group), y=!!sym(gene), color=!!sym(group), cex=2, size=2)) +
      geom_beeswarm() %>%
      return()
  }
  
  #bar plot
  bar_plot <- function(merged_data, group, gene){
    plt <- merged_data %>%
      ggplot(mapping = aes(x=!!sym(group), y=!!sym(gene)))+
      geom_bar(stat="identity")
  }
  
  #plot type
  condition <- function(type, merged_data, group, gene) {
    if(type == 'bar'){
      plot <- bar_plot(merged_data, group, gene)
    }else if(type == 'beeswarm'){
      plot <- bee_plot(merged_data, group, gene)
    }
    else if(type == 'box'){
      plot <- box_plot(merged_data, group, gene)
    }
    else if(type == 'violin'){
      plot <- violin_plot(merged_data, group, gene)
    }
    
    return(plot)
  }
  
  # Same sample info summary
  output$sum_table <- DT::renderDataTable(sample_summary(load_sample_data()))
  
  #Return standard summary table
  output$std_table <- DT::renderDataTable(load_sample_data())
  
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  
  
  output$filt_summary <- DT::renderDataTable(counts_summary(load_norm_counts(), input$var_percent, input$non_zero))
  
  #differential exp output
  output$volcano <- renderPlot(volcano_plot(load_DE_data(), input$button1, input$button2, input$slider, input$pval_col1, input$pval_col2)) # replace this NULL
  output$table <- DT::renderDataTable(draw_table(load_DE_data(), input$slider))
  
  #counts output
  output$norm_count_var_scatter <- renderPlot(diag_var_plot(load_norm_counts(), input$var_percent, input$non_zero, input$col1, input$col2))
  output$norm_count_zero_scatter <- renderPlot(diag_zero_plot(load_norm_counts(), input$var_percent, input$non_zero, input$col1, input$col2))
  output$heat <- renderPlot(heat_map(load_norm_counts()))
  output$PCA_plt <- renderPlot(plot_pca(load_norm_counts(), input$xpca, input$ypca))
  
  #individual output
  #output$value <- renderText({ input$caption })
  output$result <- renderPlot(condition(input$plot_choice, ind_expression(load_ind_counts(), load_samp_info()), input$categories, input$search))
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
