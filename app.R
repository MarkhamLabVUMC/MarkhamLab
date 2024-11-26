# This script creates the Markham Lab DESeq Dashboard. It has:
# Chatbot to fetch information about genes
# Various Plots with descriptions to analyze results from DESeq2 Analysis
# Make sure you have the 'www' folder that contains all the plots for this dashboard
# 'www' must be in same directory as app.R

# Load necessary libraries
library(shiny)
library(dplyr)

# Load your cleaned gene data
gene_data <- read.csv("D:\\MarkhamLab\\GeneScraping\\ncbi_summaries_cleaned.csv")

# Define the UI
ui <- fluidPage(
  titlePanel("Markham Lab DESeq Dashboard"),
  
  # Main tabsetPanel for the chatbot and interactive plots
  tabsetPanel(
    
    # Chatbot Tab
    tabPanel("Chatbot",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_list", "Enter Gene Symbols (comma-separated):", value = "BRCA1, TP53"),
                 checkboxGroupInput("info_type", "Select information to retrieve:", 
                                    choices = c("Description" = "description", "Summary" = "Summary"), 
                                    selected = "description"),
                 actionButton("submit", "Get Info"),
                 downloadButton("downloadData", "Download CSV")
               ),
               mainPanel(
                 tableOutput("results")
               )
             )
    ),
    
    # MA Plots grouped with dropdown to select the specific plot and option to compare two
    tabPanel("MA Plots",
             sidebarLayout(
               sidebarPanel(
                 selectInput("ma_plot_select1", "Select MA Plot:", 
                             choices = c("MA Plot" = "ma_plot_interactive.html",
                                         "MA Shrunken Plot" = "ma_shrunken_plot_interactive.html",
                                         "apeglm Plot" = "apeglm_plot_interactive.html",
                                         "Alternative Hypothesis MA" = "altHypothesis_MAplots.html")),
                 selectInput("ma_plot_select2", "Compare with:", 
                             choices = c("None" = "none", 
                                         "MA Plot" = "ma_plot_interactive.html",
                                         "MA Shrunken Plot" = "ma_shrunken_plot_interactive.html",
                                         "apeglm Plot" = "apeglm_plot_interactive.html",
                                         "Alternative Hypothesis MA" = "altHypothesis_MAplots.html")),
                 p("The MA plots visualize the relationship between the mean expression level of genes (x-axis) and the log2 fold change in expression between conditions (y-axis). Each point represents a gene, and genes that exhibit substantial expression changes and statistically significant p-values are highlighted. The standard MA plot gives an initial overview of differential expression, with red points indicating significant changes. The shrunken MA plot applies shrinkage to stabilize fold change estimates for low-count genes, reducing noise and enhancing reliability. Alternative hypothesis MA plots allow examination under specific regulatory patterns, highlighting only upregulated, downregulated, or absolute value changes. Together, these plots provide a comprehensive view of differential expression patterns, helping to pinpoint biologically meaningful gene changes.")  # Placeholder for context
               ),
               mainPanel(
                 uiOutput("ma_plot_1"),
                 uiOutput("ma_plot_2")
               )
             )
    ),
    
    # MeanSD Plots grouped with dropdown to select the specific plot and option to compare two
    tabPanel("MeanSD Plots",
             sidebarLayout(
               sidebarPanel(
                 selectInput("meansd_plot_select1", "Select MeanSD Plot:", 
                             choices = c("MeanSdPlot NTD" = "MeanSdPlot_ntd.html",
                                         "MeanSdPlot VSD" = "MeanSdPlot_vsd.html",
                                         "MeanSdPlot RLD" = "MeanSdPlot_rld.html")),
                 selectInput("meansd_plot_select2", "Compare with:", 
                             choices = c("None" = "none", 
                                         "MeanSdPlot NTD" = "MeanSdPlot_ntd.html",
                                         "MeanSdPlot VSD" = "MeanSdPlot_vsd.html",
                                         "MeanSdPlot RLD" = "MeanSdPlot_rld.html")),
                 p("The Mean-SD plot displays the relationship between the mean expression level (x-axis) and standard deviation (y-axis) for each gene across samples, under different transformations: normalized counts (NTD), variance-stabilizing transformation (VSD), and regularized log transformation (RLD). This plot allows for assessment of homoscedasticity in the data, showing whether variance remains constant across expression levels. Effective transformations stabilize variance, ensuring downstream analysis is not biased by variable dispersion across expression levels. This plot is essential for selecting the most appropriate transformation method, ultimately enhancing the reliability of differential expression results.")  # Placeholder for context
               ),
               mainPanel(
                 uiOutput("meansd_plot_1"),
                 uiOutput("meansd_plot_2")
               )
             )
    ),
    
    # Individual plot tabs with blurbs
    tabPanel("PCA Plot", 
             sidebarLayout(
               sidebarPanel(
                 p("The PCA plot displays sample relationships based on gene expression data, reducing high-dimensional data to two principal components (PC1 and PC2) on the x and y-axes, respectively. Each point represents a sample, with closer points indicating higher similarity in gene expression. Clustering of samples reveals patterns of variation driven by experimental conditions, allowing for detection of batch effects or outliers. This plot provides an overview of the main sources of variation in the dataset and helps in assessing overall data quality and the effectiveness of grouping conditions.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("pca_plot")
               )
             )
    ),
    tabPanel("Cook's Boxplot", 
             sidebarLayout(
               sidebarPanel(
                 p("This boxplot displays Cook's Distance values, a measure of each sample’s influence on the model’s fit. The y-axis shows log-transformed Cook's Distance values for each sample, helping identify influential samples that may disproportionately affect the results. Outliers or samples with high influence can distort the model, and detecting these is crucial for data quality control and accurate differential expression analysis.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("cooks_boxplot")
               )
             )
    ),
    
    tabPanel("Dispersion Plot", 
             sidebarLayout(
               sidebarPanel(
                 p("The Dispersion plot shows the relationship between gene expression levels (x-axis) and gene dispersion (y-axis), where dispersion reflects variability in gene expression across samples. Genes with high dispersion are more variable and may be prone to noise. This plot assesses how well the model accounts for expression variability, which is crucial for identifying reliable differentially expressed genes. The plot helps confirm that the model is appropriately fitted, especially for genes with varying expression levels.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("dispersion_plot")
               )
             )
    ),
    
    tabPanel("P-value Histogram", 
             sidebarLayout(
               sidebarPanel(
                 p("This histogram displays the distribution of p-values across all genes, with the x-axis showing p-value bins and the y-axis indicating frequency. A high number of low p-values suggests strong differential expression, while a uniform distribution may indicate noise or insufficient model fitting. This plot provides a visual summary of overall statistical significance in the dataset, serving as a diagnostic tool to evaluate the effectiveness of the differential expression analysis.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("pvalue_histogram")
               )
             )
    ),
    
    # Updated Heatmap and Volcano Plot with blurbs
    tabPanel("Heatmap", 
             sidebarLayout(
               sidebarPanel(
                 p("This heatmap visualizes the normalized expression levels of the top differentially expressed genes across samples (columns) and conditions (rows). Warmer colors indicate higher expression, and cooler colors indicate lower expression. The heatmap is clustered to reveal patterns and similarities in gene expression, helping to identify groups of genes or samples with shared biological functions or responses. This plot is particularly useful for spotting clusters of co-expressed genes and exploring their role in specific conditions or treatments.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("heatmap")
               )
             )
    ),
    
    tabPanel("Volcano Plot", 
             sidebarLayout(
               sidebarPanel(
                 p("The Volcano plot combines fold change (x-axis) and statistical significance (y-axis) for each gene, highlighting genes with both substantial fold changes and low p-values. Genes that meet a set significance threshold are color-coded, often with upregulated genes in one color and downregulated genes in another. This plot enables quick identification of the most strongly affected genes, allowing researchers to prioritize genes for further investigation based on both biological relevance (fold change) and statistical confidence.")  # Placeholder for context
               ),
               mainPanel(
                 htmlOutput("volcano_plot")
               )
             )
    )
  )
)

# Define the server logic
server <- function(input, output) {
  
  # Chatbot functionality with custom sorting and missing gene handling
  filtered_data <- reactive({
    req(input$gene_list)
    gene_list <- strsplit(input$gene_list, ",")[[1]]
    gene_list <- trimws(gene_list)  # Trim any white spaces
    
    # Filter the data for the selected genes and keep the order of input
    result <- gene_data %>%
      filter(Symbol %in% gene_list) %>%
      mutate(order = match(Symbol, gene_list)) %>%  # Match input order
      arrange(order) %>%
      select(-order)
    
    # Add "gene not found" for any missing genes
    missing_genes <- setdiff(gene_list, result$Symbol)
    if (length(missing_genes) > 0) {
      missing_data <- data.frame(
        Symbol = missing_genes,
        description = "gene not found",
        Summary = "gene not found"
      )
      result <- bind_rows(result, missing_data)
    }
    
    # Select the columns based on user input
    if ("description" %in% input$info_type & "Summary" %in% input$info_type) {
      result <- result
    } else if ("description" %in% input$info_type) {
      result <- result %>% select(Symbol, description)
    } else if ("Summary" %in% input$info_type) {
      result <- result %>% select(Symbol, Summary)
    }
    
    return(result)
  })
  
  # Show the filtered data in the UI
  output$results <- renderTable({
    filtered_data()
  })
  
  # Download the filtered data as a CSV
  output$downloadData <- downloadHandler(
    filename = function() { "gene_info.csv" },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  
  # Render MA plots with dropdown selections
  output$ma_plot_1 <- renderUI({
    tags$iframe(src = input$ma_plot_select1, height = 600, width = "100%")
  })
  
  output$ma_plot_2 <- renderUI({
    if (input$ma_plot_select2 != "none") {
      tags$iframe(src = input$ma_plot_select2, height = 600, width = "100%")
    }
  })
  
  # Render MeanSD plots with dropdown selections
  output$meansd_plot_1 <- renderUI({
    tags$iframe(src = input$meansd_plot_select1, height = 600, width = "100%")
  })
  
  output$meansd_plot_2 <- renderUI({
    if (input$meansd_plot_select2 != "none") {
      tags$iframe(src = input$meansd_plot_select2, height = 600, width = "100%")
    }
  })
  
  # Render individual plot tabs
  output$pca_plot <- renderUI({ tags$iframe(src = "pca_plot_interactive.html", height = 600, width = "100%") })
  output$cooks_boxplot <- renderUI({ tags$iframe(src = "Cooks_Boxplot.html", height = 600, width = "100%") })
  output$dispersion_plot <- renderUI({ tags$iframe(src = "Dispersion_Plot.html", height = 600, width = "100%") })
  output$pvalue_histogram <- renderUI({ tags$iframe(src = "pvalue_histogram.html", height = 600, width = "100%") })
  output$heatmap <- renderUI({ tags$iframe(src = "Heatmap_Highest_Genes.html", height = 600, width = "100%") })
  output$volcano_plot <- renderUI({ tags$iframe(src = "volcano_plot_labeled_interactive.html", height = 600, width = "100%") })
}

# Run the app
shinyApp(ui = ui, server = server)
