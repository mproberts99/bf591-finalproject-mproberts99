## Author: Monica Roberts
## BU BF591
## Final Project

library(shiny)
library(colourpicker) 
library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(DT)
library(pheatmap)
library(stringr)

options(shiny.maxRequestSize=30*1024^2)

ui <- fluidPage(
  titlePanel("BF591 Final Project"),
  p("Exploring RNA-Seq data from post-mortem brain tissue of patients who died from Huntington's Disease."),
  p("Reference: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810"),
  tabsetPanel(
    #---------------------------------------------------------Samples Tab------------------------------------------------------------
    tabPanel("Samples",
             p("Use this tab to explore the metadata of the samples in the study."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("samples_csvfile",
                           label = "Upload Sample Metadata:", 
                           accept = ".csv", 
                           placeholder = "metadata.csv"),
                 submitButton("Submit", 
                              width = '100%')
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            tableOutput("samples_summary")),
                   tabPanel("Table",
                            dataTableOutput("samples_table")),
                   tabPanel("Plots",
                            plotOutput("samples_plots"))
                 )
               )
             )
    ),
    #---------------------------------------------------------Counts Tab------------------------------------------------------------
    tabPanel("Counts",
             p("Use this tab to explore the counts data of the experiment and filter genes based on selected criteria."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("counts_csvfile",
                           label = "Upload counts matrix file.", 
                           accept = ".csv", 
                           placeholder = "counts.csv"),
                 sliderInput("counts_slider_var",
                             min = 0, 
                             max = 100,
                             label = "Select a minimum percentile of variance per gene:", 
                             value = 65, 
                             step = 1),
                 sliderInput("counts_slider_num",
                             min = 0, 
                             max = 70,
                             label = "Select a minimum number of non-zero samples per gene:", 
                             value = 60, 
                             step = 1),
                 submitButton("Plot", 
                              width = '100%')
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table",
                            tableOutput("counts_table")),
                   tabPanel("Filter Plots",
                            plotOutput("counts_filter_plots")),
                   tabPanel("Heatmap",
                            plotOutput("counts_heatmap")),
                   tabPanel("PCA",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("PC1_choice", 
                                                          "Choose the PC for the x-axis", 
                                                          choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                                          selected = "PC1"),
                                radioButtons("PC2_choice", 
                                             "Choose the PC for the y-axis", 
                                             choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
                                             selected = "PC2"),
                                submitButton("Plot", width='100%')
                            ),
                            mainPanel(
                              plotOutput("counts_PCA")
                            )
                            )
                   )
               )
             )
    )),
    
    #---------------------------------------------------------DE Tab------------------------------------------------------------
    tabPanel("DE",
             p("Use this tab to visualize differential expression results of the experiment."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("DE_csvfile",
                           label = "Load differential expression results", 
                           accept = ".csv", 
                           placeholder = "deseq_res.csv"),
                 radioButtons("DE_xchoice", 
                              "Choose the column for the x-axis", 
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                              selected = "log2FoldChange"),
                 radioButtons("DE_ychoice", 
                              "Choose the column for the y-axis", 
                              choices = c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj"),
                              selected = "padj"),
                 colourInput("DE_basecolor", 
                             label = "Base point color", 
                             value = "#FF57E9"),
                 colourInput("DE_highlightcolor", 
                             label = "Highlight point color", 
                             value = "#21EDCE"),
                 sliderInput(inputId = "DE_slider", 
                             min = -100, 
                             max = 0,
                             label = "Select the magnitude of the p adjusted coloring:", 
                             value = -5, 
                             step = 1),
                 submitButton("Plot", 
                              width = '100%')
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Plot",
                            plotOutput("DE_volcano")),
                   tabPanel("Table",
                            tableOutput("DE_table"))
                 )
               )
             )
    ),
    #---------------------------------------------------------GSEA Tab------------------------------------------------------------
    tabPanel("GSEA",
             p("Use this tab to visualize gene set enrichment analysis results of the data."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("GSEA_csvfile",
                           label = "Upload GSEA results table.", 
                           accept = ".csv", 
                           placeholder = "gsea_results.csv")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("NES Bar Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_Bar",
                                            min = 1, 
                                            max = 50,
                                            label = "Select number of top pathways:", 
                                            value = 40, 
                                            step = 1),
                                submitButton("Plot", 
                                             width = '100%')
                              ),
                              mainPanel(
                                plotOutput("GSEA_barplot")
                              )
                            )),
                   tabPanel("Table",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_table_pvalue",
                                            min = -20, 
                                            max = 0,
                                            label = "Select adjusted p-value to filter by:", 
                                            value = -1, 
                                            step = 1),
                                radioButtons("GSEA_pathways_choice",
                                             "Choose the types of pathways", 
                                             choices = c("All", "Positive", "Negative"),
                                             selected = "All"),
                                submitButton("Submit", width="100%")
                              ),
                              mainPanel(
                                sidebarLayout(
                                  sidebarPanel(
                                  downloadButton('download_NES_table', "Download Table")
                                ),
                                mainPanel(
                                dataTableOutput("GSEA_results_table")
                                )
                              ))
                            )),
                   tabPanel("NES Scatter Plot",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput("GSEA_scatter_pvalue",
                                            min = -20, 
                                            max = 0,
                                            label = "Select adjusted p-value to filter by:", 
                                            value = -1, 
                                            step = 1),
                                submitButton("Plot", width="100%")
                              ),
                              mainPanel(
                                plotOutput("GSEA_scatterplot")
                              )
                            ))
                 )
               )
             )
    )
  )
)

server <- function(input, output, session) {
  #---------------------------------------------------------Samples Tab------------------------------------------------------------
  samples_load_data <- reactive({
    req(input$samples_csvfile)
    file=input$samples_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
      return()
  })
  
  samples_table_summary <- function(dataf) {
    column <- colnames(dataf)
    column_type <- as.vector(sapply(dataf, typeof))
    df <- data.frame(column=column, column_type=column_type)
    num_sum <- data.frame(mean = sapply(dataf, mean), 
                          sds = sapply(dataf, sd)) %>% 
      mutate(summary = str_glue('{mean}(+/-{sds})'))
    df$numeric_summary <- num_sum$summary
    df$numeric_summary <- gsub('NA(+/-NA)', '', df$numeric_summary, fixed = TRUE)
    df_nonnum <- data.frame(sum = sapply(dataf, function(x) str_glue("Distinct Values: {length(unique(x))}")))
    df$nonnumeric_summary <- df_nonnum$sum
    return(df)
  }
  
  samples_plot <- function(dataf) {
    colnames(dataf) <- gsub('-', '_', colnames(dataf))
    p1 <- dataf %>% ggplot(aes(x=age_of_death)) +
      geom_density()
    p2 <- dataf %>% ggplot(aes(x=Bases)) +
      geom_density()
    p3 <- dataf %>% ggplot(aes(x=Bytes)) +
      geom_density()
    p4 <- dataf %>% ggplot(aes(x=mrna_seq_reads)) +
      geom_density()
    p5 <- dataf %>% ggplot(aes(x=pmi, na.rm=TRUE)) +
      geom_density()
    p6 <- dataf %>% ggplot(aes(x=rin)) +
      geom_density()
    p7 <- dataf %>% ggplot(aes(x=age_of_onset, na.rm=TRUE)) +
      geom_density()
    p8 <- dataf %>% ggplot(aes(x=cag)) +
      geom_density()
    p9 <- dataf %>% ggplot(aes(x=Duration)) +
      geom_density()
    p10 <- dataf %>% ggplot(aes(x=h_v_cortical_score, na.rm=TRUE)) +
      geom_density()
    p11 <- dataf %>% ggplot(aes(x=h_v_striatal_score, na.rm=TRUE)) +
      geom_density()
    p12 <- dataf %>% ggplot(aes(x=vonsattel_grade, na.rm=TRUE)) +
      geom_density()
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)
    return()
  }
  
  output$samples_summary <- renderTable({samples_table_summary(dataf=samples_load_data())})
  
  output$samples_table <- renderDataTable({
    samples_load_data()
  })
  
  output$samples_plots <- renderPlot({
    samples_plot(samples_load_data())
  })
  
  
  #---------------------------------------------------------Counts Tab------------------------------------------------------------
  counts_load_data <- reactive({
    req(input$counts_csvfile)
    file=input$counts_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
    return()
  })
  
  counts_filter_summary <- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    num_samples = ncol(stats_df)-4
    num_genes = nrow(stats_df)
    num_genes_filtered = nrow(filtered)
    perc_genes_filtered = (num_genes_filtered/num_genes)*100
    num_genes_out = num_genes - num_genes_filtered
    perc_genes_out = (num_genes_out/num_genes)*100
    final_stats <- tibble(num_samples=num_samples,
                          num_genes=num_genes,
                          num_genes_filtered=num_genes_filtered,
                          percent_filtered=perc_genes_filtered,
                          num_genes_filtered_out=num_genes_out,
                          percent_filtered_out=perc_genes_out)
    return(final_stats)
  }
  
  counts_filter_plots <- function(dataf, slider_var, slider_num){
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(dataf[,-1]!=0)
    stats_df$median <- apply(dataf[,-1], 1, median)
    filtered_status <- stats_df %>% 
      mutate(filter_status=ifelse(rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num, 'pass filter', 'filtered out'))
    p1 <- ggplot(filtered_status, aes(x=variance, y=median, color=filter_status)) +
      geom_point() +
      scale_y_log10() +
      scale_x_log10() +
      labs(x='Percentile of Variance', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Percentile Variance') +
      theme_dark() +
      scale_color_brewer(palette="Greens")
    p2 <- ggplot(filtered_status, aes(x=nonzero, y=median, color=filter_status)) +
      geom_point() +
      scale_y_log10() +
      labs(x='Number of Non-Zero Samples', y='Median Number of Counts', color='Filter Status', main='Median Counts vs. Number of Non-Zero Samples') +
      theme_dark() +
      scale_color_brewer(palette="Oranges")
    grid.arrange(p1, p2, nrow = 1)
    return()
  }
  
  counts_heatmap <- function(dataf, slider_var, slider_num) {
    stats_df <- dataf
    stats_df$variance <- apply(dataf[,-1], 1, var)
    stats_df <- arrange(stats_df, by=variance)
    stats_df$rank <- rank(stats_df$variance)
    stats_df$nonzero <- rowSums(stats_df!=0)
    filtered <- filter(stats_df, rank >= nrow(stats_df)*(slider_var/100) & nonzero >= slider_num)
    mat <- as.matrix(filtered[,2:70])
    rownames(mat) <- filtered$gene
    p <- pheatmap(log10(mat+1),
                 scale="row",
                 color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(11),
                 cluster_rows=TRUE,
                 cluster_cols=TRUE,
                 show_rownames=FALSE,
                 show_colnames=FALSE,
                 legend=TRUE
                 )
  }
  
  counts_PCA <- function(dataf, PC1_choice, PC2_choice) {
    expr_mat <- t(as.matrix(dataf[,-1]))
    pca <- prcomp(
      expr_mat,
      center=TRUE,
      scale=TRUE 
    )
    pca_var <- tibble(
      PC=factor(str_c("PC",1:69),str_c("PC",1:69)),
      Variance=pca$sdev**2,
      Explained_Variance=Variance/sum(Variance)*100,
      `Cumulative % Explained Variance`=cumsum(Explained_Variance)
    )
    as_tibble(pca$x) %>%
      ggplot(aes(x=!!sym(PC1_choice),y=!!sym(PC2_choice))) +
        geom_point() +
        labs(x=str_glue("{PC1_choice} Percent Variance Explained: {pca_var$Explained_Variance[pca_var$PC=={PC1_choice}]}"), 
             y=str_glue("{PC2_choice} Percent Variance Explained: {pca_var$Explained_Variance[pca_var$PC=={PC2_choice}]}")) %>% 
    return()
  }
  
  output$counts_table <- renderTable({counts_filter_summary(dataf = counts_load_data(), 
                                                            slider_var = input$counts_slider_var, 
                                                            slider_num=input$counts_slider_num)})
  
  output$counts_filter_plots <- renderPlot({counts_filter_plots(dataf = counts_load_data(),
                                                                slider_var = input$counts_slider_var, 
                                                                slider_num = input$counts_slider_num)})
  
  output$counts_heatmap <- renderPlot({counts_heatmap(dataf = counts_load_data(),
                                                      slider_var = input$counts_slider_var,
                                                      slider_num = input$counts_slider_num)})
  
  output$counts_PCA <- renderPlot({counts_PCA(dataf = counts_load_data(),
                                                      PC1_choice = input$PC1_choice,
                                                      PC2_choice = input$PC2_choice)})
  #---------------------------------------------------------DE Tab------------------------------------------------------------
  DE_load_data <- reactive({
    req(input$DE_csvfile)
    file=input$DE_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read.csv(file$datapath, header= TRUE, sep=",")}%>%
      rename(Gene = X)%>%
      return()
  })
  
  #' Volcano plot
  DE_volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    p <- ggplot(data = dataf, 
                aes(x =!!sym(x_name), y=-log10(!!sym(y_name)))) + 
      geom_point(aes(color = padj< 1*10^(slider))) +
      theme(legend.position = "bottom") +
      scale_color_manual(values = c('TRUE' = color1, 'FALSE' = color2)) +
      labs(x=x_name, y=str_glue("-log10({y_name})"), color=str_glue("{y_name} < 10^{slider}"))
    return(p)
  }
  
  #' Draw and filter table 
  DE_draw_table <- function(dataf, slider) {
    dataf %>%
      arrange(pvalue) %>% 
      filter(padj < 10^slider) %>%
      mutate(pvalue = formatC(.$pvalue, digits = 2, format = "e" ),
             padj = formatC(.$padj, digits = 2, format = "e")) %>%
      return()
  }
  
  # return volcano output 
  output$DE_volcano <- renderPlot({DE_volcano_plot(dataf = DE_load_data(), 
                                                   slider = input$DE_slider, 
                                                   x_name=input$DE_xchoice, 
                                                   y_name=input$DE_ychoice, 
                                                   color1=input$DE_basecolor, 
                                                   color2=input$DE_highlightcolor)})
  
  #return table output
  output$DE_table <- renderTable({DE_draw_table(dataf = DE_load_data(), 
                                                slider = input$DE_slider)})
  
  #---------------------------------------------------------GSEA Tab------------------------------------------------------------
  GSEA_load_data <- reactive({
    req(input$GSEA_csvfile)
    file=input$GSEA_csvfile
    if (is.null(file))
    {return(NULL)} 
    else
    {datafile=read_csv(file$datapath)} %>% 
      return()
  })
  
  GSEA_bar <- function(dataf, pathways_slider) {
    filtered <- dataf %>% arrange(padj) %>% slice_head(n=pathways_slider) %>% mutate(status= ifelse(NES > 0, 'positive', 'negative'))
    filtered$pathway <- gsub("\\_", " ", filtered$pathway)
    filtered <- filtered %>% arrange(status)
    pathways <- factor(filtered$pathway)
    filtered$pathway <- factor(filtered$pathway, levels=unique(pathways))
    status <- factor(filtered$status)
    filtered$status <- factor(filtered$status, levels=unique(status))
    plot <- filtered %>% ggplot(aes(x=stringr::str_wrap(pathway, 40), y=NES, fill=status)) +
      geom_col() +
      coord_flip() +
      theme_light() +
      theme(axis.text = element_text(size = 5)) +
      xlab("")
      return(plot)
  }
  
  GSEA_table <- function(dataf, pvalue_choice, type_pathways_choice) {
    data <- dataf %>% mutate(status = ifelse(NES > 0, 'Positive', 'Negative')) %>% filter(padj<=10**pvalue_choice)
    if (type_pathways_choice != 'All') {
      filtered <- data %>% filter(status == type_pathways_choice)
      return(filtered)
    }
    else {
      return(data)
    }
  }
  
  GSEA_scatter <- function(dataf, pvalue_choice) {
    filtered <- dataf %>% mutate(filter_status=ifelse(padj<10**pvalue_choice, "passed filter", "filtered out"))
    plot <- filtered %>% ggplot(aes(x=NES, y=-log10(padj), color=filter_status)) +
      geom_point()
    return(plot)
  }
    
  
  output$GSEA_barplot <- renderPlot({GSEA_bar(dataf = GSEA_load_data(),
                                                         pathways_slider = input$GSEA_Bar)})
  
  output$GSEA_results_table <- renderDataTable({GSEA_table(dataf=GSEA_load_data(),
                                                           pvalue_choice=input$GSEA_table_pvalue,
                                                           type_pathways_choice=input$GSEA_pathways_choice)})
  
  output$GSEA_scatterplot <- renderPlot({GSEA_scatter(dataf=GSEA_load_data(),
                                                      pvalue_choice=input$GSEA_scatter_pvalue)})
  
  output$download_NES_table <- downloadHandler(
    filename = function(){"NES_Results.csv"},
    content = function(filename){
      write.csv(GSEA_table(GSEA_load_data(), input$GSEA_table_pvalue, input$GSEA_pathways_choice), filename)
    })
}
shinyApp(ui = ui, server = server) # run app