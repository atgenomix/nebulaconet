library(shiny)
library(shinyjs)
library(DT)
library(promises)
library(future)
library(WGCNA) # Assume all WGCNA modules are loaded
library(shinycssloaders) # for withSpinner()
library(pheatmap)

NebulaCoNet_localdev <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
    plan(multisession, workers = parallel::detectCores() - 1)

    ui <- fluidPage(
        useShinyjs(),
        navbarPage(
            title = "NebulaCoNet (Local Upload)",
            tabPanel(
                title = "Data Upload & Preview",
                sidebarLayout(
                    sidebarPanel(
                        h4("Upload Files"),
                        fileInput("expression_matrix_file", "expression_matrix.csv", accept = ".csv"),
                        fileInput("gene_info_file", "gene_info.csv", accept = ".csv"),
                        fileInput("sample_info_file", "sample_info.csv", accept = ".csv"),
                        actionButton("load_data", "Load & Preview Data"),
                        width = 3
                    ),
                    mainPanel(
                        fluidRow(
                            column(
                                width = 12,
                                tabsetPanel(
                                    tabPanel(
                                        "Expression Matrix",
                                        withSpinner(DTOutput("wide_table_dt"), type = 6)
                                    ),
                                    tabPanel(
                                        "DEG Table",
                                        downloadButton("download_DEG", "Download DEG"),
                                        withSpinner(DTOutput("DEG_table"), type = 6)
                                    ),
                                    tabPanel(
                                        "Sample Info",
                                        withSpinner(DTOutput("sample_info_dt"), type = 6)
                                    )
                                )
                            )
                        )
                    )
                )
            ),
            tabPanel(
                title = "WGCNA",
                sidebarLayout(
                    sidebarPanel(
                        h4("Parallel Settings"),
                        numericInput("nThreads", "Number of Threads:",
                            value = max(1, parallel::detectCores() - 1),
                            min = 1, max = parallel::detectCores(), step = 1
                        ),
                        h4("Sample Clustering"),
                        selectInput("distMethod", "Distance Method",
                            choices = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                            selected = "euclidean"
                        ),
                        sliderInput("cutHeight", "Cutoff Height:", min = 0, max = 50, value = 15, step = 1),
                        hr(),
                        h4("Scale-Free Topology"),
                        sliderInput("powerRange", "Power Vector:", min = 1, max = 30, value = c(1, 20)),
                        numericInput("rsqCut", "RsquaredCut:", value = 0.8, min = 0, max = 1, step = 0.05),
                        numericInput("selectedPower", "Selected Power:", value = 6, min = 1, max = 30),
                        hr(),
                        h4("Module Detection"),
                        numericInput("deepSplit", "deepSplit:", value = 2, min = 0, max = 4, step = 1),
                        numericInput("minModuleSize", "Min Module Size:", value = 30, min = 5, max = 200, step = 1),
                        actionButton("runModules", "Detect Modules"),
                        width = 3
                    ),
                    mainPanel(
                        tabsetPanel(
                            tabPanel("Sample Tree", sampleClustUI("sample")),
                            tabPanel("Scale-Free Topology", sftUI("sft")),
                            tabPanel("Gene Modules", geneModuleUI("mod")),
                            tabPanel("Gene List", geneListUI("list")),
                            tabPanel(
                                "Module Heatmap",
                                fluidRow(
                                    column(
                                        width = 3,
                                        selectInput(
                                            inputId  = "heatmap_key",
                                            label    = "Sample ID 欄位 (mapping):",
                                            choices  = NULL,
                                            selected = NULL
                                        )
                                    ),
                                    column(
                                        width = 3,
                                        selectInput(
                                            "heatmap_phenotypes",
                                            "Select Phenotypes to Annotate:",
                                            choices = NULL,
                                            multiple = TRUE
                                        )
                                    ),
                                    column(
                                        width = 9,
                                        withSpinner(
                                            plotOutput("moduleHeatmap", height = "600px"),
                                            type = 6
                                        )
                                    )
                                )
                            )
                        ),
                        width = 9
                    )
                )
            )
        )
    )

    server <- function(input, output, session) {
        # Reactive containers
        results <- reactiveValues(
            normcount_data = NULL,
            exacttest_data = NULL,
            coldata        = NULL
        )
        wide_data <- reactiveVal(NULL)
        DEG_table <- reactiveVal(NULL)
        maeColData <- reactiveVal(NULL)

        observeEvent(input$load_data, {
            req(input$expression_matrix_file, input$gene_info_file, input$sample_info_file)

            withProgress(message = "Loading & filtering files...", value = 0, {
                # 1/6: Read expression_matrix.csv
                incProgress(1 / 6, detail = "Reading expression_matrix.csv")
                expr_mat <- data.table::fread(
                    input$expression_matrix_file$datapath,
                    data.table   = FALSE,
                    check.names  = FALSE
                )

                # 2/6: Read gene_info.csv
                incProgress(2 / 6, detail = "Reading gene_info.csv")
                gene_info <- data.table::fread(
                    input$gene_info_file$datapath,
                    data.table = FALSE
                )

                # 3/6: Read sample_info.csv
                incProgress(3 / 6, detail = "Reading sample_info.csv")
                sample_info <- data.table::fread(
                    input$sample_info_file$datapath,
                    data.table = FALSE
                )

                # 4/6: Check required columns
                incProgress(4 / 6, detail = "Checking columns")
                if (!"GeneSymbol" %in% colnames(expr_mat)) {
                    stop("expression_matrix.csv must contain a 'GeneSymbol' column.")
                }

                # 5/6: Gene filtering → Keep top 5000 genes by variance
                incProgress(5 / 6, detail = "Filtering top 5000 genes by variance")
                orig_n <- nrow(expr_mat)
                if (orig_n > 3000) {
                    # Only select numeric columns
                    data_mat <- expr_mat[, setdiff(colnames(expr_mat), "GeneSymbol"), drop = FALSE]
                    # Calculate variance for each gene (row)
                    gene_vars <- apply(data_mat, 1, var, na.rm = TRUE)
                    # Order by variance, keep top 3000
                    top_idx <- order(gene_vars, decreasing = TRUE)[1:3000]
                    expr_mat <- expr_mat[top_idx, , drop = FALSE]
                    showNotification(
                        sprintf("Filtered from %d genes to %d (Top 3000 by variance)", orig_n, nrow(expr_mat)),
                        type = "message"
                    )
                }

                # 6/6: Update Reactive values
                incProgress(6 / 6, detail = "Updating reactive values")
                results$normcount_data <- expr_mat
                results$exacttest_data <- gene_info
                results$coldata <- sample_info

                wide_data(expr_mat)
                DEG_table(gene_info)
                maeColData(sample_info)
            }) # end withProgress

            showNotification("Files loaded and filtered!", type = "message")
        })

        # 2. Expression preview
        output$wide_table_dt <- renderDT({
            req(wide_data())
            df <- wide_data()
            # Round numeric columns for display
            df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, 4) else x)
            datatable(df, options = list(pageLength = 20, autoWidth = TRUE))
        })

        # 3. DEG preview and download
        output$DEG_table <- renderDT(
            {
                req(DEG_table())
                df <- DEG_table()
                # Format columns
                if (all(c("logFC", "logCPM", "PValue", "FDR") %in% colnames(df))) {
                    df$logFC <- round(df$logFC, 5)
                    df$logCPM <- round(df$logCPM, 5)
                    df$PValue <- formatC(df$PValue, format = "e", digits = 5)
                    df$FDR <- formatC(df$FDR, format = "e", digits = 5)
                }
                datatable(df, options = list(pageLength = 20, autoWidth = TRUE))
            },
            server = FALSE
        )

        output$download_DEG <- downloadHandler(
            filename = function() "DEG_table.csv",
            content  = function(file) write.csv(DEG_table(), file, row.names = FALSE)
        )

        # 4. Sample Info preview
        output$sample_info_dt <- renderDT({
            req(maeColData())
            datatable(maeColData(), options = list(pageLength = 20, autoWidth = TRUE))
        })

        # 5. Enable WGCNA multithreading
        observe({
            req(input$nThreads)
            enableWGCNAThreads(input$nThreads)
        })

        # 6. When Expression Data is available, initialize WGCNA modules
        observeEvent(wide_data(), {
            req(wide_data())
            raw_df <- wide_data()
            rownames(raw_df) <- raw_df$GeneSymbol
            expr_mat <- t(as.matrix(raw_df[, setdiff(colnames(raw_df), "GeneSymbol")]))

            # Data quality check
            gsg <- goodSamplesGenes(expr_mat, verbose = 0)
            if (!gsg$allOK) expr_mat <- expr_mat[gsg$goodSamples, gsg$goodGenes]

            # Reactive wrapper for numeric conversion
            exprDataNumeric <- reactive({
                df <- as.data.frame(expr_mat)
                df[] <- lapply(df, as.numeric)
                df
            })

            # Call each Shiny module
            sampleOut <- sampleClustServer(
                "sample",
                exprData    = exprDataNumeric,
                distMethod  = reactive(input$distMethod),
                cutHeight   = reactive(input$cutHeight)
            )

            observeEvent(sampleOut$maxHeight(), {
                maxH <- sampleOut$maxHeight()
                updateSliderInput(session, "cutHeight",
                    min = 0, max = ceiling(maxH),
                    step = 1,
                    value = min(input$cutHeight, ceiling(maxH))
                )
            })

            sftServer(
                "sft",
                exprData      = sampleOut$filteredExpr,
                powerRange    = reactive(input$powerRange),
                rsqCut        = reactive(input$rsqCut),
                selectedPower = reactive(input$selectedPower)
            )

            modulesObj <- geneModuleServer(
                "mod",
                exprData   = sampleOut$filteredExpr,
                power      = reactive(input$selectedPower),
                deepSplit  = reactive(input$deepSplit),
                minSize    = reactive(input$minModuleSize),
                runTrigger = reactive(input$runModules)
            )
            observeEvent(maeColData(), {
                df <- maeColData()
                cols <- colnames(df)
                keys <- cols[sapply(cols, function(cn) length(unique(df[[cn]])) == nrow(df))]
                updateSelectInput(session, "heatmap_key",
                    choices  = keys,
                    selected = keys[1]
                )
                updateSelectInput(session, "heatmap_phenotypes",
                    choices  = setdiff(cols, keys[1]),
                    selected = character(0)
                )
            })

            output$moduleHeatmap <- renderPlot({
                req(modulesObj(), maeColData(), input$heatmap_key)

                me <- modulesObj()$MEs
                mat <- t(as.matrix(me))
                print(str(mat))
                print(input$heatmap_key)
                print(class(input$heatmap_key))

                samp <- maeColData()


                ids <- as.character(samp[[input$heatmap_key]])
                print(str(ids))
                if (any(duplicated(ids))) ids <- make.unique(ids)
                rownames(samp) <- ids
                print(str(samp))
                anno_col <- NULL
                if (length(input$heatmap_phenotypes) > 0) {
                    tmp <- samp[colnames(mat), input$heatmap_phenotypes, drop = FALSE]
                    print(str(tmp))
                    keep <- sapply(tmp, function(x) {
                        x2 <- x[!is.na(x)]
                        length(unique(x2)) > 1
                    })
                    anno_col <- tmp[, keep, drop = FALSE]
                    print(str(anno_col))
                }

                row_anno <- data.frame(Module = colnames(me))
                rownames(row_anno) <- colnames(me)

                pheatmap(
                    mat,
                    annotation_col = anno_col,
                    annotation_row = row_anno,
                    cluster_rows   = TRUE,
                    cluster_cols   = TRUE,
                    show_rownames  = TRUE,
                    show_colnames  = TRUE
                )
                print(dev.list())
            })

            geneListServer(
                "list",
                exprData    = sampleOut$filteredExpr,
                modulesObj  = modulesObj
            )
        })
    }

    # Launch App
    shinyApp(ui = ui, server = server)
}

options(shiny.maxRequestSize = 500 * 1024^2)

devtools::load_all("/Users/charleschuang/Desktop/atgenomix/ToolDev/oncoexpr")
NebulaCoNet_localdev()

dev.off()



