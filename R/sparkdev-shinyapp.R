NebulaCoNet_sparkdev <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
    plan(multisession, workers = parallel::detectCores() - 1)

    ui <- fluidPage(
        useShinyjs(),
        navbarPage(
            title = "NebulaCoNet (Spark)",
            tabPanel(
                title = "Data Upload & Preview",
                sidebarLayout(
                    sidebarPanel(
                        h4("Upload Files"),
                        #fileInput("expression_matrix_file", "expression_matrix.csv", accept = ".csv"),
                        #fileInput("gene_info_file", "gene_info.csv", accept = ".csv"),
                        #fileInput("sample_info_file", "sample_info.csv", accept = ".csv"),
                        dbBrowserUI("dbBrowser1"),
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
                                        withSpinner(DTOutput("DEG_table_view"), type = 6)
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
                            tabPanel("Sample Tree", sparksampleClustUI("sample")),
                            tabPanel("Scale-Free Topology", sparksftUI("sft")),
                            tabPanel("Gene Modules", sparkgeneModuleUI("mod")),
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
            db_info = NULL,
            table_list = NULL,
            grouplist = NULL,
            normcount_data = NULL,
            exacttest_data = NULL,
            coldata = NULL
        )
        volcano_res <- reactiveVal(NULL)
        settingMAE <- reactiveVal(NULL)
        DEG_table <- reactiveVal(NULL)
        DEG_summary <- reactiveVal(NULL)
        wide_data <- reactiveVal(NULL)
        maeColData <- reactiveVal(NULL)
        if (is.null(master)) {
            if (Sys.getenv("SPARK_K8S_ENDPOINT") != "") {
                master <- Sys.getenv("SPARK_K8S_ENDPOINT")
            } else{
                master <- "sc://172.18.0.1:15002"
            }
        } 
        sc <- sparklyr::spark_connect(master = master, method = method, version = version)
        print(sc)

        session$onSessionEnded(function() {
            if (!is.null(sc)) {
                sparklyr::spark_disconnect(sc)
                message("Spark connection disconnected.")
            }
        })

        observeEvent(sc,
            {
                req(sc)
                print("dbbrowser initialized")
                shinyjs::disable("dbBrowser1-selected_db")
                results$db_info <- dbBrowserServer("dbBrowser1", sc)
                showNotification("Waiting for initialization", type = "message", duration = 10)

                normcount_future <- trigger_cluster_query_by_pattern(
                    master, method, version,
                    pattern = "^(normcounts|expression_table)",
                    output_label = "init_tbl_normcount"
                )

                exacttest_future <- trigger_cluster_query_by_pattern(
                    master, method, version,
                    pattern = "^(exacttest|gene_table|gene_metadata)",
                    output_label = "init_tbl_exacttest"
                )

                coldata_future <- trigger_cluster_query_by_pattern(
                    master, method, version,
                    pattern = "^(coldata|sample_table|sample_metadata)",
                    output_label = "init_tbl_coldata"
                )
                all_promises <- promises::promise_all(
                    norm = normcount_future,
                    ex   = exacttest_future,
                    col  = coldata_future
                )

                all_promises %...>% (function(res_list) {
                    shinyjs::enable("dbBrowser1-selected_db")
                    showNotification("Initialization complete. Check list!", type = "message", duration = 10)
                }) %...!% (function(e) {
                    shinyjs::enable("dbBrowser1-selected_db")
                    showNotification(paste("Error:", e$message), type = "error")
                })
            },
            ignoreInit = FALSE
        )

        progressMod <- progressPopupServer("popupProgress")

        observeEvent(results$db_info$selected_db(), {
            req(results$db_info$selected_db())
            shinyjs::disable("dbBrowser1-selected_db")
            
            # output$wide_table_dt <- DT::renderDataTable({
            #     data.frame()
            # })
            # output$DEG_table_view <- DT::renderDataTable({
            #     data.frame()
            # })


            selected_db_name <- results$db_info$selected_db()
            message(sprintf("[DB Selected] %s at %s", selected_db_name, Sys.time()))

            withProgress(message = "Stage 1: Listing & filtering tables", value = 0, {
                t0 <- Sys.time()
                message(sprintf("[Stage1] Start at %s", t0))

                DBI::dbExecute(sc, paste0("USE ", selected_db_name))
                tbl_list_query <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", selected_db_name))
                tbls <- tbl_list_query$tableName
                t1 <- Sys.time()
                message(sprintf(
                    "[Stage1] Fetched %d tables at %s (%.2f sec)",
                    length(tbls), t1, as.numeric(difftime(t1, t0, "secs"))
                ))
                setProgress(value = 0.2, detail = sprintf("Fetched %d tables", length(tbls)))

                prefix <- c("^normcounts|^exacttest|^coldata|^expression_table|^gene_table|^sample_table|^gene_metadata|^sample_metadata")
                tbl_list_query_prefix <- tbl_list_query[grepl(paste(prefix, collapse = "|"), tbls), ]
                t2 <- Sys.time()
                message(sprintf(
                    "[Stage1] Prefix filter → %d tables at %s (%.2f sec)",
                    nrow(tbl_list_query_prefix), t2, as.numeric(difftime(t2, t1, "secs"))
                ))
                setProgress(value = 0.4, detail = sprintf("Prefix filter: %d tables", nrow(tbl_list_query_prefix)))

                tbls_with_prefix <- tbl_list_query_prefix$tableName
                tbls_with_time_filter <- get_latest_file_group_df(tbls_with_prefix)
                t3 <- Sys.time()
                message(sprintf(
                    "[Stage1] Time filter applied at %s (%.2f sec)",
                    t3, as.numeric(difftime(t3, t2, "secs"))
                ))
                setProgress(value = 0.6, detail = "Applied time filter")

                if (any(tbls_with_time_filter$is_latest)) {
                    sel_idx <- tbls_with_time_filter$is_latest
                    message(sprintf("[Stage1] Latest tables found at %s", Sys.time()))
                } else {
                    sel_idx <- !tbls_with_time_filter$is_latest
                    message(sprintf("[Stage1] No latest table, using all at %s", Sys.time()))
                }
                tbl_list_query_prefix_time <- tbl_list_query_prefix[sel_idx, ]
                summary_table <- tbls_with_time_filter[sel_idx, ]
                t4 <- Sys.time()
                message(sprintf(
                    "[Stage1] Selected %d tables at %s (%.2f sec)",
                    nrow(tbl_list_query_prefix_time), t4, as.numeric(difftime(t4, t3, "secs"))
                ))
                setProgress(value = 0.8, detail = sprintf("Selected %d tables", nrow(tbl_list_query_prefix_time)))

                tbls_final <- summary_table$file
                normcount_tbls <- tbl_list_query_prefix_time[grepl("^normcounts|^expression_table", tbls_final, ignore.case = TRUE), "tableName"]
                exacttest_tbls <- tbl_list_query_prefix_time[grepl("^exacttest|^gene_table", tbls_final, ignore.case = TRUE), "tableName"]
                coldata_tbls <- tbl_list_query_prefix_time[grepl("^coldata|^sample_table", tbls_final, ignore.case = TRUE), "tableName"]
                t5 <- Sys.time()
                message(sprintf(
                    "[Stage1] Categorized tables at %s (%.2f sec)",
                    t5, as.numeric(difftime(t5, t4, "secs"))
                ))

                setProgress(value = 1, detail = "Stage 1 complete")
                results$table_list <- tbl_list_query_prefix_time
            })



            req(normcount_tbls, exacttest_tbls, coldata_tbls)

            t0_norm_launch <- Sys.time()
            message(sprintf("[Stage2-normcount] Launch at %s", t0_norm_launch))
            normcount_promise <- future_promise(
                {
                    start_time <- Sys.time()
                    message(sprintf("[%s] Start querying normcounts table", start_time))

                    sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
                    on.exit(sparklyr::spark_disconnect(sc_conn))
                    DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
                    query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
                    normcount <- DBI::dbGetQuery(sc_conn, query_normcount)

                    colnames(normcount)[colnames(normcount) == "genes"] <- "GeneSymbol"
                    normcount <- normcount[, colnames(normcount) != "_c0"]

                    end_time <- Sys.time()
                    message(sprintf(
                        "[%s] Completed normcounts query (Duration: %.2f seconds)",
                        end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
                    ))

                    normcount
                },
                globals = list(
                    master = master, method = method, version = version,
                    normcount_tbls = normcount_tbls, selected_db_name = selected_db_name
                ),
                seed = TRUE
            )

            t0_exact_launch <- Sys.time()
            message(sprintf("[Stage2-exacttest] Launch at %s", t0_exact_launch))
            exacttest_promise <- future_promise(
                {
                    start_time <- Sys.time()
                    message(sprintf("[%s] Start querying exacttest table", start_time))

                    sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
                    on.exit(sparklyr::spark_disconnect(sc_conn))
                    DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
                    query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
                    exacttest <- DBI::dbGetQuery(sc_conn, query_exacttest)

                    colnames(exacttest)[colnames(exacttest) == "genes"] <- "GeneSymbol"
                    exacttest <- exacttest[, colnames(exacttest) != "_c0"]

                    end_time <- Sys.time()
                    message(sprintf(
                        "[%s] Completed exacttest query (Duration: %.2f seconds)",
                        end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
                    ))

                    exacttest
                },
                globals = list(
                    master = master, method = method, version = version,
                    exacttest_tbls = exacttest_tbls, selected_db_name = selected_db_name
                ),
                seed = TRUE
            )

            t0_coldata_launch <- Sys.time()
            message(sprintf("[Stage2-coldata] Launch at %s", t0_coldata_launch))
            coldata_promise <- if (length(coldata_tbls) > 0) {
                future_promise(
                    {
                        start_time <- Sys.time()
                        message(sprintf("[%s] Start querying coldata table", start_time))

                        sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
                        on.exit(sparklyr::spark_disconnect(sc_conn))
                        DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
                        query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
                        coldata <- DBI::dbGetQuery(sc_conn, query_coldata)

                        end_time <- Sys.time()
                        message(sprintf(
                            "[%s] Completed coldata query (Duration: %.2f seconds)",
                            end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
                        ))

                        coldata
                    },
                    globals = list(
                        master = master, method = method, version = version,
                        coldata_tbls = coldata_tbls, selected_db_name = selected_db_name
                    ),
                    seed = TRUE
                )
            } else {
                normcount_promise %...>% (function(normcount) {
                    future_promise(
                        {
                            start_time <- Sys.time()
                            message(sprintf("[%s] Generating random coldata", start_time))
                            coldata <- generate_colData_random(normcount, genecol = "GeneSymbol")
                            end_time <- Sys.time()
                            message(sprintf(
                                "[%s] Completed coldata generation (Duration: %.2f seconds)",
                                end_time, as.numeric(difftime(end_time, start_time, units = "secs"))
                            ))
                            coldata
                        },
                        seed = TRUE
                    )
                })
            }
            progressMod$addPromise(normcount_promise, label = "normcount")
            progressMod$addPromise(coldata_promise, label = "coldata")
            progressMod$addPromise(exacttest_promise, label = "exacttest")


            withProgress(message = "Stage 3: Collecting and processing data", value = 0, {
                t0_collect <- Sys.time()
                message(sprintf("[Stage3] Collection start at %s", t0_collect))

                promise_all(
                    normcount_data = normcount_promise,
                    exacttest_data = exacttest_promise,
                    coldata        = coldata_promise
                ) %...>% with({
                    t1_collect <- Sys.time()
                    message(sprintf(
                        "[Stage3] Collection completed at %s (Duration: %.2f seconds)",
                        t1_collect, as.numeric(difftime(t1_collect, t0_collect, units = "secs"))
                    ))

                    results$normcount_data <- normcount_data
                    results$exacttest_data <- exacttest_data
                    results$coldata <- coldata

                    message("=== normcount_data ===")
                    print(head(results$normcount_data))
                    message("=== exacttest_data ===")
                    print(head(results$exacttest_data))
                    message("=== coldata ===")
                    print(head(results$coldata))

                    setProgress(value = 1, detail = "Data collected")
                })
            })

            observe({
                req(results$exacttest_data, results$normcount_data, results$coldata)
                withProgress(message = "Loading & filtering files...", value = 0, {
                    # 1/6: Read expression_matrix.csv
                    incProgress(1 / 6, detail = "Reading expression_matrix.csv")
                    expr_mat <- results$normcount_data

                    # 2/6: Read gene_info.csv
                    incProgress(2 / 6, detail = "Reading gene_info.csv")
                    gene_info <- results$exacttest_data

                    # 3/6: Read sample_info.csv
                    incProgress(3 / 6, detail = "Reading sample_info.csv")
                    sample_info <- results$coldata

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
                    print("test1")
                    results$normcount_data <- expr_mat
                    results$exacttest_data <- gene_info
                    results$coldata <- sample_info
                    print("test2")
                    print(dim(expr_mat))
                    print(dim(gene_info))
                    print(dim(sample_info))
                    print(dim(wide_data()))
                    print(dim(DEG_table()))
                    print(dim(maeColData()))
                    wide_data(expr_mat)
                    DEG_table(gene_info)
                    maeColData(sample_info)
                    print(dim(wide_data()))
                    print(dim(DEG_table()))
                    print(dim(maeColData()))
                    print("test3")
                })
                # DEG_table(results$exacttest_data)
                # wide_data(results$normcount_data)
                # maeColData(results$coldata)
                message("assign reactiveVal: DEG_table, wide_data, maeColData")
            })

            # output$wide_table_dt <- DT::renderDataTable({
            #     req(wide_data())
            #     normCount_round <- as.data.frame(lapply(
            #         wide_data(),
            #         function(x) if (is.numeric(x)) round(x, 4) else x
            #     ))

            #     print("send wide data to UI")
            #     DT::datatable(
            #         normCount_round,
            #         options = list(pageLength = 20, autoWidth = TRUE)
            #     )
            # })

            message(sprintf("[Process] Completed all stages on PID %s at %s", Sys.getpid(), Sys.time()))
        })



        # observeEvent(input$load_data, {
        #     req(results$exacttest_data, results$normcount_data, results$coldata)

        #     withProgress(message = "Loading & filtering files...", value = 0, {
        #         # 1/6: Read expression_matrix.csv
        #         incProgress(1 / 6, detail = "Reading expression_matrix.csv")
        #         expr_mat <- results$normcount_data

        #         # 2/6: Read gene_info.csv
        #         incProgress(2 / 6, detail = "Reading gene_info.csv")
        #         gene_info <- results$exacttest_data

        #         # 3/6: Read sample_info.csv
        #         incProgress(3 / 6, detail = "Reading sample_info.csv")
        #         sample_info <- results$coldata

        #         # 4/6: Check required columns
        #         incProgress(4 / 6, detail = "Checking columns")
        #         if (!"GeneSymbol" %in% colnames(expr_mat)) {
        #             stop("expression_matrix.csv must contain a 'GeneSymbol' column.")
        #         }

        #         # 5/6: Gene filtering → Keep top 5000 genes by variance
        #         incProgress(5 / 6, detail = "Filtering top 5000 genes by variance")
        #         orig_n <- nrow(expr_mat)
        #         if (orig_n > 3000) {
        #             # Only select numeric columns
        #             data_mat <- expr_mat[, setdiff(colnames(expr_mat), "GeneSymbol"), drop = FALSE]
        #             # Calculate variance for each gene (row)
        #             gene_vars <- apply(data_mat, 1, var, na.rm = TRUE)
        #             # Order by variance, keep top 3000
        #             top_idx <- order(gene_vars, decreasing = TRUE)[1:3000]
        #             expr_mat <- expr_mat[top_idx, , drop = FALSE]
        #             showNotification(
        #                 sprintf("Filtered from %d genes to %d (Top 3000 by variance)", orig_n, nrow(expr_mat)),
        #                 type = "message"
        #             )
        #         }

        #         # 6/6: Update Reactive values
        #         incProgress(6 / 6, detail = "Updating reactive values")
        #         print("test1")
        #         results$normcount_data <- expr_mat
        #         results$exacttest_data <- gene_info
        #         results$coldata <- sample_info
        #         print("test2")
        #         print(dim(expr_mat))
        #         print(dim(gene_info))
        #         print(dim(sample_info))
        #         print(dim(wide_data()))
        #         print(dim(DEG_table()))
        #         print(dim(maeColData()))
        #         wide_data(expr_mat)
        #         DEG_table(gene_info)
        #         maeColData(sample_info)
        #         print(dim(wide_data()))
        #         print(dim(DEG_table()))
        #         print(dim(maeColData()))
        #         print("test3")
        #     }) # end withProgress

        #     showNotification("Files loaded and filtered!", type = "message")
        # })

        # 2. Expression preview
        output$wide_table_dt <- renderDT({
            req(wide_data())
            print("render wide_table_dt")
            df <- wide_data()[1:10,1:10]
            print(dim(df))
            # Round numeric columns for display
            #df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, 4) else x)
            datatable(df, options = list(pageLength = 20, autoWidth = TRUE))
        }, server = TRUE)

        # 3. DEG preview and download
        output$DEG_table_view <- renderDT({
                req(DEG_table())
                print("render DEG_table")
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
            server = TRUE
        )

        # output$download_DEG <- downloadHandler(
        #     filename = function() "DEG_table.csv",
        #     content  = function(file) write.csv(DEG_table(), file, row.names = FALSE)
        # )

        # 4. Sample Info preview
        output$sample_info_dt <- renderDT({
            req(maeColData())
            print("render sample_info_dt")
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
            raw_df <- wide_data()[1:100,1:100]
            print(dim(raw_df))
            print(head(raw_df))
            rownames(raw_df) <- raw_df$"GeneSymbol"
            expr_mat <- t(as.matrix(raw_df[, setdiff(colnames(raw_df), "GeneSymbol")]))
            print(dim(expr_mat))
            print(head(expr_mat))

            # Data quality check
            gsg <- goodSamplesGenes(expr_mat, verbose = 1)
            print("gsg$allOK")
            print(gsg$allOK)
            print("gsg$goodSamples")
            print(gsg$goodSamples)
            if (!gsg$allOK) expr_mat <- expr_mat[gsg$goodSamples, gsg$goodGenes]

            # Reactive wrapper for numeric conversion
            exprDataNumeric <- reactive({
                df <- as.data.frame(expr_mat)
                df[] <- lapply(df, as.numeric)
                df
            })
            print(str(exprDataNumeric()))
            # Call each Shiny module
            # 問題出在計算correlation矩陣的時候，會因為資料量太大而導致記憶體不足
            sampleOut <- sparksampleClustServer(
                "sample",
                exprData    = exprDataNumeric,
                cutHeight   = reactive(input$cutHeight)
            )
            # observeEvent(sampleOut$maxHeight(), {
            #     maxH <- sampleOut$maxHeight()
            #     updateSliderInput(session, "cutHeight",
            #         min = 0, max = ceiling(maxH),
            #         step = 1,
            #         value = min(input$cutHeight, ceiling(maxH))
            #     )
            # })

            # sparksftServer(
            #     "sft",
            #     exprData      = sampleOut$filteredExpr,
            #     powerRange    = reactive(input$powerRange),
            #     rsqCut        = reactive(input$rsqCut),
            #     selectedPower = reactive(input$selectedPower)
            # )

            # modulesObj <- sparkgeneModuleServer(
            #     "mod",
            #     exprData   = sampleOut$filteredExpr,
            #     power      = reactive(input$selectedPower),
            #     deepSplit  = reactive(input$deepSplit),
            #     minSize    = reactive(input$minModuleSize),
            #     runTrigger = reactive(input$runModules)
            # )
            # observeEvent(maeColData(), {
            #     df <- maeColData()
            #     cols <- colnames(df)
            #     keys <- cols[sapply(cols, function(cn) length(unique(df[[cn]])) == nrow(df))]
            #     updateSelectInput(session, "heatmap_key",
            #         choices  = keys,
            #         selected = keys[1]
            #     )
            #     updateSelectInput(session, "heatmap_phenotypes",
            #         choices  = setdiff(cols, keys[1]),
            #         selected = character(0)
            #     )
            # })

            # output$moduleHeatmap <- renderPlot({
            #     req(modulesObj(), maeColData(), input$heatmap_key)

            #     me <- modulesObj()$MEs
            #     mat <- t(as.matrix(me))
            #     print(str(mat))
            #     print(input$heatmap_key)
            #     print(class(input$heatmap_key))

            #     samp <- maeColData()


            #     ids <- as.character(samp[[input$heatmap_key]])
            #     print(str(ids))
            #     if (any(duplicated(ids))) ids <- make.unique(ids)
            #     rownames(samp) <- ids
            #     print(str(samp))
            #     anno_col <- NULL
            #     if (length(input$heatmap_phenotypes) > 0) {
            #         tmp <- samp[colnames(mat), input$heatmap_phenotypes, drop = FALSE]
            #         print(str(tmp))
            #         keep <- sapply(tmp, function(x) {
            #             x2 <- x[!is.na(x)]
            #             length(unique(x2)) > 1
            #         })
            #         anno_col <- tmp[, keep, drop = FALSE]
            #         print(str(anno_col))
            #     }

            #     row_anno <- data.frame(Module = colnames(me))
            #     rownames(row_anno) <- colnames(me)

            #     pheatmap(
            #         mat,
            #         annotation_col = anno_col,
            #         annotation_row = row_anno,
            #         cluster_rows   = TRUE,
            #         cluster_cols   = TRUE,
            #         show_rownames  = TRUE,
            #         show_colnames  = TRUE
            #     )
            #     print(dev.list())
            # })

            # geneListServer(
            #     "list",
            #     exprData    = sampleOut$filteredExpr,
            #     modulesObj  = modulesObj
            # )
        })
    }

    # Launch App
    shinyApp(ui = ui, server = server)
}


spark_disconnect(sc)
devtools::load_all("/Users/charleschuang/Desktop/atgenomix/ToolDev/oncoexpr")
source("/Users/charleschuang/Desktop/atgenomix/ToolDev/oncoexpr/R/mod_wgcna_spark.R")

library(sparklyr)
options(shiny.maxRequestSize = 500 * 1024^2)
NebulaCoNet_sparkdev(master = "sc://localhost:15002", method = "shell" )


#k8s
NebulaCoNet_sparkdev(master = "會拿預設的endpoint", method = "shell(預設)" )

#NebulaCoNet_localdev(master = "sc://localhost:15002")


NebulaCoNet(master = "sc://localhost:15002")
sparkCorMatrix

warnings()
