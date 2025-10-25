# --- 1. Load Libraries ---
library(data.table)
library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(DT)
library(pmd)

# --- 2. Load Static Data (if any) ---
# These are used to annotate the user's uploaded data
data("keggrall")
data("sda")
search_data_kegg <- as.data.frame(keggrall)

# --- 3. Define UI ---
ui <- dashboardPage(
    dashboardHeader(title = "Molecular Network Analysis", titleWidth = 350),
    dashboardSidebar(
        sidebarMenu(
            id = "tabs",
            menuItem("1. File Upload", tabName = "upload", icon = icon("upload")),
            menuItem("2. Network Analysis", tabName = "analysis", icon = icon("project-diagram"))
        )
    ),
    dashboardBody(
        tabItems(
            # -- Upload Tab --
            tabItem(tabName = "upload",
                    fluidRow(
                        box(width = 12, title = "Upload Your Data", status = "primary", solidHeader = TRUE,
                            p("Please provide the required files to build the molecular network. All files must be in CSV format, except for the images."),
                            hr(),
                            tags$h4("File Format Requirements:"),
                            tags$ul(
                                tags$li(tags$strong("Network Edges CSV:"), "Must contain the following columns: ", tags$code("from"), ", ", tags$code("to"), ", ", tags$code("diff"), ", and ", tags$code("cor"), ". The 'from' and 'to' columns should contain node IDs in the format 'mz_im'."),
                                tags$li(tags$strong("Node Annotations CSV:"), "Must contain the columns: ", tags$code("mz"), ", ", tags$code("im"), ", and ", tags$code("name"), ". These are used to map node IDs to metabolite names."),
                                tags$li(tags$strong("Node Images (.png):"), "Image files must be in PNG format and named after their corresponding node ID (e.g., ", tags$code("100.0765_1.135.png"), ").")
                            ),
                            hr(),
                            fileInput("network_file", "Upload Network Edges CSV", accept = ".csv"),
                            fileInput("annotation_file", "Upload Node Annotations CSV", accept = ".csv"),
                            fileInput("image_files", "Upload Node Images (.png)", accept = ".png", multiple = TRUE),
                            actionButton("process_data", "Process Uploaded Data", icon = icon("cogs"), class = "btn-success"),
                            uiOutput("upload_status")
                        )
                    )
            ),

            # -- Analysis Tab --
            tabItem(tabName = "analysis",
                    # This panel shows only after data is processed successfully
                    conditionalPanel(
                        condition = "output.data_ready",
                        fluidRow(
                            box(width = 12, title = "Molecular Network Analysis Tool", status = "primary", solidHeader = TRUE,
                                tags$ul(
                                    tags$li(tags$strong("Edge Filtering:"),
                                            "Select reactions (PMDs) from the 'Reactions (PMDs) Table' below and click 'Show Network' to display the filtered network."),
                                    tags$li(tags$strong("Network Interaction:"),
                                            "Click on a node (metabolite) to see its corresponding image. Hover over an edge to see PMD and correlation details."),
                                    tags$li(tags$strong("Metabolite Search:"),
                                            "Use the dropdown on the right to find and highlight a specific metabolite in the network.")
                                )
                            )
                        ),

                        fluidRow(
                            box(width = 12, title = "Data Exploration",
                                tabsetPanel(
                                    tabPanel("Reactions (PMDs) Table", DTOutput("pmd_table")),
                                    tabPanel("Sample PMD Search",
                                             fluidRow(
                                                 column(3, numericInput("search_number2", "Enter PMD:", value = 18.011)),
                                                 column(3, numericInput("tolerance2", "Tolerance (±):", value = 0.001)),
                                                 column(3, div(style = "margin-top: 25px;", actionButton("do_search2", "Search Sample Data", class = "btn-primary"))),
                                                 column(3, div(style = "margin-top: 25px;", downloadButton("download_results2", "Download Results", class = "btn-success")))
                                             ),
                                             hr(),
                                             DTOutput("search_results2")
                                    ),
                                    tabPanel("KEGG PMD Search",
                                             fluidRow(
                                                 column(3, numericInput("search_number_kegg", "Enter PMD:", value = 18.011)),
                                                 column(3, numericInput("tolerance_kegg", "Tolerance (±):", value = 0.001)),
                                                 column(3, div(style = "margin-top: 25px;", actionButton("do_search_kegg", "Search KEGG Data", class = "btn-primary"))),
                                                 column(3, div(style = "margin-top: 25px;", downloadButton("download_results_kegg", "Download Results", class = "btn-success")))
                                             ),
                                             hr(),
                                             DTOutput("search_results_kegg")
                                    )
                                )
                            )
                        ),

                        fluidRow(
                            box(width = 8, height = 700, status = "info",
                                actionButton("show_network", "Show Network based on PMD Selection", class = "btn-primary", width = "100%", icon = icon("eye")),
                                hr(),
                                visNetworkOutput("network", height = "550px")
                            ),
                            box(width = 4, height = 700, status = "info",
                                selectizeInput("search_dropdown", "Select or Search Node (Metabolite):", choices = NULL),
                                hr(),
                                uiOutput("image")
                            )
                        )
                    ),
                    # This panel shows if data is NOT ready
                    conditionalPanel(
                        condition = "!output.data_ready",
                        fluidRow(
                            box(width=12, status="warning", title="Pending Data",
                                h3("Please go to the 'File Upload' tab to upload and process your data first.")
                            )
                        )
                    )
            )
        )
    )
)

# --- 4. Define Server Logic ---
server <- function(input, output, session) {

    # Reactive values to store all session data
    data_store <- reactiveValues()
    filtered_edges <- reactiveVal(NULL)
    data_is_ready <- reactiveVal(FALSE)

    # This output is used by conditionalPanel to check if data is ready
    output$data_ready <- reactive({
        data_is_ready()
    })
    outputOptions(output, "data_ready", suspendWhenHidden = FALSE)

    # -- File Processing Logic --
    observeEvent(input$process_data, {
        req(input$network_file, input$annotation_file)

        withProgress(message = 'Processing your data...', value = 0, {

            # Step 1: Read uploaded data
            incProgress(0.1, detail = "Reading CSV files...")
            data_store$df <- fread(input$network_file$datapath)
            data_store$anno <- fread(input$annotation_file$datapath)

            # Step 2: Prepare data (same as original script)
            incProgress(0.2, detail = "Annotating edges...")
            data_store$df$pmd3 <- round(data_store$df$diff, 3)
            data_store$df$ms1anno <- data_store$anno$name[match(data_store$df$from, paste0(round(data_store$anno$mz,4),'_',round(data_store$anno$im)))]
            data_store$df$ms2anno <- data_store$anno$name[match(data_store$df$to, paste0(round(data_store$anno$mz,4),'_',round(data_store$anno$im)))]

            data_store$search_data2 <- as.data.frame(data_store$df[,c('from','to','diff','cor','ms1anno','ms2anno')])

            # Step 3: Create the PMD summary table
            incProgress(0.4, detail = "Summarizing PMDs...")
            dfdiff <- cbind.data.frame(pmd=as.numeric(names(table(round(data_store$df$diff,3)))),freq=as.numeric(table(round(data_store$df$diff,3))))
            sda$PMD <- round(sda$PMD,3)
            tablepmd <- merge(dfdiff, sda, by.x = 'pmd', by.y = 'PMD', all.x = T)
            tablepmd <- tablepmd[,-5]
            keggsub <- cbind.data.frame(pmd=as.numeric(names(table(round(keggrall$pmd,3)))),freq=as.numeric(table(round(keggrall$pmd,3))))
            tablepmd$freq_kegg <- keggsub$freq[match(tablepmd$pmd,keggsub$pmd)]
            data_store$tablepmd <- tablepmd

            # Step 4: Create network graph and dataframes
            incProgress(0.6, detail = "Building network graph...")
            net <- graph_from_data_frame(data_store$df, directed = F)

            # Create a vector of annotations; this may contain NAs
            node_annotations <- data_store$anno$name[match(V(net)$name, paste0(round(data_store$anno$mz,4),'_',round(data_store$anno$im)))]
            final_labels <- ifelse(is.na(node_annotations), V(net)$name, node_annotations)

            data_store$nodes <- data.frame(
                id = V(net)$name,
                label = iconv(final_labels, to = "UTF-8"), # Ensure proper encoding
                stringsAsFactors = FALSE
            )

            data_store$edges <- data.frame(
                from = as.character(ends(net, E(net))[, 1]),
                to = as.character(ends(net, E(net))[, 2]),
                title = paste('PMD:', round(E(net)$diff, 4), 'cor:', round(E(net)$cor,2)),
                diff = E(net)$diff,
                cor = E(net)$cor,
                width = abs(E(net)$cor) * 5 # Edge width proportional to correlation
            )

            # Step 5: Handle image uploads
            incProgress(0.8, detail = "Processing images...")
            if (!is.null(input$image_files)) {
                session_img_dir <- file.path(tempdir(), session$token)
                dir.create(session_img_dir, showWarnings = FALSE)

                file.copy(input$image_files$datapath, file.path(session_img_dir, input$image_files$name))

                addResourcePath("session_images", session_img_dir)
                data_store$image_path_prefix <- "session_images"

                session$onSessionEnded(function() {
                    unlink(session_img_dir, recursive = TRUE)
                })
            }

            incProgress(1, detail = "Done!")
        }) # End withProgress

        output$upload_status <- renderUI({
            tags$div(class = "alert alert-success", style = "margin-top: 15px;",
                     "Data processed successfully! You can now proceed to the 'Network Analysis' tab.")
        })

        data_is_ready(TRUE)
        updateTabItems(session, "tabs", "analysis") # Switch to analysis tab
    })

    # -- Analysis Tab Logic --

    # Render the main PMD table
    output$pmd_table <- renderDT({
        req(data_store$tablepmd)
        datatable(data_store$tablepmd, selection = 'multiple', options = list(pageLength = 10, scrollY = "200px"))
    })

    # Populate the node search dropdown
    observe({
        req(data_store$nodes)
        # Create choices that allow searching by metabolite name (label) or by m/z_im (id), like in app.r
        choices <- unique(c(
            setNames(data_store$nodes$id, data_store$nodes$label),
            setNames(data_store$nodes$id, data_store$nodes$id)
        ))
        # Clean up choices to prevent issues with NAs or empty strings
        choices <- choices[!is.na(names(choices))]
        choices <- choices[names(choices) != ""]

        updateSelectizeInput(session, "search_dropdown",
                             choices = choices,
                             server = TRUE,
                             options = list(
                                 placeholder = 'Type to search by name or m/z_im',
                                 maxOptions = 20
                             ))
    })

    # "Show Network" button logic
    observeEvent(input$show_network, {
        selected_rows <- input$pmd_table_rows_selected
        if (is.null(selected_rows) || length(selected_rows) == 0) {
            showModal(modalDialog(title = "Selection Missing", "Please select at least one PMD row from the table.", easyClose = TRUE))
            return()
        }

        selected_pmds <- data_store$tablepmd$pmd[selected_rows]

        filtered <- data_store$edges[round(data_store$edges$diff, 3) %in% round(selected_pmds, 3), ]

        if (nrow(filtered) == 0) {
            showModal(modalDialog(title = "No Data", "No network edges found for the selected PMD values.", easyClose = TRUE))
            return()
        }
        filtered_edges(filtered)
    })

    # Render the network visualization
    output$network <- renderVisNetwork({
        req(filtered_edges())

        edges_to_show <- filtered_edges()
        related_node_ids <- unique(c(edges_to_show$from, edges_to_show$to))
        nodes_to_show <- data_store$nodes[data_store$nodes$id %in% related_node_ids, ]

        visNetwork(nodes_to_show, edges_to_show) %>%
            visNodes(shape = "dot", size = 20) %>%
            visOptions(highlightNearest = TRUE, selectedBy = list(variable = "label")) %>%
            visEvents(selectNode = "function(nodes) { Shiny.onInputChange('selected_node', nodes); }") %>%
            visLayout(randomSeed = 123) %>%
            visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
    })

    # Show image when a node is clicked in the network
    observeEvent(input$selected_node, {
        node_id <- input$selected_node$nodes[[1]]
        if (!is.null(node_id)) {
            show_node_image(node_id)
        }
    })

    # Show image and select node when chosen from dropdown
    observeEvent(input$search_dropdown, {
        req(input$search_dropdown, data_store$nodes)

        # First, try to find the ID based on the label (metabolite name)
        match_id <- data_store$nodes$id[data_store$nodes$label == input$search_dropdown]

        # If no match is found by label, assume the input is already an ID
        if (length(match_id) == 0) {
            match_id <- input$search_dropdown
        }

        # Ensure we have a single, valid ID
        if (length(match_id) == 1 && !is.na(match_id)) {
            # If a network is visible, check if the node is in it and select it
            if (!is.null(filtered_edges())) {
                current_nodes <- unique(c(filtered_edges()$from, filtered_edges()$to))
                if (match_id %in% current_nodes) {
                    visNetworkProxy("network") %>% visSelectNodes(id = match_id)
                } else {
                    showNotification("Node not in current network view.", type = "warning")
                }
            }
            # Always try to show the image for the selected node
            show_node_image(match_id)
        } else {
            showNotification("Could not find the specified node.", type = "error")
        }
    })

    # **ERROR FIX HERE**: Updated function to gracefully handle missing images.
    show_node_image <- function(node_id) {
        # Check if images were uploaded in the first place
        if (is.null(data_store$image_path_prefix)) {
            output$image <- renderUI({
                tags$div(class = "alert alert-warning", "No images were uploaded.")
            })
            return()
        }

        # Check if the specific image file exists
        image_name <- paste0(node_id, ".png")
        image_disk_path <- file.path(tempdir(), session$token, image_name)

        if (file.exists(image_disk_path)) {
            image_src_path <- file.path(data_store$image_path_prefix, image_name)
            output$image <- renderUI({
                tags$div(
                    tags$h4(paste("MSI for:", node_id)),
                    tags$img(src = image_src_path, style = "width: 100%; height: auto;")
                )
            })
        } else {
            # Handle case where a specific image is not found
            output$image <- renderUI({
                tags$div(class = "alert alert-danger", paste("Image not found for node:", node_id))
            })
        }
    }

    # -- Search Logic --
    # Sample Data Search
    observeEvent(input$do_search2, {
        req(data_store$search_data2)
        results <- data_store$search_data2[abs(data_store$search_data2$diff - input$search_number2) <= input$tolerance2, ]
        output$search_results2 <- renderDT(datatable(results, options = list(pageLength = 5, scrollX=TRUE)))
    })

    # KEGG Data Search
    observeEvent(input$do_search_kegg, {
        results <- search_data_kegg[abs(search_data_kegg$pmd - input$search_number_kegg) <= input$tolerance_kegg, ]
        output$search_results_kegg <- renderDT(datatable(results, options = list(pageLength = 5, scrollX=TRUE)))
    })

    # Download Handlers
    output$download_results2 <- downloadHandler(
        filename = function() {"sample_pmd_search_results.csv"},
        content = function(file) {
            results <- data_store$search_data2[abs(data_store$search_data2$diff - input$search_number2) <= input$tolerance2, ]
            write.csv(results, file, row.names = FALSE)
        }
    )

    output$download_results_kegg <- downloadHandler(
        filename = function() {"kegg_pmd_search_results.csv"},
        content = function(file) {
            results <- search_data_kegg[abs(search_data_kegg$pmd - input$search_number_kegg) <= input$tolerance_kegg, ]
            write.csv(results, file, row.names = FALSE)
        }
    )
}

shinyApp(ui = ui, server = server)
