library(data.table)
library(shiny)
library(shinydashboard)
library(visNetwork)
library(igraph)
library(DT)
library(pmd)

# Define UI
ui <- dashboardPage(
    dashboardHeader(title = "Molecular Network Analysis", titleWidth = 350),
    dashboardSidebar(
        sidebarMenu(
            menuItem("File Upload", tabName = "upload", icon = icon("upload")),
            menuItem("Network", tabName = "network", icon = icon("network-wired")),
            menuItem("Usage Guide", tabName = "usage", icon = icon("question-circle"))
        )
    ),
    dashboardBody(
        tabItems(
            # File Upload Tab
            tabItem(tabName = "upload",
                    fluidRow(
                        box(width = 12, title = "Upload Data Files",
                            fileInput("upload_msinetsub", "Upload network csv",
                                      accept = c(".csv"),
                                      buttonLabel = "Browse...",
                                      placeholder = "No file selected"),
                            fileInput("upload_isletannoccs", "Upload annotation csv",
                                      accept = c(".csv"),
                                      buttonLabel = "Browse...",
                                      placeholder = "No file selected"),
                            fileInput("upload_images", "Upload node images ",
                                      multiple = TRUE,
                                      accept = c('.png'),
                                      buttonLabel = "Browse...",
                                      placeholder = "No images selected"),
                            actionButton("process_files", "Process Files", class = "btn-primary")
                        )
                    )
            ),

            # Network Tab
            tabItem(tabName = "network",
                    # Disable network tab until files are processed
                    conditionalPanel(
                        condition = "!input.process_files",
                        div(class = "alert alert-warning",
                            "Please upload and process files first in the 'File Upload' tab.")
                    ),
                    conditionalPanel(
                        condition = "input.process_files",
                        # First row for network and image
                        fluidRow(
                            # Molecular Network Analysis
                            box(width = 8, height = 700,
                                visNetworkOutput("network", height = "550px")
                            ),
                            # Right column with image display
                            box(width = 4, height = 700,
                                # Node dropdown and search
                                selectizeInput(
                                    "search_dropdown",
                                    "Select or Search Node:",
                                    choices = NULL,  # Will be populated dynamically
                                    options = list(
                                        placeholder = 'Type to search nodes',
                                        maxOptions = 20
                                    )
                                ),
                                hr(),
                                uiOutput("image")
                            )
                        ),
                        # Second row for PMD search
                        fluidRow(
                            box(width = 12,
                                # Data search section
                                h4("PMD Search"),
                                fluidRow(
                                    column(3,
                                           numericInput("search_number", "Enter number to search:",
                                                        value = 18.011)
                                    ),
                                    column(3,
                                           numericInput("tolerance", "Tolerance (±):",
                                                        value = 0.001)
                                    ),
                                    column(3,
                                           div(style = "margin-top: 25px;",
                                               actionButton("do_search", "Search Data",
                                                            class = "btn-primary")
                                           )
                                    ),
                                    column(3,
                                           div(style = "margin-top: 25px;",
                                               downloadButton("download_results", "Download Results",
                                                              class = "btn-success")
                                           )
                                    )
                                ),
                                hr(),
                                DTOutput("search_results")
                            )
                        )
                    )
            ),

            # Usage Guide Tab
            tabItem(tabName = "usage",
                    fluidRow(
                        box(width = 12, title = "Application Usage Guide",
                            h1("Molecular Network Analysis Tool"),
                            tags$ul(
                                tags$li(tags$strong("File Upload:"),
                                        "Upload your network, annatation and images files to start analysis."),
                                tags$li(tags$strong("Network Visualization:"),
                                        "Explore the molecular network by interacting with nodes and edges."),
                                tags$li(tags$strong("Node Search:"),
                                        "Use the dropdown to select and highlight nodes. If no image is shown, it means one annotation is linked with multiple peaks in the data."),
                                tags$li(tags$strong("PMD Search:"),
                                        "Search for PMD and related reactions by entering a PMD and specifying a tolerance range.")
                            ),
                            h4("Interaction Features:"),
                            tags$ul(
                                tags$li("Click on nodes to select and view MSI on the right panel"),
                                tags$li("Hover over edges to see PMD and spatial correlation information"),
                                tags$li("Use the node label dropdown to filter network view")
                            ),
                            p("For more information, please contact <miao.yu@jax.org>.")
                        )
                    )
            )
        )
    )
)

# Define server
server <- function(input, output, session) {
    # Reactive values to store uploaded data
    uploaded_data <- reactiveValues(
        msinetsub = NULL,
        isletannoccs = NULL,
        network = NULL,
        search_data = NULL,
        nodes = NULL,
        edges = NULL
    )

    # Process uploaded files
    observeEvent(input$process_files, {
        req(input$upload_msinetsub, input$upload_isletannoccs)

        # Load uploaded files
        uploaded_data$msinetsub <- fread(input$upload_msinetsub$datapath)
        uploaded_data$isletannoccs <- fread(input$upload_isletannoccs$datapath)

        # Use uploaded data to create network
        net <- graph_from_data_frame(uploaded_data$msinetsub, directed = F)
        V(net)$comp <- uploaded_data$isletannoccs$name[match(V(net)$name,
                                                             paste0(round(uploaded_data$isletannoccs$mz,4),'_',
                                                                    round(uploaded_data$isletannoccs$im)))]
        V(net)$comp <- iconv(V(net)$comp, to = "UTF-8")

        # Prepare network data
        uploaded_data$nodes <- data.frame(
            id = V(net)$name,
            label = V(net)$comp,
            stringsAsFactors = FALSE
        )

        uploaded_data$edges <- data.frame(
            from = as.character(ends(net, E(net))[, 1]),
            to = as.character(ends(net, E(net))[, 2]),
            title = paste('PMD:', E(net)$diff2, 'cor:', round(E(net)$cor,2))
        )

        # Prepare search data (assuming a column 'pmd' exists)
        uploaded_data$search_data <- as.data.frame(uploaded_data$msinetsub)

        # Handle image uploads
        if (!is.null(input$upload_images)) {
            # Create images directory if it doesn't exist
            dir.create("images", showWarnings = FALSE)

            # Copy uploaded images to fig directory
            for (img in input$upload_images$datapath) {
                file.copy(img, file.path("images", input$upload_images$name[input$upload_images$datapath == img]))
            }
        }

        # Notify user of successful processing
        showNotification("Files processed successfully!", type = "message")
    })

    # Download handler for search results
    output$download_results <- downloadHandler(
        filename = function() {
            paste("pmd_search_results_",
                  format(Sys.time(), "%Y%m%d_%H%M%S"),
                  ".csv",
                  sep = "")
        },
        content = function(file) {
            req(uploaded_data$search_data)
            search_value <- input$search_number
            tolerance <- input$tolerance

            # Filter data within tolerance range
            results <- uploaded_data$search_data[abs(uploaded_data$search_data$pmd - search_value) <= tolerance,]

            # Write the filtered data to the file
            write.csv(results, file, row.names = FALSE)
        }
    )

    # Add resource path for images
    addResourcePath("images", "images")

    # Dynamic update of dropdown choices
    observe({
        req(uploaded_data$nodes)

        # Combine unique node IDs and labels, create unique choices
        choices <- unique(c(
            setNames(uploaded_data$nodes$id, uploaded_data$nodes$label),  # Use labels as display, IDs as values
            setNames(uploaded_data$nodes$label, uploaded_data$nodes$label)
        ))

        updateSelectizeInput(
            session,
            "search_dropdown",
            choices = choices,
            server = TRUE
        )
    })

    # Network visualization
    output$network <- renderVisNetwork({
        req(uploaded_data$nodes, uploaded_data$edges)

        visNetwork(uploaded_data$nodes, uploaded_data$edges) %>%
            visNodes(shape = "dot", size = 30) %>%
            visEdges(font = list(align = "middle")) %>%
            visOptions(
                highlightNearest = TRUE,
                selectedBy = list(variable = "label")
            ) %>%
            visEvents(selectNode = "function(nodes) {
                  Shiny.onInputChange('selected_node', nodes);
                }") %>%
            visLayout(randomSeed = 123)
    })

    # Node selection handler
    observeEvent(input$selected_node, {
        req(input$selected_node)
        node_id <- input$selected_node$nodes[[1]]
        show_node_image(node_id)
    })

    # Dropdown search handler
    observeEvent(input$search_dropdown, {
        req(uploaded_data$nodes)

        if (input$search_dropdown != "") {
            # Find the corresponding node ID
            match_id <- uploaded_data$nodes$id[uploaded_data$nodes$label %in% input$search_dropdown]
            if (length(match_id) == 0) {
                match_id <- input$search_dropdown
            }

            # Select and highlight the node
            visNetworkProxy("network") %>%
                visSelectNodes(id = match_id)

            # Show the node's image
            show_node_image(match_id)
        }
    })

    # Function to show node image with more flexible matching
    show_node_image <- function(node_id) {
        # Try exact match first
        image_path <- paste0("images/", node_id, ".png")

        # If exact match fails, try partial matching
        if (!file.exists(image_path)) {
            # List all image files
            all_images <- list.files("images", pattern = "\\.png$")

            # Try to find images that contain the node_id
            matching_images <- all_images[grepl(gsub("[^0-9.]", "", node_id), all_images)]

            if (length(matching_images) > 0) {
                image_path <- paste0("images/", matching_images[1])
            }
        }

        if (file.exists(image_path)) {
            output$image <- renderUI({
                tags$div(
                    style = "width: 100%;",
                    tags$div(
                        style = "width: 100%;
                   height: 30px;
                   margin-bottom: 10px;
                   background: linear-gradient(to right, lightblue, red);",
                        tags$div(
                            style = "display: flex;
                     justify-content: space-between;
                     font-size: 12px;",
                            tags$span("Low"),
                            tags$span("High")
                        )
                    ),
                    tags$img(
                        src = image_path,
                        alt = paste("Image for node", node_id),
                        style = "width: 100%; height: auto; display: block;"
                    )
                )
            })
        } else {
            output$image <- renderUI({
                div(class = "alert alert-warning",
                    paste("No image found for node:", node_id))
            })
        }
    }

    # Data search functionality
    observeEvent(input$do_search, {
        req(uploaded_data$search_data)

        search_value <- input$search_number
        tolerance <- input$tolerance

        # Filter data within tolerance range
        results <- uploaded_data$search_data[abs(uploaded_data$search_data$pmd - search_value) <= tolerance,]

        output$search_results <- renderDT({
            datatable(results,
                      options = list(
                          pageLength = 10,
                          scrollX = TRUE,
                          dom = 'frtip',
                          scrollY = "400px"
                      ),
                      style = 'bootstrap'
            )
        })
    })
}

shinyApp(ui = ui, server = server)
