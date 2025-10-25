# MZCCS Annotation Tool - app.R
# This single-file Shiny application is designed for hosting.
# It converts the original R Markdown flexdashboard into a standard app structure.
# NOTE: For this app to work, ensure that 'lipidall.csv' and 'metaall.csv'
# are in the same directory as this 'app.R' file.

# 1. LOAD LIBRARIES
# ==============================================================================
# Ensure these packages are installed: install.packages(c("shiny", "DT", "data.table", "enviGCMS", "shinyjs"))
library(shiny)
library(DT)
library(data.table)
library(enviGCMS)
library(shinyjs) # Using shinyjs to conditionally show/hide UI elements

# 2. DEFINE USER INTERFACE (UI)
# ==============================================================================
# The UI is defined using fluidPage and sidebarLayout, which provides a standard
# two-column layout with inputs on the left and outputs on the right.
ui <- fluidPage(
  useShinyjs(), # Initialize shinyjs
  # Application title
  titlePanel("MZCCS Annotation Tool"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    # -- Sidebar panel for inputs --
    sidebarPanel(
      width = 3, # Adjust sidebar width
      h4("1. Select/Upload Data"),
      
      # --- Database Selection ---
      radioButtons("db_source", "Choose Database Source:",
                   choices = c("Use pre-loaded lipids database" = "lipid",
                               "Use pre-loaded metabolites database" = "metabolite",
                               "Upload my own file" = "upload"),
                   selected = "lipid"),
      
      conditionalPanel(
        condition = "input.db_source == 'upload'",
        fileInput("database", "Upload Database CSV", accept = ".csv")
      ),
      
      hr(),
      
      # --- Reference Selection ---
      radioButtons("ref_source", "Choose Experimental Data:",
                   choices = c("Use demo reference file" = "demo",
                               "Upload my own file" = "upload_ref"),
                   selected = "demo"),
      
      # Button to download the demo reference file
      downloadButton("download_demo_ref", "Download Demo Data", class = "btn-info btn-sm"),
      
      br(),br(),
      
      conditionalPanel(
        condition = "input.ref_source == 'upload_ref'",
        fileInput("reference", "Upload Reference CSV", accept = ".csv")
      ),
      
      hr(),
      h4("2. Set Parameters"),
      selectInput("ion_mode", "Ion Mode",
                  choices = c("Positive" = "positive", "Negative" = "negative"),
                  selected = "positive"),
      
      numericInput("mz_ppmshift", "M/Z PPM Shift", value = 20, min = 0),
      numericInput("ccs_shift", "CCS Shift (%)", value = 5, min = 0),
      
      hr(),
      h4("3. Run and Download"),
      # Action button to trigger the analysis
      actionButton("run_analysis", "Run Analysis", icon = icon("play"), class = "btn-primary"),
      
      br(),
      br(),
      
      # Download button for the results
      downloadButton("download_results", "Download Results", class = "btn-success")
    ),
    
    # -- Main panel for displaying outputs --
    mainPanel(
      width = 9, # Adjust main panel width
      h3("Annotation Results"),
      p("Results will appear here after clicking 'Run Analysis'."),
      # Output: Data table will be rendered here
      DT::dataTableOutput("results_table")
    )),
    hr(),
    div(
      style = "text-align: center; padding: 15px; font-size: 0.9em; color: #6c757d;",
      p(strong("Developer:"), " Miao Yu"),
      p(strong("Published for:"), " The SenNet DAC at The Jackson Laboratory (U54AG079753)"),
      p(strong("Last Updated:"), " June 17, 2025")
    )
)


# 3. DEFINE SERVER LOGIC
# ==============================================================================
# The server function contains the logic to react to user inputs and
# generate the outputs.
server <- function(input, output, session) {
  
  # --- Create Demo Reference Data ---
  # This data.table is created when the app starts and is used when the demo option is selected.
  demo_ref_data <- data.table(
    mz = c(597.4219, 532.1691, 649.4510, 484.2722, 425.3431),
    im = c(268.6109, 224.3536, 273.0371, 220.5825, 224.9311)
  )
  
  # --- Reactive expression to load the correct database ---
  # This reactive will return the data from the chosen source.
  database_data <- reactive({
    source <- input$db_source
    if (source == "lipid") {
      path <- "lipidall.csv"
    } else if (source == "metabolite") {
      path <- "metaall.csv"
    } else if (source == "upload") {
      req(input$database$datapath)
      return(fread(input$database$datapath))
    } else {
      return(NULL)
    }
    
    if (!file.exists(path)) {
      showNotification(paste("Error: Pre-loaded file '", path, "' not found in the app directory."), type = "error", duration = 10)
      return(NULL)
    }
    fread(path)
  })
  
  # --- Reactive expression to load the correct reference data ---
  reference_data <- reactive({
    source <- input$ref_source
    if (source == "demo") {
      demo_ref_data
    } else if (source == "upload_ref") {
      req(input$reference$datapath)
      fread(input$reference$datapath)
    } else {
      NULL
    }
  })
  
  # --- Reactive analysis results ---
  # Use eventReactive to ensure the analysis code runs ONLY when
  # the 'run_analysis' button is clicked.
  analysis_results <- eventReactive(input$run_analysis, {
    
    # CORRECTED: Require both database and reference data to be available from the reactive expressions.
    req(database_data(), reference_data())
    
    showNotification("Running analysis...", type = "message", duration = 5)
    
    tryCatch({
      lipid <- database_data()
      ref <- reference_data()
      
      if (input$ion_mode == "positive") {
        lipid <- lipid[adducts %in% c('[M+H]','[M+Na]','[M+NH4]','[M+H-H2O]','[M+]','[M+K]')]
      } else {
        lipid <- lipid[adducts %in% c('[M-H]','[M+Cl]','[M-H+FA]','[M-]','[M+Na-2H]')]
      }
      
      mz <- ref$mz
      im <- ref$im
      
      align <- enviGCMS::getalign(mz, lipid$mz, im, lipid$ccs,
                                  ppm = input$mz_ppmshift, deltart = input$ccs_shift)
      
      if (is.null(align) || nrow(align) == 0) {
        showNotification("No matching features found with the given parameters.", type = "warning", duration = 8)
        return(NULL)
      }
      
      anno <- cbind.data.frame(mz = mz[align$xid], im = im[align$xid],
                               db_mz = align$mz2, db_im = align$rt2)
      
      lipidanno <- merge(anno, lipid, by.x = c('db_mz','db_im'), by.y = c('mz','ccs'))
      lipidanno <- lipidanno[!duplicated(lipidanno),]
      
      showNotification("Analysis complete! Displaying results.", type = "message", duration = 5)
      
      return(lipidanno)
      
    }, error = function(e) {
      showNotification(paste("An error occurred:", e$message), type = "error", duration = 10)
      return(NULL)
    })
  })
  
  # --- Render the results table ---
  output$results_table <- DT::renderDataTable({
    req(analysis_results())
    
    DT::datatable(
      analysis_results(),
      options = list(
        pageLength = 15, scrollX = TRUE, dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE, filter = 'top', extensions = 'Buttons'
    )
  })
  
  # --- Handle demo reference file download ---
  output$download_demo_ref <- downloadHandler(
    filename = function() {
      "demo_reference_template.csv"
    },
    content = function(file) {
      fwrite(demo_ref_data, file)
    }
  )
  
  # --- Handle analysis results file download ---
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("mzccs_annotation_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(analysis_results())
      fwrite(analysis_results(), file)
    }
  )
  
}


# 4. RUN THE APPLICATION
# ==============================================================================
# This line brings the UI and Server together to launch the Shiny app.
shinyApp(ui = ui, server = server)
