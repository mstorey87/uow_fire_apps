# app.R ---------------------------------------------------------------------
# Load Required Libraries----------------------------------------------------
library(sf)
library(tidyverse)
library(shiny)
library(shinyFiles)
library(shinyWidgets)
library(leaflet)
library(mapedit)
library(wfprogression)
library(bsplus)
library(aws.s3)

Sys.setenv(
  AWS_ACCESS_KEY_ID = ,
  AWS_SECRET_ACCESS_KEY = ,
  AWS_DEFAULT_REGION = "ap-southeast-2"   # example region
)

# Sys.setenv(
#   AWS_ACCESS_KEY_ID = ,
#   AWS_SECRET_ACCESS_KEY = ,
#   AWS_DEFAULT_REGION = "ap-southeast-2"   # example region
# )

# Create Base Leaflet Map----------------------------------------------------
map <- leaflet() %>%
  addTiles() %>%
  setView(lng = 150.3119, lat = -33.7128, zoom = 9)

options(timeout = 1000)

# Define Shiny UI------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Search Fire Images and Progressions"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Search dates (YYYY-dd-mm):",
         icon("question-circle") %>% bs_embed_popover("
         Type in the first date in YYYY-dd-mm format and select number of days.
         Number of days is limited at the moment while testing the app.")
      ),
      
      textInput("fromdate", "Start date", value = '2019-12-01'),
      selectInput("ndays", "Number of days to search", choices = 1:2),
      
      tags$hr(style = "border-top: 1px solid black; margin: 3px 0;"),
      h4("2. Sources",
         icon("question-circle") %>% bs_embed_popover("
         Tick the boxes of image sources to search. 
         Images will be cropped to the search polygon.
         Only low res line scans are availabe, full res could be added later.
         MODIS and VIIRS are screen captures from NASA worldview with hotspots overlayed for time.")
      ),
      
      checkboxInput("progscheck", "Progression polygons", FALSE),
      checkboxInput("imagescheck", "Line scans", FALSE),
      checkboxInput("GAstac", "Sentinel 2 and Landsat (STAC)", FALSE),
      checkboxInput("gibs", "MODIS and VIIRS (NASA Worldview)", FALSE),
      
      h5("Himawari imagery",style = "font-weight: bold;"),
      checkboxInput("H89", "Himawari (false colour RGB)", FALSE),
      checkboxInput("H89pyro", "Himawari (band 14)", FALSE),
      selectInput("him_minutes", "Himawari timestep (minutes)",
                  choices = c(10,20,30,seq(60, 1440, by = 60)), selected = 360),
      
      tags$hr(style = "border-top: 2px solid black; margin: 5px 0;"),
      
      h4("3. Search",
         icon("question-circle") %>% bs_embed_popover("
        Run the query. 
        Box unchecked = Search only (fast).
        Box checked = Prepare a zip for download")),
      # switchInput(
      #   inputId = "gather",
      #   label = "SELECT",
      #   value = FALSE,
      #   onLabel = "Prepare ZIP",
      #   offLabel = "Search only",
      #   onStatus = "info",
      #   offStatus = "dark"
      # ),
      checkboxInput("gather", "Prepare ZIP for download", FALSE),
      
      actionButton('runquery', 'Run', class = "action-btn",width="100%",
                   style = "font-weight: bold;")
    ),
    
    mainPanel(
      tags$div(style = "font-size: 18px; font-weight: bold;",
               "4. Draw a polygon search area"),
      br(),
      editModUI("map"),
      br(),
      fileInput("polygon_file",
                "Or choose a polygon (optional, .gpkg or zipped .shp)",
                accept = c(".gpkg", ".zip")),
      #br(),
      
      # ---- SEARCH RESULTS TEXT OUTPUT ----
      fluidRow(column(12, uiOutput("display_text"))),
      br(),
      
      # ---- DOWNLOAD BUTTON SHOWN ONLY WHEN ZIP READY ----
      conditionalPanel(
        condition = "output.zipAvailable == true",
        downloadButton("downloadScans", "Download ZIP")
      )
    )
  )
)

# Define Shiny Server--------------------------------------------------------
server <- function(input, output, session) {
  
  # mapedit module
  edits <- callModule(editMod, leafmap = map, id = "map")
  
  # reactive containers
  text_list <- reactiveVal(list("Search results:"))
  zip_file_path <- reactiveVal(NULL)
  
  observeEvent(input$runquery, {
    # Reset messages and zip availability
    text_list(list("Search results:"))
    zip_file_path(NULL)
    
    # ensure a clean working output directory
    dir.out <- tempdir()
    
    # remove previous app-written files in dir.out (only common extensions we use)
    old_files <- list.files(dir.out, full.names = TRUE,
                            pattern = "\\.(gpkg|tif|shp|zip|png|jpg)$")
    if (length(old_files) > 0) file.remove(old_files)
    
    # determine polygon to search (uploaded or drawn)
    my_bbox <- NULL
    
    if (!is.null(input$polygon_file)) {
      ext <- tools::file_ext(input$polygon_file$name)
      
      if (ext == "zip") {
        unzip_dir <- tempdir()
        unzip(input$polygon_file$datapath, exdir = unzip_dir)
        shp_files <- list.files(unzip_dir, pattern = "\\.shp$", full.names = TRUE)
        
        if (length(shp_files) > 0) {
          my_bbox <- st_read(shp_files[1], quiet = TRUE) %>%
            st_union() %>% st_as_sf() %>%
            st_transform(3112) %>% st_buffer(5000) %>%
            st_transform(4326)
        } else {
          showNotification("No .shp file found in ZIP", type = "error")
          return()
        }
        
      } else if (ext == "gpkg") {
        my_bbox <- st_read(input$polygon_file$datapath, quiet = TRUE) %>%
          st_union() %>% st_as_sf() %>%
          st_transform(3112) %>% st_buffer(5000) %>%
          st_transform(4326)
      } else {
        showNotification("Unsupported file type. Upload a .gpkg or zipped .shp", type = "error")
        return()
      }
      
    } else {
      # use drawn geometry if available
      geomDraw <- edits()$finished
      if (!is.null(geomDraw)) {
        my_bbox <- geomDraw %>% st_as_sf()
      }
    }
    
    if (is.null(my_bbox)) {
      showNotification("Please draw or upload a polygon", type = "error")
      return()
    }
    
    # parse date range
    timestart <- as.POSIXct(paste0(input$fromdate," 00:00:00"), tz = "Australia/Sydney")
    timeend <- as.Date(input$fromdate) + (as.numeric(input$ndays) - 1)
    timeend <- as.POSIXct(paste0(timeend," 23:59:00"), tz = "Australia/Sydney")
    print(timestart); print(timeend)
    
    # append searched message
    text_list( append(text_list(), list(paste("Searched", timestart, "to", timeend))) )
    
    # withProgress: run searches for selected sources
    withProgress(message = "Processing...", value = 0, {
      
      # Track whether any files were written (we'll use file listing of dir.out later)
      inc <- 0.0
      
      # --- Progression polygons ---
      if (isTRUE(input$progscheck)) {
        incProgress(0.15, detail = "Progressions")
        progs <- tryCatch(
          fire_search_progs(my_bbox, timestart, timeend, return_geom = input$gather, dbpassword = "progressiondatabase",RDS=TRUE),
          error = function(e) { NULL }
        )
        
        if (is.data.frame(progs) && nrow(progs) > 0) {
          text_list( append(text_list(), list(paste0(length(unique(progs$dt_local)), " progression polygons found"))) )
          if (isTRUE(input$gather)) {
            # write progression polygons if gather requested
            outp <- file.path(dir.out, "progs.gpkg")
            sf::st_write(progs, outp, delete_dsn = TRUE, quiet = TRUE)
          }
        } else {
          text_list( append(text_list(), list("0 progression polygons found")) )
        }
        inc <- inc + 0.15
      }
      
      # --- Line scans ---
      if (isTRUE(input$imagescheck)) {
        incProgress(0.15, detail = "Line scans")
        scans <- tryCatch(
          fire_search_images(my_bbox, timestart, timeend, dbpassword = "progressiondatabase",RDS=TRUE) %>% select(-s3path),
          error = function(e) { NULL }
        )
        
        if (is.data.frame(scans) && nrow(scans) > 0) {
          text_list( append(text_list(), list(paste0(nrow(scans), " line scans found"))) )
          if (isTRUE(input$gather)) {
            sf::st_write(scans, file.path(dir.out, "scan_footprints.gpkg"), delete_dsn = TRUE, quiet = TRUE)
            
            scannames <- stringr::str_split(scans$s3name,";") %>% unlist() 
            scannames <- scannames[!is.na(scannames)]
            
            for(i in 1:length(scannames)){
              save_object(
                object = glue::glue("{scannames[i]}(lowres).tif"),
                bucket = "line-scans-low-res2",
                file = file.path(dir.out,glue::glue("{scannames[i]}(lowres).tif"))
              )
            }
            
            
            
          }
        } else {
          text_list( append(text_list(), list("0 scans found")) )
        }
        inc <- inc + 0.15
      }
      
      # --- Landsat/Sentinel (STAC) search ---
      if (isTRUE(input$GAstac)) {
        incProgress(0.15, detail = "Landsat and Sentinel 2")
        stac <- tryCatch(
          fire_stac_search(my_bbox, timestart, timeend),
          error = function(e) { NULL }
        )
        
        if (is.data.frame(stac) && nrow(stac) > 0) {
          # assuming stac has datetimelocal_chr column as in your original code
          nimgs <- if ("datetimelocal_chr" %in% names(stac)) length(unique(stac$datetimelocal_chr)) else nrow(stac)
          text_list( append(text_list(), list(paste0(nimgs, " Landsat/Sentinel images found."))) )
          if (isTRUE(input$gather)) {
            # attempt to download - wrap in tryCatch to avoid stopping the app
            tryCatch({
              fire_stac_download(my_bbox, stac, dest_folder = dir.out)
            }, error = function(e) {
              text_list( append(text_list(), list(paste0("STAC download error: ", e$message))) )
            })
          }
        } else {
          text_list( append(text_list(), list("0 Landsat or Sentinel images found")) )
        }
        inc <- inc + 0.15
      }
      
      # --- MODIS/VIIRS (GIBS) search ---
      if (isTRUE(input$gibs)) {
        incProgress(0.15, detail = "VIIRS and MODIS")
        # For GIBS we may not produce files unless gather==TRUE; we still report counts generically.
        # If fire_hotspot_map returns files (e.g., PNGs), it should write to dir.out
        text_list( append(text_list(), list("GIBS (MODIS/VIIRS) - multiple images per day possible")) )
        if (isTRUE(input$gather)) {
          tryCatch({
            fire_hotspot_map(my_bbox, timestart, timeend, dest_folder = dir.out, outwidth = 1200)
          }, error = function(e) {
            text_list( append(text_list(), list(paste0("GIBS error: ", e$message))) )
          })
        }
        inc <- inc + 0.15
      }
      
      # --- Himawari imagery search ---
      if (isTRUE(input$H89) || isTRUE(input$H89pyro)) {
        incProgress(0.15, detail = "Himawari")
        him <- tryCatch(
          fire_himawari_search(my_bbox, timestart, timeend, timestep_minutes = as.numeric(input$him_minutes)),
          error = function(e) { NULL }
        )
        
        if (is.data.frame(him) && nrow(him) > 0) {
          text_list( append(text_list(), list(paste0(length(unique(him$datetimelocal_chr)), " Himawari images")) ) )
          if (isTRUE(input$gather)) {
            # download default bands if H89 checked
            tryCatch({
              if (isTRUE(input$H89)) {
                fire_himawari_download(my_bbox, him, dest_folder = dir.out)
              }
              if (isTRUE(input$H89pyro)) {
                fire_himawari_download(my_bbox, him, dest_folder = dir.out, bands = "B14")
              }
            }, error = function(e) {
              text_list( append(text_list(), list(paste0("Himawari download error: ", e$message))) )
            })
          }
        } else {
          text_list( append(text_list(), list("0 Himawari images")) )
        }
        inc <- inc + 0.15
      }
      
      # Optionally bump progress to near completion
      setProgress(min(0.9, inc + 0.05), detail = "Finalising")
      Sys.sleep(0.1) # tiny pause to allow UI update
    }) # end withProgress
    
    # --- Prepare ZIP file for download only if gather requested ---
    if (isTRUE(input$gather)) {
      zip_file <- file.path(dirname(dir.out), "data_out.zip")
      if (file.exists(zip_file)) file.remove(zip_file)
      
      # gather files actually written by the searches
      filelist <- list.files(dir.out, full.names = TRUE,
                             pattern = "\\.(gpkg|tif|shp|png|jpg)$")
      
      # filter out temporary / system files if necessary (example)
      filelist <- filelist[!stringr::str_starts(basename(filelist), "file|spat_")]
      
      if (length(filelist) > 0) {
        # create zip (flatten paths)
        zip(zip_file, files = filelist, extras = "-j")
        zip_file_path(zip_file)
        text_list( append(text_list(), list("Data ready. Click Download ZIP")) )
      } else {
        text_list( append(text_list(), list("No files to prepare")) )
      }
    }
  }) # end observeEvent runquery
  
  # ---- DOWNLOAD HANDLER ----
  output$downloadScans <- downloadHandler(
    filename = function() { "fire_data.zip" },
    content = function(file) {
      zip_file <- zip_file_path()
      if (!is.null(zip_file) && file.exists(zip_file)) {
        file.copy(zip_file, file)
      }
    },
    contentType = "application/zip"
  )
  
  # ---- ZIP AVAILABILITY FLAG FOR CONDITIONAL PANEL ----
  output$zipAvailable <- reactive({
    !is.null(zip_file_path()) && file.exists(zip_file_path())
  })
  outputOptions(output, "zipAvailable", suspendWhenHidden = FALSE)
  
  # ---- SEARCH RESULTS TEXT ----
  output$display_text <- renderUI({
    msgs <- text_list()
    lapply(msgs, function(txt) tags$p(txt))
  })
}

# Run the Shiny application--------------------------------------------------
shinyApp(ui = ui, server = server)
