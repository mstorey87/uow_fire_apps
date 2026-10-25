library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)
library(dplyr)
library(lubridate)
library(fireobs)

# ---- helper: safe CRS handling ----
fix_crs <- function(x) {
  if (is.na(st_crs(x))) {
    st_crs(x) <- 4326
  }
  st_transform(x, 4326)
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Load Polygons and generate spread lines"),
  
  fluidRow(
    
    # ---- LEFT PANEL: controls ----
    column(3,
           wellPanel(
             fileInput("file_upload", "Upload GPKG or Zipped Shapefile",
                       accept = c(".gpkg", ".zip")),
             
             uiOutput("datetime_col_ui"),
             
             uiOutput("timezone_ui"),
             
             uiOutput("datetime_validation_ui"),
             
             tags$hr(),
             
             checkboxInput("maxlineonly", "Max line only",  TRUE),
             checkboxInput("cnvx",        "Convex hulls",  FALSE),
             
             tags$hr(),
             
             actionButton("measurespread", "Measure Spread Lines", class = "btn-primary btn-block"),
             tags$br(),
             tags$br(),
             actionButton("start_again",   "Start Over",           class = "btn-danger btn-block")
           ),
           
           verbatimTextOutput("filtered_data_info")
    ),
    
    # ---- RIGHT PANEL: map + results table ----
    column(9,
           leafletOutput("map", height = "70vh"),
           
           tags$br(),
           downloadButton("download_lines", "Download Spread Lines (GPKG)",
                          class = "btn-success",
                          style = "width: 50%;"),
           tags$br(),
           
           wellPanel(
             h4("Spread Line Results"),
             tableOutput("spread_table")
           )
    )
  )
)

# ---- SERVER ----
server <- function(input, output, session) {
  
  polygons_sf_r <- reactiveVal(NULL)
  drawn_polygon <- reactiveVal(NULL)
  spread_lines  <- reactiveVal(NULL)
  dt_col_r      <- reactiveVal(NULL)
  
  # ---- MAP INIT ----
  output$map <- renderLeaflet({
    leaflet() %>%
           addProviderTiles("OpenStreetMap",     group = "Street") %>% 
      addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
      addProviderTiles("OpenTopoMap",       group = "Topo") %>%
      addLayersControl(
        baseGroups = c("Street","Satellite",  "Topo"),
        position = "bottomleft",
        options = layersControlOptions(collapsed = FALSE)
      ) %>%
      addDrawToolbar(
        targetGroup = "draw",
        polylineOptions = FALSE,
        circleOptions = FALSE,
        circleMarkerOptions = FALSE,
        polygonOptions = drawPolygonOptions(
          shapeOptions = drawShapeOptions(color = "yellow", fillColor = "yellow", fillOpacity = 0.1),
          repeatMode = TRUE
        ),
        rectangleOptions = drawRectangleOptions(
          shapeOptions = drawShapeOptions(color = "yellow", fillColor = "yellow", fillOpacity = 0.1),
          repeatMode = TRUE
        ),
        markerOptions = FALSE,
        editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())
      ) %>%
      setView(lng = 134.0, lat = -25.0, zoom = 4)
  })
  
  # ---- helper: draw polygons with datetime hover label ----
  draw_polygons <- function(dat, dt_col) {
    
    bb <- st_bbox(dat)
    
    # Add display_label as a column so leaflet references it per feature
    if (!is.null(dt_col) && dt_col %in% names(dat)) {
      dat$display_label <- as.character(dat[[dt_col]])
      dat <- dat %>% arrange(desc(.data[[dt_col]]))
    } else {
      dat$display_label <- as.character(seq_len(nrow(dat)))
    }
    
    leafletProxy("map") %>%
      clearGroup("polygons") %>%
      clearGroup("spread_lines") %>%
      addPolygons(
        data         = dat,
        color        = "red",
        fillOpacity  = 0,
        weight       = 2.5,
        group        = "polygons",
        label        = ~display_label,
        labelOptions = labelOptions(
          style     = list("font-weight" = "normal", "font-size" = "12px"),
          direction = "bottom",
          offset    = c(0, 10)
        )
      ) %>%
      fitBounds(
        lng1 = as.numeric(bb["xmin"]),
        lat1 = as.numeric(bb["ymin"]),
        lng2 = as.numeric(bb["xmax"]),
        lat2 = as.numeric(bb["ymax"])
      )
  }
  
  # ---- TIMEZONE DROPDOWN ----
  observe({
    req(polygons_sf_r())
    output$timezone_ui <- renderUI({
      selectInput(
        "timezone",
        "Timezone of datetime column:",
        choices = c(
          "UTC"                = "UTC",
          "Australia/Canberra" = "Australia/Canberra",
          "Australia/Brisbane" = "Australia/Brisbane",
          "Australia/Adelaide" = "Australia/Adelaide",
          "Australia/Perth"    = "Australia/Perth",
          "Australia/Darwin"   = "Australia/Darwin"
        ),
        selected = if (!is.null(input$timezone)) input$timezone else "UTC"
      )
    })
  })
  
  # ---- DATETIME VALIDATION MESSAGE ----
  observe({
    req(polygons_sf_r(), input$datetime_col, input$timezone)
    
    dat    <- polygons_sf_r()
    dt_col <- input$datetime_col
    tz     <- input$timezone
    
    output$datetime_validation_ui <- renderUI({
      tryCatch({
        vals <- as.POSIXct(dat[[dt_col]], tz = tz)
        if (all(is.na(vals)) || class(dat[[dt_col]])[1] %in% c("integer", "factor", "numeric")) {
          tags$p("\u2716 Not a valid datetime column. Must be POSIXct or text e.g. 2010-12-20 12:50:00",
                 style = "color: red; font-size: 11px; margin-top: 5px;")
        } else {
          example <- as.character(vals[!is.na(vals)][1])
          tags$p(paste0("\u2714 Valid datetime. Example: ", example, " (", tz, ")"),
                 style = "color: green; font-size: 11px; margin-top: 5px;")
        }
      }, error = function(e) {
        tags$p("\u2716 Not a valid datetime column. Must be POSIXct or text e.g. 2010-12-20 12:50:00",
               style = "color: red; font-size: 11px; margin-top: 5px;")
      })
    })
  })
  
  # ---- AUTO LOAD USER DATA ----
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    file <- input$file_upload
    ext  <- tools::file_ext(file$name)
    
    tryCatch({
      
      if (ext == "gpkg") {
        dat <- st_read(file$datapath, quiet = TRUE)
      } else if (ext == "zip") {
        temp_dir <- file.path(tempdir(), paste0("shp_", as.integer(Sys.time())))
        dir.create(temp_dir)
        unzip(file$datapath, exdir = temp_dir)
        shp_files <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)
        req(length(shp_files) > 0)
        dat <- st_read(shp_files[1], quiet = TRUE)
      } else stop("Unsupported file type")
      
      # ---- FIX CRS + GEOMETRY ----
      dat <- dat %>%
        fix_crs() %>%
        st_make_valid() %>%
        filter(!st_is_empty(.))
      
      req(nrow(dat) > 0)
      
      # ---- STANDARDISE ID ----
      dat <- dat %>%
        mutate(id = if ("rowid" %in% names(.)) rowid else row_number())
      
      polygons_sf_r(dat)
      
      # ---- POPULATE DATETIME DROPDOWN ----
      all_cols       <- names(dat)
      datetime_guess <- all_cols[sapply(dat, function(x) inherits(x, c("POSIXct", "POSIXlt", "Date")))]
      default_col    <- if (length(datetime_guess) > 0) datetime_guess[1] else all_cols[1]
      
      dt_col_r(default_col)
      
      output$datetime_col_ui <- renderUI({
        selectInput(
          "datetime_col",
          "Select datetime column:",
          choices  = all_cols,
          selected = default_col
        )
      })
      
      draw_polygons(dat, default_col)
      
      output$filtered_data_info <- renderPrint({
        list(
          message    = paste("Loaded:", file$name),
          n_features = nrow(dat),
          crs        = st_crs(dat),
          columns    = all_cols
        )
      })
      
    }, error = function(e){
      showModal(modalDialog(title = "Error loading file", e$message))
    })
  })
  
  # ---- UPDATE POLYGON LABELS WHEN DATETIME COL CHANGES ----
  observeEvent(input$datetime_col, {
    req(polygons_sf_r(), input$datetime_col)
    dt_col_r(input$datetime_col)
    draw_polygons(polygons_sf_r(), input$datetime_col)
  })
  
  # ---- DRAWN POLYGON ----
  observeEvent(input$map_draw_new_feature, {
    drawn_polygon(input$map_draw_new_feature)
  })
  
  # ---- START AGAIN ----
  observeEvent(input$start_again, {
    
    req(polygons_sf_r())
    
    drawn_polygon(NULL)
    spread_lines(NULL)
    
    dat_reset <- polygons_sf_r()
    if (!is.null(dt_col_r()) && dt_col_r() %in% names(dat_reset)) {
      dat_reset$display_label <- as.character(dat_reset[[dt_col_r()]])
      dat_reset <- dat_reset %>% arrange(desc(.data[[dt_col_r()]]))
    } else {
      dat_reset$display_label <- as.character(seq_len(nrow(dat_reset)))
    }
    
    leafletProxy("map") %>%
      clearGroup("spread_lines") %>%
      clearGroup("convex_hulls") %>%
      removeDrawToolbar(clearFeatures = TRUE) %>%
      clearShapes() %>%
      addPolygons(
        data         = dat_reset,
        color        = "red",
        fillOpacity  = 0,
        weight       = 2.5,
        group        = "polygons",
        label        = ~display_label,
        labelOptions = labelOptions(
          style     = list("font-weight" = "normal", "font-size" = "12px"),
          direction = "bottom",
          offset    = c(0, 10)
        )
      ) %>%
      addDrawToolbar(
        targetGroup = "draw",
        polylineOptions = FALSE,
        circleOptions = FALSE,
        circleMarkerOptions = FALSE,
        polygonOptions = drawPolygonOptions(
          shapeOptions = drawShapeOptions(color = "yellow", fillColor = "yellow", fillOpacity = 0.1),
          repeatMode = TRUE
        ),
        rectangleOptions = drawRectangleOptions(
          shapeOptions = drawShapeOptions(color = "yellow", fillColor = "yellow", fillOpacity = 0.1),
          repeatMode = TRUE
        ),
        markerOptions = FALSE,
        editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())
      )
    
    output$spread_table       <- renderTable({ NULL })
    output$filtered_data_info <- renderPrint({ "Reloaded. Draw a polygon and measure spread." })
  })
  
  # ---- SPREAD LINES ----
  observeEvent(input$measurespread, {
    
    if (is.null(drawn_polygon())) {
      showModal(modalDialog(
        title = "No Drawn Polygon Detected",
        "Draw a box over the polygons (or section) where you want to measure spread.",
        easyClose = TRUE
      ))
      return()
    }
    
    if (is.null(input$datetime_col) || input$datetime_col == "") {
      showModal(modalDialog(
        title = "No Datetime Column Selected",
        "Please select the datetime column from the dropdown before measuring spread.",
        easyClose = TRUE
      ))
      return()
    }
    
    req(polygons_sf_r())
    
    dat       <- polygons_sf_r()
    dt_col    <- input$datetime_col
    tz_to_use <- if (!is.null(input$timezone) && input$timezone != "") input$timezone else "UTC"
    
    # ---- VALIDATE DATETIME ----
    valid <- tryCatch({
      test_vals <- as.POSIXct(dat[[dt_col]], tz = tz_to_use)
      !all(is.na(test_vals)) && !class(dat[[dt_col]])[1] %in% c("integer", "factor", "numeric")
    }, error = function(e) FALSE)
    
    if (!valid) {
      showModal(modalDialog(
        title = "Invalid Datetime Column",
        "Must be a POSIXct column or text in format e.g. 2010-12-20 12:50:00",
        easyClose = TRUE
      ))
      return()
    }
    
    # ---- BUILD dt FROM USER-SELECTED COLUMN + TIMEZONE ----
    tryCatch({
      dat <- dat %>%
        mutate(dt = as.POSIXct(.data[[dt_col]], tz = tz_to_use))
    }, error = function(e) {
      showModal(modalDialog(
        title = "Datetime Conversion Failed",
        paste("Could not convert column", dt_col, "using timezone", tz_to_use, ":", e$message),
        easyClose = TRUE
      ))
      return()
    })
    
    coords <- matrix(unlist(drawn_polygon()$geometry$coordinates[[1]]),
                     ncol = 2, byrow = TRUE)
    
    poly <- st_sfc(st_polygon(list(coords)), crs = 4326) %>%
      st_transform(st_crs(dat))
    
    filtered <- dat %>%
      st_intersection(poly) %>%
      st_transform(3112)
    
    if (nrow(filtered) < 2) {
      showModal(modalDialog(
        title = "Insufficient Data",
        "Need multiple polygons/time steps inside the drawn polygon."
      ))
      return()
    }
    
    lines <- fire_max_spread_line(
      filtered,
      "dt",
      "id",
      convex_hull = input$cnvx,
      max_only    = input$maxlineonly
    ) %>%
      st_transform(4326) %>%
    filter( !prior_intersects)  %>%
      mutate(
        layerId       = row_number(),
        line_id       = row_number(),
        timestep_mins = round(timestep_mins, 0)
      )
    
    #sf::st_write(lines,"D:\\temp\\lines_x.gpkg")
    
    spread_lines(lines)
    
    leafletProxy("map") %>%
      clearGroup("spread_lines") %>%
      clearGroup("convex_hulls") %>%
      addPolylines(
        data         = lines,
        color        = "black",
        weight       = 4,
        opacity      = 1,
        layerId      = ~layerId,
        group        = "spread_lines",
        label        = ~paste0(
          "Line: ",          line_id,
          " | Mins: ",       timestep_mins,
          " | ROS (km/h): ", round(ros_kmh, 2)
        ),
        labelOptions = labelOptions(
          style     = list("font-weight" = "bold", "font-size" = "13px"),
          textsize  = "13px",
          direction = "bottom",
          offset    = c(0, 10)
        )
      )
    
    # ---- DISPLAY CONVEX HULLS IF CHECKED ----
    if (input$cnvx) {
      hulls <- filtered %>%
        st_transform(4326) %>%
        mutate(
          display_label = as.character(dt)
        ) %>% 
        st_convex_hull()
      
      leafletProxy("map") %>%
        addPolygons(
          data         = hulls,
          color        = "orange",
          fillOpacity  = 0.1,
          weight       = 1.5,
          group        = "convex_hulls",
          label        = ~display_label,
          labelOptions = labelOptions(
            style     = list("font-weight" = "normal", "font-size" = "12px"),
            direction = "bottom",
            offset    = c(0, 10)
          )
        )
    }
    
    # ---- DISPLAY TABLE ----
    output$spread_table <- renderTable({
      lines %>%
        st_drop_geometry() %>%
        select(any_of(c("line_id", "start_time", "end_time", "timestep_mins", "line_km", "ros_kmh"))) %>%
        mutate(
          across(where(~ inherits(.x, "POSIXct")), as.character),
          timestep_mins = round(timestep_mins)
        )
    })
  })
  
  # ---- DOWNLOAD SPREAD LINES ----
  output$download_lines <- downloadHandler(
    filename = function() {
      paste0("spread_lines_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".gpkg")
    },
    content = function(file) {
      req(spread_lines())
      spread_lines() %>%
        select(any_of(c("line_id", "start_time", "end_time", "timestep_mins", "line_km", "ros_kmh", "geometry"))) %>%
        mutate(across(where(~ inherits(.x, "POSIXct")), as.character)) %>%
        st_write(file, driver = "GPKG", quiet = TRUE)
    }
  )
}

shinyApp(ui, server)