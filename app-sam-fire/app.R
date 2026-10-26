# Shiny App ----------------------------------------------------------
# Multiple TIFFs as toggleable layers
# with polygon drawing + send visible layers to second map
# and show ON/OFF checkboxes for sent layers.
# Drawn features are cleared after each SAM run.
# Now includes: Cut fire_SAM output with drawn polygons.


options(shiny.maxRequestSize = 500*1024^2)

# Sys.setenv(TMPDIR = "/srv/tmp")
# Sys.setenv(TEMP = "/srv/tmp")
# Sys.setenv(TEMPDIR = "/srv/tmp")
# terraOptions(tempdir="/srv/tmp")
# 
# reticulate::py_run_string("
# import os
# os.environ['TMPDIR'] = '/srv/tmp'
# os.environ['TEMP'] = '/srv/tmp'
# os.environ['TMP']  = '/srv/tmp'
# ")


library(shiny)
library(leaflet)
library(leaflet.extras)
library(sortable)
library(leafem)
library(terra)
library(sf)
library(dplyr)
library(smoothr)

#this one for docker build
reticulate::use_python("/srv/venv/bin/python", required = TRUE)
wk_dir <- "/srv/tmp"
chkp_dir='/srv/venv/checkpoints'

# reticulate::use_condaenv(condaenv = "my_sam2_env")
# wk_dir <- "D:\\temp\\wk_dir"
# chkp_dir="D:\\temp\\checkpoints"

# reticulate::use_python("C:/Users/mstorey/Documents/.virtualenvs/sam2_env_testa/Scripts/python.exe", required = TRUE)
# wk_dir   <- "C:/Users/mstorey/Documents/sam2/wk_dir"
# chkp_dir <- "C:/Users/mstorey/Documents/sam2/checkpoints"

sam2_models <- list(
 # "Tuned"   = list(pt = "sam2_finetuned.pt",       yaml = "sam2.1_hiera_t.yaml"),
  "Tiny"   = list(pt = "sam2.1_hiera_tiny.pt",       yaml = "sam2.1_hiera_t.yaml"),
  "Small"  = list(pt = "sam2.1_hiera_small.pt",      yaml = "sam2.1_hiera_s.yaml"),
  "Base+"  = list(pt = "sam2.1_hiera_base_plus.pt",  yaml = "sam2.1_hiera_b+.yaml"),
  "Large"  = list(pt = "sam2.1_hiera_large.pt",      yaml = "sam2.1_hiera_l.yaml")
)

fireobs:::fire_load_sam_once(chkp_dir,
  pt_name   = sam2_models[["Base+"]]$pt,
  yaml_name = sam2_models[["Base+"]]$yaml
)

 dir.create(wk_dir)

# Helper Functions----------------------------------------------------------
#  Map Creation and Raster Handling


# Create base leaflet map with default tiles and center/zoom
create_leaflet <- function() {
  leaflet() %>%
    addTiles() %>%
    setView(lng = 149.95, lat = -32.1, zoom = 9)
}

# Add a 3-band raster to an existing leaflet proxy map
# proxy: leafletProxy object
# r: terra raster object
# layer_id: string to identify layer group
# zoom_first: boolean; zoom to raster bbox if TRUE
add_raster_to_map <- function(proxy, r, layer_id, zoom_first = FALSE) {
  bbox_sf <- sf::st_as_sfc(sf::st_bbox(r)) |> sf::st_transform(4326)
  bbox <- sf::st_bbox(bbox_sf)
  proxy %>%
    addRasterRGB(r, r = 1, g = 2, b = 3, group = layer_id) %>%
    {
      if (zoom_first) fitBounds(., bbox[[1]], bbox[[2]], bbox[[3]], bbox[[4]]) else .
    }
}

# Create leaflet map with merged raster and optionally polygons overlayed
# Adds draw toolbar for polygon/rectangle/circle marker drawing
create_leaflet_with_layers <- function(r_merge, polygons = NULL) {
  map <- leaflet() %>%
    addTiles() %>%
    addRasterRGB(r_merge, r = 1, g = 2, b = 3, group = "image") %>%
    addDrawToolbar(
      targetGroup = "draw",
      polylineOptions = FALSE,
      circleOptions = FALSE,
      circleMarkerOptions = drawCircleMarkerOptions(
        repeatMode = TRUE,
        fillColor = "black",
        color = "yellow",
        opacity = 1,
        weight = 1,
        fillOpacity = 1
      ),
      polygonOptions = drawPolygonOptions(
        shapeOptions = drawShapeOptions(color = "red"),
        repeatMode = TRUE
      ),
      rectangleOptions = drawRectangleOptions(
        shapeOptions = drawShapeOptions(color = "red"),
        repeatMode = TRUE
      ),
      markerOptions = FALSE,
      editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())
    )
  
  # Add polygons if provided
  if (!is.null(polygons) && nrow(polygons) > 0) {
    map <- map %>%
      addPolygons(
        data = polygons,
        color = "yellow",
        weight = 5,
        fillColor = "lightgrey",
        fillOpacity = 0.6,
        group = "auto-polygons",
        popup = ~paste0("ID: ", id),
        layerId = ~as.character(id),
        highlightOptions = highlightOptions(color = "yellow", weight = 3, bringToFront = TRUE)
      )
  }
  
  # Determine bounding box (polygons if present, else raster)
  bbox <- if (!is.null(polygons) && nrow(polygons) > 0) {
    sf::st_bbox(polygons)
  } else {
    sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(r_merge)), 4326))
  }
  
  overlay_groups <- if (!is.null(polygons) && nrow(polygons) > 0) {
    c("image", "auto-polygons")
  } else {
    "image"
  }

  map %>%
    fitBounds(bbox[[1]], bbox[[2]], bbox[[3]], bbox[[4]]) %>%
    addLayersControl(
      overlayGroups = overlay_groups,
      options = layersControlOptions(collapsed = FALSE)
    )
}

# Clean and combine sf polygons by unioning and making valid multipolygons
# new_sf: new polygons to add
# old_sf: existing polygons to combine with (optional)
clean_combine_sf <- function(new_sf, old_sf = NULL) {
  if (is.null(old_sf)) {
    combined <- new_sf
  } else {
    st_geometry(old_sf) <- "geometry"
    st_geometry(new_sf) <- "geometry"
    old_sf <- old_sf %>% select(geometry)
    new_sf <- new_sf %>% select(geometry)
    combined <- rbind(
      st_sf(geometry = st_union(old_sf)),
      st_sf(geometry = st_union(new_sf))
    )%>% st_as_sf()
  }
  
  combined <- combined[!st_is_empty(combined), ]
  combined <- combined %>%
    st_union() %>%
    st_make_valid() %>%
    st_sfc(crs = 4326) %>%
    st_sf(geometry = .) %>%   # ✅ ensures geometry column is named
    st_cast("MULTIPOLYGON") %>%
    st_cast("POLYGON")
  combined <- combined[!st_is_empty(combined), ]
  combined$id <- seq_len(nrow(combined))
  combined
}

# Format YYYYmmddHHMMSS string to YYYY-dd-mm HH:MM for display
format_datetime_display <- function(dt_str) {
  if (is.null(dt_str) || nchar(dt_str) < 8) return(dt_str)
  year  <- substr(dt_str, 1, 4)
  month <- substr(dt_str, 5, 6)
  day   <- substr(dt_str, 7, 8)
  if (nchar(dt_str) >= 12) {
    hour <- substr(dt_str, 9, 10)
    min  <- substr(dt_str, 11, 12)
    paste0(year, "-", day, "-", month, " ", hour, ":", min)
  } else {
    paste0(year, "-", day, "-", month)
  }
}

# Help icon with CSS tooltip -----------------------------------------------

help_icon <- function(text) {
  tags$span(class = "help-icon-wrap",
    tags$span("?", class = "help-icon"),
    tags$span(text, class = "help-tooltip")
  )
}

# Button with styled hover tooltip appearing to the left
btn_tip <- function(id, label, tip_html) {
  div(class = "btn-tip-wrap",
      actionButton(id, label),
      div(class = "btn-tip", HTML(tip_html))
  )
}

help_css <- tags$style(HTML("
  .help-icon-wrap {
    position: relative;
    display: inline-block;
    margin-left: 4px;
    vertical-align: middle;
  }
  .help-icon {
    display: inline-block;
    width: 16px;
    height: 16px;
    line-height: 16px;
    text-align: center;
    border-radius: 50%;
    background: #6c757d;
    color: white;
    font-size: 11px;
    font-weight: bold;
    cursor: default;
    user-select: none;
  }
  .help-tooltip {
    visibility: hidden;
    opacity: 0;
    background-color: #333;
    color: #fff;
    text-align: left;
    padding: 8px 10px;
    border-radius: 6px;
    position: absolute;
    z-index: 9999;
    width: 280px;
    left: 22px;
    top: -4px;
    font-size: 12px;
    font-weight: normal;
    line-height: 1.4;
    transition: opacity 0.15s;
    white-space: normal;
  }
  .help-icon-wrap:hover .help-tooltip {
    visibility: visible;
    opacity: 1;
  }
  /* Button tooltip */
  .btn-tip-wrap {
    position: relative;
    width: 100%;
  }
  .btn-tip-wrap .btn {
    width: 100%;
  }
  .btn-tip {
    display: none;
    position: absolute;
    right: calc(100% + 10px);
    top: 0;
    background: #333;
    color: #fff;
    padding: 9px 12px;
    border-radius: 6px;
    width: 260px;
    font-size: 12px;
    line-height: 1.5;
    z-index: 9999;
    white-space: normal;
    box-shadow: 0 2px 8px rgba(0,0,0,0.35);
  }
  .btn-tip-wrap:hover .btn-tip {
    display: block;
  }
  .btn-tip b {
    display: block;
    margin-bottom: 4px;
    font-size: 13px;
  }
  .btn-tip ol, .btn-tip ul {
    margin: 4px 0 0 0;
    padding-left: 16px;
  }
  .btn-tip li {
    margin-bottom: 2px;
  }
  #download_gpkg {
    background-color: #28a745;
    border-color: #28a745;
    color: white;
  }
  #download_gpkg:hover {
    background-color: #218838;
    border-color: #1e7e34;
  }
  #sam2_model.shiny-input-container {
    margin-bottom: 6px !important;
  }
  .well .form-group {
    margin-bottom: 6px;
  }
  .well h4 {
    margin-top: 6px;
  }
"))

# UI Definition ----------------------------------------------------------
ui <- fluidPage(

  help_css,
  titlePanel("Semi-automated line scan mapping"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Load SAM2 model",
         help_icon("Choose which SAM2 model to use - tiny is faster but may be less accurate, large is slower but may be more accurate."),
         style = "margin-bottom: 4px;"),
      selectInput("sam2_model", label = NULL,
                  choices = names(sam2_models), selected = "Base+", width = "120px"),
      uiOutput("model_status_ui"),
      tags$hr(style = "margin: 8px 0;"),

      fileInput("tif_files", h4(
        "2. Upload 3-band TIFF(s)",
        help_icon("Click browse and select one or more geoTIFFS.
Should have start of file name be YYYYmmddHHMMSS so that time is correctly written in polygon.
We can add future ability to use jpgs/pngs.
Note that scans resolution is reduced to speed up testing.
We could use the original resolution in the future."),
        style = "margin-bottom: -10px;"
      ), accept = ".tif", multiple = TRUE),

      h4("3. Preview scan/s:",
         help_icon("Have an initial look at the scan/s then hit 'Send to mapping tool'.
For multiple scans, all tifs selected will be merged in the next step.
Click 'Send to mapping tool' to start using the auto mapping tool"),
         style = "margin-top: -18px; margin-bottom: 4px;"),

      div(style = "position: relative;",
          leafletOutput("map", height = 400),
          absolutePanel(
            top = 10, left = 44, width = "auto",
            style = "z-index: 500;",
            actionButton("send_to_second", "Send to mapping tool",
                         style = "font-size: 15px; padding: 8px 14px; white-space: nowrap;")
          )
      ),

      uiOutput("tif_order_ui")
    ),
    
    mainPanel(
      h4("4. Create polygons",
         help_icon("1. Draw multiple points or polygons over the fire area.
         One polygon and multiple points together is ok, but not mutliple polygons AND multiple polygons.
2. Click Auto-map button to create polygon.
3. Manually edit or run auto segment more times.")),
      
      # Map with floating buttons overlaid
      div(
        leafletOutput("second_map", height = 600),
        
        absolutePanel(
          top = 220, right = 20, width = "auto", draggable = FALSE,
          style = paste(
            "z-index: 500;",
            "background-color: rgba(255,255,255,0.9);",
            "padding: 10px;",
            "border-radius: 8px;",
            "display: flex;",
            "flex-direction: column;",
            "gap: 6px;"
          ),

          btn_tip("run_sam", "Auto-map", "
            <b>Auto-map</b>
            Draw points or polygons over the fire area first, then click to run the automated segmentation.
            <ol>
              <li>Draw points <em>or</em> polygons (not both)</li>
              <li>Click <b>Auto-map</b> to generate a polygon</li>
              <li>Draw more + click again to add more areas</li>
            </ol>"),
          btn_tip("remove_selected", "Remove Selected", "
            <b>Remove Selected</b>
            <ol>
              <li>Click polygon/s on the map to select them (turns blue)</li>
              <li>Click <b>Remove Selected</b> to delete them</li>
            </ol>"),
          btn_tip("keep_selected", "Keep Selected", "
            <b>Keep Selected</b>
            <ol>
              <li>Click polygon/s on the map to select them (turns blue)</li>
              <li>Click <b>Keep Selected</b> to remove all others</li>
            </ol>"),
          btn_tip("cut_polygons", "Cut", "
            <b>Cut</b>
            <ol>
              <li>Draw a polygon over the area to remove</li>
              <li>Click <b>Cut</b> to subtract it from existing polygons</li>
            </ol>
            <em>Draw over an entire polygon to delete it completely.</em>"),
          btn_tip("union_polygons", "Union", "
            <b>Union</b>
            <ol>
              <li>Draw a polygon over the area to add</li>
              <li>Click <b>Union</b> to merge it with existing polygons</li>
            </ol>
            <em>Use this to fill in areas the auto-mapper missed.</em>"),
          btn_tip("fill_holes", "Fill Holes", "
            <b>Fill Holes</b>
            Automatically fills small holes inside existing polygons.
            <br><em>For larger holes, use Draw + Union instead.</em>"),
          btn_tip("undo_last", "Undo", "
            <b>Undo</b>
            Reverts to the previous polygon state."),
          btn_tip("restart_mapping", "Restart", "
            <b>Restart</b>
            Removes all polygons and clears history so you can start again."),
        ),

      ),
      div(style = "display: flex; align-items: flex-start; gap: 0; margin-top: 10px; flex-wrap: wrap; border: 1px solid #dee2e6; border-radius: 6px; overflow: hidden;",
          div(style = "padding: 10px 16px;",
              uiOutput("smooth_ui")),
          div(style = "border-left: 1px solid #dee2e6; padding: 10px 16px; margin-left: auto;",
              uiOutput("download_ui"))
      )
    )
  )
)



# ui <- fluidPage(
#   titlePanel("Mapping from line scans"),
#   
#   sidebarLayout(
#     sidebarPanel(
#       # Upload one or multiple 3-band TIFF files
#       fileInput("tif_files", "Upload 3-band TIFF(s)", accept = ".tif", multiple = TRUE),
#       
#   
#       
#       # Show chosen directory for feedback
#       verbatimTextOutput("selected_dir"),
#       
#       #tags$hr(),
#       
#       # Button to merge selected TIF layers and send to second map
#       actionButton("send_to_second", "Merge Selected TIFS"),
#       
#       #tags$hr(),
#       
#       
#       h4("Preview and select tif/s:"),
#       
#       # Main leaflet map output
#       leafletOutput("map", height = 400),
#       # Dynamic UI to reorder TIFs
#       uiOutput("tif_order_ui")
#     ),
#     
#     mainPanel(
#       h4("Polygon mapping:"),
#       # Second map with merged raster + polygons + drawing tools
#       leafletOutput("second_map", height = 600),
#       textInput("out_dir", "Output gpkg folder", "C:/polygons"),
#       # Buttons for SAM segmentation and polygon manipulation
#       actionButton("run_sam", "Segment"),
#       actionButton("keep_selected", "Keep Selected Polygons"),
#       actionButton("cut_polygons", "Cut with drawn polygons"),
#       actionButton("union_polygons", "Union with drawn polygons"),
#       actionButton("undo_last", "Undo"),
# 
# 
#       
#  
#     )
#   )
# )

# Server Logic ----------------------------------------------------------

server <- function(input, output, session) {
  
  
  # Reactive Values Storage----
  vals <- reactiveValues(
    layers = list(),           # Stores names of layers uploaded
    rasters = list(),          # Stores terra raster objects for each layer
    selected_layers = NULL,    # Names of layers selected for merging
    r_merge = NULL,            # Merged raster for second map
    all_p = NULL,
    all_p_presmooth = NULL,    # Snapshot before smoothing; always smooth from this
    datetime = NULL,           # Datetime extracted from TIF name (original format)
    datetime_out = NULL,       # Datetime used for output (may be overridden by user)
    model_loaded = TRUE,       # Base+ loaded by default at startup
    model_name   = "Base+",

    # Combined polygons from SAM output
  )
  vals$selected_ids <- c()     # Stores IDs of selected polygons on second map

  # Load Model -----------------------------------

  output$model_status_ui <- renderUI({
    if (vals$model_loaded) {
      tags$p(paste0("\u2705 Model loaded: ", vals$model_name,
                    if (vals$model_name == "Base+") " (default)" else ""),
             style = "color: green; margin: 0; font-size: 13px;")
    } else {
      tags$p("\u26a0\ufe0f No model loaded. Click 'Load Model' to begin.",
             style = "color: #b8860b; margin: 0; font-size: 13px;")
    }
  })

  observeEvent(input$sam2_model, {
    model <- sam2_models[[input$sam2_model]]
    showNotification(paste0("Loading SAM2 ", input$sam2_model, " model..."), type = "message")
    tryCatch({
      fireobs:::fire_load_sam_once(checkpoints_dir = chkp_dir,
        pt_name   = model$pt,
        yaml_name = model$yaml,
        load_anyway = TRUE
        
      )
      vals$model_loaded <- TRUE
      vals$model_name   <- input$sam2_model
      showNotification(paste0("SAM2 ", input$sam2_model, " model loaded."), type = "message")
    }, error = function(e) {
      showNotification(paste0("Failed to load model: ", e$message), type = "error")
    })
  })

  vals$history <- list()
  push_history <- function() {
    if (!is.null(vals$all_p)) {
      vals$history[[length(vals$history) + 1]] <- vals$all_p
    }
    vals$all_p_presmooth <- NULL  # reset so next smooth uses fresh polygon
  }
  
  
  output$download_gpkg <- downloadHandler(
    filename = function() {
      req(vals$all_p, vals$ordered_layers)
      paste0(vals$datetime_out, "_polygons.gpkg")
    },
    content = function(file) {
      req(vals$all_p, vals$ordered_layers)
      tifs_used <- paste(vals$ordered_layers, collapse = ";")
      sf_to_save <- vals$all_p %>%
        dplyr::mutate(tifs_used = tifs_used, datetime = vals$datetime_out)
      sf::st_write(sf_to_save, file, delete_dsn = TRUE, quiet = TRUE)
    }
  )

  output$download_image <- downloadHandler(
    filename = function() {
      req(vals$r_merge)
      paste0(vals$datetime_out, "_image.tif")
    },
    content = function(file) {
      req(vals$r_merge)
      terra::writeRaster(vals$r_merge, file, overwrite = TRUE, datatype = "INT1U")
    }
  )
  
  # Render Main Map -----------------------------------
  # (empty base map on load)
  
  output$map <- renderLeaflet({ create_leaflet() })
  
  # Render Second Map-----------------------------------
  #  (empty base map on load)
  
  output$second_map <- renderLeaflet({
    leaflet() |> addTiles() |> setView(lng = 149.95, lat = -32.1, zoom = 9)
  })
  
  # Draw second map -----------------------------------
  # Reactive expression to read drawn features from second map
  # Converts GeoJSON features drawn by user into sf object
  
  all_drawn_sf <- reactive({
    req(input$second_map_draw_all_features)
    sf_obj <- tryCatch(
      geojsonsf::geojson_sf(jsonlite::toJSON(input$second_map_draw_all_features, auto_unbox = TRUE)),
      error = function(e) NULL
    )
    if (is.null(sf_obj) || nrow(sf_obj) == 0) {
      sf::st_sf(id = integer(0), geometry = sf::st_sfc(crs = 4326))
    } else {
      sf_obj
    }
  })
  
  
  
  # Load tiffs -----------------------------------
  # Load and display uploaded TIFFs on main map
  # Clears previous images and layers, adds each TIFF as a toggleable layer
  
  observeEvent(input$tif_files, {
    req(input$tif_files)
    vals$layers <- list()
    vals$rasters <- list()
    
    proxy <- leafletProxy("map") |> clearImages()
    
    input_files_sorted <- input$tif_files  # Initial order preserved
    
    for (i in seq_len(nrow(input_files_sorted))) {
      file <- input_files_sorted[i, ]
      tif_name <- tools::file_path_sans_ext(file$name)
      r <- terra::rast(file$datapath) |> terra::aggregate(5)
      if (terra::nlyr(r) < 3) {
        showNotification(paste("TIFF", file$name, "must have at least 3 bands!"), type = "error")
        next
      }
      layer_id <- paste0(tif_name)
      add_raster_to_map(proxy, r, layer_id, zoom_first = (i == 1))
      vals$layers[[layer_id]] <- tif_name
      vals$rasters[[layer_id]] <- r
    }
    
    proxy |> addLayersControl(
      overlayGroups = names(vals$layers),
      position = "bottomleft",
      options = layersControlOptions(collapsed = FALSE)
    )
  })
  
  # TIF order ui -----------------------------------
  # UI output for draggable TIF order list (sortable)
  # Lets user reorder TIFF layers for merging order
  
  output$tif_order_ui <- renderUI({
    req(vals$layers)
    tagList(
      tags$p(
        "Drag to set TIF order:",
        help_icon("If there are multiple scans, change the order by dragging the names - the scan at the top will be on top for the merge, and the output polygon will be assigned that scan's time.")
      ),
      rank_list(
        text = NULL,
        labels = names(vals$layers),
        input_id = "tif_order",
        options = sortable_options(multiDrag = TRUE)
      )
    )
  })
  
  # TIF order observe -----------------------------------
  # Observe changes in TIF order from draggable UI
  # Redraws the main map layers in new order and updates layers control
  
  observeEvent(input$tif_order, {
    req(vals$rasters, input$tif_order)
    
    proxy <- leafletProxy("map") %>%
      clearImages() %>%
      clearControls()
    in_tifs <- input$tif_order
    
    for (i in seq_along(in_tifs)) {
      tif_id <- in_tifs[i]
      r <- vals$rasters[[tif_id]]
      if (!is.null(r)) {
        add_raster_to_map(proxy, r, tif_id, zoom_first = (i == 1))
      }
    }
    
    # labels <- as.character(in_tifs)
    # labels[[1]] <- paste0("(top): ", labels[[1]])
    # labels[[length(labels)]] <- paste0("(bottom): ", labels[[length(labels)]])
    # 
    
    print(in_tifs)
    proxy %>% addLayersControl(
      overlayGroups = in_tifs,
      position = "bottomleft",
      options = layersControlOptions(collapsed = FALSE)
    )
  })
  
  # TIFS merge -----------------------------------
  # Merge selected TIFF layers and display on second map
  # Respects order from draggable UI, puts main TIFF first
  observeEvent(input$send_to_second, {
    
    vals$all_p <- NULL
    vals$history <- list()  # Also clear undo history if you want
    
    active_layers <- input$map_groups
    
    ordered_tifs <- input$tif_order
    main_tif <- ordered_tifs[1]
    
    vals$selected_layers <- active_layers
    
    if (!is.null(active_layers) && length(active_layers) > 0) {
      # Reorder so main_tif is first
      ordered_layers <- c(main_tif, setdiff(active_layers, main_tif))
      
      r_list <- lapply(ordered_layers, function(id) vals$rasters[[id]])
      r_list <- r_list[!sapply(r_list, is.null)]
      
      if (length(r_list) > 0) {
        r_merge <- terra::merge(terra::sprc(r_list)) |> terra::stretch()
        vals$r_merge <- r_merge
        vals$ordered_layers <- ordered_layers  # Store for saving metadata
        vals$datetime <- sub("_.*", "", ordered_layers[[1]])
        vals$datetime_out <- vals$datetime
        output$second_map <- renderLeaflet({ create_leaflet_with_layers(r_merge) })
      }
    }
  })
  
  # Restart -----------------------------------
  # Clear polygons and history but keep the merged raster on the second map

  observeEvent(input$restart_mapping, {
    req(vals$r_merge)
    vals$all_p <- NULL
    vals$history <- list()
    vals$selected_ids <- c()
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge)
    })
  })

  # SAM -----------------------------------
  # Run SAM segmentation on merged raster with drawn polygons/points as input
  # Updates polygons in vals$all_p and redraws second map with SAM output polygons
  # Clears drawn features on main map after completion
  
  observeEvent(input$run_sam, {
    req(vals$r_merge)
    sf_obj2 <- all_drawn_sf()
    if (nrow(sf_obj2) == 0) {
      showNotification("No drawn features detected!", type = "warning")
      return()
    }
    
    sf_obj2 <- sf_obj2[!sf::st_is_empty(sf_obj2), ] |> sf::st_make_valid()
    geom_types <- tolower(sf::st_geometry_type(sf_obj2))
    polygons_sf <- sf_obj2[geom_types %in% c("polygon", "multipolygon"), ]
    points_sf <- sf_obj2[geom_types %in% c("point"), ]
    
    if (nrow(polygons_sf) == 0) polygons_sf <- NULL
    if (nrow(points_sf) == 0) points_sf <- NULL
    
    vals$polygons_sf <- polygons_sf
    vals$points_sf <- points_sf
    
    # Write merged raster temporarily to disk for SAM processing
    r_merge_path <- file.path(wk_dir, "r_merge.tif")
    terra::writeRaster(vals$r_merge, r_merge_path, overwrite = TRUE, datatype = "INT1U")
    
    vals$selected_ids <- c()
    showNotification("Running segmentation", type = "message")
    
    # Run SAM (fire_SAM) - external function, returns polygons sf object
    p <- fireobs::fire_SAM(
      r_merge_path, polygons_sf, points_sf, working_dir = wk_dir
    )
    
    # Combine new SAM polygons with previously stored polygons
    combined_p <- clean_combine_sf(p, vals$all_p)
    
    push_history()
    
    vals$all_p <- combined_p
    
    # Update second map with raster + polygons
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, combined_p)
    })
    
    # Clear drawn features on second map
    leafletProxy("second_map") |> clearGroup("draw")
    
    showNotification("fire_SAM finished! Output added.", type = "message")
  })
  
  # Cut -----------------------------------
  # Cut existing fire_SAM polygons with drawn polygons on second map
  # Uses st_difference, updates polygons and redraws second map
  
  observeEvent(input$cut_polygons, {
    req(vals$all_p)
    drawn_sf <- all_drawn_sf()
    drawn_sf <- drawn_sf[!sf::st_is_empty(drawn_sf), ]
    if (nrow(drawn_sf) == 0) {
      showNotification("No drawn polygons to cut with!", type = "warning")
      return()
    }
    drawn_polygons <- drawn_sf[tolower(sf::st_geometry_type(drawn_sf)) %in% c("polygon", "multipolygon"), ]
    if (nrow(drawn_polygons) == 0) {
      showNotification("No valid polygon geometries found.", type = "warning")
      return()
    }
    diff_result <- sf::st_difference(vals$all_p, sf::st_union(drawn_polygons)) %>%
      sf::st_cast("MULTIPOLYGON") %>%
      sf::st_cast("POLYGON") %>%
      sf::st_make_valid()
    
    diff_result$id <- seq_len(nrow(diff_result))
    
    push_history()
    
    vals$all_p <- diff_result
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, diff_result)
    })
    showNotification("fire_SAM output cut by drawn polygons!", type = "message")
  })
  
  # Union -----------------------------------
  # Union existing fire_SAM polygons with drawn polygons on second map
  # Combines polygons into a single unified sf object and updates second map
  
  observeEvent(input$union_polygons, {
    req(vals$all_p)
    
    drawn_sf <- all_drawn_sf()
    drawn_sf <- drawn_sf[!sf::st_is_empty(drawn_sf), ]
    
    if (nrow(drawn_sf) == 0) {
      showNotification("No drawn polygons to union with!", type = "warning")
      return()
    }
    
    drawn_polygons <- drawn_sf[tolower(sf::st_geometry_type(drawn_sf)) %in% c("polygon", "multipolygon"), ]
    
    if (nrow(drawn_polygons) == 0) {
      showNotification("No valid polygon geometries found.", type = "warning")
      return()
    }
    
    # Make sure both sides are valid
    drawn_polygons <- sf::st_make_valid(drawn_polygons)
    all_p_valid <- sf::st_make_valid(vals$all_p)
    
    # Combine both as one sf object to safely union
    combined <- rbind(
      all_p_valid |> dplyr::select(geometry),
      drawn_polygons |> dplyr::select(geometry)
    )
    
    combined <- combined[!sf::st_is_empty(combined), ]
    
    if (nrow(combined) == 0) {
      showNotification("Combined result is empty!", type = "warning")
      return()
    }
    
    # Do the union of all geometries
    union_geom <- sf::st_union(combined)
    
    # Wrap union result back to sf object
    union_sfc <- sf::st_sfc(union_geom, crs = 4326)
    union_sfc <- union_sfc[!sf::st_is_empty(union_sfc)]
    
    union_result <- sf::st_sf(geometry = union_sfc)
    union_result$id <- seq_len(nrow(union_result))
    
    push_history()
    
    vals$all_p <- union_result
    
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, union_result)
    })
    
    showNotification("fire_SAM output unioned with drawn polygons!", type = "message")
  })
  
  # Helpers for polygon styles --------------------------------------------------

  default_poly_style <- list(
    color = "yellow", weight = 5,
    fillColor = "lightgrey", fillOpacity = 0.6
  )
  selected_poly_style <- list(
    color = "#00bfff", weight = 5,
    fillColor = "#00bfff", fillOpacity = 0.45
  )

  apply_poly_style <- function(proxy, poly, style) {
    proxy %>% addPolygons(
      data = poly,
      color = style$color, weight = style$weight,
      fillColor = style$fillColor, fillOpacity = style$fillOpacity,
      group = "auto-polygons",
      popup = ~paste0("ID: ", id),
      layerId = ~as.character(id),
      highlightOptions = highlightOptions(color = style$color, weight = 3, bringToFront = TRUE)
    )
  }

  # Click polys -----------------------------------
  # Toggle selection on click; highlight selected polygons in blue, revert deselected

  observeEvent(input$second_map_shape_click, {
    click <- input$second_map_shape_click
    req(click$id, vals$all_p)

    cid <- click$id
    poly <- vals$all_p[as.character(vals$all_p$id) == cid, ]
    if (nrow(poly) == 0) return()

    proxy <- leafletProxy("second_map")

    if (cid %in% vals$selected_ids) {
      vals$selected_ids <- setdiff(vals$selected_ids, cid)
      apply_poly_style(proxy, poly, default_poly_style)
    } else {
      vals$selected_ids <- unique(c(vals$selected_ids, cid))
      apply_poly_style(proxy, poly, selected_poly_style)
    }
  })

  # Keep -----------------------------------
  # Keep only selected polygons, discard the rest; re-render map cleanly

  observeEvent(input$keep_selected, {
    req(vals$all_p)
    if (length(vals$selected_ids) == 0) {
      showNotification("No polygons selected!", type = "warning")
      return()
    }
    kept <- vals$all_p %>% dplyr::filter(as.character(id) %in% vals$selected_ids)

    push_history()

    vals$all_p <- kept
    vals$selected_ids <- c()
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, kept)
    })
    showNotification(paste0("Kept ", nrow(kept), " polygon(s)."))
  })

  # Remove selected -----------------------------------
  # Remove selected polygons, keep the rest; re-render map cleanly

  observeEvent(input$remove_selected, {
    req(vals$all_p)
    if (length(vals$selected_ids) == 0) {
      showNotification("No polygons selected!", type = "warning")
      return()
    }
    kept <- vals$all_p %>% dplyr::filter(!as.character(id) %in% vals$selected_ids)

    push_history()

    vals$all_p <- kept
    vals$selected_ids <- c()
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, kept)
    })
    showNotification(paste0("Removed selected polygon(s). ", nrow(kept), " remaining."))
  })

  # Save gpkg-----------------------------------
  # Save current polygons to a GeoPackage file
  # Includes metadata fields from selected TIFs and datetime extracted from TIF name
  # Automatically creates output folder if missing
  # Generates unique filename by appending suffix if file exists
  
  # observeEvent(input$save_to_gpkg, {
  #   req(vals$all_p, vals$selected_layers, input$out_dir)
  #   
  #   if (!dir.exists(input$out_dir)) dir.create(input$out_dir)
  #   
  #   tif_names <- vals$ordered_layers %||% vals$selected_layers
  #   tifs_used <- paste(tif_names, collapse = ";")
  #   
  #   # Extract datetime string from first TIF name, assuming format like "20240701_xyz"
  #   datetime <- sub("_.*", "", tif_names[[1]])
  #   
  #   base_name <- paste0(datetime, "_SAM_output")
  #   save_dir <- input$out_dir
  #   ext <- ".gpkg"
  #   save_path <- file.path(save_dir, paste0(base_name, ext))
  #   
  #   # Append suffix if file exists
  #   i <- 1
  #   while (file.exists(save_path)) {
  #     suffix <- paste0("(", i, ")")
  #     save_path <- file.path(save_dir, paste0(base_name, suffix, ext))
  #     i <- i + 1
  #   }
  #   
  #   # Add metadata columns and save
  #   sf_to_save <- vals$all_p %>%
  #     dplyr::mutate(tifs_used = tifs_used, datetime = datetime)
  #   
  #   sf::st_write(sf_to_save, save_path, delete_dsn = TRUE, quiet = TRUE)
  #   
  #   showNotification(paste("GeoPackage saved to:", save_path), type = "message")
  # })
  
  
  # Undo button -------------------------------------------------------------
  
  observeEvent(input$undo_last, {
    if (length(vals$history) > 0) {
      vals$all_p <- vals$history[[length(vals$history)]]
      vals$history <- vals$history[-length(vals$history)]
      
      output$second_map <- renderLeaflet({
        create_leaflet_with_layers(vals$r_merge, vals$all_p)
      })
      
      showNotification("Undo successful!", type = "message")
    } else {
      showNotification("No previous state to undo.", type = "warning")
    }
  })
  
  
  # Download button (only shown after fire_SAM creates polygons) -----------

  output$smooth_ui <- renderUI({
    if (!is.null(vals$all_p) && nrow(vals$all_p) > 0) {
      div(
        sliderInput("smooth_amount",
                    tagList("Smooth edges:", help_icon("Test different levels of edge smoothing.")),
                    min = 0, max = 4, value = 0, step = .5, width = "220px"),
        actionButton("apply_smooth", "Apply")
      )
    }
  })

  observeEvent(input$apply_smooth, {
    req(vals$all_p)
    if (input$smooth_amount == 0) {
      showNotification("Move the slider above 0 to apply smoothing.", type = "warning")
      return()
    }
    # Use the latest unsmoothed polygon as the base.
    # all_p_presmooth is reset to NULL by push_history() whenever any real edit
    # (cut, union, auto-map, etc.) happens, so after an edit this always picks
    # up the freshly edited polygon. Repeated slider clicks reuse the same base
    # so vertices don't compound.
    if (is.null(vals$all_p_presmooth)) {
      vals$all_p_presmooth <- vals$all_p
    }
    # Push to history manually without calling push_history(), which would reset
    # all_p_presmooth and break repeated slider adjustments
    vals$history[[length(vals$history) + 1]] <- vals$all_p
    base_clean <- vals$all_p_presmooth |>
      sf::st_make_valid() |>
      sf::st_collection_extract("POLYGON") |>
      sf::st_cast("MULTIPOLYGON") |>
      sf::st_cast("POLYGON")
    smoothed <- smoothr::smooth(base_clean, method = "ksmooth",
                                smoothness = input$smooth_amount)
    smoothed <- sf::st_make_valid(smoothed)
    smoothed$id <- seq_len(nrow(smoothed))
    vals$all_p <- smoothed
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, smoothed)
    })
    showNotification("Smoothing applied.", type = "message")
  })

  output$download_ui <- renderUI({
    if (!is.null(vals$all_p) && nrow(vals$all_p) > 0) {
      div(
          div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 8px;",
              tags$label("Polygon time:", style = "margin: 0; font-weight: bold; white-space: nowrap;"),
              textInput("polygon_time", label = NULL,
                        value = format_datetime_display(vals$datetime),
                        width = "220px"),
              actionButton("apply_datetime", "Apply", style = "padding: 5px 10px;")
          ),
          div(style = "display: flex; align-items: center; gap: 10px;",
              downloadButton("download_gpkg", "Download GeoPackage"),
              downloadButton("download_image", "Download Merged Image")
          )
      )
    }
  })

  observeEvent(input$apply_datetime, {
    req(input$polygon_time)
    vals$datetime_out <- trimws(input$polygon_time)
    showNotification(paste0("Polygon time set to: ", vals$datetime_out), type = "message")
  })

  # Fill holes --------------------------------------------------------------
  
  observeEvent(input$fill_holes, {
    req(vals$all_p)
    
    # Fill holes in each polygon geometry using lwgeom::fill_holes
    # (lwgeom package provides st_fill_holes to fill polygon holes)
    
    if (!requireNamespace("nngeo", quietly = TRUE)) {
      showNotification("Please install 'nngeo' package to fill holes.", type = "error")
      return()
    }
    
    filled <- sf::st_sf(
      geometry = nngeo::st_remove_holes(vals$all_p$geometry),
      crs = sf::st_crs(vals$all_p)
    )
    
    # Preserve id column if exists
    if ("id" %in% colnames(vals$all_p)) {
      filled$id <- vals$all_p$id
    } else {
      filled$id <- seq_len(nrow(filled))
    }
    
    push_history()
    
    vals$all_p <- filled
    
    output$second_map <- renderLeaflet({
      create_leaflet_with_layers(vals$r_merge, vals$all_p)
    })
    
    showNotification("Holes filled in polygons!", type = "message")
  })
  
  
  
  
  
}

# Launch the app----------------------------------------------------------

shinyApp(ui, server)

