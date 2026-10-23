
# --- app.R ---
# Combines:
#   - Bayesian quantile ROS prediction (brms wind × FMI model)
#   - Existing fire observation map, weather, and imagery tabs

library(shiny)
library(ggplot2)
library(dplyr)
library(here)
library(brms)
library(purrr)
library(glue)
library(scales)
library(htmltools)
library(leaflet)
library(sf)
library(mapview)
library(rnaturalearth)
library(wfprogression)
library(lubridate)
library(DBI)
library(RPostgres)

# ============================================================
# --- LOAD BAYESIAN MODELS ---
# ============================================================

STAN_DIR <- "D:/UOW_projects/NHRA_ROS/DataWorking/bayesian_model/stan_models"

q_models <- list(
  fit_q50 = readRDS(file.path(STAN_DIR, "brms_quantile_q50_ros_interaction_model3.rds")),
  fit_q90 = readRDS(file.path(STAN_DIR, "brms_quantile_q90_ros_interaction_model3.rds"))
)
q_levels <- c(0.5, 0.9)

train_dat   <- q_models$fit_q50$data
pyro_levels <- levels(train_dat$pyroregion)
has_zmla    <- "zmla_mean" %in% names(train_dat)

q_cols   <- c("fit_q50" = "#2166ac", "fit_q90" = "#b2182b")
q_labels <- c("fit_q50" = "Median ROS  (q50)", "fit_q90" = "Extreme ROS  (q90)")

rng <- function(v) range(train_dat[[v]], na.rm = TRUE)
med <- function(v) round(median(train_dat[[v]], na.rm = TRUE), 2)

slider_min <- function(v, d = 1) floor(rng(v)[1]   * 10^d) / 10^d
slider_max <- function(v, d = 1) ceiling(rng(v)[2] * 10^d) / 10^d

calc_fmi <- function(temp, rh) 10 - 0.25 * (temp - rh)

fmi_data_min <- slider_min("fmi_surface")
fmi_data_max <- slider_max("fmi_surface")

# ============================================================
# --- FIRE SPREAD ELLIPSE GEOMETRY ---
# ============================================================

IGN_LON <- 150.421
IGN_LAT <- -33.517
KM_PER_DEG_LAT <- 111.0
KM_PER_DEG_LON <- 111.0 * cos(IGN_LAT * pi / 180)

fire_ellipse_pts <- function(ros_kmh, time_hr, lb_ratio, spread_dir_deg, n = 300) {
  if (is.na(ros_kmh) || ros_kmh <= 0) return(tibble(x = 0, y = 0))
  head_dist  <- ros_kmh * time_hr
  back_dist  <- head_dist / lb_ratio^2
  a          <- (head_dist + back_dist) / 2
  b          <- a / lb_ratio
  center_fwd <- (head_dist - back_dist) / 2
  dir_rad    <- spread_dir_deg * pi / 180
  fwd_x <-  sin(dir_rad); fwd_y <-  cos(dir_rad)
  rgt_x <-  cos(dir_rad); rgt_y <- -sin(dir_rad)
  cx <- center_fwd * fwd_x; cy <- center_fwd * fwd_y
  theta <- seq(0, 2 * pi, length.out = n + 1)
  tibble(
    x = a * cos(theta) * fwd_x + b * sin(theta) * rgt_x + cx,
    y = a * cos(theta) * fwd_y + b * sin(theta) * rgt_y + cy
  )
}

km_to_polygon <- function(df_xy) {
  coords <- df_xy %>%
    transmute(lon = IGN_LON + x / KM_PER_DEG_LON,
              lat = IGN_LAT + y / KM_PER_DEG_LAT) %>%
    as.matrix()
  if (!isTRUE(all.equal(coords[1, ], coords[nrow(coords), ])))
    coords <- rbind(coords, coords[1, ])
  st_polygon(list(coords))
}

# ============================================================
# --- DATABASE CONNECTION ---
# ============================================================

ensure_db_connection <- function(dbpassword) {
  if (is.null(DB) || !DBI::dbIsValid(DB)) {
    message("Connecting to database...")
    DB <<- DBI::dbConnect(
      RPostgres::Postgres(),
      dbname   = "cermb_fires",
      user     = "mstorey",
      password = dbpassword,
      host     = "charus.ad.uow.edu.au",
      port     = 5432
    )
  }
  DB
}

if (!exists("DB")) DB <- NULL
DB <- ensure_db_connection("bushfires")

# ============================================================
# --- UI ---
# ============================================================

ui <- fluidPage(
  titlePanel("Fire ROS Prediction"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Region"),
      selectInput("pyroregion", "Pyroregion",
                  choices  = pyro_levels,
                  selected = pyro_levels[1]),

      hr(),
      h4("Weather"),
      sliderInput("wind", "Wind speed — surface (km/h)",
                  min = slider_min("wind"), max = slider_max("wind"),
                  value = med("wind"), step = 0.5),

      sliderInput("temp", "Temperature (°C)",      min = 0,   max = 50,  value = 30, step = 0.5),
      sliderInput("rh",   "Relative humidity (%)", min = 1,   max = 100, value = 20, step = 1),
      helpText("FMI is calculated from temp & RH. Low FMI = hot/dry."),
      sliderInput("fmi",  "FMI (auto-calculated; drag to override)",
                  min = fmi_data_min, max = fmi_data_max,
                  value = round(calc_fmi(30, 20), 1), step = 0.1),

      hr(),
      h4("Soil & Stability"),
      sliderInput("sm_top", "Soil moisture 0–10 cm (kg/m²)",
                  min = slider_min("sm_top"), max = slider_max("sm_top"),
                  value = med("sm_top"), step = 0.5),
      sliderInput("lapse",  "Temp diff 850–700 hPa (°C)",
                  min = slider_min("lapse"), max = slider_max("lapse"),
                  value = med("lapse"), step = 0.5),

      hr(),
      h4("Terrain"),
      sliderInput("wbtsdslope", "SD slope (WBT)",
                  min = slider_min("wbtsdslope_mean"), max = slider_max("wbtsdslope_mean"),
                  value = med("wbtsdslope_mean"), step = 0.1),

      hr(),
      h4("Fire Spread (Map)"),
      sliderInput("time_hr",    "Time horizon (hours)",        min = 0.5, max = 6,   value = 1,  step = 0.5),
      sliderInput("spread_dir", "Spread direction (° from N)", min = 0,   max = 355, value = 90, step = 5),
      sliderInput("lb_ratio",   "Length-to-breadth ratio",     min = 1,   max = 5,   value = 2,  step = 0.1),

      hr(),
      helpText("200 posterior draws; debounced 700 ms.")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "mainTabs",

        tabPanel("ROS Plot",
                 plotOutput("rosPlot", height = "420px"),
                 uiOutput("rosSummary")),

        tabPanel("ROS Map",
                 leafletOutput("ellipseMap", height = "580px"),
                 uiOutput("ellipseSummary")),

        tabPanel("Australia Map",
                 h4("Observed fires with weather ±20% of slider values"),
                 actionButton("runBtn", "Map observations", icon = icon("play"), class = "btn-success"),
                 leafletOutput("ausMap", height = "600px")),

        tabPanel("Weather",
                 plotOutput("weatherPlot", height = "600px")),

        tabPanel("Imagery",
                 plotOutput("imageryPlot"))
      )
    )
  )
)

# ============================================================
# --- SERVER ---
# ============================================================

server <- function(input, output, session) {

  # Sync FMI slider to temp + rh
  observeEvent(c(input$temp, input$rh), {
    fmi_calc <- round(calc_fmi(input$temp, input$rh), 1)
    new_min  <- min(fmi_data_min, fmi_calc)
    new_max  <- max(fmi_data_max, fmi_calc)
    updateSliderInput(session, "fmi", value = fmi_calc, min = new_min, max = new_max)
  }, ignoreInit = FALSE)

  # Build newdata for model prediction
  newdata_react <- reactive({
    data.frame(
      pyroregion      = factor(input$pyroregion, levels = pyro_levels),
      wind            = input$wind,
      fmi_surface     = input$fmi,
      sm_top          = input$sm_top,
      lapse           = input$lapse,
      wbtsdslope_mean = input$wbtsdslope,
      fire_day_id     = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  newdata_d <- debounce(newdata_react, 700)

  pred_summary <- reactive({
    nd <- newdata_d()
    req(nrow(nd) == 1)

    map2_dfr(q_models, q_levels, function(fit, q) {
      d <- as.numeric(posterior_linpred(fit, newdata = nd, ndraws = 200,
                                        re_formula = NA, allow_new_levels = TRUE))
      pct_neg <- round(mean(d < 0) * 100, 1)
      tibble(
        q_key   = paste0("fit_q", q * 100),
        q_level = q,
        label   = q_labels[[paste0("fit_q", q * 100)]],
        median  = pmax(median(d), 0),
        lo      = pmax(quantile(d, 0.10), 0),
        hi      = pmax(quantile(d, 0.90), 0),
        pct_neg = pct_neg
      )
    }) %>%
      mutate(label = factor(label, levels = unname(q_labels)),
             q_key = factor(q_key,  levels = names(q_cols)))
  })

  # ── ROS Plot ──────────────────────────────────────────────
  output$rosPlot <- renderPlot({
    ps <- pred_summary()

    ggplot(ps, aes(x = label, y = median, colour = q_key)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 1.1) +
      geom_point(size = 5) +
      geom_label(aes(y = hi, label = sprintf("%.2f km/h", median)),
                 vjust = -0.4, fontface = "bold",
                 fill = "white", label.size = 0.4,
                 label.padding = unit(0.3, "lines")) +
      scale_colour_manual(values = q_cols, guide = "none") +
      scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0.05))) +
      labs(
        x = NULL, y = "Predicted ROS (km/h)",
        title    = "Predicted Rate of Spread — wind × FMI interaction model",
        subtitle = glue(
          "pyroregion={input$pyroregion}  |  ",
          "wind={input$wind} km/h  |  temp={input$temp}°C  |  RH={input$rh}%  →  ",
          "FMI={round(input$fmi, 1)}  |  ",
          "sm_top={input$sm_top}  |  lapse={input$lapse}°C  |  SD slope={input$wbtsdslope}"
        ),
        caption = "Point = posterior median  |  bars = 80% credible interval  |  200 posterior draws"
      ) +
      theme_bw(base_size = 13) +
      theme(plot.title    = element_text(face = "bold"),
            plot.subtitle = element_text(colour = "grey40", size = 10),
            plot.caption  = element_text(colour = "grey50"),
            panel.grid.major.x = element_blank())
  })

  output$rosSummary <- renderUI({
    ps   <- pred_summary()
    rows <- pmap_chr(ps, function(q_key, label, median, lo, hi, pct_neg, ...) {
      col  <- q_cols[[as.character(q_key)]]
      note <- if (pct_neg > 0)
        glue(' <span style="color:#d9a800;">({pct_neg}% draws floored at 0)</span>') else ""
      glue(
        '<span style="color:{col}; font-weight:bold;">{label}</span> &nbsp;|&nbsp; ',
        'ROS: <strong>{sprintf("%.2f", median)} km/h</strong> ',
        '<span style="color:#777;">(80% CI: {sprintf("%.2f", lo)}–{sprintf("%.2f", hi)} km/h)</span>',
        '{note}'
      )
    })
    HTML(paste0(
      '<div style="margin-top:8px; padding:10px 16px; background:#f0f4f8; ',
      'border-left:5px solid #2166ac; border-radius:4px; font-size:1.0em; line-height:2.6;">',
      paste(rows, collapse = "<br/>"), '</div>'
    ))
  })

  # ── ROS Map (fire spread ellipses) ────────────────────────
  ellipse_sf_react <- reactive({
    ps <- pred_summary()

    pmap_dfr(ps, function(q_key, median, lo, hi, ...) {
      qk       <- as.character(q_key)
      poly_hi  <- km_to_polygon(fire_ellipse_pts(hi,     input$time_hr, input$lb_ratio, input$spread_dir))
      poly_med <- km_to_polygon(fire_ellipse_pts(median, input$time_hr, input$lb_ratio, input$spread_dir))
      poly_lo  <- km_to_polygon(fire_ellipse_pts(lo,     input$time_hr, input$lb_ratio, input$spread_dir))

      ring <- tryCatch(
        st_difference(st_sfc(poly_hi, crs = 4326), st_sfc(poly_lo, crs = 4326)),
        error = function(e) st_sfc(poly_hi, crs = 4326)
      )

      bind_rows(
        st_sf(q_key = qk, band = "ring",   geometry = ring),
        st_sf(q_key = qk, band = "median", geometry = st_sfc(poly_med, crs = 4326))
      )
    }) %>%
      mutate(q_key = factor(q_key, levels = names(q_cols)))
  })

  output$ellipseMap <- renderLeaflet({
    leaflet() %>%
      addTiles(attribution = "© OpenStreetMap contributors") %>%
      setView(lng = IGN_LON, lat = IGN_LAT, zoom = 13) %>%
      addScaleBar(position = "bottomleft",
                  options = scaleBarOptions(metric = TRUE, imperial = FALSE)) %>%
      addControl(
        html = '<div style="background:white;padding:5px 8px;border-radius:4px;
                  border:1px solid #999;text-align:center;font-size:13px;line-height:1.3;">
                  <div style="font-size:20px;">&#8679;</div>
                  <div style="font-weight:bold;">N</div></div>',
        position = "topright"
      )
  })

  observe({
    ed    <- ellipse_sf_react()
    proxy <- leafletProxy("ellipseMap") %>% clearShapes() %>% clearMarkers()

    for (qk in levels(ed$q_key)) {
      col   <- q_cols[[qk]]
      proxy <- proxy %>%
        addPolygons(data = ed %>% filter(q_key == qk, band == "ring"),
                    fillColor = col, fillOpacity = 0.25,
                    color = col, weight = 1.0,
                    label = glue("{qk} — 80% CI band")) %>%
        addPolygons(data = ed %>% filter(q_key == qk, band == "median"),
                    fillOpacity = 0, color = col, weight = 2.5,
                    label = glue("{qk} — median"))
    }

    # Leader lines + labels at forward tip of each median ellipse
    ps      <- pred_summary()
    dir_rad <- input$spread_dir * pi / 180
    fwd_lon <- sin(dir_rad) / KM_PER_DEG_LON
    fwd_lat <- cos(dir_rad) / KM_PER_DEG_LAT
    rgt_lon <- cos(dir_rad) / KM_PER_DEG_LON
    rgt_lat <- -sin(dir_rad) / KM_PER_DEG_LAT

    label_offsets <- list(
      fit_q50 = list(text = "50% chance\nbeyond here", perp =  0.4, back =  0.5),
      fit_q90 = list(text = "10% chance\nbeyond here", perp =  0.4, back = -0.3)
    )

    for (qk in names(q_cols)) {
      ps_row <- ps %>% filter(q_key == qk)
      if (nrow(ps_row) == 0 || ps_row$median <= 0) next
      col     <- q_cols[[qk]]
      meta    <- label_offsets[[qk]]
      head_km <- ps_row$median * input$time_hr

      tip_lon <- IGN_LON + head_km * fwd_lon
      tip_lat <- IGN_LAT + head_km * fwd_lat
      lbl_lon <- tip_lon + meta$perp * rgt_lon - meta$back * fwd_lon
      lbl_lat <- tip_lat + meta$perp * rgt_lat - meta$back * fwd_lat

      proxy <- proxy %>%
        addCircleMarkers(lng = tip_lon, lat = tip_lat,
                         radius = 4, color = col, fillColor = col,
                         fillOpacity = 1, stroke = FALSE) %>%
        addPolylines(lng = c(tip_lon, lbl_lon), lat = c(tip_lat, lbl_lat),
                     color = col, weight = 1.5, opacity = 0.8, dashArray = "4 3") %>%
        addCircleMarkers(lng = lbl_lon, lat = lbl_lat,
                         radius = 0, opacity = 0, fillOpacity = 0,
                         label = HTML(paste0(
                           "<span style='font-weight:bold;font-size:12px;color:", col, ";",
                           "background:rgba(255,255,255,0.9);border:1.5px solid ", col, ";",
                           "border-radius:4px;padding:2px 6px;'>",
                           gsub("\n", "<br/>", meta$text), "</span>"
                         )),
                         labelOptions = labelOptions(permanent = TRUE, direction = "right",
                                                     textOnly = TRUE))
    }

    proxy <- proxy %>%
      addMarkers(lng = IGN_LON, lat = IGN_LAT, label = "Ignition",
                 popup = glue("Ignition — spread direction: {input$spread_dir}°"))

    hi_ed <- ed %>% filter(band == "ring")
    if (nrow(hi_ed) > 0 && any(!st_is_empty(hi_ed))) {
      bbox <- tryCatch(st_bbox(hi_ed), error = function(e) NULL)
      if (!is.null(bbox)) {
        pad <- max((bbox[["xmax"]] - bbox[["xmin"]]) * 0.15, 0.003)
        proxy %>% fitBounds(
          lng1 = bbox[["xmin"]] - pad, lat1 = bbox[["ymin"]] - pad,
          lng2 = bbox[["xmax"]] + pad, lat2 = bbox[["ymax"]] + pad
        )
      }
    }
  })

  output$ellipseSummary <- renderUI({
    ps   <- pred_summary()
    rows <- pmap_chr(ps, function(q_key, label, median, lo, hi, ...) {
      col     <- q_cols[[as.character(q_key)]]
      head    <- round(median * input$time_hr, 2)
      head_lo <- round(lo     * input$time_hr, 2)
      head_hi <- round(hi     * input$time_hr, 2)
      glue(
        '<span style="color:{col}; font-weight:bold;">{label}</span> &nbsp;|&nbsp; ',
        'ROS: {sprintf("%.2f", median)} km/h ',
        '<span style="color:#777;">(80% CI: {sprintf("%.2f", lo)}–{sprintf("%.2f", hi)})</span>',
        ' &nbsp;|&nbsp; Head in {input$time_hr} h: ',
        '<strong>{head} km</strong> <span style="color:#777;">({head_lo}–{head_hi} km)</span>'
      )
    })
    HTML(paste0(
      '<div style="margin-top:10px; padding:10px 16px; background:#f0f4f8; ',
      'border-left:5px solid #2166ac; border-radius:4px; font-size:1.0em; line-height:2.6;">',
      paste(rows, collapse = "<br/>"), '</div>'
    ))
  })

  # ── Australia observed fire map ───────────────────────────
  dat_ros_reactive <- eventReactive(input$runBtn, {
    fn_error <- function(x, error) c(x - x * error, x + x * error)
    isolate({
      wfprogression::fire_search_ros(
        temperature  = fn_error(input$temp, .2),
        windspeed    = fn_error(input$wind, .2),
        rh           = fn_error(input$rh,   .2),
        progressions = TRUE,
        dbpassword   = "bushfires"
      )
    })
  })

  selected_lineid <- reactiveVal(NULL)

  observeEvent(input$ausMap_shape_click, {
    click <- input$ausMap_shape_click
    if (!is.null(click$id)) {
      selected_lineid(as.numeric(click$id))
      updateTabsetPanel(session, inputId = "mainTabs", selected = "Weather")
    }
  })

  output$ausMap <- renderLeaflet({
    req(dat_ros_reactive())
    dat_ros <- dat_ros_reactive()
    d.ros   <- dat_ros$lines %>% mutate(ros = round(ros_kmh_epsg4283, 2))

    leaflet() %>%
      addTiles() %>%
      addPolygons(data = dat_ros$progressions, weight = 1, fillOpacity = 0) %>%
      addPolylines(
        data      = d.ros,
        color     = ~colorNumeric("viridis", ros)(ros),
        weight    = 4,
        opacity   = 0.8,
        layerId   = ~lineid,
        label     = ~paste("LineID:", lineid, ", ROS:", ros)
      )
  })

  # ── Weather plot for clicked line ────────────────────────
  output$weatherPlot <- renderPlot({
    req(selected_lineid(), dat_ros_reactive())
    dat_ros <- dat_ros_reactive()
    weath   <- dat_ros$weather %>% filter(lineid == selected_lineid())
    lines   <- dat_ros$lines   %>% filter(lineid == selected_lineid())

    combined <- weath %>%
      select(time_round, sfcwind_mean, tas_mean, hurs_mean, lineid) %>%
      st_drop_geometry() %>%
      left_join(lines %>% st_drop_geometry(), by = "lineid") %>%
      filter(time_round > start_time_utc - hours(2)) %>%
      mutate(
        time_round     = with_tz(time_round,     "Australia/Sydney"),
        start_time_utc = with_tz(start_time_utc, "Australia/Sydney"),
        end_time_utc   = with_tz(end_time_utc,   "Australia/Sydney")
      )

    ggplot(combined, aes(x = time_round)) +
      geom_line(aes(y = sfcwind_mean * 3.6, col = "wind (km/h)"),    linewidth = 1) +
      geom_line(aes(y = tas_mean - 273.15,  col = "temperature (°C)"), linewidth = 1) +
      geom_line(aes(y = hurs_mean,          col = "RH"),               linewidth = 1, linetype = "dashed") +
      geom_vline(aes(xintercept = start_time_utc), colour = "black", linewidth = 1, linetype = "dotted") +
      geom_vline(aes(xintercept = end_time_utc),   colour = "black", linewidth = 1, linetype = "dotted") +
      scale_color_manual(
        name   = "Variables",
        values = c("wind (km/h)" = "black", "temperature (°C)" = "red", "RH" = "blue")
      ) +
      scale_x_datetime(name = "Sydney time", date_labels = "%H:%M") +
      coord_cartesian(ylim = c(0, 100)) +
      theme_bw(base_size = 14) +
      labs(y = "Weather variables", x = "Time",
           title = paste("Weather for LineID:", selected_lineid()))
  })

  # ── Imagery ───────────────────────────────────────────────
  output$imageryPlot <- renderPlot({
    req(selected_lineid())
    cowplot::plot_grid(wfprogression::fire_read_plot(selected_lineid(), "bushfires"))
  })
}

# ============================================================
# --- RUN ---
# ============================================================

shinyApp(ui, server) |>
  shiny::runApp(host = "127.0.0.1", port = 5826)
