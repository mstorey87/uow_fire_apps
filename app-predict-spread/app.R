
# --- app.R ---
# Bayesian quantile ROS prediction (brms wind × FMI model)

library(shiny)
library(ggplot2)
library(dplyr)
library(purrr)
library(glue)
library(scales)
library(htmltools)
library(leaflet)
library(sf)
library(bsplus)

# Hover-triggered (?) popover helper, in the style of app-sam-4/app.R
help_icon <- function(content, placement = "right") {
  icon("question-circle") %>%
    bs_embed_popover(content = content, title = NULL, placement = placement,
                      trigger = "hover", html = TRUE)
}

# ============================================================
# --- LOAD PRECOMPUTED PREDICTION LOOKUP ---
# ============================================================
# This app predicts from a precomputed grid of the brms models' additive terms
# (built by precompute_lookup.R) rather than from the brmsfit objects directly,
# so it needs no brms/rstan/Stan at runtime. The model is additive:
#   mu ~ 0 + pyroregion + t2(wind, fmi_surface) + s(sm_top) + s(lapse) + s(wbtsdslope_mean)
# so each term is interpolated off its own grid, per posterior draw, and summed -
# this reproduces posterior_linpred() to within the grid's interpolation error
# (verified < 0.005 km/h against the fitted models; see precompute_lookup.R).

STAN_DIR <- Sys.getenv("STAN_DIR", "stan_models")

lookup   <- readRDS(file.path(STAN_DIR, "prediction_lookup.rds"))
q_levels <- c(0.5, 0.9)

pyro_levels <- lookup$fit_q50$pyro_levels
var_stats   <- lookup$fit_q50$var_stats

q_cols   <- c("fit_q50" = "#2166ac", "fit_q90" = "#b2182b")
q_labels <- c("fit_q50" = "Median ROS  (q50)", "fit_q90" = "Extreme ROS  (q90)")

rng <- function(v) c(var_stats[[v]]$min, var_stats[[v]]$max)
med <- function(v) round(var_stats[[v]]$median, 2)

slider_min <- function(v, d = 1) floor(rng(v)[1]   * 10^d) / 10^d
slider_max <- function(v, d = 1) ceiling(rng(v)[2] * 10^d) / 10^d

calc_fmi <- function(temp, rh) 10 - 0.25 * (temp - rh)

fmi_data_min <- slider_min("fmi_surface")
fmi_data_max <- slider_max("fmi_surface")

# --- Interpolation helpers (vectorized across posterior draws, clamped at grid edges) ---

interp1_draws <- function(grid, mat, x) {
  x  <- min(max(x, grid[1]), grid[length(grid)])
  ih <- min(max(findInterval(x, grid, all.inside = TRUE) + 1, 2), length(grid))
  il <- ih - 1
  w  <- (x - grid[il]) / (grid[ih] - grid[il])
  mat[, il] * (1 - w) + mat[, ih] * w
}

interp2_draws <- function(grid_x, grid_y, arr, x, y) {
  x <- min(max(x, grid_x[1]), grid_x[length(grid_x)])
  y <- min(max(y, grid_y[1]), grid_y[length(grid_y)])
  ixh <- min(max(findInterval(x, grid_x, all.inside = TRUE) + 1, 2), length(grid_x)); ixl <- ixh - 1
  iyh <- min(max(findInterval(y, grid_y, all.inside = TRUE) + 1, 2), length(grid_y)); iyl <- iyh - 1
  wx <- (x - grid_x[ixl]) / (grid_x[ixh] - grid_x[ixl])
  wy <- (y - grid_y[iyl]) / (grid_y[iyh] - grid_y[iyl])
  (1 - wx) * (1 - wy) * arr[, ixl, iyl] + wx * (1 - wy) * arr[, ixh, iyl] +
    (1 - wx) * wy       * arr[, ixl, iyh] + wx * wy       * arr[, ixh, iyh]
}

# Reconstructs posterior draws of the linear predictor for one quantile model
# by summing its interpolated additive terms - the lookup-table equivalent of
# posterior_linpred(fit, newdata, re_formula = NA, allow_new_levels = TRUE).
predict_draws <- function(lk, pyroregion, wind, fmi_surface, sm_top, lapse, wbtsdslope_mean) {
  lk$pyro_mat[, as.character(pyroregion)] +
    interp2_draws(lk$grid_wind, lk$grid_fmi, lk$t2_arr, wind, fmi_surface) +
    interp1_draws(lk$grid_sm_top, lk$sm_top_mat, sm_top) +
    interp1_draws(lk$grid_lapse, lk$lapse_mat, lapse) +
    interp1_draws(lk$grid_wbt, lk$wbt_mat, wbtsdslope_mean)
}

# ============================================================
# --- MCARTHUR MK5 & VESTA MK2 FIRE SPREAD MODELS ---
# ============================================================
# Independent empirical fire-behaviour models, included for comparison against
# the Bayesian model above. Equations from CSIRO SPARK's model library:
#   https://research.csiro.au/spark/resources/model-library/mcarthur/
#   https://research.csiro.au/spark/resources/model-library/dry-eucalypt-forest-vesta-mk-2/

# McArthur Mark 5 forest fire danger meter: FFDI, then head fire rate of spread.
# wind: km/h at 10 m in the open. drought_factor: 1-10. fuel_load: t/ha.
calc_mcarthur <- function(temp, rh, drought_factor, wind, fuel_load) {
  ffdi <- 2 * exp(-0.450 + 0.987 * log(drought_factor) - 0.0345 * rh + 0.0338 * temp + 0.0234 * wind)
  ros_kmh <- 0.0012 * ffdi * fuel_load
  list(ffdi = ffdi, ros = ros_kmh)
}

# Vesta Mk2 dry eucalypt forest model (Cheney/Gould/McCaw/Cruz).
# wind10: km/h at 10 m in the open. waf: wind adjustment factor (open-10m -> 2m under canopy).
# ws: surface + near-surface fuel load (t/ha). hu: understorey height (m).
# drought_factor: 0-10. slope_deg: slope in the spread direction, degrees (+ve = upslope).
# hour: 0-23, selects which fuel-moisture-content time-of-day equation applies.
calc_vesta <- function(wind10, temp, rh, waf, ws, hu, drought_factor, slope_deg, hour) {
  wind2 <- wind10 / waf

  mc <- if (hour >= 12 && hour < 17) {
    2.76 + 0.124 * rh - 0.0187 * temp
  } else if ((hour >= 9 && hour < 12) || (hour >= 17 && hour < 20)) {
    3.6 + 0.169 * rh - 0.045 * temp
  } else {
    3.08 + 0.198 * rh - 0.0483 * temp
  }

  moisture_coeff <- 0.9082 + 0.1206 * mc - 0.03106 * mc^2 + 0.001853 * mc^3 - 0.00003467 * mc^4
  moisture_coeff <- if (mc < 4.1) 1 else if (mc > 24) 0 else moisture_coeff

  fuel_availability <- 1.008 / (1 + 104.9 * exp(-0.9306 * drought_factor))
  moisture_coeff <- moisture_coeff * fuel_availability

  slope_coeff  <- 2^(0.1 * abs(slope_deg))
  slope_factor <- if (slope_deg >= 0) slope_coeff else slope_coeff / (2 * slope_coeff - 1)

  ros1_head <- if (wind2 < 2) {
    0.03 * moisture_coeff * slope_factor
  } else {
    (0.03 + 0.05024 * (wind2 - 1)^0.92628 * (ws * 0.1)^0.79928) * moisture_coeff * slope_factor
  }
  ros2_head <- 0.19591 * wind2^0.8257 * (ws * 0.1)^0.4672 * hu^0.495 * moisture_coeff * slope_factor
  # Phase III (Eq. 16 in Cruz 2021, "The Vesta Mk 2 rate of fire spread model: a
  # user's guide") has no slope term - only Phases I and II (Eq. 14a/15) do.
  ros3_head <- 0.05235 * wind10^1.19128 * moisture_coeff

  g2 <- -23.9315 + 1.7033 * wind2 + 12.0822 * moisture_coeff + 0.9524 * ws
  p2 <- if (ws < 1) 0 else 1 / (1 + exp(-g2))

  # P(III) = 0 if R(II) < 0.3 km/h is Eq. 11 in the user's guide - a deliberate,
  # documented threshold (not a bug): "a Phase II rate of fire spread constraint...
  # if the Phase II rate of fire spread is below 0.3 km/h the likelihood of fire
  # spread to Phase III is nil." At typical dry-eucalypt fuel loads this doesn't
  # create a visible jump because R(II) only crosses 0.3 once the logistic below is
  # already near zero; the earlier "fix" here (dropping this guard entirely) was
  # incorrect relative to the source and has been reverted.
  g3 <- -32.3074 + 0.2951 * wind10 + 26.8734 * moisture_coeff
  p3 <- if (ros2_head < 0.3) 0 else 1 / (1 + exp(-g3))

  ros_kmh <- if (p2 < 0.5) {
    ros1_head * (1 - p2) + ros2_head * p2
  } else {
    (ros1_head * (1 - p2) + ros2_head * p2) * (1 - p3) + ros3_head * p3
  }

  list(mc = mc, ros = ros_kmh)
}

# ============================================================
# --- FIRE SPREAD ELLIPSE GEOMETRY ---
# ============================================================

IGN_LON_DEFAULT <- 150.421
IGN_LAT_DEFAULT <- -33.517
KM_PER_DEG_LAT <- 111.0
km_per_deg_lon <- function(lat) 111.0 * cos(lat * pi / 180)

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

km_to_polygon <- function(df_xy, ign_lon, ign_lat) {
  kpd_lon <- km_per_deg_lon(ign_lat)
  coords <- df_xy %>%
    transmute(lon = ign_lon + x / kpd_lon,
              lat = ign_lat + y / KM_PER_DEG_LAT) %>%
    as.matrix()
  if (!isTRUE(all.equal(coords[1, ], coords[nrow(coords), ])))
    coords <- rbind(coords, coords[1, ])
  st_polygon(list(coords))
}

# ============================================================
# --- VALUE BOX (card) HELPER ---
# ============================================================
# Self-contained equivalent of the value/info-box style shown at
# https://www.listendata.com/2021/03/value-info-box-in-shiny-rmarkdown.html
# (colored left border, big value, faded icon accent) without adding the
# GitHub-only summaryBox package as a dependency.

# Perceived luminance (0-1) of a #RRGGBB/#RGB hex color, or NA if not a hex color
# (e.g. a CSS name like "white"). Used to decide whether a theme color needs a
# dark outline to stay legible against the card's white background - bright
# colors like white or yellow would otherwise wash out.
hex_luminance <- function(color) {
  if (!grepl("^#[0-9A-Fa-f]{3}([0-9A-Fa-f]{3})?$", color)) {
    return(if (toupper(color) == "WHITE") 1 else NA_real_)
  }
  hex <- sub("^#", "", color)
  if (nchar(hex) == 3) hex <- paste0(rep(strsplit(hex, "")[[1]], each = 2), collapse = "")
  rgb <- vapply(c(1, 3, 5), function(i) strtoi(substr(hex, i, i + 1), 16L) / 255, numeric(1))
  sum(rgb * c(0.2126, 0.7152, 0.0722))
}

value_box <- function(title, value, subtitle = NULL, note = NULL, icon_name = "chart-line", letter = NULL, color = "#2166ac", width = 6) {
  # Bright theme colors (white, yellow, ...) would otherwise be low-contrast or
  # invisible against the card's white background - add a dark outline around the
  # card and a dark text stroke so the border, value and badge stay legible while
  # still reading as the intended bright color.
  is_light      <- isTRUE(hex_luminance(color) > 0.75)
  card_outline  <- if (is_light) "border:1.5px solid #333; " else ""
  text_outline  <- if (is_light) "-webkit-text-stroke: 1px #333; text-shadow: 0 0 2px rgba(0,0,0,0.6); " else ""
  badge_opacity <- if (is_light) "0.55" else "0.2"

  badge <- if (!is.null(letter)) {
    div(letter, style = paste0(
      "position:absolute; right:14px; top:6px; font-size:48px; font-weight:800; ",
      "font-family:sans-serif; line-height:1; color:", color, "; opacity:", badge_opacity, "; ",
      text_outline
    ))
  } else {
    icon(icon_name, style = paste0(
      "position:absolute; right:16px; top:14px; font-size:42px; color:", color, "; opacity:0.18;"
    ))
  }

  column(
    width = width,
    div(
      style = paste0(
        "position:relative; overflow:hidden; background:#ffffff; ",
        card_outline,
        "border-left:6px solid ", color, "; border-radius:8px; ",
        "box-shadow:0 1px 4px rgba(0,0,0,0.12); padding:18px 20px; margin-bottom:20px;"
      ),
      badge,
      div(style = "font-size:0.95em; color:#555; font-weight:600; text-transform:uppercase; letter-spacing:0.03em;",
          title),
      div(style = paste0("font-size:2.1em; font-weight:700; color:", color, "; margin-top:4px; ", text_outline),
          value),
      if (!is.null(subtitle))
        div(style = "font-size:1.05em; color:#777; margin-top:4px;", subtitle),
      if (!is.null(note))
        div(style = "font-size:1em; color:#999; margin-top:10px; padding-top:10px; border-top:1px solid #eee;",
            note)
    )
  )
}

# Colors for the combined comparison map (Bayesian + McArthur + Vesta)
q_cols_compare <- c(q_cols, mcarthur = "#8E44AD", vesta = "#000000")

# Line pattern per model on the combined comparison map (SVG dash arrays)
q_dash_compare <- c(mcarthur = "1 12", vesta = "1 12")

# ============================================================
# --- UI ---
# ============================================================

ui <- fluidPage(
  bsplus::use_bs_popover(),
  titlePanel("Fire ROS Prediction"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Region"),
      selectInput("pyroregion",
                  label = tagList("Pyroregion", help_icon(paste0(
                    "Describes broad regions of similar fire characteristics, from the paper ",
                    "“Pyrogeography in flux: Reorganization of Australian fire regimes in a hotter ",
                    "world”. The default region (a), along with region d, covers most of the ",
                    "forested east of NSW and gives similar predictions."
                  ))),
                  choices  = pyro_levels,
                  selected = pyro_levels[1]),

      hr(),
      h4("Weather"),
      sliderInput("wind",
                  label = tagList("Wind speed — surface (km/h)", help_icon(paste0(
                    "Based on modelled/gridded weather data, not on-ground observations. This often ",
                    "underestimates the wind speed actually reaching the fire ground, by around 5–15 km/h."
                  ))),
                  min = slider_min("wind"), max = slider_max("wind"),
                  value = med("wind"), step = 0.5),

      sliderInput("temp",
                  label = tagList("Temperature (°C)", help_icon(
                    "Used, together with relative humidity, to calculate the Fuel Moisture Index (FMI) below."
                  )),
                  min = 0, max = 50, value = 30, step = 0.5),
      sliderInput("rh",
                  label = tagList("Relative humidity (%)", help_icon(
                    "Used, together with temperature, to calculate the Fuel Moisture Index (FMI) below."
                  )),
                  min = 1, max = 100, value = 20, step = 1),
      helpText("FMI is calculated from temp & RH. Low FMI = hot/dry."),
      sliderInput("fmi",
                  label = tagList("FMI (auto-calculated; drag to override)", help_icon(
                    "Sharples Fuel Moisture Index. Calculated as:<br><em>FMI = 10 − 0.25 × (Temperature − Relative Humidity)</em><br>Lower FMI means hotter, drier fuel."
                  )),
                  min = fmi_data_min, max = fmi_data_max,
                  value = round(calc_fmi(30, 20), 1), step = 0.1),

      hr(),
      h4("Soil, Stability & Terrain", help_icon(
        "Soil moisture, atmospheric stability and terrain ruggedness generally have a more minor effect on the predicted rate of spread than wind and FMI."
      )),
      sliderInput("sm_top", "Soil moisture 0–10 cm (kg/m²)",
                  min = slider_min("sm_top"), max = slider_max("sm_top"),
                  value = med("sm_top"), step = 0.5),
      sliderInput("lapse",
                  label = tagList("Temp diff 850–700 hPa (°C)", help_icon(
                    "Higher values indicate greater atmospheric instability. Highly correlated with (but simpler than) the C-Haines index. Most effect on the model occurs above around 12°C."
                  )),
                  min = slider_min("lapse"), max = slider_max("lapse"),
                  value = med("lapse"), step = 0.5),
      sliderInput("wbtsdslope",
                  label = tagList("SD slope (WBT)", help_icon(
                    "Standard deviation of slope across the general area (about 1 km windows). Higher values = more rugged terrain; lower values = flatter plains."
                  )),
                  min = slider_min("wbtsdslope_mean"), max = slider_max("wbtsdslope_mean"),
                  value = med("wbtsdslope_mean"), step = 0.1),

      hr(),
      helpText("200 posterior draws; debounced 700 ms.",
               br(), "Spread-period/direction/shape sliders for the ROS Map(s) are below the map on those tabs."),

      hr(),
      tags$details(
        tags$summary(
          style = "cursor:pointer; font-weight:bold; color:#555;",
          "McArthur & Vesta Inputs", help_icon(
            "Extra inputs used only by the McArthur Mk5 and Vesta Mk2 models (Model Comparison tabs). Wind, temperature and RH above are shared with these models."
          )
        ),
        div(
          style = "margin-top:10px;",
          sliderInput("drought_factor",
                      label = tagList("Drought Factor (1–10)", help_icon(
                        "Soil/fuel dryness index used by both models - higher values mean drier conditions and more fuel available to burn."
                      )),
                      min = 1, max = 10, value = 10, step = 0.5),
          sliderInput("fuel_load",
                      label = tagList("Fine fuel load (t/ha)", help_icon(
                        "Surface fuel load. Used as McArthur's fuel load and Vesta's surface fuel load (Ws)."
                      )),
                      min = 5, max = 35, value = 20, step = 0.5),
          sliderInput("waf",
                      label = tagList("Wind adjustment factor (WAF)", help_icon(
                        "Vesta only. Converts the open 10 m wind speed above down to the near-surface (2 m) wind under the forest canopy; typically 2-3 for forests."
                      )),
                      min = 1, max = 5, value = 2.5, step = 0.1),
          sliderInput("understorey_ht",
                      label = tagList("Understorey height (m)", help_icon(
                        "Vesta only. Height of the near-surface/elevated fuel layer."
                      )),
                      min = 0.1, max = 3, value = 0.5, step = 0.1),
          sliderInput("slope_deg",
                      label = tagList("Slope (° in spread direction)", help_icon(
                        "Vesta only. Positive = upslope, negative = downslope. Fire spreads faster upslope."
                      )),
                      min = -20, max = 20, value = 0, step = 1),
          sliderInput("fire_hour",
                      label = tagList("Time of day (hour)", help_icon(
                        "Vesta only. Used to pick the right fuel moisture equation for the time of day."
                      )),
                      min = 0, max = 23, value = 14, step = 1)
        )
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "mainTabs",

        tabPanel("ROS Plot",
                 plotOutput("rosPlot", height = "420px"),
                 uiOutput("rosSummary")),

        tabPanel("ROS Map",
                 helpText("Click anywhere on the map to set a new ignition point and redraw the ROS polygons there.",
                          br(), "Note: after switching to this tab, move a slider slightly to trigger the initial draw."),
                 leafletOutput("ellipseMap", height = "580px"),
                 hr(),
                 h4("Fire Spread Settings"),
                 fluidRow(
                   column(4, sliderInput("time_hr",    "Spread period (hours)",        min = 0.5, max = 6,   value = 1,  step = 0.5, width = "100%")),
                   column(4, sliderInput("spread_dir", "Spread direction (° from N)", min = 0,   max = 355, value = 90, step = 5, width = "100%")),
                   column(4, sliderInput("lb_ratio",   "Length-to-breadth ratio",     min = 1,   max = 5,   value = 2,  step = 0.1, width = "100%"))
                 ),
                 uiOutput("ellipseSummary")),

        tabPanel("ROS Description",
                 uiOutput("rosDescription")),

        tabPanel("Model Comparison",
                 br(),
                 uiOutput("modelComparisonCards"),
                 div(style = "margin-top:6px; color:#666; font-size:0.85em;",
                     "McArthur Mk5 and Vesta Mk2 are independent empirical fire-behaviour models shown here for ",
                     "comparison with the Bayesian model above. Their model-specific inputs are in the ",
                     "“McArthur & Vesta Inputs” drop-down at the bottom of the sidebar.")),

        tabPanel("Model Comparison Map",
                 helpText("Click anywhere on the map to set a new ignition point and redraw all four models' polygons there.",
                          br(), "Note: after switching to this tab, move a slider slightly to trigger the initial draw.",
                          br(), "Uses the same Spread Period / Direction / Length-to-breadth ratio sliders as the ROS Map tab."),
                 leafletOutput("ellipseMapCompare", height = "580px"),
                 uiOutput("ellipseSummaryCompare"))
      )
    )
  )
)

# ============================================================
# --- SERVER ---
# ============================================================

server <- function(input, output, session) {

  # Ignition point - defaults to the original fixed location, but the user can
  # click anywhere on either ROS map to move it (shared by both maps below).
  ignition <- reactiveVal(list(lon = IGN_LON_DEFAULT, lat = IGN_LAT_DEFAULT))

  observeEvent(input$ellipseMap_click, {
    click <- input$ellipseMap_click
    ignition(list(lon = click$lng, lat = click$lat))
  })

  observeEvent(input$ellipseMapCompare_click, {
    click <- input$ellipseMapCompare_click
    ignition(list(lon = click$lng, lat = click$lat))
  })

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

    map2_dfr(lookup, q_levels, function(lk, q) {
      d <- predict_draws(lk, nd$pyroregion, nd$wind, nd$fmi_surface,
                          nd$sm_top, nd$lapse, nd$wbtsdslope_mean)
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

  # ── McArthur Mk5 & Vesta Mk2 (shared by cards + combined map) ─
  model_compare_inputs <- reactive({
    list(
      wind = input$wind, temp = input$temp, rh = input$rh,
      drought_factor = req(input$drought_factor), fuel_load = req(input$fuel_load),
      waf = req(input$waf), understorey_ht = req(input$understorey_ht),
      slope_deg = req(input$slope_deg), fire_hour = req(input$fire_hour)
    )
  })

  model_compare_d <- debounce(model_compare_inputs, 700)

  model_compare_calc <- reactive({
    mi <- model_compare_d()
    list(
      mcarthur = calc_mcarthur(mi$temp, mi$rh, mi$drought_factor, mi$wind, mi$fuel_load),
      vesta    = calc_vesta(mi$wind, mi$temp, mi$rh, mi$waf, mi$fuel_load, mi$understorey_ht,
                             mi$drought_factor, mi$slope_deg, mi$fire_hour)
    )
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

  # ── ROS Description (plain-English text) ──────────────────
  output$rosDescription <- renderUI({
    ps  <- pred_summary()
    q50 <- ps %>% filter(q_key == "fit_q50")
    q90 <- ps %>% filter(q_key == "fit_q90")

    html <- paste0(
      '<p>Given a wind speed of <strong>', input$wind, ' km/h</strong>, ',
      'a relative humidity of <strong>', input$rh, '%</strong>, ',
      'and a temperature of <strong>', input$temp, '&deg;C</strong>, we would expect:</p>',
      '<ul style="padding-left:20px;">',
      '<li style="margin-bottom:10px;"><strong>50% of fires</strong> under similar conditions to spread faster than ',
      '<strong>', sprintf("%.2f", q50$median), ' km/h</strong> ',
      '<span style="color:#777;">(with some uncertainty - this could be as low as ',
      sprintf("%.2f", q50$lo), ' or as high as ', sprintf("%.2f", q50$hi), ' km/h)</span></li>',
      '<li><strong>10% of fires</strong> under similar conditions to spread faster than ',
      '<strong>', sprintf("%.2f", q90$median), ' km/h</strong> ',
      '<span style="color:#777;">(with some uncertainty - this could be as low as ',
      sprintf("%.2f", q90$lo), ' or as high as ', sprintf("%.2f", q90$hi), ' km/h)</span></li>',
      '</ul>'
    )

    HTML(paste0(
      '<div style="margin-top:8px; padding:14px 18px; background:#f0f4f8; ',
      'border-left:5px solid #2166ac; border-radius:4px; font-size:1.05em; line-height:1.6;">',
      html, '</div>'
    ))
  })

  # ── Model Comparison cards ──────────────────────────────────
  output$modelComparisonCards <- renderUI({
    ps  <- pred_summary()
    q50 <- ps %>% filter(q_key == "fit_q50")
    q90 <- ps %>% filter(q_key == "fit_q90")
    mc  <- model_compare_calc()

    fluidRow(
      value_box(
        title    = "Median ROS (q50)",
        value    = paste0(sprintf("%.2f", q50$median), " km/h"),
        subtitle = paste0("80% CI: ", sprintf("%.2f", q50$lo), "–", sprintf("%.2f", q50$hi), " km/h"),
        note     = "Bayesian model - 50% of fires expected to spread faster than this.",
        icon_name = "fire",
        color    = q_cols[["fit_q50"]],
        width    = 3
      ),
      value_box(
        title    = "Extreme ROS (q90)",
        value    = paste0(sprintf("%.2f", q90$median), " km/h"),
        subtitle = paste0("80% CI: ", sprintf("%.2f", q90$lo), "–", sprintf("%.2f", q90$hi), " km/h"),
        note     = "Bayesian model - 10% of fires expected to spread faster than this.",
        icon_name = "triangle-exclamation",
        color    = q_cols[["fit_q90"]],
        width    = 3
      ),
      value_box(
        title    = "McArthur Mk5",
        value    = paste0(sprintf("%.2f", mc$mcarthur$ros), " km/h"),
        subtitle = paste0("FFDI: ", sprintf("%.1f", mc$mcarthur$ffdi)),
        note     = "Forest fire danger meter (McArthur 1973/CSIRO).",
        letter   = "M",
        color    = q_cols_compare[["mcarthur"]],
        width    = 3
      ),
      value_box(
        title    = "Vesta Mk2",
        value    = paste0(sprintf("%.2f", mc$vesta$ros), " km/h"),
        subtitle = paste0("Fuel moisture: ", sprintf("%.1f", mc$vesta$mc), "%"),
        note     = "Dry eucalypt forest model",
        letter   = "V",
        color    = q_cols_compare[["vesta"]],
        width    = 3
      )
    )
  })

  # ── ROS Map (Bayesian q50/q90 ellipses) ────────────────────
  ellipse_sf_react <- reactive({
    ps  <- pred_summary()
    ign <- ignition()

    pmap_dfr(ps, function(q_key, median, lo, hi, ...) {
      qk       <- as.character(q_key)
      poly_hi  <- km_to_polygon(fire_ellipse_pts(hi,     input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)
      poly_med <- km_to_polygon(fire_ellipse_pts(median, input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)
      poly_lo  <- km_to_polygon(fire_ellipse_pts(lo,     input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)

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
      setView(lng = IGN_LON_DEFAULT, lat = IGN_LAT_DEFAULT, zoom = 13) %>%
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
    ign     <- ignition()
    kpdlon  <- km_per_deg_lon(ign$lat)
    dir_rad <- input$spread_dir * pi / 180
    fwd_lon <- sin(dir_rad) / kpdlon
    fwd_lat <- cos(dir_rad) / KM_PER_DEG_LAT
    rgt_lon <- cos(dir_rad) / kpdlon
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

      tip_lon <- ign$lon + head_km * fwd_lon
      tip_lat <- ign$lat + head_km * fwd_lat
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
      addMarkers(lng = ign$lon, lat = ign$lat, label = "Ignition",
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

  # ── Model Comparison Map (Bayesian + McArthur + Vesta ellipses) ─
  ellipse_sf_compare_react <- reactive({
    ps  <- pred_summary()
    ign <- ignition()
    mc  <- model_compare_calc()

    bayes_sf <- pmap_dfr(ps, function(q_key, median, lo, hi, ...) {
      qk       <- as.character(q_key)
      poly_hi  <- km_to_polygon(fire_ellipse_pts(hi,     input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)
      poly_med <- km_to_polygon(fire_ellipse_pts(median, input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)
      poly_lo  <- km_to_polygon(fire_ellipse_pts(lo,     input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat)

      ring <- tryCatch(
        st_difference(st_sfc(poly_hi, crs = 4326), st_sfc(poly_lo, crs = 4326)),
        error = function(e) st_sfc(poly_hi, crs = 4326)
      )

      bind_rows(
        st_sf(q_key = qk, band = "ring",   geometry = ring),
        st_sf(q_key = qk, band = "median", geometry = st_sfc(poly_med, crs = 4326))
      )
    })

    model_sf <- bind_rows(
      st_sf(q_key = "mcarthur", band = "single",
            geometry = st_sfc(
              km_to_polygon(fire_ellipse_pts(mc$mcarthur$ros, input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat),
              crs = 4326)),
      st_sf(q_key = "vesta", band = "single",
            geometry = st_sfc(
              km_to_polygon(fire_ellipse_pts(mc$vesta$ros, input$time_hr, input$lb_ratio, input$spread_dir), ign$lon, ign$lat),
              crs = 4326))
    )

    bind_rows(bayes_sf, model_sf) %>%
      mutate(q_key = factor(q_key, levels = names(q_cols_compare)))
  })

  output$ellipseMapCompare <- renderLeaflet({
    leaflet() %>%
      addTiles(attribution = "© OpenStreetMap contributors") %>%
      setView(lng = IGN_LON_DEFAULT, lat = IGN_LAT_DEFAULT, zoom = 13) %>%
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
    ed    <- ellipse_sf_compare_react()
    ign   <- ignition()
    mc    <- model_compare_calc()
    proxy <- leafletProxy("ellipseMapCompare") %>% clearShapes() %>% clearMarkers()

    for (qk in c("fit_q50", "fit_q90")) {
      col   <- q_cols_compare[[qk]]
      proxy <- proxy %>%
        addPolygons(data = ed %>% filter(q_key == qk, band == "ring"),
                    fillColor = col, fillOpacity = 0.25,
                    color = col, weight = 1.0,
                    label = glue("{qk} — 80% CI band")) %>%
        addPolygons(data = ed %>% filter(q_key == qk, band == "median"),
                    fillOpacity = 0, color = col, weight = 2.5,
                    label = glue("{qk} — median"))
    }

    model_labels <- list(
      mcarthur = glue("McArthur Mk5 — {sprintf('%.2f', mc$mcarthur$ros)} km/h"),
      vesta    = glue("Vesta Mk2 — {sprintf('%.2f', mc$vesta$ros)} km/h")
    )

    for (qk in c("mcarthur", "vesta")) {
      col   <- q_cols_compare[[qk]]
      proxy <- proxy %>%
        addPolygons(data = ed %>% filter(q_key == qk, band == "single"),
                    fillOpacity = 0,
                    color = col, weight = 5, dashArray = q_dash_compare[[qk]],
                    label = model_labels[[qk]])
    }

    proxy <- proxy %>%
      addMarkers(lng = ign$lon, lat = ign$lat, label = "Ignition",
                 popup = glue("Ignition — spread direction: {input$spread_dir}°"))

    if (nrow(ed) > 0 && any(!st_is_empty(ed))) {
      bbox <- tryCatch(st_bbox(ed), error = function(e) NULL)
      if (!is.null(bbox)) {
        pad <- max((bbox[["xmax"]] - bbox[["xmin"]]) * 0.15, 0.003)
        proxy %>% fitBounds(
          lng1 = bbox[["xmin"]] - pad, lat1 = bbox[["ymin"]] - pad,
          lng2 = bbox[["xmax"]] + pad, lat2 = bbox[["ymax"]] + pad
        )
      }
    }
  })

  output$ellipseSummaryCompare <- renderUI({
    ps  <- pred_summary()
    mc  <- model_compare_calc()
    q50 <- ps %>% filter(q_key == "fit_q50")
    q90 <- ps %>% filter(q_key == "fit_q90")

    bayes_rows <- pmap_chr(ps, function(q_key, label, median, lo, hi, ...) {
      col     <- q_cols_compare[[as.character(q_key)]]
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

    model_rows <- c(
      glue(
        '<span style="color:{q_cols_compare[["mcarthur"]]}; font-weight:bold;">McArthur Mk5</span> &nbsp;|&nbsp; ',
        'ROS: {sprintf("%.2f", mc$mcarthur$ros)} km/h ',
        '&nbsp;|&nbsp; Head in {input$time_hr} h: <strong>{round(mc$mcarthur$ros * input$time_hr, 2)} km</strong>'
      ),
      glue(
        '<span style="color:{q_cols_compare[["vesta"]]}; font-weight:bold;">Vesta Mk2</span> &nbsp;|&nbsp; ',
        'ROS: {sprintf("%.2f", mc$vesta$ros)} km/h ',
        '&nbsp;|&nbsp; Head in {input$time_hr} h: <strong>{round(mc$vesta$ros * input$time_hr, 2)} km</strong>'
      )
    )

    HTML(paste0(
      '<div style="margin-top:10px; padding:10px 16px; background:#f0f4f8; ',
      'border-left:5px solid #2166ac; border-radius:4px; font-size:1.0em; line-height:2.6;">',
      paste(c(bayes_rows, model_rows), collapse = "<br/>"), '</div>'
    ))
  })

}

# ============================================================
# --- RUN ---
# ============================================================

shinyApp(ui, server) |>
  shiny::runApp(host = Sys.getenv("SHINY_HOST", "127.0.0.1"),
                port = as.integer(Sys.getenv("SHINY_PORT", "5826")))
