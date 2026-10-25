# UOW Fire Apps

A collection of [R Shiny](https://shiny.posit.co/) apps for mapping and analysing bushfire spread, developed at the University of Wollongong. The apps cover three aspects of rate of spread analysis:

1. **Map** fire perimeters from airborne line scan imagery (`app-sam-fire`)
2. **Measure** rate of spread (ROS) between successive perimeters (`app-measure-spread`)
3. **Predict** rate of spread from weather and landscape conditions using a Bayesian quantile model (`app-predict-spread`)

| App | What it does | Main inputs | Output |
|---|---|---|---|
| [`app-sam-fire`](app-sam-fire) | Semi-automated fire perimeter mapping from line scans using SAM2 image segmentation | 3-band GeoTIFF line scans | Fire perimeter polygons (GeoPackage) |
| [`app-measure-spread`](app-measure-spread) | Measures spread lines and rates of spread between timed fire perimeters | Timestamped perimeter polygons (GeoPackage or zipped shapefile) | Spread lines with distance, time step and ROS (GeoPackage) |
| [`app-predict-spread`](app-predict-spread) | Predicts ROS from a Bayesian quantile model and compares it with McArthur Mk5 and Vesta Mk2 | Weather, soil moisture, atmospheric stability and fuel inputs (sliders) | ROS estimates with uncertainty, plots and spread ellipses on a map |

---

## The apps

### `app-sam-fire`: Semi-automated line scan mapping

Maps fire perimeters from airborne line scan imagery. It uses Meta's [Segment Anything Model 2 (SAM2)](https://github.com/facebookresearch/sam2) through the [`fireobs`](https://github.com/mstorey87/fireobs) R package.
The tool demonstrates how the Segment Anything model can be applied specifically to the issue of mapping active fires from images, to save manual mapping time.
The process can captured the edges of actively burning area and burnt area.
The process requires the combined automated component with manual oversight and correction. 
Some scans may only need minimal correction, some may require more depending on if the image is washed out, how bright the burning area is, the contrast between burnt and unburnt area, cloud cover etc.
Meta releases new versions of the Segment Anything model (version 3 is available as at Oct 2026).
The accuracy and speed could be improved by adopting new versions, using a GPU enabled computer and conducting a specific model training process on line scan data.

Workflow:

1. **Upload** one or more 3-band GeoTIFF line scans. File names should start with a timestamp in `YYYYmmddHHMMSS` format so the polygon gets the right date and time, although the process should still work if that is not done.
2. **Preview** the scans on a map. In cases where multiple scans are uploaded to cover a single fire area, drag them to set the merge order (the top scan is on top, and its time and name are used for the polygon). Untick any scans you don't need, then click **Send to mapping tool**.
3. **Create polygons.** Draw a series of points *or* polygons (not both) over the fire area, then click **Auto-map** to run the sam2 model, which produces the polygon. Refine the result with:
   - **Keep Selected**: keep only the polygons you click on
   - **Cut**: remove the area under a manually drawn polygon (ie a polygon you draw over the auto-generated polygon)
   - **Union**: add a drawn area to the fire polygon
   - **Fill Holes**: fill small gaps automatically
   - **Undo** / **Restart**
4. **Download** the result as a GeoPackage, including the date/time and the scans used.

Scans are downsampled (5× aggregation) to keep segmentation fast. Uploads of up to 500 MB are allowed.

### `app-measure-spread`: Measure spread lines and ROS

Measures how far and how fast a fire spread between successive perimeters.

1. **Upload** fire progression polygons as a GeoPackage (`.gpkg`) or zipped shapefile (`.zip`). The layer needs a posix date/time column for each perimeter.
2. **Choose** the date/time column and its time zone. The app checks that the values parse.
3. **Draw** a polygon on the map around the perimeters you want to analyse (at least two time steps are needed). This is so the user can focus on a specific location in the context of a larger fire.
4. Choose whether to keep only the **maximum spread line** per time step, and whether to use **convex hulls** of the perimeters.
5. Click **Measure Spread Lines**. The lines are drawn on the map and listed in a results table.
6. **Download** the spread lines as a GeoPackage.

Spread lines are calculated with `fireobs::fire_max_spread_line()` in GDA94 / Geoscience Australia Lambert (EPSG:3112). This package is available on github. Contact Michael Storey.

### `app-predict-spread`: Fire ROS prediction

Interactive predictions from a Bayesian quantile regression model of forward rate of spread, fitted with `brms` to observed Australian fire spread data. The model structure is:

```
ROS ~ 0 + pyroregion + t2(wind, FMI) + s(soil moisture) + s(lapse rate) + s(WBT SD slope)
```

It is fitted for two quantiles: the **median (q50)** and an **extreme (q90)** rate of spread.

Inputs:

- **Region:** pyroregion
- **Weather:** surface wind speed, temperature and relative humidity. These give the Fuel Moisture Index (`FMI = 10 − 0.25 × (T − RH)`), which can also be set by hand.
- **Soil, stability and terrain:** 0–10 cm soil moisture, 850–700 hPa temperature difference (lapse rate) and wet-bulb temperature SD slope
- **Comparison models:** drought factor, fine fuel load, wind adjustment factor, understorey height, slope and time of day

Tabs:

- **ROS Plot:** predicted ROS across wind speeds for both quantiles, with credible intervals
- **ROS Map:** spread ellipses for a chosen ignition point, spread period, direction and length-to-breadth ratio
- **ROS Description:** summary of the prediction
- **Model Comparison:** the Bayesian predictions next to the [McArthur Mk5](https://research.csiro.au/spark/resources/model-library/mcarthur/) and [Vesta Mk2](https://research.csiro.au/spark/resources/model-library/dry-eucalypt-forest-vesta-mk-2/) models

The deployed app does **not** need Stan or `brms`. `precompute_lookup.R` evaluates each additive model term on a fine grid for 200 posterior draws and saves it to `stan_models/prediction_lookup.rds`. The app rebuilds predictions by interpolating and summing these terms. `prediction_app_explainer.Rmd` walks through the method step by step.

---

## Installing and running locally

All apps need **R (≥ 4.3 recommended)**. The spatial packages (`sf`, `terra`) need GDAL, GEOS and PROJ. On Windows and macOS these come with the CRAN binary packages. On Linux, install them first, for example on Ubuntu/Debian:

```bash
sudo apt install libgdal-dev gdal-bin libproj-dev libgeos-dev libudunits2-dev libcurl4-openssl-dev libssl-dev libxml2-dev
```

Clone the repository:

```bash
git clone https://github.com/mstorey87/uow_fire_apps.git
```

### `app-predict-spread` (simplest)

Install the R packages:

```r
install.packages(c("shiny", "ggplot2", "dplyr", "purrr", "glue", "scales",
                   "htmltools", "leaflet", "sf", "bsplus"))
```

Then run from the app folder so the lookup table in `stan_models/` is found:

```r
setwd("uow_fire_apps/app-predict-spread")
shiny::runApp()
```

Settings can be changed with environment variables:

| Variable | Default | Purpose |
|---|---|---|
| `STAN_DIR` | `stan_models` | Folder containing `prediction_lookup.rds` |
| `SHINY_HOST` | `127.0.0.1` | Host to bind to (only used when run with `Rscript app.R`) |
| `SHINY_PORT` | `5826` | Port (only used when run with `Rscript app.R`) |

Rebuilding the lookup table (only needed if the fitted models change) also needs `brms`/`rstan` and the fitted model objects. Set `SOURCE_STAN_DIR` at the top of `precompute_lookup.R` to point to them.

### `app-measure-spread`

```r
install.packages(c("shiny", "leaflet", "leaflet.extras", "sf", "dplyr", "lubridate", "remotes"))
remotes::install_github("mstorey87/fireobs")

shiny::runApp("uow_fire_apps/app-measure-spread")
```

### `app-sam-fire` (needs Python + SAM2)

This app calls SAM2 from Python through `reticulate`, so it needs a Python environment as well as R.

**1. R packages**

```r
install.packages(c("shiny", "leaflet", "leaflet.extras", "sortable", "leafem", "terra",
                   "sf", "dplyr", "bsplus", "reticulate", "geojsonsf", "jsonlite",
                   "lwgeom", "nngeo", "smoothr", "units", "purrr", "magick", "remotes"))
remotes::install_github("mstorey87/fireobs")
```

**2. Python environment with SAM2** (Python ≥ 3.10; conda or venv both work)

```bash
python -m venv sam2_env
```

Activate the environment, then:

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cpu
pip install git+https://github.com/facebookresearch/sam2.git
pip install opencv-python-headless
```

Use the CUDA build of PyTorch instead of the CPU build if you have a supported GPU.

**3. SAM2 checkpoint and config**

Download the model weights and config into a checkpoints folder, keeping this layout:

```
checkpoints/
├── sam2.1_hiera_base_plus.pt      # https://dl.fbaipublicfiles.com/segment_anything_2/092824/sam2.1_hiera_base_plus.pt
└── sam2.1/
    └── sam2.1_hiera_b+.yaml       # from sam2/configs/sam2.1/ in the SAM2 GitHub repo
```

**4. Point the app at your local paths**

At the top of `app-sam-fire/app.R`, the paths are set for the Docker deployment (`/srv/venv/...`, `/srv/tmp`). For a local install, swap them for your own Python environment, checkpoints folder and a writable working directory. The file has a commented-out example of this:

```r
reticulate::use_condaenv(condaenv = "my_sam2_env")   # or reticulate::use_virtualenv("path/to/sam2_env")
fireobs:::fire_load_sam_once("C:/path/to/checkpoints",
                                   pt_name   = "sam2.1_hiera_base_plus.pt",
                                   yaml_name = "sam2.1_hiera_b+.yaml")
wk_dir <- "C:/temp/wk_dir"
```

**5. Run**

```r
shiny::runApp("uow_fire_apps/app-sam-fire")
```

Segmentation runs on the CPU by default and can take a while on large scans. On a multi-core machine, setting `OMP_NUM_THREADS` / `TORCH_NUM_THREADS` can make it more responsive.

---

## Repository layout

```
uow_fire_apps/
├── app-measure-spread/
│   └── app.R
├── app-predict-spread/
│   ├── app.R
│   ├── precompute_lookup.R              # builds the lookup table from fitted brms models
│   ├── prediction_app_explainer.Rmd     # walkthrough of the prediction method
│   └── stan_models/prediction_lookup.rds
├── app-sam-fire/
│   └── app.R
└── uow_fire_apps.Rproj
```

## Related

- [`fireobs`](https://github.com/mstorey87/fireobs): R package with the fire progression mapping, SAM2 wrapper and spread-line functions these apps use.

## Contact

**M Storey**, University of Wollongong
GitHub: [@mstorey87](https://github.com/mstorey87)

For questions, bug reports or collaboration, please open an [issue](https://github.com/mstorey87/uow_fire_apps/issues) or contact M Storey directly.
