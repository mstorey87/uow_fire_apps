# UOW Fire Apps



A collection of [R Shiny](https://shiny.posit.co/) apps for mapping and analysing bushfire spread, developed at the University of Wollongong. The apps cover three aspects of rate of spread analysis:

1. **Map** fire perimeters from airborne line scan imagery (`app-sam-fire`)
2. **Measure** rate of spread (ROS) between successive perimeters (`app-measure-spread`)
3. **Predict** rate of spread from weather and landscape conditions using a Bayesian quantile model (`app-predict-spread`)

| App | What it does | Main inputs | Output |
|---|---|---|---|
| [`app-sam-fire`](app-sam-fire) | Semi-automated fire perimeter mapping from line scans using SAM2 image segmentation | 3-band GeoTIFF line scans | Fire perimeter polygons (GeoPackage) |
| [`app-measure-spread`](app-measure-spread) | Measures spread lines and rates of spread between timed fire perimeters | Timestamped perimeter polygons (GeoPackage or zipped shapefile) | Spread lines with distance, time step and ROS (GeoPackage) |
| [`app-predict-spread`](app-predict-spread) | Predicts ROS from a Bayesian quantile model and compares it with McArthur Mk5 and Vesta Mk2 | Weather, soil moisture, atmospheric stability and fuel inputs (sliders) | ROS estimates with uncertainty, plots and spread ellipses on a map | This is a work in progress |

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

1. **Upload** one or more 3-band GeoTIFF line scans. File names should start with a timestamp in `YYYYmmddHHMMSS` format so the polygon gets the right date and time, although the process should still work if that is not done - time can be entered manually.
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

Measures how far and how fast a fire spread between successive perimeters. The methods is the same that is used in "Derivation of a Bayesian fire spread model using large-scale wildfire observations" (https://doi.org/10.1016/j.envsoft.2021.105127)

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

Get a copy of the repository, either with git:

```bash
git clone https://github.com/mstorey87/uow_fire_apps.git
```

or without git: on the GitHub page click **Code → Download ZIP** and unzip it.

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

This app runs SAM2 in Python, called from R through `reticulate`. Everything it needs apart from R, RStudio and Python is installed by [`app-sam-fire/install_sam_app.Rmd`](app-sam-fire/install_sam_app.Rmd).

**1. Install R, RStudio and Python yourself first**

| Software | Where to get it | Notes |
|---|---|---|
| R (≥ 4.3) | <https://cran.r-project.org/> | |
| RStudio Desktop | <https://posit.co/download/rstudio-desktop/> | Used to run the install document and the app |
| Python 3.10–3.12 | <https://www.python.org/downloads/> | 3.11 or 3.12 recommended. Avoid the newest Python release, because PyTorch and other packages can take a few months to support it |

**You need to know the path to your Python.** The install document asks for it. Typical locations:

- **Windows:** `C:/Users/<you>/AppData/Local/Programs/Python/Python311/python.exe`. If Python is on your PATH, `where python` in a Command Prompt shows it.
- **Mac/Linux:** `/usr/bin/python3` or `/usr/local/bin/python3`. `which python3` in a terminal shows it.

Git is not needed.

**2. Run the install document**

1. Open `app-sam-fire/install_sam_app.Rmd` in RStudio.
2. Restart R (**Session → Restart R**). The setup only works in a session that hasn't used Python yet.
3. In the **Settings** chunk at the top:
   - set `python_path` to your Python path, using forward slashes (or leave it as `NULL` to let it search for one);
   - optionally change `torch_build` to a CUDA build (e.g. `"cu126"`) if you have an NVIDIA GPU;
   - optionally change where the model weights (`checkpoints_dir`) and temporary files (`wk_dir`) go.
4. Click **Run All** (Ctrl+Alt+R), or **Knit**.

It takes 10–30 minutes the first time. It:

- installs the R packages, including `fireobs` from GitHub;
- creates a Python virtual environment (`sam2_env`) and installs PyTorch, OpenCV and SAM2 into it;
- checks that Python and SAM2 load, and reports whether a GPU was found;
- downloads the SAM2 model weights (about 300 MB) and loads the model as a final test.

**3. Point the app at your setup**

At the end, the install document prints three lines (`reticulate::use_python(...)`, `wk_dir` and `chkp_dir`). Paste them over the matching lines near the top of `app-sam-fire/app.R`. The ones there now are for the Docker deployment (`/srv/venv/...`, `/srv/tmp`).

**4. Run the app**

Open `app-sam-fire/app.R` in RStudio and click **Run App**, or run the `shiny::runApp(...)` command printed at the end of the install document. It uses the full path to the app folder, so it works whatever R's working directory is.

Segmentation runs on the CPU unless you installed a CUDA build of PyTorch and have an NVIDIA GPU. On the CPU it can take a while on large scans.

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
│   ├── app.R
│   └── install_sam_app.Rmd              # installs everything the app needs
└── uow_fire_apps.Rproj
```

## Related

- [`fireobs`](https://github.com/mstorey87/fireobs): R package with the fire progression mapping, SAM2 wrapper and spread-line functions these apps use.

## Contact

**M Storey**, University of Wollongong
GitHub: [@mstorey87](https://github.com/mstorey87)

For questions, bug reports or collaboration, please open an [issue](https://github.com/mstorey87/uow_fire_apps/issues) or contact M Storey directly.
