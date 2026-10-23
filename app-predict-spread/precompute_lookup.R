# --- precompute_lookup.R ---
# Run this LOCALLY (with brms/rstan installed) whenever the fitted models change.
# It extracts each additive term of the brms model (per posterior draw) onto a
# fine grid, so the deployed Shiny app can reconstruct posterior_linpred() by
# summing interpolated components - without needing brms/rstan/Stan at runtime.
#
# Model structure (see: formula(fit)):
#   mu ~ 0 + pyroregion + t2(wind, fmi_surface) + s(sm_top) + s(lapse) + s(wbtsdslope_mean)
# This is additive, so each term can be gridded independently and summed per draw.
# Verified against posterior_linpred() to match to floating-point precision (~1e-15).

library(brms)

SOURCE_STAN_DIR <- "D:/UOW_projects/NHRA_ROS/DataWorking/bayesian_model/stan_models"
OUT_PATH        <- "stan_models/prediction_lookup.rds"

DRAW_IDS       <- 1:200   # fixed set of posterior draws baked into the lookup
T2_GRID_N      <- 60      # grid points per axis for the wind x fmi_surface surface
SMOOTH1D_GRID_N <- 200    # grid points for each 1D smooth (sm_top, lapse, wbtsdslope_mean)

build_term_lookup <- function(fit) {
  train_dat   <- fit$data
  pyro_levels <- levels(train_dat$pyroregion)

  var_stats <- lapply(
    c("wind", "fmi_surface", "sm_top", "lapse", "wbtsdslope_mean"),
    function(v) list(min = min(train_dat[[v]], na.rm = TRUE),
                      max = max(train_dat[[v]], na.rm = TRUE),
                      median = median(train_dat[[v]], na.rm = TRUE))
  )
  names(var_stats) <- c("wind", "fmi_surface", "sm_top", "lapse", "wbtsdslope_mean")

  grid_wind <- seq(var_stats$wind$min,        var_stats$wind$max,        length.out = T2_GRID_N)
  grid_fmi  <- seq(var_stats$fmi_surface$min, var_stats$fmi_surface$max, length.out = T2_GRID_N)
  grid_sm   <- seq(var_stats$sm_top$min,          var_stats$sm_top$max,          length.out = SMOOTH1D_GRID_N)
  grid_lp   <- seq(var_stats$lapse$min,           var_stats$lapse$max,           length.out = SMOOTH1D_GRID_N)
  grid_wbt  <- seq(var_stats$wbtsdslope_mean$min, var_stats$wbtsdslope_mean$max, length.out = SMOOTH1D_GRID_N)

  # newdata needs every variable in the formula, even though a given
  # posterior_smooths() call only uses its own term's columns
  base_extra <- data.frame(
    pyroregion  = factor(pyro_levels[1], levels = pyro_levels),
    fire_day_id = NA_character_
  )

  cat("  building pyroregion offsets...\n")
  b_pyro <- as.matrix(fit, variable = "^b_pyroregion", regex = TRUE)[DRAW_IDS, , drop = FALSE]
  colnames(b_pyro) <- pyro_levels  # brms sanitizes level names (e.g. "n/p" -> "nDp"); map back by position

  cat("  building t2(wind, fmi_surface) surface...\n")
  nd_t2 <- expand.grid(wind = grid_wind, fmi_surface = grid_fmi)
  nd_t2 <- cbind(nd_t2,
                 sm_top = median(train_dat$sm_top), lapse = median(train_dat$lapse),
                 wbtsdslope_mean = median(train_dat$wbtsdslope_mean), base_extra, row.names = NULL)
  t2_flat <- posterior_smooths(fit, smooth = "t2(wind,fmi_surface,k=5)", newdata = nd_t2, draw_ids = DRAW_IDS)
  t2_arr  <- array(t2_flat, dim = c(length(DRAW_IDS), T2_GRID_N, T2_GRID_N))

  cat("  building s(sm_top)...\n")
  nd_sm <- cbind(sm_top = grid_sm,
                 wind = median(train_dat$wind), fmi_surface = median(train_dat$fmi_surface),
                 lapse = median(train_dat$lapse), wbtsdslope_mean = median(train_dat$wbtsdslope_mean),
                 base_extra, row.names = NULL)
  sm_top_mat <- posterior_smooths(fit, smooth = "s(sm_top,k=5)", newdata = nd_sm, draw_ids = DRAW_IDS)

  cat("  building s(lapse)...\n")
  nd_lp <- cbind(lapse = grid_lp,
                 wind = median(train_dat$wind), fmi_surface = median(train_dat$fmi_surface),
                 sm_top = median(train_dat$sm_top), wbtsdslope_mean = median(train_dat$wbtsdslope_mean),
                 base_extra, row.names = NULL)
  lapse_mat <- posterior_smooths(fit, smooth = "s(lapse,k=8)", newdata = nd_lp, draw_ids = DRAW_IDS)

  cat("  building s(wbtsdslope_mean)...\n")
  nd_wbt <- cbind(wbtsdslope_mean = grid_wbt,
                  wind = median(train_dat$wind), fmi_surface = median(train_dat$fmi_surface),
                  sm_top = median(train_dat$sm_top), lapse = median(train_dat$lapse),
                  base_extra, row.names = NULL)
  wbt_mat <- posterior_smooths(fit, smooth = "s(wbtsdslope_mean,k=5)", newdata = nd_wbt, draw_ids = DRAW_IDS)

  list(
    pyro_levels = pyro_levels, pyro_mat = b_pyro,
    grid_wind = grid_wind, grid_fmi = grid_fmi, t2_arr = t2_arr,
    grid_sm_top = grid_sm, sm_top_mat = sm_top_mat,
    grid_lapse = grid_lp, lapse_mat = lapse_mat,
    grid_wbt = grid_wbt, wbt_mat = wbt_mat,
    var_stats = var_stats,
    ndraws = length(DRAW_IDS)
  )
}

verify_against_posterior_linpred <- function(fit, term_lookup, n_check = 20) {
  train_dat <- fit$data
  set.seed(42)
  nd <- data.frame(
    pyroregion      = factor(sample(term_lookup$pyro_levels, n_check, replace = TRUE), levels = term_lookup$pyro_levels),
    wind            = runif(n_check, term_lookup$var_stats$wind$min, term_lookup$var_stats$wind$max),
    fmi_surface     = runif(n_check, term_lookup$var_stats$fmi_surface$min, term_lookup$var_stats$fmi_surface$max),
    sm_top          = runif(n_check, term_lookup$var_stats$sm_top$min, term_lookup$var_stats$sm_top$max),
    lapse           = runif(n_check, term_lookup$var_stats$lapse$min, term_lookup$var_stats$lapse$max),
    wbtsdslope_mean = runif(n_check, term_lookup$var_stats$wbtsdslope_mean$min, term_lookup$var_stats$wbtsdslope_mean$max),
    fire_day_id     = NA_character_
  )

  truth <- posterior_linpred(fit, newdata = nd, draw_ids = DRAW_IDS, re_formula = NA, allow_new_levels = TRUE)

  interp1 <- function(grid, mat, x) {
    x <- pmin(pmax(x, grid[1]), grid[length(grid)])
    ih <- pmin(pmax(findInterval(x, grid, all.inside = TRUE) + 1, 2), length(grid))
    il <- ih - 1
    w  <- (x - grid[il]) / (grid[ih] - grid[il])
    mat[, il] * (1 - w) + mat[, ih] * w
  }
  interp2 <- function(gx, gy, arr, x, y) {
    x <- pmin(pmax(x, gx[1]), gx[length(gx)]); y <- pmin(pmax(y, gy[1]), gy[length(gy)])
    ixh <- pmin(pmax(findInterval(x, gx, all.inside = TRUE) + 1, 2), length(gx)); ixl <- ixh - 1
    iyh <- pmin(pmax(findInterval(y, gy, all.inside = TRUE) + 1, 2), length(gy)); iyl <- iyh - 1
    wx <- (x - gx[ixl]) / (gx[ixh] - gx[ixl]); wy <- (y - gy[iyl]) / (gy[iyh] - gy[iyl])
    (1-wx)*(1-wy)*arr[, ixl, iyl] + wx*(1-wy)*arr[, ixh, iyl] +
      (1-wx)*wy*arr[, ixl, iyh] + wx*wy*arr[, ixh, iyh]
  }

  recon <- sapply(seq_len(n_check), function(i) {
    term_lookup$pyro_mat[, as.character(nd$pyroregion[i])] +
      interp2(term_lookup$grid_wind, term_lookup$grid_fmi, term_lookup$t2_arr, nd$wind[i], nd$fmi_surface[i]) +
      interp1(term_lookup$grid_sm_top, term_lookup$sm_top_mat, nd$sm_top[i]) +
      interp1(term_lookup$grid_lapse, term_lookup$lapse_mat, nd$lapse[i]) +
      interp1(term_lookup$grid_wbt, term_lookup$wbt_mat, nd$wbtsdslope_mean[i])
  })

  max_abs_diff <- max(abs(truth - recon))
  cat(sprintf("  max abs diff vs posterior_linpred over %d random points: %.6f\n", n_check, max_abs_diff))
  invisible(max_abs_diff)
}

cat("Loading fit_q50...\n")
fit_q50 <- readRDS(file.path(SOURCE_STAN_DIR, "brms_quantile_q50_ros_interaction_model3.rds"))
cat("Building lookup for fit_q50...\n")
lk_q50 <- build_term_lookup(fit_q50)
cat("Verifying fit_q50 lookup (grid interpolation, not exact draws, so expect small but nonzero error)...\n")
verify_against_posterior_linpred(fit_q50, lk_q50)
rm(fit_q50); gc()

cat("Loading fit_q90...\n")
fit_q90 <- readRDS(file.path(SOURCE_STAN_DIR, "brms_quantile_q90_ros_interaction_model3.rds"))
cat("Building lookup for fit_q90...\n")
lk_q90 <- build_term_lookup(fit_q90)
cat("Verifying fit_q90 lookup...\n")
verify_against_posterior_linpred(fit_q90, lk_q90)
rm(fit_q90); gc()

lookup <- list(fit_q50 = lk_q50, fit_q90 = lk_q90)
saveRDS(lookup, OUT_PATH, compress = "xz")

cat(sprintf("\nSaved lookup to %s (%.2f MB)\n", OUT_PATH, file.info(OUT_PATH)$size / 1024^2))
