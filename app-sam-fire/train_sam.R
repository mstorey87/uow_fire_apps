library(sf)
library(dplyr)
library(terra)
library(purrr)
library(png)
library(reticulate)

# =============================================================================
# STEP 1: Load & filter fire progression data
# =============================================================================
# All fire progression polygons from the curated geopackage
dat.progs <- st_read(
  "D:\\UOW_projects\\NHRA_ROS\\DataWorking\\package_development\\wfprogression_otherdata\\vector_data_export_v31.gpkg"
)

# Hand-picked best linescan scans for training (filenames without extension)
list.best <- list.files("C:\\app_dev\\train_sam\\best_scans\\") %>%
  tibble(nm = basename(.)) %>%
  mutate(nm = purrr::map_chr(stringr::str_split(nm, "___"), 1))

# Keep only single-scan linescan records that match the curated list
dat.progs.filt <- dat.progs %>%
  filter(
    sourcetype == "linescan",
    !stringr::str_detect(s3name, ";"),   # exclude multi-scan entries
    s3name %in% list.best$nm
  ) %>%
  mutate(
    s3name, sourcetype, dt_local,
    tif_path = paste0("D:\\ROS_database_all_noreproject_v2\\", s3name, ".tif"),
    geom,
    .keep = "none"
  )

# Each element of poly_list = one fire-front polygon group (one scan x timestamp)
poly_list <- dat.progs.filt %>%
  group_by(dt_local, s3name) %>%
  group_split()# %>%
  #head(10)   # <- remove or increase for a full training run

scans_list <- purrr::map(poly_list, ~ .x$dt_local[1])

# =============================================================================
# STEP 2: Crop TIF tiles and export as matched PNG pairs (image + mask)
# =============================================================================
# PNG is the simplest format for Python/PIL - no geospatial library needed.
# Images: 3-channel, normalised to [0, 1] before writing (PIL reads back as uint8).
# Masks : grayscale, 0 = background / 255 = fire.

out_img  <- "C:/app_dev/train_sam/images"
out_mask <- "C:/app_dev/train_sam/masks"
dir.create(out_img,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_mask, recursive = TRUE, showWarnings = FALSE)

# Robust per-band min-max stretch to [0, 1], tolerates NA / Inf
norm01 <- function(x) {
  mm   <- terra::minmax(x)
  vmin <- mm[1, ];  vmax <- mm[2, ]
  if (all(is.na(vmin)) || all(is.na(vmax))) return(x * 0)
  rng <- vmax - vmin
  rng[rng == 0 | is.na(rng)] <- 1          # avoid divide-by-zero on flat bands
  y <- (x - vmin) / rng
  y[y < 0] <- 0;  y[y > 1] <- 1
  y
}

for (i in seq_along(poly_list)) {

  p_sf  <- poly_list[[i]]
  r_pth <- p_sf$tif_path[1]

  if (!file.exists(r_pth)) {
    message(sprintf("[%d/%d] Skipped - TIF not found: %s", i, length(poly_list), basename(r_pth)))
    next
  }

  r <- rast(r_pth)

  # Use first 3 bands; if single-band thermal scan, replicate across RGB
  img <- if (nlyr(r) >= 3) r[[1:3]] else c(r, r, r)
  if (all(is.na(values(img)))) next         # skip blank rasters

  img <- norm01(img)                         # scale each band to [0, 1]

  # Rasterise the fire-front polygon(s) onto the image grid -> binary mask
  p_proj <- st_transform(p_sf, crs(r))
  v      <- vect(p_proj)
  mask   <- rasterize(v, r[[1]], field = 1, background = 0)
  mask[is.na(mask)] <- 0

  if (terra::global(mask, "sum", na.rm = TRUE)[1, 1] == 0) next  # skip empty masks

  # Unique filename: scan name + timestamp
  name <- paste0(
    p_sf$s3name[1], "_",
    format(as.POSIXct(p_sf$dt_local[1]), "%Y%m%d_%H%M%S")
  )

  nr <- nrow(img);  nc <- ncol(img)

  # terra::values() returns cells row-by-row (row 1 left to right, then row 2 ...)
  # matrix(..., byrow = TRUE) reshapes correctly to H x W
  vals    <- terra::values(img, mat = TRUE)
  img_arr <- array(0, dim = c(nr, nc, 3))
  for (b in 1:3) {
    img_arr[, , b] <- matrix(vals[, b], nrow = nr, ncol = nc, byrow = TRUE)
  }
  img_arr[is.na(img_arr)] <- 0             # replace any residual NAs with black

  mask_vals <- terra::values(mask, mat = FALSE)
  mask_arr  <- matrix(mask_vals, nrow = nr, ncol = nc, byrow = TRUE)
  mask_arr[is.na(mask_arr)] <- 0

  # writePNG expects [0, 1]; mask values 0/1 map to black/white in the PNG file
  png::writePNG(img_arr,  file.path(out_img,  paste0(name, ".png")))
  png::writePNG(mask_arr, file.path(out_mask, paste0(name, ".png")))

  message(sprintf("[%d/%d] Saved: %s", i, length(poly_list), name))
}

message(sprintf(
  "PNG export complete.\n  Images : %s\n  Masks  : %s",
  out_img, out_mask
))

# =============================================================================
# STEP 3: Configure Python environment
# =============================================================================
reticulate::use_condaenv(condaenv = "my_sam2_env")

chkp_dir  <- "D:/temp/checkpoints"
yaml_name <- "configs/sam2.1/sam2.1_hiera_t.yaml"  # Tiny - fastest to fine-tune
pt_name   <- "sam2.1_hiera_tiny.pt"                 # swap for _small/_base_plus as needed
save_path <- file.path(chkp_dir, "sam2_finetuned.pt")

# =============================================================================
# STEP 4: Define Python helper functions
# =============================================================================
reticulate::py_run_string('
import os, torch, numpy as np
import torch.nn.functional as F
from PIL import Image
from torch.utils.data import Dataset, DataLoader


def get_box_from_mask(mask):
    """Return [[x1, y1, x2, y2]] (shape 1x4) covering the True region, or None."""
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    if not rows.any() or not cols.any():
        return None
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]
    return np.array([[cmin, rmin, cmax, rmax]], dtype=np.float32)  # [1, 4] xyxy


# -----------------------------------------------------------------------------
# Load SAM2.1 from a local checkpoint
# -----------------------------------------------------------------------------
def load_model(chkp_dir, yaml_name, pt_name):
    from sam2.build_sam import build_sam2
    ckpt_path = os.path.join(chkp_dir, pt_name)
    device    = "cuda" if torch.cuda.is_available() else "cpu"
    model = build_sam2(yaml_name, ckpt_path, device=device,
                       mode="eval", apply_postprocessing=False)
    print(f"Loaded {pt_name} on {device}")
    return model, device


# =============================================================================
# STEP 4b: Pre-compute and cache image embeddings (CPU speed-up)
# =============================================================================
# The image encoder is the slowest part (~3-5 sec/image on CPU).
# Since it is frozen, its output is identical every epoch.
# Running it once and saving to disk means each training epoch only
# runs the tiny decoder (~0.3 sec/image) - a ~10-20x speed-up.
# -----------------------------------------------------------------------------
def cache_embeddings(model, img_dir, cache_dir):
    """
    Run every image through the frozen encoder once and save the result.
    Skips images that already have a cached .pt file (safe to re-run).
    """
    from sam2.sam2_image_predictor import SAM2ImagePredictor
    os.makedirs(cache_dir, exist_ok=True)
    predictor = SAM2ImagePredictor(model)
    files = sorted(f for f in os.listdir(img_dir) if f.endswith(".png"))
    print(f"Caching embeddings for {len(files)} images -> {cache_dir}")

    for i, fname in enumerate(files):
        out_path = os.path.join(cache_dir, fname.replace(".png", ".pt"))
        if os.path.exists(out_path):
            continue   # already cached - skip
        img_np = np.array(Image.open(os.path.join(img_dir, fname)).convert("RGB"))
        with torch.no_grad():
            predictor.set_image(img_np)
        # Save embedding + high-res features + original image dimensions.
        # orig_hw saved as a 1-D tensor [H, W] so DataLoader stacks it
        # into shape [1, 2] and [0].tolist() correctly returns [H, W].
        # (A plain Python tuple would be collated element-by-element,
        #  causing an "unpack" error at training time.)
        torch.save({
            "image_embed"   : predictor._features["image_embed"],
            "high_res_feats": predictor._features["high_res_feats"],
            "orig_hw"       : torch.tensor(list(predictor._orig_hw[0])),  # shape [2]
        }, out_path)
        if (i + 1) % 10 == 0 or (i + 1) == len(files):
            print(f"  {i+1}/{len(files)} cached")

    print("Embedding cache complete.")


# -----------------------------------------------------------------------------
# Dataset: loads pre-cached embeddings + masks (no image encoding at train time)
# -----------------------------------------------------------------------------
class FireDatasetCached(Dataset):
    def __init__(self, cache_dir, mask_dir):
        cache_files = set(f.replace(".pt", ".png") for f in os.listdir(cache_dir) if f.endswith(".pt"))
        mask_files  = set(os.listdir(mask_dir))
        self.files     = sorted(cache_files & mask_files)
        self.cache_dir = cache_dir
        self.mask_dir  = mask_dir

    def __len__(self):
        return len(self.files)

    def __getitem__(self, idx):
        fname      = self.files[idx]
        cache_path = os.path.join(self.cache_dir, fname.replace(".png", ".pt"))
        mask_path  = os.path.join(self.mask_dir,  fname)
        cached = torch.load(cache_path, map_location="cpu", weights_only=False)
        # writePNG stored 0->0 and 1->255, threshold at 127 recovers binary mask
        mask   = np.array(Image.open(mask_path).convert("L")) > 127
        return cached, mask, fname


# =============================================================================
# STEP 5: Training loop (uses cached embeddings - encoder never called again)
# =============================================================================
def train_sam2(model, device, cache_dir, mask_dir,
               num_epochs=20, lr=1e-5, save_path=None):

    dataset = FireDatasetCached(cache_dir, mask_dir)
    # batch_size=1: embeddings can differ in spatial dims across images
    loader  = DataLoader(dataset, batch_size=1, shuffle=True)
    print(f"Training on {len(dataset)} cached samples for {num_epochs} epochs")

    # Freeze everything; unfreeze only the two small trainable heads
    for p in model.parameters():
        p.requires_grad = False
    for p in model.sam_mask_decoder.parameters():
        p.requires_grad = True
    for p in model.sam_prompt_encoder.parameters():
        p.requires_grad = True

    optimizer = torch.optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()), lr=lr
    )
    # Halve LR every 5 epochs to stabilise late training
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.5)

    for epoch in range(num_epochs):
        total_loss, n_ok = 0.0, 0

        for cached_batch, masks, fnames in loader:
            mask_np = masks[0].numpy()   # HxW bool

            box = get_box_from_mask(mask_np)
            if box is None:
                continue

            # Load pre-computed encoder outputs directly - no set_image() needed.
            # DataLoader adds a leading batch dim when stacking, so index [0] to remove it.
            # image_embed    : [1, 1, C, 64, 64]  -> [0] -> [1, C, 64, 64]
            # high_res_feats : list of [1, 1, C, H, W] -> [0] on each -> [1, C, H, W]
            # orig_hw        : [1, 2]              -> [0].tolist() -> [H, W]
            image_embed    = cached_batch["image_embed"][0].to(device)
            high_res_feats = [lvl[0].to(device)
                              for lvl in cached_batch["high_res_feats"]]
            orig_hw        = tuple(cached_batch["orig_hw"][0].tolist())   # (H, W)

            # Convert pixel box to SAM2 1024-space [1, 2, 2]
            box_t = torch.tensor(box, dtype=torch.float32, device=device)
            # Manual normalisation matching SAM2Transforms.transform_boxes:
            #   x / W * 1024,  y / H * 1024
            h, w  = orig_hw
            scale = torch.tensor([[w, h, w, h]], dtype=torch.float32, device=device)
            unnorm_box = (box_t / scale * 1024).reshape(1, 2, 2)
            box_labels  = torch.tensor([[2, 3]], dtype=torch.int, device=device)

            model.sam_mask_decoder.train()
            model.sam_prompt_encoder.train()

            with torch.enable_grad():
                sparse_emb, dense_emb = model.sam_prompt_encoder(
                    points=(unnorm_box, box_labels),
                    boxes=None,
                    masks=None,
                )

                low_res_masks, iou_pred, _, _ = model.sam_mask_decoder(
                    image_embeddings=image_embed,
                    image_pe=model.sam_prompt_encoder.get_dense_pe(),
                    sparse_prompt_embeddings=sparse_emb,
                    dense_prompt_embeddings=dense_emb,
                    multimask_output=False,
                    repeat_image=False,
                    high_res_features=high_res_feats,  # already [1, C, H, W] - no unsqueeze needed
                )

                # Upsample logits to original image size
                pred = F.interpolate(
                    low_res_masks, size=orig_hw,
                    mode="bilinear", align_corners=False
                ).squeeze()   # [H, W]

                gt   = torch.tensor(mask_np.astype(np.float32), device=device)

                # BCE + Dice loss
                bce  = F.binary_cross_entropy_with_logits(pred, gt)
                prob = torch.sigmoid(pred)
                dice = 1.0 - (2.0 * (prob * gt).sum()) / (prob.sum() + gt.sum() + 1e-6)
                loss = bce + dice

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            total_loss += loss.item()
            n_ok += 1

        scheduler.step()
        avg = total_loss / max(n_ok, 1)
        print(f"Epoch {epoch+1}/{num_epochs}  "
              f"loss={avg:.4f}  samples={n_ok}  "
              f"lr={scheduler.get_last_lr()[0]:.2e}")

    if save_path:
        # Wrap in {"model": ...} to match the format build_sam2/_load_checkpoint expects
        torch.save({"model": model.state_dict()}, save_path)
        print(f"Saved fine-tuned weights -> {save_path}")

    return model
')

# =============================================================================
# STEP 5: Load model, cache embeddings once, then train
# =============================================================================
out_cache <- "C:/app_dev/train_sam/embeddings"   # cached .pt files live here

message("Loading SAM2.1 model...")
reticulate::py_run_string(sprintf(
  'model, device = load_model(r"%s", "%s", "%s")',
  chkp_dir, yaml_name, pt_name
))

# Run every image through the encoder ONCE and save to disk.
# On subsequent runs this step is skipped for already-cached files.
# NOTE: if you see an "unpack" error, delete the embeddings folder and re-run
#       this step to regenerate the cache in the corrected format.
message("Caching image embeddings (slow first run, instant after)...")
reticulate::py_run_string(sprintf(
  'cache_embeddings(model, r"%s", r"%s")',
  out_img, out_cache
))

# Training only runs the decoder (~0.3 sec/image vs ~4 sec with encoder)
message(sprintf("Training on cached embeddings for %d samples...", length(poly_list)))
reticulate::py_run_string(sprintf(
  'train_sam2(model, device, r"%s", r"%s", num_epochs=20, lr=1e-5, save_path=r"%s")',
  out_cache, out_mask, save_path
))

message("Training complete. Fine-tuned weights saved to: ", save_path)

# =============================================================================
# STEP 6: Verify the fine-tuned weights load correctly
# =============================================================================
# Loads the saved .pt back into a fresh model instance and runs a test image
# through the full pipeline to confirm nothing is broken before using in app.R

reticulate::py_run_string(sprintf('
import os, torch, numpy as np
from PIL import Image
from sam2.build_sam import build_sam2
from sam2.sam2_image_predictor import SAM2ImagePredictor

print("--- Loading fine-tuned model ---")

# Build a fresh model with the same architecture (no pretrained weights yet)
ft_model = build_sam2(
    r"%s", os.path.join(r"%s", r"%s"),
    device=device, mode="eval", apply_postprocessing=False
)

# Load fine-tuned weights - unwrap the {"model": ...} wrapper that build_sam2 expects
state = torch.load(r"%s", map_location=device, weights_only=True)["model"]
ft_model.load_state_dict(state)
ft_model.eval()
print("Weights loaded OK")

# --- Quick inference test on the first training image ---
test_imgs = sorted(f for f in os.listdir(r"%s") if f.endswith(".png"))
if test_imgs:
    test_img_path  = os.path.join(r"%s", test_imgs[0])
    test_mask_path = os.path.join(r"%s", test_imgs[0])

    img_np  = np.array(Image.open(test_img_path).convert("RGB"))
    mask_np = np.array(Image.open(test_mask_path).convert("L")) > 127

    # Get GT bounding box as prompt
    rows = np.any(mask_np, axis=1); cols = np.any(mask_np, axis=0)
    rmin, rmax = np.where(rows)[0][[0,-1]]; cmin, cmax = np.where(cols)[0][[0,-1]]
    box = np.array([cmin, rmin, cmax, rmax], dtype=float)

    # Run predictor
    predictor = SAM2ImagePredictor(ft_model)
    with torch.no_grad():
        predictor.set_image(img_np)
        masks, scores, _ = predictor.predict(box=box, multimask_output=False)

    pred_mask = masks[0] > 0                          # HxW bool (threshold at 0 for logits)
    gt_mask   = mask_np.astype(bool)
    intersect = (pred_mask & gt_mask).sum()
    union     = (pred_mask | gt_mask).sum()
    iou       = intersect / union if union > 0 else 0.0

    print(f"Test image : {test_imgs[0]}")
    print(f"Pred mask  : {pred_mask.sum()} px  |  GT mask: {gt_mask.sum()} px")
    print(f"IoU score  : {iou:.3f}  (1.0 = perfect, >0.7 = good)")
    print(f"SAM score  : {scores[0]:.3f}")
else:
    print("No test images found - skipping inference check")

print("--- Fine-tuned model check complete ---")
', yaml_name, chkp_dir, pt_name, save_path, out_img, out_img, out_mask))
