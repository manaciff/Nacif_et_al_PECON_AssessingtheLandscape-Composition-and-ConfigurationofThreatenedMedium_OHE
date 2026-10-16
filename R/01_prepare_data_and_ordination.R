# ==============================================================================
# 01_prepare_data_and_ordination.R
# Data preparation and unconstrained ordination
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   1. Assembles the analytical dataset of the 67 sampling units: canonical
#      taxon labels, the response at full precision, the metric gaps filled
#      from script 00, structural zeros recoded only where a class is absent,
#      and the polygon centroids.
#   2. Describes the landscape inside each occupied habitat envelope (OHE).
#   3. Ordinates the land cover composition by non-metric multidimensional
#      scaling (NMDS) on Bray-Curtis dissimilarities, fits the configuration
#      metrics to the ordination (envfit), and tests the grouping factors with
#      betadisper before PERMANOVA, because the design is unbalanced.
#
# Input
#   Dados/Processados/Data_Raw_FINAL.csv               consolidated table
#   Dados/Processados/Tabela_Metricas_Recomputada.csv  written by script 00
#   Dados/Processados/AOH_area_full_precision.csv      optional; written by
#                                                      arcpy_AOH_area_full_precision.py
#   08_Dados_Especies/Dados_geo_especies/Sp_data_singlepart/*.shp
#
# Output
#   Dados/Processados/Data_Raw_WithCoords.csv   the input of scripts 02 to 08
#   Outputs/Manuscrito/PartI/                   (docs/TABLE_MAP.md lists every file)
#     TableS4_PERMANOVA_Global.csv              (Table 1)
#     Figure2_Ridgeline_Composition.tiff        (Figure 2)
#     Figure3_Ridgeline_Configuration.tiff      (Figure 3)
#     Figure1_NMDS_Composition_Genus.tiff       (Figure 4)
#     Tables S3 to S8 and S11 to S14, and TableS_Sample_Sizes_PartI.csv (S27)
#
# Guards
#   Stops unless the input has 67 rows, all 67 carry coordinates, and all 67
#   are complete over the analysis pool.
#
# Figures
#   Titles belong in the captions, not inside the images.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)          # project-relative paths
  library(sf)            # spatial vector data (centroid extraction)
  library(tidyverse)     # data wrangling and ggplot2
  library(vegan)         # NMDS, PERMANOVA, betadisper, envfit
  library(ggridges)      # ridgeline density plots
  library(ggrepel)       # non-overlapping labels in the NMDS figure
  library(scales)        # axis formatting
})

GLOBAL_SEED <- 123
N_PERM      <- 9999

# The full sampling design. Every analysis must run on this many units; the
# guards below enforce it.
N_UNITS_EXPECTED <- 67

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartI")
geo_path    <- here("08_Dados_Especies", "Dados_geo_especies")
shp_path    <- file.path(geo_path, "Sp_data_singlepart")

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Okabe-Ito colors and a fixed symbol for each genus, the same as in script 02.
colors_genus <- c(
  "Puma"         = "#D55E00",
  "Alouatta"     = "#009E73",
  "Mazama"       = "#CC79A7",
  "Brachyteles"  = "#E69F00",
  "Tayassu"      = "#0072B2",
  "Myrmecophaga" = "#F0E442",
  "Leopardus"    = "#56B4E9",
  "Bradypus"     = "#999999"
)
shapes_genus <- c(Alouatta = 1, Brachyteles = 2, Bradypus = 3, Leopardus = 4,
                  Myrmecophaga = 5, Mazama = 6, Puma = 7, Tayassu = 8)


# 2. Helper functions ----------------------------------------------------------

# Export a ggplot as a TIFF (LZW compression), at the size given in mm.
save_publication_plot <- function(plot_obj, file_name,
                                  width_mm  = 190,
                                  height_mm = 130,
                                  dpi       = 500) {
  if (is.null(plot_obj)) return(invisible(NULL))
  ggsave(filename    = file_name,
         plot        = plot_obj,
         path        = output_path,
         width       = width_mm,
         height      = height_mm,
         units       = "mm",
         dpi         = dpi,
         device      = "tiff",
         compression = "lzw",
         bg          = "white")
  message(sprintf("   Figure exported: %s (%d x %d mm, %d dpi)",
                  file_name, width_mm, height_mm, dpi))
}

# canonical_taxon(): one label per taxon, whatever spelling reaches it.
#
# The project files spell some taxa in more than one way:
#   aguriba        raster ud6_aguriba1.tif, missing the second "a" of guariba
#   baracnoides    missing the "h" of arachnoides, in an intermediate table
#   concolor       Puma rasters ud10 to ud12 drop the "p" prefix of ud1 to ud9
#   B_hypoxantus   invalid spelling of Brachyteles hypoxanthus (Kuhl, 1820)
#
# Every non-letter is removed and the name is lowercased before matching, so
# the separator never matters. An unknown name stops the run, so that no taxon
# can drop out of a join silently. Scripts 00, 02 and 03 carry the same
# function; keep the four copies identical.
canonical_taxon <- function(x) {
  key <- gsub("[^a-z]", "", tolower(as.character(x)))
  out <- dplyr::case_when(
    key %in% c("aguariba", "aguriba", "alouattaguariba")           ~ "A_guariba",
    key %in% c("barach", "baracnoides", "barachnoides",
               "brachytelesarachnoides")                           ~ "B_arachnoides",
    key %in% c("bcrinitus", "btorquatus", "bradypuscrinitus")       ~ "B_crinitus",
    key %in% c("bhypox", "bhypoxantus", "bhypoxanthus",
               "brachyteleshypoxanthus")                           ~ "B_hypoxanthus",
    key %in% c("lwiedii", "leoparduswiedii")                        ~ "L_wiedii",
    key %in% c("mazama", "mazamasp")                                ~ "Mazama",
    key %in% c("mtrida", "mtridactyla", "myrmecophagatridactyla")   ~ "M_tridactyla",
    key %in% c("pconcolor", "concolor", "pumaconcolor")             ~ "P_concolor",
    key %in% c("tpecari", "tayassupecari")                          ~ "T_pecari",
    TRUE ~ NA_character_
  )
  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognized taxon name(s): ",
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "),
         "\nAdd them to canonical_taxon() in scripts 00 to 03 before rerunning.")
  }
  out
}

# One row of a global PERMANOVA table.
tidy_adonis <- function(adonis_obj, group_label) {
  row <- as.data.frame(adonis_obj)[1L, , drop = FALSE]
  data.frame(
    Grouping = group_label,
    Df       = row$Df,
    SumOfSqs = round(row$SumOfSqs, 4),
    R2       = round(row$R2, 4),
    F_value  = round(row$F, 4),
    p_value  = row[["Pr(>F)"]]
  )
}

# Pairwise PERMANOVA between every pair of groups, with p-values adjusted for
# multiple comparisons (false discovery rate by default).
pairwise.adonis <- function(x, factors,
                            sim.method = "bray",
                            p.adjust.m = "fdr",
                            n_perm     = N_PERM,
                            seed       = GLOBAL_SEED) {
  co <- combn(unique(as.character(factors)), 2)
  pairs <- character(0); F.Model <- numeric(0)
  R2 <- numeric(0); p.value <- numeric(0)

  for (elem in seq_len(ncol(co))) {
    sub_idx <- factors %in% c(co[1, elem], co[2, elem])
    x1 <- vegdist(x[sub_idx, , drop = FALSE], method = sim.method)
    set.seed(seed)
    ad <- adonis2(x1 ~ factors[sub_idx], permutations = n_perm)
    pairs   <- c(pairs, paste(co[1, elem], "vs", co[2, elem]))
    F.Model <- c(F.Model, ad$F[1])
    R2      <- c(R2, ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }

  data.frame(
    pairs      = pairs,
    F.Model    = round(F.Model, 4),
    R2         = round(R2, 4),
    p.value    = round(p.value, 4),
    p.adjusted = round(p.adjust(p.value, method = p.adjust.m), 4),
    sig        = ifelse(p.adjust(p.value, method = p.adjust.m) < 0.05, "*", "")
  )
}

# Ridgelines with the density bounded to the possible range of each metric.
# Without bounds, the kernel tails cross into impossible values, such as
# negative densities or a percentage of landscape above 100%. ggridges
# evaluates the density only between `from` and `to`, so the curve is truncated
# at the bounds. One layer is drawn per metric, because each metric has its own
# bounds; NA means no bound on that side. The bandwidth of each metric is the
# one ggridges picks by default for a panel: the mean of the bw.nrd0()
# bandwidths of the genera with more than one value.
bounded_ridges <- function(df, bounds, ...) {
  lapply(names(bounds), function(m) {
    b  <- bounds[[m]]
    dm <- df[df$Metric == m, , drop = FALSE]
    xs <- split(dm$Value, droplevels(factor(dm$GENUS)))
    bw <- mean(vapply(xs[lengths(xs) > 1], stats::bw.nrd0, numeric(1)))
    geom_density_ridges(data = dm, bandwidth = bw,
                        from = if (is.na(b[1])) NULL else b[1],
                        to   = if (is.na(b[2])) NULL else b[2],
                        ...)
  })
}

# Because each metric is a separate layer, ggridges also visits the panels of the
# other metrics, where the layer has no data, and warns that max() or min()
# received no values. Those warnings are harmless and are muffled here; any
# other warning still reaches the console.
save_ridges <- function(...) {
  withCallingHandlers(
    save_publication_plot(...),
    warning = function(w) {
      if (grepl("no non-missing arguments to (max|min)", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
}


# 3. Data [1/10] ---------------------------------------------------------------
# Data_Raw_FINAL.csv is the consolidated table of the 67 sampling units, one
# row per unit, keyed by (SPECIES, UA_ID). UA_ID numbers the units within each
# taxon.
message("\n[1/10] Loading and formatting data...")

data_raw <- read_csv(file.path(data_path, "Data_Raw_FINAL.csv"),
                     show_col_types = FALSE)

# Data_Raw_FINAL.csv names the area column UA_Area_ha; the chain calls it
# AOH_unit_area_ha.
if ("UA_Area_ha" %in% names(data_raw) &&
    !("AOH_unit_area_ha" %in% names(data_raw))) {
  data_raw <- data_raw %>% rename(AOH_unit_area_ha = UA_Area_ha)
  message("   Column UA_Area_ha renamed to AOH_unit_area_ha.")
}

# Canonical taxon labels, applied before anything is joined or counted.
data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
data_raw$UA_ID   <- as.integer(data_raw$UA_ID)

if (nrow(data_raw) != N_UNITS_EXPECTED) {
  stop(sprintf("Data_Raw_FINAL.csv holds %d rows, expected %d.",
               nrow(data_raw), N_UNITS_EXPECTED))
}


# 3.1 The response at full precision, when available ---------------------------
# The area column of Data_Raw_FINAL.csv is stored as whole numbers, rounded
# somewhere between ArcGIS and the table. arcpy_AOH_area_full_precision.py, in
# the root of the repository, exports the area of every sampling unit from the
# same polygons at full precision. If the file it writes is present, it is used
# here; if not, the run continues on the stored values and says so. Nothing is
# imputed either way. The precision matters for the identity
# log(OHE) = log(NP) + log(100) - log(PD) checked in scripts 04 and 05.
precision_file <- file.path(data_path, "AOH_area_full_precision.csv")
if (file.exists(precision_file)) {
  prec <- read_csv(precision_file, show_col_types = FALSE) %>%
    mutate(SPECIES = canonical_taxon(SPECIES),
           UA_ID   = as.integer(UA_ID)) %>%
    dplyr::select(SPECIES, UA_ID, area_planar_ha)

  if (anyDuplicated(prec[, c("SPECIES", "UA_ID")]) > 0)
    stop("AOH_area_full_precision.csv carries duplicated (SPECIES, UA_ID) keys.")

  joined <- data_raw %>%
    dplyr::select(SPECIES, UA_ID, stored = AOH_unit_area_ha) %>%
    left_join(prec, by = c("SPECIES", "UA_ID"))

  missing_units <- sum(is.na(joined$area_planar_ha))
  if (missing_units > 0) {
    print(as.data.frame(joined[is.na(joined$area_planar_ha),
                               c("SPECIES", "UA_ID")]), row.names = FALSE)
    stop(sprintf(paste0("%d unit(s) have no area in AOH_area_full_precision.csv. ",
                        "Rerun arcpy_AOH_area_full_precision.py over all nine ",
                        "shapefiles, or delete the file to run on the stored ",
                        "values."), missing_units))
  }

  rel <- abs(joined$area_planar_ha - joined$stored) / joined$stored
  # A departure of more than five percent is not rounding; it means the polygons
  # exported are not the ones the stored table describes.
  if (max(rel) > 0.05) {
    worst <- joined[which.max(rel), ]
    stop(sprintf(paste0("Full-precision area departs from the stored value by ",
                        "%.1f%% on %s UA %d (%.1f against %.1f ha). That is too ",
                        "large to be rounding. Check that the shapefiles are the ",
                        "ones the analysis used before replacing the response."),
                 100 * max(rel), worst$SPECIES, worst$UA_ID,
                 worst$area_planar_ha, worst$stored))
  }

  data_raw <- data_raw %>%
    left_join(prec, by = c("SPECIES", "UA_ID")) %>%
    mutate(AOH_unit_area_ha = area_planar_ha) %>%
    dplyr::select(-area_planar_ha)

  message(sprintf(paste0("   Response replaced by the full-precision area: ",
                         "median change %.3f%%, maximum %.3f%%."),
                  100 * median(rel), 100 * max(rel)))
} else {
  message(paste0("   AOH_area_full_precision.csv not found; using the stored ",
                 "response. See arcpy_AOH_area_full_precision.py in the root of ",
                 "the repository."))
}


# 3.2 Metric repair and conditional structural zeros [2/10] --------------------
# Two sampling units, A_guariba UA 6 and B_arachnoides UA 1, entered
# Data_Raw_FINAL.csv with empty cells in ED and LPI of every class, and in NP, PD
# and PROX_MN of the farming class. The cells are filled from
# Tabela_Metricas_Recomputada.csv, written by script 00. Nothing is imputed:
# each value is measured again from the same per-unit raster, with the same
# reclassification and the same 8-neighbor rule used for the other 65 units.
#
# FRAGSTATS omits the row of a class that is absent from a unit, and that
# omission became NA in the consolidated table. A unit with no herbaceous patch
# has a class area of zero, a percentage of landscape of zero, no patch and no
# edge: its values are structural zeros, not unknown. A unit with forest cover
# and an empty ED_Forest cell, however, has a missing value, not a zero. A cell
# is therefore set to zero only when the class area of its class is zero. If an
# NA remains in a class that is present, the run stops and names the cells.
message("\n[2/10] Repairing metric gaps and recoding structural zeros...")

recomp_file <- file.path(data_path, "Tabela_Metricas_Recomputada.csv")
if (!file.exists(recomp_file)) {
  stop("Tabela_Metricas_Recomputada.csv not found in: ", data_path, "\n",
       "Run script 00 first. Without it, two sampling units carry empty ED, ",
       "LPI and farming-class cells, and the analyses would drop them.")
}

recomp <- read_csv(recomp_file, show_col_types = FALSE) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

repair_cols <- intersect(
  setdiff(names(recomp), c("SPECIES", "UA_ID", "UA_Area_ha_recomp")),
  names(data_raw)
)

n_before <- nrow(data_raw)
data_raw <- data_raw %>%
  left_join(
    recomp %>%
      select(SPECIES, UA_ID, all_of(repair_cols)) %>%
      rename_with(~paste0(.x, "__recomp"), all_of(repair_cols)) %>%
      # Explicit marker of a successful match. Testing a repair column instead
      # would be unreliable, because a column such as LPI_Water is legitimately
      # NA wherever the class is absent.
      mutate(matched__recomp = TRUE),
    by = c("SPECIES", "UA_ID")
  )
stopifnot(nrow(data_raw) == n_before)

n_unmatched <- sum(is.na(data_raw$matched__recomp))
if (n_unmatched > 0) {
  print(as.data.frame(
    data_raw[is.na(data_raw$matched__recomp), c("SPECIES", "UA_ID")]),
    row.names = FALSE)
  stop(sprintf(paste0("%d unit(s) of Data_Raw_FINAL.csv have no counterpart in ",
                      "the recomputed table. The join key is (SPECIES, UA_ID) ",
                      "after canonical_taxon(); check script 00."), n_unmatched))
}
data_raw$matched__recomp <- NULL

message(sprintf("   Recomputed table joined: %d column(s) available for repair.",
                length(repair_cols)))

repair_log <- list()
n_filled   <- 0L
for (cl in repair_cols) {
  src <- paste0(cl, "__recomp")
  idx <- is.na(data_raw[[cl]]) & !is.na(data_raw[[src]])
  if (any(idx)) {
    repair_log[[length(repair_log) + 1L]] <- data.frame(
      SPECIES = data_raw$SPECIES[idx],
      UA_ID   = data_raw$UA_ID[idx],
      Column  = cl,
      Value   = data_raw[[src]][idx]
    )
    data_raw[[cl]][idx] <- data_raw[[src]][idx]
    n_filled <- n_filled + sum(idx)
  }
}
data_raw <- data_raw %>% select(-ends_with("__recomp"))

if (n_filled > 0) {
  repair_df <- bind_rows(repair_log)
  write.csv(repair_df,
            file.path(output_path, "TableS_Metric_Repair_Log.csv"),
            row.names = FALSE)
  message(sprintf("   %d cell(s) filled from the recomputed table, in %d unit(s):",
                  n_filled, n_distinct(paste(repair_df$SPECIES, repair_df$UA_ID))))
  print(as.data.frame(repair_df %>% count(SPECIES, UA_ID, name = "cells_filled")),
        row.names = FALSE)
} else {
  message("   No gap to fill: Data_Raw_FINAL.csv is already complete.")
}

# Conditional structural zeros
CLASSES_LULC <- c("Forest", "Herbaceous", "Agropecuaria", "Water", "Non_Vegetated")
ZERO_METRICS <- c("CA", "PLAND", "NP", "PD", "ED", "LPI")

n_zeroed  <- 0L
bad_cells <- list()

for (cls in CLASSES_LULC) {
  ca_col <- paste0("CA_", cls)
  if (!ca_col %in% names(data_raw)) next

  ca_val    <- data_raw[[ca_col]]
  is_absent <- is.na(ca_val) | ca_val <= 0

  for (mt in ZERO_METRICS) {
    col <- paste0(mt, "_", cls)
    if (!col %in% names(data_raw)) next

    to_zero <- is_absent & is.na(data_raw[[col]])
    if (any(to_zero)) {
      data_raw[[col]][to_zero] <- 0
      n_zeroed <- n_zeroed + sum(to_zero)
    }

    still_na <- !is_absent & is.na(data_raw[[col]])
    if (any(still_na)) {
      bad_cells[[length(bad_cells) + 1L]] <- data.frame(
        SPECIES = data_raw$SPECIES[still_na],
        UA_ID   = data_raw$UA_ID[still_na],
        Class   = cls,
        Column  = col,
        CA      = ca_val[still_na]
      )
    }
  }
}

message(sprintf("   %d structural zeros recoded (class genuinely absent).", n_zeroed))

# The stop applies to the three focal classes. Water and the non-vegetated class
# were excluded at the kernel stage and enter no analysis, so an unresolved cell
# in them is reported but does not stop the run.
FOCAL_CLASSES <- c("Forest", "Herbaceous", "Agropecuaria")

if (length(bad_cells) > 0) {
  bad_df   <- bind_rows(bad_cells)
  bad_focal <- bad_df %>% filter(Class %in% FOCAL_CLASSES)
  bad_other <- bad_df %>% filter(!Class %in% FOCAL_CLASSES)

  if (nrow(bad_other) > 0) {
    message(sprintf("   %d cell(s) unresolved in non-focal classes (not used in any analysis):",
                    nrow(bad_other)))
    print(as.data.frame(bad_other %>% count(Column, name = "n_units")),
          row.names = FALSE)
  }

  if (nrow(bad_focal) > 0) {
    print(as.data.frame(bad_focal), row.names = FALSE)
    stop("The cells above are NA in a focal class that is present in the unit. ",
         "They are missing values, not structural zeros, and must not be ",
         "recoded to zero. Run script 00.")
  }
}

# Per-patch means stay NA where they are undefined: a mean over zero patches, or
# a nearest-neighbor distance for a single patch. They are not zeros and are not
# filled; their count is reported.
undefined_cols <- grep("^(PROX_MN|Mean_ENN|Mean_FRAC|Mean_AREA)_",
                       names(data_raw), value = TRUE)
n_undef <- sum(is.na(data_raw[, undefined_cols]))
if (n_undef > 0) {
  message(sprintf("   %d per-patch mean(s) remain undefined (no patch of the class).",
                  n_undef))
  undef_report <- data_raw %>%
    select(SPECIES, UA_ID, all_of(undefined_cols)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value)) %>%
    count(Column, name = "n_units")
  print(as.data.frame(undef_report), row.names = FALSE)
}


# 3.3 Guard on the analytical sample size --------------------------------------
# Every analysis of scripts 01 to 03 draws its predictors from the pool below.
# If a unit is incomplete over that pool, complete.cases() or drop_na() would
# remove it later, and the manuscript would report a sample smaller than the
# design. The run stops here instead and names the cells at fault.
ANALYSIS_POOL <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest", "ED_Forest", "Mean_FRAC_Forest", "Mean_ENN_Forest_m",
  "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria", "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)
pool_present <- intersect(ANALYSIS_POOL, names(data_raw))
n_complete   <- sum(complete.cases(data_raw[, pool_present]))

message(sprintf("   Units complete over the %d-predictor analysis pool: %d of %d.",
                length(pool_present), n_complete, nrow(data_raw)))

if (n_complete != N_UNITS_EXPECTED) {
  offenders <- data_raw %>%
    filter(!complete.cases(data_raw[, pool_present])) %>%
    select(SPECIES, UA_ID, all_of(pool_present)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  print(as.data.frame(offenders), row.names = FALSE)
  stop(sprintf(paste0("Only %d of %d units are complete over the analysis pool. ",
                      "Scripts 02 and 03 would drop the rest. Resolve ",
                      "the cells listed above before proceeding."),
               n_complete, N_UNITS_EXPECTED))
}

# Factor conversions
data_raw <- data_raw %>%
  mutate(across(c(ORDER, GENUS, SPECIES, DIET, LOCOMOTION), as.factor),
         UA_ID = as.integer(UA_ID))

message(sprintf("   Data loaded: %d observations across %d species and %d genera.",
                nrow(data_raw),
                n_distinct(data_raw$SPECIES),
                n_distinct(data_raw$GENUS)))

# (SPECIES, UA_ID) must be unique before the spatial join.
dup_check <- data_raw %>%
  group_by(SPECIES, UA_ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Duplicated (SPECIES, UA_ID) pairs detected. Spatial join will be ambiguous.")
  print(dup_check)
} else {
  message("   Integrity OK: every (SPECIES, UA_ID) pair is unique.")
}


# 4. Centroids of the sampling units [3/10] ------------------------------------
# Reads every .shp file in shp_path and extracts the centroid of each polygon.
# All layers are reprojected to the coordinate system of the first file.
#
# Fields:
#   species        SP_ID (preferred) or SPECIES; if neither exists, the name is
#                  taken from the file name, without the projection suffix
#   sampling unit  ORIG_FID (required; the run stops if it is absent)
#
# Join key: (SPECIES, UA_ID) in the table against (SP_ID, ORIG_FID) in the
# shapefiles, which gives a one-to-one match.
message("\n[3/10] Extracting polygon centroids from per-species shapefiles...")

shp_files <- list.files(shp_path, pattern = "\\.shp$", full.names = TRUE)
if (length(shp_files) == 0) {
  stop(sprintf(
    "No .shp files found in: %s\nCheck shp_path in Section 1.", shp_path
  ))
}

ref_crs <- st_crs(st_read(shp_files[1], quiet = TRUE))

extract_centroids <- function(file) {
  layer <- st_read(file, quiet = TRUE)
  if (st_crs(layer) != ref_crs) layer <- st_transform(layer, ref_crs)

  nm <- names(layer)

  # Species identifier: SP_ID preferred, SPECIES as fallback
  sp_field   <- nm[grep("^(sp_id|species)$", nm, ignore.case = TRUE)][1]

  # Sampling-unit identifier: ORIG_FID required
  orig_field <- nm[grep("orig_?fid", nm, ignore.case = TRUE)][1]
  if (is.na(orig_field)) {
    stop(sprintf(
      "ORIG_FID field not found in: %s\nAvailable fields: %s",
      basename(file), paste(nm, collapse = ", ")
    ))
  }

  # Species vector, then canonical labels. The shapefile of Brachyteles
  # hypoxanthus carries SP_ID = "B_hypoxanthus" while the table carried
  # "B_hypoxantus"; without this step the key would not match.
  sp_id <- if (!is.na(sp_field)) {
    as.character(layer[[sp_field]])
  } else {
    rep(
      sub("_(utm|albers|wgs).*$", "",
          tools::file_path_sans_ext(basename(file)),
          ignore.case = TRUE),
      nrow(layer)
    )
  }
  sp_id <- canonical_taxon(sp_id)

  # Invalid rings would make st_centroid() return an empty geometry, so the
  # geometry is repaired first. For a concave polygon the centroid can fall
  # outside the polygon; the point on the surface is then used instead, which
  # keeps the coordinate inside the sampling unit it describes.
  geom <- st_make_valid(st_geometry(layer))
  cent <- suppressWarnings(st_centroid(geom))

  outside <- !suppressMessages(
    as.logical(sf::st_within(cent, geom, sparse = FALSE)[cbind(seq_along(cent),
                                                               seq_along(cent))])
  )
  outside[is.na(outside)] <- TRUE
  if (any(outside)) {
    cent[outside] <- suppressWarnings(st_point_on_surface(geom[outside]))
    message(sprintf("   %s: %d centroid(s) fell outside the polygon; replaced by st_point_on_surface().",
                    basename(file), sum(outside)))
  }

  xy <- st_coordinates(cent)

  data.frame(
    SPECIES = sp_id,
    UA_ID   = as.integer(layer[[orig_field]]),
    coord_x = xy[, 1],
    coord_y = xy[, 2],
    stringsAsFactors = FALSE
  )
}

ua_coords <- bind_rows(lapply(shp_files, extract_centroids)) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

message(sprintf("   Centroid table built: %d records from %d shapefiles.",
                nrow(ua_coords), length(shp_files)))

# A duplicated key would make the join multiply rows instead of matching them.
if (anyDuplicated(ua_coords[, c("SPECIES", "UA_ID")]) > 0) {
  print(ua_coords[duplicated(ua_coords[, c("SPECIES", "UA_ID")]),
                  c("SPECIES", "UA_ID")])
  stop("Duplicated (SPECIES, UA_ID) in the centroid table.")
}

# Join the centroids to the table by (SPECIES, UA_ID).
data_raw <- data_raw %>%
  mutate(SPECIES_char = as.character(SPECIES)) %>%
  left_join(ua_coords,
            by = c("SPECIES_char" = "SPECIES", "UA_ID" = "UA_ID")) %>%
  select(-SPECIES_char)

n_matched <- sum(!is.na(data_raw$coord_x))
message(sprintf("   Centroids matched: %d of %d (%.1f%%).",
                n_matched, nrow(data_raw),
                100 * n_matched / nrow(data_raw)))

# A missing centroid would remove the unit from every spatial analysis of
# script 02, so it stops the run.
if (n_matched < nrow(data_raw)) {
  unmatched <- data_raw %>%
    filter(is.na(coord_x)) %>%
    distinct(SPECIES, UA_ID)
  print(as.data.frame(unmatched), row.names = FALSE)
  stop(sprintf(paste0("%d sampling unit(s) have no centroid. Check that the ",
                      "shapefile of the taxon exists in %s and that its SP_ID ",
                      "and ORIG_FID match (SPECIES, UA_ID) after ",
                      "canonical_taxon()."),
               nrow(unmatched), shp_path))
}

stopifnot(n_matched == N_UNITS_EXPECTED)
message(sprintf("   All %d sampling units carry centroid coordinates.", n_matched))


# 5. Multivariate analyses of composition [4/10] -------------------------------
# Composition matrix: the class area of the three focal classes (forest,
# farming, and herbaceous and shrub vegetation). Water and the non-vegetated
# class are excluded, as at the kernel stage.
message("\n[4/10] Multivariate analyses (NMDS, betadisper, PERMANOVA, envfit)...")

comp_vars <- data_raw %>%
  select(CA_Forest, CA_Herbaceous, CA_Agropecuaria) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
names(comp_vars) <- c("Forest", "Shrub/Herb", "Farming")

# NMDS on Bray-Curtis distances of the square root of the class area.
LULC_sqrt <- sqrt(comp_vars)
LULC_dist <- vegdist(LULC_sqrt, method = "bray", na.rm = TRUE)

set.seed(GLOBAL_SEED)
nmds_result <- metaMDS(LULC_dist, k = 2, trymax = 200, trace = FALSE)
message(sprintf("   NMDS stress (k = 2): %.3f", nmds_result$stress))

# The stress is quoted in the Results and in the caption of Figure 4.
write.csv(
  data.frame(
    Metric = c("NMDS stress (k = 2)", "Dimensions", "N sampling units",
               "Distance", "Transformation", "Permutations"),
    Value  = c(round(nmds_result$stress, 4), 2, nrow(LULC_sqrt),
               "Bray-Curtis", "square root of Class Area", N_PERM)
  ),
  file.path(output_path, "TableS_NMDS_Summary.csv"), row.names = FALSE
)
message("   TableS_NMDS_Summary.csv exported (stress for the Results and captions).")

# envfit: forest cover and six configuration metrics of forest. LPI_Forest is
# not used, because it is redundant with PLAND_Forest.
config_metrics <- data_raw %>%
  select(Mean_FRAC_Forest, Mean_ENN_Forest_m, PROX_MN_Forest,
         Mean_AREA_Forest_ha, PLAND_Forest, PD_Forest, ED_Forest)

# envfit(na.rm = TRUE) would drop any unit with an NA without reporting it, so
# the number of complete units is checked first.
n_envfit <- sum(complete.cases(config_metrics))
message(sprintf("   envfit fitted on %d of %d units.", n_envfit, nrow(config_metrics)))
if (n_envfit != N_UNITS_EXPECTED) {
  print(as.data.frame(
    data_raw[!complete.cases(config_metrics), c("SPECIES", "UA_ID")]),
    row.names = FALSE)
  stop(sprintf("envfit would run on %d units instead of %d.",
               n_envfit, N_UNITS_EXPECTED))
}

set.seed(GLOBAL_SEED)
envfit_result <- envfit(nmds_result, config_metrics,
                        permutations = N_PERM, na.rm = TRUE)

envfit_vectors <- as.data.frame(scores(envfit_result, display = "vectors"))
envfit_vectors$Metric  <- rownames(envfit_vectors)
envfit_vectors$R2      <- envfit_result$vectors$r
envfit_vectors$p_value <- envfit_result$vectors$pvals
write.csv(envfit_vectors,
          file.path(output_path, "TableS2_Envfit_Configuration.csv"),
          row.names = FALSE)

# Significant vectors, lengthened by 1.2 for display, with short labels.
envfit_codes <- c(PLAND_Forest = "PLAND", PD_Forest = "PD", ED_Forest = "ED",
                  Mean_FRAC_Forest = "FRAC_MN", Mean_ENN_Forest_m = "ENN_MN",
                  PROX_MN_Forest = "PROX_MN", Mean_AREA_Forest_ha = "AREA_MN")
envfit_sig <- subset(envfit_vectors, p_value < 0.05)
envfit_sig$NMDS1  <- envfit_sig$NMDS1 * 1.2
envfit_sig$NMDS2  <- envfit_sig$NMDS2 * 1.2
envfit_sig$Label  <- unname(envfit_codes[envfit_sig$Metric])
message(sprintf("   envfit: %d significant predictors (p < 0.05).",
                nrow(envfit_sig)))

# betadisper: homogeneity of multivariate dispersion. Significance is assessed
# with permutest() rather than anova(), because the distances to the centroid
# are not normally distributed (Anderson & Walsh 2013). The observed F is the
# one anova() would return; only the p-value changes. A fixed seed is set before
# each call. betadisper() warns that some squared distances are negative and
# changed to zero: Bray-Curtis dissimilarities are not Euclidean, and the
# warning reports the correction that vegan applies.
disp_genus   <- betadisper(LULC_dist, data_raw$GENUS)
disp_species <- betadisper(LULC_dist, data_raw$SPECIES)
disp_diet    <- betadisper(LULC_dist, data_raw$DIET)
disp_loc     <- betadisper(LULC_dist, data_raw$LOCOMOTION)

set.seed(GLOBAL_SEED); permd_genus   <- permutest(disp_genus,   permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_species <- permutest(disp_species, permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_diet    <- permutest(disp_diet,    permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_loc     <- permutest(disp_loc,     permutations = N_PERM)

betadisper_results <- bind_rows(
  as.data.frame(permd_genus$tab)   %>% mutate(Group = "Genus",      Source = rownames(.)),
  as.data.frame(permd_species$tab) %>% mutate(Group = "Species",    Source = rownames(.)),
  as.data.frame(permd_diet$tab)    %>% mutate(Group = "Diet",       Source = rownames(.)),
  as.data.frame(permd_loc$tab)     %>% mutate(Group = "Locomotion", Source = rownames(.))
) %>% select(Group, Source, everything())

message(sprintf(
  "   betadisper (permutest, %d perm): Genus p = %.4f | Species p = %.4f | Diet p = %.4f | Locomotion p = %.4f",
  N_PERM,
  permd_genus$tab[1, "Pr(>F)"],   permd_species$tab[1, "Pr(>F)"],
  permd_diet$tab[1, "Pr(>F)"],    permd_loc$tab[1, "Pr(>F)"]))

write.csv(betadisper_results,
          file.path(output_path, "TableS3_Betadisper.csv"),
          row.names = FALSE)

# PERMANOVA
set.seed(GLOBAL_SEED)
adonis_genus   <- adonis2(LULC_dist ~ GENUS,      data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_species <- adonis2(LULC_dist ~ SPECIES,    data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_diet    <- adonis2(LULC_dist ~ DIET,       data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_loc     <- adonis2(LULC_dist ~ LOCOMOTION, data = data_raw, permutations = N_PERM)

permanova_global <- bind_rows(
  tidy_adonis(adonis_genus,   "Genus"),
  tidy_adonis(adonis_species, "Species"),
  tidy_adonis(adonis_diet,    "Diet"),
  tidy_adonis(adonis_loc,     "Locomotion")
)
write.csv(permanova_global,
          file.path(output_path, "TableS4_PERMANOVA_Global.csv"),
          row.names = FALSE)

write.csv(pairwise.adonis(LULC_sqrt, data_raw$GENUS),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Genus.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$SPECIES),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Species.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$DIET),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Diet.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$LOCOMOTION),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Locomotion.csv"),
          row.names = FALSE)
message("   Pairwise PERMANOVA tables exported.")


# 6. Descriptive tables by species [5/10] --------------------------------------
message("\n[5/10] Building descriptive tables by species...")

# Shannon diversity of the three focal classes.
comp_vars_pland <- data_raw %>%
  select(PLAND_Forest, PLAND_Agropecuaria, PLAND_Herbaceous) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
data_raw$Shannon_Index <- diversity(comp_vars_pland, index = "shannon")

data_raw$MPS_Forest_ha <- data_raw$CA_Forest / data_raw$NP_Forest
data_raw$MPS_Forest_ha[is.infinite(data_raw$MPS_Forest_ha)] <- NA

# Composition by species
table_composition <- data_raw %>%
  group_by(SPECIES) %>%
  summarise(
    N_SUs                  = n(),
    PLAND_Forest_Mean      = mean(PLAND_Forest,         na.rm = TRUE),
    PLAND_Forest_SD        = sd(PLAND_Forest,           na.rm = TRUE),
    PLAND_Forest_Median    = median(PLAND_Forest,       na.rm = TRUE),
    PLAND_Forest_Min       = min(PLAND_Forest,          na.rm = TRUE),
    PLAND_Forest_Max       = max(PLAND_Forest,          na.rm = TRUE),
    PLAND_Farming_Mean     = mean(PLAND_Agropecuaria,   na.rm = TRUE),
    PLAND_Farming_SD       = sd(PLAND_Agropecuaria,     na.rm = TRUE),
    PLAND_Farming_Median   = median(PLAND_Agropecuaria, na.rm = TRUE),
    PLAND_ShrubHerb_Mean   = mean(PLAND_Herbaceous,     na.rm = TRUE),
    PLAND_ShrubHerb_SD     = sd(PLAND_Herbaceous,       na.rm = TRUE),
    PLAND_ShrubHerb_Median = median(PLAND_Herbaceous,   na.rm = TRUE),
    Shannon_Mean           = mean(Shannon_Index,        na.rm = TRUE),
    Shannon_SD             = sd(Shannon_Index,          na.rm = TRUE),
    Shannon_Median         = median(Shannon_Index,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(PLAND_Forest_Mean))

write.csv(table_composition,
          file.path(output_path, "Table1_Composition_By_Species.csv"),
          row.names = FALSE)

# Configuration by species
table_configuration <- data_raw %>%
  group_by(SPECIES) %>%
  summarise(
    N_SUs                  = n(),
    NP_Forest_Mean         = mean(NP_Forest,             na.rm = TRUE),
    NP_Forest_Median       = median(NP_Forest,           na.rm = TRUE),
    PD_Forest_Mean         = mean(PD_Forest,             na.rm = TRUE),
    PD_Forest_SD           = sd(PD_Forest,               na.rm = TRUE),
    PD_Forest_Median       = median(PD_Forest,           na.rm = TRUE),
    AREA_MN_Forest_Mean    = mean(Mean_AREA_Forest_ha,   na.rm = TRUE),
    AREA_MN_Forest_SD      = sd(Mean_AREA_Forest_ha,     na.rm = TRUE),
    AREA_MN_Forest_Median  = median(Mean_AREA_Forest_ha, na.rm = TRUE),
    MPS_Forest_Mean        = mean(MPS_Forest_ha,         na.rm = TRUE),
    MPS_Forest_Median      = median(MPS_Forest_ha,       na.rm = TRUE),
    ED_Forest_Mean         = mean(ED_Forest,             na.rm = TRUE),
    ED_Forest_Median       = median(ED_Forest,           na.rm = TRUE),
    FRAC_Forest_Mean       = mean(Mean_FRAC_Forest,      na.rm = TRUE),
    FRAC_Forest_Median     = median(Mean_FRAC_Forest,    na.rm = TRUE),
    ENN_Forest_Mean        = mean(Mean_ENN_Forest_m,     na.rm = TRUE),
    ENN_Forest_Median      = median(Mean_ENN_Forest_m,   na.rm = TRUE),
    PROX_Forest_Mean       = mean(PROX_MN_Forest,        na.rm = TRUE),
    PROX_Forest_Median     = median(PROX_MN_Forest,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(AREA_MN_Forest_Mean))

write.csv(table_configuration,
          file.path(output_path, "Table2_Configuration_By_Species.csv"),
          row.names = FALSE)
message("   Composition and Configuration tables exported.")


# 7. Figure 4: NMDS with envfit vectors [6/10] ---------------------------------
# Polygons enclose the sampling units of each genus; arrows are the significant
# envfit vectors (forest cover and the configuration metrics of forest).
message("\n[6/10] Building Figure 4 (NMDS + envfit)...")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Genus <- data_raw$GENUS

hulls_genus <- nmds_scores %>%
  group_by(Genus) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()

genera_present      <- as.character(unique(data_raw$GENUS))
colors_genus_subset <- colors_genus[names(colors_genus) %in% genera_present]

# Vector labels start just beyond the arrow tips. The site points enter the same
# layer with empty labels, so that the labels are repelled from them as well
# (ggrepel does not draw empty labels).
nmds_labels <- rbind(
  data.frame(x = envfit_sig$NMDS1, y = envfit_sig$NMDS2, label = envfit_sig$Label),
  data.frame(x = nmds_scores$NMDS1, y = nmds_scores$NMDS2, label = "")
)
nmds_nudge_x <- ifelse(nmds_labels$label == "", 0, 0.08 * nmds_labels$x)
nmds_nudge_y <- ifelse(nmds_labels$label == "", 0, 0.08 * nmds_labels$y)

plot_nmds_genus <- ggplot() +
  geom_polygon(data  = hulls_genus,
               aes(x = NMDS1, y = NMDS2,
                   fill = Genus, group = Genus),
               alpha = 0.30, colour = "grey30", linewidth = 0.3) +
  geom_point(data  = nmds_scores,
             aes(x = NMDS1, y = NMDS2,
                 colour = Genus, shape = Genus),
             size = 2.2, stroke = 0.8, alpha = 0.9) +
  geom_segment(data = envfit_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.18, "cm")),
               colour = "grey15", linewidth = 0.7) +
  ggrepel::geom_text_repel(data = nmds_labels,
                           aes(x = x, y = y, label = label),
                           nudge_x = nmds_nudge_x, nudge_y = nmds_nudge_y,
                           colour = "black", size = 2.8, fontface = "bold",
                           segment.colour = "grey60", min.segment.length = 0.2,
                           box.padding = 0.2, max.overlaps = Inf,
                           seed = GLOBAL_SEED) +
  scale_fill_manual(values  = colors_genus_subset, na.value = "grey90") +
  scale_colour_manual(values = colors_genus_subset, na.value = "grey50") +
  scale_shape_manual(values  = shapes_genus[genera_present]) +
  coord_equal() +
  labs(x = "NMDS 1", y = "NMDS 2") +
  theme_classic(base_size = 9) +
  theme(legend.text     = element_text(face = "italic"),
        legend.position = "right")

save_publication_plot(plot_nmds_genus,
                      file_name = "Figure1_NMDS_Composition_Genus.tiff",
                      width_mm  = 190, height_mm = 150)


# 8. Figures 2 and 3: the landscape inside the envelopes [7/10] ----------------
message("\n[7/10] Building Figure 2 (composition ridgeline) and Figure 3 (configuration ridgeline)...")

# The Metric codes below are also the column names of Table3_Medians_*.csv and
# Table4_Medians_*.csv; the panels show the descriptive labels.

# Figure 2: the three focal classes
df_long_comp <- data_raw %>%
  select(GENUS,
         `Forest`     = PLAND_Forest,
         `Farming`    = PLAND_Agropecuaria,
         `Shrub/Herb` = PLAND_Herbaceous) %>%
  pivot_longer(cols      = -GENUS,
               names_to  = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Metric = factor(Metric,
                         levels = c("Forest", "Farming", "Shrub/Herb")))

median_comp <- df_long_comp %>%
  group_by(Metric, GENUS) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")

# Genera ordered by their median forest cover, the most forested at the top.
genus_order_comp <- df_long_comp %>%
  filter(Metric == "Forest") %>%
  group_by(GENUS) %>%
  summarise(med = median(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(GENUS) %>% as.character()

df_long_comp$GENUS <- factor(df_long_comp$GENUS, levels = genus_order_comp)
median_comp$GENUS  <- factor(median_comp$GENUS,  levels = genus_order_comp)

comp_strip_labels <- c("Forest"     = "Forest",
                       "Farming"    = "Farming",
                       "Shrub/Herb" = "Herbaceous and shrub vegetation")

# Percentages are bounded by 0 and 100; the herbaceous panel keeps its own
# upper limit, so that its narrow range stays readable.
comp_bounds <- list("Forest"     = c(0, 100),
                    "Farming"    = c(0, 100),
                    "Shrub/Herb" = c(0, NA))

plot_ridge_comp <- ggplot(df_long_comp,
                          aes(x = Value, y = GENUS, fill = GENUS)) +
  bounded_ridges(df_long_comp, comp_bounds,
                 scale            = 1.15,
                 rel_min_height   = 0.01,
                 alpha            = 0.75,
                 colour           = "grey25",
                 linewidth        = 0.3,
                 quantile_lines   = TRUE,
                 quantiles        = 0.5,
                 vline_colour     = "grey10",
                 vline_width      = 0.7) +
  scale_fill_manual(values = colors_genus[genus_order_comp]) +
  scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1,
             labeller = as_labeller(comp_strip_labels)) +
  labs(x = "Percentage of landscape (PLAND)", y = NULL) +
  theme_ridges(font_size = 9, grid = TRUE, center_axis_labels = TRUE) +
  theme(legend.position    = "none",
        strip.text         = element_text(face = "bold", size = 8,
                                          colour = "white",
                                          margin = margin(3, 3, 3, 3)),
        strip.background   = element_rect(fill = "grey25", colour = NA),
        axis.text.y        = element_text(face = "italic", size = 8),
        axis.text.x        = element_text(size = 7),
        axis.title.x       = element_text(face = "bold", size = 9,
                                          margin = margin(t = 6)),
        panel.spacing      = unit(1, "lines"),
        plot.margin        = margin(6, 8, 6, 6))

save_ridges(plot_ridge_comp,
            file_name = "Figure2_Ridgeline_Composition.tiff",
            width_mm  = 190, height_mm = 95)

# Figure 3: configuration metrics of forest
df_long_conf <- data_raw %>%
  transmute(
    GENUS,
    `PD (n / 100 ha)` = PD_Forest,
    `AREA_MN (ha)`    = Mean_AREA_Forest_ha,
    `ED (m / ha)`     = ED_Forest,
    `FRAC_MN`         = Mean_FRAC_Forest,
    `ENN_MN (m)`      = Mean_ENN_Forest_m,
    `PROX_MN`         = PROX_MN_Forest
  ) %>%
  pivot_longer(cols      = -GENUS,
               names_to  = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value), is.finite(Value)) %>%
  mutate(Metric = factor(Metric,
                         levels = c("PD (n / 100 ha)", "AREA_MN (ha)",
                                    "ED (m / ha)",     "FRAC_MN",
                                    "ENN_MN (m)",      "PROX_MN")))

df_long_conf$GENUS <- factor(df_long_conf$GENUS, levels = genus_order_comp)

conf_strip_labels <- c(
  "PD (n / 100 ha)" = "Patch density\n(PD, patches/100 ha)",
  "AREA_MN (ha)"    = "Mean patch area\n(AREA_MN, ha)",
  "ED (m / ha)"     = "Edge density\n(ED, m/ha)",
  "FRAC_MN"         = "Shape complexity\n(FRAC_MN)",
  "ENN_MN (m)"      = "Nearest-neighbor distance\n(ENN_MN, m)",
  "PROX_MN"         = "Proximity index\n(PROX_MN)"
)

# Every metric is non-negative; FRAC_MN lies between 1 and 2 and needs no bound.
conf_bounds <- list("PD (n / 100 ha)" = c(0, NA),
                    "AREA_MN (ha)"    = c(0, NA),
                    "ED (m / ha)"     = c(0, NA),
                    "FRAC_MN"         = c(NA, NA),
                    "ENN_MN (m)"      = c(0, NA),
                    "PROX_MN"         = c(0, NA))

plot_ridge_conf <- ggplot(df_long_conf,
                          aes(x = Value, y = GENUS, fill = GENUS)) +
  bounded_ridges(df_long_conf, conf_bounds,
                 scale            = 1.10,
                 rel_min_height   = 0.01,
                 alpha            = 0.75,
                 colour           = "grey25",
                 linewidth        = 0.3,
                 quantile_lines   = TRUE,
                 quantiles        = 0.5,
                 vline_colour     = "grey10",
                 vline_width      = 0.7) +
  scale_fill_manual(values = colors_genus[genus_order_comp]) +
  scale_x_continuous(labels = label_number(big.mark = ","),
                     expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 3,
             labeller = as_labeller(conf_strip_labels)) +
  labs(x = "Metric value",
       y = NULL) +
  theme_ridges(font_size = 9, grid = TRUE, center_axis_labels = TRUE) +
  theme(legend.position    = "none",
        strip.text         = element_text(face = "bold", size = 8,
                                          colour = "white", lineheight = 0.9,
                                          margin = margin(3, 3, 3, 3)),
        strip.background   = element_rect(fill = "grey25", colour = NA),
        axis.text.y        = element_text(face = "italic", size = 8),
        axis.text.x        = element_text(size = 7),
        axis.title.x       = element_text(face = "bold", size = 9,
                                          margin = margin(t = 6)),
        panel.spacing.x    = unit(1.0, "lines"),
        panel.spacing.y    = unit(0.8, "lines"),
        plot.margin        = margin(6, 8, 6, 6))

save_ridges(plot_ridge_conf,
            file_name = "Figure3_Ridgeline_Configuration.tiff",
            width_mm  = 190, height_mm = 200)


# 9. Median tables [8/10] ------------------------------------------------------
message("\n[8/10] Exporting median summary tables (for captions and Results)...")

median_comp_wide <- median_comp %>%
  pivot_wider(names_from = Metric, values_from = Median)
write.csv(median_comp_wide,
          file.path(output_path, "Table3_Medians_Composition_PLAND.csv"),
          row.names = FALSE)

median_conf_wide <- df_long_conf %>%
  group_by(Metric, GENUS) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metric, values_from = Median)
write.csv(median_conf_wide,
          file.path(output_path, "Table4_Medians_Configuration.csv"),
          row.names = FALSE)

# Statistics of the minor classes and of the largest patch index, quoted in
# Supplementary Tables S3, S4 and S6.
write.csv(
  data_raw %>%
    group_by(GENUS) %>%
    summarise(PLAND_NonVegetated_Median = median(PLAND_Non_Vegetated, na.rm = TRUE),
              PLAND_Water_Median        = median(PLAND_Water, na.rm = TRUE),
              .groups = "drop"),
  file.path(output_path, "Table3b_Medians_MinorClasses.csv"), row.names = FALSE)

write.csv(
  data_raw %>%
    group_by(SPECIES) %>%
    summarise(
      N_SUs                     = n(),
      PLAND_Herbaceous_Mean     = mean(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_Herbaceous_SD       = sd(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_Herbaceous_Median   = median(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_NonVegetated_Mean   = mean(PLAND_Non_Vegetated, na.rm = TRUE),
      PLAND_NonVegetated_SD     = sd(PLAND_Non_Vegetated, na.rm = TRUE),
      PLAND_NonVegetated_Median = median(PLAND_Non_Vegetated, na.rm = TRUE),
      LPI_Forest_Mean           = mean(LPI_Forest, na.rm = TRUE),
      LPI_Forest_Median         = median(LPI_Forest, na.rm = TRUE),
      PLAND_Forest_Min          = min(PLAND_Forest, na.rm = TRUE),
      PLAND_Forest_Max          = max(PLAND_Forest, na.rm = TRUE),
      .groups = "drop"),
  file.path(output_path, "Table1b_MinorClasses_And_LPI_By_Species.csv"), row.names = FALSE)

message("   Table3b and Table1b exported.")


# 9.1 Sample size used by each analysis (Table S27) ----------------------------
sample_sizes_I <- data.frame(
  Analysis = c("Rows loaded",
               "NMDS on Bray-Curtis of sqrt(CA)",
               "PERMANOVA and betadisper",
               "envfit of configuration metrics",
               "Units with centroid coordinates",
               "Complete over the analysis pool"),
  N        = c(nrow(data_raw),
               nrow(LULC_sqrt),
               nrow(as.matrix(LULC_dist)),
               n_envfit,
               sum(!is.na(data_raw$coord_x)),
               sum(complete.cases(data_raw[, pool_present])))
)
write.csv(sample_sizes_I,
          file.path(output_path, "TableS_Sample_Sizes_PartI.csv"),
          row.names = FALSE)
message("\n   Sample size per analysis:")
print(sample_sizes_I, row.names = FALSE)

if (any(sample_sizes_I$N != N_UNITS_EXPECTED)) {
  stop("At least one analysis of script 01 did not use all ", N_UNITS_EXPECTED,
       " sampling units. See the table above.")
}


# 10. Export the data for scripts 02 to 08 [9/10] ------------------------------
message("\n[9/10] Exporting Data_Raw_WithCoords.csv...")

if (!all(c("coord_x", "coord_y") %in% names(data_raw))) {
  stop(
    "coord_x / coord_y columns are missing from data_raw.\n",
    "Check Section 4 (centroid extraction) for errors before proceeding."
  )
}

pct_with_coords <- 100 * mean(!is.na(data_raw$coord_x))
if (pct_with_coords < 100) {
  print(as.data.frame(
    data_raw[is.na(data_raw$coord_x), c("SPECIES", "UA_ID")]), row.names = FALSE)
  stop(sprintf(paste0("Only %.1f%% of rows have coordinates. Script 02 would ",
                      "drop the rest from the MEM screening, the variation ",
                      "partitioning and the RDA. Resolve them in Section 4."),
               pct_with_coords))
}

# Final guard: scripts 02 to 08 all read this file.
stopifnot(
  nrow(data_raw) == N_UNITS_EXPECTED,
  sum(complete.cases(data_raw[, pool_present])) == N_UNITS_EXPECTED,
  all(!is.na(data_raw$coord_x)), all(!is.na(data_raw$coord_y))
)

write.csv(data_raw,
          file.path(data_path, "Data_Raw_WithCoords.csv"),
          row.names = FALSE)
message(sprintf(
  "   Data_Raw_WithCoords.csv exported (%d rows, %.1f%% with coordinates).",
  nrow(data_raw), pct_with_coords
))


# 11. Output inventory [10/10] -------------------------------------------------
message("\n[10/10] All analyses of script 01 completed.")
message("\nOutputs in: ", output_path)
message("  Figures:")
message("    Figure1_NMDS_Composition_Genus.tiff      (Figure 4)")
message("    Figure2_Ridgeline_Composition.tiff       (Figure 2)")
message("    Figure3_Ridgeline_Configuration.tiff     (Figure 3)")
message("  Tables:")
message("    Table1_Composition_By_Species.csv")
message("    Table2_Configuration_By_Species.csv")
message("    Table3_Medians_Composition_PLAND.csv")
message("    Table4_Medians_Configuration.csv")
message("    TableS2_Envfit_Configuration.csv")
message("    TableS3_Betadisper.csv")
message("    TableS4_PERMANOVA_Global.csv")
message("    TableS5_PERMANOVA_Pairwise_{Genus,Species,Diet,Locomotion}.csv")
message("    TableS_Metric_Repair_Log.csv            (cells filled from script 00)")
message("    TableS_Sample_Sizes_PartI.csv           (n per analysis; all = 67)")
message("\nData written for scripts 02 to 08:")
message("    ", file.path(data_path, "Data_Raw_WithCoords.csv"))


# 12. Session information ------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part01.txt"))
