# ==============================================================================
# 00_recompute_class_metrics.R
# Recompute the class-level landscape metrics from the per-unit rasters
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# What this script does
#   Recomputes every class-level metric of the 67 sampling units directly from
#   the per-unit land cover rasters and compares the result with the values
#   stored in the consolidated dataset. Two units entered that dataset with
#   empty cells; script 01 fills them from the table written here. Nothing is
#   imputed: every filled value is measured again from the same raster.
#
# Input
#   Dados/Rasters/UA_RASTER/*.tif     67 per-unit rasters (not in the repository;
#                                     see README.md)
#   Dados/FRAGSTATS_RESULT/*Patch.CSV FRAGSTATS patch tables (proximity index)
#   Dados/Processados/Data_Raw_FINAL.csv
#
# Output
#   Dados/Processados/Tabela_Metricas_Recomputada.csv
#   Outputs/Manuscrito/Part0/TableS_Recompute_Agreement.csv   (Table S32)
#   Outputs/Manuscrito/Part0/TableS_Repaired_Cells.csv        (Table S33)
#
# Run time
#   A few minutes. More than 99% of each raster is empty padding, so every
#   raster is cropped with terra::trim() before the metrics are computed. The
#   result of each unit is cached, so an interrupted run resumes where it
#   stopped.
# ==============================================================================


# 1. Setup ---------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(landscapemetrics)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tidyr)
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

# RECOMPUTE_PATCH_MEANS: FALSE computes only the metrics needed to repair the
# dataset (CA, PLAND, NP, PD, ED and LPI). TRUE also computes AREA_MN, FRAC_MN
# and ENN_MN, which never had empty cells and take most of the run time.
RECOMPUTE_PATCH_MEANS <- FALSE

# USE_CACHE: TRUE stores the result of each unit in cache_path and skips units
# already computed. Delete that folder to force a full recomputation.
USE_CACHE <- TRUE

raster_path    <- here("Dados", "Rasters", "UA_RASTER")
fragstats_path <- here("Dados", "FRAGSTATS_RESULT")
out_path       <- here("Dados", "Processados")
report_path    <- here("Outputs", "Manuscrito", "Part0")
cache_path     <- here("Dados", "Processados", "_part0_cache")

if (!dir.exists(report_path)) dir.create(report_path, recursive = TRUE)
if (USE_CACHE && !dir.exists(cache_path)) dir.create(cache_path, recursive = TRUE)

N_UNITS_EXPECTED <- 67

# Class codes of the MapBiomas aggregation used in the study. Farming is the
# union of pasture (code 3) and agriculture (code 4), so code 4 is reclassified
# into code 3 BEFORE the metrics are computed. Adjacent pasture and agriculture
# pixels then form a single patch, which is what makes NP, PD and LPI correct
# for the farming class.
CLASS_CODES    <- c(Forest = 1, Herbaceous = 2, Agropecuaria = 3,
                    Water = 5, Non_Vegetated = 6)
FARMING_FROM   <- 4
FARMING_TO     <- 3


# 2. Taxon names ---------------------------------------------------------------
# The project files spell some taxa in more than one way. Every variant is
# listed here, so the correction is explicit:
#
#   aguriba       raster ud6_aguriba1.tif, missing the second "a" of guariba
#   baracnoides   missing the "h" of arachnoides, in an intermediate table
#   concolor      Puma rasters ud10 to ud12 drop the "p" prefix of ud1 to ud9
#   B_hypoxantus  invalid spelling of Brachyteles hypoxanthus (Kuhl, 1820)
#
# canonical_taxon() stops on an unknown name instead of returning NA, so a
# misspelled taxon cannot drop out of a join silently. Scripts 01, 02 and 03
# carry the same function; a new variant must be added to all four copies.

canonical_taxon <- function(x) {
  key <- x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z]", "")   # drops underscores, spaces, digits, dots

  dplyr::case_when(
    key %in% c("aguariba", "aguriba", "alouattaguariba")        ~ "A_guariba",
    key %in% c("barach", "baracnoides", "barachnoides",
               "brachytelesarachnoides")                        ~ "B_arachnoides",
    key %in% c("bcrinitus", "btorquatus", "bradypuscrinitus")    ~ "B_crinitus",
    key %in% c("bhypox", "bhypoxantus", "bhypoxanthus",
               "brachyteleshypoxanthus")                        ~ "B_hypoxanthus",
    key %in% c("lwiedii", "leoparduswiedii")                     ~ "L_wiedii",
    key %in% c("mazama", "mazamasp")                             ~ "Mazama",
    key %in% c("mtrida", "mtridactyla",
               "myrmecophagatridactyla")                         ~ "M_tridactyla",
    key %in% c("pconcolor", "concolor", "pumaconcolor")          ~ "P_concolor",
    key %in% c("tpecari", "tayassupecari")                       ~ "T_pecari",
    TRUE ~ NA_character_
  ) -> out

  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognized taxon name(s): ",
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "),
         ". Add them to canonical_taxon() in Parts 0, I, II and III before ",
         "rerunning.", call. = FALSE)
  }
  out
}


# 3. Index the rasters [1/5] ---------------------------------------------------
message("\n[1/5] Indexing the per-unit rasters...")

files <- list.files(raster_path, pattern = "\\.tif$", full.names = TRUE)
message(sprintf("   %d rasters found in %s", length(files), raster_path))

if (length(files) != N_UNITS_EXPECTED) {
  stop(sprintf(paste0("Expected %d rasters, found %d. The recomputation must ",
                      "cover every sampling unit before it is used to repair ",
                      "the dataset."),
               N_UNITS_EXPECTED, length(files)))
}

keys <- tibble(
  path    = files,
  file    = basename(files),
  # File name pattern: ud<unit>_<taxon><n>.tif
  UA_ID   = as.integer(str_match(basename(files), "^ud(\\d+)_")[, 2]),
  raw_tax = str_match(basename(files), "^ud\\d+_([a-z]+)\\d*\\.tif$")[, 2]
) %>%
  mutate(SPECIES = canonical_taxon(raw_tax))

if (any(is.na(keys$SPECIES)) || any(is.na(keys$UA_ID))) {
  print(keys %>% filter(is.na(SPECIES) | is.na(UA_ID)), n = Inf)
  stop("Unparsed file names. Add the pattern to canonical_taxon() before rerunning.")
}
if (anyDuplicated(keys[, c("SPECIES", "UA_ID")]) > 0) {
  print(keys[duplicated(keys[, c("SPECIES", "UA_ID")]), ], n = Inf)
  stop("Duplicated (SPECIES, UA_ID) after normalization.")
}

message(sprintf("   %d unique (taxon, unit) pairs across %d taxa.",
                nrow(keys), n_distinct(keys$SPECIES)))
print(keys %>% count(SPECIES, name = "n_units"), n = Inf)


# 4. Class-level metrics from the rasters [2/5] --------------------------------
# directions = 8 and the default count_boundary = FALSE of lsm_c_ed() are the
# settings of the original extraction. Changing them here would require
# recomputing the whole column for all 67 units.
message("\n[2/5] Computing class-level metrics with landscapemetrics...")

# ED, LPI, NP and PD are the columns that are empty in the two incomplete
# units, so they must be part of the core set. The three per-patch means are
# optional (see RECOMPUTE_PATCH_MEANS).
WHAT_CORE <- c("lsm_c_ca", "lsm_c_pland", "lsm_c_np", "lsm_c_pd",
               "lsm_c_ed", "lsm_c_lpi")
WHAT_SLOW <- c("lsm_c_area_mn", "lsm_c_frac_mn", "lsm_c_enn_mn")
WHAT_USE  <- if (RECOMPUTE_PATCH_MEANS) c(WHAT_CORE, WHAT_SLOW) else WHAT_CORE

message(sprintf("   Metrics requested: %s", paste(WHAT_USE, collapse = ", ")))
if (!RECOMPUTE_PATCH_MEANS) {
  message("   AREA_MN, FRAC_MN and ENN_MN skipped (RECOMPUTE_PATCH_MEANS = FALSE).")
  message("   They are not repaired by this script and ENN_MN dominates the run time.")
}

unit_metrics <- function(path, sp, ua) {
  r <- terra::rast(path)
  n_before <- terra::ncell(r)

  # Merge agriculture into pasture, so that farming is a single class.
  r <- terra::classify(r, cbind(FARMING_FROM, FARMING_TO))

  # Crop to the bounding box of the valid pixels. Each per-unit raster covers
  # the extent of every unit of the taxon, so most of the file is NoData.
  r <- terra::trim(r)

  # Add back a one-cell NoData margin, so that every valid pixel keeps the
  # neighborhood it had in the original file and edge density does not react
  # to the crop. Section 6 checks this against the stored values.
  r <- terra::extend(r, 1)

  n_after <- terra::ncell(r)

  # freq() counts the valid cells without loading the raster into memory.
  fr      <- terra::freq(r)
  area_ha <- sum(fr$count) * prod(terra::res(r)) / 10000

  res <- landscapemetrics::calculate_lsm(
    r,
    level      = "class",
    what       = WHAT_USE,
    directions = 8,
    verbose    = FALSE
  )

  attr(res, "n_before") <- n_before
  attr(res, "n_after")  <- n_after

  res %>%
    mutate(SPECIES = sp, UA_ID = ua, UA_Area_ha_recomp = area_ha) %>%
    select(SPECIES, UA_ID, UA_Area_ha_recomp, class, metric, value)
}

# The cache file name includes the metric set, so switching
# RECOMPUTE_PATCH_MEANS never reuses a result computed with the other set.
cache_file <- function(sp, ua) {
  file.path(cache_path,
            sprintf("%s_UA%02d_%s.rds", sp, ua,
                    if (RECOMPUTE_PATCH_MEANS) "full" else "core"))
}

n_units  <- nrow(keys)
t_start  <- Sys.time()
out_list <- vector("list", n_units)
n_cached <- 0L

for (i in seq_len(n_units)) {
  sp <- keys$SPECIES[i]; ua <- keys$UA_ID[i]; pth <- keys$path[i]
  cf <- cache_file(sp, ua)

  if (USE_CACHE && file.exists(cf)) {
    out_list[[i]] <- readRDS(cf)
    n_cached <- n_cached + 1L
    message(sprintf("   [%2d/%2d] %-16s UA %-3d  cached", i, n_units, sp, ua))
    next
  }

  t0  <- Sys.time()
  res <- unit_metrics(pth, sp, ua)
  dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (USE_CACHE) saveRDS(res, cf)
  out_list[[i]] <- res

  done      <- i - n_cached
  elapsed   <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  remaining <- if (done > 0) (elapsed / done) * (n_units - i) else NA_real_

  message(sprintf(
    "   [%2d/%2d] %-16s UA %-3d  %6.1f s  |  %d -> %d cells (%.0fx smaller)  |  ETA %s",
    i, n_units, sp, ua, dt,
    attr(res, "n_before"), attr(res, "n_after"),
    attr(res, "n_before") / attr(res, "n_after"),
    if (is.na(remaining)) "-" else format(round(remaining / 60, 1), nsmall = 1)
  ))
}

long_metrics <- bind_rows(out_list)

message(sprintf("   Done in %.1f min (%d unit(s) taken from cache).",
                as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                n_cached))

# Map the landscapemetrics names onto the column names of the project. The
# forest columns must match Data_Raw_FINAL.csv exactly, because scripts 01 to
# 05 refer to them by name.
metric_prefix <- c(
  ca      = "CA",
  pland   = "PLAND",
  np      = "NP",
  pd      = "PD",
  ed      = "ED",
  lpi     = "LPI",
  area_mn = "Mean_AREA",
  frac_mn = "Mean_FRAC",
  enn_mn  = "Mean_ENN"
)
metric_suffix <- c(area_mn = "_ha", enn_mn = "_m")   # unit suffixes

# Only the metrics actually computed are mapped.
metric_prefix <- metric_prefix[names(metric_prefix) %in% unique(long_metrics$metric)]

class_labels <- tibble(class = unname(CLASS_CODES), class_name = names(CLASS_CODES))

wide_raster <- long_metrics %>%
  inner_join(class_labels, by = "class") %>%
  mutate(
    column = paste0(metric_prefix[metric], "_", class_name,
                    dplyr::coalesce(metric_suffix[metric], ""))
  ) %>%
  select(SPECIES, UA_ID, UA_Area_ha_recomp, column, value) %>%
  pivot_wider(names_from = column, values_from = value)

message(sprintf("   %d units by %d metric columns.",
                nrow(wide_raster), ncol(wide_raster) - 3))


# 5. Proximity index from the FRAGSTATS patch tables [3/5] ---------------------
# PROX_MN is the mean of the patch-level proximity index over the patches of
# a class. For the farming class, the pasture and agriculture patches are
# pooled before the mean is taken; the search neighborhoods themselves were
# evaluated on the map with the two classes separate. This operational
# definition, stated in the supplementary material, is applied identically to
# all 67 units.
message("\n[3/5] Recovering PROX_MN from the FRAGSTATS patch tables...")

patch_files <- list.files(fragstats_path, pattern = "Patch\\.CSV$",
                          ignore.case = TRUE, full.names = TRUE)
if (length(patch_files) == 0) {
  stop("No *_Patch.CSV found in ", fragstats_path)
}

read_fragstats <- function(f) {
  readr::read_csv(f, col_types = readr::cols(.default = "c"),
                  name_repair = "unique_quiet") %>%
    rename_with(str_squish) %>%
    mutate(across(where(is.character), str_squish))
}

prox_table <- patch_files %>%
  map_dfr(read_fragstats) %>%
  select(any_of(c("SPECIES", "UA", "TYPE", "PROX"))) %>%
  mutate(
    SPECIES = canonical_taxon(SPECIES),
    UA_ID   = as.integer(UA),
    PROX    = as.numeric(PROX),
    class_name = dplyr::case_when(
      TYPE == "cls_1"               ~ "Forest",
      TYPE == "cls_2"               ~ "Herbaceous",
      TYPE %in% c("cls_3", "cls_4") ~ "Agropecuaria",
      TYPE == "cls_5"               ~ "Water",
      TYPE == "cls_6"               ~ "Non_Vegetated",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SPECIES), !is.na(UA_ID), !is.na(class_name), !is.na(PROX)) %>%
  group_by(SPECIES, UA_ID, class_name) %>%
  summarise(PROX_MN = mean(PROX), .groups = "drop") %>%
  mutate(column = paste0("PROX_MN_", class_name)) %>%
  select(SPECIES, UA_ID, column, PROX_MN) %>%
  pivot_wider(names_from = column, values_from = PROX_MN)

message(sprintf("   PROX_MN recovered for %d units.", nrow(prox_table)))

recomputed <- wide_raster %>%
  left_join(prox_table, by = c("SPECIES", "UA_ID")) %>%
  arrange(SPECIES, UA_ID)

if (nrow(recomputed) != N_UNITS_EXPECTED) {
  stop(sprintf("Recomputed table has %d rows, expected %d.",
               nrow(recomputed), N_UNITS_EXPECTED))
}

write_csv(recomputed, file.path(out_path, "Tabela_Metricas_Recomputada.csv"))
message(sprintf("   Tabela_Metricas_Recomputada.csv exported (%d units).",
                nrow(recomputed)))


# 6. Agreement with the stored values [4/5] ------------------------------------
# Compares every recomputed cell with the stored cell of Data_Raw_FINAL.csv,
# wherever a stored value exists. Cells that are empty in the stored file are
# the target of the repair and are listed in Section 7.
message("\n[4/5] Comparing recomputed against stored values...")

stored <- read_csv(file.path(out_path, "Data_Raw_FINAL.csv"), show_col_types = FALSE) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

shared_cols <- intersect(
  setdiff(names(recomputed), c("SPECIES", "UA_ID", "UA_Area_ha_recomp")),
  names(stored)
)
message(sprintf("   %d metric columns present in both tables.", length(shared_cols)))

comparison <- shared_cols %>%
  map_dfr(function(cl) {
    j <- stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(!is.na(stored_value), !is.na(new_value))

    if (nrow(j) < 3) {
      return(tibble(Column = cl, N = nrow(j), Pearson_r = NA_real_,
                    Mean_Bias = NA_real_, Median_Rel_Error = NA_real_,
                    Max_Rel_Error = NA_real_))
    }

    rel <- abs(j$new_value - j$stored_value) /
      ifelse(abs(j$stored_value) < 1e-9, NA_real_, abs(j$stored_value))

    tibble(
      Column           = cl,
      N                = nrow(j),
      Pearson_r        = suppressWarnings(cor(j$stored_value, j$new_value)),
      Mean_Bias        = mean(j$new_value - j$stored_value),
      Median_Rel_Error = median(rel, na.rm = TRUE),
      Max_Rel_Error    = suppressWarnings(max(rel, na.rm = TRUE))
    )
  }) %>%
  mutate(
    Verdict = dplyr::case_when(
      is.na(Median_Rel_Error)     ~ "too few pairs to judge",
      Median_Rel_Error < 0.01     ~ "reproduces the stored pipeline; fill gaps only",
      Median_Rel_Error < 0.05     ~ "small systematic offset; inspect before use",
      TRUE                        ~ "not interchangeable; replace the whole column"
    )
  ) %>%
  arrange(desc(Median_Rel_Error))

write_csv(comparison, file.path(report_path, "TableS_Recompute_Agreement.csv"))

message("\n   Agreement, worst ten metrics first:")
print(as.data.frame(comparison %>% slice_head(n = 10)), row.names = FALSE, digits = 4)

n_bad <- sum(comparison$Median_Rel_Error >= 0.01, na.rm = TRUE)
if (n_bad > 0) {
  message(sprintf(
    "\n   WARNING. %d column(s) do not reproduce the stored values to within 1%%.",
    n_bad))
  message("   For those columns, replace the values of ALL 67 units with the")
  message("   recomputed ones rather than filling only the two gaps, so that a")
  message("   single measurement procedure underlies the column, and record the")
  message("   change in the Methods.")
} else {
  message("\n   All shared columns reproduce the stored values to within 1%.")
  message("   Filling only the empty cells is therefore methodologically safe.")
}


# 7. Cells that script 01 will fill [5/5] --------------------------------------
message("\n[5/5] Listing the cells that Part I will fill...")

repaired <- shared_cols %>%
  map_dfr(function(cl) {
    stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(is.na(stored_value), !is.na(new_value)) %>%
      transmute(SPECIES, UA_ID, Column = cl, Filled_value = new_value)
  }) %>%
  arrange(SPECIES, UA_ID, Column)

write_csv(repaired, file.path(report_path, "TableS_Repaired_Cells.csv"))

message(sprintf("   %d empty cells will be filled, in %d sampling unit(s):",
                nrow(repaired), n_distinct(paste(repaired$SPECIES, repaired$UA_ID))))
print(as.data.frame(repaired %>% count(SPECIES, UA_ID, name = "cells_filled")),
      row.names = FALSE)

still_empty <- shared_cols %>%
  map_dfr(function(cl) {
    stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(is.na(stored_value), is.na(new_value)) %>%
      transmute(SPECIES, UA_ID, Column = cl)
  })

if (nrow(still_empty) > 0) {
  message(sprintf(
    "\n   %d cell(s) remain empty after recomputation. These are genuinely",
    nrow(still_empty)))
  message("   undefined, not missing: a per-patch mean over a class that has no")
  message("   patch in the unit, or ENN over a single patch. They are structural")
  message("   zeros only for CA, PLAND, NP, PD, ED and LPI, and Part I recodes")
  message("   them accordingly. They are NOT zeros for AREA_MN, FRAC_MN, ENN_MN")
  message("   or PROX_MN.")
  print(as.data.frame(still_empty %>% count(Column, name = "n_units")),
        row.names = FALSE)
}

message("\n============================================================")
message(sprintf("  Part 0 complete. %d units recomputed.", nrow(recomputed)))
if (!RECOMPUTE_PATCH_MEANS) {
  message("  AREA_MN, FRAC_MN and ENN_MN were not recomputed. They are complete in")
  message("  Data_Raw_FINAL.csv for all 67 units and need no repair. To extend the")
  message("  agreement table to them, set RECOMPUTE_PATCH_MEANS <- TRUE and rerun;")
  message("  the cache keeps the two metric sets separate, so nothing is lost.")
}
message("  Next: run Part I, which reads Tabela_Metricas_Recomputada.csv and")
message("  repairs Data_Raw_FINAL.csv before any analysis.")
message("============================================================\n")


# 8. Session information -------------------------------------------------------
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(report_path, "session_info_part00.txt"))
